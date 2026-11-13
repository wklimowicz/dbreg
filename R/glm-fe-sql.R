# Keep coefficient state and projected rows in owned SQL tables. Ordinary CG
# transfers scalar diagnostics; weak systems use bounded exact working cells
# for QR recovery, then upload residuals and coefficients for certification.
dbglm_fe_sql_step = function(conn, working, x, fe, nobs, tol, manager,
                               cluster, cell, control, ops = NULL, score_keys = NULL, warm = NULL) {
  limits = projection_limits()
  budget = list2env(list(rows = 0, cells = 0, upload_bytes = 0), parent = emptyenv())
  # Projection retries must release only tables created by this attempt.
  # The caller's coefficient state and immutable source remain borrowed.
  attempt = list2env(list(tables = character()), parent = emptyenv())
  success = FALSE
  create = function(sql) {
    table = manager[["create"]](sql)
    attempt[["tables"]] = c(attempt[["tables"]], table)
    table
  }
  drop = function(table) {
    manager[["drop"]](table)
    attempt[["tables"]] = setdiff(attempt[["tables"]], table)
  }
  # CG overwrites its own working states in place where the manager can.
  replace = if (!is.null(manager[["replace"]])) function(table, sql) {
    value = manager[["replace"]](table, sql)
    attempt[["tables"]] = c(setdiff(attempt[["tables"]], table), value)
    value
  }
  on.exit(if (!success) for (table in rev(attempt[["tables"]])) drop(table))
  columns = paste0("v", seq_len(length(x) + 1L))
  varying = !is.null(ops)
  if (!varying) ops = lapply(fe, function(key) list(f = key, v = "1.0"))
  axes = setdiff(unique(vapply(ops, `[[`, character(1), "v")), "1.0")
  ids = unique(c(fe, axes, "score", if (cluster) "__dbglm_cluster", if (cell) "score_sq", score_keys))
  # The coordinate and normalization steps stay inline: each statement that
  # creates a table carries a fixed cost on remote engines.
  input = paste0("(SELECT ", paste(c(ids, "w", "score/w AS v1",
    if (length(x)) paste(x, "AS", columns[-1L])), collapse = ","), " FROM ", working, ") glm_input")
  stats = fe_cg_fetch(conn, paste("SELECT",
    paste(c(paste0("MAX(ABS(", columns, ")) AS ", columns),
      "MIN(w) AS minimum_weight", "MAX(w) AS maximum_weight"), collapse = ","), "FROM", input),
    length(columns) + 2L, budget, limits)
  scales = pmax(as.numeric(stats[1, columns, drop = FALSE]), .Machine[["double.xmin"]])
  weak = stats[["minimum_weight"]] < 1e-8 * stats[["maximum_weight"]]
  if (any(!is.finite(scales))) projection_abort("singular", "Invalid GLM SQL projection scales.")
  normalized = paste("SELECT", paste(c(ids, "w",
    paste0(columns, "/(", dbglm_number(scales), ") AS ", columns)), collapse = ","),
    "FROM", input)
  # Structural columns keep their units across IRLS iterations, so their
  # earlier FE coefficients seed CG. The score column always starts at zero.
  start = if (!is.null(warm) && length(x) && length(warm[["scales"]]) == length(columns)) {
    list(table = warm[["coefficients"]], ratio = c(NA, warm[["scales"]][-1L] / scales[-1L]))
  }
  # Exact QR recovery when its dense design fits the bounded budget, else NULL
  # after releasing any table the declined attempt created.
  bounded_qr = function(source) {
    before = attempt[["tables"]]
    tryCatch(dbglm_fe_qr_projection(conn, source, columns, ids, ops, create, drop, budget, limits),
      dbreg_projection_capacity = function(e) {
        for (table in rev(setdiff(attempt[["tables"]], before))) drop(table)
        NULL
      })
  }
  # Window and CG projections return rows already multiplied back by scales.
  exact_unavailable = FALSE
  if (!varying && length(fe) == 1L) {
    result = dbglm_fe_window_projection(conn, paste0("(", normalized, ") glm_normalized"),
      columns, ids, fe, create, drop, unscale = scales)
    rows = result[["table"]]
  } else {
    normalized = create(normalized)
    # Weak weights prefer exact QR recovery. When its dense design cannot fit
    # (large FE systems), CG runs instead under the same certificate below.
    result = if (weak) bounded_qr(normalized)
    exact = !is.null(result)
    exact_unavailable = weak && !exact
    if (!exact) result = fe_sql_cg(conn, normalized, columns, ids, ops, rep(1, length(columns)),
      tol / 16, control, create, drop, budget, limits, retain_coefficients = TRUE, allow_zero = varying,
      replace = replace, start = start, output_scales = scales)
    drop(normalized)
    rows = if (exact) {
      scaled = create(paste("SELECT", paste(c(ids, "w",
        paste0(columns, "*(", dbglm_number(scales), ") AS ", columns)), collapse = ","),
        "FROM", result[["table"]]))
      drop(result[["table"]])
      scaled
    } else result[["table"]]
  }
  for (table in result[["intermediates"]]) drop(table)
  # Plain FE certificates also return the working rows' group scores and the
  # global information moments, saving two scalar queries per projection.
  certify = function(rows) if (varying) fe_coordinate_certificate(conn, rows, ops, columns, nobs, scales) else
    fe_projection_certificate(conn, rows, fe, columns, nobs, scales,
      grouped = c(group_score = dbglm_fe_group_score_sql(conn)), global = fe_information_terms(conn, columns))
  certificate = certify(rows)
  valid_certificate = function(value) is.finite(value[["max_scaled"]]) && value[["max_scaled"]] <= tol &&
    is.finite(value[["max_relative"]]) && value[["max_relative"]] <= tol
  remaining = control[["fe_max_iter"]] - if (length(result[["cycles"]])) result[["cycles"]] else 0L
  for (refinement in seq_len(8L)) {
    if (valid_certificate(certificate) || remaining <= 0L) break
    # Keep the QR residual itself. Reconstructing it from nearly equal raw
    # coordinates would reintroduce the roundoff this refinement removes.
    correction = if (!exact_unavailable) bounded_qr(rows)
    unit = scales
    if (is.null(correction)) {
      # Large FE systems cannot hold the dense QR design. Refine the retained
      # residual by CG instead, in normalized units, with per-column targets
      # relative to each residual RMS (the certificate's relative criterion).
      exact_unavailable = TRUE
      rms = as.numeric(fe_cg_fetch(conn, paste("SELECT", paste0("SQRT(",
        sql_sum(conn, paste0("w*POWER(", columns, "/(", dbglm_number(scales), "),2)")), "/",
        sql_sum(conn, "w"), ") AS r", seq_along(columns), collapse = ","), "FROM", rows),
        length(columns), budget, limits)[1, ])
      if (any(!is.finite(rms))) projection_abort("singular", "Invalid GLM SQL residual scales.")
      relative = ifelse(rms > 128 * .Machine[["double.eps"]], pmin(1, rms), 1)
      normal = create(paste("SELECT", paste(c(ids, "w",
        paste0(columns, "/(", dbglm_number(scales), ") AS ", columns)), collapse = ","), "FROM", rows))
      refine_control = control
      refine_control[["fe_max_iter"]] = remaining
      correction = fe_sql_cg(conn, normal, columns, ids, ops, relative, tol / 16, refine_control,
        create, drop, budget, limits, retain_coefficients = TRUE, allow_zero = varying,
        replace = replace, output_scales = scales)
      drop(normal)
      unit = rep(1, length(columns))
    }
    adjusted = create(paste("SELECT c.__cg_op,c.__cg_group,",
      paste0("c.c", seq_along(columns), "+COALESCE(d.c", seq_along(columns), ",0.0)/(",
        dbglm_number(unit), ") AS c", seq_along(columns), collapse = ","),
      "FROM (SELECT * FROM", result[["coefficients"]], ") c LEFT JOIN", correction[["coefficients"]],
      "d USING (__cg_op,__cg_group)"))
    drop(rows)
    drop(if (is.null(result[["coefficient_table"]])) result[["coefficients"]] else result[["coefficient_table"]])
    result[["coefficient_table"]] = NULL
    drop(correction[["coefficients"]])
    rows = correction[["table"]]
    result[["coefficients"]] = adjusted
    if (identical(correction[["method"]], "bounded_qr")) {
      result[["method"]] = "bounded_qr"
      remaining = remaining - 1L
    } else remaining = remaining - correction[["cycles"]]
    certificate = certify(rows)
  }
  if (!valid_certificate(certificate)) {
    if (remaining <= 0L) projection_abort("nonconverged", "GLM SQL projection did not converge within its iteration cap.")
    projection_abort("certificate_failed", "GLM SQL working rows failed the projection certificate.")
  }
  information = if (is.null(certificate[["global"]])) fe_information_query(conn, rows, columns, budget, limits) else {
    fe_information_reserve(columns, budget, limits)
    fe_information_read(certificate[["global"]], columns, budget)
  }
  if (information[["n"]] != nobs) {
    projection_abort("certificate_failed", "GLM SQL projection changed the row count.")
  }
  M = information[["cross"]]
  result = list(H = M[-1L, -1L, drop = FALSE], score = M[-1L, 1L],
    coefficients = result[["coefficients"]], coefficient_table = result[["coefficient_table"]], scales = scales,
    residual_source = rows, certificate = certificate,
    fe_group_score = if (isTRUE(is.finite(certificate[["grouped"]][["group_score"]]))) certificate[["grouped"]][["group_score"]],
    method = if (!varying && length(fe) == 1L) "sql_window" else
      if (identical(result[["method"]], "bounded_qr")) "sql_qr" else "sql_cg")
  success = TRUE
  result
}

# Release a projection's owned tables. Inline coefficients are an aggregate
# over `coefficient_table`, which stays owned while they are in use.
dbglm_fe_sql_release = function(manager, projection) {
  if (!is.null(projection[["residual_source"]])) manager[["drop"]](projection[["residual_source"]])
  owner = projection[["coefficient_table"]]
  manager[["drop"]](if (is.null(owner)) projection[["coefficients"]] else owner)
}

dbglm_fe_sql_start = function(base, fe, alpha, manager, ops = NULL) {
  keys = if (is.null(ops)) fe else vapply(ops, `[[`, character(1), "f")
  sql = vapply(seq_along(keys), function(i) paste("SELECT DISTINCT", i,
    "AS __cg_op,", keys[i], "AS __cg_group, CAST(", dbglm_number(if (i == 1L) alpha else 0),
    "AS DOUBLE) AS c1 FROM", base), character(1))
  manager[["create"]](paste(sql, collapse = " UNION ALL "))
}

# Dampen toward the same constant predictor used by public initialization.
# No fitted rows cross the client boundary, and the caller releases the old state.
dbglm_fe_sql_damp = function(effects, intercept, manager) {
  manager[["create"]](paste("SELECT __cg_op,__cg_group,c1/2.0 + CASE WHEN __cg_op=1 THEN",
    dbglm_number(intercept / 2), "ELSE 0.0 END AS c1 FROM", effects))
}

# With `inline`, return an unaliased subquery instead of an owned table: the
# step trials re-read the same expression without a separate CREATE.
dbglm_fe_sql_delta = function(direction, correction, manager, column = "delta", inline = FALSE) {
  scales = direction[["scales"]]
  terms = c(paste0("c1*(", dbglm_number(scales[1L]), ")"),
    if (length(correction)) paste0("- c", seq_along(correction) + 1L, "*(",
      dbglm_number(scales[-1L] * correction), ")"))
  sql = paste("SELECT __cg_op,__cg_group,",
    paste(terms, collapse = " "), "AS", column, "FROM", direction[["coefficients"]])
  if (inline) paste0("(", sql, ")") else manager[["create"]](sql)
}

dbglm_fe_sql_update = function(effects, delta, step, manager) {
  manager[["create"]](paste0("SELECT __cg_op,__cg_group,e.c1+(",
    dbglm_number(step), ")*COALESCE(d.delta,0) AS c1 FROM ", effects, " e LEFT JOIN ", delta,
    " d USING (__cg_op,__cg_group)"))
}

dbglm_varying_basis = function(conn, sample, manager) {
  create = manager[["create"]]
  drop = manager[["drop"]]
  base = sample[["base"]]
  trends = sample[["trend_columns"]]
  fields = names(DBI::dbGetQuery(conn, paste("SELECT * FROM", base, "LIMIT 0")))
  # Positive prior weights define a fixed coordinate change. Iteration-specific
  # working weights still define the projections and score certificates.
  weighted = create(paste("SELECT *,__dbglm_a AS w FROM", base))
  row = DBI::dbGetQuery(conn, paste("SELECT",
    paste(sprintf("MAX(ABS(%s)) AS %s", trends, trends), collapse = ","), "FROM", weighted))
  scales = setNames(pmax(as.numeric(row[1, ]), .Machine[["double.xmin"]]), names(trends))
  if (any(!is.finite(scales))) projection_abort("singular", "Invalid GLM varying trend scales.")
  basis = fe_varying_basis(conn, weighted, sample[["fe"]], sample[["slopes"]],
    trends, scales, c(fields, "w"), create, drop)
  derived = basis[["table"]]
  basis[["table"]] = create(paste("SELECT", paste(c(fields, basis[["basis"]]), collapse = ","),
    "FROM", basis[["table"]]))
  # The coordinate builder owns its derived intermediate, not the source.
  drop(derived)
  drop(weighted)
  basis
}

# A single additive group has an exact weighted-mean projection. Keep its
# coefficients and residual coordinates in SQL without a client dictionary
# or an iterative coefficient solve; use the same certificate as general CG.
dbglm_fe_window_projection = function(conn, source, columns, ids, fe, create, drop, unscale = NULL) {
  mean = function(column, window = NULL) paste0("(",
    sql_stable_sum(conn, paste0("w*", column), window = window), ")/(",
    sql_stable_sum(conn, "w", window = window), ")")
  window = paste("PARTITION BY", fe)
  averages = vapply(columns, mean, character(1), window)
  saved = safe_identifiers(paste0("dbglm_window_mean_", seq_along(columns)), taken = c(ids, "w", columns))
  first = create(paste("SELECT", paste(c(ids, "w", paste(averages, "AS", saved),
    paste0(columns, "-(", averages, ") AS ", columns)), collapse = ","), "FROM", source))
  # Subtract the remainder's mean as well. A first subtraction can leave
  # roundoff larger than the true within variation in nearly absorbed cells.
  residual = paste0(columns, "-(", vapply(columns, mean, character(1), window), ")")
  # Optionally multiply back by the caller's column scales in the same pass.
  residual = if (length(unscale)) paste0("(", residual, ")*(", dbglm_number(unscale), ") AS ", columns) else
    paste(residual, "AS", columns)
  table = create(paste("SELECT", paste(c(ids, "w", residual), collapse = ","), "FROM", first))
  # Materialize the compact G-row coefficients table and drop `first`
  # immediately so step guards and halving trials do not re-aggregate all rows.
  coefficients = create(paste("SELECT 1 AS __cg_op,", fe, "AS __cg_group,",
    paste(paste0("MAX(", saved, ") + (", vapply(columns, mean, character(1)), ") AS c", seq_along(columns)), collapse = ","),
    "FROM", first, "GROUP BY", fe))
  drop(first)
  list(table = table, coefficients = coefficients, coefficient_table = NULL)
}

# A weak working metric can hide FE directions in normal equations or a CG
# stopping floor. Recover all directions from a bounded, exact compression of
# the actual working coordinates. Householder QR never forms cross-products or
# a complete Q; its residual is retained separately from rounded coefficients.
dbglm_fe_qr_projection = function(conn, source, columns, ids, ops, create, drop, budget, limits) {
  families = unique(vapply(ops, `[[`, character(1), "f"))
  axes = setdiff(unique(vapply(ops, `[[`, character(1), "v")), "1.0")
  keys = unique(c(families, axes, columns))
  sizes = fe_cg_fetch(conn, paste("SELECT",
    paste0("COUNT(DISTINCT ", families, ") AS n", seq_along(families), collapse = ","),
    "FROM", source), length(families), budget, limits)
  sizes = vapply(sizes, function(value) as.numeric(value)[1L], numeric(1), USE.NAMES = FALSE)[
    match(vapply(ops, `[[`, character(1), "f"), families)]
  width = sum(sizes)
  fields = length(keys) + 2L
  # Include the coefficient state and upload buffers separately: a redundant
  # incidence design can have more coordinate columns than exact working cells.
  state_entries = width * length(columns)
  available_bytes = limits[["memory_bytes"]] - 32 * state_entries
  cap = floor(min(limits[["summary_rows"]] - budget[["rows"]],
    (limits[["transfer_cells"]] - budget[["cells"]]) / fields,
    limits[["design_entries"]] / max(1, width),
    available_bytes / (64 * max(1, width + length(columns)))))
  unavailable = function() projection_abort("system_limit",
    paste("GLM fixed-effect projection could not be certified: exact QR recovery",
      "exceeds the bounded QR recovery budget (weak or near-separated fixed-effect information)."))
  if (!is.finite(width) || width < 1 || state_entries > limits[["design_entries"]] ||
      cap < 1 || max(sizes) > cap) unavailable()
  count = as.numeric(fe_cg_fetch(conn, paste0("SELECT COUNT(*) AS n FROM (SELECT DISTINCT ",
    paste(keys, collapse = ","), " FROM ", source, " LIMIT ",
    format(cap + 1, scientific = FALSE, trim = TRUE), ") glm_qr_keys"), 1L, budget, limits)[["n"]])
  if (!is.finite(count) || count > cap || count < 1) unavailable()
  # Account for the count query before reserving the bounded cell transfer.
  if (budget[["rows"]] + count > limits[["summary_rows"]] ||
      budget[["cells"]] + count * fields > limits[["transfer_cells"]]) unavailable()
  cells = create(paste0("SELECT ROW_NUMBER() OVER (ORDER BY ", paste(keys, collapse = ","),
    ") AS __dbglm_qr_cell, q.* FROM (SELECT ", paste(keys, collapse = ","), ",",
    sql_stable_sum(conn, "w"), " AS w FROM ", source, " GROUP BY ",
    paste(keys, collapse = ","), ") q"))
  values = DBI::dbGetQuery(conn, paste("SELECT * FROM", cells, "ORDER BY __dbglm_qr_cell"))
  budget[["rows"]] = budget[["rows"]] + nrow(values)
  budget[["cells"]] = budget[["cells"]] + nrow(values) * ncol(values)
  # These are bounded generated ranks and numeric working coordinates, not
  # public keys. Native drivers may return their BIGINT ranks as decimal text.
  values[] = lapply(values, as.numeric)
  if (nrow(values) != count || ncol(values) != fields ||
      any(!is.finite(as.matrix(values))) || any(values[["w"]] <= 0)) {
    projection_abort("certificate_failed", "Invalid bounded GLM working cells.")
  }
  A = matrix(0, nrow(values), width)
  domains = lapply(ops, function(op) sort(unique(as.numeric(values[[op[["f"]]]]))))
  starts = c(0, cumsum(lengths(domains)))
  if (!identical(as.numeric(lengths(domains)), sizes)) {
    projection_abort("certificate_failed", "GLM QR coordinates changed their group counts.")
  }
  for (i in seq_along(ops)) {
    group = match(as.numeric(values[[ops[[i]][["f"]]]]), domains[[i]])
    axis = if (ops[[i]][["v"]] == "1.0") rep(1, nrow(values)) else values[[ops[[i]][["v"]]]]
    A[cbind(seq_len(nrow(values)), starts[i] + group)] = axis
  }
  # Select the structural basis independently of the changing working weights.
  # Integer incidence columns have their exact gauge dependencies here; varying
  # axes have already been normalized by the prior-weight coordinate builder.
  structural = qr(A, tol = 1e-12, LAPACK = FALSE)
  if (structural[["rank"]] < 1L) projection_abort("singular", "GLM QR design has no identified coordinates.")
  keep = structural[["pivot"]][seq_len(structural[["rank"]])]
  root = sqrt(values[["w"]])
  factor = qr(A[, keep, drop = FALSE] * root, LAPACK = TRUE)
  response = as.matrix(values[columns]) * root
  solved = qr.coef(factor, response)
  coefficients = matrix(0, width, length(columns))
  coefficients[keep, ] = solved
  transformed = qr.qty(factor, response)
  transformed[seq_along(keep), ] = 0
  residual = qr.qy(factor, transformed) / root
  if (any(!is.finite(coefficients)) || any(!is.finite(residual))) {
    projection_abort("singular", "Bounded GLM QR lost numerical information.")
  }
  pieces = vapply(seq_along(ops), function(i) {
    at = seq.int(starts[i] + 1L, starts[i + 1L])
    lookup = fe_coefficient_sql(conn, cbind(domains[[i]], coefficients[at, , drop = FALSE]), budget, limits)
    paste0("SELECT ", i, " AS __cg_op,c1 AS __cg_group,",
      paste0("c", seq_along(columns) + 1L, " AS c", seq_along(columns), collapse = ","),
      " FROM (", lookup, ") qr_coefficients")
  }, character(1))
  state = create(paste(pieces, collapse = " UNION ALL "))
  lookup = fe_coefficient_sql(conn, residual, budget, limits)
  joined = paste0("s.", keys, "=c.", keys, collapse = " AND ")
  rows = create(paste0("SELECT ", paste(c(paste0("s.", ids), "s.w",
    paste0("r.c", seq_along(columns), " AS ", columns)), collapse = ","), " FROM ", source,
    " s JOIN ", cells, " c ON ", joined, " JOIN (", lookup, ") r ON c.__dbglm_qr_cell=r.id"))
  drop(cells)
  list(table = rows, coefficients = state, cycles = 0L, method = "bounded_qr")
}
