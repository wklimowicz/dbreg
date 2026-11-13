# Project each slope and the score correction onto the same weighted FE space.
# The correction response is the score divided by information weight, so
# subtracting two absolute working-response solutions is unnecessary.
dbglm_fe_refine = function(system, rhs, solver, scales, tol, max_iter) {
  if (length(scales) != ncol(rhs) || any(!is.finite(scales)) || any(scales <= 0)) {
    projection_abort("singular", "Invalid GLM working-column norms.")
  }
  projection = solver[["solve"]](rhs)
  diagonal = sqrt(system[["diag"]])
  previous = Inf
  floor = FALSE
  refinements = min(8L, max_iter)
  for (refine in seq_len(refinements)) {
    residual = rhs - fe_multiply(system, projection)
    # Both denominators are weighted norms: changing a working column's units
    # or replicating the weighted sample cannot change the refinement decision.
    error = max(abs(sweep(residual / diagonal, 2L, scales, "/")))
    if (is.finite(error) && error <= tol / 8) break
    if (is.finite(error) && error <= 1024 * .Machine[["double.eps"]] && error > previous / 2) {
      floor = TRUE
      break
    }
    if (refine == refinements) {
      projection_abort("nonconverged", "GLM fixed-effect moment refinement did not converge.",
        list(residual = error, target = tol / 8, units = "weighted-column norm"))
    }
    previous = error
    projection = projection + solver[["solve"]](residual)
  }
  list(projection = projection, refinement = list(residual = error, target = tol / 8,
    corrections = refine - 1L, roundoff_floor = floor))
}

# Uncertain subtraction must not decide structural rank. The score is always
# evaluated on the actual projected rows; only the structural Gram uses FWL.
dbglm_fe_information = function(conn, source, columns, information, budget, limits, nobs) {
  raw = information[["cross"]][-1L, -1L, drop = FALSE]
  H = information[["within"]][-1L, -1L, drop = FALSE]
  recompute = any(!is.finite(H))
  if (length(H) && !recompute) {
    diagonal = diag(H)
    recompute = any(diagonal <= 0) || any(diagonal < 1e-8 * diag(raw))
    if (!recompute) {
      units = outer(sqrt(diagonal), sqrt(diagonal))
      # A small eigenvalue is uncertain relative to the operands of the
      # subtraction, not only to its result. This conservative heuristic is
      # not a total forward-error bound. Recompute before deciding rank.
      normalized = H / units
      uncertainty = 64 * .Machine[["double.eps"]] *
        norm((abs(raw) + abs(raw - H)) / units, type = "I")
      recompute = any(!is.finite(units) | units <= 0) ||
        any(!is.finite(normalized)) || !is.finite(uncertainty)
      if (!recompute) {
        values = eigen(normalized, symmetric = TRUE, only.values = TRUE)[["values"]]
        recompute = min(values) <= max(1e-8 * max(values), uncertainty)
      }
    }
  }
  if (recompute) {
    actual = fe_information_query(conn, source, columns[-1L], budget, limits)
    if (actual[["n"]] != nobs) projection_abort("certificate_failed", "GLM projected information changed the row count.")
    H = actual[["cross"]]
  }
  cells = length(columns)
  if (budget[["cells"]] + cells > limits[["transfer_cells"]] ||
      budget[["rows"]] + 1 > limits[["summary_rows"]]) {
    projection_abort("transfer_limit", "GLM projected scores exceed the transfer budget.")
  }
  terms = if (length(columns) > 1L) paste(sql_stable_sum(conn, paste0("w*(", columns[1L], ")*(", columns[-1L], ")")),
    "AS", paste0("s", seq_len(length(columns) - 1L))) else character()
  row = DBI::dbGetQuery(conn, paste("SELECT", paste(c("COUNT(*) AS n", terms), collapse = ","), "FROM", source))
  if (nrow(row) != 1L || !isTRUE(as.numeric(row[["n"]]) == nobs) || any(!is.finite(as.numeric(row[1, ])))) {
    projection_abort("certificate_failed", "GLM projected score certificate is invalid.")
  }
  budget[["cells"]] = budget[["cells"]] + cells
  budget[["rows"]] = budget[["rows"]] + 1L
  list(H = H, score = as.numeric(row[1L, -1L, drop = FALSE]), recomputed = recompute)
}

dbglm_fe_step = function(conn, working, x, fe, sizes, nobs, tol, cluster = FALSE, cell = FALSE,
                         score_keys = NULL, max_iter = 10000L, lookup_manager = NULL, raw_state = NULL,
                         certify_weak = FALSE) {
  if (!is.null(raw_state) &&
      (!identical(raw_state[["working"]], working) ||
       !identical(raw_state[["working_columns"]], x) ||
       !identical(raw_state[["information_type"]], "fisher") ||
       !is.matrix(raw_state[["H"]]) || !is.numeric(raw_state[["H"]]) ||
       !identical(dim(raw_state[["H"]]), rep(length(x), 2L)) ||
       !is.numeric(raw_state[["score"]]) || length(raw_state[["score"]]) != length(x) ||
       any(!is.finite(c(raw_state[["H"]], raw_state[["score"]]))))) {
    projection_abort("singular", "Invalid matching GLM working moments.")
  }
  limits = projection_limits()
  width = length(x) + 1L
  if (sum(sizes) * width > limits[["design_entries"]] ||
      !fe_dense_feasible(sizes, width, limits)) {
    projection_abort("system_limit", "GLM fixed-effect information exceeds its projection budget.")
  }
  budget = list2env(list(rows = 0, cells = 0, upload_bytes = 0), parent = emptyenv())
  columns = paste0("v", seq_len(width))
  keys = unique(c(fe, score_keys))
  source_sql = paste0("SELECT w,score,", if (cell) "score_sq,", if (cluster) "__dbglm_cluster,", paste(keys, collapse = ","),
    ",score/w AS ", columns[1L],
    if (length(x)) paste0(",", paste(paste(x, "AS", columns[-1L]), collapse = ",")),
    " FROM ", working)
  source = paste0("(", source_sql, ") glm_score_coordinates")
  # Reuse this source-scale query for the final SQL certificate. Reserve its
  # scalar transfer before dispatch, including the added weighted column norms.
  bounded = is.null(raw_state) && identical(detect_backend(conn)[["name"]], "spark")
  extra_scale = if (bounded) c("MAX(ABS(w)) AS max_weight",
    paste0("SUM(CASE WHEN ", paste(paste0(c("w", columns), " IS NULL"), collapse = " OR "),
      " THEN 1 ELSE 0 END) AS missing")) else character()
  weak_scale = if (length(fe) > 1L) c("MIN(w) AS information_min_weight", "MAX(w) AS information_max_weight") else character()
  scale_cells = 2L * width + length(extra_scale) + length(weak_scale)
  if (budget[["cells"]] + scale_cells > limits[["transfer_cells"]] ||
      budget[["rows"]] + 1 > limits[["summary_rows"]]) {
    projection_abort("transfer_limit", "GLM working-column norms exceed the transfer budget.")
  }
  scale = DBI::dbGetQuery(conn, paste0("SELECT ", paste(c(
    paste0("MAX(ABS(", columns, ")) AS ", columns),
    paste0(sql_stable_sum(conn, paste0("w*", columns, "*", columns)), " AS norm", seq_len(width)), extra_scale, weak_scale),
    collapse = ","), " FROM ", source))
  if (nrow(scale) != 1L || any(!is.finite(as.numeric(scale[1, ])))) {
    projection_abort("singular", "Invalid GLM working-column norms.")
  }
  # DuckDB can attempt the bounded solve on weak rows and use its existing
  # projection certificate to decide whether SQL recovery is needed.
  if (!certify_weak && length(weak_scale) && scale[["information_min_weight"]] < 1e-8 * scale[["information_max_weight"]]) {
    projection_abort("weak_information", "Weak GLM fixed-effect information requires bounded QR recovery.")
  }
  if (bounded && (scale[["missing"]] != 0 || scale[["max_weight"]] <= 0)) {
    projection_abort("singular", "Invalid GLM working-column bounds.")
  }
  bounds = if (bounded) list(weight = as.numeric(scale[["max_weight"]]),
    columns = as.numeric(scale[1, columns, drop = FALSE])) else NULL
  squared_norms = as.numeric(scale[1, paste0("norm", seq_len(width)), drop = FALSE])
  if (any(squared_norms < 0)) projection_abort("singular", "Invalid GLM working-column norms.")
  norms = sqrt(squared_norms)
  norms[norms == 0] = 1
  budget[["cells"]] = budget[["cells"]] + scale_cells
  budget[["rows"]] = budget[["rows"]] + 1L
  summaries = query_fe_summaries(conn, source, fe, columns, sizes, budget, limits)
  info = summaries[["info"]]
  system = fe_system(info, summaries[["pairs"]], limits)
  if (is.null(system)) projection_abort("system_limit", "GLM fixed-effect system exceeds its projection budget.")
  if (any(!is.finite(system[["diag"]]) | system[["diag"]] <= 0)) {
    projection_abort("singular", "Invalid GLM fixed-effect information weights.")
  }
  cache = list2env(list(projection_limits = limits, pcg_max_iter = max_iter), parent = emptyenv())
  solver = fe_solver(system, cache)
  B = do.call(rbind, lapply(info, `[[`, "b"))
  rhs = solver[["project"]](B)
  refined = dbglm_fe_refine(system, rhs, solver, norms, tol, max_iter)
  projection = refined[["projection"]]
  global = NULL
  if (!is.null(raw_state)) {
    H = raw_state[["H"]]
    score = raw_state[["score"]]
    global = list(cross = rbind(c(squared_norms[1L], score), cbind(score, H)),
      n = nobs, sum_w = sum(info[[1L]][["h"]]))
  }
  information = fe_fwl_information(conn, source, columns, B, projection, system, budget, limits, global, bounds)
  pieces = lapply(seq_along(fe), function(i) {
    positions = seq.int(system[["starts"]][i] + 1L, system[["starts"]][i + 1L])
    projection[positions, , drop = FALSE]
  })
  lookups = fe_coefficient_lookups(conn, pieces, budget, limits, lookup_manager)
  success = FALSE
  on.exit(if (!success) for (table in rev(lookups[["tables"]])) lookup_manager[["drop"]](table))
  joins = vapply(seq_along(fe), function(i) {
    alias = paste0("p", i)
    paste0(" JOIN ", lookups[["sources"]][i], " ", alias, " ON s.", fe[i], "=", alias, ".id")
  }, character(1))
  projected = vapply(seq_len(width), function(j) {
    adjustment = paste0("p", seq_along(fe), ".c", j, collapse = "+")
    paste0("(s.", columns[j], "-(", adjustment, ")) AS ", columns[j])
  }, character(1))
  broadcast = if (length(lookups[["tables"]])) paste0("/*+ BROADCAST(",
    paste0("p", seq_along(fe), collapse = ","), ") */ ") else ""
  residual_source = paste0("(SELECT ", broadcast, "s.w,s.score,", if (cell) "s.score_sq,", if (cluster) "s.__dbglm_cluster,", paste0("s.", keys, collapse = ","), ",",
    paste(projected, collapse = ","), " FROM (", source_sql, ") s",
    paste(joins, collapse = ""), ") glm_within_coordinates")
  if (information[["n"]] != nobs) {
    projection_abort("certificate_failed", "GLM fixed-effect projection changed the row count.")
  }
  within = dbglm_fe_information(conn, residual_source, columns, information, budget, limits, nobs)
  certificate = fe_projection_certificate(conn, residual_source, fe, columns, nobs,
    as.numeric(scale[1, columns, drop = FALSE]),
    grouped = c(group_score = dbglm_fe_group_score_sql(conn)))
  if (!is.finite(certificate[["max_scaled"]]) || certificate[["max_scaled"]] > tol ||
      !is.finite(certificate[["max_relative"]]) || certificate[["max_relative"]] > tol) {
    projection_abort("certificate_failed", "GLM fixed-effect working rows failed the projection certificate.",
      c(certificate, list(target = tol)))
  }
  result = list(H = within[["H"]], score = within[["score"]],
    projection = pieces, residual_source = residual_source, certificate = certificate,
    fe_group_score = if (isTRUE(is.finite(certificate[["grouped"]][["group_score"]]))) certificate[["grouped"]][["group_score"]],
    transfer_cells = budget[["cells"]], upload_bytes = budget[["upload_bytes"]],
    refinement = refined[["refinement"]], lookup_tables = lookups[["tables"]],
    information_recomputed = within[["recomputed"]], reused_working_moments = !is.null(raw_state))
  success = TRUE
  result
}
