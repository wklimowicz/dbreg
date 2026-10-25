# Bounded direct nuisance solves for moderate FE dimensions. The database sends
# only FE cross-products and FE-by-design moments; no observation vectors enter R.
# The largest FE block is diagonal and is absorbed analytically, so only the
# remaining levels form a dense system (Schur complement).
fe_dense_feasible = function(op_sizes, columns = 1L, limits = projection_limits()) {
  op_sizes = as.numeric(op_sizes)
  if (!length(op_sizes) || any(!is.finite(op_sizes))) return(FALSE)
  rest = sum(op_sizes) - max(op_sizes)
  rest <= limits[["dense_levels"]] &&
    max(op_sizes) * rest <= limits[["cross_entries"]] &&
    projection_memory(op_sizes, columns) <= limits[["memory_bytes"]]
}

# info: per-op list(id, h); pairs: list of list(i, j, ir, jr, h) with local
# level indices into info[[i]]$id and info[[j]]$id.
fe_system = function(info, pairs) {
  sizes = vapply(info, function(x) length(x[["id"]]), integer(1))
  if (!fe_dense_feasible(sizes)) return(NULL)
  starts = cumsum(c(0L, sizes))
  n = starts[length(info) + 1L]
  a = which.max(sizes)
  big = seq.int(starts[a] + 1L, starts[a + 1L])
  rest = setdiff(seq_len(n), big)
  pos = integer(n)
  pos[rest] = seq_along(rest)
  diagonal = as.numeric(unlist(lapply(info, `[[`, "h"), use.names = FALSE))
  Cm = matrix(0, length(big), length(rest))
  E = matrix(0, length(rest), length(rest))
  E[cbind(seq_along(rest), seq_along(rest))] = diagonal[rest]
  for (p in pairs) {
    gi = starts[p[["i"]]] + p[["ir"]]
    gj = starts[p[["j"]]] + p[["jr"]]
    if (p[["i"]] == a) {
      Cm[cbind(gi - starts[a], pos[gj])] = p[["h"]]
    } else if (p[["j"]] == a) {
      Cm[cbind(gj - starts[a], pos[gi])] = p[["h"]]
    } else {
      E[cbind(pos[gi], pos[gj])] = p[["h"]]
      E[cbind(pos[gj], pos[gi])] = p[["h"]]
    }
  }
  list(n = n, starts = starts, big = big, rest = rest, Dg = diagonal[big], Cm = Cm, E = E, diag = diagonal)
}

fe_multiply = function(sys, X) {
  x1 = X[sys[["big"]], , drop = FALSE]
  x2 = X[sys[["rest"]], , drop = FALSE]
  out = matrix(0, sys[["n"]], ncol(X))
  out[sys[["big"]], ] = sys[["Dg"]] * x1 + sys[["Cm"]] %*% x2
  out[sys[["rest"]], ] = crossprod(sys[["Cm"]], x1) + sys[["E"]] %*% x2
  out
}

# Factor the Schur complement on the non-absorbed levels; the identified basis
# is chosen once per fit so FE normalisations stay fixed across iterations.
fe_solver = function(sys, cache, retain_null = FALSE) {
  if (any(!is.finite(sys[["Dg"]])) || any(sys[["Dg"]] <= 0)) projection_abort("singular", "Fixed-effect information is singular.")
  scaled = sys[["Cm"]] / sys[["Dg"]]
  S = sys[["E"]] - crossprod(sys[["Cm"]], scaled)
  sdiag = sqrt(pmax(diag(S), 0))
  if (is.null(cache[["moment_basis"]])) {
    positive = which(sdiag > 0 & is.finite(sdiag))
    keep = integer()
    if (length(positive)) {
      C = S[positive, positive, drop = FALSE] / outer(sdiag[positive], sdiag[positive])
      pivoted = suppressWarnings(chol(C, pivot = TRUE, tol = 1e-11))
      keep = positive[attr(pivoted, "pivot")[seq_len(attr(pivoted, "rank"))]]
    }
    if (!length(keep) && !length(sys[["big"]])) stop("No identified fixed-effect span remains.", call. = FALSE)
    cache[["moment_basis"]] = keep
  }
  keep = cache[["moment_basis"]]
  scale = sdiag[keep]
  factor = if (length(keep)) tryCatch(chol(S[keep, keep, drop = FALSE] / outer(scale, scale)), error = function(e) NULL) else matrix(0, 0, 0)
  if (is.null(factor)) projection_abort("singular", "Fixed-effect information is singular.")
  # B lies in range(H) exactly, but plain SQL sums (e.g. Spark) leave roundoff
  # in the redundant FE directions that no solve can remove. Null vectors of
  # H follow from the dropped Schur pivots: v2[keep] = -S_kk^{-1} S_kd.
  dropped = setdiff(seq_along(sys[["rest"]]), keep)
  null_q = NULL
  scaled_null_q = NULL
  information_scale = sqrt(pmax(sys[["diag"]], .Machine[["double.xmin"]]))
  if (length(dropped)) {
    v2 = matrix(0, length(sys[["rest"]]), length(dropped))
    v2[cbind(dropped, seq_along(dropped))] = 1
    if (length(keep)) {
      v2[keep, ] = -backsolve(factor, backsolve(factor, S[keep, dropped, drop = FALSE] / scale, transpose = TRUE)) / scale
    }
    V = matrix(0, sys[["n"]], length(dropped))
    V[sys[["rest"]], ] = v2
    V[sys[["big"]], ] = -scaled %*% v2
    if (retain_null) null_q = qr.Q(qr(V))
    # Remove SQL roundoff in information-scaled coordinates. An unscaled
    # projection can move O(eps) from a large group into an O(1e-20) equation,
    # corrupting that group's learned effect by many orders of magnitude.
    scaled_null_q = qr.Q(qr(information_scale * V))
  }
  solve = function(rhs) {
    r1 = rhs[sys[["big"]], , drop = FALSE]
    r2 = rhs[sys[["rest"]], , drop = FALSE] - crossprod(scaled, r1)
    x2 = matrix(0, length(sys[["rest"]]), ncol(rhs))
    if (length(keep)) {
      x2[keep, ] = backsolve(factor, backsolve(factor, r2[keep, , drop = FALSE] / scale, transpose = TRUE)) / scale
    }
    out = matrix(0, sys[["n"]], ncol(rhs))
    out[sys[["rest"]], ] = x2
    out[sys[["big"]], ] = r1 / sys[["Dg"]] - scaled %*% x2
    out
  }
  project = function(rhs) if (is.null(scaled_null_q)) rhs else
    rhs - information_scale * (scaled_null_q %*% crossprod(scaled_null_q, rhs / information_scale))
  list(solve = solve, project = project, null_basis = null_q)
}

# One statement returns the global row (__set = 0), per-FE marginals
# (__set = i) and pairwise FE information (__set = K + pair index). Pairs need
# only the weight sum, so each is a narrow GROUP BY over one shared source;
# GROUPING SETS would evaluate every aggregate for every set.
fe_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL, compact = FALSE) {
  if (identical(detect_backend(conn)[["name"]], "spark")) {
    sql = dbglm_databricks_moments_sql(source, fe, value_terms, global_terms, limit)
    if (!is.null(sql)) return(if (compact) pack_fe_moments_sql(conn, sql, length(fe), value_terms, global_terms, globals_packed = TRUE) else sql)
  }
  K = length(fe)
  pairs = utils::combn(K, 2L, simplify = FALSE)
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  null = function(terms) if (length(terms)) paste0("CAST(NULL AS DOUBLE) AS ", alias(terms)) else character()
  id = function(f) if (is.null(f)) "CAST(NULL AS BIGINT)" else paste0("CAST(", f, " AS BIGINT)")
  branch = function(set, i, j, columns, group) paste0("SELECT ", set, " AS __set,", id(i), " AS __i,", id(j), " AS __j,",
    paste(columns, collapse = ","), " FROM __dbglm_src", if (length(group)) paste(" GROUP BY", paste(group, collapse = ",")))
  branches = c(
    if (length(global_terms)) branch(0L, NULL, NULL, c(global_terms, null(value_terms)), character()),
    vapply(seq_len(K), function(i) branch(i, fe[i], NULL, c(null(global_terms), value_terms), fe[i]), character(1)),
    vapply(seq_along(pairs), function(k) {
      p = pairs[[k]]
      branch(K + k, fe[p[1]], fe[p[2]], c(null(global_terms), value_terms[1L], null(value_terms[-1L])), fe[p])
    }, character(1)))
  materialized = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "" else "MATERIALIZED "
  sql = paste0("WITH __dbglm_src AS ", materialized, "(SELECT * FROM ", source, ") ", paste(branches, collapse = " UNION ALL "))
  if (!is.null(limit)) sql = paste0("SELECT * FROM (", sql, ") __dbglm_m ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  if (compact) pack_fe_moments_sql(conn, sql, K, value_terms, global_terms) else sql
}

# Pair rows carry only indices and weight. Only marginal rows carry RHS
# vectors; only the global row carries global moments. Packing before DBI
# transfer avoids recreating a wide mostly-NULL rectangle in client memory.
pack_fe_moments_sql = function(conn, sql, K, value_terms, global_terms = character(), globals_packed = FALSE) {
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  values = alias(value_terms[-1L])
  native = identical(detect_backend(conn)[["name"]], "spark")
  payload = if (length(values)) paste0("CASE WHEN __set BETWEEN 1 AND ", K,
    " THEN ", if (native) "to_json(array(" else "list_value(", paste(values, collapse = ","),
    if (native) ")) END" else ") END") else if (native) "CAST(NULL AS VARCHAR)" else "CAST(NULL AS DOUBLE[])"
  fields = c("__set", "__i", "__j", "__h", paste(payload, "AS __values"))
  if (length(global_terms)) {
    globals = alias(global_terms)
    encoded = if (globals_packed) "__globals" else paste0("CASE WHEN __set=0 THEN list_value(",
      paste(globals, collapse = ","), ") END")
    fields = c(fields, paste(encoded, "AS __globals"))
  }
  result = paste0("SELECT ", paste(fields, collapse = ","), " FROM (", sql, ") compact_fe_moments")
  if (length(global_terms)) attr(result, "fe_global_names") = alias(global_terms)
  result
}

read_fe_moments = function(conn, sql, max_rows = Inf) {
  result = backend_get_query(conn, sql)
  if (nrow(result) > max_rows) return(NULL)
  if ("__values" %in% names(result)) {
    payload = result[["__values"]]
    rows = if (is.list(payload)) which(lengths(payload) > 0L) else which(!is.na(payload))
    values = lapply(payload[rows], function(x) if (is.list(payload)) as.numeric(x) else as.numeric(jsonlite::fromJSON(x)))
    if (length(unique(lengths(values))) > 1L) stop("Inconsistent fixed-effect marginal widths.", call. = FALSE)
    result[["__values"]] = NULL
    attr(result, "fe_values") = list(rows = rows, values = if (length(rows)) do.call(rbind, values) else matrix(numeric(), 0L, 0L))
  }
  if (!"__globals" %in% names(result)) return(result)
  global = which(as.integer(result[["__set"]]) == 0L)
  payload = result[["__globals"]]
  if (length(global) != 1L || (is.list(payload) && !length(payload[[global]])) ||
      (!is.list(payload) && is.na(payload[global]))) {
    stop("GLM moment query did not return its global summary.", call. = FALSE)
  }
  values = if (is.list(payload)) stats::setNames(as.list(as.numeric(payload[[global]])), attr(sql, "fe_global_names")) else
    jsonlite::fromJSON(payload[global])
  # Do not recreate mostly-NULL global columns in R: their number grows as
  # p^2, and padding them across every FE-pair row can exhaust client memory.
  globals = as.data.frame(lapply(values, as.numeric), check.names = FALSE)
  result[["__globals"]] = NULL
  attr(result, "dbglm_globals") = globals
  result
}

# Value columns are __h, __b1.., __ex1..; see fe_moments_sql().
unpack_fe_moments = function(res, K, n_cols, n_extra) {
  sets = as.integer(res[["__set"]])
  b_cols = if (n_cols) paste0("__b", seq_len(n_cols)) else character()
  ex_cols = if (n_extra) paste0("__ex", seq_len(n_extra)) else character()
  compact = attr(res, "fe_values", exact = TRUE)
  info = vector("list", K)
  for (i in seq_len(K)) {
    rows = which(sets == i)
    if (!length(rows)) return(NULL)
    rows = rows[order(res[["__i"]][rows])]
    h = as.numeric(res[["__h"]][rows])
    if (is.null(compact)) {
      b = as.matrix(res[rows, b_cols, drop = FALSE])
      ex = as.matrix(res[rows, ex_cols, drop = FALSE])
    } else {
      positions = match(rows, compact[["rows"]])
      if ((n_cols + n_extra) && (anyNA(positions) || ncol(compact[["values"]]) != n_cols + n_extra)) return(NULL)
      b = if (n_cols) compact[["values"]][positions, seq_len(n_cols), drop = FALSE] else matrix(numeric(), length(rows), 0L)
      ex = if (n_extra) compact[["values"]][positions, n_cols + seq_len(n_extra), drop = FALSE] else matrix(numeric(), length(rows), 0L)
    }
    if (anyNA(res[["__i"]][rows]) || any(!is.finite(h)) || any(h <= 0) || any(!is.finite(b)) || any(!is.finite(ex))) return(NULL)
    info[[i]] = list(id = res[["__i"]][rows], h = h, b = unname(b), ex = unname(ex))
  }
  combos = utils::combn(K, 2L, simplify = FALSE)
  pairs = vector("list", length(combos))
  for (k in seq_along(combos)) {
    p = combos[[k]]
    rows = which(sets == K + k)
    if (!length(rows)) return(NULL)
    h = as.numeric(res[["__h"]][rows])
    ir = match(res[["__j"]][rows], info[[p[2]]][["id"]])
    jr = match(res[["__i"]][rows], info[[p[1]]][["id"]])
    if (any(!is.finite(h)) || anyNA(ir) || anyNA(jr)) return(NULL)
    pairs[[k]] = list(i = p[2], j = p[1], ir = ir, jr = jr, h = h)
  }
  sys = fe_system(info, pairs)
  if (is.null(sys)) return(NULL)
  B = do.call(rbind, lapply(info, `[[`, "b"))
  B_extra = do.call(rbind, lapply(info, `[[`, "ex"))
  fe_score_max = if (n_extra >= 1L) max(abs(B_extra[, 1L]) / sqrt(pmax(sys[["diag"]], .Machine[["double.xmin"]]))) else 0
  list(info = info, sys = sys, B = B, B_extra = B_extra, starts = sys[["starts"]], n = sys[["n"]],
    rows = sum(sets != 0L), fe_score_max = fe_score_max)
}

fe_value_terms = function(conn, columns, extra_sums = character()) {
  c("FSUM(w) AS __h",
    if (length(columns)) sprintf("FSUM(w*(%s)) AS __b%d", columns, seq_along(columns)),
    if (length(extra_sums)) sprintf("%s AS __ex%d", vapply(extra_sums, function(e) sql_stable_sum(conn, e), character(1)), seq_along(extra_sums)))
}

query_fe_moments = function(conn, table, fe, columns, extra_sums = character(),
                                   max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6),
                                   max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6),
                                   max_cells = NULL, fe_sizes = NULL) {
  validate_moment_limits(max_edges, max_entries)
  width = length(columns) + length(extra_sums)
  if (!is.null(max_cells)) {
    validate_moment_limits(max_cells, max_cells)
    # Reserve the actual marginal payload when sizes are already known, or
    # conservatively reserve the allowed FE-by-design bound otherwise.
    payload = if (is.null(fe_sizes)) max_entries * width / max(1L, length(columns)) else sum(fe_sizes)*width
    max_edges = min(max_edges, floor((max_cells - payload) / 5L) - 1L)
    if (max_edges < 1L) return(NULL)
  }
  sql = fe_moments_sql(conn, table, fe, fe_value_terms(conn, columns, extra_sums), limit = max_edges + 1, compact = TRUE)
  res = read_fe_moments(conn, sql, max_rows = max_edges)
  if (is.null(res) || !nrow(res)) return(NULL)
  unpacked = unpack_fe_moments(res, length(fe), length(columns), length(extra_sums))
  if (is.null(unpacked)) return(NULL)
  if (unpacked[["n"]] * max(1L, length(columns)) > max_entries) return(NULL)
  unpacked[["transfer_cells"]] = nrow(res)*5 + unpacked[["n"]]*width
  unpacked
}

validate_moment_limits = function(max_edges, max_entries) {
  valid = function(x) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 1 && x <= 2^53 - 2 && x == floor(x)
  if (!valid(max_edges) || !valid(max_entries)) {
    stop("GLM moment transfer limits must be positive whole numbers below 2^53-1.", call. = FALSE)
  }
  invisible(NULL)
}

project_fe_moments = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop, cache,
                                  materialize = TRUE, extra_sums = character(), need_table = TRUE,
                                  precomputed = NULL, key = "__dbglm_key") {
  all_intercepts = all(vapply(ops, function(op) identical(op[["v"]], "1.0"), logical(1)))
  n_cols = length(columns)
  n_extra = length(extra_sums)
  if (!is.null(precomputed)) {
    info = lapply(precomputed[["info"]], function(x) list(id = x[["id"]], h = x[["h"]]))
    sys = precomputed[["sys"]]
    B = if (n_cols) precomputed[["B"]][, seq_len(n_cols), drop = FALSE] else matrix(0, precomputed[["n"]], 0L)
    B_extra = if (n_extra) precomputed[["B_extra"]][, seq_len(n_extra), drop = FALSE] else matrix(0, precomputed[["n"]], 0L)
    n = precomputed[["n"]]
    moment_rows = precomputed[["rows"]]
  } else if (all_intercepts && length(ops) >= 2L) {
    fe = vapply(ops, `[[`, character(1), "f")
    unpacked = query_fe_moments(conn, table, fe, columns, extra_sums,
      max_cells = projection_limits()[["transfer_cells"]], fe_sizes = cache[["fe_sizes"]])
    if (is.null(unpacked)) return(NULL)
    info = lapply(unpacked[["info"]], function(x) list(id = x[["id"]], h = x[["h"]]))
    sys = unpacked[["sys"]]
    B = unpacked[["B"]]
    B_extra = unpacked[["B_extra"]]
    n = unpacked[["n"]]
    moment_rows = unpacked[["rows"]]
  } else {
    limits = cache[["projection_limits"]]
    if (is.null(limits)) limits = projection_limits()
    budget = list2env(list(rows = 0, cells = 0), parent = emptyenv())
    info = lapply(ops, function(op) projection_fetch(conn, paste0("SELECT ", op[["f"]], " AS id,FSUM(w*(", op[["v"]], ")*(", op[["v"]], ")) AS h",
      if (n_cols) paste0(",", paste(sprintf("FSUM(w*(%s)*%s) AS b%d", op[["v"]], columns, seq_len(n_cols)), collapse = ",")),
      if (n_extra) paste0(",", paste(sprintf("FSUM((%s)*%s) AS ex%d", op[["v"]], extra_sums, seq_len(n_extra)), collapse = ",")),
      " FROM ", table, " GROUP BY ", op[["f"]], " ORDER BY ", op[["f"]]),
      budget, 2L+n_cols+n_extra, limits))
    pairs = list()
    for (i in seq_along(ops)) if (i > 1) for (j in seq_len(i - 1L)) {
      oi = ops[[i]]
      oj = ops[[j]]
      edge = projection_fetch(conn, paste0("SELECT ", oi[["f"]], " AS i,", oj[["f"]], " AS j,FSUM(w*(", oi[["v"]], ")*(", oj[["v"]], ")) AS h FROM ",
        table, " GROUP BY ", oi[["f"]], ",", oj[["f"]]), budget, 3L, limits)
      pairs[[length(pairs) + 1L]] = list(i = i, j = j, h = as.numeric(edge[["h"]]),
        ir = match(edge[["i"]], info[[i]][["id"]]), jr = match(edge[["j"]], info[[j]][["id"]]))
    }
    sys = fe_system(lapply(info, function(x) list(id = x[["id"]], h = as.numeric(x[["h"]]))), pairs)
    if (is.null(sys)) return(NULL)
    n = sys[["n"]]
    moment_rows = NULL
    B = if (n_cols) do.call(rbind, lapply(info, function(x) as.matrix(x[, paste0("b", seq_len(n_cols)), drop = FALSE]))) else matrix(0, n, 0L)
    B_extra = if (n_extra) do.call(rbind, lapply(info, function(x) as.matrix(x[, paste0("ex", seq_len(n_extra)), drop = FALSE]))) else matrix(0, n, 0L)
  }
  ids = unlist(lapply(seq_along(ops), function(i) paste(i, info[[i]][["id"]], sep = "/")), use.names = FALSE)
  if (!is.null(cache[["moment_keys"]]) && !identical(ids, cache[["moment_keys"]])) stop("The fixed-effect sample changed within fitting.", call. = FALSE)
  cache[["moment_keys"]] = ids
  if (all_intercepts && length(ops) >= 2L && !is.null(moment_rows)) {
    # The owned base snapshot is fixed before IRLS; working expressions only
    # change values, never FE keys or row membership. Thus group/pair support
    # has a fixed, already bounded cardinality for subsequent evaluations.
    cache[["moment_rows"]] = moment_rows
    cache[["moment_fe"]] = vapply(ops, `[[`, character(1), "f")
  }
  diagonal = sqrt(pmax(sys[["diag"]], 0))
  solver = fe_solver(sys, cache)
  solve_rhs = solver[["solve"]]
  if (n_cols) B = solver[["project"]](B)
  coefficients = matrix(0, n, n_cols)
  if (n_cols) coefficients = solve_rhs(B)
  residual = 0
  residual_rhs = matrix(0, n, n_cols)
  if (n_cols) {
    max_cycles = min(control[["fe_max_iter"]], 8L)
    for (cycle in seq_len(max_cycles)) {
      residual_rhs = B - fe_multiply(sys, coefficients)
      normalized = residual_rhs / pmax(diagonal, .Machine[["double.xmin"]])
      residual = max(abs(sweep(normalized, 2, scales, "/")))
      if (is.finite(residual) && residual <= target) break
      if (cycle == max_cycles) {
        projection_abort("nonconverged", "GLM fixed-effect projection did not converge to its group-equation tolerance.", list(residual = residual, target = target))
      }
      coefficients = coefficients + solve_rhs(residual_rhs)
    }
  } else {
    cycle = 1L
  }
  cross_reduction = if (n_cols) {
    cr = crossprod(B, coefficients) + crossprod(residual_rhs, coefficients)
    0.5 * (cr + t(cr))
  } else matrix(0, 0L, 0L)
  extra_reduction = if (n_extra && n_cols) crossprod(coefficients, B_extra) else NULL
  if (!need_table) {
    return(list(table = NULL, raw_projected = NULL, cycles = cycle, residual = residual, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = character()))
  }
  if (!n_cols) {
    projected = create(paste("SELECT", paste(c(key, identifiers, "w"), collapse = ","), "FROM", table), materialize = materialize)
    return(list(table = projected, raw_projected = NULL, cycles = cycle, residual = 0, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = character()))
  }
  frame = dbglm_coefficient_frame(info, coefficients)
  is_dbx = isTRUE(detect_backend(conn)[["is_databricks"]])
  inline_intermediate = !materialize && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))
  dictionary = dbglm_coefficient_lookup(conn, frame, create, inline = inline_intermediate)
  lookup = dictionary[["source"]]
  lookup_sub = dictionary[["tables"]]
  ctes = dbglm_coefficient_ctes(lookup, length(ops), "proj")
  joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN proj_g", i, " ON s.", ops[[i]][["f"]], "=proj_g", i, ".__glm_group"), character(1))
  broadcast = if (is_dbx) paste0("/*+ BROADCAST(", paste0("proj_g", seq_along(ops), collapse = ","), ") */ ") else ""
  if (inline_intermediate) {
    tilde_exprs = vapply(seq_len(n_cols), function(j) {
      terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
        sprintf("(%s)*COALESCE(proj_g%d.c%d,0)", v, i, j)
      }, character(1))
      paste0("CAST((s.", columns[j], "-(", paste(terms, collapse = " + "), ")) AS DOUBLE) AS ", columns[j], "_tilde")
    }, character(1))
    raw_sql = paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast, "s.*, s.w AS __w, ",
      paste(tilde_exprs, collapse = ","), " FROM ", table, " s", paste(joins, collapse = ""))
    raw_projected = paste0("(", backend_sql(conn, raw_sql), ")")
    return(list(table = NULL, raw_projected = raw_projected, cycles = cycle, residual = residual, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = lookup_sub))
  }
  expressions = vapply(seq_len(n_cols), function(j) {
    terms = vapply(seq_along(ops), function(i) {
      v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
      sprintf("(%s)*COALESCE(proj_g%d.c%d,0)", v, i, j)
    }, character(1))
    paste0("CAST((s.", columns[j], "-(", paste(terms, collapse = " + "), ")) AS DOUBLE) AS ", columns[j])
  }, character(1))
  projected = create(paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast,
    paste(paste0("s.", c(key, identifiers, "w")), collapse = ","), ",",
    paste(expressions, collapse = ","), " FROM ", table, " s", paste(joins, collapse = "")), materialize = TRUE)
  for (st in lookup_sub) drop(st)
  list(table = projected, raw_projected = NULL, cycles = cycle, residual = residual, method = "moments",
    cross_reduction = cross_reduction, extra_reduction = extra_reduction,
    coefficients = coefficients, info = info, ops = ops, sub_tables = character())
}
