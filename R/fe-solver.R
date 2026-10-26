# Shared projection policy. Numerical tolerances remain estimator controls;
# these limits bound representations and allocations, not statistical accuracy.
projection_limits = function() {
  limits = list(
    summary_rows = getOption("dbreg.duckdb_ap_max_edges", 1e6),
    design_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6),
    transfer_cells = getOption("dbreg.ap_moment_max_cells", 8e6),
    dense_levels = getOption("dbreg.glm_fe_dense_max", 4096),
    cross_entries = getOption("dbreg.glm_fe_cross_max", 5e7),
    memory_bytes = getOption("dbreg.projection_max_bytes", 512 * 1024^2)
  )
  for (name in names(limits)) {
    value = limits[[name]]
    lower = if (name %in% c("dense_levels", "cross_entries")) 0 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < lower || value > 2^53 - 2 || value != floor(value)) {
      stop("Projection limits must be positive whole numbers (dense limits may be zero): ", name, ".", call. = FALSE)
    }
  }
  limits
}

projection_memory = function(sizes, columns = 1L, n_edges = NULL) {
  rest = sum(sizes) - max(sizes)
  if (is.null(n_edges)) {
    possible = (sum(sizes)^2 - sum(sizes^2)) / 2
    n_edges = min(possible, getOption("dbreg.duckdb_ap_max_edges", 1e6))
  }
  # Sparse cross blocks and scaled copies, plus a conservative dense Schur
  # factor/basis allowance and RHS buffers. This plans allocations, not RSS.
  64 * n_edges + 8 * (6 * rest^2 + 4 * sum(sizes) * max(columns, 1L))
}

projection_abort = function(reason, message, details = list()) {
  kind = if (reason %in% c("singular", "nonconverged", "certificate_failed")) "numerical" else "capacity"
  stop(structure(list(message = message, call = NULL, reason = reason, details = details),
    class = c(paste0("dbreg_projection_", kind), "dbreg_projection_error", "error", "condition")
  ))
}

projection_attempt = function(method, expr) {
  tryCatch(
    {
      result = force(expr)
      details = if (is.null(result)) list() else result[intersect(names(result), c("transfer_cells", "estimated_bytes", "iterations", "max_abs", "observation_tables"))]
      list(
        method = method, status = if (is.null(result)) "declined" else "completed",
        reason = if (is.null(result)) "not_applicable" else NULL, details = details, result = result
      )
    },
    dbreg_projection_error = function(e) list(
      method = method, status = "declined",
      reason = e[["reason"]], details = e[["details"]], result = NULL
    )
  )
}

record_projection_attempt = function(context, attempt) {
  summary = attempt[setdiff(names(attempt), "result")]
  previous = context[["projection_diagnostics"]]
  # Record method transitions/rejections, not one duplicate per IRLS update.
  if (!length(previous) || !identical(utils::tail(previous, 1L)[[1L]], summary)) {
    context[["projection_diagnostics"]] = c(previous, list(summary))
  }
  invisible(attempt[["result"]])
}

# Account across the separate marginal/pair statements needed by varying slopes.
# Reserve an overflow probe before each DBI fetch and propagate typed capacity
# rejection through the same solver policy as packed intercept-only summaries.
projection_fetch = function(conn, sql, budget, width, limits) {
  rows = min(
    limits[["summary_rows"]] - budget[["rows"]],
    floor((limits[["transfer_cells"]] - budget[["cells"]]) / width) - 1L
  )
  if (rows < 1L) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  result = backend_get_query(conn, paste(sql, "LIMIT", sprintf("%.0f", rows + 1L)))
  if (nrow(result) > rows) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  budget[["rows"]] = budget[["rows"]] + nrow(result)
  budget[["cells"]] = budget[["cells"]] + width * (nrow(result) + 1L)
  result
}

# Bounded direct nuisance solves for moderate FE dimensions. The database sends
# only FE cross-products and FE-by-design moments; no observation vectors enter R.
# The largest FE block is diagonal and is absorbed analytically, so only the
# remaining levels form the smaller Schur system; both cross blocks stay sparse.
fe_dense_feasible = function(op_sizes, columns = 1L, limits = projection_limits(), n_edges = NULL) {
  op_sizes = as.numeric(op_sizes)
  if (!length(op_sizes) || any(!is.finite(op_sizes))) return(FALSE)
  rest = sum(op_sizes) - max(op_sizes)
  if (is.null(n_edges)) n_edges = min((sum(op_sizes)^2 - sum(op_sizes^2)) / 2, limits[["summary_rows"]])
  rest <= limits[["dense_levels"]] &&
    (is.null(n_edges) || n_edges <= limits[["cross_entries"]]) &&
    projection_memory(op_sizes, columns, n_edges) <= limits[["memory_bytes"]]
}

# info: per-op list(id, h); pairs: list of list(i, j, ir, jr, h) with local
# level indices into info[[i]]$id and info[[j]]$id.
fe_system = function(info, pairs) {
  sizes = vapply(info, function(x) length(x[["id"]]), integer(1))
  n_edges = sum(vapply(pairs, function(pair) length(pair[["h"]]), integer(1)))
  if (!fe_dense_feasible(sizes, n_edges = n_edges)) return(NULL)
  starts = cumsum(c(0L, sizes))
  n = starts[length(info) + 1L]
  a = which.max(sizes)
  big = seq.int(starts[a] + 1L, starts[a + 1L])
  rest = setdiff(seq_len(n), big)
  pos = integer(n)
  pos[rest] = seq_along(rest)
  diagonal = as.numeric(unlist(lapply(info, `[[`, "h"), use.names = FALSE))
  cross_i = cross_j = cross_h = list()
  rest_i = rest_j = rest_h = list()
  for (k in seq_along(pairs)) {
    pair = pairs[[k]]
    gi = starts[pair[["i"]]] + pair[["ir"]]
    gj = starts[pair[["j"]]] + pair[["jr"]]
    if (pair[["i"]] == a || pair[["j"]] == a) {
      cross_i[[k]] = if (pair[["i"]] == a) gi - starts[a] else gj - starts[a]
      cross_j[[k]] = if (pair[["i"]] == a) pos[gj] else pos[gi]
      cross_h[[k]] = pair[["h"]]
    } else {
      rest_i[[k]] = c(pos[gi], pos[gj])
      rest_j[[k]] = c(pos[gj], pos[gi])
      rest_h[[k]] = rep(pair[["h"]], 2L)
    }
  }
  Cm = Matrix::sparseMatrix(
    i = as.integer(unlist(cross_i)), j = as.integer(unlist(cross_j)),
    x = as.numeric(unlist(cross_h)), dims = c(length(big), length(rest))
  )
  E = Matrix::sparseMatrix(
    i = c(seq_along(rest), as.integer(unlist(rest_i))),
    j = c(seq_along(rest), as.integer(unlist(rest_j))),
    x = c(diagonal[rest], as.numeric(unlist(rest_h))), dims = c(length(rest), length(rest))
  )

  list(n = n, starts = starts, big = big, rest = rest, Dg = diagonal[big], Cm = Cm, E = E, diag = diagonal, sizes = sizes)
}

fe_multiply = function(sys, X) {
  x1 = X[sys[["big"]], , drop = FALSE]
  x2 = X[sys[["rest"]], , drop = FALSE]
  out = matrix(0, sys[["n"]], ncol(X))
  out[sys[["big"]], ] = sys[["Dg"]] * x1 + as.matrix(sys[["Cm"]] %*% x2)
  out[sys[["rest"]], ] = as.matrix(Matrix::crossprod(sys[["Cm"]], x1) + sys[["E"]] %*% x2)
  out
}

# Two-way intercept information has one signed constant null direction per
# connected component. Pin one smaller-family level in each component.
fe_two_components = function(sys) {
  C = sys[["Cm"]]
  if (length(sys[["sizes"]]) != 2L || any(C < 0)) return(NULL)
  if (any(abs(as.numeric(Matrix::rowSums(C)) - sys[["Dg"]]) > 1e-12 * sys[["Dg"]]) ||
    any(abs(as.numeric(Matrix::colSums(C)) - sys[["diag"]][sys[["rest"]]]) >
      1e-12 * sys[["diag"]][sys[["rest"]]])) return(NULL)
  edges = Matrix::summary(C)
  edge_i = edges[["i"]]
  edge_j = edges[["j"]]
  edge_h = edges[["x"]]
  parent = seq_len(sys[["n"]])
  rank = integer(sys[["n"]])
  for (k in seq_along(edge_h)) {
    if (edge_h[[k]] <= 0) next
    i = sys[["big"]][edge_i[[k]]]
    j = sys[["rest"]][edge_j[[k]]]
    while (parent[[i]] != i) {
      parent[[i]] = parent[[parent[[i]]]]
      i = parent[[i]]
    }
    while (parent[[j]] != j) {
      parent[[j]] = parent[[parent[[j]]]]
      j = parent[[j]]
    }
    if (i == j) next
    if (rank[[i]] < rank[[j]]) parent[[i]] = j else {
      parent[[j]] = i
      if (rank[[i]] == rank[[j]]) rank[[i]] = rank[[i]] + 1L
    }
  }
  roots = integer(length(parent))
  for (k in seq_along(parent)) {
    i = k
    while (parent[[i]] != i) i = parent[[i]]
    roots[[k]] = i
  }
  match(roots, unique(roots))
}

fe_sparse_solver = function(sys, S, scaled, cache, retain_null) {
  if (length(sys[["sizes"]]) != 2L) return(NULL)
  nr = length(sys[["rest"]])
  sparse = nr > 512L && Matrix::nnzero(S) / nr^2 < 0.15
  components = cache[["sparse_components"]]
  pattern = list(sys[["Cm"]]@i, sys[["Cm"]]@p, sys[["Cm"]]@x > 0)
  if (is.null(components) || !identical(pattern, cache[["sparse_pattern"]])) {
    components = fe_two_components(sys)
    if (is.null(components)) return(NULL)
    cache[["sparse_components"]] = components
    cache[["sparse_pattern"]] = pattern
    cache[["sparse_factor"]] = NULL
  }
  rest_components = components[sys[["rest"]]]
  keep = which(duplicated(rest_components))
  scale = sqrt(pmax(as.numeric(Matrix::diag(S))[keep], 0))
  if (any(!is.finite(scale)) || any(scale <= 0)) projection_abort("singular", "Fixed-effect information is singular.")
  factor = NULL
  if (length(keep)) {
    factor = tryCatch(
      {
        if (sparse) {
          D = Matrix::Diagonal(x = 1 / scale)
          matrix = Matrix::forceSymmetric(D %*% S[keep, keep, drop = FALSE] %*% D)
          previous = cache[["sparse_factor"]]
          if (is.null(previous)) Matrix::Cholesky(matrix, perm = TRUE, LDL = TRUE) else Matrix::update(previous, parent = matrix)
        } else chol(as.matrix(S[keep, keep, drop = FALSE]) / outer(scale, scale))
      },
      error = function(error) NULL
    )
    if (is.null(factor)) projection_abort("singular", "Fixed-effect factorization failed.")
    if (sparse) cache[["sparse_factor"]] = factor
  }
  solve = function(rhs) {
    r1 = rhs[sys[["big"]], , drop = FALSE]
    r2 = rhs[sys[["rest"]], , drop = FALSE] - as.matrix(Matrix::crossprod(scaled, r1))
    x2 = matrix(0, nr, ncol(rhs))
    if (length(keep)) {
      rhs2 = r2[keep, , drop = FALSE] / scale
      x2[keep, ] = if (sparse) as.matrix(Matrix::solve(factor, rhs2)) / scale else
        backsolve(factor, backsolve(factor, rhs2, transpose = TRUE)) / scale
    }
    out = matrix(0, sys[["n"]], ncol(rhs))
    out[sys[["rest"]], ] = x2
    out[sys[["big"]], ] = r1 / sys[["Dg"]] - as.matrix(scaled %*% x2)
    out
  }
  sign = rep(-1, sys[["n"]])
  sign[sys[["big"]]] = 1
  component_weight = as.numeric(rowsum(sys[["diag"]], components, reorder = FALSE))
  project = function(rhs) {
    totals = rowsum(sign * rhs, components, reorder = FALSE)
    rhs - (sign * (sys[["diag"]] / component_weight[components])) * totals[components, , drop = FALSE]
  }
  null_q = NULL
  if (retain_null) {
    count = max(components)
    if (8 * sys[["n"]] * count > projection_limits()[["memory_bytes"]]) {
      projection_abort("system_limit", "Fixed-effect null basis exceeds the prediction memory budget.")
    }
    null_q = matrix(0, sys[["n"]], count)
    sizes = tabulate(components)
    null_q[cbind(seq_len(sys[["n"]]), components)] = sign / sqrt(sizes[components])
  }
  list(solve = solve, project = project, null_basis = null_q, method = if (sparse) "sparse_cholesky" else "dense_cholesky")
}

# Factor the Schur complement on the non-absorbed levels; the identified basis
# is chosen once per fit so FE normalisations stay fixed across iterations.
fe_solver = function(sys, cache, retain_null = FALSE) {
  if (any(!is.finite(sys[["Dg"]])) || any(sys[["Dg"]] <= 0)) projection_abort("singular", "Fixed-effect information is singular.")
  scaled = sys[["Cm"]]
  if (inherits(scaled, "CsparseMatrix")) {
    scaled@x = scaled@x / sys[["Dg"]][scaled@i + 1L]
  } else scaled = scaled / sys[["Dg"]]
  S = Matrix::forceSymmetric(sys[["E"]] - Matrix::crossprod(sys[["Cm"]], scaled))
  sparse = fe_sparse_solver(sys, S, scaled, cache, retain_null)
  if (!is.null(sparse)) return(sparse)
  S = as.matrix(S)
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
    if (32 * sys[["n"]] * length(dropped) > projection_limits()[["memory_bytes"]]) {
      projection_abort("system_limit", "Fixed-effect null projection exceeds its memory budget.")
    }
    v2 = matrix(0, length(sys[["rest"]]), length(dropped))
    v2[cbind(dropped, seq_along(dropped))] = 1
    if (length(keep)) {
      v2[keep, ] = -backsolve(factor, backsolve(factor, S[keep, dropped, drop = FALSE] / scale, transpose = TRUE)) / scale
    }
    V = matrix(0, sys[["n"]], length(dropped))
    V[sys[["rest"]], ] = v2
    V[sys[["big"]], ] = -as.matrix(scaled %*% v2)
    if (retain_null) null_q = qr.Q(qr(V))
    # Remove SQL roundoff in information-scaled coordinates. An unscaled
    # projection can move O(eps) from a large group into an O(1e-20) equation,
    # corrupting that group's learned effect by many orders of magnitude.
    scaled_null_q = qr.Q(qr(information_scale * V))
  }
  solve = function(rhs) {
    r1 = rhs[sys[["big"]], , drop = FALSE]
    r2 = rhs[sys[["rest"]], , drop = FALSE] - as.matrix(Matrix::crossprod(scaled, r1))
    x2 = matrix(0, length(sys[["rest"]]), ncol(rhs))
    if (length(keep)) {
      x2[keep, ] = backsolve(factor, backsolve(factor, r2[keep, , drop = FALSE] / scale, transpose = TRUE)) / scale
    }
    out = matrix(0, sys[["n"]], ncol(rhs))
    out[sys[["rest"]], ] = x2
    out[sys[["big"]], ] = r1 / sys[["Dg"]] - as.matrix(scaled %*% x2)
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
fe_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL, compact = FALSE, text_keys = FALSE) {
  if (identical(detect_backend(conn)[["name"]], "spark")) {
    sql = dbglm_databricks_moments_sql(source, fe, value_terms, global_terms, limit, text_keys = text_keys)
    if (!is.null(sql)) return(if (compact) pack_fe_moments_sql(conn, sql, length(fe), value_terms, global_terms, globals_packed = TRUE) else sql)
  }
  K = length(fe)
  pairs = utils::combn(K, 2L, simplify = FALSE)
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  null = function(terms) if (length(terms)) paste0("CAST(NULL AS DOUBLE) AS ", alias(terms)) else character()
  id_type = if (text_keys) "VARCHAR" else "BIGINT"
  id = function(f) if (is.null(f)) paste0("CAST(NULL AS ", id_type, ")") else paste0("CAST(", f, " AS ", id_type, ")")
  branch = function(set, i, j, columns, group) paste0(
    "SELECT ", set, " AS __set,", id(i), " AS __i,", id(j), " AS __j,",
    paste(columns, collapse = ","), " FROM __dbglm_src", if (length(group)) paste(" GROUP BY", paste(group, collapse = ","))
  )
  branches = c(
    if (length(global_terms)) branch(0L, NULL, NULL, c(global_terms, null(value_terms)), character()),
    vapply(seq_len(K), function(i) branch(i, fe[i], NULL, c(null(global_terms), value_terms), fe[i]), character(1)),
    vapply(seq_along(pairs), function(k) {
      p = pairs[[k]]
      branch(K + k, fe[p[1]], fe[p[2]], c(null(global_terms), value_terms[1L], null(value_terms[-1L])), fe[p])
    }, character(1))
  )
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
  payload = if (length(values)) paste0(
    "CASE WHEN __set BETWEEN 1 AND ", K,
    " THEN ", if (native) "to_json(array(" else "list_value(", paste(values, collapse = ","),
    if (native) ")) END" else ") END"
  ) else if (native) "CAST(NULL AS VARCHAR)" else "CAST(NULL AS DOUBLE[])"
  fields = c("__set", "__i", "__j", "__h", paste(payload, "AS __values"))
  if (length(global_terms)) {
    globals = alias(global_terms)
    encoded = if (globals_packed) "__globals" else paste0(
      "CASE WHEN __set=0 THEN list_value(",
      paste(globals, collapse = ","), ") END"
    )
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
    values = lapply(payload[rows], function(x) if (is.list(payload)) as.numeric(x) else as.numeric(read_numeric_json(x)))
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
    read_numeric_json(payload[global])
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
    if (anyNA(res[["__i"]][rows]) || anyDuplicated(res[["__i"]][rows]) || any(!is.finite(h)) || any(h <= 0) || any(!is.finite(b)) || any(!is.finite(ex))) return(NULL)
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
  list(
    info = info, sys = sys, B = B, B_extra = B_extra, starts = sys[["starts"]], n = sys[["n"]],
    rows = sum(sets != 0L), fe_score_max = fe_score_max
  )
}

fe_value_terms = function(conn, columns, extra_sums = character()) {
  c(
    "FSUM(w) AS __h",
    if (length(columns)) sprintf("FSUM(w*(%s)) AS __b%d", columns, seq_along(columns)),
    if (length(extra_sums)) sprintf("%s AS __ex%d", vapply(extra_sums, function(e) sql_stable_sum(conn, e), character(1)), seq_along(extra_sums))
  )
}

query_fe_moments = function(conn, table, fe, columns, extra_sums = character(),
                            max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6),
                            max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6),
                            max_cells = NULL, fe_sizes = NULL, global_terms = character(), text_keys = FALSE) {
  validate_moment_limits(max_edges, max_entries)
  width = length(columns) + length(extra_sums)
  if (!is.null(max_cells)) {
    validate_moment_limits(max_cells, max_cells)
    # Reserve the actual marginal payload when sizes are already known, or
    # conservatively reserve the allowed FE-by-design bound otherwise.
    payload = if (is.null(fe_sizes)) max_entries * width / max(1L, length(columns)) else sum(fe_sizes) * width
    max_edges = min(max_edges, floor((max_cells - payload - length(global_terms)) / 5L) - 1L)
    if (max_edges < 1L) return(NULL)
  }
  sql = fe_moments_sql(conn, table, fe, fe_value_terms(conn, columns, extra_sums),
    global_terms = global_terms,
    limit = max_edges + 1, compact = TRUE, text_keys = text_keys
  )
  res = read_fe_moments(conn, sql, max_rows = max_edges)
  if (is.null(res) || !nrow(res)) return(NULL)
  unpacked = unpack_fe_moments(res, length(fe), length(columns), length(extra_sums))
  if (is.null(unpacked)) return(NULL)
  if (unpacked[["n"]] * max(1L, length(columns)) > max_entries) return(NULL)
  unpacked[["globals"]] = attr(res, "dbglm_globals", exact = TRUE)
  unpacked[["transfer_cells"]] = nrow(res) * 5 + unpacked[["n"]] * width + length(global_terms)
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
      max_cells = projection_limits()[["transfer_cells"]], fe_sizes = cache[["fe_sizes"]]
    )
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
    info = lapply(ops, function(op) projection_fetch(
      conn, paste0(
        "SELECT ", op[["f"]], " AS id,FSUM(w*(", op[["v"]], ")*(", op[["v"]], ")) AS h",
        if (n_cols) paste0(",", paste(sprintf("FSUM(w*(%s)*%s) AS b%d", op[["v"]], columns, seq_len(n_cols)), collapse = ",")),
        if (n_extra) paste0(",", paste(sprintf("FSUM((%s)*%s) AS ex%d", op[["v"]], extra_sums, seq_len(n_extra)), collapse = ",")),
        " FROM ", table, " GROUP BY ", op[["f"]], " ORDER BY ", op[["f"]]
      ),
      budget, 2L + n_cols + n_extra, limits
    ))
    pairs = list()
    for (i in seq_along(ops)) if (i > 1) for (j in seq_len(i - 1L)) {
      oi = ops[[i]]
      oj = ops[[j]]
      edge = projection_fetch(conn, paste0(
        "SELECT ", oi[["f"]], " AS i,", oj[["f"]], " AS j,FSUM(w*(", oi[["v"]], ")*(", oj[["v"]], ")) AS h FROM ",
        table, " GROUP BY ", oi[["f"]], ",", oj[["f"]]
      ), budget, 3L, limits)
      pairs[[length(pairs) + 1L]] = list(
        i = i, j = j, h = as.numeric(edge[["h"]]),
        ir = match(edge[["i"]], info[[i]][["id"]]), jr = match(edge[["j"]], info[[j]][["id"]])
      )
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
    return(list(
      table = NULL, raw_projected = NULL, cycles = cycle, residual = residual, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = character()
    ))
  }
  if (!n_cols) {
    projected = create(paste("SELECT", paste(c(key, identifiers, "w"), collapse = ","), "FROM", table), materialize = materialize)
    return(list(
      table = projected, raw_projected = NULL, cycles = cycle, residual = 0, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = character()
    ))
  }
  frame = dbglm_coefficient_frame(info, coefficients)
  is_dbx = isTRUE(detect_backend(conn)[["is_databricks"]])
  inline_intermediate = !materialize && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))
  dictionary = dbglm_coefficient_lookup(conn, frame, create, inline = inline_intermediate)
  lookup = dictionary[["source"]]
  lookup_sub = dictionary[["tables"]]
  ctes = dbglm_coefficient_ctes(lookup, length(ops), "proj")
  text_keys = is.character(frame[["__glm_group"]])
  joins = vapply(seq_along(ops), function(i) {
    key = paste0("s.", ops[[i]][["f"]])
    if (text_keys) key = paste0("CAST(", key, " AS VARCHAR)")
    paste0(" LEFT JOIN proj_g", i, " ON ", key, "=proj_g", i, ".__glm_group")
  }, character(1))
  broadcast = if (is_dbx) paste0("/*+ BROADCAST(", paste0("proj_g", seq_along(ops), collapse = ","), ") */ ") else ""
  if (inline_intermediate) {
    tilde_exprs = vapply(seq_len(n_cols), function(j) {
      terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
        sprintf("(%s)*COALESCE(proj_g%d.c%d,0)", v, i, j)
      }, character(1))
      paste0("CAST((s.", columns[j], "-(", paste(terms, collapse = " + "), ")) AS DOUBLE) AS ", columns[j], "_tilde")
    }, character(1))
    raw_sql = paste0(
      "WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast, "s.*, s.w AS __w, ",
      paste(tilde_exprs, collapse = ","), " FROM ", table, " s", paste(joins, collapse = "")
    )
    raw_projected = paste0("(", backend_sql(conn, raw_sql), ")")
    return(list(
      table = NULL, raw_projected = raw_projected, cycles = cycle, residual = residual, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = lookup_sub
    ))
  }
  expressions = vapply(seq_len(n_cols), function(j) {
    terms = vapply(seq_along(ops), function(i) {
      v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
      sprintf("(%s)*COALESCE(proj_g%d.c%d,0)", v, i, j)
    }, character(1))
    paste0("CAST((s.", columns[j], "-(", paste(terms, collapse = " + "), ")) AS DOUBLE) AS ", columns[j])
  }, character(1))
  projected = create(paste0(
    "WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast,
    paste(paste0("s.", c(key, identifiers, "w")), collapse = ","), ",",
    paste(expressions, collapse = ","), " FROM ", table, " s", paste(joins, collapse = "")
  ), materialize = TRUE)
  for (st in lookup_sub) drop(st)
  list(
    table = projected, raw_projected = NULL, cycles = cycle, residual = residual, method = "moments",
    cross_reduction = cross_reduction, extra_reduction = extra_reduction,
    coefficients = coefficients, info = info, ops = ops, sub_tables = character()
  )
}

# Work in D2-scaled coordinates so component deflation remains stable for
# unequal group weights. Bounded Schur fill permits symmetric Gauss-Seidel;
# otherwise use an additive block SGS/Jacobi preconditioner without that fill.
two_way_cg_preconditioner = function(C, d1, d2, limits, available_bytes) {
  n1 = nrow(C)
  n2 = ncol(C)
  system = list(
    n = n1 + n2, big = seq_len(n1), rest = n1 + seq_len(n2),
    Cm = C, Dg = d1, diag = c(d1, d2), sizes = c(n1, n2)
  )
  components = fe_two_components(system)
  if (is.null(components)) return(NULL)
  components = components[n1 + seq_len(n2)]
  total = as.numeric(rowsum(d2, components, reorder = TRUE))
  root = sqrt(d2)
  null = root / sqrt(total[components])
  project = function(value) {
    sums = rowsum(null * value, components, reorder = TRUE)
    value - null * sums[components, , drop = FALSE]
  }
  row_scaled = C
  row_scaled@x = row_scaled@x / d1[row_scaled@i + 1L]
  diagonal = pmax(1 - as.numeric(Matrix::colSums(C * row_scaled)) / d2, .Machine[["double.eps"]])
  walk = function(value) {
    as.matrix(Matrix::crossprod(C, as.matrix(C %*% (value / root)) / d1)) / root
  }
  degree = as.numeric(Matrix::rowSums(C != 0))
  fill = n2 + sum(degree^2)
  if (fill <= limits[["cross_entries"]] && 64 * fill <= available_bytes) {
    scaled = Matrix::Diagonal(x = 1 / sqrt(d1)) %*% C %*% Matrix::Diagonal(x = 1 / root)
    S = Matrix::Diagonal(n2) - Matrix::crossprod(scaled)
    S = S + Matrix::Diagonal(x = diagonal - as.numeric(Matrix::diag(S)))
    lower = Matrix::tril(S)
    upper = Matrix::t(lower)
    precondition = function(value) as.matrix(Matrix::solve(
      upper,
      diagonal * as.matrix(Matrix::solve(lower, value))
    ))
    method = "schur_sgs"
  } else {
    precondition = function(value) value / diagonal + value + walk(value)
    method = "block_sgs"
  }
  list(
    project = project, precondition = precondition, walk = walk,
    root = root, method = method, components = max(components)
  )
}

#' Bounded two-way fixed-effect projection for DuckDB and Databricks
#'
#' The observation relation remains in the database. Weighted FE-pair counts and
#' FE-level sums define a sparse Schur complement over the smaller FE family.
#' Only these bounded sufficient statistics are transferred to R.
#'
#' Returns NULL when the bounds or numerical checks require the join/window fallback.
#' @keywords internal
#' @noRd
dbreg_duckdb_cg_projections = function(
  conn, from_statement, fe, yvar, xvars_sql, xvar_names, weights,
  cluster_var = NULL, verbose = FALSE, max_iter = 100L, tol = 1e-10,
  max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6),
  max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6), scales = NULL,
  materialize_result = TRUE, inline_result = FALSE, extra_sums = character()
) {
  valid_cap = function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  }
  if (!valid_cap(max_edges) || !valid_cap(max_entries)) {
    stop("DuckDB AP transfer limits must be positive whole numbers.", call. = FALSE)
  }
  stopifnot(length(fe) == 2L)
  backend_info = detect_backend(conn)
  vars = c(yvar, xvar_names)
  limits = projection_limits()
  width = length(vars) + length(extra_sums)
  key_bytes = function(f) {
    value = paste0("CAST(", f, " AS VARCHAR)")
    if (backend_info[["name"]] == "duckdb") paste0("OCTET_LENGTH(encode(", value, "))") else
      paste0("OCTET_LENGTH(", value, ")")
  }
  w = if (is.null(weights)) "1.0" else weights
  raw_vars = sql_cast_double(c(yvar, xvars_sql), conn)
  w_double = sql_cast_double(w, conn)
  center_terms = paste0(
    vapply(
      paste0("(", w_double, ")*(", raw_vars, ")"),
      function(value) sql_stable_sum(conn, value), character(1)
    ), "/(", sql_stable_sum(conn, w_double),
    ") AS center", seq_along(raw_vars)
  )
  preflight = backend_get_query(conn, paste(
    "SELECT COUNT(*) AS source_n,",
    paste(c(
      sprintf("COUNT(DISTINCT %s) AS n%d", fe, 1:2),
      sprintf("MAX(%s) AS bytes%d", vapply(fe, key_bytes, character(1)), 1:2), center_terms
    ), collapse = ","), from_statement
  ))
  sizes = as.numeric(preflight[c("n1", "n2")][1L, ])
  key_width = as.numeric(preflight[c("bytes1", "bytes2")][1L, ])
  if (any(!is.finite(c(sizes, key_width))) || any(sizes < 1L)) return(NULL)
  centers = vapply(seq_along(vars), function(i) as.numeric(preflight[[paste0("center", i)]]), numeric(1))
  if (any(!is.finite(centers))) return(NULL)
  vars_sql = paste0("(", raw_vars, "-(", dbglm_number(centers), "))")
  if (sum(sizes) * length(vars) > max_entries) projection_abort("system_limit", "Sparse FE system exceeds its design budget.")
  row_bytes = 4 * (sum(key_width) + 128 + 8 * (width + 3L))
  remaining = limits[["memory_bytes"]] - 8 * sum(sizes) * max(1L, width) * 12
  max_edges = min(
    max_edges, floor(limits[["transfer_cells"]] / (width + 3L)) - 1L,
    floor(remaining / row_bytes) - 1L
  )
  if (max_edges < 1L) projection_abort("system_limit", "Sparse FE system exceeds its memory or transfer budget.")
  equation_metric = !is.null(scales)
  if (equation_metric && (length(scales) != length(vars) ||
    any(!is.finite(scales)) || any(scales <= 0))) {
    stop("Projection scales must be positive and match the projected columns.", call. = FALSE)
  }
  aggregate = if (equation_metric) "FSUM" else "SUM"
  ids = setdiff(unique(c(fe, cluster_var)), vars)
  if ("__w" %in% c(ids, vars)) {
    stop("Model column conflicts with the internal AP column __w.", call. = FALSE)
  }
  seed = paste0(Sys.getpid(), "_", sprintf("%08d", sample.int(1e8, 1L)))
  tables = setNames(
    paste0("dbreg_ap_cg_", seed, "_", c("source", "edges", "a", "b", "result")),
    c("source", "edges", "a", "b", "result")
  )
  views = if (isTRUE(backend_info[["is_databricks"]]) &&
    !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
    unname(tables[c("source", "edges", "a", "b", if (!materialize_result) "result")])
  } else character()
  state = new.env(parent = emptyenv())
  state[["created"]] = character()
  drop_state = function(table) {
    if (table %in% views) try(DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", table)), silent = TRUE) else
      drop_ap_table_if_exists(conn, table, backend_info)
  }
  success = FALSE
  on.exit(
    {
      if (!success) {
        for (table in rev(state[["created"]])) drop_state(table)
      }
    },
    add = TRUE
  )
  create = function(name, sql) {
    state[["created"]] = unique(c(state[["created"]], tables[[name]]))
    if (isTRUE(backend_info[["is_databricks"]])) {
      DBI::dbExecute(conn, backend_sql(conn, paste(if (tables[[name]] %in% views) "CREATE TEMP VIEW" else "CREATE TEMP TABLE", tables[[name]], "AS", sql)))
    } else create_temp_table_as(conn, tables[[name]], sql, "duckdb")
  }
  w = if (is.null(weights)) "1.0" else weights
  # String keys preserve BIGINT identifiers beyond R's exact double range.
  keys = c("dbreg_cg_key_1", "dbreg_cg_key_2")
  key_sql = sprintf("CAST(%s AS VARCHAR)", fe)
  sum_terms = c(
    sprintf("%s((%s) * (%s)) AS s%d", aggregate, w, vars_sql, seq_along(vars)),
    if (length(extra_sums)) sprintf("%s(%s) AS ex%d", aggregate, extra_sums, seq_along(extra_sums))
  )
  edges_sql = paste0(
    "SELECT ", paste(sprintf("%s AS %s", key_sql, keys), collapse = ", "),
    ", ", aggregate, "(", w, ") AS ew, ",
    paste(sum_terms, collapse = ", "), " ", from_statement,
    " GROUP BY ", paste(fe, collapse = ", "),
    " LIMIT ", sprintf("%.0f", max_edges + 1)
  )
  edges = backend_get_query(conn, edges_sql)
  if (!nrow(edges) || nrow(edges) > max_edges) {
    if (verbose) message("[dbreg] FE summaries exceed DuckDB transfer limits; using join AP")
    projection_abort("transfer_limit", "FE pair summaries exceed sparse projection limits.", list(summary_rows = max_edges))
  }
  u1 = unique(edges[[keys[1]]])
  u2 = unique(edges[[keys[2]]])
  entries = (as.numeric(length(u1)) + as.numeric(length(u2))) * length(vars)
  if (!is.finite(entries) || entries > max_entries) {
    if (verbose) message("[dbreg] FE summaries exceed DuckDB transfer limits; using join AP")
    return(NULL)
  }
  scols = paste0("s", seq_along(vars))
  s_mat = as.matrix(edges[scols])
  ew = as.numeric(edges[["ew"]])
  if (any(!is.finite(ew)) || any(ew <= 0) || any(!is.finite(s_mat))) return(NULL)
  idx1 = match(edges[[keys[1]]], u1)
  idx2 = match(edges[[keys[2]]], u2)
  d1 = as.numeric(rowsum(ew, idx1, reorder = FALSE))
  d2 = as.numeric(rowsum(ew, idx2, reorder = FALSE))
  s1 = rowsum(s_mat, idx1, reorder = FALSE)
  s2 = rowsum(s_mat, idx2, reorder = FALSE)
  ex1 = ex2 = NULL
  if (length(extra_sums)) {
    ex_mat = as.matrix(edges[paste0("ex", seq_along(extra_sums))])
    ex1 = rowsum(ex_mat, idx1, reorder = FALSE)
    ex2 = rowsum(ex_mat, idx2, reorder = FALSE)
  }
  if (any(!is.finite(c(d1, d2))) || any(c(d1, d2) <= 0)) return(NULL)
  y1 = s1 / d1
  y2 = s2 / d2
  means = list(
    local({
      df = data.frame(k = u1, gw = d1, stringsAsFactors = FALSE)
      names(df)[1] = keys[1]
      df
    }),
    local({
      df = data.frame(k = u2, gw = d2, stringsAsFactors = FALSE)
      names(df)[1] = keys[2]
      df
    })
  )
  swapped = FALSE
  if (nrow(means[[1]]) < nrow(means[[2]])) {
    fe = rev(fe)
    keys = rev(keys)
    key_sql = rev(key_sql)
    means = rev(means)
    tmp_idx = idx1
    idx1 = idx2
    idx2 = tmp_idx
    tmp_d = d1
    d1 = d2
    d2 = tmp_d
    tmp_s = s1
    s1 = s2
    s2 = tmp_s
    tmp_y = y1
    y1 = y2
    y2 = tmp_y
    if (length(extra_sums)) {
      tmp_ex = ex1
      ex1 = ex2
      ex2 = tmp_ex
    }
    swapped = TRUE
  }
  C = Matrix::sparseMatrix(
    i = idx1,
    j = idx2, x = ew,
    dims = c(nrow(means[[1]]), nrow(means[[2]]))
  )

  # S = D2 - C' D1^-1 C. The iterative operator never constructs a
  # dense Schur matrix; every vector stays orthogonal to component nulls.
  preconditioner = two_way_cg_preconditioner(
    C, d1, d2, limits,
    max(0, remaining - nrow(edges) * row_bytes)
  )
  if (is.null(preconditioner)) return(NULL)
  project = preconditioner[["project"]]
  precondition = preconditioner[["precondition"]]
  root = preconditioner[["root"]]
  rhs = s2 - as.matrix(Matrix::crossprod(C, y1))
  rhs_scaled = project(rhs / root)
  apply_s = function(value) value - preconditioner[["walk"]](value)
  normalize = function(value) {
    if (equation_metric) sweep(value, 2L, scales, "/") else value / root
  }
  scaled_b = matrix(0, nrow(y2), ncol(y2))
  r = rhs_scaled
  z = project(precondition(r))
  direction = z
  rz = colSums(r * z)
  residual = max(abs(normalize(r)))
  iterations = 0L
  for (iter in seq_len(max_iter)) {
    if (is.finite(residual) && residual <= tol) break
    active = apply(abs(normalize(r)), 2L, max) > tol
    product = project(apply_s(direction))
    denominator = colSums(direction * product)
    if (any(!is.finite(denominator[active])) || any(denominator[active] <= 0)) return(NULL)
    alpha = numeric(ncol(y2))
    alpha[active] = rz[active] / denominator[active]
    if (any(!is.finite(alpha))) return(NULL)
    scaled_b = project(scaled_b + sweep(direction, 2L, alpha, "*"))
    # Recompute the true residual rather than trusting the CG recurrence.
    r = project(rhs_scaled - apply_s(scaled_b))
    residual = max(abs(normalize(r)))
    iterations = iter
    if (verbose) message("[AP] iter ", iter, ": max abs mean = ", sprintf("%.4e", residual))
    if (is.finite(residual) && residual <= tol) break
    z = project(precondition(r))
    new_rz = colSums(r * z)
    beta = numeric(ncol(y2))
    beta[active] = new_rz[active] / rz[active]
    if (any(!is.finite(beta))) return(NULL)
    direction = project(z + sweep(direction, 2L, beta, "*"))
    rz = new_rz
  }
  if (!is.finite(residual) || residual > tol) return(NULL)
  b = scaled_b / root
  a = y1 - as.matrix(C %*% b) / d1
  # Exact cell-level residual verification in R: r1 = s1 - d1*a - C*b, r2 = s2 - C'*a - d2*b
  r1 = s1 - a * d1 - as.matrix(C %*% b)
  r2 = s2 - as.matrix(Matrix::crossprod(C, a)) - b * d2
  norm1 = if (equation_metric) sweep(abs(r1) / sqrt(d1), 2L, scales, "/") else abs(r1) / d1
  norm2 = if (equation_metric) sweep(abs(r2) / sqrt(d2), 2L, scales, "/") else abs(r2) / d2
  actual = max(norm1, norm2)
  if (!is.finite(actual) || actual > tol) return(NULL)
  # Cross-product reduction from FE projection: V' W V - V_tilde' W V_tilde
  # Since V_tilde = V - D1*a - D2*b, V' W V - V_tilde' W V_tilde = 2*(S1'a + S2'b) - (a' D1' D1 a + b' D2' D2 b + 2 a' C b)
  cross_reduction = crossprod(s1, a) + crossprod(s2, b) + crossprod(r1, a) + crossprod(r2, b)
  cross_reduction = 0.5 * (cross_reduction + t(cross_reduction))
  extra_reduction = if (length(extra_sums)) crossprod(a, ex1) + crossprod(b, ex2) else NULL
  # Public reductions refer to the original columns, while the streamed
  # within rows use centered coordinates. Restore the constant contribution.
  total_centered = colSums(s_mat)
  cross_reduction = cross_reduction + outer(centers, total_centered) +
    outer(total_centered, centers) + sum(ew) * tcrossprod(centers)
  if (length(extra_sums)) extra_reduction = extra_reduction + outer(centers, colSums(ex1))
  raw_a = sweep(a, 2L, centers, "+")
  raw_s1 = s1 + outer(d1, centers)
  raw_s2 = s2 + outer(d2, centers)

  if (inline_result && isTRUE(backend_info[["is_databricks"]]) &&
    !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
    lookup_sqls = character(2)
    for (k in 1:2) {
      adjustment = means[[k]][keys[k]]
      effect = if (k == 1L) a else b
      for (j in seq_along(vars)) adjustment[[paste0("a", j)]] = effect[, j]
      lsql = sql_values_frame(conn, adjustment)
      if (is.null(lsql)) return(NULL)
      lookup_sqls[k] = lsql
    }
    w_select = if (is.null(weights)) "1.0 AS __w" else paste0("s.", weights, " AS __w")
    result_sql = paste0(
      "(SELECT /*+ BROADCAST(a, b) */ s.*, ", paste(c(
        w_select,
        sprintf("(%s) - a.a%d - b.a%d AS %s_tilde", vars_sql, seq_along(vars), seq_along(vars), vars)
      ), collapse = ", "),
      " FROM (SELECT * ", from_statement, ") s JOIN (", lookup_sqls[1],
      ") a ON CAST(s.", fe[1], " AS VARCHAR) = a.", keys[1],
      " JOIN (", lookup_sqls[2], ") b ON CAST(s.", fe[2], " AS VARCHAR) = b.", keys[2], ")"
    )
    success = TRUE
    return(list(
      table = result_sql, base_table = NULL, all_tables = character(), view_tables = character(),
      iterations = iterations, max_abs = actual, method = "bounded_pcg", backend_info = backend_info, preconditioner = preconditioner[["method"]],
      components = preconditioner[["components"]],
      cross_reduction = cross_reduction, extra_reduction = extra_reduction, s1 = raw_s1, s2 = raw_s2, a = raw_a, b = b,
      transfer_cells = nrow(edges) * (width + 3L), estimated_bytes = nrow(edges) * row_bytes + 8 * sum(sizes) * max(1L, width) * 12
    ))
  }
  for (k in 1:2) {
    adjustment = means[[k]][keys[k]]
    effect = if (k == 1L) a else b
    for (j in seq_along(vars)) adjustment[[paste0("a", j)]] = effect[, j]
    if (isTRUE(backend_info[["is_databricks"]])) {
      lookup_sql = sql_values_frame(conn, adjustment)
      if (is.null(lookup_sql)) return(NULL)
      create(if (k == 1L) "a" else "b", lookup_sql)
    } else dbWriteTable(conn, tables[[if (k == 1L) "a" else "b"]], adjustment,
      temporary = TRUE, overwrite = TRUE
    )
    state[["created"]] = unique(c(state[["created"]], tables[[if (k == 1L) "a" else "b"]]))
  }
  broadcast_hint = if (isTRUE(backend_info[["is_databricks"]])) "/*+ BROADCAST(a, b) */ " else ""
  w_select = if (is.null(weights)) "1.0 AS __w" else paste0("s.", weights, " AS __w")
  create("result", paste0(
    "SELECT ", broadcast_hint, paste(c(
      paste0("s.", ids), w_select,
      sprintf("(%s) - a.a%d - b.a%d AS %s_tilde", vars_sql, seq_along(vars), seq_along(vars), vars)
    ), collapse = ", "),
    " FROM (SELECT * ", from_statement, ") s JOIN ", tables[["a"]],
    " a ON CAST(s.", fe[1], " AS VARCHAR) = a.", keys[1],
    " JOIN ", tables[["b"]], " b ON CAST(s.", fe[2], " AS VARCHAR) = b.", keys[2]
  ))
  if (!isTRUE(backend_info[["is_databricks"]]) || isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
    checks = vapply(fe, function(f) {
      terms = if (equation_metric) sprintf(
        "FSUM(__w * %s_tilde) / SQRT(FSUM(__w)) / %.17g AS m%d", vars, scales, seq_along(vars)
      ) else sprintf("SUM(__w * %s_tilde) / SUM(__w) AS m%d", vars, seq_along(vars))
      paste0(
        "SELECT MAX(GREATEST(", paste(sprintf("ABS(m%d)", seq_along(vars)), collapse = ", "),
        ")) AS max_abs FROM (SELECT ",
        paste(terms, collapse = ", "),
        " FROM ", tables[["result"]], " GROUP BY ", f, ") q"
      )
    }, character(1))
    actual = backend_get_query(conn, paste0(
      "SELECT MAX(max_abs) AS max_abs FROM (",
      paste(checks, collapse = " UNION ALL "), ") q"
    ))[["max_abs"]]
    n_result = backend_get_query(conn, paste("SELECT COUNT(*) AS n FROM", tables[["result"]]))[["n"]]
    n_source = as.numeric(preflight[["source_n"]])
    if (!is.finite(actual) || actual > tol || n_result != n_source) return(NULL)
  }
  success = TRUE
  list(
    table = tables[["result"]], base_table = NULL, all_tables = state[["created"]], view_tables = intersect(state[["created"]], views),
    iterations = iterations, max_abs = actual, method = "bounded_pcg", backend_info = backend_info, preconditioner = preconditioner[["method"]],
    components = preconditioner[["components"]],
    cross_reduction = cross_reduction, extra_reduction = extra_reduction, s1 = raw_s1, s2 = raw_s2, a = raw_a, b = b,
    transfer_cells = nrow(edges) * (width + 3L), estimated_bytes = nrow(edges) * row_bytes + 8 * sum(sizes) * max(1L, width) * 12
  )
}

# Exact weighted within transformation using the bounded FE solver shared with
# GLMs. Raw rows and key dictionaries stay in SQL, including string/BIGINT keys.
# The result obeys the existing AP table contract used by OLS and IV frames.
dbreg_moment_projections = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights", "cluster_var", "verbose"))
  backend = detect_backend(conn)
  if (!(backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) ||
    length(fe) < 2L || length(fe) > 30L) return(NULL)
  if (isTRUE(backend[["is_databricks"]])) validate_databricks_ap_backend(backend)

  if (!is.null(inputs[["retain"]])) {
    fast = dbreg_fwl_moments(inputs, x_sql, x_names, max_iter, tol)
    if (!is.null(fast)) return(fast)
  }
  limits = inputs[["projection_limits"]]
  if (is.null(limits)) limits = projection_limits()
  max_edges = limits[["summary_rows"]]
  max_entries = limits[["design_entries"]]
  max_cells = limits[["transfer_cells"]]
  validate_moment_limits(max_edges, max_entries)
  validate_moment_limits(max_cells, max_cells)
  vars = c(yvar, x_names)
  columns = paste0("dbreg_m_v", seq_along(vars))
  ids = unique(c(fe, cluster_var, inputs[["retained_key_safe"]]))
  id_columns = paste0("dbreg_m_id", seq_along(ids))
  fe_columns = paste0("dbreg_m_f", seq_along(fe))

  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources))
  manager = fit_table_manager(conn, resources, backend)
  create = manager[["create"]]
  drop = manager[["drop"]]
  w = if (is.null(weights)) "1.0" else sql_cast_double(weights, conn)
  source_sql = paste("SELECT", paste(c(
    sprintf("%s AS %s", ids, id_columns),
    sprintf("%s AS %s", sql_cast_double(c(yvar, x_sql), conn), columns),
    paste(w, "AS w")
  ), collapse = ","), from_statement)
  fe_ids = id_columns[match(fe, ids)]
  stats = backend_get_query(conn, paste("SELECT", paste(c(
    "COUNT(*) AS n", "SUM(w) AS sw",
    sprintf("MAX(ABS(%s)) AS original%d", columns, seq_along(columns)),
    sprintf("COUNT(DISTINCT %s) AS n%d", fe_ids, seq_along(fe)),
    sprintf("SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS null%d", fe_ids, seq_along(fe)),
    sprintf("SUM(w*%s)/SUM(w) AS c%d", columns, seq_along(columns))
  ), collapse = ","), "FROM (", source_sql, ") projection_preflight"))
  sizes = vapply(seq_along(fe), function(i) as.numeric(stats[[paste0("n", i)]]), numeric(1))
  centers = vapply(seq_along(columns), function(i) as.numeric(stats[[paste0("c", i)]]), numeric(1))
  nobs = as.numeric(stats[["n"]])
  if (nobs == 0 || any(!is.finite(centers)) || any(sizes < 1) ||
    any(vapply(seq_along(fe), function(i) as.numeric(stats[[paste0("null", i)]]) > 0, logical(1))) ||
    sum(sizes) * length(columns) > max_entries || !fe_dense_feasible(sizes, length(columns), limits)) {
    projection_abort(
      "system_limit", "Fixed-effect system exceeds its projection budget.",
      list(fe_sizes = sizes, columns = length(columns), estimated_bytes = projection_memory(sizes, length(columns)), limits = limits)
    )
  }
  source = create(source_sql)

  # Rank distinct keys, never the N-row source once per FE. Comparisons happen
  # in their SQL types, so IDs above 2^53 and arbitrary text remain distinct.
  dictionaries = vapply(seq_along(fe), function(i) {
    f = fe_ids[i]
    paste0(
      " JOIN (SELECT ", f, ",ROW_NUMBER() OVER (ORDER BY ", f, ") AS ", fe_columns[i],
      " FROM (SELECT DISTINCT ", f, " FROM ", source, ") keys) d", i, " ON s.", f, "=d", i, ".", f
    )
  }, character(1))
  fields = c(
    paste0("s.", id_columns), "s.w",
    sprintf("(s.%s-(%s)) AS %s", columns, dbglm_number(centers), columns),
    sprintf("d%d.%s", seq_along(fe), fe_columns)
  )
  encoded = create(paste0("SELECT ", paste(fields, collapse = ","), " FROM ", source, " s", paste(dictionaries, collapse = "")))
  moments = query_fe_moments(conn, encoded, fe_columns, columns,
    max_edges = max_edges, max_entries = max_entries, max_cells = max_cells, fe_sizes = sizes
  )
  if (is.null(moments)) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.", list(limits = limits))
  # A stricter solve than the requested absolute group-mean tolerance leaves
  # room for SQL reconstruction roundoff. The SQL certificate below is final.
  target = min(1e-10, tol * sqrt(min(moments[["sys"]][["diag"]])) / 4)
  projected = project_fe_moments(conn, encoded, columns,
    c(id_columns, fe_columns), lapply(fe_columns, function(f) list(f = f, v = "1.0")),
    rep(1, length(columns)), target, list(fe_max_iter = max_iter), create, drop,
    new.env(parent = emptyenv()),
    precomputed = moments, key = NULL
  )
  if (is.null(projected)) return(NULL)

  certificate = dbreg_projection_certificate(conn, projected[["table"]], fe_columns, columns, nobs)
  if (!is.finite(certificate) || certificate > tol) {
    if (verbose) message("[dbreg] Bounded FE projection failed row certification; using iterative AP")
    projection_abort("certificate_failed", "Bounded FE projection failed row certification.", list(residual = certificate, target = tol))
  }
  result = create(paste("SELECT", paste(c(
    sprintf("%s AS %s", id_columns, ids), "w AS __w",
    sprintf("%s AS %s_tilde", columns, vars)
  ), collapse = ","), "FROM", projected[["table"]]))
  # The returned physical table is now owned by the frame. Release every
  # intermediate here, even on Databricks, before the frame begins its passes.
  name = gsub("[\"`]", "", result)
  resources[["tables"]] = setdiff(resources[["tables"]], name)
  if (verbose) message(
    "[dbreg] Bounded FE moments: ", projected[["cycles"]],
    " solve/refinement steps; certified max abs mean = ", sprintf("%.4e", certificate)
  )
  list(
    table = result, base_table = NULL, all_tables = result,
    iterations = projected[["cycles"]], max_abs = certificate, method = "bounded_moments",
    transfer_cells = moments[["transfer_cells"]], estimated_bytes = projection_memory(sizes, length(columns)),
    nobs = nobs, n_fe_levels = sizes,
    original_scale = setNames(vapply(seq_along(columns), function(i) as.numeric(stats[[paste0("original", i)]]), numeric(1)), vars)
  )
}

# Certify the actual transformed SQL rows, rather than only the small system.
# One grouping-sets pass checks group equations, finite values and row count.
dbreg_projection_certificate = function(conn, table, fe, columns, nobs) {
  means = sprintf("ABS(SUM(w*%s)/SUM(w)) AS r%d", columns, seq_along(columns))
  values = c("w", columns)
  invalid = paste(c("w <= 0", paste0(values, " IS NULL OR NOT ", sql_is_finite(conn, values))), collapse = " OR ")
  sql = paste0(
    "SELECT ", paste(c(
      sprintf("MAX(r%d) AS r%d", seq_along(columns), seq_along(columns)),
      "MAX(bad) AS bad", "SUM(n) AS count_sum"
    ), collapse = ","),
    " FROM (SELECT ", paste(c(means, "COUNT(*) AS n", paste0("SUM(CASE WHEN ", invalid, " THEN 1 ELSE 0 END) AS bad")), collapse = ","),
    " FROM ", table, " GROUP BY GROUPING SETS (", paste0("(", fe, ")", collapse = ","), ")) certified"
  )
  row = backend_get_query(conn, sql)
  if (is.na(row[["bad"]]) || row[["bad"]] > 0 ||
    !isTRUE(as.numeric(row[["count_sum"]]) == nobs * length(fe))) return(Inf)
  max(vapply(row[paste0("r", seq_along(columns))], as.numeric, numeric(1)))
}

# General linear nuisance projection, including group-specific slopes. The
# shared projector only solves weighted group equations; no IRLS/family is used.
dbreg_varying_projections = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights", "cluster_var"))
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources), add = TRUE)
  manager = fit_table_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  ids = unique(c(fe, cluster_var, inputs[["retained_key_safe"]]))
  id_names = paste0("__dbreg_vary_id", seq_along(ids))
  columns = paste0("__dbreg_vary_v", seq_along(c(yvar, x_names)))
  trends = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  trend_names = stats::setNames(paste0("__dbreg_vary_t", seq_along(trends)), trends)
  fields = c(
    "ROW_NUMBER() OVER () AS __dbglm_key", paste(ids, "AS", id_names),
    paste(sql_cast_double(c(yvar, x_sql), conn), "AS", columns),
    paste(sql_cast_double(trends, conn), "AS", unname(trend_names)),
    paste(if (is.null(weights)) "1.0" else sql_cast_double(weights, conn), "AS w")
  )
  source = create(paste("SELECT", paste(fields, collapse = ","), from_statement))
  weight_max = as.numeric(backend_get_query(conn, paste("SELECT MAX(w) AS m FROM", source))[["m"]])
  if (!is.finite(weight_max) || weight_max <= 0) stop("No positive-weight observations remain.", call. = FALSE)
  weight_scale = if (weight_max < 1 || weight_max > 2) 2^floor(log2(weight_max)) else 1
  w = paste0("w/(", dbglm_number(weight_scale), ")")
  centers = as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf(
    "SUM((%s)*%s)/SUM(%s) AS c%d", w, columns, w, seq_along(columns)
  ), collapse = ","), "FROM", source))[1, ])
  coded_fe = paste0("__dbreg_vary_f", seq_along(fe))
  attr(coded_fe, "slopes") = lapply(inputs[["fe_slopes"]], function(vars) unname(trend_names[vars]))
  ranked = paste0("(", sql_rank_fixed_effects(conn, source, as.list(id_names[match(fe, ids)]), coded_fe), ")")
  encoded = create(paste("SELECT", paste(c(
    "__dbglm_key", coded_fe, unname(trend_names),
    sprintf("(%s-(%s)) AS %s", columns, dbglm_number(centers), columns), paste(w, "AS w")
  ), collapse = ","), "FROM", ranked))
  control = dbglm_control(fe_tol = min(tol, 1e-12), fe_max_iter = max_iter)
  projected = solve_partial_wls(conn, encoded, columns, coded_fe, control, create, drop, response = FALSE)
  result = create(paste0(
    "SELECT ", paste(c(
      paste0("s.", id_names, " AS ", ids),
      paste0("p.", columns, " AS ", c(yvar, x_names), "_tilde"), "s.w AS __w"
    ), collapse = ","),
    " FROM ", projected[["table"]], " p JOIN ", source, " s ON p.__dbglm_key=s.__dbglm_key"
  ))
  name = gsub("[\"`]", "", result)
  resources[["tables"]] = setdiff(resources[["tables"]], name)
  list(
    table = result, base_table = NULL, all_tables = result,
    iterations = projected[["cycles"]], max_abs = projected[["residual"]], method = projected[["method"]]
  )
}

# Read-only FWL: aggregate raw centered cross-products and FE information,
# solve the compact group equations, and subtract their projection in R.
dbreg_fwl_moments = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights"))
  limits = inputs[["projection_limits"]]
  if (is.null(limits)) limits = projection_limits()
  vars = c(yvar, x_names)
  taken = unique(c(fe, vars, weights, inputs[["cluster_var"]], inputs[["retained_key_safe"]], "w", "__w"))
  columns = safe_identifiers(paste0("dbreg_fwl_v", seq_along(vars)), taken = taken)
  keys = safe_identifiers(paste0("dbreg_fwl_f", seq_along(fe)), taken = c(taken, columns, paste0(columns, "_tilde")))
  w = if (is.null(weights)) "1.0" else sql_cast_double(weights, conn)
  passthrough = unique(c(fe, inputs[["cluster_var"]], inputs[["retained_key_safe"]]))
  id_columns = safe_identifiers(paste0("dbreg_fwl_id", seq_along(passthrough)),
    taken = c(taken, columns, paste0(columns, "_tilde"), keys)
  )
  fe_ids = id_columns[match(fe, passthrough)]
  output_columns = safe_identifiers(paste0(vars, "_tilde"), taken = passthrough)
  weight_column = safe_identifiers("dbreg_fwl_weight__", taken = c(passthrough, output_columns))
  raw = paste("SELECT", paste(c(
    paste(passthrough, "AS", id_columns), paste(sql_cast_double(c(yvar, x_sql), conn), "AS", columns),
    paste(w, "AS w")
  ), collapse = ","), from_statement)
  preflight = backend_get_query(conn, paste("SELECT", paste(c(
    "COUNT(*) AS n", "MIN(w) AS wmin", "MAX(w) AS wmax",
    sprintf("COUNT(DISTINCT %s) AS size%d", fe_ids, seq_along(fe)),
    sprintf("MIN(typeof(%s)) AS type%d", fe_ids, seq_along(fe)),
    sprintf(
      paste0(
        "MAX(CASE WHEN UPPER(typeof(%1$s)) IN ('TINYINT','SMALLINT','INT','INTEGER','BIGINT',",
        "'UTINYINT','USMALLINT','UINTEGER','UBIGINT','DOUBLE','FLOAT','DATE','BOOLEAN') ",
        "THEN 64 ELSE %2$s END) AS bytes%3$d"
      ), fe_ids,
      vapply(fe_ids, function(key) sql_text_bytes(conn, key), character(1)), seq_along(fe)
    ),
    sprintf("%s AS center%d", vapply(columns, function(column) sql_projection_mean(conn, column), character(1)), seq_along(columns)),
    sprintf("MAX(ABS(%s)) AS maximum%d", columns, seq_along(columns))
  ), collapse = ","), "FROM (", raw, ") raw"))
  nobs = as.numeric(preflight[["n"]])
  sizes = vapply(seq_along(fe), function(i) as.numeric(preflight[[paste0("size", i)]]), numeric(1))
  centers = vapply(seq_along(columns), function(i) as.numeric(preflight[[paste0("center", i)]]), numeric(1))
  maximum = vapply(seq_along(columns), function(i) as.numeric(preflight[[paste0("maximum", i)]]), numeric(1))
  if (!nobs || any(!is.finite(c(centers, maximum))) || any(sizes < 1) ||
    sum(sizes) * length(columns) > limits[["design_entries"]] ||
    !fe_dense_feasible(sizes, length(columns), limits)) return(NULL)
  wmax = as.numeric(preflight[["wmax"]])
  weight_scale = 2^floor(log2(wmax))
  if (!is.finite(weight_scale) || weight_scale <= 0 || as.numeric(preflight[["wmin"]]) / weight_scale <= 0) return(NULL)
  key_bytes = vapply(seq_along(fe), function(i) as.numeric(preflight[[paste0("bytes", i)]]), numeric(1))
  row_bytes = 4 * (sum(sort(key_bytes, decreasing = TRUE)[1:2]) + 128 + 8 * (length(columns) + 3L))
  system_bytes = projection_memory(sizes, length(columns))
  max_rows = min(limits[["summary_rows"]], floor((limits[["memory_bytes"]] - system_bytes) / row_bytes))
  if (!is.finite(max_rows) || max_rows < 1) return(NULL)
  key_sql = vapply(seq_along(fe), function(i) {
    type = tolower(as.character(preflight[[paste0("type", i)]]))
    if (type %in% c("double", "float", "real")) {
      paste0("CASE WHEN ", fe_ids[[i]], "=0 THEN CAST(0.0 AS DOUBLE) ELSE ", fe_ids[[i]], " END")
    } else fe_ids[[i]]
  }, character(1))
  source = paste0("(SELECT ", paste(c(
    id_columns, paste(key_sql, "AS", keys),
    sprintf("(%s-(%s)) AS %s", columns, dbglm_number(centers), columns),
    paste0("w/(", dbglm_number(weight_scale), ") AS w")
  ), collapse = ","), " FROM (", raw, ") raw)")
  global_terms = fwl_information_terms(conn, columns)
  moments = query_fe_moments(conn, source, keys, columns,
    max_edges = max_rows,
    max_entries = limits[["design_entries"]], max_cells = limits[["transfer_cells"]], fe_sizes = sizes,
    global_terms = global_terms, text_keys = TRUE
  )
  if (is.null(moments)) return(NULL)
  global = lapply(moments[["globals"]], as.numeric)
  if (!identical(as.numeric(global[["fwl_n"]]), nobs) ||
    !identical(as.numeric(vapply(moments[["info"]], function(info) length(info[["id"]]), integer(1))), sizes)) return(NULL)
  target = min(1e-10, tol * sqrt(min(moments[["sys"]][["diag"]])) / 4)
  cache = new.env(parent = emptyenv())
  projected = project_fe_moments(conn, source, columns, keys,
    lapply(keys, function(key) list(f = key, v = "1.0")), rep(1, length(columns)), target,
    list(fe_max_iter = max_iter),
    create = NULL, drop = NULL, cache = cache,
    precomputed = moments, need_table = FALSE
  )
  raw_moments = dbglm_read_information(global, length(columns))
  M = raw_moments - projected[["cross_reduction"]]
  M = (M + t(M)) / 2
  positive = diag(raw_moments) > 0
  recompute = any(diag(M)[positive] < 1e-8 * diag(raw_moments)[positive])
  if (any(!is.finite(M))) return(NULL)
  need_rows = retain_sample(inputs[["retain"]]) || recompute || !identical(inputs[["vcov_type_req"]], "iid")
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources))
  table = NULL
  if (need_rows) {
    manager = fit_table_manager(conn, resources)
    rows = project_fe_moments(conn, source, columns, id_columns,
      lapply(keys, function(key) list(f = key, v = "1.0")), rep(1, length(columns)), target,
      list(fe_max_iter = max_iter), manager[["create"]], manager[["drop"]], cache,
      precomputed = moments, materialize = FALSE
    )
    projected_source = rows[["raw_projected"]]
    if (recompute) {
      actual = backend_get_query(conn, paste(
        "SELECT",
        paste(fwl_information_terms(conn, paste0(columns, "_tilde")), collapse = ","),
        "FROM", projected_source
      ))
      M = dbglm_read_information(lapply(actual, as.numeric), length(columns))
      if (any(!is.finite(M))) return(NULL)
    }
    table = paste0(
      "(SELECT ", paste(c(
        paste(id_columns, "AS", passthrough),
        paste0("w*(", dbglm_number(weight_scale), ") AS ", weight_column),
        paste(paste0(columns, "_tilde"), "AS", output_columns)
      ), collapse = ","),
      " FROM ", projected_source, ")"
    )
  }
  row = list(n_total = nobs, sum_w = global[["fwl_sw"]], sum_wy = 0, sum_wy_sq = M[1L, 1L])
  for (i in seq_along(fe)) row[[paste0("n_fe", i)]] = sizes[[i]]
  for (i in seq_along(x_names)) {
    row[[sprintf("sum_w%d_y", i)]] = M[i + 1L, 1L]
    for (j in seq_len(i)) row[[sprintf("sum_w%d_%d", j, i)]] = M[j + 1L, i + 1L]
  }
  for (i in seq_along(c(x_names, yvar))) {
    row[[paste0("dbreg_moment_center_", i)]] = 0
    row[[paste0("dbreg_moment_original_", i)]] = maximum[match(c(x_names, yvar)[[i]], vars)]
  }
  row[["dbreg_moment_weight_scale"]] = weight_scale
  if (retain_prediction(inputs[["retain"]])) {
    total = sum(sizes) + length(x_names)
    if (total * max(1, total - max(sizes), length(columns)) <= prediction_state_limit()) {
      inputs[["linear_prediction_moments"]] = list(
        info = lapply(moments[["info"]], function(group) group[c("id", "h")]),
        sys = moments[["sys"]], B = moments[["B"]], alpha = projected[["coefficients"]],
        cross = raw_moments, within = M, centers = centers, x = x_names,
        sizes = sizes, nobs = nobs, sum_w = global[["fwl_sw"]], weight_scale = weight_scale,
        min_weight = as.numeric(preflight[["wmin"]]) / weight_scale,
        types = vapply(seq_along(fe), function(i) as.character(preflight[[paste0("type", i)]]), character(1))
      )
    }
  }
  tables = resources[["tables"]]
  resources[["tables"]] = character()
  list(
    table = table, base_table = NULL, all_tables = tables, inline = TRUE, moments = as.data.frame(row),
    cols = setNames(output_columns, vars), weight_column = weight_column,
    iterations = projected[["cycles"]], max_abs = projected[["residual"]] / sqrt(min(moments[["sys"]][["diag"]])),
    method = "bounded_moments", observation_tables = 0L, transfer_cells = moments[["transfer_cells"]],
    estimated_bytes = system_bytes + moments[["rows"]] * row_bytes, nobs = nobs, n_fe_levels = sizes
  )
}

# Raw and projected moment queries use the same stable scalar aggregations.
fwl_information_terms = function(conn, columns) {
  terms = c("COUNT(*) AS fwl_n", paste(sql_stable_sum(conn, "w"), "AS fwl_sw"))
  for (i in seq_along(columns)) for (j in seq_len(i)) {
    terms = c(terms, paste(sql_stable_sum(
      conn,
      paste0("w*(", columns[[i]], ")*(", columns[[j]], ")")
    ), "AS", sprintf("h_%d_%d", i, j)))
  }
  terms
}

# Shared by iterative execution and the inspectable SQL plan.
sql_ap_group_means = function(table, vars, fe) {
  terms = sprintf("SUM(__w * %s_tilde) / SUM(__w) AS %s_mean", vars, vars)
  paste0(
    "SELECT ", paste(c(fe, terms), collapse = ", "),
    " FROM ", table, " GROUP BY ", fe
  )
}

sql_ap_join_projection = function(table, means, ids, vars, fe) {
  cols = c(
    paste0("t.", ids), "t.__w",
    sprintf("t.%s_tilde - m.%s_mean AS %s_tilde", vars, vars, vars)
  )
  paste0(
    "SELECT ", paste(cols, collapse = ", "), " FROM ", table,
    " t JOIN ", means, " m ON t.", fe, " = m.", fe
  )
}

sql_ap_convergence = function(table, vars, fe) {
  means = paste(sprintf("SUM(__w * %s_tilde) / SUM(__w) AS %s_mean", vars, vars), collapse = ", ")
  inner = paste(sprintf("SELECT %s FROM %s GROUP BY %s", means, table, fe), collapse = " UNION ALL ")
  outer = paste(sprintf("MAX(ABS(%s_mean)) AS max_%s", vars, vars), collapse = ", ")
  paste0("SELECT ", outer, " FROM (", inner, ") dbreg_convergence")
}

sql_ap_original_scale = function(table, vars, conn) {
  paste("SELECT", paste(sprintf("MAX(ABS(%s)) AS original_%d", sql_cast_double(vars, conn), seq_along(vars)), collapse = ","), "FROM", table)
}

sql_ap_rms = function(table, vars, conn) {
  values = sql_cast_double(paste0(vars, "_tilde"), conn)
  paste("SELECT", paste(sprintf("SQRT(SUM(__w*(%s)*(%s))/SUM(__w)) AS rms_%d", values, values, seq_along(vars)), collapse = ","), "FROM", table)
}

dbreg_ap_query_plan = function(inputs, x_sql, x_names) {
  conn = inputs[["conn"]]
  fe = inputs[["fe"]]
  y = inputs[["yvar"]]
  vars = c(y, x_names)
  ids = unique(c(fe, inputs[["cluster_var"]]))
  w = sql_weight_expr(inputs[["weights"]])
  if (is.null(w)) w = "1.0"
  backend = detect_backend(conn)[["name"]]
  method = ap_method(backend)
  base = paste0(
    "SELECT ", paste(c(frame_select(y, x_sql, x_names,
      id_cols = setdiff(ids, vars)
    ), paste0(w, " AS __w")), collapse = ", "),
    " ", inputs[["from_statement"]]
  )
  tilde = paste0(vars, "_tilde")
  initial = paste0(
    "SELECT ", paste(c(ids, "__w", sprintf(
      "%s AS %s",
      if (method == "window") sql_cast_double(vars, conn) else vars, tilde
    )), collapse = ", "),
    " FROM (<model_frame_sql>) dbreg_ap_base"
  )
  projection = setNames(vapply(fe, function(f) {
    if (method == "window") return(sql_ap_window_sweeps("<ap_current>", c(ids, "__w"), tilde, f, 1L))
    paste0(
      "WITH dbreg_means AS (", sql_ap_group_means("<ap_current>", vars, f), ") ",
      sql_ap_join_projection("<ap_current>", "dbreg_means", ids, vars, f)
    )
  }, character(1)), fe)
  frame = list(
    cte_sql = "WITH demeaned AS (SELECT * FROM <ap_converged>)", cte_name = "demeaned",
    cols = setNames(tilde, vars), weights_expr = "__w", intercept = FALSE, cast = TRUE,
    rank_source = paste0("(", base, ") dbreg_ap_rank_source"),
    prefix_terms = c(sql_count(conn, "n_total"), vapply(seq_along(fe), function(k)
      sql_count(conn, sprintf("n_fe%d", k), fe[[k]], distinct = TRUE), character(1)))
  )
  structure(list(
    strategy = "demean", backend = backend, method = method,
    tolerance = inputs[["fe_tol"]], max_iterations = inputs[["fe_max_iter"]],
    model_frame_sql = base, initialization_sql = initial,
    design_matrix = data.frame(term = restore_names(standardize_coef_names(x_names), inputs[["name_map"]]), sql = x_sql),
    projection_sql = projection, convergence_sql = sql_ap_convergence("<ap_current>", vars, fe),
    original_scale_sql = sql_ap_original_scale("(<model_frame_sql>) dbreg_ap_original", vars, conn),
    residual_scale_sql = sql_ap_rms("<ap_current>", vars, conn),
    final_moments_sql = frame_moments_sql(frame, y, x_names, conn),
    note = paste(
      "Templates for iterative fixed-effect projection; replace the angle-bracket placeholders.",
      "Repeat the projections until convergence. The automatic solver may instead use bounded FE moments.",
      "Also require group means divided by projected RMS to meet tolerance for columns above their original-value roundoff bound.",
      "No projection or coefficient estimation is executed by sql_only."
    )
  ), class = "dbreg_query_plan")
}

#' Inspect a fixed-effect projection SQL plan
#' @param x A plan returned by `dbreg()` or `dbivreg()` with `sql_only = TRUE`.
#' @param ... Unused.
#' @export
print.dbreg_query_plan = function(x, ...) {
  cat("Fixed-effect SQL plan (", x[["backend"]], ", ", x[["method"]], ")\n", sep = "")
  for (name in c("model_frame_sql", "initialization_sql", "projection_sql", "convergence_sql", "original_scale_sql", "residual_scale_sql", "final_moments_sql")) {
    cat("\n", name, ":\n", paste(x[[name]], collapse = "\n\n"), "\n", sep = "")
  }
  cat("\n", x[["note"]], "\n", sep = "")
  invisible(x)
}

# SQL owns row state; bounded FE moments or FE-level CG determine projections.
solve_partial_wls = function(conn, table, x, fe, control, create, drop, response = TRUE, cache = NULL, materialize = TRUE, scales = NULL, extra_sums = character(), need_table = TRUE, precomputed = NULL) {
  columns = c(x, if (response) "z")
  backend_info = detect_backend(conn)
  slopes = attr(fe, "slopes")
  if (is.null(cache)) cache = new.env(parent = emptyenv())
  if (length(fe) > 1) {
    dbglm_moments_enabled(conn, table, fe, control, cache, length(columns))
    if (!cache[["moment_enabled"]] && length(fe) == 2L && !length(unlist(slopes)) &&
      identical(control[["projection_method"]], "auto") && length(columns)) {
      identifiers = unique(c(fe, unlist(slopes)))
      if (is.null(scales)) {
        norm = as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf(
          "FSUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)
        ), collapse = ","), "FROM", table))[1, ])
        scales = sqrt(pmax(norm, .Machine[["double.xmin"]]))
      }
      target = control[["fe_tol"]]
      if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
      first = utils::tail(columns, 1)
      rest = utils::head(columns, -1)
      order = match(c(first, rest), columns)
      attempt = projection_attempt("bounded_pcg", dbreg_duckdb_cg_projections(conn, paste("FROM", table), fe, first, rest, rest, "w",
        cluster_var = "__dbglm_key", max_iter = control[["fe_max_iter"]], tol = target, scales = scales[order],
        materialize_result = FALSE, inline_result = !materialize, extra_sums = extra_sums
      ))
      projected = record_projection_attempt(cache, attempt)
      if (!is.null(projected)) {
        on.exit(for (t in rev(projected[["all_tables"]])) {
          if (!materialize && isTRUE(backend_info[["is_databricks"]]) && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) next
          if (t %in% projected[["view_tables"]]) try(DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", t)), silent = TRUE) else
            drop_ap_table_if_exists(conn, t, backend_info)
        }, add = TRUE)
        if (!materialize && isTRUE(backend_info[["is_databricks"]]) && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
          cache[["views"]] = unique(c(cache[["views"]], projected[["view_tables"]]))
          cache[["tables"]] = unique(c(cache[["tables"]], setdiff(projected[["all_tables"]], projected[["view_tables"]])))
        }
        result = create(paste("SELECT __dbglm_key,", paste(c(
          identifiers, "__w AS w",
          sprintf("%s_tilde AS %s", columns, columns)
        ), collapse = ","), "FROM", projected[["table"]]), materialize = materialize)
        inv_order = match(columns, c(first, rest))
        reduction = if (!is.null(projected[["cross_reduction"]])) projected[["cross_reduction"]][inv_order, inv_order, drop = FALSE] else NULL
        extra_red = if (!is.null(projected[["extra_reduction"]])) projected[["extra_reduction"]][inv_order, , drop = FALSE] else NULL
        raw_proj = if (!materialize && isTRUE(backend_info[["is_databricks"]]) &&
          !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) projected[["table"]] else NULL
        return(list(
          table = result, raw_projected = raw_proj, cycles = projected[["iterations"]], residual = projected[["max_abs"]], method = "bounded_pcg",
          cross_reduction = reduction, extra_reduction = extra_red, sub_tables = if (!materialize) projected[["all_tables"]] else character(),
          s1 = projected[["s1"]][, inv_order, drop = FALSE], s2 = projected[["s2"]][, inv_order, drop = FALSE],
          a = projected[["a"]][, inv_order, drop = FALSE], b = projected[["b"]][, inv_order, drop = FALSE]
        ))
      }
    }
    if (cache[["moment_enabled"]]) {
      ops = lapply(fe, function(f) list(f = f, v = "1.0"))
      for (i in seq_along(fe)) for (v in slopes[[i]]) ops[[length(ops) + 1]] = list(f = fe[i], v = v)
      identifiers = unique(c(fe, unlist(slopes)))
      if (is.null(scales)) {
        norm = if (length(columns)) as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf("FSUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)), collapse = ","), "FROM", table))[1, ]) else numeric()
        scales = sqrt(pmax(norm, .Machine[["double.xmin"]]))
      }
      target = control[["fe_tol"]]
      if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
      attempt = projection_attempt("bounded_moments", project_fe_moments(conn, table, columns, identifiers, ops, scales, target, control, create, drop, cache,
        materialize = materialize, extra_sums = extra_sums, need_table = need_table, precomputed = precomputed
      ))
      ans = record_projection_attempt(cache, attempt)
      if (!is.null(ans)) return(ans)
      cache[["moment_enabled"]] = FALSE
    }
    if (isTRUE(backend_info[["is_databricks"]]) && !length(unlist(slopes))) {
      return(dbglm_project_window(conn, table, columns, fe, control, create, drop, response))
    }
  }
  projection = create(paste("SELECT __dbglm_key,", paste(unique(c(fe, unlist(slopes), "w", columns)), collapse = ","), "FROM", table))
  ops = lapply(fe, function(f) list(f = f, v = "1.0"))
  basis = character()
  for (i in seq_along(fe)) {
    local_basis = character()
    for (j in seq_along(slopes[[i]])) {
      v = slopes[[i]][j]
      name = paste0("__dbglm_v", i, "_", j)
      partition = paste0(" OVER (PARTITION BY ", fe[i], ")")
      centered = paste0(v, " - SUM(w*", v, ")", partition, "/SUM(w)", partition)
      next_projection = create(paste("SELECT *,", centered, "AS", name, "FROM", projection))
      drop(projection)
      projection = next_projection
      for (prev in local_basis) {
        adjustment = paste0(
          name, " - ", prev, "*COALESCE(SUM(w*", prev, "*", name, ")", partition,
          "/NULLIF(SUM(w*", prev, "*", prev, ")", partition, ",0),0)"
        )
        next_projection = create(paste("SELECT * EXCLUDE(", name, "),", adjustment, "AS", name, "FROM", projection))
        drop(projection)
        projection = next_projection
      }
      local_basis = c(local_basis, name)
      basis = c(basis, name)
      ops[[length(ops) + 1]] = list(f = fe[i], v = name)
    }
  }
  identifiers = unique(c(fe, unlist(slopes), basis))
  if (!length(columns)) return(list(table = projection, cycles = 0L, residual = 0, method = "group"))
  norm = backend_get_query(conn, paste("SELECT", paste(sprintf("SUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)), collapse = ","), "FROM", projection))
  scales = sqrt(pmax(as.numeric(norm[1, ]), .Machine[["double.xmin"]]))
  target = control[["fe_tol"]]
  if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
  if (length(fe) > 1) return(dbglm_project_cg(conn, projection, columns, identifiers, ops, scales, target, control, create, drop))
  residual = 0
  for (cycle in seq_len(control[["fe_max_iter"]])) {
    for (op in ops) {
      f = op[["f"]]
      v = op[["v"]]
      values = sprintf("(%1$s - (%3$s)*COALESCE(SUM(w*(%3$s)*%1$s) OVER (PARTITION BY %2$s)/NULLIF(SUM(w*(%3$s)*(%3$s)) OVER (PARTITION BY %2$s),0),0)) AS %1$s", columns, f, v)
      next_projection = create(paste("SELECT __dbglm_key,", paste(c(identifiers, "w", values), collapse = ","), "FROM", projection))
      drop(projection)
      projection = next_projection
    }
    residual = 0
    for (op in ops) {
      f = op[["f"]]
      v = op[["v"]]
      terms = sprintf("COALESCE(ABS(SUM(w*(%3$s)*%1$s))/NULLIF(SQRT(SUM(w*(%3$s)*(%3$s))),0),0)/%2$s AS r%4$d", columns, dbglm_number(scales), v, seq_along(columns))
      row = backend_get_query(conn, paste0(
        "SELECT ", paste(sprintf("MAX(r%d) AS r%d", seq_along(columns), seq_along(columns)), collapse = ","),
        " FROM (SELECT ", paste(terms, collapse = ","), " FROM ", projection, " GROUP BY ", f, ") t"
      ))
      residual = max(residual, as.numeric(row[1, ]))
    }
    if (residual <= target) return(list(table = projection, cycles = cycle, residual = residual, method = "group"))
  }
  stop("GLM fixed-effect projection did not converge.", call. = FALSE)
}

# Direct FE moments need the non-absorbed levels to fit a dense solve in R;
# see fe_dense_feasible(). The decision is cached for the whole fit.
dbglm_moments_enabled = function(conn, table, fe, control, cache, columns = 1L) {
  if (!is.null(cache[["moment_enabled"]])) return(cache[["moment_enabled"]])
  if (identical(control[["projection_method"]], "sql")) {
    cache[["moment_enabled"]] = FALSE
    return(FALSE)
  }
  if (is.null(cache[["projection_limits"]])) cache[["projection_limits"]] = projection_limits()
  sizes = cache[["fe_sizes"]]
  if (is.null(sizes)) {
    count = backend_get_query(conn, paste("SELECT", paste(sprintf("COUNT(DISTINCT %s) AS n%d", fe, seq_along(fe)), collapse = ","), "FROM", table))
    sizes = vapply(count, as.numeric, numeric(1))
  }
  cache[["fe_sizes"]] = sizes
  op_sizes = rep(sizes, 1 + lengths(attr(fe, "slopes")))
  cache[["moment_enabled"]] = sum(op_sizes) * max(1L, columns) <= cache[["projection_limits"]][["design_entries"]] &&
    fe_dense_feasible(op_sizes, columns, limits = cache[["projection_limits"]])
  if (!cache[["moment_enabled"]]) record_projection_attempt(
    cache,
    list(
      method = "bounded_moments", status = "declined", reason = "system_limit",
      details = list(fe_sizes = op_sizes, estimated_bytes = projection_memory(op_sizes, columns), limits = cache[["projection_limits"]])
    )
  )
  cache[["moment_enabled"]]
}

# Preconditioned CG on FE coefficients. Every row correction is reconstructed
# from group lookups, so roundoff cannot introduce arbitrary row-space effects.
# State size is O(number of FE coefficients * number of projected columns).
dbglm_project_cg = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop) {
  index = seq_along(columns)
  r = paste0("__cg_r", index)
  d = paste0("__cg_d", index)
  g = paste0("__cg_g", index)
  a = paste0("__cg_a", index)
  direction_names = paste0("__cg_value", index)
  key = "__cg_op,__cg_group"
  moments = function(source, inputs, outputs, info = NULL) {
    pieces = vapply(seq_along(ops), function(i) {
      op = ops[[i]]
      terms = c(
        if (is.null(info)) sprintf("SQRT(FSUM(w*(%1$s)*(%1$s))) AS __cg_scale", op[["v"]]),
        sprintf("FSUM(w*(%s)*(%s)) AS u%d", op[["v"]], inputs, index)
      )
      paste0(
        "SELECT ", i, " AS __cg_op,", op[["f"]], " AS __cg_group,", paste(terms, collapse = ","),
        " FROM ", source, " GROUP BY ", op[["f"]]
      )
    }, character(1))
    denominator = if (is.null(info)) "u.__cg_scale" else "i.__cg_scale"
    fields = sprintf("COALESCE(u.u%d/NULLIF(%s,0),0) AS %s", index, denominator, outputs)
    create(paste0(
      "SELECT u.__cg_op,u.__cg_group,", denominator, " AS __cg_scale,", paste(fields, collapse = ","),
      " FROM (", paste(pieces, collapse = " UNION ALL "), ") u",
      if (!is.null(info)) paste0(" JOIN ", info, " i USING (", key, ")")
    ))
  }
  rows = function(state, coefficients, outputs, subtract = FALSE) {
    ctes = vapply(seq_along(ops), function(i) paste0(
      "cg", i, " AS (SELECT __cg_group,__cg_scale,",
      paste(coefficients, collapse = ","), " FROM ", state, " WHERE __cg_op=", i, ")"
    ), character(1))
    fields = vapply(index, function(j) {
      terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
        sprintf("(%s)*COALESCE(cg%d.%s/NULLIF(cg%d.__cg_scale,0),0)", v, i, coefficients[j], i)
      }, character(1))
      value = paste(terms, collapse = " + ")
      if (subtract) value = paste0("s.", columns[j], "-(", value, ")")
      paste0("(", value, ") AS ", outputs[j])
    }, character(1))
    joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN cg", i, " ON s.", ops[[i]][["f"]], "=cg", i, ".__cg_group"), character(1))
    create(paste0(
      "WITH ", paste(ctes, collapse = ","), " SELECT ",
      paste(paste0("s.", c("__dbglm_key", identifiers, "w")), collapse = ","), ",", paste(fields, collapse = ","),
      " FROM ", table, " s", paste(joins, collapse = "")
    ))
  }
  info = moments(table, columns, r)
  state = create(paste("SELECT *,", paste(c(paste(r, "AS", d), paste("0.0 AS", g)), collapse = ","), "FROM", info))
  active = rep(TRUE, length(columns))
  floor_rr = (32 * .Machine[["double.eps"]] * scales)^2
  for (cycle in seq_len(control[["fe_max_iter"]])) {
    direction = rows(state, d, direction_names)
    applied = moments(direction, direction_names, a, info)
    drop(direction)
    joined = create(paste0("SELECT s.*,", paste(paste0("a.", a), collapse = ","), " FROM ", state, " s JOIN ", applied, " a USING (", key, ")"))
    moment = backend_get_query(conn, paste("SELECT", paste(c(
      sprintf("FSUM(%1$s*%1$s) AS rr%2$d", r, index),
      sprintf("FSUM(%s*%s) AS da%d", d, a, index)
    ), collapse = ","), "FROM", joined))
    rr = vapply(index, function(j) as.numeric(moment[[paste0("rr", j)]]), numeric(1))
    da = vapply(index, function(j) as.numeric(moment[[paste0("da", j)]]), numeric(1))
    if (any(!is.finite(rr)) || any(!is.finite(da)) || any(active & da <= 0 & rr > pmax((target * scales)^2, floor_rr))) {
      stop("GLM fixed-effect projection lost numerical information.", call. = FALSE)
    }
    alpha = ifelse(active & da > 0 & rr > floor_rr, rr / da, 0)
    updated = create(paste("SELECT", paste(c(
      key, "__cg_scale",
      sprintf("(%s+(%s)*%s) AS %s", g, dbglm_number(alpha), d, g),
      sprintf("(%s-(%s)*%s) AS %s", r, dbglm_number(alpha), a, r), d
    ), collapse = ","), "FROM", joined))
    drop(joined)
    drop(applied)
    drop(state)
    state = updated
    restart = cycle %% 20L == 0L || all(alpha == 0)
    if (cycle %% 4L == 0L || cycle == 1L || cycle == control[["fe_max_iter"]] || restart) {
      projected = rows(state, g, columns, subtract = TRUE)
      certificate = moments(projected, columns, r, info)
      row = backend_get_query(conn, paste("SELECT", paste(sprintf("MAX(ABS(%s))/%s AS r%d", r, dbglm_number(scales), index), collapse = ","), "FROM", certificate))
      residual_columns = as.numeric(row[1, ])
      residual = max(residual_columns)
      active = residual_columns > target
      if (isTRUE(getOption("dbreg.glm_trace")) && cycle %% 100L == 0L) message("[projection] cycle ", cycle, " residual ", residual, " target ", target)
      if (residual <= target) {
        drop(certificate)
        drop(state)
        drop(info)
        drop(table)
        return(list(table = projected, cycles = cycle, residual = residual, method = "cg"))
      }
      drop(projected)
      if (restart) {
        restarted = create(paste0(
          "SELECT s.__cg_op,s.__cg_group,s.__cg_scale,", paste(paste0("s.", g), collapse = ","), ",",
          paste(c(paste0("c.", r, " AS ", r), paste0("c.", r, " AS ", d)), collapse = ","),
          " FROM ", state, " s JOIN ", certificate, " c USING (", key, ")"
        ))
        drop(state)
        state = restarted
      }
      drop(certificate)
    }
    if (!restart) {
      new_rr = as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf("FSUM(%1$s*%1$s) AS r%2$d", r, index), collapse = ","), "FROM", state))[1, ])
      beta = ifelse(active & rr > floor_rr, new_rr / rr, 0)
      next_state = create(paste(
        "SELECT * EXCLUDE(", paste(d, collapse = ","), "),",
        paste(sprintf("((%s)*%s+(%s)*%s) AS %s", dbglm_number(as.numeric(active)), r, dbglm_number(beta), d, d), collapse = ","), "FROM", state
      ))
      drop(state)
      state = next_state
    }
  }
  stop("GLM fixed-effect projection did not converge.", call. = FALSE)
}

# Reuse Spark's bounded-plan window AP chunks with the GLM's scaled equation
# certificate and its tighter response-dependent inner tolerance.
dbglm_project_window = function(conn, table, columns, fe, control, create, drop, response) {
  if (!length(columns)) return(list(table = create(paste(
    "SELECT __dbglm_key,",
    paste(c(fe, "w"), collapse = ","), "FROM", table
  )), cycles = 0L, residual = 0, method = "window"))
  norm = backend_get_query(conn, paste("SELECT", paste(sprintf(
    "SUM(w*POWER(%s,2)) AS n%d",
    columns, seq_along(columns)
  ), collapse = ","), "FROM", table))
  scales = sqrt(pmax(as.numeric(norm[1, ]), .Machine[["double.xmin"]]))
  target = control[["fe_tol"]]
  if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
  first = utils::tail(columns, 1)
  rest = utils::head(columns, -1)
  vars_all = c(first, rest)
  tilde_names = paste0(vars_all, "_tilde")
  id_cols = unique(c(fe, "__dbglm_key"))
  tilde_cols = sprintf("CAST((%s) AS DOUBLE) AS %s", vars_all, tilde_names)
  cur_table = create(paste0("SELECT ", paste(c(id_cols, "w AS __w", tilde_cols), collapse = ", "), " FROM ", table))
  chunk_cycles = getOption("dbreg.spark_ap_chunk_cycles", 3L)
  iter = 0L
  max_abs = Inf
  while (iter < control[["fe_max_iter"]] && max_abs > target) {
    cycles = min(chunk_cycles, control[["fe_max_iter"]] - iter)
    iter = iter + cycles
    window_sql = sql_ap_window_sweeps(cur_table, c(id_cols, "__w"), tilde_names, fe, cycles)
    alt_table = create(window_sql)
    drop(cur_table)
    cur_table = alt_table
    diagnostics = vapply(fe, function(fe_k) {
      aliases = sprintf("dbreg_mean_%04d", seq_along(tilde_names))
      means = sprintf("ABS(SUM(__w * %s))/SQRT(SUM(__w))/%.17g AS %s", tilde_names, scales[match(vars_all, columns)], aliases)
      paste0(
        "SELECT MAX(GREATEST(", paste(aliases, collapse = ", "),
        ")) AS max_abs FROM (SELECT ", fe_k, ", ", paste(means, collapse = ", "),
        " FROM ", cur_table, " GROUP BY ", fe_k, ") dbreg_ap_conv"
      )
    }, character(1))
    convergence_query = paste0("SELECT GREATEST(", paste(sprintf("(%s)", diagnostics), collapse = ", "), ") AS max_abs")
    max_abs = as.numeric(backend_get_query(conn, convergence_query)[1, 1L])
  }
  result = create(paste("SELECT __dbglm_key,", paste(c(
    fe, "__w AS w",
    sprintf("%s_tilde AS %s", columns, columns)
  ), collapse = ","), "FROM", cur_table))
  drop(cur_table)
  list(table = result, cycles = iter, residual = max_abs, method = "window")
}

# Per-fit controls stay in inputs rather than mutating global R options.
validate_fit_diagnostics = function(inputs) {
  tol = inputs[["fe_tol"]]
  iterations = inputs[["fe_max_iter"]]
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`fe_tol` must be one finite, positive number.", call. = FALSE)
  }
  if (!is.numeric(iterations) || length(iterations) != 1L || !is.finite(iterations) ||
    iterations < 1 || iterations > .Machine[["integer.max"]] || iterations != floor(iterations)) {
    stop("`fe_max_iter` must be a positive integer.", call. = FALSE)
  }
  diagnostics = inputs[["diagnostics"]]
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  inputs[["fe_max_iter"]] = as.integer(iterations)
}

# This optional scan is separate from estimation; missingness counts can overlap.
sample_diagnostics = function(inputs, from_statement, variables) {
  if (!isTRUE(inputs[["diagnostics"]]) || isTRUE(inputs[["sql_only"]])) return(NULL)
  roles = rep(names(variables), lengths(variables))
  variables = unlist(variables, use.names = FALSE)
  aliases = paste0("n_missing_", seq_along(variables))
  terms = c(
    sql_count(inputs[["conn"]], "n_input"),
    sprintf("SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS %s", variables, aliases)
  )
  weights = inputs[["weights"]]
  if (!is.null(weights)) terms = c(
    terms,
    sprintf("SUM(CASE WHEN %s <= 0 THEN 1 ELSE 0 END) AS n_nonpositive_weight", weights)
  )
  sql = paste0(
    "SELECT ", paste(terms, collapse = ", "),
    " FROM (SELECT * ", from_statement, ") dbreg_audit"
  )
  counts = DBI::dbGetQuery(inputs[["conn"]], sql)
  list(
    n_input = as.numeric(counts[["n_input"]]),
    n_nonpositive_weight = if (is.null(weights)) 0 else as.numeric(counts[["n_nonpositive_weight"]]),
    detail = data.frame(
      role = roles,
      variable = restore_names(variables, inputs[["source_map"]]),
      n_missing = vapply(aliases, function(a) as.numeric(counts[[a]]), numeric(1)), row.names = NULL
    ),
    query_string = sql
  )
}

finalize_fit_diagnostics = function(result, inputs) {
  estimates = setNames(result[["coeftable"]][, "estimate"], rownames(result[["coeftable"]]))
  dropped = setdiff(result[["collin.var"]], names(estimates))
  result[["collin.coef"]] = c(estimates, setNames(rep(NA_real_, length(dropped)), dropped))
  if (is.null(result[["coef_names_all"]])) {
    visible = result[["coef_names"]]
    if (is.null(visible)) visible = result[["xvars"]]
    result[["coef_names_all"]] = unique(c(visible, dropped))
  }
  result[["vcov_type"]] = inputs[["vcov_type_req"]]
  result[["cluster_var"]] = restore_names(inputs[["cluster_var"]], inputs[["source_map"]])
  result[["ssc"]] = if (inputs[["vcov_type_req"]] == "cluster") inputs[["ssc"]] else NULL
  result[["n_clusters"]] = attr(result[["vcov"]], "n_clusters")
  result[["n_params"]] = attr(result[["vcov"]], "n_params")
  result[["sample_diagnostics"]] = inputs[["sample_diagnostics"]]
  result[["nobs_input"]] = if (is.null(result[["sample_diagnostics"]])) NA_real_ else result[["sample_diagnostics"]][["n_input"]]
  result[["nobs_dropped"]] = result[["nobs_input"]] - result[["nobs_orig"]]
  projected = !is.null(result[["ap_iterations"]])
  result[["ap_converged"]] = if (projected) TRUE else NA
  result[["ap_tolerance"]] = if (projected) inputs[["fe_tol"]] else NA_real_
  result[["ap_max_iter"]] = if (projected) inputs[["fe_max_iter"]] else NA_integer_
  result[["projection_diagnostics"]] = inputs[["projection_diagnostics"]]
  result
}
