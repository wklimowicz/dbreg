# Bounded algebra policy and shared additive fixed-effect system solver.
projection_limits = function() {
  memory_bytes = getOption("dbreg.projection_max_bytes", 1024^3)
  # By default compact FE summaries may use an eighth of the memory budget at
  # 64 planned bytes per pair row (projection_memory()), i.e. about 2.1M rows
  # at 1 GiB; a fixed 1e6 declined real-like scaled designs (Finding J).
  default_rows = if (is.numeric(memory_bytes) && length(memory_bytes) == 1L && is.finite(memory_bytes) && memory_bytes >= 512) floor(memory_bytes / 512) else 1e6
  limits = list(
    summary_rows = getOption("dbreg.projection_max_rows", getOption("dbreg.duckdb_ap_max_edges", default_rows)),
    design_entries = getOption("dbreg.projection_max_entries", getOption("dbreg.duckdb_ap_max_entries", 2e6)),
    transfer_cells = getOption("dbreg.projection_max_cells", getOption("dbreg.ap_moment_max_cells", 16e6)),
    dense_levels = getOption("dbreg.fe_dense_max", getOption("dbreg.glm_fe_dense_max", 4096)),
    cross_entries = getOption("dbreg.fe_cross_max", getOption("dbreg.glm_fe_cross_max", 5e7)),
    upload_bytes = getOption("dbreg.projection_max_upload_bytes", 8 * 1024^2),
    memory_bytes = memory_bytes
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

# Two families whose smaller one exceeds the dense limit use an iterative
# Schur solve (fe_sparse_solver) that never forms the Schur matrix.
fe_iterative_schur = function(sizes, limits = projection_limits()) {
  length(sizes) == 2L && sum(sizes) - max(sizes) > limits[["dense_levels"]]
}

projection_memory = function(sizes, columns = 1L, n_edges = NULL, limits = projection_limits()) {
  rest = sum(sizes) - max(sizes)
  if (is.null(n_edges)) {
    possible = (sum(sizes)^2 - sum(sizes^2)) / 2
    n_edges = min(possible, limits[["summary_rows"]])
  }
  # Sparse cross blocks and scaled copies, plus a conservative dense Schur
  # factor/basis allowance (or CG vectors) and RHS buffers. This plans
  # allocations, not RSS.
  schur = if (fe_iterative_schur(sizes, limits)) 16 * rest * max(columns, 1L) else 6 * rest^2
  64 * n_edges + 8 * (schur + 4 * sum(sizes) * max(columns, 1L))
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
    if (isTRUE(context[["verbose"]]) && identical(summary[["status"]], "declined")) {
      message(projection_decline_message(summary, context[["verbose_prefix"]]))
    }
  }
  invisible(attempt[["result"]])
}

# One line per declined projection method: its reason and the limits or sizes
# behind it (a decline is otherwise silent and can cost hours of fallback).
projection_decline_message = function(summary, prefix = NULL) {
  details = summary[["details"]]
  scalar = function(x) is.numeric(x) && length(x) >= 1L && length(x) <= 8L && all(is.finite(x))
  fmt = function(x) paste(vapply(x, function(v) format(signif(v, 4), scientific = abs(v) >= 1e7 || (v != 0 && abs(v) < 1e-3), trim = TRUE), character(1)), collapse = "/")
  parts = character()
  for (name in setdiff(names(details), "limits")) {
    if (scalar(details[[name]])) parts = c(parts, paste0(name, "=", fmt(details[[name]])))
  }
  limits = details[["limits"]]
  if (is.list(limits) && length(limits)) {
    parts = c(parts, paste0("limits: ", paste(paste0(names(limits), "=", vapply(limits, fmt, character(1))), collapse = ", ")))
  }
  paste0(if (is.null(prefix)) "[dbreg]" else prefix, " ", summary[["method"]], " declined (", summary[["reason"]], ")",
    if (length(parts)) paste0(": ", paste(parts, collapse = "; ")))
}


fe_dense_feasible = function(op_sizes, columns = 1L, limits = projection_limits(), n_edges = NULL) {
  op_sizes = as.numeric(op_sizes)
  if (!length(op_sizes) || any(!is.finite(op_sizes))) return(FALSE)
  rest = sum(op_sizes) - max(op_sizes)
  if (is.null(n_edges)) n_edges = min((sum(op_sizes)^2 - sum(op_sizes^2)) / 2, limits[["summary_rows"]])
  (rest <= limits[["dense_levels"]] || fe_iterative_schur(op_sizes, limits)) &&
    (is.null(n_edges) || n_edges <= limits[["cross_entries"]]) &&
    projection_memory(op_sizes, columns, n_edges, limits) <= limits[["memory_bytes"]]
}

# info: per-op list(id, h); pairs: list of list(i, j, ir, jr, h) with local
# level indices into info[[i]]$id and info[[j]]$id.
fe_system = function(info, pairs, limits = projection_limits()) {
  sizes = vapply(info, function(x) length(x[["id"]]), integer(1))
  n_edges = sum(vapply(pairs, function(pair) length(pair[["h"]]), integer(1)))
  if (!fe_dense_feasible(sizes, n_edges = n_edges, limits = limits)) return(NULL)
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

# Two families: pin one smaller-family level per connected component, then
# factor the Schur complement (sparse or dense Cholesky), or with S = NULL
# solve it by preconditioned CG without ever forming it.
fe_sparse_solver = function(sys, S, scaled, cache, retain_null) {
  if (length(sys[["sizes"]]) != 2L) return(NULL)
  nr = length(sys[["rest"]])
  iterative = is.null(S)
  sparse = !iterative && nr > 512L && Matrix::nnzero(S) / nr^2 < 0.15
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
  if (iterative) {
    limits = if (is.null(cache[["projection_limits"]])) projection_limits() else cache[["projection_limits"]]
    d2 = sys[["diag"]][sys[["rest"]]]
    pre = two_way_cg_preconditioner(sys[["Cm"]], sys[["Dg"]], d2, limits, limits[["memory_bytes"]], components = rest_components)
    if (is.null(pre)) projection_abort("singular", "Fixed-effect information is singular.")
    stats = new.env(parent = emptyenv())
    stats[["iterations"]] = 0L
    solve = function(rhs) {
      r1 = rhs[sys[["big"]], , drop = FALSE]
      r2 = rhs[sys[["rest"]], , drop = FALSE] - as.matrix(Matrix::crossprod(scaled, r1))
      # Scaled Schur system (I - walk) u = r2 / sqrt(d2), u = sqrt(d2) * x2.
      b = pre[["project"]](r2 / pre[["root"]])
      tol = 1e-12 * pmax(apply(abs(b), 2L, max), .Machine[["double.xmin"]])
      solved = fe_pcg(function(v) v - pre[["walk"]](v), b, pre[["precondition"]], pre[["project"]],
        function(r) apply(abs(r), 2L, max) / tol, 1,
        if (is.null(cache[["pcg_max_iter"]])) 2000L else cache[["pcg_max_iter"]])
      if (is.null(solved)) projection_abort("nonconverged", "Fixed-effect CG solve broke down or exhausted its iteration budget.")
      stats[["iterations"]] = stats[["iterations"]] + solved[["iterations"]]
      x2 = solved[["x"]] / pre[["root"]]
      out = matrix(0, sys[["n"]], ncol(rhs))
      out[sys[["rest"]], ] = x2
      out[sys[["big"]], ] = r1 / sys[["Dg"]] - as.matrix(scaled %*% x2)
      out
    }
  } else {
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
    mem_limit = if (is.null(cache[["projection_limits"]])) projection_limits()[["memory_bytes"]] else cache[["projection_limits"]][["memory_bytes"]]
    if (8 * sys[["n"]] * count > mem_limit) {
      projection_abort("system_limit", "Fixed-effect null basis exceeds the prediction memory budget.")
    }
    null_q = matrix(0, sys[["n"]], count)
    sizes = tabulate(components)
    null_q[cbind(seq_len(sys[["n"]]), components)] = sign / sqrt(sizes[components])
  }
  list(solve = solve, project = project, null_basis = null_q,
    stats = if (iterative) stats else NULL, preconditioner = if (iterative) pre[["method"]] else NULL,
    method = if (iterative) "schur_pcg" else if (sparse) "sparse_cholesky" else "dense_cholesky")
}

# Numerical rank from unit-scaled pivoted-Cholesky pivots. Exact FE
# dependencies (e.g. a band nested in one FE and determined by another) should
# pivot to zero, but plain SQL sums over millions of rows (Spark) leave them at
# ~1e-10, above the fixed tolerance. Such a pivot sits many orders of magnitude
# below the identified ones, so a drop of >= 1e6 below sqrt(eps) ends the rank.
fe_pivot_rank = function(pivots) {
  rank = length(pivots)
  if (rank < 2L) return(rank)
  gap = which(pivots[-1L] <= sqrt(.Machine[["double.eps"]]) & pivots[-1L] * 1e6 <= pivots[-rank])
  if (length(gap)) gap[[1L]] else rank
}

# Factor the Schur complement on the non-absorbed levels; the identified basis
# is chosen once per fit so FE normalisations stay fixed across iterations.
fe_solver = function(sys, cache, retain_null = FALSE) {
  if (any(!is.finite(sys[["Dg"]])) || any(sys[["Dg"]] <= 0)) projection_abort("singular", "Fixed-effect information is singular.")
  scaled = sys[["Cm"]]
  if (inherits(scaled, "CsparseMatrix")) {
    scaled@x = scaled@x / sys[["Dg"]][scaled@i + 1L]
  } else scaled = scaled / sys[["Dg"]]
  if (fe_iterative_schur(sys[["sizes"]], if (is.null(cache[["projection_limits"]])) projection_limits() else cache[["projection_limits"]])) {
    return(fe_sparse_solver(sys, NULL, scaled, cache, retain_null))
  }
  S = Matrix::forceSymmetric(sys[["E"]] - Matrix::crossprod(sys[["Cm"]], scaled))
  sparse = fe_sparse_solver(sys, S, scaled, cache, retain_null)
  if (!is.null(sparse)) return(sparse)
  S = as.matrix(S)
  sdiag = sqrt(pmax(diag(S), 0))
  if (is.null(cache[["moment_basis"]])) {
    # A level the big family absorbs has a Schur diagonal that is roundoff
    # relative to its own information; unit scaling below would otherwise
    # promote that roundoff to a full pivot.
    positive = which(is.finite(sdiag) & sdiag^2 > 1e-11 * sys[["diag"]][sys[["rest"]]])
    keep = integer()
    if (length(positive)) {
      C = S[positive, positive, drop = FALSE] / outer(sdiag[positive], sdiag[positive])
      pivoted = suppressWarnings(chol(C, pivot = TRUE, tol = 1e-11))
      rank = fe_pivot_rank(diag(pivoted)[seq_len(attr(pivoted, "rank"))]^2)
      keep = positive[attr(pivoted, "pivot")[seq_len(rank)]]
    }
    if (!length(keep) && !length(sys[["big"]])) stop("No identified fixed-effect span remains.", call. = FALSE)
    cache[["moment_basis"]] = keep
  }
  keep = cache[["moment_basis"]]
  scale = sdiag[keep]
  factor = if (length(keep)) tryCatch(chol(S[keep, keep, drop = FALSE] / outer(scale, scale)), error = function(e) NULL) else matrix(0, 0, 0)
  if (is.null(factor) || any(!is.finite(factor))) projection_abort("singular", "Fixed-effect information is singular.")
  # B lies in range(H) exactly, but plain SQL sums (e.g. Spark) leave roundoff
  # in the redundant FE directions that no solve can remove. Null vectors of
  # H follow from the dropped Schur pivots: v2[keep] = -S_kk^{-1} S_kd.
  dropped = setdiff(seq_along(sys[["rest"]]), keep)
  null_q = NULL
  scaled_null_q = NULL
  information_scale = sqrt(pmax(sys[["diag"]], .Machine[["double.xmin"]]))
  if (length(dropped)) {
    mem_limit = if (is.null(cache[["projection_limits"]])) projection_limits()[["memory_bytes"]] else cache[["projection_limits"]][["memory_bytes"]]
    if (32 * sys[["n"]] * length(dropped) > mem_limit) {
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
    if (retain_null) null_q = qr.Q(qr(V, LAPACK = TRUE))
    # Remove SQL roundoff in information-scaled coordinates. An unscaled
    # projection can move O(eps) from a large group into an O(1e-20) equation,
    # corrupting that group's learned effect by many orders of magnitude.
    scaled_null_q = qr.Q(qr(information_scale * V, LAPACK = TRUE))
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
  list(solve = solve, project = project, null_basis = null_q, method = "dense_schur")
}


fe_pcg = function(apply_op, rhs, precondition, project, residuals, tol, max_iter, trace = NULL) {
  x = matrix(0, nrow(rhs), ncol(rhs))
  r = rhs
  z = project(precondition(r))
  direction = z
  rz = colSums(r * z)
  current = residuals(r)
  iterations = 0L
  for (iter in seq_len(max_iter)) {
    if (all(is.finite(current)) && max(current) <= tol) break
    active = current > tol
    product = project(apply_op(direction))
    denominator = colSums(direction * product)
    if (any(!is.finite(denominator[active])) || any(denominator[active] <= 0)) return(NULL)
    alpha = numeric(ncol(rhs))
    alpha[active] = rz[active] / denominator[active]
    if (any(!is.finite(alpha))) return(NULL)
    x = project(x + sweep(direction, 2L, alpha, "*"))
    r = project(r - sweep(product, 2L, alpha, "*"))
    current = residuals(r)
    iterations = iter
    if (!is.null(trace)) trace(iter, max(current))
    refreshed = FALSE
    if (all(is.finite(current)) && max(current) <= tol) {
      r = project(rhs - apply_op(x))
      current = residuals(r)
      if (all(is.finite(current)) && max(current) <= tol) break
      refreshed = TRUE
    }
    z = project(precondition(r))
    new_rz = colSums(r * z)
    if (refreshed) {
      direction = z
    } else {
      beta = numeric(ncol(rhs))
      beta[active] = new_rz[active] / rz[active]
      if (any(!is.finite(beta))) return(NULL)
      direction = project(z + sweep(direction, 2L, beta, "*"))
    }
    rz = new_rz
  }
  if (!all(is.finite(current)) || max(current) > tol) return(NULL)
  list(x = x, r = r, residual = max(current), iterations = iterations)
}

# Work in D2-scaled coordinates so component deflation remains stable for
# unequal group weights. Bounded Schur fill permits symmetric Gauss-Seidel;
# otherwise use an additive block SGS/Jacobi preconditioner without that fill.
two_way_cg_preconditioner = function(C, d1, d2, limits, available_bytes, components = NULL) {
  n1 = nrow(C)
  n2 = ncol(C)
  if (is.null(components)) {
    system = list(
      n = n1 + n2, big = seq_len(n1), rest = n1 + seq_len(n2),
      Cm = C, Dg = d1, diag = c(d1, d2), sizes = c(n1, n2)
    )
    components = fe_two_components(system)
    if (is.null(components)) return(NULL)
    components = components[n1 + seq_len(n2)]
  }
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
