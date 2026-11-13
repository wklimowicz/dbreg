# Aggregate only compact sufficient statistics; every read uses the owned base.
fe_information_query = function(conn, source, columns, budget, limits, bounds = NULL) {
  fe_information_reserve(columns, budget, limits)
  terms = fe_information_terms(conn, columns, bounds, detect_backend(conn))
  fe_information_read(DBI::dbGetQuery(conn, paste("SELECT", paste(terms, collapse = ","), "FROM", source)), columns, budget)
}

fe_information_reserve = function(columns, budget, limits) {
  width = length(columns)
  cells = width * (width + 1) / 2 + 2
  if (budget[["cells"]] + cells > limits[["transfer_cells"]] || budget[["rows"]] + 1 > limits[["summary_rows"]]) {
    projection_abort("transfer_limit", "Global FWL moments exceed the transfer budget.")
  }
  invisible(cells)
}

# Global moment terms; callers may evaluate them together with other scalar
# aggregates of the same source and read them back with fe_information_read.
fe_information_terms = function(conn, columns, bounds = NULL, backend = detect_backend(conn)) {
  width = length(columns)
  bounded = is.list(bounds) && is.numeric(bounds[["weight"]]) &&
    length(bounds[["weight"]]) == 1L && is.finite(bounds[["weight"]]) &&
    bounds[["weight"]] > 0 && is.numeric(bounds[["columns"]]) &&
    length(bounds[["columns"]]) == width &&
    all(is.finite(bounds[["columns"]]) & bounds[["columns"]] >= 0)
  terms = c("COUNT(*) AS n", paste(sql_stable_sum(conn, "w", backend_info = backend,
    max_abs = if (bounded) bounds[["weight"]] else NULL), "AS sw"))
  for (i in seq_len(width)) for (j in seq_len(i)) {
    bound = if (bounded) bounds[["weight"]] * bounds[["columns"]][i] * bounds[["columns"]][j] else NULL
    terms = c(terms, paste(sql_stable_sum(conn, paste0("w*(", columns[i], ")*(", columns[j], ")"),
      backend_info = backend, max_abs = bound), "AS", paste0("h", i, "_", j)))
  }
  terms
}

fe_information_read = function(row, columns, budget) {
  width = length(columns)
  cells = width * (width + 1) / 2 + 2
  if (nrow(row) != 1L) stop("Incomplete fixed-effect global moments.", call. = FALSE)
  row = lapply(row[c("n", "sw", unlist(lapply(seq_len(width), function(i) paste0("h", i, "_", seq_len(i)))))], as.numeric)
  if (any(lengths(row) != 1L) || any(!is.finite(unlist(row))) || row[["n"]] < 1 || row[["sw"]] <= 0) {
    projection_abort("singular", "Nonfinite fixed-effect global moments.")
  }
  cross = matrix(0, width, width)
  for (i in seq_len(width)) for (j in seq_len(i)) cross[i, j] = cross[j, i] = row[[paste0("h", i, "_", j)]]
  budget[["cells"]] = budget[["cells"]] + cells
  budget[["rows"]] = budget[["rows"]] + 1
  list(cross = cross, n = row[["n"]], sum_w = row[["sw"]])
}

fe_fwl_information = function(conn, source, columns, B, coefficients, sys, budget, limits, global = NULL, bounds = NULL) {
  # Include all small dense cross-product buffers in the existing memory plan.
  if (projection_memory(sys[["sizes"]], length(columns), limits = limits) + 64 * length(columns)^2 > limits[["memory_bytes"]]) {
    projection_abort("system_limit", "FWL cross-products exceed the projection memory budget.")
  }
  if (is.null(global)) global = fe_information_query(conn, source, columns, budget, limits, bounds)
  residual = B - fe_multiply(sys, coefficients)
  reduction = crossprod(B, coefficients) + crossprod(residual, coefficients)
  reduction = (reduction + t(reduction)) / 2
  global[["within"]] = global[["cross"]] - reduction
  global
}

fe_fwl_moments = function(conn, projected, columns, information, scale, original_scale, weight_scale, nobs, sizes, budget, limits) {
  if (!identical(information[["n"]], nobs)) stop("Fixed-effect source count changed.", call. = FALSE)
  M = information[["within"]]
  raw = information[["cross"]]
  if (any(!is.finite(M))) projection_abort("singular", "Nonfinite FWL projection information.")
  positive = diag(raw) > 0
  # Subtracting nearly equal cross-products can destroy a small residual.
  # Recompute the actual projected products in SQL, without writing rows.
  recompute = any(diag(M) < 0) || any(diag(M)[positive] < 1e-8 * diag(raw)[positive])
  if (recompute) {
    expressions = paste0(columns, "/(", sprintf("%.17e", scale), ")")
    actual = fe_information_query(conn, projected, expressions, budget, limits)
    if (!identical(actual[["n"]], nobs)) stop("Fixed-effect projected count changed.", call. = FALSE)
    M = actual[["cross"]]
  }
  M = M * outer(scale, scale)
  if (any(!is.finite(M))) projection_abort("singular", "Nonfinite original-scale FWL moments.")
  row = list(n_total = nobs, sum_w = information[["sum_w"]], sum_wy = 0, sum_wy_sq = M[1L, 1L])
  for (i in seq_along(sizes)) row[[paste0("n_fe", i)]] = sizes[i]
  for (i in seq_len(length(columns) - 1L)) {
    row[[paste0("sum_w", i, "_y")]] = M[i + 1L, 1L]
    row[[paste0("dbreg_moment_original_", i)]] = original_scale[i + 1L]
    for (j in seq_len(i)) row[[paste0("sum_w", j, "_", i)]] = M[j + 1L, i + 1L]
  }
  for (i in seq_along(columns)) row[[paste0("dbreg_moment_center_", i)]] = 0
  row[["dbreg_moment_weight_scale"]] = weight_scale
  out = as.data.frame(row)
  attr(out, "fwl_recomputed") = recompute
  out
}
