# Shared projection policy. Numerical tolerances remain estimator controls;
# these limits bound representations and allocations, not statistical accuracy.
projection_limits = function() {
  limits = list(
    summary_rows = getOption("dbreg.duckdb_ap_max_edges", 1e6),
    design_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6),
    transfer_cells = getOption("dbreg.ap_moment_max_cells", 8e6),
    dense_levels = getOption("dbreg.glm_fe_dense_max", 4096),
    cross_entries = getOption("dbreg.glm_fe_cross_max", 5e7),
    memory_bytes = getOption("dbreg.projection_max_bytes", 512*1024^2))
  for (name in names(limits)) {
    value = limits[[name]]
    lower = if (name %in% c("dense_levels", "cross_entries")) 0 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < lower || value > 2^53-2 || value != floor(value)) {
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
    class = c(paste0("dbreg_projection_", kind), "dbreg_projection_error", "error", "condition")))
}

projection_attempt = function(method, expr) {
  tryCatch({
    result = force(expr)
    details = if (is.null(result)) list() else result[intersect(names(result), c("transfer_cells", "estimated_bytes", "iterations", "max_abs"))]
    list(method = method, status = if (is.null(result)) "declined" else "completed",
      reason = if (is.null(result)) "not_applicable" else NULL, details = details, result = result)
  }, dbreg_projection_error = function(e) list(method = method, status = "declined",
    reason = e[["reason"]], details = e[["details"]], result = NULL))
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
  rows = min(limits[["summary_rows"]]-budget[["rows"]],
    floor((limits[["transfer_cells"]]-budget[["cells"]])/width)-1L)
  if (rows < 1L) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  result = backend_get_query(conn, paste(sql, "LIMIT", sprintf("%.0f", rows+1L)))
  if (nrow(result) > rows) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  budget[["rows"]] = budget[["rows"]]+nrow(result)
  budget[["cells"]] = budget[["cells"]]+width*(nrow(result)+1L)
  result
}
