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
  terms = c(sql_count(inputs[["conn"]], "n_input"),
    sprintf("SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS %s", variables, aliases))
  weights = inputs[["weights"]]
  if (!is.null(weights)) terms = c(terms,
    sprintf("SUM(CASE WHEN %s <= 0 THEN 1 ELSE 0 END) AS n_nonpositive_weight", weights))
  sql = paste0("SELECT ", paste(terms, collapse = ", "),
    " FROM (SELECT * ", from_statement, ") dbreg_audit")
  counts = DBI::dbGetQuery(inputs[["conn"]], sql)
  list(
    n_input = as.numeric(counts[["n_input"]]),
    n_nonpositive_weight = if (is.null(weights)) 0 else as.numeric(counts[["n_nonpositive_weight"]]),
    detail = data.frame(role = roles,
      variable = restore_names(variables, inputs[["source_map"]]),
      n_missing = vapply(aliases, function(a) as.numeric(counts[[a]]), numeric(1)), row.names = NULL),
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
  result
}
