dbglm_likelihood = function(conn, sample, x, fit, family) {
  n = sample[["n"]]
  base = sample[["base"]]
  weight_scale = sample[["weight_scale"]]
  deviance = fit[["deviance"]]
  weighted = switch(family[["family"]],
    gaussian = dbglm_gaussian_loglik(n, deviance, weight_scale,
      sample[["moments"]][["log_weights"]] + n * log(weight_scale)),
    poisson = (dbglm_saturated_loglik(conn, sample[["original_base"]]) - deviance / 2) * weight_scale,
    binomial = -deviance * weight_scale / 2)
  reported = switch(family[["family"]],
    gaussian = weighted - 1,
    poisson = weighted,
    binomial = {
      working = fit[["working"]]
      if (is.null(working)) {
        rows = dbglm_working_sql(base, x, fit[["beta"]], family[["link"]], conn = conn, cell = sample[["cell"]])
        working = paste0("(", rows, "SELECT * FROM working) likelihood")
      }
      dbglm_binomial_loglik(conn, working, weight_scale, cell = sample[["cell"]])
    })
  list(reported = reported, weighted = weighted)
}

dbglm_null_fit = function(conn, sample, result, intercept) {
  family = result[["family"]]
  control = result[["control"]]
  base = sample[["original_base"]]
  if (!intercept && length(sample[["y_center"]]) && sample[["y_center"]] != 0) {
    base = paste0("(SELECT ", sql_select_except(conn, "__dbglm_y"),
      ", (__dbglm_y+(", dbglm_number(sample[["y_center"]]), ")) AS __dbglm_y FROM ", base, ") null_uncentered")
  }
  sample[["base"]] = base
  sample[["cell"]] = FALSE
  moments = sample[["moments"]]
  x = if (intercept) "__dbglm_intercept" else character()
  if (intercept && length(sample[["fe"]])) {
    base = paste0("(SELECT *,1.0 AS __dbglm_intercept FROM ", base, ") null_sample")
    sample[["base"]] = base
  }
  full_is_null = identical(result[["coef_names"]], "(Intercept)")
  record = list(definition = "global_same_sample_weights_offset_structural_constant",
    nobs = result[["nobs_orig"]], constant = intercept, offset = result[["offset"]],
    family = family[["family"]], link = family[["link"]])
  state = if (!result[["finite_mle"]]) {
    list(available = FALSE, reason = "nonfinite_full_mle")
  } else if (full_is_null && result[["converged"]]) {
    list(available = TRUE, method = "full_is_null", iterations = 0L, evaluations = 0L,
      likelihood = list(reported = result[["loglik"]], weighted = result[["weighted_loglik"]]))
  } else tryCatch({
    if (intercept) {
      alpha = dbglm_initial_intercept(conn, base, family, moments,
        offset_present = isTRUE(sample[["offset_present"]]))
      control[["tolerate_max_iter"]] = FALSE
      fit = dbglm_irls(conn, base, x, alpha, family, control, anchor = alpha)
      method = "scoring_constant"
      iterations = fit[["iterations"]]
      evaluations = 1L + iterations + fit[["step_halvings"]]
    } else {
      fit = dbglm_evaluate(conn, base, x, numeric(), family[["link"]])
      if (fit[["bad"]] || !is.finite(fit[["deviance"]]) ||
          dbglm_boundary(fit, family[["link"]])) {
        stop("Null model reaches a numerical boundary.", call. = FALSE)
      }
      fit[["beta"]] = numeric()
      method = "fixed_offset"
      iterations = 0L
      evaluations = 1L
    }
    list(available = TRUE, method = method, iterations = iterations,
      evaluations = evaluations, likelihood = dbglm_likelihood(conn, sample, x, fit, family))
  }, error = function(error) {
    message = conditionMessage(error)
    warning("Requested null model unavailable: ", message, call. = FALSE)
    reason = if (grepl("boundary", message)) "null_boundary" else
      if (grepl("did not converge within|step halving failed", message)) "budget_exceeded" else
        "numerical_error"
    list(available = FALSE, reason = reason, message = message)
  })
  record[["fit"]] = state[setdiff(names(state), "likelihood")]
  for (kind in c("reported", "weighted")) {
    ll = state[["likelihood"]][[kind]]
    if (is.null(ll)) ll = NA_real_
    available = isTRUE(state[["available"]]) && is.finite(ll)
    full = result[[if (kind == "reported") "loglik" else "weighted_loglik"]]
    reason = if (!isTRUE(state[["available"]])) state[["reason"]] else
      if (!is.finite(ll)) "nonfinite_null_likelihood" else NULL
    ratio_reason = if (!available) reason else
      if (!is.finite(full)) "nonfinite_full_likelihood" else
        if (ll == 0) "zero_null_likelihood" else NULL
    ratio = if (is.null(ratio_reason)) 1 - full / ll else NA_real_
    if (is.null(ratio_reason) && !is.finite(ratio)) ratio_reason = "nonfinite_ratio"
    record[[kind]] = list(logLik = ll, available = available, reason = reason,
      ratio = ratio, ratio_available = is.null(ratio_reason), ratio_reason = ratio_reason)
  }
  record
}
