# Apply identical parity gates to saved fits without repeating database work.
# label identifies the independent reference; sparse checks are not fixest runs.
teacher_wide_compare = function(bundle, ref, out, label = "fixest") {
  fit = bundle[["fit"]]
  keys = bundle[["keys"]]
  n = bundle[["n_input"]]
  elapsed = bundle[["elapsed"]]
  mode = bundle[["engine"]]
  stopifnot(n == ref[["n_input"]], bundle[["exits"]] == ref[["exits"]])
  checks = c(names = identical(names(coef(fit)), names(ref[["coef"]])),
    nobs = nobs(fit) == ref[["nobs"]], sample = identical(keys, sort(ref[["keys"]])),
    coef = all(abs(coef(fit) - ref[["coef"]]) <= 1e-7 + 1e-6 * abs(ref[["coef"]])),
    vcov = all(abs(vcov(fit) - ref[["vcov"]]) <= 1e-8 + 1e-5 * abs(ref[["vcov"]])),
    se = all(abs(sqrt(diag(vcov(fit))) - sqrt(diag(ref[["vcov"]]))) <= 1e-8 + 1e-6 * sqrt(diag(ref[["vcov"]]))),
    likelihood = abs(as.numeric(logLik(fit)) - ref[["logLik"]]) / nobs(fit) < 1e-9,
    model_df = fit[["n_params_model"]] == ref[["model_df"]],
    covariance_df = fit[["n_params_vcov"]] == ref[["covariance_df"]],
    score = fit[["score_max"]] <= fit[["control"]][["score_tol"]],
    projection = fit[["projection_residual"]] <= fit[["control"]][["fe_tol"]], converged = fit[["converged"]])
  write.csv(data.frame(check = names(checks), passed = unname(checks)), file.path(out, paste0(mode, "_", label, "_checks.csv")), row.names = FALSE)
  metrics = data.frame(engine = mode, reference = label, n_input = n, nobs = nobs(fit), exits = ref[["exits"]],
    coefficients = length(coef(fit)), elapsed = elapsed, coef_error = max(abs(coef(fit) - ref[["coef"]])),
    vcov_error = max(abs(vcov(fit) - ref[["vcov"]])),
    se_error = max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(ref[["vcov"]])))),
    likelihood_error = abs(as.numeric(logLik(fit)) - ref[["logLik"]]),
    score = fit[["score_max"]], projection = fit[["projection_residual"]],
    iterations = fit[["iterations"]], passed = all(checks))
  write.csv(metrics, file.path(out, paste0(mode, "_", label, "_metrics.csv")), row.names = FALSE)
  print(metrics)
  stopifnot(all(checks))
  invisible(metrics)
}
