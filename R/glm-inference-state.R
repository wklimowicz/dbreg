retain_glm_score_frame = function(result, fit, conn, base, x, fe, original, resources,
  components, row_id, weights, weight_scale, counts, ssc) {
  if (isFALSE(result[["finite_mle"]])) stop("Full inference retention requires a finite GLM fit.", call. = FALSE)
  prefix = ""
  if (!length(fe)) {
    prefix = dbglm_working_sql(base, x, fit[["beta"]], result[["family"]][["link"]])
    from = "FROM working b"
    projected = paste0("b.", x)
  } else if (!is.null(fit[["raw_projected"]])) {
    from = paste("FROM", fit[["raw_projected"]], "b")
    projected = paste0("b.", x, "_tilde")
  } else {
    from = paste("FROM", fit[["working"]], "b", if (length(x)) paste("JOIN", fit[["projected"]], "p USING (__dbglm_key)") else "")
    projected = paste0("p.", x)
  }
  score_names = if (length(x)) paste0("dbreg_fit_score", seq_along(x)) else character()
  terms = c("b.__dbglm_key AS dbreg_fit_key__", "b.eta AS dbreg_fit_eta__", "b.mu AS dbreg_fit_mu__",
    "(b.__dbglm_y-b.mu) AS dbreg_fit_residual__", if (length(x)) paste0(projected, "*b.score AS ", score_names))
  query = paste(prefix, "SELECT", paste(terms, collapse = ","), from)
  state = retain_score_sample(conn, original, resources, query, "dbreg_fit_key__",
    if (length(x)) paste0("s.", score_names) else character(), "s.dbreg_fit_eta__", "s.dbreg_fit_residual__",
    mu = "s.dbreg_fit_mu__", row_id = row_id, fe_components = components, weights = weights)
  if (state[["nobs"]] != result[["nobs_orig"]]) stop("Retained GLM scores differ from the fitted sample.", call. = FALSE)
  state[["bread"]] = if (length(x)) dbglm_solve(fit[["H"]])[["inverse"]] else matrix(numeric(), 0L, 0L)
  state[["coefficient_names"]] = names(stats::coef(structure(result, class = c("dbglm", "dbreg"))))
  dimnames(state[["bread"]]) = rep(list(state[["coefficient_names"]]), 2L)
  state[["kind"]] = if (result[["family"]][["family"]] == "gaussian") "gaussian" else "glm"
  state[["dispersion_scaled"]] = result[["dispersion_scaled"]]
  state[["weight_scale"]] = weight_scale
  state[["K_model"]] = counts[["model"]]
  state[["fe_df"]] = fe_intercept_df(state[["fe_sizes"]])
  state[["default_ssc"]] = ssc
  state[["default_vcov"]] = if (result[["vcov_type"]] == "cluster") stats::reformulate(result[["cluster_var"]], env = baseenv()) else result[["vcov_type"]]
  attach_inference_state(result, state)
}
