retain_linear_score_frame = function(fit, frame, y, inputs) {
  if (!retain_sample(inputs[["retain"]])) return(invisible(NULL))
  conn = inputs[["conn"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  residual = fit_residual_sql(fit, frame, y, conn)
  x = vapply(fit[["x"]], function(name) if (name == "(Intercept)") "1.0" else
    fit_column_sql(fit, frame, name, conn), character(1))
  w = frame[["weights_expr"]]
  if (is.null(w)) w = "1.0"
  scale = fit[["weight_scale"]]
  if (is.null(scale)) scale = 1
  w = paste0("(", w, ")/(", dbglm_number(scale), ")")
  scores = sprintf("CAST((%s) AS DOUBLE)*(%s)*(%s)", w, residual, x)
  score_names = if (length(x)) paste0("dbreg_fit_score", seq_along(x)) else character()
  select = c(paste(inputs[["retained_key_safe"]], "AS dbreg_fit_key__"),
    paste(residual, "AS dbreg_fit_residual__"), if (length(scores)) paste(scores, "AS", score_names))
  query = paste(frame[["cte_sql"]], "SELECT", paste(select, collapse = ","), "FROM", frame[["cte_name"]])
  original = inputs[["retained_source"]]
  public_y = restore_names(y, inputs[["source_map"]])
  weights = restore_names(inputs[["weights"]], inputs[["source_map"]])
  components = lapply(inputs[["prediction_fe_components"]], function(columns) restore_names(unique(columns), inputs[["source_map"]]))
  state = retain_score_sample(conn, original, inputs[["sample_resources"]], query, "dbreg_fit_key__",
    if (length(scores)) paste0("s.", score_names) else character(),
    paste0("CAST(r.", q(public_y), " AS DOUBLE)-s.dbreg_fit_residual__"), "s.dbreg_fit_residual__",
    row_id = inputs[["row_id"]], fe_components = components, weights = weights)
  if (!is.finite(state[["rss"]]) || abs(state[["rss"]] - fit[["rss"]]) > 1e-6*max(1, abs(fit[["rss"]]))) {
    stop("Retained OLS scores do not reproduce the fitted RSS.", call. = FALSE)
  }
  state[["bread"]] = unname(fit[["bread"]])*scale
  if (!is.null(fit[["sandwich_bread"]])) state[["score_projection"]] = t(unname(fit[["sandwich_bread"]])*scale)
  state[["coefficient_names"]] = fit[["x"]]
  state[["kind"]] = "ols"
  state[["weight_scale"]] = scale
  inputs[["inference_state"]] = state
  invisible(NULL)
}

finalize_linear_inference = function(result, inputs) {
  if (!retain_sample(inputs[["retain"]])) return(result)
  if (is.null(inputs[["inference_state"]]) && identical(result[["strategy"]], "compress")) {
    reconstruct_compressed_score_frame(result, inputs)
  }
  state = inputs[["inference_state"]]
  if (is.null(state)) stop("This strategy has not retained its inference frame.", call. = FALSE)
  if (state[["nobs"]] != result[["nobs_orig"]]) stop("Retained score rows differ from the fitted sample.", call. = FALSE)
  map = inputs[["name_map"]]
  variable = names(map) %in% names(inputs[["source_map"]])
  map[variable] = coef_var_name(map[variable])
  names = standardize_coef_names(restore_names(state[["coefficient_names"]], map))
  if (!identical(names, names(coef.dbreg(result)))) stop("Retained inference coefficients do not match the fitted design.", call. = FALSE)
  state[["coefficient_names"]] = names
  dimnames(state[["bread"]]) = list(names, names)
  state[["K_model"]] = length(names) + if (length(state[["fe_sizes"]])) sum(state[["fe_sizes"]]) - length(state[["fe_sizes"]]) + 1 else 0
  state[["default_ssc"]] = dbglm_ssc(list(K.fixef = if (inputs[["ssc"]] == "full") "full" else "nonnested"))
  state[["default_vcov"]] = if (inputs[["vcov_type_req"]] == "cluster")
    stats::reformulate(restore_names(inputs[["cluster_var"]], inputs[["source_map"]]), env = baseenv()) else inputs[["vcov_type_req"]]
  result = attach_inference_state(result, state)
  result
}

# Compression keeps the original slopes fixed. Reuse the shared within frame
# to reconstruct original-row scores; no cell expansion is fetched into R.
reconstruct_compressed_score_frame = function(result, inputs) {
  conn = inputs[["conn"]]
  design = sql_model_matrix(inputs[["fml"]], conn,
    sub("^FROM\\s+", "", inputs[["from_statement"]], ignore.case = TRUE),
    fe_vars = inputs[["fe"]], intercept = inputs[["intercept"]], factor_refs = inputs[["factor_refs"]])
  register_factor_display(inputs, design[["factor_levels"]], design[["level_aliases"]], design[["col_names"]],
    design[["aliases"]], factor_classes = design[["factor_classes"]])
  map = inputs[["name_map"]]
  variable = names(map) %in% names(inputs[["source_map"]])
  map[variable] = coef_var_name(map[variable])
  public = standardize_coef_names(restore_names(design[["aliases"]], map))
  beta = coef.dbreg(result)
  has_intercept = "(Intercept)" %in% names(beta)
  selected = match(setdiff(names(beta), "(Intercept)"), public)
  if (anyNA(selected)) stop("Could not align compressed fitted scores.", call. = FALSE)
  x_sql = design[["select_exprs"]][selected]
  x = design[["aliases"]][selected]
  frame = if (length(inputs[["fe"]])) build_demean_frame(inputs, x_sql, x) else build_base_frame(inputs, x_sql, x)
  on.exit(frame[["cleanup"]](), add = TRUE)
  row = backend_get_query(conn, frame_moments_sql(frame, inputs[["yvar"]], x, conn))
  S = cross_moments(row, x, inputs[["yvar"]], intercept = has_intercept)
  names = c(if (has_intercept) "(Intercept)", x)
  fixed = matrix(unname(beta), ncol = 1L, dimnames = list(names, NULL))
  H = S[names, names, drop = FALSE]
  fit = list(coef = fixed, x = names,
    bread = if (length(names)) {
      A = moment_coordinates(S, names)
      A %*% solve_with_fallback(H, matrix(0, length(names), 1L))[["XtX_inv"]] %*% t(A)/attr(S, "weight_scale")
    } else matrix(numeric(), 0L, 0L),
    weight_scale = attr(S, "weight_scale"), centers = attr(S, "centers"),
    sandwich_bread = if (length(names)) moment_coordinates(S, names) %*% solve_with_fallback(H, matrix(0, length(names), 1L))[["XtX_inv"]]/attr(S, "weight_scale") else matrix(numeric(), 0L, 0L),
    rss = attr(result[["vcov"]], "rss"))
  retain_linear_score_frame(fit, frame, inputs[["yvar"]], inputs)
}

attach_inference_state = function(result, inference) {
  state = result[["fitted_state"]]
  if (is.null(state)) state = list(version = 1L, specification = NULL, coefficients = NULL,
    sample = list(nobs = inference[["nobs"]]), inference = NULL, resources = NULL, prediction = NULL)
  state[["inference"]] = inference
  state[["resources"]] = list(storage = "database", version = 1L,
    tables = unique(c(state[["resources"]][["tables"]], inference[["resources"]][["tables"]])))
  result[["fitted_state"]] = state
  result
}
