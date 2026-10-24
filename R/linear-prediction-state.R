# Use the original linear fit's slopes and RSS; only nuisance coefficients are
# recovered. Both model classes share the bounded FE kernel and prediction state.
finalize_linear_prediction = function(result, inputs) {
  if (!identical(inputs[["retain"]], "prediction") || !length(inputs[["fe"]])) return(result)
  conn = inputs[["conn"]]
  components = lapply(inputs[["prediction_fe_components"]], unique)
  from = inputs[["from_statement"]]
  source = paste0("(SELECT * ", from, ")")
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  design = sql_model_matrix(inputs[["fml"]], conn, sub("^FROM\\s+", "", from, ignore.case = TRUE),
    expand = "all", fe_vars = inputs[["fe"]], intercept = inputs[["intercept"]], factor_refs = inputs[["factor_refs"]])
  register_factor_display(inputs, design[["factor_levels"]], design[["level_aliases"]], design[["col_names"]], design[["aliases"]], factor_classes = design[["factor_classes"]])
  map = inputs[["name_map"]]
  variable = names(map) %in% names(inputs[["source_map"]])
  map[variable] = coef_var_name(map[variable])
  design_names = standardize_coef_names(restore_names(design[["aliases"]], map))
  beta = coef.dbreg(result)
  if (length(setdiff(names(beta), design_names))) stop("Could not align fitted OLS prediction slopes.", call. = FALSE)
  full_beta = stats::setNames(rep(0, length(design_names)), design_names)
  full_beta[names(beta)] = beta
  taken = unique(c(inputs[["yvar"]], inputs[["xvars"]], inputs[["fe"]], unlist(components),
    inputs[["weights"]], inputs[["cluster_var"]], inputs[["row_id_safe"]]))
  fe = safe_identifiers(paste0("__dbglm_fe", seq_along(components)), taken = taken, fold_case = folds_identifier_case(conn))
  attr(fe, "slopes") = rep(list(NULL), length(fe))
  ranked = paste0("(", sql_rank_fixed_effects(conn, source, components, fe), ")")
  x = if (length(design_names)) paste0("__dbglm_prediction_x", seq_along(design_names)) else character()
  weights = if (is.null(inputs[["weights"]])) "1.0" else q(inputs[["weights"]])
  summary_terms = c("COUNT(*) AS n", paste0("MAX(", weights, ") AS weight_max"),
    paste0("AVG(CAST(", q(inputs[["yvar"]]), " AS DOUBLE)) AS y_center"),
    sprintf("AVG(CAST((%s) AS DOUBLE)) AS center%d", design[["select_exprs"]], seq_along(x)),
    sprintf("COUNT(DISTINCT %s) AS size%d", fe, seq_along(fe)))
  row = dbglm_get_query(conn, paste("SELECT", paste(summary_terms, collapse = ","), "FROM", ranked))
  n = as.numeric(row[["n"]])
  if (n != result[["nobs_orig"]]) stop("The OLS sample changed before prediction recovery.", call. = FALSE)
  center = vapply(seq_along(x), function(i) as.numeric(row[[paste0("center", i)]]), numeric(1))
  y_center = as.numeric(row[["y_center"]])
  if (any(!is.finite(c(center, y_center)))) stop("OLS prediction centering is nonfinite.", call. = FALSE)
  sizes = vapply(seq_along(fe), function(i) as.numeric(row[[paste0("size", i)]]), numeric(1))
  weight_max = as.numeric(row[["weight_max"]])
  weight_scale = if (weight_max < 1) 2^floor(log2(weight_max)) else 1
  keys = prediction_key_fields(conn, components)
  fields = c(fe, keys[["fields"]],
    paste0("(CAST(", q(inputs[["yvar"]]), " AS DOUBLE)-(", dbglm_number(y_center), ")) AS eta"),
    "CAST(0.0 AS DOUBLE) AS __dbglm_o",
    paste0("CAST(", weights, " AS DOUBLE)/(", dbglm_number(weight_scale), ") AS __dbglm_a"),
    "CAST(1.0 AS DOUBLE) AS w",
    sprintf("(CAST((%s) AS DOUBLE)-(%s)) AS %s", design[["select_exprs"]], dbglm_number(center), x))
  frame = paste0("(", dbglm_sql(conn, paste("SELECT", paste(fields, collapse = ","), "FROM", ranked)), ")")
  manager = dbglm_table_manager(conn, inputs[["sample_resources"]])
  public_components = lapply(components, restore_names, map = inputs[["source_map"]])
  result[["fitted_state"]] = build_prediction_state(conn, frame, frame, fe, public_components,
    x, design_names, center, character(), sizes, manager[["create"]], manager[["drop"]], n,
    dictionary_components = keys[["columns"]], beta = unname(full_beta),
    expected_rss = attr(result[["vcov"]], "rss") / weight_scale, constant = y_center,
    recovery_control = list(tol = inputs[["fe_tol"]], max_iter = inputs[["fe_max_iter"]]))
  levels = inputs[["factor_levels"]]
  if (length(levels)) names(levels) = restore_names(names(levels), inputs[["source_map"]])
  result[["factor_levels"]] = levels
  result
}

predict_linear_state = function(object, newdata, interval) {
  if (is.null(newdata)) stop("newdata is required; prediction retention does not keep input rows.", call. = FALSE)
  if (interval != "none") stop("Intervals with learned effects are not implemented; use interval='none'.", call. = FALSE)
  data = as.data.frame(newdata)
  rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
  state = object[["fitted_state"]]
  database = identical(state[["prediction"]][["storage"]], "database")
  if (database) fitted_resource_connection(object[["resource_state"]])
  required = unique(c(all.vars(rhs), unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns"))))
  absent = setdiff(required, names(data))
  if (length(absent)) stop("newdata is missing required columns: ", paste(absent, collapse = ", "), call. = FALSE)
  valid = rep(TRUE, nrow(data))
  for (name in required) {
    valid = valid & !is.na(data[[name]])
    if (is.numeric(data[[name]]) && !inherits(data[[name]], "integer64")) valid = valid & is.finite(data[[name]])
  }
  prediction = rep(NA_real_, nrow(data))
  if (!any(valid)) return(prediction)
  rows = data[valid, , drop = FALSE]
  positions = which(valid)
  lookup = if (database) prediction_database_lookup(state, object[["resource_state"]], rows) else NULL
  groups = if (!database) prediction_match_groups(state, rows) else NULL
  known = if (database) !is.na(lookup[["gamma"]]) else if (length(groups))
    Reduce(`&`, lapply(groups, function(group) !is.na(group))) else rep(TRUE, nrow(rows))
  if (!any(known)) return(prediction)
  if (!all(known)) {
    rows = rows[known, , drop = FALSE]
    groups = lapply(groups, function(group) group[known])
    if (database) lookup = lookup[known, , drop = FALSE]
    positions = positions[known]
  }
  levels = object[["factor_levels"]]
  levels = levels[intersect(names(levels), all.vars(rhs))]
  coded = prediction_factors(list(factor_levels = levels), rows)
  design = as.matrix(fitted_model_matrix(rhs, coded))
  absent = setdiff(state[["specification"]][["design_names"]], colnames(design))
  if (length(absent)) stop("Could not reconstruct fitted OLS design: ", paste(absent, collapse = ", "), call. = FALSE)
  prediction[positions] = predict_from_fitted_state(state, design, rows, groups = groups, lookup = lookup)
  prediction
}
