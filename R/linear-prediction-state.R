# Use the original linear fit's slopes and RSS; only nuisance coefficients are
# recovered. Both model classes share the bounded FE kernel and prediction state.
finalize_linear_prediction = function(result, inputs) {
  if (!retain_prediction(inputs[["retain"]]) || !length(inputs[["fe"]])) return(result)
  reused = tryCatch(finalize_fwl_prediction(result, inputs),
    dbreg_prediction_overflow = function(error) NULL)
  if (!is.null(reused)) return(reused)
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
  trend_sources = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  trend_columns = if (length(trend_sources)) safe_identifiers(
    paste0("__dbglm_prediction_t", seq_along(trend_sources)), taken = c(taken, fe),
    fold_case = folds_identifier_case(conn)) else character()
  trend_aliases = stats::setNames(trend_columns, restore_names(trend_sources, inputs[["source_map"]]))
  attr(fe, "slopes") = lapply(inputs[["fe_slopes"]], function(vars) trend_columns[match(vars, trend_sources)])
  ranked = paste0("(", sql_rank_fixed_effects(conn, source, components, fe), ")")
  x = if (length(design_names)) paste0("__dbglm_prediction_x", seq_along(design_names)) else character()
  weights = if (is.null(inputs[["weights"]])) "1.0" else q(inputs[["weights"]])
  summary_terms = c("COUNT(*) AS n", paste0("MAX(", weights, ") AS weight_max"),
    paste0("AVG(CAST(", q(inputs[["yvar"]]), " AS DOUBLE)) AS y_center"),
    sprintf("AVG(CAST((%s) AS DOUBLE)) AS center%d", design[["select_exprs"]], seq_along(x)),
    sprintf("COUNT(DISTINCT %s) AS size%d", fe, seq_along(fe)),
    if (length(trend_sources)) sprintf("MIN(%1$s)/2.0+MAX(%1$s)/2.0 AS trend_center%2$d",
      sql_cast_double(q(trend_sources), conn), seq_along(trend_sources)))
  row = backend_get_query(conn, paste("SELECT", paste(summary_terms, collapse = ","), "FROM", ranked))
  n = as.numeric(row[["n"]])
  if (n != result[["nobs_orig"]]) stop("The OLS sample changed before prediction recovery.", call. = FALSE)
  center = vapply(seq_along(x), function(i) as.numeric(row[[paste0("center", i)]]), numeric(1))
  y_center = as.numeric(row[["y_center"]])
  if (any(!is.finite(c(center, y_center)))) stop("OLS prediction centering is nonfinite.", call. = FALSE)
  sizes = vapply(seq_along(fe), function(i) as.numeric(row[[paste0("size", i)]]), numeric(1))
  weight_max = as.numeric(row[["weight_max"]])
  weight_scale = if (weight_max < 1) 2^floor(log2(weight_max)) else 1
  trend_centers = vapply(seq_along(trend_sources), function(i) as.numeric(row[[paste0("trend_center", i)]]), numeric(1))
  attr(trend_aliases, "centers") = stats::setNames(trend_centers, names(trend_aliases))
  keys = prediction_key_fields(conn, components)
  fields = c(fe, keys[["fields"]],
    paste0("(CAST(", q(inputs[["yvar"]]), " AS DOUBLE)-(", dbglm_number(y_center), ")) AS eta"),
    "CAST(0.0 AS DOUBLE) AS __dbglm_o",
    paste0("CAST(", weights, " AS DOUBLE)/(", dbglm_number(weight_scale), ") AS __dbglm_a"),
    "CAST(1.0 AS DOUBLE) AS w",
    if (length(trend_sources)) paste0("(", sql_cast_double(q(trend_sources), conn), "-(",
      dbglm_number(trend_centers), ")) AS ", trend_columns),
    sprintf("(CAST((%s) AS DOUBLE)-(%s)) AS %s", design[["select_exprs"]], dbglm_number(center), x))
  frame = paste0("(", backend_sql(conn, paste("SELECT", paste(fields, collapse = ","), "FROM", ranked)), ")")
  manager = fit_table_manager(conn, inputs[["sample_resources"]])
  public_components = lapply(components, restore_names, map = inputs[["source_map"]])
  result[["fitted_state"]] = build_prediction_state(conn, frame, frame, fe, public_components,
    x, design_names, center, trend_aliases, sizes, manager[["create"]], manager[["drop"]], n,
    dictionary_components = keys[["columns"]], beta = unname(full_beta),
    expected_rss = attr(result[["vcov"]], "rss") / weight_scale, constant = y_center,
    recovery_control = list(tol = inputs[["fe_tol"]], max_iter = inputs[["fe_max_iter"]]))
  levels = inputs[["factor_levels"]]
  if (length(levels)) names(levels) = restore_term_names(names(levels), inputs[["source_map"]])
  result[["factor_levels"]] = levels
  result
}

# Native typed IDs already arrived in compact FE summaries. Reuse them as
# prediction dictionaries, without scanning rows to rank the same keys again.
fwl_prediction_dictionary = function(ids, sql_type, column) {
  type = prediction_key_type(sql_type)
  if (is.na(type)) return(NULL)
  values = switch(type, numeric = as.numeric(ids), float = as.numeric(ids),
    logical = ifelse(toupper(ids) == "TRUE", TRUE, ifelse(toupper(ids) == "FALSE", FALSE, NA)), ids)
  keys = prediction_tuple_key(stats::setNames(list(values), column), type)
  if (anyNA(keys) || anyDuplicated(keys)) return(NULL)
  list(columns = column, types = type, key = keys, id = seq_along(ids), sql_types = sql_type)
}

finalize_fwl_prediction = function(result, inputs) {
  moments = inputs[["linear_prediction_moments"]]
  if (is.null(moments) || moments[["nobs"]] != result[["nobs_orig"]]) return(NULL)
  x = moments[["x"]]
  p = length(x)
  map = inputs[["name_map"]]
  variables = names(map) %in% names(inputs[["source_map"]])
  map[variables] = coef_var_name(map[variables])
  design_names = standardize_coef_names(restore_names(x, map))
  fitted = coef.dbreg(result)
  if (length(setdiff(names(fitted), design_names))) return(NULL)
  beta = stats::setNames(rep(0, p), design_names)
  beta[names(fitted)] = fitted
  alpha = moments[["alpha"]]
  gamma = as.numeric(alpha[, 1L] - alpha[, seq_len(p) + 1L, drop = FALSE] %*% unname(beta))
  starts = cumsum(c(0, moments[["sizes"]]))
  nuisance_bound = sum(vapply(seq_along(moments[["sizes"]]), function(i)
    max(abs(gamma[starts[[i]] + seq_len(moments[["sizes"]][[i]])])), numeric(1)))
  x_bound = sqrt(pmax(diag(moments[["cross"]])[-1L], 0) / moments[["min_weight"]])
  prediction_bound = nuisance_bound + sum(abs(beta) * x_bound)
  response_scale = max(1, sqrt(max(moments[["cross"]][1L, 1L], 0) / moments[["sum_w"]]))
  # Keep the existing row certificate for numerically risky cancellation.
  if (!is.finite(prediction_bound) || prediction_bound > 1e8 * response_scale) return(NULL)
  components = lapply(inputs[["prediction_fe_components"]], unique)
  public = lapply(components, restore_names, map = inputs[["source_map"]])
  conn = inputs[["conn"]]
  source = paste0("(SELECT * ", inputs[["from_statement"]], ")")
  dictionaries = vector("list", length(components))
  info = moments[["info"]]
  for (i in seq_along(components)) {
    ids = info[[i]][["id"]]
    dictionary = if (length(components[[i]]) == 1L) {
      fwl_prediction_dictionary(ids, moments[["types"]][[i]], public[[i]])
    } else {
      key = paste0("CAST(", DBI::dbQuoteIdentifier(conn, inputs[["fe"]][[i]]), " AS VARCHAR)")
      dictionary = prediction_fe_dictionary(conn, source, key, components[[i]],
        floor(prediction_state_limit() / length(components)), public_components = public[[i]])
      dictionary[["id"]] = match(as.character(dictionary[["id"]]), ids)
      if (length(dictionary[["id"]]) != length(ids) || anyNA(dictionary[["id"]])) return(NULL)
      dictionary
    }
    if (is.null(dictionary)) return(NULL)
    dictionaries[[i]] = dictionary
    info[[i]][["id"]] = seq_along(ids)
  }
  order = c(seq_len(p) + 1L, 1L)
  precomputed = list(info = info, sys = moments[["sys"]],
    B = moments[["B"]][, order, drop = FALSE], alpha = alpha[, order, drop = FALSE],
    cross = moments[["cross"]][order, order, drop = FALSE],
    within = moments[["within"]][order, order, drop = FALSE], dictionaries = dictionaries)
  state = tryCatch(build_portable_prediction_state(conn, NULL, NULL, inputs[["fe"]], public,
    x, design_names, moments[["centers"]][-1L], character(), moments[["sizes"]], NULL, NULL,
    moments[["nobs"]], beta = unname(beta),
    expected_rss = attr(result[["vcov"]], "rss") / moments[["weight_scale"]],
    constant = moments[["centers"]][[1L]], precomputed = precomputed),
    dbreg_prediction_overflow = function(error) NULL)
  if (is.null(state)) return(NULL)
  result[["fitted_state"]] = state
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
  required = unique(c(all.vars(rhs), object[["offset"]], unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns")),
    unlist(lapply(state[["specification"]][["operations"]], `[[`, "variable"))))
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
  levels = levels[intersect(names(levels), c(all.vars(rhs), names(explicit_factor_sources(attr(stats::terms(rhs), "term.labels")))))]
  design = as.matrix(build_newdata_design(rhs, rows, levels, object[["i_columns"]], object[["i_specs"]]))
  absent = setdiff(state[["specification"]][["design_names"]], colnames(design))
  if (length(absent)) stop("Could not reconstruct fitted OLS design: ", paste(absent, collapse = ", "), call. = FALSE)
  prediction[positions] = predict_from_fitted_state(state, design, rows, groups = groups, lookup = lookup) + linear_prediction_offset(object, rows)
  prediction
}
