# Strategy execution functions for dbreg
# Each function takes an inputs environment and returns a result list.

#' Execute moments strategy (no fixed effects)
#' @keywords internal
execute_moments_strategy = function(inputs) {
  env2env(
    inputs,
    environment(), c(
      "from_statement", "conn", "sql_only", "verbose",
      "data_only", "vcov_type_req", "fml", "fe", "intercept", "xvars", "yvar",
      "weights", "compression_ratio_est"
    )
  )
  # Get SQL expressions for design matrix terms
  # For interactions/factors, this expands to CASE WHEN expressions
  # Expand factors and interactions to SQL expressions (numeric terms pass
  # through unchanged)
  table_ref = sub("^FROM\\s+", "", from_statement, ignore.case = TRUE)
  sql_design = sql_model_matrix(
    fml,
    conn,
    table_ref,
    expand = "all",
    fe_vars = fe,
    reserved_names = sql_design_reserved(inputs),
    factor_refs = inputs[["factor_refs"]],
    intercept = intercept
  )
  xvars_sql = sql_design[["select_exprs"]]
  xvar_names = sql_design[["aliases"]]
  register_sql_design(inputs, sql_design)
  
  frame = build_base_frame(inputs, xvars_sql, xvar_names)
  moments_sql = frame_moments_sql(frame, yvar, xvar_names, conn)

  if (sql_only) {
    return(moments_sql)
  }
  if (verbose) {
    message(if (!is.null(weights)) "[dbreg] Executing weighted moments SQL\n" else "[dbreg] Executing moments SQL\n")
  }
  moments_df = dbGetQuery(conn, moments_sql)
  if (data_only) {
    return(moments_df)
  }
  n_total = as.numeric(moments_df[["n_total"]])

  S = cross_moments(moments_df, xvar_names, yvar, intercept = frame[["intercept"]])
  vars_all = c(if (frame[["intercept"]]) "(Intercept)", xvar_names)
  p = length(vars_all)
  fit = fit_ols(S, yvar, vars_all, drop_collinear = TRUE, verbose = verbose)
  betahat = fit[["coef"]]
  p = length(fit[["x"]])
  rss = fit[["rss"]]
  tss = moment_tss(S, yvar)
  df_res = n_total - p
  
  # Compute meat matrix if needed (HC1 or cluster)
  meat = fit_meat(fit, frame, yvar, inputs)
  
  vcov_mat = compute_vcov(
    ssc = inputs[["ssc"]],
    vcov_type = vcov_type_req,
    strategy = "moments",
    XtX_inv = fit[["bread"]],
    sandwich_bread = fit[["sandwich_bread"]],
    rss = rss,
    df_res = df_res,
    nobs_orig = n_total,
    n_params = p,
    meat = meat
  )
  attr(vcov_mat, "rss") = rss
  attr(vcov_mat, "tss") = tss

  coeftable = gen_coeftable(betahat, vcov_mat, df_res)
  layout = linear_fit_layout(inputs[["design_columns"]], fit[["x"]])

  list(
    coeftable = coeftable,
    linear_design = layout,
    vcov = vcov_mat,
    fml = fml,
    yvar = yvar,
    xvars = standardize_coef_names(xvars),
    collin.var = standardize_coef_names(fit[["dropped"]]),
    fe = NULL,
    weights = weights,
    query_string = moments_sql,
    nobs = 1L,
    nobs_orig = n_total,
    strategy = "moments",
    prediction_coordinates = fit_prediction_state(fit, yvar, vcov_mat, rss, df_res, n_total, meat),
    compression_ratio_est = compression_ratio_est,
    df_residual = df_res
  )
}

#' Execute demean strategy (1+ fixed effects)
#' 
#' Double demeaning / within estimator. Gives identical coefficients to 
#' fixed effects regression.
#' 
#' @keywords internal
execute_demean_strategy = function(inputs) {
  env2env(
    inputs,
    environment(),
    c(
      "from_statement", "conn", "sql_only", "data_only",
      "vcov_type_req", "fml", "fe", "xvars", "yvar", "verbose", "weights",
      "compression_ratio_est", "ssc"
    )
  )
  # Expand factors and interactions to SQL expressions (numeric terms pass
  # through unchanged)
  table_ref = sub("^FROM\\s+", "", from_statement, ignore.case = TRUE)
  sql_design = sql_model_matrix(
    fml,
    conn,
    table_ref,
    expand = "all",
    fe_vars = fe,
    reserved_names = sql_design_reserved(inputs),
    factor_refs = inputs[["factor_refs"]]
  )
  xvars_sql = sql_design[["select_exprs"]]
  xvar_names = sql_design[["aliases"]]
  register_sql_design(inputs, sql_design)
  
  frame = build_demean_frame(inputs, xvars_sql, xvar_names)
  if (inherits(frame, "dbreg_query_plan")) return(frame)
  on.exit(frame[["cleanup"]](), add = TRUE)
  demean_sql = frame_moments_sql(frame, yvar, xvar_names, conn)
  if (sql_only) {
    return(demean_sql)
  }

  # Execute SQL and build matrices
  if (verbose) {
    message(if (!is.null(weights)) "[dbreg] Executing weighted demean SQL\n" else "[dbreg] Executing demean SQL\n")
  }
  demean_df = if (!is.null(frame[["moments"]])) frame[["moments"]] else dbGetQuery(conn, demean_sql)
  if (data_only) {
    return(demean_df)
  }
  n_total = as.numeric(demean_df[["n_total"]])
  n_fe_levels = vapply(seq_along(fe), function(k) {
    val = demean_df[[sprintf("n_fe%d", k)]]
    if (is.null(val)) 1 else as.numeric(val)
  }, numeric(1))

  S = cross_moments(demean_df, xvar_names, yvar, intercept = FALSE)
  fit = fit_ols(S, yvar, xvar_names, drop_collinear = TRUE, verbose = verbose)
  betahat = fit[["coef"]]
  xvar_names_kept = fit[["x"]]
  collin_vars = fit[["dropped"]]
  p_kept = length(xvar_names_kept)
  rss = fit[["rss"]]
  n_fe = length(fe)
  df_fe = sum(n_fe_levels) - (n_fe - 1)
  df_res = n_total - p_kept - df_fe
  
  # Compute meat matrix if needed (HC1 or cluster)
  meat = fit_meat(fit, frame, yvar, inputs)
  n_params_cluster = linear_parameter_count(inputs, p_kept + df_fe, df_fe)

  vcov_mat = compute_vcov(
    ssc = inputs[["ssc"]],
    vcov_type = vcov_type_req,
    strategy = "demean",
    XtX_inv = fit[["bread"]],
    sandwich_bread = fit[["sandwich_bread"]],
    rss = rss,
    df_res = df_res,
    nobs_orig = n_total,
    n_params = n_params_cluster,
    meat = meat
  )
  attr(vcov_mat, "rss") = rss
  attr(vcov_mat, "tss") = source_tss(inputs)

  coeftable = gen_coeftable(betahat, vcov_mat, df_res)
  layout = linear_fit_layout(inputs[["design_columns"]], fit[["x"]])

  list(
    coeftable = coeftable,
    linear_design = layout,
    vcov = vcov_mat,
    fml = fml,
    yvar = yvar,
    xvars = standardize_coef_names(xvar_names_kept),
    collin.var = standardize_coef_names(collin_vars),
    fe = fe,
    weights = weights,
    query_string = demean_sql,
    nobs = 1L,
    nobs_orig = n_total,
    strategy = "demean",
    compression_ratio_est = compression_ratio_est,
    df_residual = df_res,
    n_fe_levels = n_fe_levels
  )
}

#' Execute true Mundlak/CRE strategy
#'
#' Regresses Y on X plus group means of X for each fixed effect.
#' Y is NOT demeaned - predictions are on the original scale.
#'
#' @keywords internal
execute_mundlak_strategy = function(inputs) {
  env2env(
    inputs,
    environment(), c(
      "from_statement", "conn", "cluster_var", "sql_only",
      "verbose", "data_only", "vcov_type_req", "yvar", "fe", "fml", "xvars",
      "weights", "compression_ratio_est"
    )
  )
  n_fe = length(fe)

  if (n_fe == 0) {
    stop("mundlak strategy requires at least one fixed effect")
  }

  # Expand factors and interactions to SQL expressions (numeric terms pass
  # through unchanged)
  table_ref = sub("^FROM\\s+", "", from_statement, ignore.case = TRUE)
  sql_design = sql_model_matrix(
    fml,
    conn,
    table_ref,
    expand = "all",
    fe_vars = fe,
    reserved_names = sql_design_reserved(inputs),
    factor_refs = inputs[["factor_refs"]]
  )
  xvars_sql = sql_design[["select_exprs"]]
  xvar_names = sql_design[["aliases"]]
  register_sql_design(inputs, sql_design)
  if (!length(xvar_names)) stop("Mundlak requires at least one expanded regressor.", call. = FALSE)
  
  source_roles = unique(c(fe, yvar, xvars, xvar_names, weights, cluster_var))
  private = moment_sql_names(c("base", "augmented", "dbreg_mundlak_location", paste0("fe", seq_along(fe), "_means")),
    paste(from_statement, paste(source_roles, collapse = " ")), conn)
  base_name = private[[1L]]
  augmented_name = private[[2L]]
  location_name = private[[3L]]
  mean_tables = private[-seq_len(3L)]
  reserved = c(source_roles, private)
  weight_alias = safe_identifiers("dbreg_mundlak_weight__", taken = reserved, fold_case = folds_identifier_case(conn))
  scale_alias = safe_identifiers("dbreg_mundlak_scale__", taken = c(reserved, weight_alias), fold_case = folds_identifier_case(conn))
  weights_expr_base = if (!is.null(weights)) paste0("(", sql_cast_double(weights, conn), ")/", scale_alias)
  weights_expr_aug = if (!is.null(weights)) sql_cast_double(weight_alias, conn)
  mean_aliases = safe_identifiers(unlist(lapply(seq_along(fe), function(k) sprintf("dbreg_bar%d_%d", k, seq_along(xvar_names)))),
    taken = c(reserved, weight_alias, scale_alias), fold_case = folds_identifier_case(conn))
  anchors = safe_identifiers(paste0("dbreg_mundlak_center_", seq_along(xvar_names)),
    taken = c(reserved, weight_alias, scale_alias, mean_aliases), fold_case = folds_identifier_case(conn))
  names(anchors) = xvar_names


  # Build base CTE with expanded columns AND original xvars (for group means)
  base_select = c(fe, yvar, xvars)
  for (i in seq_along(xvar_names)) {
    # Only add expanded terms that aren't already in original xvars
    if (!xvar_names[i] %in% xvars) {
      base_select = c(base_select, sprintf("%s AS %s", xvars_sql[i], xvar_names[i]))
    }
  }
  if (!is.null(weights) && !weights %in% c(fe, yvar, xvars, xvar_names)) {
    base_select = c(base_select, weights)
  }
  if (!is.null(cluster_var) && !cluster_var %in% c(fe, yvar, xvars, xvar_names, weights)) {
    base_select = c(base_select, cluster_var)
  }

  # Group means of every regressor column (incl. factor dummies and
  # interactions), one set per FE. Anything less isn't FE-equivalent (#77).
  # Means get opaque internal names, with display names registered below.
  cte_parts = character(0)
  join_parts = character(0)
  aug_bar_parts = character(0)
  xbar_all = character(0)
  for (k in seq_along(fe)) {
    fe_k = fe[k]
    xbar_k = mean_aliases[(k - 1L) * length(xvar_names) + seq_along(xvar_names)]
    xbar_all = c(xbar_all, xbar_k)
    means_cols = vapply(
      seq_along(xvar_names),
      function(i) sql_weighted_mean(paste0("(", sql_cast_double(xvar_names[i], conn), "-l.", anchors[i], ")"), weights_expr_base, xbar_k[i]),
      character(1)
    )
    cte_parts = c(cte_parts, if (length(means_cols)) {
      sprintf("%s AS (SELECT %s, %s FROM %s CROSS JOIN %s l GROUP BY %s)", mean_tables[k], fe_k, paste(means_cols, collapse = ", "), base_name, location_name, fe_k)
    } else {
      sprintf("%s AS (SELECT DISTINCT %s FROM %s)", mean_tables[k], fe_k, base_name)
    })
    join_parts = c(join_parts, sprintf("JOIN %s m%d ON b.%s = m%d.%s", mean_tables[k], k, fe_k, k, fe_k))
    aug_bar_parts = c(aug_bar_parts, paste0("m", k, ".", xbar_k))
  }

  # Select columns for augmented table (include FE for counting)
  aug_select_parts = c(sprintf("b.%s", fe), sprintf("b.%s", yvar), sprintf("b.%s", xvar_names))
  if (!is.null(cluster_var) && !cluster_var %in% c(fe, yvar, xvar_names)) {
    aug_select_parts = c(aug_select_parts, sprintf("b.%s AS %s", cluster_var, cluster_var))
  }
  if (!is.null(weights)) {
    aug_select_parts = c(aug_select_parts, sprintf("b.%s AS %s", weights, weight_alias))
  }
  aug_select = paste(c(aug_select_parts, aug_bar_parts, paste0("l.", anchors)), collapse = ", ")

  # Display names for the means: "<coef name>_bar_<fe>"
  map = inputs[["name_map"]]
  coef_map = map
  is_var = names(coef_map) %in% names(inputs[["source_map"]])
  coef_map[is_var] = coef_var_name(coef_map[is_var])
  term_display = standardize_coef_names(restore_names(xvar_names, coef_map))
  fe_display = restore_names(fe, map)
  mundlak_means = data.frame(
    term = rep(term_display, length(fe)),
    fe = rep(fe_display, each = length(xvar_names)),
    name = paste0(rep(term_display, length(fe)), "_bar_", rep(fe_display, each = length(xvar_names)))
  )
  inputs[["name_map"]] = c(map, setNames(mundlak_means[["name"]], xbar_all))

  mean_columns = inputs[["design_columns"]][rep(seq_len(nrow(inputs[["design_columns"]])), length(fe)), , drop = FALSE]
  mean_columns[["id"]] = paste0("mean", seq_along(xbar_all))
  mean_columns[["alias"]] = xbar_all
  mean_columns[["role"]] = "nuisance"
  mean_columns[["recipe"]] = Map(function(recipe, group) list(kind = "mean", source = recipe, variable = group),
    mean_columns[["recipe"]], rep(fe, each = length(xvar_names)))
  inputs[["design_columns"]] = rbind(inputs[["design_columns"]], mean_columns)

  # Regressors before means, so collinearity (e.g. a regressor that's constant
  # within groups) drops the mean and keeps the regressor, as in CRE
  all_regressors = c(xvar_names, xbar_all)

  frame = list(
    cte_sql = paste0(
      "WITH ", base_name, " AS (SELECT ", paste(unique(base_select), collapse = ", "), " ", from_statement, "),\n",
      location_name, " AS (SELECT ", paste(c(
        if (!is.null(weights)) paste0("MAX(", sql_cast_double(weights, conn), ") AS ", scale_alias),
        paste0("MIN(", sql_cast_double(xvar_names, conn), ")/2+MAX(", sql_cast_double(xvar_names, conn), ")/2 AS ", anchors)), collapse = ","), " FROM ", base_name, "),\n",
      paste(cte_parts, collapse = ",\n"), ",\n",
      augmented_name, " AS (SELECT ", aug_select, " FROM ", base_name, " b CROSS JOIN ", location_name, " l ", paste(join_parts, collapse = " "), ")"
    ),
    cte_name = augmented_name,
    cols = setNames(c(yvar, all_regressors), c(yvar, all_regressors)),
    origins = setNames(rep(unname(anchors), length(fe)), xbar_all),
    weights_expr = weights_expr_aug,
    intercept = TRUE,
    cast = TRUE,
    prefix_terms = c(
      sql_count(conn, "n_total"),
      if (n_fe >= 1) sql_count(conn, "n_fe1", fe[1], distinct = TRUE) else "1 AS n_fe1",
      if (n_fe >= 2) sql_count(conn, "n_fe2", fe[2], distinct = TRUE) else "1 AS n_fe2"
    )
  )
  mundlak_sql = frame_moments_sql(frame, yvar, all_regressors, conn)

  if (sql_only) {
    return(mundlak_sql)
  }

  if (verbose) {
    message(if (!is.null(weights)) "[dbreg] Executing weighted mundlak SQL\n" else "[dbreg] Executing mundlak SQL\n")
  }
  mundlak_df = dbGetQuery(conn, mundlak_sql)
  if (data_only) {
    return(mundlak_df)
  }

  n_total = as.numeric(mundlak_df[["n_total"]])
  n_fe1 = as.numeric(mundlak_df[["n_fe1"]])
  n_fe2 = if (is.null(mundlak_df[["n_fe2"]])) NULL else as.numeric(mundlak_df[["n_fe2"]])

  S = cross_moments(mundlak_df, all_regressors, yvar)
  # Prediction computes the same centered group means, before restoring their
  # raw coordinate origins for the public coefficients and covariance.
  mean_origins = attr(S, "origins")[xbar_all]
  mean_columns[["recipe"]] = Map(function(recipe, origin) {
    recipe[["center"]] = origin
    recipe
  }, mean_columns[["recipe"]], as.list(mean_origins))
  inputs[["design_columns"]][["recipe"]][match(xbar_all, inputs[["design_columns"]][["alias"]])] = mean_columns[["recipe"]]
  vars_all = c("(Intercept)", all_regressors)
  p = length(vars_all)
  fit = fit_ols(S, yvar, vars_all, drop_collinear = TRUE, verbose = verbose)
  betahat = fit[["coef"]]
  p = length(fit[["x"]])
  rss = fit[["rss"]]
  tss = moment_tss(S, yvar)
  df_res = n_total - p

  # Compute meat matrix if needed (HC1 or cluster)
  meat = fit_meat(fit, frame, yvar, inputs)

  vcov_mat = compute_vcov(
    ssc = inputs[["ssc"]],
    vcov_type = vcov_type_req,
    strategy = "mundlak",
    XtX_inv = fit[["bread"]],
    sandwich_bread = fit[["sandwich_bread"]],
    rss = rss,
    df_res = df_res,
    nobs_orig = n_total,
    n_params = p,
    meat = meat
  )
  attr(vcov_mat, "rss") = rss
  attr(vcov_mat, "tss") = tss

  coeftable = gen_coeftable(betahat, vcov_mat, df_res)
  layout = linear_fit_layout(inputs[["design_columns"]], fit[["x"]])

  list(
    coeftable = coeftable,
    linear_design = layout,
    vcov = vcov_mat,
    fml = fml,
    yvar = yvar,
    xvars = standardize_coef_names(xvar_names),
    collin.var = fit[["dropped"]],
    fe = fe,
    weights = weights,
    query_string = mundlak_sql,
    nobs = 1L,
    nobs_orig = n_total,
    strategy = "mundlak",
    compression_ratio_est = compression_ratio_est,
    df_residual = df_res,
    mundlak_means = mundlak_means,
    prediction_coordinates = fit_prediction_state(fit, yvar, vcov_mat, rss, df_res, n_total, meat),
    n_fe1 = n_fe1,
    n_fe2 = n_fe2
  )
}

#' Execute compress strategy (groupby compression)
#' @keywords internal
execute_compress_strategy = function(inputs) {
  env2env(
    inputs,
    environment(),
    c(
      "has_interactions", "conn", "sql_only", "data_only", "vcov_type_req",
      "fml", "fe", "xvars", "weights", "yvar", "verbose", "intercept",
      "compression_ratio_est", "ssc"
    )
  )
  from_statement = inputs[["from_statement"]]
  # catch for sampled (limited) queries
  if (grepl("LIMIT\\s+\\d+\\s*$", from_statement, ignore.case = TRUE)) {
    from_statement = glue("FROM (SELECT * {from_statement})")
  }

  component_vars = unique(unlist(inputs[["composite_fes"]], use.names = FALSE))
  component_types = if (length(component_vars)) category_key_metadata(conn, sub("^FROM\\s+", "", from_statement), component_vars)[["types"]] else character()
  inputs[["composite_integer_keys"]] = component_vars[component_types == "integer"]
  overlap = intersect(xvars, component_vars[component_types %in% c("integer", "numeric", "float")])
  moment_fields = c("n", "sum_w", "sum_wy", "sum_wy_sq", "sum_w2", "sum_w2y", "sum_w2y_sq", "mean_Y", "wts", "dbreg_cell")
  protected = component_vars[tolower(component_vars) %in% tolower(moment_fields)]
  cell_sources = if (length(protected)) setNames(protected,
    safe_identifiers(paste0("dbreg_joint_component_", seq_along(protected)), taken = c(sql_design_reserved(inputs), moment_fields), fold_case = folds_identifier_case(conn))) else character()
  inputs[["composite_cell_sources"]] = cell_sources
  component_columns = component_vars
  component_columns[match(unname(cell_sources), component_vars)] = names(cell_sources)
  inputs[["composite_integer_keys"]] = component_columns[component_types == "integer"]
  if (length(cell_sources)) has_interactions = TRUE
  if (length(overlap)) has_interactions = TRUE

  # Handle interactions: expand to SQL expressions
  if (isTRUE(has_interactions)) {
    # Extract table name from FROM statement for sql_model_matrix
    table_ref = sub("^FROM\\s+", "", from_statement, ignore.case = TRUE)
    
    # Get SQL expansions for RHS terms (expand interactions only, keep main effects as-is)
    sql_design = sql_model_matrix(
      fml,
      conn,
      table_ref,
      expand = "interactions",
      fe_vars = fe,
      reserved_names = sql_design_reserved(inputs),
      factor_refs = inputs[["factor_refs"]],
      intercept = intercept
    )
    
    # Build SELECT expressions with aliases
    select_exprs = if (length(sql_design[["aliases"]])) paste0(sql_design[["select_exprs"]], " AS ", sql_design[["aliases"]]) else character()
    xvars_sql = paste(select_exprs, collapse = ", ")
    xvar_names = sql_design[["aliases"]]
    register_sql_design(inputs, sql_design)
  } else {
    xvars_sql = paste(xvars, collapse = ", ")
    xvar_names = xvars
  }
  
  rhs_aliases = character()
  if (length(overlap)) {
    positions = match(overlap, xvar_names)
    positions = positions[!is.na(positions)]
    if (length(positions)) {
      original = xvar_names[positions]
      aliases = safe_identifiers(paste0("dbreg_joint_rhs_", seq_along(positions)), taken = c(sql_design_reserved(inputs), xvar_names), fold_case = folds_identifier_case(conn))
      select_exprs[positions] = paste(sql_cast_double(original, conn), "AS", aliases)
      xvar_names[positions] = aliases
      columns = inputs[["design_columns"]]
      columns[["alias"]][match(original, columns[["alias"]])] = aliases
      for (position in match(original, inputs[["design_columns"]][["alias"]])) columns[["recipe"]][[position]][["kind"]] = "numeric"
      inputs[["design_columns"]] = columns
      rhs_aliases = aliases
      inputs[["composite_rhs_aliases"]] = setNames(original, aliases)
      inputs[["name_map"]] = c(inputs[["name_map"]], setNames(restore_names(original, inputs[["source_map"]]), aliases))
    }
  }
  if (length(cell_sources)) select_exprs = c(select_exprs, paste(unname(cell_sources), "AS", names(cell_sources)))
  emitted_names = c(xvar_names, names(cell_sources))
  group_cols = unique(c(xvar_names, fe, component_columns))
  expanded_from = if (isTRUE(has_interactions)) paste0("FROM (SELECT ", paste(c(select_exprs,
    setdiff(unique(c(fe, yvar, weights, component_columns)), emitted_names)), collapse = ","), " ", from_statement, ") dbreg_compress_source") else from_statement
  packet = compressed_moment_query(conn, expanded_from, group_cols, yvar, weights, intercept)
  query_string = packet[["sql"]]

  if (sql_only) {
    return(query_string)
  }
  if (verbose) {
    message(if (!is.null(weights)) "[dbreg] Executing weighted compress strategy SQL\n" else "[dbreg] Executing compress strategy SQL\n")
  }
  validate_composite_key_transport(inputs)
  if (!data_only) inputs[["composite_dictionaries"]] = composite_prediction_dictionaries(inputs)
  compressed_dat = read_compressed_moments(conn, packet, integer_text = inputs[["composite_integer_keys"]])
  nobs_orig = sum(compressed_dat[["n"]])
  nobs_comp = nrow(compressed_dat)
  compression_ratio = nobs_comp / max(nobs_orig, 1)

  if (verbose && compression_ratio > 0.8) {
    warning(paste0(
      sprintf(
        "[dbreg] compression ineffective (%.1f%% of original rows). ",
        100 * compression_ratio
      ),
      sprintf("Consider strategy = '%s'.", if (length(fe)) "demean" else "moments")
    ))
  }

  if (length(fe)) {
    for (f in fe) {
      compressed_dat[[f]] = factor(compressed_dat[[f]])
    }
  }
  if (data_only) {
    compressed_dat[c("dbreg_cell", rhs_aliases)] = NULL
    return(restore_compressed_data(compressed_dat))
  }

  # Build design matrix
  # Use expanded column names if interactions were present
  design_vars = if (isTRUE(has_interactions)) xvar_names else xvars
  # Factor dummies (incl. FE) built in R are named after their (safe) variable
  # and level; register display names for them
  is_cat = vapply(c(design_vars, fe), function(v) {
    v %in% names(compressed_dat) && (is.factor(compressed_dat[[v]]) || is.character(compressed_dat[[v]]))
  }, logical(1))
  register_factor_display(inputs, lapply(setNames(nm = c(design_vars, fe)[is_cat]), function(v) unique(compressed_dat[[v]])))
  matrix_formula = if (length(c(design_vars, fe))) reformulate(c(design_vars, fe), intercept = intercept) else ~1
  categorical = unique(c(design_vars, fe))[vapply(compressed_dat[unique(c(design_vars, fe))], function(value)
    is.factor(value) || is.character(value) || is.logical(value), logical(1))]
  contrasts = lapply(compressed_dat[categorical], function(value) {
    levels = if (is.logical(value)) c("FALSE", "TRUE") else if (is.factor(value)) levels(value) else levels(factor(value))
    stats::contr.treatment(levels)
  })
  if (!length(contrasts)) contrasts = NULL
  matrix_data = linear_matrix_data(compressed_dat)
  X = sparse.model.matrix(matrix_formula, matrix_data, contrasts.arg = contrasts)
  columns = matrix_column_layout(X, matrix_formula, matrix_data,
    inputs[["design_columns"]], slope_terms = seq_along(design_vars))
  absorbed = composite_absorbed_columns(columns, inputs[["composite_fes"]])
  if (length(absorbed)) X[, absorbed] = 0
  generated = vapply(columns[["recipe"]], function(recipe)
    !recipe[["kind"]] %in% c("raw", "numeric") && !is.null(recipe[["variable"]]), logical(1))
  columns[["alias"]] = allocate_sql_aliases(columns[["alias"]], generated,
    sql_design_reserved(inputs), folds_identifier_case(conn))
  colnames(X) = columns[["alias"]]
  if (ncol(X) == 0) {
    stop("Design matrix has zero columns.")
  }
  numeric_terms = which(vapply(c(design_vars, fe), function(name)
    is.numeric(compressed_dat[[name]]) && !is.factor(compressed_dat[[name]]), logical(1)))
  S = compressed_cross_moments(X, compressed_dat, yvar,
    center_columns = which(attr(X, "assign") %in% numeric_terms))
  fit = fit_ols(S, yvar, colnames(X), drop_collinear = TRUE, verbose = verbose)
  layout = list(version = 1L, columns = columns, fitted_ids = columns[["id"]][match(fit[["x"]], columns[["alias"]])])
  compressed_design = list(formula = matrix_formula, data_names = names(compressed_dat), contrasts = contrasts,
    positions = match(fit[["x"]], columns[["alias"]]))
  collin_vars = fit[["dropped"]]
  X = X[, fit[["x"]], drop = FALSE]
  betahat = fit[["coef"]]
  XtX_inv = fit[["bread"]]
  weight_scale = fit[["weight_scale"]]
  score_X = center_compressed_design(X, c("(Intercept)" = 0, fit[["centers"]])[colnames(X)])
  yhat = as.numeric(score_X %*% fit[["moment_coef"]])

  sum_w = compressed_dat[["sum_w"]]
  sum_wy = compressed_dat[["sum_wy"]]
  sum_wy_sq = compressed_dat[["sum_wy_sq"]]
  rss_g = sum_wy_sq - 2 * yhat * sum_wy + sum_w * (yhat^2)
  rss_total = sum(rss_g)*weight_scale
  df_fe = fe_intercept_df(vapply(compressed_dat[fe], nlevels, integer(1)))
  n_params_model = ncol(X)
  if (length(fe)) {
    roles = columns[["role"]][match(fit[["x"]], columns[["alias"]])]
    n_params_model = sum(roles == "slope" & fit[["x"]] != "(Intercept)") + df_fe
  }
  df_res = nobs_orig - n_params_model

  # Calculate TSS for R2
  sum_wy_total = sum(compressed_dat[["sum_wy"]])
  sum_wy_sq_total = sum(compressed_dat[["sum_wy_sq"]])
  sum_w_total = sum(compressed_dat[["sum_w"]])
  tss = (sum_wy_sq_total - (sum_wy_total^2 / sum_w_total))*weight_scale
  
  # For clustered SEs, need to query cluster-by-cell stats
  meat = NULL
  n_params_cluster = linear_parameter_count(inputs, n_params_model, df_fe)
  if (vcov_type_req == "hc1" && !is.null(weights)) {
    sum_w2 = compressed_dat[["sum_w2"]]
    sum_w2y = compressed_dat[["sum_w2y"]]
    sum_w2y_sq = compressed_dat[["sum_w2y_sq"]]
    rss_g_w2 = sum_w2y_sq - 2 * yhat * sum_w2y + sum_w2 * (yhat^2)
    meat = crossprod(score_X, Diagonal(x = as.numeric(rss_g_w2)) %*% score_X)
  }
  if (vcov_type_req == "cluster") {
    # Cells are keyed on the expanded columns, so expose them (and only the
    # other columns the query needs) to it
    cells_from = if (isTRUE(has_interactions)) {
      cols = setdiff(unique(c(inputs[["cluster_var"]], yvar, weights, fe, component_columns)), emitted_names)
      paste0("FROM (SELECT ", paste(c(cols, select_exprs), collapse = ", "), " ", from_statement, ") AS dbreg_cells")
    } else {
      from_statement
    }
    # Many clusters but few cells: sum the scores by cluster in the database
    # instead of fetching one row per cluster and cell
    meat = if (use_cluster_sql(nrow(compressed_dat), ncol(X), inputs)) {
      compute_meat_cluster_compress_sql(inputs, cells_from, group_cols, compressed_dat[["dbreg_cell"]], score_X, list(yhat),
        y_center = attr(compressed_dat, "moment_center"), weight_scale = weight_scale)[[1]]
    } else {
      compute_meat_cluster_compress(inputs, cells_from, group_cols, compressed_dat, score_X, yhat,
        y_center = attr(compressed_dat, "moment_center"), weight_scale = weight_scale)
    }
  }
  
  if (vcov_type_req == "hc1" && is.null(meat)) {
    meat = crossprod(score_X, Diagonal(x = as.numeric(rss_g)) %*% score_X)
  }
  if (!is.null(meat)) attr(meat, "weight_scale") = weight_scale
  vcov_mat = compute_vcov(
    ssc = inputs[["ssc"]],
    vcov_type = vcov_type_req,
    strategy = "compress",
    XtX_inv = XtX_inv,
    sandwich_bread = fit[["sandwich_bread"]],
    rss = rss_total,
    df_res = df_res,
    nobs_orig = nobs_orig,
    n_params = n_params_cluster,
    X = X,
    rss_g = rss_g,
    meat = meat
  )
  attr(vcov_mat, "rss") = rss_total
  attr(vcov_mat, "tss") = tss

  coeftable = gen_coeftable(betahat, vcov_mat, df_res)

  # Coefficient names for xvars (excluding intercept and FE dummies)
  # Match coefficients that start with any design_var name
  all_coef_names = rownames(coeftable)
  all_coef_names = all_coef_names[all_coef_names != "(Intercept)"]
  design_pattern = paste0("^(", paste(standardize_coef_names(design_vars), collapse = "|"), ")")
  coef_names = all_coef_names[grepl(design_pattern, all_coef_names)]

  return(
    list(
      coeftable = coeftable,
      data = restore_compressed_data(compressed_dat),
      linear_design = layout,
      compressed_design = compressed_design,
      vcov = vcov_mat,
      fml = fml,
      yvar = yvar,
      xvars = standardize_coef_names(xvars),
      collin.var = standardize_coef_names(collin_vars),
      coef_names = coef_names,
      fe = fe,
      weights = weights,
      query_string = query_string,
      nobs = nobs_comp,
      nobs_orig = nobs_orig,
      strategy = "compress",
      prediction_coordinates = fit_prediction_state(fit, yvar, vcov_mat, rss_total, df_res, nobs_orig, meat),
      compression_ratio = compression_ratio,
      compression_ratio_est = compression_ratio_est,
      df_residual = df_res
    )
  )
}

#' Degrees of freedom of FE variables nested within cluster variable
#' 
#' For ssc = "nested", we exclude an FE's parameters from K if the FE is nested
#' within the cluster variable (i.e., each FE value belongs to exactly one
#' cluster). Matching fixest, each nested FE contributes G - 1 parameters (its
#' levels less one for the absorbed intercept), so the result is
#' sum(G_j - 1) over nested FEs.
#' 
#' @keywords internal
nested_fe_df = function(inputs, from_statement = inputs[["from_statement"]]) {
  env2env(inputs, environment(), c("conn", "fe", "cluster_var"))
  if (is.null(fe) || length(fe) == 0 || is.null(cluster_var)) {
    return(0L)
  }
  backend = detect_backend(conn)[["name"]]
  nested_df = 0L
  for (f in fe) {
    # Match K.exact=FALSE's declared-FE shortcut when clustering on an FE.
    # Otherwise establish nesting from functional dependencies on fitted rows.
    if (cluster_var %in% fe && f != cluster_var) next
    n_levels = tryCatch({
      if (identical(f, cluster_var)) {
        count_sql = glue("SELECT CAST(COUNT(DISTINCT {f}) AS BIGINT) AS n FROM (SELECT * {from_statement}) t")
        dbGetQuery(conn, count_sql)[["n"]]
      } else {
        # Compute group count and nesting violation in a single scan.
        check_sql = glue(
          "SELECT CAST(COUNT(*) AS BIGINT) AS n, ",
          "SUM(CASE WHEN c > 1 THEN 1 ELSE 0 END) AS bad ",
          "FROM (SELECT {f}, COUNT(DISTINCT {cluster_var}) AS c ",
          "FROM (SELECT * {from_statement}) t GROUP BY {f}) g"
        )
        res = dbGetQuery(conn, check_sql)
        if (nrow(res) == 1L && isTRUE(as.numeric(res[["bad"]]) == 0)) res[["n"]] else 0L
      }
    }, error = function(e) {
      # Fall back to keeping the FE in K (larger, conservative SEs), but say so
      warning(
        "[dbreg] Could not check whether fixed effect '", f, "' is nested in '",
        cluster_var, "' (", conditionMessage(e), "). Treating it as not nested ",
        "for the ssc = \"nested\" correction.",
        call. = FALSE
      )
      0L
    })
    nested_df = nested_df + max(as.numeric(n_levels) - 1, 0)
  }
  nested_df
}
