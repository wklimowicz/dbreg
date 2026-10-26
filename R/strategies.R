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
    intercept = intercept,
    factor_refs = inputs[["factor_refs"]]
  )
  xvars_sql = sql_design[["select_exprs"]]
  xvar_names = sql_design[["aliases"]]
  register_factor_display(inputs, sql_design[["factor_levels"]], sql_design[["level_aliases"]], sql_design[["col_names"]], sql_design[["aliases"]], i_specs = sql_design[["i_specs"]], i_columns = sql_design[["i_columns"]], factor_classes = sql_design[["factor_classes"]])

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

  list(
    linear_design = native_linear_layout(sql_design, fit[["x"]], frame[["intercept"]]),
    coeftable = coeftable,
    vcov = vcov_mat,
    fml = fml,
    yvar = yvar,
    xvars = standardize_coef_names(xvars),
    fe = NULL,
    weights = weights,
    query_string = moments_sql,
    nobs = 1L,
    nobs_orig = n_total,
    strategy = "moments",
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
    factor_refs = inputs[["factor_refs"]]
  )
  xvars_sql = sql_design[["select_exprs"]]
  xvar_names = sql_design[["aliases"]]
  register_factor_display(inputs, sql_design[["factor_levels"]], sql_design[["level_aliases"]], sql_design[["col_names"]], sql_design[["aliases"]], i_specs = sql_design[["i_specs"]], i_columns = sql_design[["i_columns"]], factor_classes = sql_design[["factor_classes"]])

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
  demean_df = dbGetQuery(conn, demean_sql)
  if (data_only) {
    return(demean_df)
  }
  n_total = as.numeric(demean_df[["n_total"]])
  n_fe_levels = vapply(seq_along(fe), function(k) {
    val = demean_df[[sprintf("n_fe%d", k)]]
    if (is.null(val)) 1L else as.integer(val)
  }, integer(1))

  S = cross_moments(demean_df, xvar_names, yvar, intercept = FALSE)
  fit = fit_ols(S, yvar, xvar_names, drop_collinear = TRUE, verbose = verbose)
  betahat = fit[["coef"]]
  xvar_names_kept = fit[["x"]]
  collin_vars = fit[["dropped"]]
  p_kept = length(xvar_names_kept)
  rss = fit[["rss"]]
  n_fe = length(fe)
  df_fe = linear_fe_df(inputs, n_fe_levels)
  df_res = n_total - p_kept - df_fe

  # Compute meat matrix if needed (HC1 or cluster)
  meat = fit_meat(fit, frame, yvar, inputs)
  n_params_cluster = p_kept + df_fe - if (ssc[["K.fixef"]] == "none") fe_intercept_df(n_fe_levels) else 0 # K for CR1 correction
  if (vcov_type_req == "cluster") {
    # For ssc = "nested", exclude nested FE levels from K
    if (ssc[["K.fixef"]] == "nonnested") {
      n_params_cluster = p_kept + df_fe - nested_fe_df(inputs)
    }
  }

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
  attr(vcov_mat, "tss") = moment_tss(S, yvar)

  coeftable = gen_coeftable(betahat, vcov_mat, df_res)

  list(
    coeftable = coeftable,
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
    n_fe_levels = n_fe_levels,
    ap_iterations = frame[["ap_res"]][["iterations"]],
    ap_max_abs = frame[["ap_res"]][["max_abs"]],
    ap_method = frame[["ap_res"]][["method"]]
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
    factor_refs = inputs[["factor_refs"]]
  )
  xvars_sql = sql_design[["select_exprs"]]
  xvar_names = sql_design[["aliases"]]
  register_factor_display(inputs, sql_design[["factor_levels"]], sql_design[["level_aliases"]], sql_design[["col_names"]], sql_design[["aliases"]], i_specs = sql_design[["i_specs"]], i_columns = sql_design[["i_columns"]], factor_classes = sql_design[["factor_classes"]])

  weights_expr_base = sql_weight_expr(weights)
  weights_expr_aug = if (is.null(weights)) NULL else sql_weight_expr("weights")


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
    xbar_k = if (length(xvar_names)) sprintf("dbreg_bar%d_%d", k, seq_along(xvar_names)) else character(0)
    xbar_all = c(xbar_all, xbar_k)
    means_cols = vapply(
      seq_along(xvar_names),
      function(i) sql_weighted_mean(xvar_names[i], weights_expr_base, xbar_k[i]),
      character(1)
    )
    cte_parts = c(cte_parts, if (length(means_cols)) {
      sprintf("fe%d_means AS (SELECT %s, %s FROM base GROUP BY %s)", k, fe_k, paste(means_cols, collapse = ", "), fe_k)
    } else {
      sprintf("fe%d_means AS (SELECT DISTINCT %s FROM base)", k, fe_k)
    })
    join_parts = c(join_parts, sprintf("JOIN fe%d_means m%d ON b.%s = m%d.%s", k, k, fe_k, k, fe_k))
    aug_bar_parts = c(aug_bar_parts, paste0("m", k, ".", xbar_k))
  }

  # Select columns for augmented table (include FE for counting)
  aug_select_parts = c(sprintf("b.%s", fe), sprintf("b.%s", yvar), sprintf("b.%s", xvar_names))
  if (!is.null(cluster_var) && !cluster_var %in% c(fe, yvar, xvar_names)) {
    aug_select_parts = c(aug_select_parts, sprintf("b.%s AS %s", cluster_var, cluster_var))
  }
  if (!is.null(weights)) {
    aug_select_parts = c(aug_select_parts, sprintf("b.%s AS weights", weights))
  }
  aug_select = paste(c(aug_select_parts, aug_bar_parts), collapse = ", ")

  # Display names for the means: "<coef name>_bar_<fe>"
  map = inputs[["name_map"]]
  coef_map = coefficient_name_map(map, inputs[["source_map"]])
  term_display = standardize_coef_names(restore_names(xvar_names, coef_map))
  fe_display = restore_names(fe, map)
  mundlak_means = data.frame(
    term = rep(term_display, length(fe)),
    fe = rep(fe_display, each = length(xvar_names)),
    name = paste0(rep(term_display, length(fe)), "_bar_", rep(fe_display, each = length(xvar_names)))
  )
  inputs[["name_map"]] = c(map, setNames(mundlak_means[["name"]], xbar_all))

  # Regressors before means, so collinearity (e.g. a regressor that's constant
  # within groups) drops the mean and keeps the regressor, as in CRE
  all_regressors = c(xvar_names, xbar_all)

  frame = list(
    cte_sql = paste0(
      "WITH base AS (SELECT ", paste(base_select, collapse = ", "), " ", from_statement, "),\n",
      paste(cte_parts, collapse = ",\n"), ",\n",
      "augmented AS (SELECT ", aug_select, " FROM base b ", paste(join_parts, collapse = " "), ")"
    ),
    cte_name = "augmented",
    cols = setNames(c(yvar, all_regressors), c(yvar, all_regressors)),
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
  n_fe1 = mundlak_df[["n_fe1"]]
  n_fe2 = mundlak_df[["n_fe2"]]

  S = cross_moments(mundlak_df, all_regressors, yvar)
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

  list(
    coeftable = coeftable,
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
  sql_design = NULL
  compressed_layout = NULL
  from_statement = inputs[["from_statement"]]
  composite_columns = unique(unlist(inputs[["composite_fes"]], use.names = FALSE))
  trends = as.character(unique(unlist(inputs[["fe_slopes"]], use.names = FALSE)))
  trend_aliases = if (length(trends)) safe_identifiers(paste0("dbreg_compress_trend", seq_along(trends), "__"),
    taken = c(
      xvars, fe, yvar, weights, inputs[["cluster_var"]], composite_columns,
      trends, inputs[["row_id_safe"]], inputs[["retained_key_safe"]]
    )
  ) else character()
  trend_sources = stats::setNames(trends, trend_aliases)
  nuisance = unlist(lapply(seq_along(fe), function(i) {
    slopes = inputs[["fe_slopes"]][[i]]
    if (length(slopes)) paste(fe[i], trend_aliases[match(slopes, trends)], sep = ":") else character()
  }), use.names = FALSE)
  if (length(trends)) inputs[["name_map"]] = c(inputs[["name_map"]], stats::setNames(
    coef_var_name(restore_names(trends, inputs[["source_map"]])), trend_aliases
  ))
  # catch for sampled (limited) queries
  if (grepl("LIMIT\\s+\\d+\\s*$", from_statement, ignore.case = TRUE)) {
    from_statement = glue("FROM (SELECT * {from_statement})")
  }

  # Handle interactions: expand to SQL expressions
  if (isTRUE(has_interactions) || length(inputs[["factor_refs"]])) {
    # Extract table name from FROM statement for sql_model_matrix
    table_ref = sub("^FROM\\s+", "", from_statement, ignore.case = TRUE)

    # Get SQL expansions for RHS terms (expand interactions only, keep main effects as-is)
    sql_design = sql_model_matrix(
      fml,
      conn,
      table_ref,
      expand = "interactions",
      fe_vars = fe,
      intercept = intercept,
      factor_refs = inputs[["factor_refs"]]
    )

    # Build SELECT expressions with aliases
    select_exprs = paste0(sql_design[["select_exprs"]], " AS ", sql_design[["aliases"]])
    xvars_sql = paste(select_exprs, collapse = ", ")
    xvar_names = sql_design[["aliases"]]
    register_factor_display(inputs, sql_design[["factor_levels"]], sql_design[["level_aliases"]], sql_design[["col_names"]], sql_design[["aliases"]], i_specs = sql_design[["i_specs"]], i_columns = sql_design[["i_columns"]], factor_classes = sql_design[["factor_classes"]])
    # Both the compression and clustered score passes must read the expanded
    # columns. Aliases from a SELECT list do not exist in the original source.
    projection = frame_select(yvar, sql_design[["select_exprs"]], xvar_names,
      id_cols = setdiff(fe, c(yvar, xvar_names)), extra = c(weights, inputs[["cluster_var"]], composite_columns, trends, inputs[["retained_key_safe"]])
    )
    from_statement = paste0(
      "FROM (SELECT ", paste(projection, collapse = ", "),
      " ", from_statement, ") dbreg_compress_design"
    )
    xvars_sql = paste(xvar_names, collapse = ", ")
  } else {
    xvars_sql = paste(xvars, collapse = ", ")
    xvar_names = xvars
  }

  # FE columns (no expansion needed - used for grouping)
  fe_sql = if (length(fe)) paste(fe, collapse = ", ") else NULL

  # Combined columns for SELECT and GROUP BY
  all_cols_sql = if (!is.null(fe_sql)) paste(xvars_sql, fe_sql, sep = ", ") else xvars_sql
  extra_components = setdiff(composite_columns, c(xvar_names, fe))
  if (length(extra_components)) all_cols_sql = paste(all_cols_sql, paste(extra_components, collapse = ", "), sep = ", ")
  group_cols = unique(c(xvar_names, fe, extra_components, trends))
  group_cols_sql = paste(group_cols, collapse = ", ")
  query_string = compressed_cell_moments_sql(conn, from_statement, yvar, group_cols, weights, intercept)
  if (is.null(query_string)) {
    cell = compression_moment_source(conn, from_statement, yvar, weights, intercept)
    weights_expr = cell[["weights"]]
    moment_terms = build_weighted_moment_terms(
      y_sql = cell[["y"]],
      weights_expr = weights_expr,
      prefix_terms = c("COUNT(*) AS n", cell[["prefix"]]),
      include_w_sq = TRUE
    )

    query_string = paste0(
      "WITH cte AS (\n    SELECT\n        ",
      paste(c(group_cols, moment_terms), collapse = ",\n        "),
      "\n    ",
      cell[["from"]],
      if (length(group_cols)) paste0("\n    GROUP BY ", group_cols_sql) else "",
      "\n    )\n    SELECT\n    *,\n    sum_wy / sum_w AS mean_Y,\n    sqrt(sum_w) AS wts\n    FROM cte"
    )
  }

  if (sql_only) {
    return(query_string)
  }
  if (verbose) {
    message(if (!is.null(weights)) "[dbreg] Executing weighted compress strategy SQL\n" else "[dbreg] Executing compress strategy SQL\n")
  }
  compressed_frame = read_compressed_frame(inputs, query_string, from_statement)
  compressed_dat = compressed_frame[["data"]]
  compressed_classes = vapply(compressed_dat, function(column) class(column)[[1L]], character(1))
  query_string = compressed_frame[["sql"]]
  nobs_orig = as.numeric(sum(compressed_dat[["n"]]))
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

  if (length(inputs[["composite_fes"]])) {
    inputs[["composite_levels"]] = lapply(names(inputs[["composite_fes"]]), function(f) {
      keys = tuple_key(compressed_dat, inputs[["composite_fes"]][[f]])
      setNames(as.character(compressed_dat[[f]][!duplicated(keys)]), keys[!duplicated(keys)])
    })
    names(inputs[["composite_levels"]]) = names(inputs[["composite_fes"]])
  }
  if (length(fe)) {
    for (f in fe) {
      compressed_dat[[f]] = if (f %in% compressed_frame[["integers"]]) exact_integer_factor(compressed_dat[[f]]) else
        exact_numeric_factor(compressed_dat[[f]])
    }
  }
  numeric_i = vapply(Filter(
    function(spec) identical(spec[["source_type"]], "numeric"),
    inputs[["i_specs"]]
  ), `[[`, character(1), "first")
  for (nm in intersect(names(inputs[["factor_levels"]]), names(compressed_dat))) {
    if (!nm %in% fe && restore_names(nm, inputs[["source_map"]]) %in% numeric_i) next
    previous = compressed_dat[[nm]]
    compressed_dat[[nm]] = factor(previous, levels = inputs[["factor_levels"]][[nm]])
    if (!is.null(attr(previous, "dbreg_key_values"))) {
      attr(compressed_dat[[nm]], "dbreg_key_values") = attr(previous, "dbreg_key_values")[match(levels(compressed_dat[[nm]]), levels(previous))]
      attr(compressed_dat[[nm]], "dbreg_key_type") = attr(previous, "dbreg_key_type")
    }
  }
  if (data_only) {
    return(restore_compressed_data(compressed_dat))
  }

  # Build design matrix
  # Use expanded column names if interactions were present
  design_vars = if (isTRUE(has_interactions)) xvar_names else xvars
  # Factor dummies (incl. FE) built in R are named after their (safe) variable
  # and level; register display names for them
  is_cat = vapply(c(design_vars, fe), function(v) {
    v %in% names(compressed_dat) && (is.factor(compressed_dat[[v]]) || is.character(compressed_dat[[v]]) || is.logical(compressed_dat[[v]]))
  }, logical(1))
  categorical = unique(c(design_vars, fe)[is_cat])
  factor_levels = lapply(setNames(nm = categorical), function(v) {
    saved = inputs[["factor_levels"]][[v]]
    if (!is.null(saved)) return(saved)
    if (is.logical(compressed_dat[[v]])) c(FALSE, TRUE) else levels(factor(compressed_dat[[v]]))
  })
  register_factor_display(inputs, factor_levels, factor_classes = compressed_classes[categorical])
  matrix_formula = reformulate(if (length(c(design_vars, fe, nuisance))) c(design_vars, fe, nuisance) else "1", intercept = intercept)
  matrix_data = compressed_trend_data(compressed_dat, trend_sources)
  # Group intercepts absorb any common trend location. Remove it before
  # forming the nuisance Gram matrix, keeping large calendar origins stable.
  trend_centers = vapply(matrix_data[trend_aliases], function(value) min(value) / 2 + max(value) / 2, numeric(1))
  for (alias in trend_aliases) matrix_data[[alias]] = matrix_data[[alias]] - trend_centers[[alias]]
  X = fitted_model_matrix(matrix_formula, matrix_data)
  if (!length(fe)) {
    compressed_layout = compressed_linear_layout(X, matrix_formula, matrix_data,
      sql_design, fml, inputs)
    colnames(X) = compressed_layout[["aliases"]]
  }
  numeric_sql = if (length(trends)) stats::setNames(paste0(
    "(", sql_cast_double(trends, conn), "-(",
    dbglm_number(trend_centers), "))"
  ), trend_aliases) else character()
  score_expressions = if (vcov_type_req == "cluster") compressed_score_expressions(
    conn, matrix_formula,
    matrix_data, X, c(compressed_classes, stats::setNames(rep("numeric", length(trends)), trend_aliases)), numeric_sql
  )
  slope_terms = match(design_vars, attr(terms(matrix_formula), "term.labels"))
  slope_names = standardize_coef_names(colnames(X)[attr(X, "assign") %in% slope_terms])
  if (ncol(X) == 0) {
    stop("Design matrix has zero columns.")
  }
  Y = compressed_dat[, "mean_Y"]
  wts = compressed_dat[["wts"]]
  # Center numeric terms without expanding sparse factor indicators to N*G.
  numeric_terms = which(vapply(c(design_vars, fe), function(name)
    is.numeric(compressed_dat[[name]]) && !is.factor(compressed_dat[[name]]), logical(1)))
  S = compressed_cross_moments(X, compressed_dat, yvar,
    center_columns = which(attr(X, "assign") %in% numeric_terms)
  )
  ordinary = colnames(X)[attr(X, "assign") %in% slope_terms]
  absorbed = compressed_slope_rank(X, matrix_data, inputs, ordinary, trend_sources)
  fit = fit_ols(S, yvar, setdiff(colnames(X), absorbed), drop_collinear = TRUE, verbose = verbose)
  fit[["dropped"]] = c(absorbed, fit[["dropped"]])
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
  rss_total = sum(rss_g) * weight_scale
  # Use the same declared FE-count policy as the within estimator. A dummy
  # matrix can have extra rank deficiencies on disconnected/nested graphs;
  # its numerical rank must not silently select a different SSC convention.
  n_params_full = ncol(X)
  if (length(fe)) {
    n_slopes = length(intersect(slope_names, standardize_coef_names(colnames(X))))
    n_params_full = n_slopes + linear_fe_df(inputs, vapply(compressed_dat[fe], nlevels, integer(1)))
  }
  df_res = nobs_orig - n_params_full

  # Calculate TSS for R2
  sum_wy_total = sum(compressed_dat[["sum_wy"]])
  sum_wy_sq_total = sum(compressed_dat[["sum_wy_sq"]])
  sum_w_total = sum(compressed_dat[["sum_w"]])
  tss = (sum_wy_sq_total - (sum_wy_total^2 / sum_w_total)) * weight_scale

  # Cluster cross-products stay in SQL, even for observation-sized groups.
  meat = NULL
  n_params_cluster = n_params_full - if (ssc[["K.fixef"]] == "none" && length(fe))
    fe_intercept_df(vapply(compressed_dat[fe], nlevels, integer(1))) else 0 # K for CR1 correction
  if (vcov_type_req == "hc1" && !is.null(weights)) {
    sum_w2 = compressed_dat[["sum_w2"]]
    sum_w2y = compressed_dat[["sum_w2y"]]
    sum_w2y_sq = compressed_dat[["sum_w2y_sq"]]
    rss_g_w2 = sum_w2y_sq - 2 * yhat * sum_w2y + sum_w2 * (yhat^2)
    meat = crossprod(score_X, Diagonal(x = as.numeric(rss_g_w2)) %*% score_X)
  }
  if (vcov_type_req == "cluster") {
    score_frame = list(
      cte_sql = paste0("WITH dbreg_compress_scores AS (SELECT * ", from_statement, ")"),
      cte_name = "dbreg_compress_scores", weights_expr = sql_weight_expr(weights),
      cols = c(setNames(sql_cast_double(yvar, conn), yvar), score_expressions)
    )
    meat = fit_meats(list(fit), score_frame, yvar, inputs)[[1L]]
    # For ssc = "nested", exclude nested FE levels from K
    if (ssc[["K.fixef"]] == "nonnested") {
      n_params_cluster = n_params_full - nested_fe_df(inputs, from_statement)
    }
  }

  if (!is.null(meat)) attr(meat, "weight_scale") = weight_scale
  if (vcov_type_req == "hc1" && is.null(meat)) {
    meat = crossprod(score_X, Diagonal(x = as.numeric(rss_g)) %*% score_X)
    attr(meat, "weight_scale") = weight_scale
  }
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

  if (retain_sample(inputs[["retain"]])) {
    retain_compressed_score_sample(
      fit, score_X, yhat, compressed_dat,
      group_cols, from_statement, slope_names, inputs
    )
  }

  coeftable = gen_coeftable(betahat, vcov_mat, df_res)

  # Use term assignments, not name prefixes (an FE can start with a slope name).
  coef_names = intersect(slope_names, rownames(coeftable))

  return(
    list(
      linear_design = if (!is.null(compressed_layout)) native_linear_layout(
        list(columns = compressed_layout[["columns"]]), fit[["x"]], FALSE),
      coeftable = coeftable,
      compressed_design = list(
        formula = matrix_formula, trend_sources = trend_sources, trend_centers = trend_centers,
        source_names = names(compressed_dat), columns = colnames(X)
      ),
      data = restore_compressed_data(compressed_dat),
      vcov = vcov_mat,
      fml = fml,
      yvar = yvar,
      xvars = standardize_coef_names(xvars),
      collin.var = standardize_coef_names(collin_vars),
      coef_names = coef_names,
      coef_names_all = slope_names,
      fe = fe,
      weights = weights,
      query_string = query_string,
      nobs = nobs_comp,
      nobs_orig = nobs_orig,
      strategy = "compress",
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
    # Nested iff no FE value spans more than one cluster
    nested_sql = sql_limit(
      glue(
        "SELECT {f} FROM (SELECT * {from_statement}) t ",
        "GROUP BY {f} HAVING COUNT(DISTINCT {cluster_var}) > 1"
      ),
      1,
      backend
    )
    count_sql = glue("SELECT COUNT(DISTINCT {f}) AS n FROM (SELECT * {from_statement}) t")
    n_levels = tryCatch(
      if (nrow(dbGetQuery(conn, nested_sql)) == 0) dbGetQuery(conn, count_sql)[["n"]] else 0L,
      error = function(e) {
        # Fall back to keeping the FE in K (larger, conservative SEs), but say so
        warning(
          "[dbreg] Could not check whether fixed effect '", f, "' is nested in '",
          cluster_var, "' (", conditionMessage(e), "). Treating it as not nested ",
          "for the ssc = \"nested\" correction.",
          call. = FALSE
        )
        0L
      }
    )
    nested_df = nested_df + max(as.integer(n_levels) - 1L, 0L)
  }
  nested_df
}
