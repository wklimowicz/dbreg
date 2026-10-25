# Estimation frames
#
# A frame is the CTE that both the moments pass and the score pass read from.
# Each build_*_frame() returns a list with:
#   cte_sql, cte_name  CTE text and the name of its final table
#   cols               named character: variable name -> SQL column in cte_name
#   weights_expr       SQL weight expression in cte_name (NULL if unweighted)
#   intercept          whether the frame keeps a constant (FALSE once demeaned)
#   cast               whether moments CAST columns before aggregating
#   prefix_terms       count aggregates (n_total, n_fe1, ...)
#   cleanup            function dropping any temp tables the frame created
# build_base_frame() also returns base_sql, the bare SELECT behind the CTE.

# SELECT list: id columns, the outcome, aliased regressors, then any extra
# (weights, cluster) columns not already selected
frame_select = function(y, x_sql, x_names, id_cols = NULL, extra = NULL) {
  c(
    id_cols,
    y,
    sprintf("%s AS %s", x_sql, x_names),
    setdiff(unique(extra), c(id_cols, y, x_names))
  )
}

# Global location/weight normalization shared by compressed OLS and IV.
# Grouped sums then retain the within-cell variation instead of subtracting
# raw outcome squares at the scale of a possibly enormous intercept.
compression_moment_source = function(conn, from, y, weights, intercept = TRUE) {
  w = sql_weight_expr(weights)
  location = paste0("SELECT ", if (intercept) paste0("AVG(", sql_cast_double(y, conn), ")") else "MAX(0.0)",
    " AS dbreg_cell_center__, ", if (is.null(w)) "MAX(1.0)" else paste0("MAX(", w, ")"),
    " AS dbreg_cell_scale__ ", from)
  list(from = paste0("FROM (SELECT * ", from, ") dbreg_cell_source CROSS JOIN (", location, ") dbreg_cell_location"),
    y = paste0("(", sql_cast_double(y, conn), "-dbreg_cell_center__)"),
    weights = if (!is.null(w)) paste0("(", w, ")/dbreg_cell_scale__"),
    prefix = c("MAX(dbreg_cell_center__) AS dbreg_cell_center__", "MAX(dbreg_cell_scale__) AS dbreg_cell_scale__"))
}

restore_compressed_data = function(dat) {
  c = dat[["dbreg_cell_center__"]]
  scale = dat[["dbreg_cell_scale__"]]
  if (is.null(c)) return(dat)
  for (suffix in c("", "2")) {
    w = paste0("sum_w", suffix)
    wy = paste0(w, "y")
    yy = paste0(wy, "_sq")
    if (!all(c(w, wy, yy) %in% names(dat))) next
    dat[[yy]] = (dat[[yy]]+2*c*dat[[wy]]+c*c*dat[[w]])*scale^(if (suffix == "2") 2 else 1)
    dat[[wy]] = (dat[[wy]]+c*dat[[w]])*scale^(if (suffix == "2") 2 else 1)
    dat[[w]] = dat[[w]]*scale^(if (suffix == "2") 2 else 1)
  }
  if ("mean_Y" %in% names(dat)) dat[["mean_Y"]] = dat[["mean_Y"]]+c
  if ("wts" %in% names(dat)) dat[["wts"]] = sqrt(dat[["sum_w"]])
  dat
}

# Translate the cell design's actual contrasts into SQL. This also handles
# numeric FE labels and exact integer dictionaries without string round-trips.
compressed_score_expressions = function(conn, formula, dat, X, classes) {
  assign = attr(X, "assign")
  terms = attr(stats::terms(formula), "term.labels")
  out = vapply(seq_len(ncol(X)), function(j) {
    if (assign[[j]] == 0L) return("1.0")
    variable = terms[[assign[[j]]]]
    values = dat[[variable]]
    if (!is.factor(values) && !is.character(values) && !is.logical(values)) return(variable)
    levels = if (is.factor(values)) levels(values) else sort(unique(values))
    first = match(as.character(levels), as.character(values))
    raw = attr(values, "dbreg_key_values")
    if (is.null(raw)) raw = levels
    type = if (identical(attr(values, "dbreg_key_type"), "integer")) "integer64" else classes[[variable]]
    literals = sql_factor_literals(conn, setNames(list(raw), variable), setNames(type, variable))[[variable]]
    loads = as.numeric(X[first, j])
    keep = which(!is.na(first) & is.finite(loads) & loads != 0)
    if (!length(keep)) return("0.0")
    paste0("(CASE ", paste0("WHEN ", variable, "=", literals[keep], " THEN ", dbglm_number(loads[keep]), collapse = " "), " ELSE 0.0 END)")
  }, character(1))
  setNames(out, colnames(X))
}

# Rank distinct typed FE tuples, then join their IDs to the rows. Sorting the
# bounded level dictionaries avoids a window sort over the full observation set.
sql_rank_fixed_effects = function(conn, source, components, fe) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  joins = vapply(seq_along(fe), function(i) {
    columns = q(unique(components[[i]]))
    keys = paste(columns, collapse = ",")
    alias = paste0("__dbglm_d", i)
    paste0(" JOIN (SELECT ", keys, ", DENSE_RANK() OVER (ORDER BY ", keys,
      ") AS ", fe[i], " FROM (SELECT DISTINCT ", keys, " FROM ", source, ") __dbglm_k) ", alias, " ON ",
      paste(sprintf("__dbglm_r.%s = %s.%s", columns, alias, columns), collapse = " AND "))
  }, character(1))
  paste0("SELECT __dbglm_r.*, ", paste(sprintf("__dbglm_d%d.%s", seq_along(fe), fe), collapse = ", "),
    " FROM ", source, " __dbglm_r", paste(joins, collapse = ""))
}

#' Untransformed estimation frame
#'
#' Selects the outcome, regressors (aliased to `x_names`), weights and cluster
#' columns into a `base` CTE.
#'
#' @param inputs Inputs environment.
#' @param x_sql,x_names SQL expressions and column names for the regressors.
#' @keywords internal
build_base_frame = function(inputs, x_sql, x_names) {
  env2env(
    inputs,
    environment(),
    c("yvar", "weights", "cluster_var", "from_statement", "conn", "intercept")
  )
  y = yvar
  sel = frame_select(y, x_sql, x_names, extra = c(weights, cluster_var, inputs[["retained_key_safe"]]))
  base_sql = paste0("SELECT ", paste(sel, collapse = ", "), " ", from_statement)
  list(
    cte_sql = paste0("WITH base AS (\n  ", base_sql, "\n)"),
    cte_name = "base",
    base_sql = base_sql,
    cols = setNames(c(y, x_names), c(y, x_names)),
    weights_expr = sql_weight_expr(weights),
    intercept = !isFALSE(intercept),
    cast = TRUE,
    prefix_terms = sql_count(conn, "n_total"),
    cleanup = function() invisible(NULL)
  )
}

#' Demeaned (within) estimation frame
#'
#' Sweeps out the fixed effects from the outcome and regressors, producing
#' `<var>_tilde` columns in a `demeaned` CTE. One FE and balanced, unweighted
#' two-way panels are demeaned analytically in SQL. Otherwise (3+ FE, or two
#' FE that are weighted or unbalanced) fixed effects are projected into a
#' temp table, which `cleanup()` drops.
#'
#' @inheritParams build_base_frame
#' @param tag Prefix for messages, e.g. `"dbreg"`.
#' @keywords internal
build_demean_frame = function(inputs, x_sql, x_names, tag = "dbreg") {
  env2env(
    inputs,
    environment(),
    c(
      "conn", "fe", "yvar", "weights", "cluster_var", "verbose", "sql_only",
      "from_statement"
    )
  )
  y = yvar
  vars = c(y, x_names)
  w_base = sql_weight_expr(weights)

  use_ap = length(fe) > 2
  if (length(fe) == 2) {
    is_balanced = inputs[["is_balanced"]]
    if (is.null(is_balanced)) {
      is_balanced = dbreg_is_balanced_panel(inputs)
    }
    use_ap = !is.null(weights) || !isTRUE(is_balanced)
  }
  if (use_ap && verbose) {
    message(sprintf("[%s] Using fixed-effect projection for FE demeaning", tag))
  }
  if (use_ap && sql_only) {
    return(dbreg_ap_query_plan(inputs, x_sql, x_names))
  }

  cleanup = function() invisible(NULL)
  if (use_ap) {
    ap_res = dbreg_alternating_projections(inputs, x_sql, x_names)
    ap_tables = unique(c(ap_res[["table"]], ap_res[["base_table"]], ap_res[["all_tables"]]))
    cleanup = function() {
      backend = detect_backend(conn)[["name"]]
      for (tbl in ap_tables) {
        drop_table_if_exists(conn, tbl, backend)
      }
    }
    cte_sql = paste0("WITH demeaned AS (SELECT * FROM ", ap_res[["table"]], ")")
    weights_expr = "__w"
  } else {
    # Analytic demeaning: subtract group means for each FE (and add back the
    # overall mean for two-way panels)
    means = function(alias, suffix, group = NULL) {
      cols = vapply(vars, function(v) sql_weighted_mean(v, w_base, paste0(v, suffix)), character(1))
      sprintf(
        "%s AS (SELECT %s FROM base%s)",
        alias,
        paste(c(group, cols), collapse = ", "),
        if (is.null(group)) "" else paste0(" GROUP BY ", group)
      )
    }
    if (length(fe) == 1) {
      mean_ctes = means("group_means", "_mean", fe[1])
      tilde = sprintf("(b.%s - gm.%s_mean) AS %s_tilde", vars, vars, vars)
      joins = sprintf("JOIN group_means gm ON b.%s = gm.%s", fe[1], fe[1])
    } else {
      mean_ctes = c(
        means("unit_means", "_u", fe[1]),
        means("time_means", "_t", fe[2]),
        means("overall", "_o")
      )
      tilde = sprintf("(b.%s - um.%s_u - tm.%s_t + o.%s_o) AS %s_tilde", vars, vars, vars, vars, vars)
      joins = c(
        sprintf("JOIN unit_means um ON b.%s = um.%s", fe[1], fe[1]),
        sprintf("JOIN time_means tm ON b.%s = tm.%s", fe[2], fe[2]),
        "CROSS JOIN overall o"
      )
    }
    passthrough = c(
      if (!is.null(weights)) sprintf("b.%s AS weights", weights),
      if (!is.null(cluster_var) && !cluster_var %in% fe) sprintf("b.%s AS %s", cluster_var, cluster_var),
      if (length(inputs[["retained_key_safe"]])) sprintf("b.%s AS %s", inputs[["retained_key_safe"]], inputs[["retained_key_safe"]])
    )
    sel = frame_select(y, x_sql, x_names, id_cols = fe, extra = c(weights, cluster_var, inputs[["retained_key_safe"]]))
    cte_sql = paste0(
      "WITH base AS (\n  SELECT ", paste(sel, collapse = ", "), " ", from_statement, "\n),\n",
      paste(mean_ctes, collapse = ",\n"), ",\n",
      "demeaned AS (\n  SELECT ", paste(c(paste0("b.", fe), tilde, passthrough), collapse = ",\n    "),
      "\n  FROM base b\n  ", paste(joins, collapse = "\n  "), "\n)"
    )
    weights_expr = if (is.null(weights)) NULL else sql_weight_expr("weights")
  }

  list(
    cte_sql = cte_sql,
    cte_name = "demeaned",
    cols = setNames(paste0(vars, "_tilde"), vars),
    weights_expr = weights_expr,
    intercept = FALSE,
    cast = TRUE,
    prefix_terms = c(
      sql_count(conn, "n_total"),
      vapply(seq_along(fe), function(k) sql_count(conn, sprintf("n_fe%d", k), fe[k], distinct = TRUE), character(1))
    ),
    cleanup = cleanup,
    ap_res = if (use_ap) ap_res else NULL,
    rank_source = paste0("(SELECT ", paste(frame_select(yvar, x_sql, x_names), collapse = ","), " ", from_statement, ") dbreg_rank_source")
  )
}

#' Moments query over an estimation frame
#'
#' Appends the cross-moment aggregates for `y` and `x` (names in `frame$cols`)
#' to the frame's CTE. Read the result back with [cross_moments()].
#'
#' @param frame Result of a `build_*_frame()` function.
#' @param y,x Outcome and regressor names.
#' @param conn Database connection (for the cast type).
#' @keywords internal
frame_moments_sql = function(frame, y, x, conn) {
  col = function(v) {
    sql = unname(frame[["cols"]][v])
    if (frame[["cast"]]) sql_cast_double(sql, conn) else sql
  }
  # Center before multiplying: reconstructing a variance from raw squares loses
  # all useful digits when a variable has a large location and small variation.
  variables = c(x, y)
  center_fields = paste0("dbreg_moment_center_", seq_along(variables))
  rank_fields = paste0("dbreg_moment_original_", seq_along(variables))
  tokens = c(frame[["cols"]], frame[["prefix_terms"]], frame[["weights_expr"]])
  taken = unique(unlist(regmatches(tokens, gregexpr("[A-Za-z_][A-Za-z_0-9]*", tokens))))
  private = safe_identifiers(c(center_fields, "dbreg_moment_weight_scale", rank_fields), taken = taken, fold_case = TRUE)
  centers = private[seq_along(variables)]
  weight_field = private[[length(variables)+1L]]
  rank_names = utils::tail(private, length(variables))
  center_terms = if (isTRUE(frame[["intercept"]])) paste0("AVG(", vapply(variables, col, character(1)), ")") else rep("MAX(0.0)", length(variables))
  w = frame[["weights_expr"]]
  scale = if (is.null(w)) "MAX(1.0)" else paste0("MAX(", sql_cast_double(w, conn), ")")
  stats = paste0(",\ndbreg_moment_location AS (SELECT ", paste(c(paste(center_terms, "AS", centers),
    paste(scale, "AS", weight_field)), collapse = ","), " FROM ", frame[["cte_name"]], ")")
  centered = setNames(paste0("(", vapply(variables, col, character(1)), "-", centers, ")"), variables)
  rank_terms = character()
  rank_join = ""
  if (!is.null(frame[["rank_source"]])) {
    # Projection noise must not be normalized into an identified column. Bound
    # that noise in each variable's original units, rather than using the scale
    # of an unrelated regressor as an absolute rank threshold.
    stats = paste0(stats, ",\ndbreg_moment_original AS (SELECT ",
      paste(paste0("MAX(ABS(", sql_cast_double(variables, conn), ")) AS ", rank_names), collapse = ","),
      " FROM ", frame[["rank_source"]], ")")
    rank_terms = paste0("MAX(", rank_names, ") AS ", rank_fields)
    rank_join = " CROSS JOIN dbreg_moment_original"
  }
  terms = build_weighted_moment_terms(
    y_sql = centered[[y]],
    x_sql = unname(centered[x]),
    weights_expr = if (!is.null(w)) paste0("(", sql_cast_double(w, conn), ")/", weight_field),
    intercept = frame[["intercept"]],
    prefix_terms = frame[["prefix_terms"]]
  )
  paste0(
    frame[["cte_sql"]], stats, ",\nmoments AS (\n  SELECT\n    ",
    paste(c(terms, rank_terms, paste0("MAX(", c(centers, weight_field), ") AS ", c(center_fields, "dbreg_moment_weight_scale"))), collapse = ",\n    "),
    "\n  FROM ", frame[["cte_name"]], " CROSS JOIN dbreg_moment_location", rank_join, "\n)\nSELECT * FROM moments"
  )
}
