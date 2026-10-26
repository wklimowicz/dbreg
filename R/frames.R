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

# Unweighted cells need one N-row aggregation. Center the cell means in a
# second, small aggregation; within-cell covariance preserves outcome variation.
compressed_cell_moments_sql = function(conn, from, y, groups, weights = NULL,
                                       intercept = TRUE, fit_columns = TRUE) {
  if (!is.null(weights) || !identical(detect_backend(conn)[["name"]], "duckdb")) return(NULL)
  private = safe_identifiers(c("dbreg_cell_mean__", "dbreg_cell_within__", "dbreg_cell_anchor__"),
    taken = c(groups, y), fold_case = folds_identifier_case(conn))
  mean = private[[1L]]
  within = private[[2L]]
  anchor = private[[3L]]
  value = sql_cast_double(y, conn)
  cells = paste0("dbreg_cells AS (SELECT ", paste(c(groups, "COUNT(*) AS n",
    paste0("FAVG(", value, ") AS ", mean),
    paste0("COVAR_POP(", value, ",", value, ")*COUNT(*) AS ", within)), collapse = ","), " ", from,
    if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ",")), ")")
  location = if (intercept) paste0(
    "dbreg_cell_anchor AS (SELECT MIN(", mean, ")/2.0+MAX(", mean, ")/2.0 AS ", anchor, " FROM dbreg_cells),",
    "dbreg_cell_location AS (SELECT MAX(", anchor, ")+SUM(n*(", mean, "-", anchor,
    "))/SUM(n) AS dbreg_cell_center__ FROM dbreg_cells CROSS JOIN dbreg_cell_anchor)") else
      "dbreg_cell_location AS (SELECT MAX(0.0) AS dbreg_cell_center__ FROM dbreg_cells)"
  centered = paste0("(", mean, "-dbreg_cell_center__)")
  wy = paste0("n*", centered)
  yy = paste0("(", within, "+n*", centered, "*", centered, ")")
  fields = c(groups, "n", "n AS sum_w", paste(wy, "AS sum_wy"), paste(yy, "AS sum_wy_sq"),
    "n AS sum_w2", paste(wy, "AS sum_w2y"), paste(yy, "AS sum_w2y_sq"),
    "dbreg_cell_center__", "1.0 AS dbreg_cell_scale__",
    if (fit_columns) c(paste(centered, "AS mean_Y"), "SQRT(n) AS wts"))
  paste0("WITH ", cells, ",", location, " SELECT ", paste(fields, collapse = ","),
    " FROM dbreg_cells CROSS JOIN dbreg_cell_location")
}

# Translate the cell design's actual contrasts into SQL. This also handles
# numeric FE labels and exact integer dictionaries without string round-trips.
compressed_score_expressions = function(conn, formula, dat, X, classes, numeric_sql = character()) {
  assign = attr(X, "assign")
  terms = attr(stats::terms(formula), "term.labels")
  out = vapply(seq_len(ncol(X)), function(j) {
    if (assign[[j]] == 0L) return("1.0")
    variables = split_interaction_term(terms[[assign[[j]]]])
    categorical = variables[vapply(dat[variables], function(value) {
      is.factor(value) || is.character(value) || is.logical(value)
    }, logical(1))]
    numeric = setdiff(variables, categorical)
    product_sql = if (length(numeric)) paste(vapply(numeric, function(name) {
      if (name %in% names(numeric_sql)) numeric_sql[[name]] else sql_cast_double(name, conn)
    }, character(1)), collapse = "*") else "1.0"
    if (!length(categorical)) return(product_sql)
    product = Reduce(`*`, lapply(dat[numeric], as.numeric), init = rep(1, nrow(dat)))
    column = as.numeric(X[, j])
    rows = which(is.finite(product) & product != 0 & column != 0)
    if (!length(rows)) return("0.0")
    keys = tuple_key(dat[rows, , drop = FALSE], categorical)
    first = rows[!duplicated(keys)]
    loads = column[first] / product[first]
    conditions = lapply(categorical, function(variable) {
      values = dat[[variable]]
      levels = if (is.factor(values)) levels(values) else sort(unique(values))
      raw = attr(values, "dbreg_key_values")
      if (is.null(raw)) raw = levels
      type = if (identical(attr(values, "dbreg_key_type"), "integer")) "integer64" else classes[[variable]]
      literals = sql_factor_literals(conn, setNames(list(raw), variable), setNames(type, variable))[[variable]]
      paste0(variable, "=", literals[match(as.character(values[first]), as.character(levels))])
    })
    condition = Reduce(function(a, b) paste(a, b, sep = " AND "), conditions)
    paste0("(", product_sql, ")*(CASE ", paste0("WHEN ", condition, " THEN ",
      dbglm_number(loads), collapse = " "), " ELSE 0.0 END)")
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

  varying = length(unlist(inputs[["fe_slopes"]])) > 0L
  use_ap = length(fe) > 2 || varying
  if (length(fe) == 2 && !varying) {
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
    if (varying) stop("SQL-only projection plans for varying slopes are not implemented.", call. = FALSE)
    return(dbreg_ap_query_plan(inputs, x_sql, x_names))
  }

  cleanup = function() invisible(NULL)
  rank_columns = rank_values = NULL
  rank_names = safe_identifiers(paste0("dbreg_frame_original_", seq_along(vars)),
    taken = c(vars, fe, weights, cluster_var, inputs[["retained_key_safe"]]), fold_case = folds_identifier_case(conn))
  if (use_ap) {
    ap_res = dbreg_alternating_projections(inputs, x_sql, x_names)
    ap_tables = unique(c(if (!isTRUE(ap_res[["inline"]])) ap_res[["table"]], ap_res[["base_table"]], ap_res[["all_tables"]]))
    cleanup = function() {
      backend = detect_backend(conn)[["name"]]
      for (tbl in ap_tables) {
        drop_table_if_exists(conn, tbl, backend)
      }
    }
    cte_sql = if (!is.null(ap_res[["table"]])) paste0("WITH demeaned AS (SELECT * FROM ", ap_res[["table"]], ")") else NULL
    weights_expr = if (is.null(weights)) NULL else if (is.null(ap_res[["weight_column"]])) "__w" else ap_res[["weight_column"]]
    rank_values = ap_res[["original_scale"]]
  } else {
    # Analytic demeaning: subtract group means for each FE (and add back the
    # overall mean for two-way panels)
    means = function(alias, suffix, group = NULL) {
      cols = vapply(vars, function(v) sql_weighted_mean(v, w_base, paste0(v, suffix)), character(1))
      sprintf(
        "%s AS (SELECT %s FROM base%s)",
        alias,
        paste(c(group, cols, if (alias %in% c("group_means", "unit_means"))
          paste0("MAX(ABS(", vapply(vars, sql_cast_double, character(1), conn = conn), ")) AS ", rank_names)), collapse = ", "),
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
    rank_columns = setNames(rank_names, vars)
    passthrough = c(
      paste0(if (length(fe) == 1L) "gm." else "um.", rank_names, " AS ", rank_names),
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
    cols = if (use_ap && !is.null(ap_res[["cols"]])) ap_res[["cols"]] else setNames(paste0(vars, "_tilde"), vars),
    weights_expr = weights_expr,
    intercept = FALSE,
    cast = TRUE,
    prefix_terms = c(
      if (use_ap && !is.null(ap_res[["nobs"]])) paste(dbglm_number(ap_res[["nobs"]]), "AS n_total") else sql_count(conn, "n_total"),
      vapply(seq_along(fe), function(k) {
        if (use_ap && !is.null(ap_res[["n_fe_levels"]])) paste(dbglm_number(ap_res[["n_fe_levels"]][[k]]), "AS", sprintf("n_fe%d", k)) else
          sql_count(conn, sprintf("n_fe%d", k), fe[k], distinct = TRUE)
      }, character(1))
    ),
    cleanup = cleanup,
    ap_res = if (use_ap) ap_res else NULL,
    rank_columns = rank_columns, rank_values = rank_values,
    moments = if (use_ap) ap_res[["moments"]] else NULL,
    rank_source = paste0("(SELECT ", paste(frame_select(yvar, x_sql, x_names), collapse = ","), " ", from_statement, ") dbreg_rank_source")
  )
}

# Stable unweighted centered moments in one aggregation pass. Built-in
# covariance accumulators avoid subtracting raw squares at large locations.
build_centered_moment_terms = function(y, x, prefix = character()) {
  product = function(a, b, name) paste0("COVAR_POP(", a, ",", b, ")*COUNT(*) AS ", name)
  terms = c(prefix, "COUNT(*) AS sum_w", "MAX(0.0) AS sum_wy", product(y, y, "sum_wy_sq"))
  for (i in seq_along(x)) {
    terms = c(terms, paste0("MAX(0.0) AS sum_w", i), product(x[[i]], y, paste0("sum_w", i, "_y")))
    for (j in seq_len(i)) terms = c(terms, product(x[[j]], x[[i]], paste0("sum_w", j, "_", i)))
  }
  terms
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
  if (!is.null(frame[["moments"]])) {
    row = frame[["moments"]]
    return(paste("SELECT", paste(paste(dbglm_number(as.numeric(row[1L, ])), "AS",
      as.character(DBI::dbQuoteIdentifier(conn, names(row)))), collapse = ",")))
  }
  col = function(v) {
    sql = unname(frame[["cols"]][v])
    if (frame[["cast"]]) sql_cast_double(sql, conn) else sql
  }
  # The common unweighted OLS path needs neither a location scan nor a
  # location CROSS JOIN when the backend has stable covariance aggregates.
  w = frame[["weights_expr"]]
  backend = detect_backend(conn)[["name"]]
  variables = c(x, y)
  if (is.null(w) && isTRUE(frame[["intercept"]]) && is.null(frame[["rank_source"]]) &&
      identical(backend, "duckdb")) {
    terms = build_centered_moment_terms(col(y), vapply(x, col, character(1)), frame[["prefix_terms"]])
    location = c(paste0(if (backend == "duckdb") "FAVG(" else "AVG(", vapply(variables, col, character(1)), ") AS dbreg_moment_center_", seq_along(variables)),
      "MAX(1.0) AS dbreg_moment_weight_scale")
    return(paste0(frame[["cte_sql"]], " SELECT ", paste(c(terms, location), collapse = ","), " FROM ", frame[["cte_name"]]))
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
  location_join = " CROSS JOIN dbreg_moment_location"
  location_terms = paste0("MAX(", c(centers, weight_field), ") AS ", c(center_fields, "dbreg_moment_weight_scale"))
  if (!isTRUE(frame[["intercept"]]) && is.null(w)) {
    stats = ""
    centered = setNames(vapply(variables, col, character(1)), variables)
    location_join = ""
    location_terms = c(paste0("MAX(0.0) AS ", center_fields), "MAX(1.0) AS dbreg_moment_weight_scale")
  }
  rank_terms = character()
  rank_join = ""
  if (!is.null(frame[["rank_values"]])) {
    rank_terms = paste0("MAX(", dbglm_number(unname(frame[["rank_values"]][variables])), ") AS ", rank_fields)
  } else if (!is.null(frame[["rank_columns"]])) {
    rank_terms = paste0("MAX(", unname(frame[["rank_columns"]][variables]), ") AS ", rank_fields)
  } else if (!is.null(frame[["rank_source"]])) {
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
    paste(c(terms, rank_terms, location_terms), collapse = ",\n    "),
    "\n  FROM ", frame[["cte_name"]], location_join, rank_join, "\n)\nSELECT * FROM moments"
  )
}
