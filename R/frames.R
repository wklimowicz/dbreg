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
  sel = frame_select(y, x_sql, x_names, extra = c(weights, cluster_var))
  base_sql = paste0("SELECT ", paste(sel, collapse = ", "), " ", from_statement)
  cte_qual = if (identical(detect_backend(conn)[["name"]], "duckdb")) "NOT MATERIALIZED " else ""
  list(
    cte_sql = paste0("WITH base AS ", cte_qual, "(\n  ", base_sql, "\n)"),
    cte_name = "base",
    base_sql = base_sql,
    cols = setNames(c(y, x_names), c(y, x_names)),
    weights_expr = sql_weight_expr(weights),
    intercept = !isFALSE(intercept),
    cast = FALSE,
    prefix_terms = sql_count(conn, "n_total"),
    cleanup = function() invisible(NULL)
  )
}

#' Demeaned (within) estimation frame
#'
#' Sweeps out the fixed effects from the outcome and regressors, producing
#' `<var>_tilde` columns in a `demeaned` CTE. One FE and balanced, unweighted
#' two-way panels are demeaned analytically in SQL. Otherwise (3+ FE, or two
#' FE that are weighted or unbalanced) alternating projections are run into a
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
    if (!is.null(weights)) {
      use_ap = TRUE
    } else {
      is_balanced = inputs[["is_balanced"]]
      if (is.null(is_balanced)) {
        is_balanced = dbreg_is_balanced_panel(inputs)
      }
      use_ap = !isTRUE(is_balanced)
    }
  }
  if (use_ap && verbose) {
    message(sprintf("[%s] Using alternating projections for FE demeaning", tag))
  }
  if (use_ap && sql_only) {
    return(dbreg_ap_query_plan(inputs, x_sql, x_names))
  }

  cleanup = function() invisible(NULL)
  if (use_ap) {
    ap_res = dbreg_alternating_projections(inputs, x_sql, x_names)
    inputs[["ap_diagnostics"]] = ap_res[c("iterations", "accelerated_steps", "max_abs", "max_relative", "method")]
    ap_tables = c(ap_res[["table"]], ap_res[["base_table"]])
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
      if (!is.null(cluster_var) && !cluster_var %in% fe) sprintf("b.%s AS %s", cluster_var, cluster_var)
    )
    sel = frame_select(y, x_sql, x_names, id_cols = fe, extra = c(weights, cluster_var))
    cte_sql = paste0(
      "WITH base AS (\n  SELECT ", paste(sel, collapse = ", "), " ", from_statement, "\n),\n",
      paste(mean_ctes, collapse = ",\n"), ",\n",
      "demeaned AS (\n  SELECT ", paste(c(paste0("b.", fe), tilde, passthrough), collapse = ",\n    "),
      "\n  FROM base b\n  ", paste(joins, collapse = "\n  "), "\n)"
    )
    weights_expr = if (is.null(weights)) NULL else sql_weight_expr("weights")
  }

  fe_count_terms = if (use_ap) {
    vapply(seq_along(fe), function(k) sql_count(conn, sprintf("n_fe%d", k), fe[k], distinct = TRUE), character(1))
  } else {
    count_ctes = if (length(fe) == 1) "group_means" else c("unit_means", "time_means")
    vapply(seq_along(count_ctes), function(k) {
      sprintf("(SELECT %s FROM %s) AS n_fe%d", sql_count_expr(detect_backend(conn)[["name"]]), count_ctes[k], k)
    }, character(1))
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
      fe_count_terms
    ),
    cleanup = cleanup,
    rank_values = if (use_ap) ap_res[["original_scale"]] else NULL
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
  col = function(v) sql_cast_double(unname(frame[["cols"]][v]), conn)
  variables = c(x, y)
  fields = paste0("dbreg_moment_center_", seq_along(variables))
  w = frame[["weights_expr"]]
  rank_sql = rank_join = ""
  rank_terms = character()
  if (length(x) && !is.null(frame[["rank_values"]])) {
    original = sql_cast_double(sprintf("%.17g", unname(frame[["rank_values"]][x])), conn)
    rank_terms = sprintf("MAX(%s) AS dbreg_moment_original_%d", original, seq_along(x))
  } else if (length(x) && !is.null(frame[["rank_source"]])) {
    rank_fields = paste0("dbreg_moment_original_", seq_along(x))
    rank_cte = moment_sql_names("dbreg_original", frame[["cte_sql"]], conn)[[1L]]
    rank_sql = paste0(",\n", rank_cte, " AS (SELECT ",
      paste(sprintf("(MAX(%s)/2-MIN(%s)/2) AS %s", sql_cast_double(x, conn), sql_cast_double(x, conn), rank_fields), collapse = ","),
      " FROM ", frame[["rank_source"]], ")")
    rank_terms = paste0("MAX(", rank_fields, ") AS ", rank_fields)
    rank_join = paste0(" CROSS JOIN ", rank_cte)
  }
  private = moment_sql_names(c("dbreg_moment_location", fields, "dbreg_moment_weight_scale"),
    frame[["cte_sql"]], conn)
  location_name = private[[1L]]
  centers = private[seq_along(variables) + 1L]
  weight_field = private[[length(variables) + 2L]]
  center_terms = if (isTRUE(frame[["intercept"]])) paste0("MIN(", vapply(variables, col, character(1)), ")/2+MAX(", vapply(variables, col, character(1)), ")/2") else rep("MAX(0.0)", length(variables))
  scale = if (is.null(w)) "MAX(1.0)" else paste0("MAX(", sql_cast_double(w, conn), ")")
  stats = paste0(",\n", location_name, " AS (SELECT ", paste(c(paste(center_terms, "AS", centers),
    paste(scale, "AS", weight_field)), collapse = ","), " FROM ", frame[["cte_name"]], ")")
  centered = setNames(paste0("(", vapply(variables, col, character(1)), "-", centers, ")"), variables)
  location_join = paste(" CROSS JOIN", location_name)
  location_terms = paste0("MAX(", c(centers, weight_field), ") AS ", c(fields, "dbreg_moment_weight_scale"))
  if (!isTRUE(frame[["intercept"]]) && is.null(w)) {
    stats = ""
    location_join = ""
    centered = setNames(vapply(variables, col, character(1)), variables)
    location_terms = c(paste0("MAX(0.0) AS ", fields), "MAX(1.0) AS dbreg_moment_weight_scale")
  }

  terms = build_weighted_moment_terms(conn = conn, y_sql = centered[[y]], x_sql = unname(centered[x]),
    weights_expr = if (!is.null(w)) paste0("(", sql_cast_double(w, conn), ")/", weight_field),
    intercept = frame[["intercept"]], prefix_terms = frame[["prefix_terms"]])
  paste0(frame[["cte_sql"]], rank_sql, stats, ",\nmoments AS (SELECT ",
    paste(c(terms, location_terms, rank_terms), collapse = ",\n"), " FROM ", frame[["cte_name"]],
    rank_join, location_join, ") SELECT * FROM moments")
}

