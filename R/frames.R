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
  sel = frame_select(y, x_sql, x_names, extra = c(weights, cluster_var, if (!is.null(inputs[["retained_original"]])) inputs[["row_id_safe"]]))
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
#' FE that are weighted or unbalanced) a bounded FE sufficient-statistic solve
#' keeps transformed rows inline over an owned source snapshot. An expected
#' decline falls back to observation-table AP. `cleanup()` releases owned state.
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

  varying = length(unlist(inputs[["fe_slopes"]])) > 0L
  use_ap = length(fe) > 2 || varying
  if (length(fe) == 2 && !varying) {
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
    if (varying) stop("SQL-only projection plans for varying slopes are not implemented.", call. = FALSE)
    return(dbreg_ap_query_plan(inputs, x_sql, x_names))
  }

  cleanup = function() invisible(NULL)
  if (use_ap) {
    bound = projection_frame_inputs(inputs, x_sql, x_names)
    ap_res = dbreg_alternating_projections(bound[["inputs"]], bound[["x"]], bound[["x"]])
    inputs[["projection_diagnostics"]] = bound[["inputs"]][["projection_diagnostics"]]
    inputs[["projection_limits"]] = bound[["inputs"]][["projection_limits"]]
    names(ap_res[["original_scale"]]) = vars
    inputs[["ap_diagnostics"]] = ap_res[intersect(c("iterations", "accelerated_steps", "max_abs", "max_relative", "method", "transfer_cells", "summary_rows", "upload_bytes", "estimated_bytes", "solver", "solver_iterations", "preconditioner", "transport", "summary_statements", "wire_rows", "wire_bytes", "observation_tables", "transformed_observation_tables", "moment_method", "fwl_recomputed", "peak_tables", "cg_restarts", "cg_retries"), names(ap_res))]
    ap_tables = c(ap_res[["table"]], ap_res[["base_table"]])
    cleanup = if (!is.null(ap_res[["resources"]])) function() release_fit_resources(ap_res[["resources"]]) else function() {
      backend = detect_backend(conn)[["name"]]
      for (tbl in ap_tables) {
        drop_table_if_exists(conn, tbl, backend)
      }
    }
    mapped = projection_frame_columns(conn, ap_res[["table"]], bound, vars)
    cte_sql = paste0("WITH demeaned AS (", mapped[["sql"]], ")")
    projected_columns = mapped[["columns"]]
    weights_expr = mapped[["weights"]]
  } else {
    # Analytic demeaning: subtract group means for each FE (and add back the
    # overall mean for two-way panels).
    keys = if (!is.null(inputs[["retained_original"]])) inputs[["row_id_safe"]]
    sel = frame_select(y, x_sql, x_names, id_cols = fe, extra = c(weights, cluster_var, keys))
    base_sql = paste0("SELECT ", paste(sel, collapse = ", "), " ", from_statement)
    centered = function(v, prefix = "b.") sql_cast_double(paste0(prefix, v), conn)
    local_from = "base b"
    # Demeaning roundoff scales with the raw magnitudes, so collinearity is
    # judged against the largest absolute value, read off the same GROUP BY
    # as the means and carried through the existing join (no extra scan).
    first = if (length(fe) == 1) "gm" else "um"
    rank_names = if (length(x_names)) setNames(safe_identifiers(paste0(x_names, "_half"),
      taken = c(fe, keys, cluster_var, weights, paste0(vars, "_tilde"))), x_names) else NULL
    analytic_weights = if (is.null(weights)) NULL else sql_weight_expr(paste0("b.", weights))
    means = function(alias, suffix, group = NULL, ranges = FALSE) {
      cols = vapply(vars, function(v) sql_weighted_mean(centered(v), analytic_weights, paste0(v, suffix)), character(1))
      if (ranges && length(x_names)) cols = c(cols, paste0("MAX(ABS(", centered(x_names), ")) AS ", rank_names))
      sprintf(
        "%s AS (SELECT %s FROM %s%s)",
        alias,
        paste(c(if (length(group)) paste0("b.", group), cols), collapse = ", "),
        local_from,
        if (is.null(group)) "" else paste0(" GROUP BY b.", group)
      )
    }
    if (length(fe) == 1) {
      mean_ctes = means("group_means", "_mean", fe[1], ranges = TRUE)
      tilde = sprintf("(%s - gm.%s_mean) AS %s_tilde", centered(vars, "b."), vars, vars)
      joins = sprintf("JOIN group_means gm ON b.%s = gm.%s", fe[1], fe[1])
    } else {
      mean_ctes = c(
        means("unit_means", "_u", fe[1], ranges = TRUE),
        means("time_means", "_t", fe[2]),
        means("overall", "_o")
      )
      tilde = sprintf("(%s - um.%s_u - tm.%s_t + o.%s_o) AS %s_tilde", centered(vars, "b."), vars, vars, vars, vars)
      joins = c(
        sprintf("JOIN unit_means um ON b.%s = um.%s", fe[1], fe[1]),
        sprintf("JOIN time_means tm ON b.%s = tm.%s", fe[2], fe[2]),
        "CROSS JOIN overall o"
      )
    }
    weight_alias = safe_identifiers("weights", taken = c(fe, keys, cluster_var))
    passthrough = c(
      if (!is.null(weights)) sprintf("b.%s AS %s", weights, weight_alias),
      if (!is.null(cluster_var) && !cluster_var %in% fe) sprintf("b.%s AS %s", cluster_var, cluster_var)
    )
    key_columns = setdiff(keys, c(fe, cluster_var))
    if (length(key_columns)) passthrough = c(passthrough, paste0("b.", key_columns))
    cte_sql = paste0(
      "WITH base AS (\n  ", base_sql, "\n),\n",
      paste(mean_ctes, collapse = ",\n"), ",\n",
      "demeaned AS (\n  SELECT ", paste(c(paste0("b.", fe), tilde, passthrough), collapse = ",\n    "),
      "\n  FROM ", local_from, "\n  ", paste(joins, collapse = "\n  "), "\n)"
    )
    weights_expr = if (is.null(weights)) NULL else sql_weight_expr(weight_alias)
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
    cols = if (use_ap) projected_columns else setNames(paste0(vars, "_tilde"), vars),
    weights_expr = weights_expr,
    intercept = FALSE,
    cast = TRUE,
    prefix_terms = c(
      sql_count(conn, "n_total"),
      fe_count_terms
    ),
    cleanup = cleanup,
    moments = if (use_ap) ap_res[["moments"]] else NULL,
    rank_values = if (use_ap) ap_res[["original_scale"]] else NULL,
    rank_columns = if (!use_ap && length(x_names)) rank_names else NULL,
    rank_cte = if (!use_ap && length(x_names)) (if (length(fe) == 1) "group_means" else "unit_means") else NULL
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
  private = moment_sql_names(c("dbreg_moment_location", fields, "dbreg_moment_weight_scale",
    "dbreg_moment_original", paste0("dbreg_moment_original_", seq_along(x))), frame[["cte_sql"]], conn)
  location_name = private[[1L]]
  centers = private[seq_along(variables)+1L]
  weight_field = private[[length(variables)+2L]]
  original_name = private[[length(variables)+3L]]
  original_fields = utils::tail(private, length(x))
  center_terms = if (isTRUE(frame[["intercept"]])) paste0("MIN(", vapply(variables, col, character(1)), ")/2+MAX(", vapply(variables, col, character(1)), ")/2") else rep("MAX(0.0)", length(variables))
  scale = if (is.null(w)) "MAX(1.0)" else paste0("MAX(", sql_cast_double(w, conn), ")")
  stats = paste0(",\n", location_name, " AS (SELECT ", paste(c(paste(center_terms, "AS", centers),
    paste(scale, "AS", weight_field)), collapse = ","), " FROM ", frame[["cte_name"]], ")")
  centered = setNames(paste0("(", vapply(variables, col, character(1)), "-", centers, ")"), variables)
  location_join = paste(" CROSS JOIN", location_name)
  location_terms = paste0("MAX(", c(centers, weight_field), ") AS ", c(fields, "dbreg_moment_weight_scale"))
  # Weight rescaling only guards overflow; demeaned FE frames use the stored weights.
  fe_frame = identical(frame[["cte_name"]], "demeaned")
  if (!isTRUE(frame[["intercept"]]) && (is.null(w) || fe_frame)) {
    stats = ""
    location_join = ""
    centered = setNames(vapply(variables, col, character(1)), variables)
    location_terms = c(paste0("MAX(0.0) AS ", fields), "MAX(1.0) AS dbreg_moment_weight_scale")
  }
  rank_terms = character()
  rank_join = ""
  rank_labels = paste0("dbreg_moment_original_", seq_along(x))
  if (length(x) && !is.null(frame[["rank_values"]])) {
    original = sql_cast_double(sprintf("%.17e", unname(frame[["rank_values"]][x])), conn)
    rank_terms = paste0("MAX(", original, ") AS ", rank_labels)
  } else if (length(x) && !is.null(frame[["rank_columns"]])) {
    if (!is.null(frame[["rank_cte"]])) {
      rank_terms = paste0("(SELECT MAX(", sql_cast_double(unname(frame[["rank_columns"]][x]), conn),
        ") FROM ", frame[["rank_cte"]], ") AS ", rank_labels)
    } else {
      rank_terms = paste0("MAX(", sql_cast_double(unname(frame[["rank_columns"]][x]), conn), ") AS ", rank_labels)
    }
  } else if (length(x) && !is.null(frame[["rank_source"]])) {
    values = sql_cast_double(x, conn)
    stats = paste0(stats, ",\n", original_name, " AS (SELECT ", paste(paste0("(MAX(",
      values, ")/2-MIN(", values, ")/2) AS ", original_fields), collapse = ","), " FROM ", frame[["rank_source"]], ")")
    rank_terms = paste0("MAX(", original_fields, ") AS ", rank_labels)
    rank_join = paste(" CROSS JOIN", original_name)
  }
  terms = build_weighted_moment_terms(conn = conn, y_sql = centered[[y]], x_sql = unname(centered[x]),
    weights_expr = if (!is.null(w)) paste0("(", sql_cast_double(w, conn), ")",
      if (!fe_frame || isTRUE(frame[["intercept"]])) paste0("/", weight_field)),
    intercept = frame[["intercept"]], prefix_terms = frame[["prefix_terms"]])
  if (length(frame[["origins"]])) {
    origins = setNames(rep("0.0", length(variables)), variables)
    origins[names(frame[["origins"]])] = frame[["origins"]]
    location_terms = c(location_terms, paste0("MAX(", sql_cast_double(origins, conn),
      ") AS dbreg_moment_origin_", seq_along(variables)))
  }
  paste0(frame[["cte_sql"]], stats, ",\nmoments AS (SELECT ",
    paste(c(terms, rank_terms, location_terms), collapse = ",\n"), " FROM ", frame[["cte_name"]],
    location_join, rank_join, ") SELECT * FROM moments")
}
