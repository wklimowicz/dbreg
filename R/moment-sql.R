# Shared SQL location/scale metadata. Allocate against the entire selected frame
# so a user column cannot become ambiguous after a location CROSS JOIN.
moment_sql_names = function(names, sql, conn) {
  tokens = unique(unlist(regmatches(sql, gregexpr("[A-Za-z_][A-Za-z_0-9]*", sql))))
  safe_identifiers(names, taken = tokens, fold_case = folds_identifier_case(conn))
}

# One row per design cell, with centered outcome statistics and normalized
# weights. Grouping keys keep their SQL types; only numeric aggregates change.
compressed_moment_query = function(conn, from, groups, y, weights, intercept = TRUE, fit_columns = TRUE) {
  private = moment_sql_names(c("dbreg_cells", "dbreg_cell_location", "dbreg_cell_anchor",
    "dbreg_cell_source", "dbreg_cell_center__", "dbreg_cell_scale__",
    "dbreg_cell_mean__", "dbreg_cell_within__", "dbreg_cell_anchor__"),
    paste(from, paste(c(groups, y, weights), collapse = " ")), conn)
  cells = private[[1L]]
  location = private[[2L]]
  anchor_table = private[[3L]]
  source = private[[4L]]
  center = private[[5L]]
  scale = private[[6L]]
  mean = private[[7L]]
  within = private[[8L]]
  anchor = private[[9L]]
  group_by = if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ",")) else ""
  row_id = if (length(groups)) paste0("ROW_NUMBER() OVER (ORDER BY ", paste(groups, collapse = ","), ")") else "1"
  value = sql_cast_double(y, conn)
  backend_name = detect_backend(conn)[["name"]]
  if (is.null(weights) && backend_name %in% c("duckdb", "spark")) {
    # One source aggregation: within-cell covariance avoids subtracting large
    # raw squares, then the smaller cell result supplies the global location.
    avg_fn = if (backend_name == "duckdb") "FAVG" else "AVG"
    cells_sql = paste0(cells, " AS (SELECT ", paste(c(groups, "COUNT(*) AS n",
      paste0(avg_fn, "(", value, ") AS ", mean),
      paste0("COVAR_POP(", value, ",", value, ")*COUNT(*) AS ", within)), collapse = ","),
      " ", from, group_by, ")")
    location_sql = if (intercept) paste0(anchor_table, " AS (SELECT MIN(", mean,
      ")/2.0+MAX(", mean, ")/2.0 AS ", anchor, " FROM ", cells, "),",
      location, " AS (SELECT MAX(", anchor, ")+SUM(n*(", mean, "-", anchor,
      "))/SUM(n) AS ", center, " FROM ", cells, " CROSS JOIN ", anchor_table, ")") else
      paste0(location, " AS (SELECT MAX(0.0) AS ", center, " FROM ", cells, ")")
    centered = paste0("(", mean, "-", center, ")")
    wy = paste0("n*", centered)
    yy = paste0("(", within, "+n*", centered, "*", centered, ")")
    fields = c(groups, "n", "n AS sum_w", paste(wy, "AS sum_wy"), paste(yy, "AS sum_wy_sq"),
      "n AS sum_w2", paste(wy, "AS sum_w2y"), paste(yy, "AS sum_w2y_sq"),
      center, paste("1.0 AS", scale), if (fit_columns) c(paste(centered, "AS mean_Y"), "SQRT(n) AS wts"),
      paste(row_id, "AS dbreg_cell"))
    sql = paste0("WITH ", cells_sql, ",", location_sql, " SELECT ", paste(fields, collapse = ","),
      " FROM ", cells, " CROSS JOIN ", location)
  } else {
    w = if (!is.null(weights)) sql_cast_double(weights, conn)
    location_sql = paste0("SELECT ", if (intercept) paste0("AVG(", value, ")") else "MAX(0.0)",
      " AS ", center, ",", if (is.null(w)) "MAX(1.0)" else paste0("MAX(", w, ")"), " AS ", scale, " ", from)
    projected = paste(unique(c(groups, y, weights)), collapse = ",")
    normalized_from = paste0("FROM (SELECT ", projected, " ", from, ") ", source,
      " CROSS JOIN (", location_sql, ") ", location)
    terms = build_weighted_moment_terms(conn = conn, y_sql = paste0("(", value, "-", center, ")"),
      weights_expr = if (!is.null(w)) paste0("(", w, ")/", scale), include_w_sq = TRUE,
      prefix_terms = c("COUNT(*) AS n", paste0("MAX(", center, ") AS ", center), paste0("MAX(", scale, ") AS ", scale)))
    sql = paste0("WITH ", cells, " AS (SELECT ", paste(c(groups, terms), collapse = ","), " ",
      normalized_from, group_by, ") SELECT *,",
      if (fit_columns) "sum_wy/sum_w AS mean_Y,SQRT(sum_w) AS wts,", row_id, " AS dbreg_cell FROM ", cells)
  }
  list(sql = sql, center_field = center, scale_field = scale)
}

read_compressed_moments = function(conn, query) {
  data = decode_compressed_moments(dbGetQuery(conn, query[["sql"]]))
  if (!nrow(data)) stop("No observations in the compressed estimation sample.", call. = FALSE)
  attr(data, "moment_center") = as.numeric(data[[query[["center_field"]]]][[1L]])
  attr(data, "weight_scale") = as.numeric(data[[query[["scale_field"]]]][[1L]])
  data[c(query[["center_field"]], query[["scale_field"]])] = NULL
  data
}

# Preserve the existing raw data_only/object$data contract. Estimation and
# inference use normalized statistics before this conversion.
restore_compressed_data = function(data) {
  center = attr(data, "moment_center")
  scale = attr(data, "weight_scale")
  if (is.null(center)) return(data)
  for (suffix in c("", "2")) {
    w = paste0("sum_w", suffix)
    wy = paste0(w, "y")
    yy = paste0(wy, "_sq")
    multiplier = scale^(if (suffix == "2") 2 else 1)
    data[[yy]] = (data[[yy]]+2*center*data[[wy]]+center^2*data[[w]])*multiplier
    data[[wy]] = (data[[wy]]+center*data[[w]])*multiplier
    data[[w]] = data[[w]]*multiplier
  }
  if ("mean_Y" %in% names(data)) data[["mean_Y"]] = data[["mean_Y"]]+center
  if ("wts" %in% names(data)) data[["wts"]] = sqrt(data[["sum_w"]])
  attr(data, "moment_center") = NULL
  attr(data, "weight_scale") = NULL
  data
}
