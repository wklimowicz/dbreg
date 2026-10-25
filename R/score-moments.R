# Aggregate cluster score cross-products in SQL. Even observation-sized cluster
# dictionaries return one matrix row; intersections are typed GROUP BY tuples.
sql_score_moments = function(from_statement, score_sql, cluster_sql = character(), projection = NULL) {
  columns = if (length(score_sql)) paste0("dbreg_score_", seq_along(score_sql)) else character()
  values = if (length(columns)) paste(if (length(cluster_sql)) paste0("SUM(", score_sql, ")") else score_sql,
    "AS", columns) else "1 AS dbreg_score_empty__"
  source = paste0("(SELECT ", paste(values, collapse = ","), " ", from_statement,
    if (length(cluster_sql)) paste(" GROUP BY", paste(cluster_sql, collapse = ",")) else "", ") dbreg_scores")
  if (!is.null(projection) && length(columns)) {
    stopifnot(nrow(projection) == length(columns), ncol(projection) == length(columns))
    expressions = vapply(seq_len(ncol(projection)), function(i)
      sql_linear_predictor(columns, projection[, i]), character(1))
    projected = paste0("dbreg_influence_", seq_along(columns))
    source = paste0("(SELECT ", paste(expressions, "AS", projected, collapse = ","), " FROM ", source, ") dbreg_influences")
    columns = projected
  }
  paste("SELECT", paste(c("COUNT(*) AS n_groups", dbglm_moment_terms(columns, "1.0")), collapse = ","), "FROM", source)
}

score_moments = function(conn, from_statement, score_sql, cluster_sql = character(), projection = NULL) {
  row = backend_get_query(conn, sql_score_moments(from_statement, score_sql, cluster_sql, projection))
  if (nrow(row) != 1L) stop("Score moment aggregation must return one row.", call. = FALSE)
  list(matrix = dbglm_read_information(lapply(row, as.numeric), length(score_sql)), groups = as.numeric(row[["n_groups"]]))
}
