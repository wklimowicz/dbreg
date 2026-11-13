# Ordinary generated aggregates preserve DuckDB compensation without changing
# Spark's existing double-SUM semantics. Stable working sums remain separate.
sql_sum = function(conn, expr, backend_info = detect_backend(conn)) {
  if (!length(expr)) return(character())
  aggregate = if (identical(backend_info[["name"]], "duckdb")) "FSUM" else "SUM"
  paste0(aggregate, "(", expr, ")")
}

# Emit the complete finite predicate where its owning expression is built.
# NULL stays unknown; NaN/infinities are false on DuckDB and Spark.
sql_is_finite = function(conn, expr, backend_info = detect_backend(conn)) {
  if (!length(expr)) return(character())
  backend = backend_info[["name"]]
  if (identical(backend, "duckdb")) return(paste0("isfinite(", expr, ")"))
  if (identical(backend, "spark")) return(paste0(
    "(NOT isnan(CAST((", expr, ") AS DOUBLE)) AND ABS(CAST((", expr,
    ") AS DOUBLE)) <= 1.7976931348623157e308)"
  ))
  stop("Finite SQL predicates support DuckDB and Spark only.", call. = FALSE)
}

# Spark sums doubles in a nondeterministic order; at millions of rows that
# noise exceeds the line-search deviance tolerance and the score tolerance
# (sum(x*score) cancels heavily). Fixed-point DECIMAL addition is exact, so
# these sums no longer depend on summation order.
sql_stable_sum = function(conn, expr, backend_info = detect_backend(conn)) {
  if (!identical(backend_info[["name"]], "spark")) return(sql_sum(conn, expr, backend_info))
  value = paste0("(", expr, ")")
  ok = paste0(sql_is_finite(conn, value, backend_info), " AND ABS(", value, ") < 1e10")
  # DECIMAL(38,18) stabilizes distributed cancellation but rounds very small
  # scores to zero. When every term is small, ordinary double aggregation
  # preserves their relative scale and avoids that absolute quantization floor.
  sprintf(paste0(
    "(CASE WHEN MAX(ABS(%2$s)) < 1e-9 THEN SUM(%2$s) ",
    "WHEN SUM(CASE WHEN %1$s THEN 0 ELSE 1 END) = 0 THEN ",
    "CAST(SUM(CAST(CASE WHEN %1$s THEN (%2$s) END AS DECIMAL(38,18))) AS DOUBLE) ELSE SUM(%2$s) END)"
  ), ok, expr)
}


