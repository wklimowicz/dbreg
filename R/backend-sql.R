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
sql_stable_sum = function(conn, expr, backend_info = detect_backend(conn), window = NULL, max_abs = NULL) {
  over = if (is.null(window)) "" else paste0(" OVER (", window, ")")
  if (!identical(backend_info[["name"]], "spark")) return(paste0(sql_sum(conn, expr, backend_info), over))
  aggregate = function(name, value) paste0(name, "(", value, ")", over)
  value = paste0("(", expr, ")")
  ok = paste0(sql_is_finite(conn, value, backend_info), " AND ABS(", value, ") < 1e10")
  # Preserve very small scores rather than quantizing them to decimal zero.
  ordinary = aggregate("SUM", expr)
  # Only callers that verified complete finite rows may supply this bound.
  # The tenfold margin also avoids relying on rounded bounds near the cast cap.
  if (is.numeric(max_abs) && length(max_abs) == 1L && is.finite(max_abs) &&
      max_abs >= 0 && max_abs < 1e9) {
    return(paste0("(CASE WHEN ", aggregate("MAX", paste0("ABS(", expr, ")")),
      " < 1e-9 THEN ", ordinary, " ELSE CAST(",
      aggregate("SUM", paste0("CAST(", value, " AS DECIMAL(38,18))")),
      " AS DOUBLE) END)"))
  }
  paste0("(CASE WHEN ", aggregate("MAX", paste0("ABS(", expr, ")")),
    " < 1e-9 THEN ", ordinary, " WHEN ",
    aggregate("SUM", paste0("CASE WHEN ", ok, " THEN 0 ELSE 1 END")),
    " = 0 THEN CAST(", aggregate("SUM", paste0("CAST(CASE WHEN ", ok,
      " THEN (", expr, ") END AS DECIMAL(38,18))")),
    " AS DOUBLE) ELSE ", ordinary, " END)")
}


# These feature paths already require native star exclusion. Never discover a
# schema inside a working-state iteration to emulate an unsupported dialect.
sql_select_except = function(conn, columns, backend_info = detect_backend(conn)) {
  if (!length(columns)) return("*")
  keyword = switch(backend_info[["name"]], duckdb = "EXCLUDE", spark = "EXCEPT",
    stop("Column exclusion supports DuckDB and Spark only.", call. = FALSE))
  paste0("* ", keyword, "(", paste(DBI::dbQuoteIdentifier(conn, columns), collapse = ","), ")")
}
