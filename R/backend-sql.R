# Translate only generated SQL tokens, preserving quoted source identifiers and
# literals. DuckDB retains its compensated aggregates; Spark uses double SUM.
backend_sql = function(conn, sql) {
  if (!identical(detect_backend(conn)[["name"]], "spark")) return(sql)
  # Consume quoted runs in chunks, without backtracking over each character.
  # FE coefficient lookups can contain megabytes of base64 in one literal;
  # character-by-character alternatives exhaust PCRE's match limit there.
  pattern = "'(?:''|[^']++)*+'|\"(?:\"\"|[^\"]++)*+\"|`(?:``|[^`]++)*+`|[A-Za-z_][A-Za-z_0-9]*|[[:space:]]+|."
  tokens = regmatches(sql, gregexpr(pattern, sql, perl = TRUE))[[1]]
  if (sum(nchar(tokens, type = "bytes")) != nchar(sql, type = "bytes")) {
    stop("Could not tokenize the complete generated GLM SQL statement.", call. = FALSE)
  }
  tokens[tokens == "FSUM"] = "SUM"
  tokens[tokens == "EXCLUDE"] = "EXCEPT"
  tokens[tokens == "VARCHAR"] = "STRING"
  repeat {
    positions = which(tokens == "isfinite")
    if (!length(positions)) break
    start = positions[1]
    open = start + 1L
    while (grepl("^[[:space:]]+$", tokens[open])) open = open + 1L
    stopifnot(tokens[open] == "(")
    depth = 1L
    end = open
    while (depth > 0L) {
      end = end + 1L
      if (tokens[end] == "(") depth = depth + 1L
      if (tokens[end] == ")") depth = depth - 1L
    }
    arg = paste0(tokens[seq.int(open + 1L, end - 1L)], collapse = "")
    value = paste0("(NOT isnan(CAST((", arg, ") AS DOUBLE)) AND ABS(CAST((", arg, ") AS DOUBLE)) <= 1.7976931348623157e308)")
    tokens = c(if (start > 1L) tokens[seq_len(start - 1L)], value, if (end < length(tokens)) tokens[seq.int(end + 1L, length(tokens))])
  }
  paste0(tokens, collapse = "")
}

# Spark sums doubles in a nondeterministic order; at millions of rows that
# noise exceeds the line-search deviance tolerance and the score tolerance
# (sum(x*score) cancels heavily). Fixed-point DECIMAL addition is exact, so
# these sums no longer depend on summation order.
sql_stable_sum = function(conn, expr) {
  if (!identical(detect_backend(conn)[["name"]], "spark")) return(sprintf("FSUM(%s)", expr))
  ok = sprintf("isfinite(%1$s) AND ABS(%1$s) < 1e10", paste0("(", expr, ")"))
  # DECIMAL(38,18) stabilizes distributed cancellation but rounds very small
  # scores to zero. When every term is small, ordinary double aggregation
  # preserves their relative scale and avoids that absolute quantization floor.
  sprintf(paste0("(CASE WHEN MAX(ABS(%2$s)) < 1e-9 THEN SUM(%2$s) ",
    "WHEN SUM(CASE WHEN %1$s THEN 0 ELSE 1 END) = 0 THEN ",
    "CAST(SUM(CAST(CASE WHEN %1$s THEN (%2$s) END AS DECIMAL(38,18))) AS DOUBLE) ELSE SUM(%2$s) END)"), ok, expr)
}


backend_get_query = function(conn, sql) DBI::dbGetQuery(conn, backend_sql(conn, sql))

# Bounded lookup literals hold FE coefficients or requested prediction keys,
# never training observation state. A size limit bounds SQL text. The caller
# owns materialization and cleanup just as for every other projector state.
sql_values_frame = function(conn, frame, max_bytes = 8 * 1024^2,
                             format = getOption("dbreg.glm_databricks_lookup", "json")) {
  if (!nrow(frame)) stop("A coefficient lookup cannot be empty.", call. = FALSE)
  if (!is.character(format) || length(format) != 1 || !format %in% c("json", "values")) {
    stop("GLM FE lookup format must be json or values.", call. = FALSE)
  }
  if (any(vapply(frame, function(x) is.numeric(x) && any(!is.finite(x)), logical(1)))) {
    stop("Nonfinite projection coefficients.", call. = FALSE)
  }
  if (format == "json") {
    # Arrays avoid repeating field names for every FE. Seventeen digits preserve
    # double round trips (digits=NA rounds some extremes to infinity). Base64
    # keeps SQL string escaping independent of arbitrary key text/Unicode.
    payload = jsonlite::toJSON(frame, dataframe = "values", digits = 17, na = "null")
    encoded = gsub("[\r\n]", "", jsonlite::base64_enc(charToRaw(enc2utf8(payload))))
    if (nchar(encoded, type = "bytes") + 1024 > max_bytes) return(NULL)
    fields = sprintf("CAST(item[%d] AS %s) AS %s", seq_along(frame) - 1L,
      ifelse(vapply(frame, is.numeric, logical(1)), "DOUBLE", "STRING"),
      DBI::dbQuoteIdentifier(conn, names(frame)))
    return(paste0("SELECT ", paste(fields, collapse = ","),
      " FROM explode(from_json(CAST(unbase64('", encoded, "') AS STRING), ",
      "'ARRAY<ARRAY<STRING>>')) AS lookup(item)"))
  }
  fields = lapply(frame, function(x) {
    if (is.numeric(x)) {
      if (any(!is.finite(x))) stop("Nonfinite projection coefficients.", call. = FALSE)
      ifelse(x == 0, "0.0", sprintf("%.17ge0", x))
    } else as.character(DBI::dbQuoteString(conn, as.character(x)))
  })
  # sprintf can already emit an exponent; append e0 only to decimal literals.
  for (j in which(vapply(frame, is.numeric, logical(1)))) {
    fields[[j]] = sub("(e[+-]?[0-9]+)e0$", "\\1", fields[[j]])
  }
  rows = do.call(paste, c(fields, sep = ","))
  if (sum(nchar(rows, type = "bytes")) + 3 * nrow(frame) > max_bytes) return(NULL)
  names = paste(DBI::dbQuoteIdentifier(conn, names(frame)), collapse = ",")
  paste0("SELECT * FROM VALUES ", paste0("(", rows, ")", collapse = ","), " AS lookup(", names, ")")
}
