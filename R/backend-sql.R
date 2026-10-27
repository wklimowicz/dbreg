# Translate only generated SQL tokens, preserving quoted source identifiers and
# literals. Finite predicates are emitted at their owning builders. This
# temporary adapter still owns exclusion and text-type token families.
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
  tokens[tokens == "EXCLUDE"] = "EXCEPT"
  tokens[tokens == "VARCHAR"] = "STRING"
  paste0(tokens, collapse = "")
}

# Text spelling stays at construction; generic dialects retain legacy VARCHAR.
sql_text_type = function(conn, backend_info = detect_backend(conn)) {
  if (identical(backend_info[["name"]], "spark")) "STRING" else "VARCHAR"
}

sql_cast_text = function(conn, expr, backend_info = detect_backend(conn)) {
  if (!length(expr)) return(character())
  paste0("CAST(", expr, " AS ", sql_text_type(conn, backend_info), ")")
}

# These feature paths already require native star exclusion. Never discover a
# schema inside a working-state iteration to emulate an unsupported dialect.
sql_select_except = function(conn, columns, backend_info = detect_backend(conn)) {
  if (!length(columns)) return("*")
  keyword = switch(backend_info[["name"]], duckdb = "EXCLUDE", spark = "EXCEPT",
    stop("Column exclusion supports DuckDB and Spark only.", call. = FALSE))
  paste0("* ", keyword, "(", paste(DBI::dbQuoteIdentifier(conn, columns), collapse = ","), ")")
}

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


backend_get_query = function(conn, sql) DBI::dbGetQuery(conn, backend_sql(conn, sql))

# Projection preflights share stable locations and UTF-8 byte bounds across
# DuckDB and native Spark, without fetching an observation-sized key column.
sql_projection_mean = function(conn, expr) {
  if (identical(detect_backend(conn)[["name"]], "duckdb")) return(paste0("FAVG(", expr, ")"))
  paste0("(", sql_stable_sum(conn, expr), ")/COUNT(*)")
}

sql_text_bytes = function(conn, expr) {
  backend_info = detect_backend(conn)
  value = sql_cast_text(conn, expr, backend_info)
  if (identical(backend_info[["name"]], "duckdb")) value = paste0("encode(", value, ")")
  paste0("OCTET_LENGTH(", value, ")")
}

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
  backend_info = detect_backend(conn)
  if (format == "json") {
    # Arrays avoid repeating field names for every FE. Seventeen digits preserve
    # double round trips (digits=NA rounds some extremes to infinity). Base64
    # keeps SQL string escaping independent of arbitrary key text/Unicode.
    payload = json_lookup_rows(frame)
    encoded = raw_to_base64(charToRaw(enc2utf8(payload)))
    if (nchar(encoded, type = "bytes") + 1024 > max_bytes) return(NULL)
    fields = sprintf(
      "CAST(item[%d] AS %s) AS %s", seq_along(frame) - 1L,
      ifelse(vapply(frame, is.numeric, logical(1)), "DOUBLE", sql_text_type(conn, backend_info)),
      DBI::dbQuoteIdentifier(conn, names(frame))
    )
    return(paste0(
      "SELECT ", paste(fields, collapse = ","),
      " FROM explode(from_json(", sql_cast_text(conn, paste0("unbase64('", encoded, "')"), backend_info), ", ",
      "'ARRAY<ARRAY<STRING>>')) AS lookup(item)"
    ))
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

# Native Spark string parsing differs from DBI's standard quote doubling and
# backslash handling. Base64 keeps arbitrary UTF-8 values out of that parser.
# Other backends retain their DBI literal representation and native value type.
sql_value_literals = function(conn, values) {
  if (!length(values)) return(character())
  # ODBC can quote logicals as 1/0, which Spark cannot use as booleans.
  if (is.logical(values)) {
    result = rep("NULL", length(values))
    result[!is.na(values)] = ifelse(values[!is.na(values)], "TRUE", "FALSE")
    return(result)
  }
  if (is.factor(values)) values = as.character(values)
  if (is.numeric(values) && !inherits(values, "integer64")) {
    valid = !is.na(values)
    if (any(!is.finite(values[valid]))) stop("SQL factor levels must be finite or missing.", call. = FALSE)
    result = rep("NULL", length(values))
    result[valid] = dbglm_number(values[valid])
    return(result)
  }
  backend_info = if (is.character(values)) detect_backend(conn) else NULL
  if (is.character(values) && isTRUE(backend_info[["is_databricks"]])) {
    return(vapply(values, function(value) {
      if (is.na(value)) return(sql_cast_text(conn, "NULL", backend_info))
      encoded = raw_to_base64(charToRaw(enc2utf8(value)))
      sql_cast_text(conn, paste0("unbase64('", encoded, "')"), backend_info)
    }, character(1), USE.NAMES = FALSE))
  }
  as.character(DBI::dbQuoteLiteral(conn, values))
}

sql_factor_literals = function(conn, levels, classes = character()) {
  result = lapply(names(levels), function(name) {
    value = sql_value_literals(conn, levels[[name]])
    # Spark's stacked category query returns exact integer text. Restore the
    # numeric SQL type explicitly so comparison never promotes BIGINT to double.
    if (is.character(levels[[name]])) {
      type = switch(unname(classes[name]),
        integer64 = "BIGINT",
        integer = "DOUBLE",
        numeric = "DOUBLE",
        logical = "BOOLEAN",
        Date = "DATE",
        NULL
      )
      if (!is.null(type)) value = paste0("CAST(", value, " AS ", type, ")")
    }
    value
  })
  stats::setNames(result, names(levels))
}

# Bounded SQL lookup payloads need only UTF-8 strings and numeric arrays.
raw_to_base64 = function(value) {
  n = length(value)
  if (!n) return("")
  padding = (3L - n %% 3L) %% 3L
  bytes = matrix(as.integer(c(value, raw(padding))), ncol = 3L, byrow = TRUE)
  words = bytes[, 1L] * 65536L + bytes[, 2L] * 256L + bytes[, 3L]
  alphabet = strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "", fixed = TRUE)[[1L]]
  indices = vapply(c(18L, 12L, 6L, 0L), function(shift) bitwAnd(bitwShiftR(words, shift), 63L) + 1L, integer(length(words)))
  output = alphabet[as.vector(t(matrix(indices, ncol = 4L)))]
  if (padding) output[length(output) - seq_len(padding) + 1L] = "="
  paste0(output, collapse = "")
}

json_string = function(value) {
  value = enc2utf8(value)
  value = gsub("\\", "\\\\", value, fixed = TRUE)
  value = gsub("\"", "\\\"", value, fixed = TRUE)
  for (code in seq_len(31L)) {
    value = gsub(intToUtf8(code), sprintf("\\u%04x", code), value, fixed = TRUE)
  }
  ifelse(is.na(value), "null", paste0("\"", value, "\""))
}

json_lookup_rows = function(frame) {
  fields = lapply(frame, function(value) {
    if (is.numeric(value) && !inherits(value, "integer64")) {
      ifelse(is.na(value), "null", sprintf("%.17g", value))
    } else if (is.logical(value)) {
      ifelse(is.na(value), "null", ifelse(value, "true", "false"))
    } else json_string(as.character(value))
  })
  rows = do.call(paste, c(fields, sep = ","))
  paste0("[", paste0("[", rows, "]", collapse = ","), "]")
}

# Spark generates flat numeric arrays/objects, never nested user JSON. Parse
# that deliberately narrow schema without adding a package dependency.
read_numeric_json = function(value) {
  value = trimws(value)
  object = startsWith(value, "{") && endsWith(value, "}")
  array = startsWith(value, "[") && endsWith(value, "]")
  if (!object && !array) stop("Invalid numeric SQL payload.", call. = FALSE)
  body = trimws(substr(value, 2L, nchar(value) - 1L))
  if (!nzchar(body)) return(if (object) list() else numeric())
  if (endsWith(body, ",")) stop("Invalid numeric SQL payload.", call. = FALSE)
  fields = trimws(strsplit(body, ",", fixed = TRUE)[[1L]])
  keys = NULL
  if (object) {
    valid = grepl("^\"[A-Za-z_][A-Za-z_0-9]*\"[[:space:]]*:", fields)
    if (!all(valid)) stop("Invalid numeric SQL payload keys.", call. = FALSE)
    keys = sub("^\"([A-Za-z_][A-Za-z_0-9]*)\"[[:space:]]*:.*$", "\\1", fields)
    fields = trimws(sub("^\"[A-Za-z_][A-Za-z_0-9]*\"[[:space:]]*:", "", fields))
  }
  valid = grepl("^(null|-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?)$", fields)
  if (!all(valid) || anyDuplicated(keys)) stop("Invalid numeric SQL payload values.", call. = FALSE)
  result = rep(NA_real_, length(fields))
  nonnull = fields != "null"
  result[nonnull] = as.numeric(fields[nonnull])
  if (object) stats::setNames(as.list(result), keys) else result
}
