prediction_frame_types = function(frame, conn) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  vapply(frame, function(value) {
    if (inherits(value, "integer64")) "BIGINT" else
    if (inherits(value, "Date")) "DATE" else
    if (inherits(value, "POSIXct")) if (native) "TIMESTAMP" else "TIMESTAMPTZ" else
    if (is.factor(value) || is.character(value)) if (native) "STRING" else "VARCHAR" else
    if (is.logical(value)) "BOOLEAN" else
    if (is.integer(value)) "INTEGER" else
    if (is.numeric(value) && !is.complex(value)) "DOUBLE" else
      stop("Unsupported prediction input column type: ", class(value)[[1L]], ".", call. = FALSE)
  }, character(1))
}

prediction_frame_text = function(value) {
  if (inherits(value, "integer64")) return(as.character(value))
  if (inherits(value, "POSIXct")) return(format(value, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
  if (is.numeric(value) && !inherits(value, "Date")) {
    text = dbglm_number(value)
    text[is.na(value)] = NA_character_
    text[!is.na(value) & value == Inf] = "Infinity"
    text[!is.na(value) & value == -Inf] = "-Infinity"
    return(text)
  }
  as.character(value)
}

# Bounded incoming batches/portable dictionaries, never a collected SQL panel.
prediction_upload_frame = function(conn, frame, resources, sql_types = NULL) {
  if (is.null(sql_types)) sql_types = prediction_frame_types(frame, conn)
  stopifnot(length(sql_types) == ncol(frame), ncol(frame) > 0L)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  name = basename(tempfile("dbreg_prediction_input_"))
  resources[["tables"]] = c(resources[["tables"]], name)
  fields = sprintf("CAST(%s AS %s) AS %s", q(names(frame)), sql_types, q(names(frame)))
  empty = paste("SELECT", paste(sprintf("CAST(NULL AS %s) AS %s", sql_types, q(names(frame))), collapse = ","), "WHERE FALSE")
  DBI::dbExecute(conn, sql_sample_table(conn, name, empty, names(frame)))
  if (!nrow(frame)) return(q(name))
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  batch = max(1L, min(4096L, floor(prediction_state_limit() / ncol(frame))))
  first = 1L
  while (first <= nrow(frame)) {
    last = min(nrow(frame), first + batch - 1L)
    chunk = as.data.frame(lapply(frame[first:last, , drop = FALSE], prediction_frame_text), stringsAsFactors = FALSE,
      check.names = FALSE)
    if (native) {
      values = dbglm_values_sql(conn, chunk, max_bytes = 8 * 1024^2, format = "json")
      if (is.null(values)) {
        if (batch == 1L) stop("A prediction input row exceeds the SQL transfer limit.", call. = FALSE)
        batch = max(1L, floor(batch / 2L))
        next
      }
      source = paste0("(", values, ") incoming")
      DBI::dbExecute(conn, paste("INSERT INTO", q(name), "SELECT", paste(fields, collapse = ","), "FROM", source))
    } else {
      registered = basename(tempfile("dbreg_prediction_chunk_"))
      tryCatch({
        duckdb::duckdb_register(conn, registered, chunk)
        DBI::dbExecute(conn, paste("INSERT INTO", q(name), "SELECT", paste(fields, collapse = ","), "FROM", q(registered)))
      }, finally = try(duckdb::duckdb_unregister(conn, registered), silent = TRUE))
    }
    first = last + 1L
  }
  q(name)
}

# Keys are length-prefixed UTF-8 byte strings. Decode only a bounded portable
# dictionary already held in R; no SQL observation data is collected here.
prediction_decode_keys = function(dictionary) {
  width = length(dictionary[["columns"]])
  columns = stats::setNames(rep(list(character(length(dictionary[["key"]]))), width), paste0("k", seq_len(width)))
  for (i in seq_along(dictionary[["key"]])) {
    raw = charToRaw(enc2utf8(dictionary[["key"]][[i]]))
    position = 1L
    for (j in seq_len(width)) {
      remaining = seq.int(position, length(raw))
      colon = remaining[which(raw[remaining] == as.raw(58))[[1L]]]
      size = suppressWarnings(as.numeric(rawToChar(raw[seq.int(position, colon - 1L)])))
      if (!is.finite(size) || size < 0 || size != floor(size) || colon + size > length(raw)) {
        stop("Malformed portable prediction key.", call. = FALSE)
      }
      value = if (size) rawToChar(raw[seq.int(colon + 1L, colon + size)]) else ""
      Encoding(value) = "UTF-8"
      columns[[j]][[i]] = value
      position = colon + size + 1L
    }
    if (position != length(raw) + 1L) stop("Malformed portable prediction key.", call. = FALSE)
  }
  as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
}

prediction_integer_sql_type = function(values, conn) {
  negative = startsWith(values, "-")
  digits = sub("^-", "", values)
  long_limit = ifelse(negative, "9223372036854775808", "9223372036854775807")
  if (all(nchar(digits) < 19L | nchar(digits) == 19L & digits <= long_limit)) return("BIGINT")
  if (max(nchar(digits)) <= 38L) return("DECIMAL(38,0)")
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    stop("Portable integer keys exceed this backend's exact decimal range.", call. = FALSE)
  }
  if (any(negative)) "HUGEINT" else "UHUGEINT"
}
