# Native Spark string parsing differs from DBI's standard quote doubling and
# backslash handling. Base64 keeps arbitrary UTF-8 values out of that parser.
# Other backends retain their DBI literal representation and native value type.
sql_value_literals = function(conn, values) {
  if (!length(values)) return(character())
  # ODBC can quote logicals as 1/0, which Spark cannot use as booleans.
  if (is.logical(values)) {
    result = rep("NULL",length(values))
    result[!is.na(values)] = ifelse(values[!is.na(values)],"TRUE","FALSE")
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
  if (is.character(values) && isTRUE(detect_backend(conn)[["is_databricks"]])) {
    return(vapply(values, function(value) {
      if (is.na(value)) return("CAST(NULL AS STRING)")
      encoded = gsub("[\r\n]", "", jsonlite::base64_enc(charToRaw(enc2utf8(value))))
      paste0("CAST(unbase64('", encoded, "') AS STRING)")
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
      type = switch(unname(classes[name]),integer64="BIGINT",integer="DOUBLE",numeric="DOUBLE",logical="BOOLEAN",Date="DATE",NULL)
      if (!is.null(type)) value = paste0("CAST(", value, " AS ",type,")")
    }
    value
  })
  stats::setNames(result, names(levels))
}
