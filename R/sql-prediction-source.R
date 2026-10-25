prediction_sql_kind = function(type) {
  type = tolower(gsub(" ", "", type, fixed = TRUE))
  if (grepl("^(u?(tinyint|smallint|int|integer|bigint|hugeint)|decimal\\([0-9]+,0\\))$", type)) "integer" else
  if (grepl("^(double|float|real|decimal)", type)) "numeric" else
  if (grepl("^(varchar|char|string|enum)", type)) "character" else
  if (type %in% c("boolean", "bool")) "logical" else
  if (type == "date") "date" else
  if (startsWith(type, "timestamp")) "datetime" else "unsupported"
}

prediction_factor_class = function(design, name) {
  class = unname(design[["factor_classes"]][name])
  if (!length(class) || is.na(class)) class = class(design[["factor_levels"]][[name]])[[1L]]
  class
}

sql_prediction_source = function(object, context, conn, from_statement, row_id, raw_required, dataframe_labels = NULL) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  dictionaries = context[["state"]][["specification"]][["dictionaries"]]
  key_sources = lapply(dictionaries, function(dictionary) {
    columns = dictionary[["columns"]]
    for (j in seq_along(columns)) {
      if (dictionary[["types"]][[j]] == "character" && !is.null(dataframe_labels)) columns[[j]] = dataframe_labels[[columns[[j]]]]
    }
    columns
  })
  needed = unique(c(raw_required, unname(context[["design"]][["sources"]]), unlist(key_sources), row_id))
  schema = DBI::dbGetQuery(conn, paste("SELECT *", from_statement, "LIMIT 0"))
  absent = setdiff(needed, names(schema))
  if (length(absent)) stop("Prediction source is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  type_row = DBI::dbGetQuery(conn, paste("SELECT", paste0("typeof(", q(needed), ") AS t", seq_along(needed), collapse = ","), from_statement, "LIMIT 1"))
  types = if (nrow(type_row)) stats::setNames(vapply(type_row, function(value) value[[1L]], character(1)), needed) else
    prediction_frame_types(schema[needed], conn)
  kinds = vapply(types, prediction_sql_kind, character(1))
  classes = vapply(schema[needed], function(value) class(value)[[1L]], character(1))
  reserved = c(names(context[["design"]][["sources"]]), "__dbglm_y", "__dbglm_a", "__dbglm_o", "__prediction_eta")
  raw_aliases = safe_identifiers(paste0("dbreg_prediction_raw", seq_along(needed), "__"), taken = reserved,
    fold_case = folds_identifier_case(conn))
  names(raw_aliases) = needed
  key_aliases = safe_identifiers(paste0("dbreg_prediction_id", seq_along(row_id), "__"), taken = c(reserved, raw_aliases), fold_case = folds_identifier_case(conn))
  x = if (length(context[["beta"]])) safe_identifiers(paste0("dbreg_prediction_x", seq_along(context[["beta"]]), "__"),
    taken = c(reserved, raw_aliases, key_aliases), fold_case = folds_identifier_case(conn)) else character()
  varying = context[["varying"]]
  varying_aliases = if (length(varying)) safe_identifiers(paste0("dbreg_prediction_v", seq_along(varying), "__"),
    taken = c(reserved, raw_aliases, key_aliases, x), fold_case = folds_identifier_case(conn)) else character()
  names(varying_aliases) = varying
  offset_alias = safe_identifiers("dbreg_prediction_offset__", taken = c(reserved, raw_aliases, key_aliases, x, varying_aliases), fold_case = folds_identifier_case(conn))
  raw = paste("SELECT", paste(q(needed), "AS", raw_aliases, collapse = ","), from_statement)
  precision = character()
  source_terms = vapply(names(context[["design"]][["sources"]]), function(name) {
    source = context[["design"]][["sources"]][[name]]
    value = paste0("raw.", raw_aliases[[source]])
    if (name %in% context[["design"]][["numeric_sources"]]) {
      value = sql_prediction_numeric(conn, value, classes[[source]])
    } else {
      class = prediction_factor_class(context[["design"]], name)
      kind = kinds[[source]]
      compatible = switch(class, character = kind == "character", factor = kind == "character",
        ordered = kind == "character", logical = kind == "logical", Date = kind == "date",
        integer64 = kind %in% c("integer", "numeric"), integer = kind %in% c("integer", "numeric"),
        numeric = kind %in% c("integer", "numeric"), FALSE)
      if (!compatible) stop("Prediction factor `", source, "` requires its fitted column type.", call. = FALSE)
      if (class == "integer64") value = paste0("TRY_CAST(", value, " AS BIGINT)") else
        if (class %in% c("numeric", "integer")) value = paste0("CAST(", value, " AS DOUBLE)")
    }
    paste(value, "AS", q(name))
  }, character(1))
  # Precision checks use original fields, before a categorical numeric cast.
  for (name in names(context[["design"]][["factor_levels"]])) {
    source = context[["design"]][["sources"]][[name]]
    class = prediction_factor_class(context[["design"]], name)
    value = paste0("r.", raw_aliases[[source]])
    if (class == "integer64" && kinds[[source]] == "numeric") precision = c(precision,
      paste0("(ABS(", value, ")>9007199254740992 OR ", value, "<>FLOOR(", value, "))"))
    if (class %in% c("numeric", "integer") && kinds[[source]] == "integer") precision = c(precision,
      paste0("(", value, ">9007199254740992 OR ", value, "< -9007199254740992)"))
  }
  for (i in seq_along(dictionaries)) for (j in seq_along(key_sources[[i]])) {
    source = key_sources[[i]][[j]]
    type = dictionaries[[i]][["types"]][[j]]
    kind = kinds[[source]]
    compatible = switch(type, integer = kind %in% c("integer", "numeric", "character"),
      numeric = kind %in% c("integer", "numeric"), float = kind %in% c("integer", "numeric"),
      character = kind == "character", logical = kind == "logical", date = kind %in% c("date", "character"), FALSE)
    if (!compatible) stop("Prediction FE key `", source, "` requires its fitted column type.", call. = FALSE)
    value = paste0("r.", raw_aliases[[source]])
    if (type == "integer" && kind == "numeric") precision = c(precision,
      paste0("(ABS(", value, ")>9007199254740992 OR ", value, "<>FLOOR(", value, "))"))
    if (type %in% c("numeric", "float") && kind == "integer") precision = c(precision,
      paste0("(", value, ">9007199254740992 OR ", value, "< -9007199254740992)"))
  }
  output_keys = vapply(seq_along(row_id), function(i) {
    value = paste0("raw.", raw_aliases[[row_id[[i]]]])
    # SQL NaN is non-NULL on DuckDB/Spark, but is a missing key in R.
    if (kinds[[row_id[[i]]]] == "numeric") value = paste0("CASE WHEN isnan(CAST(", value,
      " AS DOUBLE)) THEN NULL ELSE ", value, " END")
    paste(value, "AS", key_aliases[[i]])
  }, character(1))
  extras = c(output_keys,
    vapply(varying, function(name) paste(sql_prediction_numeric(conn, paste0("raw.", raw_aliases[[name]]), classes[[name]]), "AS", varying_aliases[[name]]), character(1)),
    paste(if (is.null(context[["offset"]])) "CAST(0.0 AS DOUBLE)" else
      sql_prediction_numeric(conn, paste0("raw.", raw_aliases[[context[["offset"]]]]), classes[[context[["offset"]]]]), "AS", offset_alias))
  prepared = paste0("(SELECT raw.*,", paste(c(source_terms, extras), collapse = ","), " FROM (", raw, ") raw) r")
  fields = c("r.*", sprintf("(CAST((%s) AS DOUBLE)-(%s)) AS %s", context[["design"]][["expressions"]], dbglm_number(context[["center"]]), x))
  frame = paste0("(SELECT ", paste(fields, collapse = ","), " FROM ", prepared, ") r")
  raw_values = if (length(raw_required)) paste0("r.", raw_aliases[raw_required]) else character()
  valid = if (length(raw_values)) paste0(raw_values, " IS NOT NULL") else character()
  numeric = raw_required[classes[raw_required] %in% c("numeric", "integer", "integer64")]
  if (length(numeric)) {
    values = paste0("r.", raw_aliases[numeric])
    valid = c(valid, paste0("isfinite(", values, ")"))
  }
  list(frame = frame, x = x, varying = varying_aliases, offset = offset_alias,
    key_aliases = key_aliases, raw_aliases = raw_aliases, key_sources = key_sources,
    types = types, classes = classes, valid = if (length(valid)) paste(valid, collapse = " AND ") else "TRUE",
    precision = if (length(precision)) paste(precision, collapse = " OR ") else NULL)
}
