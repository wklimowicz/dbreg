sql_prediction_values_table = function(conn, fields, values, resources, budget) {
  rows = paste0("(", do.call(paste, c(values, list(sep = ","))), ")")
  bytes = nchar(rows, type = "bytes") + 1
  budget[["used"]] = budget[["used"]] + sum(bytes) + sum(nchar(fields))
  if (budget[["used"]] > budget[["limit"]]) {
    projection_abort("transfer_limit", "Prediction input or dictionaries exceed the SQL upload budget.")
  }
  create = fe_projection_tables(conn, resources)
  aliases = as.character(DBI::dbQuoteIdentifier(conn, fields))
  table = NULL
  start = 1L
  while (start <= length(rows)) {
    remaining = seq.int(start, min(start + 511L, length(rows)))
    accepted = remaining[cumsum(bytes[remaining]) <= 128 * 1024]
    end = if (length(accepted)) max(accepted) else start
    query = paste0("SELECT * FROM VALUES ", paste(rows[seq.int(start, end)], collapse = ","),
      " AS uploaded(", paste(aliases, collapse = ","), ")")
    if (is.null(table)) table = create(query, fields)
    else DBI::dbExecute(conn, paste("INSERT INTO", table, query))
    start = end + 1L
  }
  table
}

sql_prediction_input_literals = function(conn, value) {
  if (is.factor(value)) value = as.character(value)
  kind = if (inherits(value, "integer64") || is.integer(value)) "BIGINT" else
    if (inherits(value, "Date")) "DATE" else if (inherits(value, "POSIXt")) "TIMESTAMP" else
      if (is.character(value)) "STRING" else if (is.logical(value)) "BOOLEAN" else
        if (is.double(value)) "DOUBLE" else NULL
  if (is.null(kind) || !is.null(dim(value))) {
    stop("Unsupported R column type for database prediction.", call. = FALSE)
  }
  present = !is.na(value)
  literal = rep(paste0("CAST(NULL AS ", kind, ")"), length(value))
  val = value[present]
  literal[present] = switch(kind,
    BIGINT = category_sql_literals(conn, as.character(val), "integer"),
    DATE = category_sql_literals(conn, val, "date"),
    TIMESTAMP = paste0("CAST(", category_string_literals(conn,
      format(val, "%Y-%m-%dT%H:%M:%OS6+00:00", tz = "UTC")), " AS TIMESTAMP)"),
    STRING = category_string_literals(conn, val),
    BOOLEAN = category_sql_literals(conn, val, "logical"),
    DOUBLE = ifelse(is.finite(val), sql_cast_double(sprintf("%.17e", val), conn),
      paste0("CAST('", ifelse(val < 0, "-", ""), "Infinity' AS DOUBLE)")))
  list(literals = literal, kind = kind)
}

sql_native_prediction_input = function(conn, data, needed, resources) {
  if (!is.data.frame(data)) stop("newdata must be a data frame.", call. = FALSE)
  if (!all(needed %in% names(data))) {
    stop("Missing prediction source columns: ", paste(setdiff(needed, names(data)), collapse = ", "), call. = FALSE)
  }
  limits = projection_limits()
  limit = min(limits[["upload_bytes"]], limits[["memory_bytes"]] / 8)
  if (as.numeric(utils::object.size(data[needed])) > 4 * limit) {
    projection_abort("transfer_limit", "Prediction input exceeds the SQL upload budget.")
  }
  literals = lapply(data[needed], function(value) sql_prediction_input_literals(conn, value))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  if (!nrow(data)) {
    select = paste0("CAST(NULL AS ", vapply(literals, `[[`, character(1), "kind"), ") AS ", q(needed))
    return(fe_projection_tables(conn, resources)(paste("SELECT", paste(select, collapse = ","), "WHERE FALSE"), needed))
  }
  values = lapply(literals, `[[`, "literals")
  budget = list2env(list(used = 0, limit = limit), parent = emptyenv())
  sql_prediction_values_table(conn, needed, values, resources, budget)
}
