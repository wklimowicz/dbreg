# Bounded typed category dictionaries and selector semantics for SQL i().
read_i_schema = function(conn, query) {
  result = DBI::dbSendQuery(conn, query)
  on.exit(DBI::dbClearResult(result))
  columns = DBI::dbColumnInfo(result)
  schema = DBI::dbFetch(result, n = 0)
  driver = if (inherits(result, "OdbcResult")) "odbc" else
    if (inherits(result, "MySQLResult")) "rmysql" else "unknown"
  list(schema = schema, columns = columns, driver = driver)
}

# Native type evidence classifies i() keys independently of R transport.
# Literal construction remains separate from this precision/type probe.
i_native_metadata = function(conn, source, column, schema, backend_info = detect_backend(conn)) {
  q = as.character(DBI::dbQuoteIdentifier(conn, column))
  if (is.null(schema)) {
    schema = read_i_schema(conn, sql_limit(paste("SELECT", q, "FROM", source), 0, backend_info[["name"]]))
  }
  unavailable = function(detail) stop("Cannot preserve i() key `", column,
    "` exactly: ", detail, ".", call. = FALSE)
  info = schema[["columns"]]
  if (!is.data.frame(info) || !all(c("name", "type") %in% names(info)) ||
      sum(info[["name"]] == column, na.rm = TRUE) != 1L) unavailable("native column metadata is unavailable")
  value = as.character(info[["type"]][match(column, info[["name"]])])
  if (length(value) != 1L || is.na(value)) unavailable("native column type is unavailable")
  backend = backend_info[["name"]]
  driver = schema[["driver"]]
  if (!is.character(driver) || length(driver) != 1L || is.na(driver)) driver = "unknown"
  native_type = NULL
  if (driver == "odbc" && backend %in% c("sqlserver", "mysql")) {
    # SQL type constants exposed by odbc::dbColumnInfo(OdbcResult).
    native_type = c("-6" = "TINYINT", "5" = "SMALLINT", "4" = "INTEGER", "-5" = "BIGINT",
      "8" = "DOUBLE", "7" = "REAL", "6" = "FLOAT", "2" = "NUMERIC", "3" = "DECIMAL",
      "-7" = "BIT", "91" = "DATE", "1" = "CHAR", "12" = "VARCHAR",
      "-1" = "VARCHAR", "-8" = "CHAR", "-9" = "VARCHAR", "-10" = "VARCHAR")[value]
  } else if (driver == "rmysql" && backend == "mysql") {
    # RMySQL's public native labels omit NEWDECIMAL on some versions. Unknown
    # labels stay unknown: the accompanying Sclass is not native type evidence.
    labels = c("TINYINT", "SMALLINT", "MEDIUMINT", "INTEGER", "BIGINT", "DOUBLE", "FLOAT",
      "DECIMAL", "DATE", "ENUM", "SET", "VAR_STRING", "STRING")
    if (value %in% labels) native_type = value
  }
  if (is.null(native_type) || is.na(native_type)) unavailable(paste0(
    "unsupported native metadata from ", driver, " (", value, ")"))
  native_type = unname(native_type)
  type = if (native_type %in% c("TINYINT", "SMALLINT", "MEDIUMINT", "INTEGER", "BIGINT")) "integer" else
    if (native_type == "DOUBLE") "numeric" else
    if (native_type %in% c("CHAR", "VARCHAR", "ENUM", "SET", "VAR_STRING", "STRING")) "character" else
    if (native_type == "DATE") "date" else NA_character_
  if (native_type == "BIT" && backend == "sqlserver") type = "logical"
  if (backend == "sqlserver" && native_type %in% c("FLOAT", "NUMERIC", "DECIMAL")) {
    # Only known numeric native families are variant-compatible. In particular,
    # never apply this probe to (max) strings, XML, LOBs or mixed sql_variant.
    property = function(name, cast) paste0("CAST(SQL_VARIANT_PROPERTY(CAST(", q,
      " AS sql_variant), '", name, "') AS ", cast, ")")
    query = paste0("SELECT ", property("BaseType", "VARCHAR(32)"), " AS base_type, ",
      property("Precision", "INT"), " AS precision, ", property("Scale", "INT"), " AS scale",
      " FROM (SELECT ", q, " FROM ", source, ") AS dbreg_key_source WHERE ", q, " IS NOT NULL")
    metadata = DBI::dbGetQuery(conn, sql_limit(query, 1, backend))
    if (nrow(metadata) != 1L || !all(c("base_type", "precision", "scale") %in% names(metadata))) {
      unavailable("native numeric precision/scale is unavailable")
    }
    base = tolower(as.character(metadata[["base_type"]]))
    precision = metadata[["precision"]]
    scale = metadata[["scale"]]
    if (anyNA(c(base, precision, scale)) || !is.numeric(precision) || !is.numeric(scale) ||
        precision != floor(precision) || scale != floor(scale)) unavailable("invalid native numeric precision/scale")
    if (native_type == "FLOAT" && base %in% c("float", "real") &&
        precision %in% c(24, 53) && scale == 0) {
      type = if (precision == 53) "numeric" else "float"
      native_type = if (precision == 53) "DOUBLE" else "REAL"
    } else if (native_type %in% c("NUMERIC", "DECIMAL") && base %in% c("numeric", "decimal") &&
               precision >= 1 && precision <= 38 && scale >= 0 && scale <= precision) {
      if (scale != 0) unavailable("fractional DECIMAL keys are unsupported")
      type = "integer"
      native_type = paste0("DECIMAL(", precision, ",0)")
    } else unavailable("inconsistent native numeric precision/scale")
  } else if (backend == "mysql" && native_type %in% c("NUMERIC", "DECIMAL")) {
    # MySQL's native decimal CHAR representation preserves declared fractional
    # width, even for zero/trailing zeros. This does not recover precision.
    server = DBI::dbGetInfo(conn)[["dbms.name"]]
    if (!length(server) || is.na(server[[1L]]) || !grepl("mysql", server[[1L]], ignore.case = TRUE) ||
        grepl("mariadb", server[[1L]], ignore.case = TRUE)) {
      unavailable("this server's native DECIMAL text-scale contract is unverified")
    }
    query = paste0("SELECT CAST(", q, " AS CHAR) AS native_decimal_text FROM (SELECT ",
      q, " FROM ", source, ") AS dbreg_key_source WHERE ", q, " IS NOT NULL")
    text = DBI::dbGetQuery(conn, sql_limit(query, 1, backend))
    text = if (ncol(text) == 1L) text[[1L]] else NULL
    if (!is.character(text) || length(text) != 1L || is.na(text) ||
        !grepl("^[+-]?[0-9]+(\\.[0-9]+)?$", text)) unavailable("native DECIMAL text scale is unavailable")
    if (grepl(".", text, fixed = TRUE)) unavailable("fractional DECIMAL keys are unsupported")
    type = "integer"
  } else if (native_type %in% c("FLOAT", "REAL")) type = "float"
  if (is.na(type)) unavailable(paste0("unsupported native key type ", native_type))
  list(types = type, sql_types = native_type)
}

#' Parse a `fixest`-style `i()` term
#'
#' @param term A term label.
#' @param env Environment in which to evaluate `ref` and `keep` (e.g. the
#'   formula's), or `NULL` to only validate the term.
#' @return `NULL` if `term` is not an `i()` call, else a list with the factor
#'   `f`, interacted variable `var` and the `ref` and `keep` levels (each
#'   `NULL` if absent or not evaluated).
#' @keywords internal
parse_i_term = function(term, env = NULL) {
  e = if (is.call(term)) term else tryCatch(str2lang(term), error = function(e) NULL)
  if (!is.call(e)) return(NULL)
  head = e[[1L]]
  qualified = is.call(head) && identical(head[[1L]], as.name("::")) &&
    as.character(head[[2L]]) %in% c("dbreg", "fixest") && identical(head[[3L]], as.name("i"))
  if (!identical(head, as.name("i")) && !qualified) return(NULL)
  e[[1L]] = as.name("i")
  label = paste(deparse(e, width.cutoff = 500L), collapse = "")
  e = tryCatch(match.call(i_template, e), error = function(err) {
    stop("Invalid i() term: ", label, " (", conditionMessage(err), ")", call. = FALSE)
  })
  unsupported = intersect(names(e), c("bin", "ref2", "keep2", "bin2"))
  if (length(unsupported)) {
    stop("i() arguments not supported yet: ", paste(unsupported, collapse = ", "), call. = FALSE)
  }
  # An expression in the second position denotes reference levels; a symbol
  # denotes a numeric/logical source column. Named ref/keep use the formula env.
  if (!is.null(e[["var"]]) && !is.symbol(e[["var"]])) {
    if ("ref" %in% names(e)) stop("Specify positional or named i() references, not both.", call. = FALSE)
    e[["ref"]] = e[["var"]]
    e[["var"]] = NULL
  }
  args = list(e[["factor_var"]], e[["var"]])
  if (any(vapply(args, function(x) is.symbol(x) && startsWith(as.character(x), "i."), logical(1)))) {
    stop("The i. prefix in i() terms is not supported yet: ", label, call. = FALSE)
  }
  if (!is.symbol(args[[1L]]) || !(is.null(args[[2L]]) || is.symbol(args[[2L]]))) {
    stop("i() terms must take plain variable names, e.g. i(period, treat, ref = 5): ", label, call. = FALSE)
  }
  value = function(arg) {
    expression = e[[arg]]
    if (is.null(expression) || is.null(env)) return(NULL)
    result = tryCatch(eval(expression, env), error = function(err) {
      stop("Could not evaluate `", arg, "` in ", label, ": ", conditionMessage(err), call. = FALSE)
    })
    if (!is.atomic(result) || anyNA(result) ||
      (is.numeric(result) && any(!is.finite(result)))) {
      stop("i() references/keep values must be finite, nonmissing vectors.", call. = FALSE)
    }
    result
  }
  ref = value("ref")
  keep = value("keep")
  if (!is.null(env)) {
    # Freeze evaluated selectors before compacting the formula environment.
    for (arg in intersect(c("ref", "keep"), names(e))) {
      e[arg] = list(if (arg == "ref") ref else keep)
    }
  }
  first = as.character(e[["factor_var"]])
  second = if (!is.null(e[["var"]])) as.character(e[["var"]])
  list(
    f = first, var = second, first = first, second = second,
    ref = ref, keep = keep, explicit = any(c("ref", "keep") %in% names(e)),
    call = e, term = paste(deparse(e, width.cutoff = 500L), collapse = "")
  )
}

# Selection is shared by native SQL expansion and retained typed dictionaries.
i_level_indices = function(value, keys, labels = as.character(keys), allow_partial = is.character(value)) {
  if (is.character(value)) {
    found = integer()
    for (item in value) {
      indices = if (allow_partial && startsWith(item, "@!")) {
        which(!grepl(substring(item, 3L), labels))
      } else if (allow_partial && startsWith(item, "@")) {
        grep(substring(item, 2L), labels)
      } else {
        exact = match(item, labels)
        if (is.na(exact) && allow_partial) pmatch(tolower(item), tolower(labels), duplicates.ok = TRUE) else exact
      }
      found = union(found, indices[!is.na(indices)])
    }
  } else {
    found = unique(match(value, keys))
    found = found[!is.na(found)]
  }
  if (!length(found)) stop("The i() reference/keep value(s) match no level.", call. = FALSE)
  found
}

i_column_limit = function() {
  value = getOption("dbreg.i_max_columns", 1000L)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) stop("dbreg.i_max_columns must be a positive integer.", call. = FALSE)
  value
}

i_selected_items = function(value, keys, labels, type) {
  allow_partial = is.character(value)
  if (is.logical(value) && type %in% c("integer", "numeric")) value = as.numeric(value)
  if (is.numeric(value) && type == "logical") value = as.logical(value)
  if (type == "integer" && !is.character(value)) value = sub("^[0-9]+:", "", category_key_part(value, "integer"))
  i_level_indices(value, keys, labels, allow_partial = allow_partial)
}

make_i_spec = function(declaration, keys, type, second_class = NULL, intercept = TRUE) {
  labels = as.character(keys)
  if (type == "numeric" && anyDuplicated(labels)) labels = sprintf("%.17g", keys)
  encoded = category_key_part(keys, type)
  if (anyNA(encoded) || anyDuplicated(encoded) || anyDuplicated(labels)) stop("Could not preserve i() levels exactly.", call. = FALSE)
  selected = seq_along(keys)
  if (!is.null(declaration[["second"]]) && !second_class %in% c("numeric", "integer", "integer64", "logical")) stop("The second i() variable must be numeric or logical; categorical interactions use ordinary factor syntax.", call. = FALSE)
  ref = declaration[["ref"]]
  if (isTRUE(ref) || !declaration[["explicit"]] && intercept && is.null(declaration[["second"]])) {
    selected = setdiff(selected, 1L)
  } else if (!is.null(ref)) selected = setdiff(selected, i_selected_items(ref, keys, labels, type))
  if (!is.null(declaration[["keep"]])) selected = intersect(selected, i_selected_items(declaration[["keep"]], keys, labels, type))
  names = if (length(selected)) paste0(
    declaration[["first"]], "::", labels[selected],
    if (!is.null(declaration[["second"]])) paste0(":", declaration[["second"]])
  ) else character()
  list(
    first = declaration[["first"]], second = declaration[["second"]], second_class = second_class,
    values = keys, keys = encoded, type = type, labels = labels, selected = selected, names = names, term = declaration[["term"]]
  )
}

# Only bounded, typed dictionaries cross the SQL/R boundary. Integer levels are
# transferred as decimal strings, never through the driver's double conversion.
i_dictionaries = function(conn, from_sql, declarations, classes, intercept, native_schema = NULL) {
  result = list()
  if (!length(declarations)) return(result)
  limit = i_column_limit()
  for (declaration in declarations) {
    first = declaration[["first"]]
    dictionary = read_category_dictionary(conn, from_sql, first, native_schema, limit)
    keys = dictionary[["values"]]
    type = dictionary[["type"]]
    second = declaration[["second"]]
    spec = make_i_spec(declaration, keys, type, if (!is.null(second)) classes[[second]], intercept || length(result) > 0L)
    key = make.unique(c(names(result), first))[[length(result) + 1L]]
    result[[key]] = spec
    if (sum(vapply(result, function(x) length(x[["selected"]]), integer(1))) > limit) stop("i() columns exceed dbreg.i_max_columns.", call. = FALSE)
  }
  result
}

category_key_part = function(value, type) {
  if (type == "integer") {
    if (is.numeric(value) && !inherits(value, "integer64")) {
      bad = !is.na(value) & (!is.finite(value) | value != floor(value) | abs(value) > 2^53)
      if (any(bad)) stop("Integer category keys must be exact; use integer64 or decimal strings for large keys.", call. = FALSE)
      value = format(value, scientific = FALSE, trim = TRUE, digits = 22)
    } else value = as.character(value)
    valid = !is.na(value) & grepl("^[+-]?[0-9]+$", value)
    value[!valid] = NA_character_
    value = sub("^\\+", "", value)
    value = sub("^(-?)0+([0-9])", "\\1\\2", value)
    value[value %in% c("-0", "+0")] = "0"
  } else if (type %in% c("numeric", "float")) {
    if (!is.numeric(value) || inherits(value, "integer64")) {
      stop("Numeric category keys require numeric newdata columns.", call. = FALSE)
    }
    if (type == "float") value = readBin(writeBin(as.double(value), raw(), size = 4L), "double", n = length(value), size = 4L)
    value[value == 0 & !is.na(value)] = 0
    value = ifelse(is.finite(value), sprintf("%.17g", value), NA_character_)
  } else if (type == "logical") {
    if (!is.logical(value)) stop("Logical category keys require logical newdata columns.", call. = FALSE)
    value = ifelse(is.na(value), NA_character_, ifelse(value, "TRUE", "FALSE"))
  } else if (type == "date") {
    if (!inherits(value, "Date") && !is.character(value)) stop("Date category keys require Date or ISO date-string columns.", call. = FALSE)
    value = as.character(value)
  } else value = as.character(value)
  ifelse(is.na(value), NA_character_, paste0(nchar(enc2utf8(value), type = "bytes"), ":", enc2utf8(value)))
}

category_key_type = function(type, backend = NULL) {
  type = tolower(type)
  if (identical(backend, "sqlite") && type == "real") return("numeric")
  if (grepl("^(tinyint|smallint|int|integer|bigint|hugeint|u(tinyint|smallint|integer|bigint|hugeint)|decimal\\([0-9]+,0\\))$", type)) "integer" else if (type %in% c("double", "double precision")) "numeric" else if (type %in% c("float", "real")) "float" else if (type %in% c("boolean", "bool")) "logical" else if (grepl("^(varchar|char|string|text|enum)", type)) "character" else if (type == "date") "date" else NA_character_
}

category_key_metadata = function(conn, retained, components) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  backend = detect_backend(conn)[["name"]]
  if (identical(backend, "sqlite")) {
    sql_types = vapply(components, function(component) {
      col = q(component)
      row = DBI::dbGetQuery(conn, paste0("SELECT typeof(", col, ") AS t FROM (SELECT ", col, " FROM ", retained, ") dbreg_meta_src WHERE ", col, " IS NOT NULL LIMIT 1"))
      if (nrow(row) == 1L) return(as.character(row[["t"]][[1L]]))
      fallback = DBI::dbGetQuery(conn, paste0("SELECT typeof(", col, ") AS t FROM ", retained, " LIMIT 1"))
      if (nrow(fallback) != 1L) stop("Cannot discover i() levels from an empty source.", call. = FALSE)
      as.character(fallback[["t"]][[1L]])
    }, character(1), USE.NAMES = FALSE)
    types = vapply(sql_types, function(type) category_key_type(type, backend), character(1))
  } else {
    type_fn = if (identical(backend, "postgres")) "pg_typeof" else "typeof"
    sql_types_df = DBI::dbGetQuery(conn, paste0(
      "SELECT ", paste0(type_fn, "(", q(components), ") AS t", seq_along(components), collapse = ","),
      " FROM ", retained, " LIMIT 1"
    ))
    if (nrow(sql_types_df) != 1L) stop("Cannot discover i() levels from an empty source.", call. = FALSE)
    sql_types = vapply(sql_types_df, function(type) type[[1L]], character(1), USE.NAMES = FALSE)
    types = vapply(sql_types, function(type) category_key_type(type, backend), character(1))
  }
  if (anyNA(types)) {
    unknown = tolower(sql_types[[which(is.na(types))[[1L]]]])
    stop("Category dictionaries do not support key SQL type '", unknown, "'.", call. = FALSE)
  }
  list(types = unname(types), sql_types = unname(sql_types))
}

category_integer_in_range = function(values) {
  negative = startsWith(values, "-")
  digits = sub("^-", "", values)
  limit = ifelse(negative, "9223372036854775808", "9223372036854775807")
  !is.na(values) & grepl("^-?[0-9]+$", values) &
    (nchar(digits) < 19L | nchar(digits) == 19L & digits <= limit)
}

sql_cast_text = function(conn, expr, backend_info = detect_backend(conn)) {
  type = switch(backend_info[["name"]], spark = "STRING", mysql = "CHAR", sqlserver = "VARCHAR(MAX)", "VARCHAR")
  paste0("CAST(", expr, " AS ", type, ")")
}

category_sql_literals = function(conn, keys, type) {
  if (type == "integer") {
    target = if (detect_backend(conn)[["name"]] == "mysql") "SIGNED" else "BIGINT"
    return(paste0("CAST(", as.character(DBI::dbQuoteString(conn, as.character(keys))), " AS ", target, ")"))
  }
  if (type == "numeric") return(sql_cast_double(sprintf("%.17e", keys), conn))
  if (type == "logical") {
    bit = detect_backend(conn)[["name"]] == "sqlserver"
    return(ifelse(keys, if (bit) "1" else "TRUE", if (bit) "0" else "FALSE"))
  }
  value = category_string_literals(conn, as.character(keys))
  if (type == "date") value = paste0("CAST(", value, " AS DATE)")
  value
}

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

# Spark's parser does not share DBI's quote/backslash rules. Encode key bytes
# outside that parser; the native comparison still uses the original string.
category_string_literals = function(conn, values) {
  if (detect_backend(conn)[["name"]] != "spark") {
    return(as.character(DBI::dbQuoteString(conn, values)))
  }
  vapply(values, function(value) {
    if (is.na(value)) return("CAST(NULL AS STRING)")
    paste0("CAST(unbase64('", raw_to_base64(charToRaw(enc2utf8(value))), "') AS STRING)")
  }, character(1), USE.NAMES = FALSE)
}

# Read typed keys once for either ordinary factors or bounded i() expansion.
# Ordinary factors keep their existing unbounded discovery contract; only i()
# reads use its option-controlled limit.
read_category_dictionary = function(conn, from_sql, column_name, native_schema = NULL, limit = NULL) {
  backend_info = detect_backend(conn)
  source = sub("^FROM\\s+", "", from_sql, ignore.case = TRUE)
  adapter = backend_info[["name"]] %in% c("sqlserver", "mysql")
  metadata = if (adapter) i_native_metadata(conn, source, column_name, native_schema, backend_info) else
    category_key_metadata(conn, source, column_name)
  type = metadata[["types"]][[1L]]
  if (type == "float") stop("i() floating-point categories require a DOUBLE source column.", call. = FALSE)
  column = as.character(DBI::dbQuoteIdentifier(conn, column_name))
  private_key = as.character(DBI::dbQuoteIdentifier(conn, "dbreg_i_key__"))
  key_sql = paste0("dbreg_i_levels.", private_key)
  value = if (type == "integer") sql_cast_text(conn, key_sql, backend_info) else key_sql
  # Qualify the native key so a textual output alias cannot change its order.
  # The one-column source filters missing levels without appending to a WHERE.
  query = paste0(
    "SELECT ", value, " AS value FROM (SELECT ", column, " AS ", private_key,
    " ", from_sql, ") dbreg_i_levels WHERE ", key_sql, " IS NOT NULL GROUP BY ", key_sql,
    " ORDER BY ", key_sql
  )
  if (!is.null(limit)) query = sql_limit(query, sprintf("%.0f", limit + 2L), backend_info[["name"]])
  keys = DBI::dbGetQuery(conn, query)[[1L]]
  if (adapter && type == "logical" && !is.logical(keys)) {
    if (!is.numeric(keys) || anyNA(keys) || any(!keys %in% c(0, 1))) {
      stop("Native BIT i() levels must be 0 or 1.", call. = FALSE)
    }
    keys = as.logical(keys)
  }
  if (!is.null(limit) && length(keys) > limit + 1L) stop("i() level dictionary exceeds dbreg.i_max_columns.", call. = FALSE)
  if (adapter && type == "integer") {
    if (!is.character(keys)) stop("Exact integer i() levels require native character transport.", call. = FALSE)
    keys = sub("^[0-9]+:", "", category_key_part(keys, "integer"))
    if (anyNA(keys)) stop("Could not preserve i() integer levels exactly.", call. = FALSE)
  }
  # All other backends keep their existing text collation. The two adapters
  # retain native SQL order, including ENUMs whose driver flags are absent.
  if (!adapter && type == "character" && !startsWith(tolower(metadata[["sql_types"]][[1L]]), "enum")) keys = sort(keys)
  if (type == "integer" && !all(category_integer_in_range(keys))) stop("i() integer levels must fit signed BIGINT.", call. = FALSE)
  list(values = keys, type = type)
}
