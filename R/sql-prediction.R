# Compile fitted structural coordinates against a private, keyed source.
# Every input row reaches the output, including missing predictor rows.
sql_prediction_recipe_inputs = function(recipe) {
  children = Filter(Negate(is.null), recipe[c("left", "right")])
  c(recipe[["variable"]], recipe[["second"]],
    unlist(lapply(children, sql_prediction_recipe_inputs), use.names = FALSE))
}

# Fallback matrix recipes use character labels for both character and logical
# factors. Their R evaluator compares as.character() labels, independently of
# the type of the prediction input.
sql_prediction_key_type = function(recipe, levels) {
  type = recipe[["key_type"]]
  if (!is.null(type)) return(type)
  if (inherits(levels, "integer64")) "integer" else if (inherits(levels, "Date")) "date" else
    if (is.numeric(levels)) "numeric" else "character"
}

sql_prediction_character_columns = function(object, state) {
  dictionaries = object[["linear_design"]][["dictionaries"]]
  fields = function(recipe) {
    levels = recipe[["levels"]]
    if (!is.null(recipe[["dictionary"]])) levels = dictionaries[[recipe[["dictionary"]]]]
    own = if (recipe[["kind"]] %in% c("constant", "indicator", "i") &&
        !is.null(levels) && sql_prediction_key_type(recipe, levels) == "character") recipe[["variable"]]
    children = Filter(Negate(is.null), recipe[c("left", "right")])
    c(own, unlist(lapply(children, fields), use.names = FALSE))
  }
  unique(c(unlist(lapply(sql_prediction_design(object, state), fields), use.names = FALSE),
    unlist(lapply(state[["dictionaries"]], function(dictionary) {
      dictionary[["columns"]][dictionary[["types"]] == "character"]
    }), use.names = FALSE)))
}

# R's double and logical formatting is not SQL CAST formatting (1 versus 1.0,
# TRUE versus true). Keep native row keys intact and upload separate labels.
sql_prediction_character_input = function(data, fields) {
  fields = intersect(fields, names(data))
  if (!length(fields)) return(list(data = data, map = NULL))
  aliases = safe_identifiers(paste0("dbreg_character_key", seq_along(fields)),
    taken = names(data), fold_case = TRUE)
  map = stats::setNames(aliases, fields)
  for (name in fields) data[[map[[name]]]] = as.character(data[[name]])
  list(data = data, map = map)
}

validate_sql_character_keys = function(conn, input, fields) {
  if (!length(fields)) return(invisible(NULL))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  types = DBI::dbGetQuery(conn, paste("SELECT", paste0("typeof(", q(fields),
    ") AS k", seq_along(fields), collapse = ","), "FROM", input, "LIMIT 1"))
  if (!nrow(types)) return(invisible(NULL))
  valid = vapply(types, function(value) category_key_type(value[[1L]]) %in% c("character", "integer", "logical"), logical(1))
  if (any(!valid)) stop("Character prediction keys require character, integer or logical database source columns; ",
    "convert these columns explicitly to character: ", paste(fields[!valid], collapse = ", "), ".", call. = FALSE)
  invisible(NULL)
}

sql_prediction_category_value = function(conn, value, type, column) {
  if (type == "character") {
    if (is.logical(column)) return(paste0("(CASE WHEN ", value, " IS NULL THEN ",
      sql_cast_text(conn, "NULL"), " WHEN ", value, " THEN 'TRUE' ELSE 'FALSE' END)"))
    return(sql_cast_text(conn, value))
  }
  if (type == "date" && is.character(column)) {
    converted = paste0("TRY_CAST(", value, " AS DATE)")
    return(paste0("(CASE WHEN ", sql_cast_text(conn, converted), "=", value,
      " THEN ", converted, " END)"))
  }
  value
}

# Learned FE dictionaries include native integer domains wider than BIGINT.
# Prediction comparisons retain the backend's full exact integer range.
sql_prediction_integer_type = function(conn) {
  backend = detect_backend(conn)
  if (backend[["name"]] == "duckdb") "HUGEINT" else
    if (isTRUE(backend[["is_databricks"]])) "DECIMAL(38,0)" else "BIGINT"
}

sql_prediction_category_literals = function(conn, keys, type) {
  if (type != "integer") return(category_sql_literals(conn, keys, type))
  paste0("CAST(", category_string_literals(conn, as.character(keys)), " AS ",
    sql_prediction_integer_type(conn), ")")
}

sql_prediction_integer_key = function(conn, value, column) {
  precision = valid = NULL
  if (is.numeric(column) && !inherits(column, "integer64")) {
    kind = paste0("LOWER(typeof(", value, "))")
    floating = paste0(kind, " IN ('double','float','real')")
    number = sql_cast_double(value, conn)
    bound = paste0("(CASE WHEN ", kind, " IN ('float','real') THEN 16777216.0 ELSE 9007199254740992.0 END)")
    precision = paste0("NOT (", floating, " AND ABS(", number, ")>", bound, ")")
    # Check native DECIMAL values before casting: integer casts round their
    # fractional parts, and converting to double can erase a small fraction.
    valid = paste0("(", sql_is_finite(conn, number), " AND ", value,
      "=FLOOR(", value, "))")
  }
  if (is.character(column)) {
    valid = if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
      paste0("(", value, " RLIKE '^[+-]?[0-9]+$')")
    } else paste0("regexp_full_match(", value, ",'[+-]?[0-9]+')")
  }
  list(value = paste0("TRY_CAST(", value, " AS ", sql_prediction_integer_type(conn), ")"),
    precision = precision, valid = valid)
}

sql_prediction_recipe_rules = function(recipe, conn, column, schema, dictionaries) {
  kind = recipe[["kind"]]
  if (kind == "product") {
    return(c(sql_prediction_recipe_rules(recipe[["left"]], conn, column, schema, dictionaries),
      sql_prediction_recipe_rules(recipe[["right"]], conn, column, schema, dictionaries)))
  }
  name = recipe[["variable"]]
  if (!length(name)) return(list())
  value = column(name)
  numeric = function(field) {
    if (!(is.numeric(schema[[field]]) || is.logical(schema[[field]])) ||
        inherits(schema[[field]], c("Date", "POSIXt"))) {
      stop("Prediction variable must be numeric: ", field, call. = FALSE)
    }
    paste0("(", column(field), " IS NOT NULL AND ",
      sql_is_finite(conn, sql_cast_double(column(field), conn)), ")")
  }
  if (kind %in% c("numeric", "raw", "power")) return(list(list(missing = numeric(name))))
  levels = recipe[["levels"]]
  if (!is.null(recipe[["dictionary"]])) levels = dictionaries[[recipe[["dictionary"]]]]
  missing = c(paste(value, "IS NOT NULL"), if (length(recipe[["second"]])) numeric(recipe[["second"]]))
  if (is.null(levels)) return(list(list(missing = missing)))
  type = sql_prediction_key_type(recipe, levels)
  value = sql_prediction_category_value(conn, value, type, schema[[name]])
  if (type == "date" && is.character(schema[[name]])) missing = c(missing, paste(value, "IS NOT NULL"))
  integer = if (type == "integer") sql_prediction_integer_key(conn, value, schema[[name]])
  if (length(integer)) value = integer[["value"]]
  missing = c(missing, integer[["valid"]])
  domain = paste0("COALESCE(", value, " IN (",
    paste(sql_prediction_category_literals(conn, levels, type), collapse = ","), "),FALSE)")
  list(list(missing = missing, precision = integer[["precision"]],
    indicator = if (kind == "i") domain else NULL,
    factor = if (kind != "i") domain else NULL))
}

sql_prediction_design = function(object, state) {
  columns = object[["linear_design"]][["columns"]]
  if (!is.null(state[["design_ids"]])) {
    columns = columns[match(state[["design_ids"]], columns[["id"]]), , drop = FALSE]
  }
  columns[["recipe"]]
}

sql_prediction_columns = function(object, row_id, state) {
  recipes = sql_prediction_design(object, state)
  groups = unlist(lapply(state[["dictionaries"]], `[[`, "columns"), use.names = FALSE)
  trends = unlist(lapply(state[["ops"]], `[[`, "trend"), use.names = FALSE)
  unique(c(row_id, object[["offset"]], groups, trends,
    unlist(lapply(recipes, sql_prediction_recipe_inputs), use.names = FALSE)))
}

sql_fitted_predictor = function(object, conn, input, row_id, prediction, state, resources, character_map = NULL) {
  validate_linear_layout(object)
  layout = object[["linear_design"]]
  fitted = match(layout[["fitted_ids"]], layout[["columns"]][["id"]])
  recipes = sql_prediction_design(object, state)
  dictionaries = layout[["dictionaries"]]
  offset = object[["offset"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  raw_column = function(x) paste0("r.", q(x))
  column = function(x) {
    mapped = character_map[x]
    if (length(mapped)) x[!is.na(mapped)] = unname(mapped[!is.na(mapped)])
    raw_column(x)
  }
  schema = DBI::dbGetQuery(conn, paste("SELECT * FROM", input, "LIMIT 0"))
  character_fields = sql_prediction_character_columns(object, state)
  mapped = character_map[character_fields]
  if (length(mapped)) character_fields[!is.na(mapped)] = unname(mapped[!is.na(mapped)])
  validate_sql_character_keys(conn, input, character_fields)
  if (length(character_map)) schema[names(character_map)] = schema[unname(character_map)]
  rules = unique(unlist(lapply(recipes, sql_prediction_recipe_rules,
    conn = conn, column = column, schema = schema, dictionaries = dictionaries), recursive = FALSE))
  if (length(offset)) {
    rules = c(rules, sql_prediction_recipe_rules(list(kind = "numeric", variable = offset),
      conn, column, schema, dictionaries))
  }
  conditions = function(name) unique(unlist(lapply(rules, `[[`, name), use.names = FALSE))
  expr = vapply(if (is.null(state)) recipes[fitted] else recipes, learned_linear_recipe_sql, character(1),
    conn = conn, source = column, dictionaries = dictionaries, schema = schema)
  effects = NULL
  beta = coef(object, fe = TRUE)
  if (!is.null(state)) {
    effects = if (is_database_effect_state(state)) sql_database_effects(conn, state, column, schema)
      else sql_prediction_effects(conn, state, column, schema, resources)
    beta = state[["beta"]]
    if (length(expr)) expr = paste0("((", expr, ")-(", dbglm_number(state[["center"]]), "))")
  } else if (!is.null(object[["prediction_coordinates"]])) {
    sc = object[["prediction_coordinates"]]
    beta = sc[["beta"]]
    if (length(expr)) {
      for (j in seq_along(expr)) {
        if (sc[["centers"]][j] != 0 || sc[["scales"]][j] != 1) {
          expr[j] = paste0("(((", expr[j], ")-(", dbglm_number(sc[["centers"]][j]),
            "))/(", dbglm_number(sc[["scales"]][j]), "))")
        }
      }
    }
  }
  constant_val = if (!is.null(state[["constant"]])) state[["constant"]]
    else if (!is.null(object[["prediction_coordinates"]])) object[["prediction_coordinates"]][["y_center"]]
    else 0
  constant = if (length(constant_val) && constant_val != 0) dbglm_number(constant_val) else "0.0"
  terms = c(constant, if (length(beta)) paste0("(", dbglm_number(beta), ")*(", expr, ")"), effects[["terms"]],
    if (length(offset)) sql_cast_double(column(offset), conn))
  eta = paste(terms, collapse = "+")
  check = function(parts, action) {
    if (length(parts)) paste0("WHEN NOT (", paste(parts, collapse = " AND "), ") THEN ", action)
  }
  backend = detect_backend(conn)
  error = function(message) paste0("CAST(",
    if (isTRUE(backend[["is_databricks"]])) "raise_error" else "error",
    "(", as.character(DBI::dbQuoteString(conn, message)), ") AS DOUBLE)")
  absent = "CAST(NULL AS DOUBLE)"
  branches = c(check(conditions("missing"), absent),
    check(c(conditions("precision"), effects[["precision"]]), error("Prediction integer keys require exact values.")),
    check(effects[["valid"]], absent),
    check(conditions("indicator"), absent),
    check(conditions("factor"), if (isTRUE(state[["r_prediction"]])) absent else error("newdata contains unseen factor levels.")),
    if (is_database_effect_state(state)) check(sql_database_estimability(state, expr, effects), absent) else
      if (!is.null(state)) check(sql_learned_estimability(state, expr, effects), absent),
    check(sql_is_finite(conn, eta), absent))
  value = paste("CASE", paste(branches, collapse = " "), "ELSE", eta, "END")
  keys = paste(raw_column(row_id), "AS", q(row_id))
  query = paste("SELECT", paste(c(keys, paste(value, "AS", q(prediction))), collapse = ","),
    "FROM", input, "r", paste(effects[["joins"]], collapse = " "))
  if (!is.null(state) && nchar(query, type = "bytes") > projection_limits()[["upload_bytes"]]) {
    projection_abort("transfer_limit", "Portable SQL prediction equation exceeds the upload budget.")
  }
  query
}

sql_prediction_requested = function(conn, table, path, row_id, prediction_table) {
  !all(vapply(list(conn, table, path, row_id, prediction_table), is.null, logical(1)))
}

predict_keyed_sql = function(object, newdata, conn, table, path, row_id, prediction_table, type = "link", for_r = FALSE) {
  if (is.null(row_id)) stop("Database prediction requires explicit row_id columns.", call. = FALSE)
  validate_sample_args(row_id, prediction_table, conn, table, FALSE, FALSE)
  if (sum(!vapply(list(newdata, table, path), is.null, logical(1))) != 1L) {
    stop("Provide exactly one of newdata, table, or path for database prediction.", call. = FALSE)
  }
  if (inherits(object, "dbivreg") || (!for_r && identical(object[["strategy"]], "mundlak"))) {
    stop("Database prediction does not support IV or Mundlak fits.", call. = FALSE)
  }
  state = if (length(object[["fe"]])) object[["fitted_state"]]
  if (length(object[["fe"]]) && is.null(state)) {
    stop("Database prediction with fixed effects requires a retained learned-effect state.", call. = FALSE)
  }
  if (is_database_effect_state(state)) {
    owner = retained_fit_connection(object[["resource_state"]])
    if (is.null(conn) && inherits(table, "tbl_lazy")) conn = dbplyr::remote_con(table)
    if (!is.null(conn) && !resource_same_connection(conn, owner)) {
      stop("Database-backed effects require the fitted owning session.", call. = FALSE)
    }
    conn = owner
    attr(conn, "dbreg_resource_lease") = object[["resource_state"]]
    state[["r_prediction"]] = for_r
  }
  if (is.null(conn) && !is.null(table) && !inherits(table, "tbl_lazy")) {
    stop("A named table requires a database connection.", call. = FALSE)
  }
  inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) retained_fit_connection(inherited)
  resources = fit_resources()
  on.exit(release_fit_resources(resources))
  needed = sql_prediction_columns(object, row_id, state)
  character_map = NULL
  if (is.data.frame(newdata)) {
    normalized = sql_prediction_character_input(as.data.frame(newdata), sql_prediction_character_columns(object, state))
    newdata = normalized[["data"]]
    character_map = normalized[["map"]]
    needed = c(needed, unname(character_map))
  }
  if (!is.null(newdata) && !is.null(conn) && isTRUE(detect_backend(conn)[["is_databricks"]])) {
    resources[["conn"]] = conn
    validate_native_temp_query(conn)
    setup = list(conn = conn, from_statement = paste("FROM",
      sql_native_prediction_input(conn, newdata, needed, resources)))
  } else {
    setup = setup_db_connection(conn, table, newdata, path, caller = "prediction", resources = resources)
  }
  conn = setup[["conn"]]
  if (is.null(inherited)) inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) retained_fit_connection(inherited)
  if (!is.null(prediction_table) && !is.null(inherited) && isTRUE(inherited[["session"]][["own_conn"]])) {
    stop("Prediction exports require a caller-owned connection.", call. = FALSE)
  }
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Database prediction requires DuckDB or native Databricks.", call. = FALSE)
  }
  validate_native_temp_query(conn, backend)
  if (!is.null(prediction_table) && sample_table_exists(conn, prediction_table)) {
    stop("Prediction export table already exists.", call. = FALSE)
  }
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  prediction = safe_identifiers("prediction", taken = row_id, fold_case = folds_identifier_case(conn))
  fields = c(row_id, prediction)
  create = fe_projection_tables(conn, resources)
  # Snapshot the source before checking the keys and evaluating fitted recipes.
  input_schema = DBI::dbGetQuery(conn, paste("SELECT *", setup[["from_statement"]], "LIMIT 0"))
  if (!all(needed %in% names(input_schema))) {
    stop("Missing prediction source columns: ", paste(setdiff(needed, names(input_schema)), collapse = ", "), call. = FALSE)
  }
  input = create(paste("SELECT", paste(q(needed), collapse = ","), setup[["from_statement"]]), needed)
  from = paste("FROM", input)
  validate_source_row_id(conn, from, row_id)
  output = create(sql_fitted_predictor(object, conn, input, row_id, prediction, state, resources, character_map), fields)
  if (!is.null(state)) validate_source_row_id(conn, paste("FROM", output), row_id)
  if (type == "response" && inherits(object, "dbglm")) {
    eta = q(prediction)
    link = object[["family"]][["link"]]
    mean = if (link == "log") paste0("CASE WHEN ", eta, ">",
      dbglm_number(log(.Machine[["double.xmax"]])), " THEN CAST('Infinity' AS DOUBLE) ELSE EXP(",
      eta, ") END") else dbglm_family_sql(link, eta)[["mu"]]
    value = paste0("CASE WHEN ", eta, " IS NULL OR isnan(", eta, ") THEN CAST(NULL AS DOUBLE) ELSE ", mean, " END")
    output = create(paste("SELECT", paste(c(q(row_id), paste(value, "AS", eta)), collapse = ","),
      "FROM", output), fields)
  }
  if (is.null(prediction_table)) {
    name = resources[["tables"]][match(output, q(resources[["tables"]]))]
    lease = retain_fit_resources(resources, name, session = inherited[["session"]])
  } else {
    DBI::dbExecute(conn, sql_sample_table(conn, prediction_table, paste("SELECT * FROM", output), fields))
    name = prediction_table
    lease = retain_fit_resources(resources, character(), session = inherited[["session"]])
  }
  attr(conn, "dbreg_resource_lease") = lease
  structure(list(conn = conn, table = name, row_id = row_id, prediction_column = prediction,
    resource_state = lease), class = "dbreg_prediction")
}

#' Keyed database predictions
#'
#' A prediction result owns a temporary table in its database session. Use
#' `head()` for a bounded preview and `as.data.frame()` for collection. Close
#' the result to release its table. Caller-named exports remain available.
#' Collection and previews use ascending native-key order, rather than the
#' original order of an R upload. Match rows using the explicit key columns.
#' Learned character categories use R's \code{as.character()} labels for R
#' prediction inputs, while preserving the original row keys. Database sources
#' may supply character, exact integer or logical columns for those categories;
#' logical values use the labels \code{"TRUE"} and \code{"FALSE"}. Floating-point
#' database columns require an explicit conversion to the intended character
#' labels. Date categories accept Date columns or canonical ISO date strings;
#' malformed strings produce missing predictions.
#' @param x,con A keyed database prediction result.
#' @param n Maximum rows to preview.
#' @param ... Unused.
#' @return A data frame from collection, invisibly `NULL` from `close()`,
#'   or the input invisibly from `print()`.
#' @name dbreg-prediction
NULL

collect_sql_prediction = function(x, n = NULL) {
  conn = retained_fit_connection(x[["resource_state"]])
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  keys = x[["row_id"]]
  fields = c(keys, x[["prediction_column"]])
  table = q(x[["table"]])
  select = q(fields)
  backend = detect_backend(conn)
  if (isTRUE(backend[["is_databricks"]])) {
    types = DBI::dbGetQuery(conn, paste("SELECT",
      paste0("typeof(", q(keys), ") AS k", seq_along(keys), collapse = ","),
      "FROM", table, "LIMIT 1"))
    exact = which(vapply(types, function(value) length(value) && grepl("^(bigint$|decimal)", tolower(value[1L])), logical(1)))
    select[exact] = paste(sql_cast_text(conn, select[exact]), "AS", q(keys[exact]))
  }
  sql = paste("SELECT", paste(select, collapse = ","), "FROM", table,
    "ORDER BY", paste(paste0(table, ".", q(keys)), collapse = ","))
  if (!is.null(n)) sql = sql_limit(sql, n, backend[["name"]])
  read_keyed_query(conn, sql, keys, exact_keys = TRUE)
}

#' @rdname dbreg-prediction
#' @export
as.data.frame.dbreg_prediction = function(x, ...) collect_sql_prediction(x)

#' @rdname dbreg-prediction
#' @export
head.dbreg_prediction = function(x, n = 6L, ...) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) {
    stop("n must be one nonnegative whole number.", call. = FALSE)
  }
  collect_sql_prediction(x, n)
}

#' @rdname dbreg-prediction
#' @export
close.dbreg_prediction = function(con, ...) close_retained_lease(con[["resource_state"]])

#' @rdname dbreg-prediction
#' @export
print.dbreg_prediction = function(x, ...) {
  cat("Database predictions; keys: ", paste(x[["row_id"]], collapse = ", "),
    "; value: ", x[["prediction_column"]], ". Use head() to preview.\n", sep = "")
  invisible(x)
}
