prediction_required_columns = function(object) {
  if (inherits(object, "dbglm")) {
    variables = object[["model_spec"]][["source_vars"]]
    for (name in intersect(variables, names(object[["factor_sources"]]))) variables[variables == name] = object[["factor_sources"]][[name]]
  } else variables = all.vars(stats::formula(object[["fml"]], lhs = 0, rhs = 1))
  state = object[["fitted_state"]]
  unique(c(variables, object[["offset"]],
    unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns")),
    unlist(lapply(state[["specification"]][["operations"]], `[[`, "variable"))))
}

prediction_factor_sources = function(object) {
  if (inherits(object, "dbglm")) {
    variables = names(object[["model_spec"]][["factor_levels"]])
    for (name in intersect(variables, names(object[["factor_sources"]]))) variables[variables == name] = object[["factor_sources"]][[name]]
    unique(variables)
  } else intersect(names(object[["factor_levels"]]), all.vars(stats::formula(object[["fml"]], lhs = 0, rhs = 1)))
}

prediction_dataframe_source = function(object, data, required) {
  data = as.data.frame(data)
  absent = setdiff(required, names(data))
  if (length(absent)) stop("newdata is missing required columns: ", paste(absent, collapse = ", "), call. = FALSE)
  state = object[["fitted_state"]]
  character_keys = unlist(lapply(state[["specification"]][["dictionaries"]], function(dictionary) {
    dictionary[["columns"]][dictionary[["types"]] == "character"]
  }), use.names = FALSE)
  factors = unique(c(prediction_factor_sources(object), character_keys))
  labels = safe_identifiers(paste0("dbreg_prediction_factor", seq_along(factors), "__"), taken = names(data))
  if (!length(factors)) labels = character()
  for (i in seq_along(factors)) data[[labels[[i]]]] = as.character(data[[factors[[i]]]])
  list(data = data[unique(c(required, labels))], labels = stats::setNames(labels, factors))
}

sql_prediction_context = function(object, conn, dataframe_labels = NULL) {
  state = object[["fitted_state"]]
  if (is.null(state[["prediction"]])) state = NULL
  if (length(object[["fe"]]) && is.null(state)) {
    stop("Database prediction with fixed effects requires fitting with retain='prediction'.", call. = FALSE)
  }
  if (!is.null(state) && !identical(state[["version"]], 1L)) stop("Unsupported fitted-state version.", call. = FALSE)
  design = sql_fitted_prediction_design(object, conn, factor_labels = !is.null(dataframe_labels))
  for (name in names(design[["factor_levels"]])) {
    public = design[["sources"]][[name]]
    if (!is.null(dataframe_labels)) design[["sources"]][name] = dataframe_labels[[public]]
  }
  if (is.null(state)) {
    beta = stats::coef(object)
    if (inherits(object, "dbglm")) names(beta) = object[["design_names"]]
    index = match(names(beta), design[["names"]])
    if (anyNA(index)) stop("Could not align fitted SQL prediction coefficients.", call. = FALSE)
    design[["expressions"]] = design[["expressions"]][index]
    design[["names"]] = design[["names"]][index]
    return(list(kind = "none", design = design, beta = unname(beta), center = rep(0, length(beta)),
      constant = 0, state = NULL, offset = object[["offset"]], varying = character()))
  }
  expected = state[["specification"]][["design_names"]]
  index = match(expected, design[["names"]])
  if (anyNA(index)) stop("Could not align the fitted learned-effect SQL design.", call. = FALSE)
  design[["expressions"]] = design[["expressions"]][index]
  design[["names"]] = expected
  database = identical(state[["prediction"]][["storage"]], "database")
  beta = if (database) state[["coefficients"]][["slopes"]] else
    state[["coefficients"]][["joint"]][state[["coefficients"]][["fe_n"]] + seq_along(expected)]
  constant = state[["prediction"]][["constant"]]
  if (is.null(constant)) constant = 0
  list(kind = if (database) "database" else "portable", design = design, beta = beta,
    center = state[["specification"]][["design_center"]], constant = constant,
    state = state, offset = object[["offset"]], varying = unique(unlist(lapply(state[["specification"]][["operations"]], `[[`, "variable"))))
}

prediction_portable_tables = function(context, conn, resources) {
  state = context[["state"]]
  if (context[["kind"]] == "none") return(list())
  if (context[["kind"]] == "database") return(state[["specification"]][["dictionaries"]])
  dictionaries = state[["specification"]][["dictionaries"]]
  coefficients = state[["coefficients"]]
  operations = state[["specification"]][["operations"]]
  null = state[["prediction"]][["null_basis"]]
  lapply(seq_along(dictionaries), function(i) {
    dictionary = dictionaries[[i]]
    frame = prediction_decode_keys(dictionary)
    key_types = vapply(seq_along(dictionary[["types"]]), function(j) {
      switch(dictionary[["types"]][[j]], integer = prediction_integer_sql_type(frame[[j]], conn),
        numeric = "DOUBLE", float = "REAL", character = "STRING", logical = "BOOLEAN", date = "DATE")
    }, character(1))
    types = stats::setNames(key_types, names(frame))
    frame[["id"]] = as.numeric(dictionary[["id"]])
    types[["id"]] = "DOUBLE"
    members = which(vapply(operations, `[[`, integer(1), "group") == i)
    for (op in members) {
      indices = coefficients[["starts"]][[op]] + match(dictionary[["id"]], coefficients[["ids"]][[op]])
      frame[[paste0("gamma", op)]] = coefficients[["joint"]][indices]
      types[[paste0("gamma", op)]] = "DOUBLE"
      if (!is.null(null)) for (j in seq_len(ncol(null))) {
        name = paste0("null", op, "_", j)
        frame[[name]] = null[indices, j]
        types[[name]] = "DOUBLE"
      }
    }
    dictionary[["table"]] = prediction_upload_frame(conn, frame, resources, types)
    dictionary[["keys"]] = paste0("k", seq_along(key_types))
    dictionary[["sql_types"]] = key_types
    dictionary[["operations"]] = members
    dictionary
  })
}

sql_prediction_numeric = function(conn, column, class) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  if (class == "Date") {
    return(if (native) paste0("DATEDIFF(", column, ",DATE '1970-01-01')") else
      paste0("DATE_DIFF('day',DATE '1970-01-01',", column, ")"))
  }
  if (class == "POSIXct") {
    return(if (native) paste0("unix_micros(", column, ")/1e6") else paste0("epoch(", column, ")"))
  }
  if (!class %in% c("numeric", "integer", "integer64")) stop("Numeric prediction inputs require numeric or date/time SQL columns.", call. = FALSE)
  paste0("CAST(", column, " AS DOUBLE)")
}

sql_prediction_key_cast = function(conn, value, type, sql_type, input_type = "STRING") {
  cast = paste0("TRY_CAST(", value, " AS ", sql_type, ")")
  if (type == "integer" && grepl("^(varchar|char|string|enum)", tolower(input_type))) {
    valid = if (isTRUE(detect_backend(conn)[["is_databricks"]])) paste0(value, " RLIKE '^[+-]?[0-9]+$'") else
      paste0("regexp_full_match(", value, ",'[+-]?[0-9]+')")
    return(paste0("CASE WHEN ", valid, " THEN ", cast, " ELSE CAST(NULL AS ", sql_type, ") END"))
  }
  if (type == "date" && grepl("^(varchar|char|string|enum)", tolower(input_type))) {
    # Date strings match the fitted textual date, rather than silently accepting
    # a timestamp or a differently formatted string for the same calendar day.
    text = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "STRING" else "VARCHAR"
    valid = paste0("((CAST(", cast, " AS ", text, ")=", value, " AND ", value,
      " NOT IN ('infinity','-infinity')) OR ", value, " IN ('Inf','-Inf'))")
    return(paste0("CASE WHEN ", valid, " THEN ", cast, " ELSE CAST(NULL AS ", sql_type, ") END"))
  }
  cast
}
