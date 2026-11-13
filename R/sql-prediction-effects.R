# Upload only the bounded portable effect dictionaries, never prediction rows.
sql_prediction_effect_table = function(conn, dictionary, ops, indices, state, resources, budget) {
  keys = dictionary[["keys"]]
  types = dictionary[["types"]]
  width = length(types)
  decoded = vapply(keys, decode_learned_key, character(width), width = width, USE.NAMES = FALSE)
  dim(decoded) = c(width, length(keys))
  fields = paste0("k", seq_len(width))
  values = lapply(seq_len(width), function(i) {
    type = types[i]
    value = decoded[i, ]
    if (type %in% c("numeric", "float")) value = as.numeric(value)
    if (type == "logical") value = value == "TRUE"
    literal = sql_prediction_category_literals(conn, value, if (type == "float") "numeric" else type)
    if (type == "float") paste0("CAST(", literal, " AS FLOAT)") else literal
  })
  null = state[["null"]]
  number = function(x) sql_cast_double(dbglm_number(x), conn)
  for (i in indices) {
    op = ops[[i]]
    index = op[["index"]]
    present = !is.na(index)
    coefficient = numeric(length(keys))
    coefficient[present] = state[["alpha"]][index[present]]
    fields = c(fields, paste0("a", i), paste0("ok", i))
    values = c(values, list(number(coefficient), ifelse(present, "TRUE", "FALSE")))
    if (!is.null(op[["trend"]])) {
      fields = c(fields, paste0("center", i), paste0("scale", i))
      values = c(values, list(number(op[["center"]]), number(op[["scale"]])))
    }
    if (length(null)) {
      coordinates = matrix(0, nrow = length(keys), ncol = ncol(null))
      coordinates[present, ] = null[index[present], , drop = FALSE]
      fields = c(fields, paste0("n", i, "_", seq_len(ncol(null))))
      values = c(values, lapply(seq_len(ncol(null)), function(j) number(coordinates[, j])))
    }
  }
  sql_prediction_values_table(conn, fields, values, resources, budget)
}

sql_effect_key_rules = function(conn, dictionary, column, schema) {
  keys = dictionary[["columns"]]
  values = column(keys)
  types = dictionary[["types"]]
  precision = valid = character()
  for (i in seq_along(keys)) values[i] = sql_prediction_category_value(conn, values[i], types[i], schema[[keys[i]]])
  for (i in which(types == "integer")) {
    integer = sql_prediction_integer_key(conn, values[i], schema[[keys[i]]])
    precision = c(precision, integer[["precision"]])
    valid = c(valid, integer[["valid"]])
    values[i] = integer[["value"]]
  }
  for (i in which(types == "logical")) {
    if (!is.logical(schema[[keys[i]]])) stop("Logical prediction keys require logical source columns.", call. = FALSE)
  }
  floats = which(types == "float")
  values[floats] = if (length(floats)) paste0("CAST(", values[floats], " AS FLOAT)") else character()
  list(values = values, valid = valid, precision = precision)
}

sql_prediction_effects = function(conn, state, column, schema, resources) {
  dictionaries = state[["dictionaries"]]
  ops = learned_effect_operations(state)
  limits = projection_limits()
  budget = list2env(list(used = 0, limit = min(limits[["upload_bytes"]], limits[["memory_bytes"]] / 8)),
    parent = emptyenv())
  groups = vapply(ops, `[[`, integer(1), "group")
  aliases = paste0("g", seq_along(dictionaries))
  group_alias = aliases[groups]
  joins = precision = valid = character()
  for (j in seq_along(dictionaries)) {
    dictionary = dictionaries[[j]]
    keys = dictionary[["columns"]]
    key_rules = sql_effect_key_rules(conn, dictionary, column, schema)
    values = key_rules[["values"]]
    valid = c(valid, key_rules[["valid"]])
    precision = c(precision, key_rules[["precision"]])
    table = sql_prediction_effect_table(conn, dictionary, ops, which(groups == j), state, resources, budget)
    equal = paste0(values, "=", aliases[j], ".k", seq_along(keys))
    joins = c(joins, paste("LEFT JOIN", table, aliases[j], "ON", paste(equal, collapse = " AND ")))
    valid = c(valid, paste0(aliases[j], ".k1 IS NOT NULL"))
  }
  multiplier = vapply(seq_along(ops), function(i) {
    trend = ops[[i]][["trend"]]
    if (is.null(trend)) return("1.0")
    if (!(is.numeric(schema[[trend]]) || is.logical(schema[[trend]]))) {
      stop("Prediction requires numeric trend column: ", trend, call. = FALSE)
    }
    paste0("((", sql_cast_double(column(trend), conn), "-", group_alias[i], ".center", i,
      ")/", group_alias[i], ".scale", i, ")")
  }, character(1))
  valid = c(valid, sql_is_finite(conn, multiplier),
    paste0("(", group_alias, ".ok", seq_along(ops), " OR ABS(", multiplier, ")<=1e-10)"))
  list(joins = joins, valid = valid, precision = precision, multiplier = multiplier,
    terms = paste0(multiplier, "*", group_alias, ".a", seq_along(ops)), aliases = group_alias)
}

sql_learned_estimability = function(state, x, effects) {
  null = state[["null"]]
  if (!length(null)) return(character())
  axis = seq_along(effects[["multiplier"]])
  score = vapply(seq_len(ncol(null)), function(k) {
    structural = if (length(x)) paste0("(", dbglm_number(null[length(state[["alpha"]]) + seq_along(x), k]),
      ")*", x)
    grouped = paste0(effects[["multiplier"]], "*", effects[["aliases"]], ".n", axis, "_", k)
    paste0("POWER((", paste(c(structural, grouped), collapse = "+"), "),2)")
  }, character(1))
  norm = paste0("POWER(", c(x, effects[["multiplier"]]), ",2)")
  paste0("(", paste(score, collapse = "+"), ")<=1e-14*(", paste(norm, collapse = "+"), ")")
}
