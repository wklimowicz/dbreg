# Reproduce the fitted recipe in SQL over training rows. Only the structural
# columns are compiled here. The original fixed effects and Mundlak means are
# represented by the learned group coordinates.
learned_linear_recipe_sql = function(recipe, conn, source, dictionaries, schema) {
  kind = recipe[["kind"]]
  if (kind == "constant") return(sql_cast_double("1.0", conn))
  if (kind == "product") {
    left = learned_linear_recipe_sql(recipe[["left"]], conn, source, dictionaries, schema)
    right = learned_linear_recipe_sql(recipe[["right"]], conn, source, dictionaries, schema)
    return(paste0("(", left, ")*(", right, ")"))
  }
  variable = source(recipe[["variable"]])
  if (kind %in% c("numeric", "raw")) {
    value = sql_cast_double(variable, conn)
    center = recipe[["center"]]
    return(if (is.null(center)) value else paste0("(", value, "-(", dbglm_number(center), "))"))
  }
  if (kind == "power") {
    return(paste0("POWER(", sql_cast_double(variable, conn), ",", dbglm_number(recipe[["exponent"]]), ")"))
  }
  if (kind %in% c("indicator", "i")) {
    levels = recipe[["levels"]]
    if (!is.null(recipe[["dictionary"]])) levels = dictionaries[[recipe[["dictionary"]]]]
    value = recipe[["value"]]
    type = sql_prediction_key_type(recipe, levels)
    variable = sql_prediction_category_value(conn, variable, type, schema[[recipe[["variable"]]]])
    if (type == "integer") variable = paste0("TRY_CAST(", variable, " AS ", sql_prediction_integer_type(conn), ")")
    target = sql_prediction_category_literals(conn, value, type)
    # DOUBLE indicators: DuckDB types 1.0/0.0 as DECIMAL, whose products overflow.
    out = sql_cast_double(paste0("(CASE WHEN ", variable, "=", target, " THEN 1.0 ELSE 0.0 END)"), conn)
    if (length(recipe[["second"]])) out = paste0(out, "*", sql_cast_double(source(recipe[["second"]]), conn))
    return(out)
  }
  stop("Unsupported learned OLS design recipe: ", kind, call. = FALSE)
}

linear_learned_state = function(object, inputs, rows) {
  conn = inputs[["conn"]]
  layout = object[["linear_design"]]
  columns = layout[["columns"]]
  columns = columns[columns[["role"]] == "slope" & columns[["alias"]] != "(Intercept)", , drop = FALSE]
  groups = object[["fe_components"]]
  slopes = object[["varying_slopes"]]
  if (!length(slopes)) slopes = rep(list(character()), length(groups))
  trend_names = unique(unlist(slopes, use.names = FALSE))
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  source = function(value) paste0("b.", q(value))
  fe = paste0("f", seq_along(groups))
  key_columns = lapply(seq_along(groups), function(i) paste0("k", i, "_", seq_along(groups[[i]])))
  key_fields = unlist(Map(function(g, k) paste(source(g), "AS", k), groups, key_columns), use.names = FALSE)
  ranks = vapply(groups, function(g) paste0("DENSE_RANK() OVER (ORDER BY ", paste(source(g), collapse = ","), ")"), character(1))
  x = if (nrow(columns)) paste0("x", seq_len(nrow(columns)))
  original = inputs[["retained_original"]]
  schema = DBI::dbGetQuery(conn, paste("SELECT * FROM", original, "LIMIT 0"))
  expr = vapply(columns[["recipe"]], learned_linear_recipe_sql, character(1),
    conn = conn, source = source, dictionaries = layout[["dictionaries"]], schema = schema)
  trend_fields = if (length(trend_names)) stats::setNames(paste0("trend", seq_along(trend_names)), trend_names)
  offset = if (length(object[["offset"]])) sql_cast_double(source(object[["offset"]]), conn) else "0.0"
  fields = c(paste(ranks, "AS", fe), key_fields,
    if (length(x)) paste(expr, "AS", x),
    if (length(trend_names)) paste(sql_cast_double(source(trend_names), conn), "AS", trend_fields),
    paste0("(s.__dbreg_saved_fit-", offset, ") AS eta"),
    "s.__dbreg_saved_weight AS prior_weight")
  join = paste0("s.", rows[["keys"]], "=", source(inputs[["row_id"]]), collapse = " AND ")
  create = fe_projection_tables(conn, inputs[["resources"]])
  base = create(paste("SELECT", paste(fields, collapse = ","), "FROM", rows[["table"]],
    "s JOIN", original, "b ON", join, "WHERE s.__dbreg_estimated"),
    c(fe, unlist(key_columns, use.names = FALSE), x, unname(trend_fields), "eta", "prior_weight"))
  centers_sql = function(value, name) paste0("MIN(", value, ")/2+MAX(", value, ")/2 AS ", name)
  moments = DBI::dbGetQuery(conn, paste("SELECT", paste(c("COUNT(*) AS n", "MAX(prior_weight) AS max_weight",
    centers_sql("eta", "eta"), if (length(x)) centers_sql(x, x), paste("MAX(", fe, ") AS", fe)),
    collapse = ","), "FROM", base))
  sizes = as.numeric(moments[1, fe, drop = TRUE])
  centers = if (length(x)) as.numeric(moments[1, x, drop = TRUE]) else numeric()
  constant = moments[["eta"]]
  scale = 2^floor(log2(as.numeric(moments[["max_weight"]])))
  if (!isTRUE(as.numeric(moments[["n"]]) == nobs(object)) ||
      any(!is.finite(c(centers, constant, scale))) || scale <= 0) {
    projection_abort("certificate_failed", "The retained OLS sample has invalid prediction coordinates.")
  }
  limits = projection_limits()
  check_learned_budget(sizes, slopes, length(x), groups, limits)
  dictionaries = lapply(seq_along(groups), function(i) learned_key_dictionary(conn, base,
    fe[i], key_columns[[i]], groups[[i]], sizes[i], limits[["memory_bytes"]] / (8 * length(groups))))
  coefficients = coef(object, fe = TRUE)
  positions = match(columns[["id"]], layout[["fitted_ids"]])
  beta = unname(coefficients[positions])
  beta[is.na(positions)] = 0
  design = if (length(x)) paste0("(", x, "-(", dbglm_number(centers), ")) AS ", x)
  working = paste0("(SELECT ", paste(c(fe, unname(trend_fields), design,
    paste0("(eta-(", dbglm_number(constant), ")) AS eta"),
    paste0("prior_weight/(", dbglm_number(scale), ") AS __dbglm_a"),
    "0.0 AS __dbglm_o"), collapse = ","), " FROM ", base, ") learned_ols")
  budget = list2env(list(rows = 0, cells = 0, summary_bytes = 0, upload_bytes = 0), parent = emptyenv())
  state = recover_learned_effects(conn, working, fe, slopes, trend_fields, sizes, x, beta,
    centers, nobs(object), budget, limits)
  state[["dictionaries"]] = dictionaries
  state[["design_ids"]] = columns[["id"]]
  state[["constant"]] = constant
  state
}
