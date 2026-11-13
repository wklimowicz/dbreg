is_database_effect_state = function(state) identical(state[["storage"]], "database")

# For a single family the nuisance projection is a weighted group mean.
# The dictionary stays in SQL. Only structural information crosses DBI.
database_one_effect = function(conn, working, key_source, fe, keys, group,
                                x, beta, center, n, resources) {
  p = length(x)
  entries = getOption("dbreg.prediction_max_entries", 2e6)
  if ((p + 1)^2 > entries) {
    projection_abort("system_limit", "The database prediction information exceeds its resource budget.")
  }
  key_type = category_key_metadata(conn, key_source, keys)[["types"]]
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  aliases = paste0("k", seq_along(keys))
  means = if (p) paste0("m", seq_len(p))
  slope = if (p) paste(paste0("(", dbglm_number(beta), ")*s.", x), collapse = "+") else "0.0"
  target = "(s.eta-s.__dbglm_o)"
  weight = "s.__dbglm_a"
  wsum = sql_stable_sum(conn, weight)
  weighted = function(expr) paste0(sql_stable_sum(conn, paste0(weight, "*(", expr, ")")), "/", wsum)
  columns = c(paste0("s.", fe, " AS id"),
    paste(weighted(paste0(target, "-(", slope, ")")), "AS gamma"),
    if (p) paste(weighted(paste0("s.", x)), "AS", means))
  aggregate = paste("SELECT", paste(columns, collapse = ","), "FROM (SELECT * FROM", working, ") s GROUP BY", paste0("s.", fe))
  domain = paste("SELECT DISTINCT", paste(c(paste(fe, "AS id"), paste(q(keys), "AS", aliases)),
    collapse = ","), "FROM", key_source)
  create = fe_projection_tables(conn, resources)
  dictionary = create(paste("SELECT a.*,", paste(paste0("d.", aliases), collapse = ","),
    "FROM (", aggregate, ") a JOIN (", domain, ") d ON a.id=d.id"), c("id", "gamma", means, aliases))
  from = paste("FROM (SELECT * FROM", working, ") s JOIN", dictionary, "g ON", paste0("s.", fe, "=g.id"))
  error = paste0("ABS((", slope, ")+g.gamma-", target, ")")
  certificate = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n, MAX(", error,
    ") AS error, MAX(ABS(", target, ")) AS scale", from))
  if (!isTRUE(as.numeric(certificate[["n"]]) == n) || !is.finite(certificate[["error"]]) ||
      certificate[["error"]] > 1e-8 * max(1, certificate[["scale"]])) {
    projection_abort("certificate_failed", "The database dictionary does not reproduce the fitted predictor.")
  }
  null = NULL
  if (p) {
    within = paste0("(SELECT s.__dbglm_a AS w,", paste(paste0("(s.", x, "-g.", means, ") AS ", x),
      collapse = ","), " ", from, ") database_within")
    limits = projection_limits()
    budget = list2env(list(rows = 0, cells = 0, summary_bytes = 0), parent = emptyenv())
    H = fe_information_query(conn, within, x, budget, limits)[["cross"]]
    raw = fe_information_query(conn, paste0("(SELECT __dbglm_a AS w,",
      paste(x, collapse = ","), " FROM ", working, ") database_design"), x, budget, limits)[["cross"]]
    null = database_within_null(H, raw)
  }
  owned = resources[["tables"]][match(dictionary, q(resources[["tables"]]))]
  list(storage = "database", tables = owned,
    dictionaries = list(list(table = dictionary, columns = group, types = key_type, keys = aliases, means = means)),
    beta = beta, center = center, null = null)
}

database_within_null = function(H, raw) {
  p = nrow(H)
  absorbed = diag(H) <= 1e-12 * pmax(diag(raw), .Machine[["double.xmin"]])
  H[absorbed, ] = 0
  H[, absorbed] = 0
  system = list(n = p + 1L, big = 1L, rest = 1L + seq_len(p), Dg = 1,
    Cm = matrix(0, 1L, p), E = H, diag = c(1, diag(raw)))
  solver = fe_solver(system, new.env(parent = emptyenv()), retain_null = TRUE)
  basis = solver[["null_basis"]]
  if (length(basis)) basis[-1L, , drop = FALSE]
}

learned_effect_attempt = function(portable, database = NULL) {
  attempt = projection_attempt("portable_prediction", portable)
  if (is.null(attempt[["result"]]) && !missing(database) && !identical(attempt[["reason"]], "certificate_failed")) {
    fallback = projection_attempt("database_prediction", database)
    if (!identical(fallback[["reason"]], "not_applicable")) {
      fallback[["portable_reason"]] = attempt[["reason"]]
      attempt = fallback
    }
  }
  if (is.null(attempt[["result"]])) {
    warning("Learned fixed effects could not be retained (", attempt[["reason"]],
      "); predict() on new data is unavailable for this fit. Inspect prediction_diagnostics.", call. = FALSE)
  }
  attempt
}

predict_database_effects_r = function(object, data) {
  if (!nrow(data)) return(numeric())
  lease = object[["resource_state"]]
  conn = retained_fit_connection(lease)
  attr(conn, "dbreg_resource_lease") = lease
  id = safe_identifiers("dbreg_prediction_position", taken = names(data), fold_case = folds_identifier_case(conn))
  data[[id]] = seq_len(nrow(data))
  result = predict_keyed_sql(object, data, conn, NULL, NULL, id, NULL, for_r = TRUE)
  on.exit(close(result))
  values = as.data.frame(result)
  values[[result[["prediction_column"]]]]
}

sql_database_effects = function(conn, state, column, schema) {
  dictionaries = state[["dictionaries"]]
  aliases = paste0("g", seq_along(dictionaries))
  joins = valid = precision = character()
  means = vector("list", length(dictionaries))
  for (i in seq_along(dictionaries)) {
    dictionary = dictionaries[[i]]
    alias = aliases[i]
    rules = sql_effect_key_rules(conn, dictionary, column, schema)
    condition = paste0(rules[["values"]], "=", alias, ".", dictionary[["keys"]], collapse = " AND ")
    joins = c(joins, paste("LEFT JOIN", dictionary[["table"]], alias, "ON", condition))
    valid = c(valid, rules[["valid"]], paste0(alias, ".id IS NOT NULL"))
    precision = c(precision, rules[["precision"]])
    means[[i]] = if (length(dictionary[["means"]])) paste0(alias, ".", dictionary[["means"]])
  }
  if (length(dictionaries) == 2L) valid = c(valid, "g1.component=g2.component")
  group_means = if (length(means[[1L]])) {
    vapply(seq_along(means[[1L]]), function(j) paste0("(",
      paste(vapply(means, `[[`, character(1), j), collapse = "+"), ")"), character(1))
  }
  list(joins = joins, valid = valid, precision = precision,
    terms = paste0(aliases, ".gamma"), means = group_means)
}

sql_database_estimability = function(state, x, effects) {
  null = state[["null"]]
  if (!length(null)) return(character())
  if (length(state[["design_scale"]])) x = paste0("(", x, "/(", dbglm_number(state[["design_scale"]]), "))")
  difference = paste0("(", x, "-", effects[["means"]], ")")
  score = vapply(seq_len(ncol(null)), function(i) {
    dot = paste0("(", dbglm_number(null[, i]), ")*", difference, collapse = "+")
    paste0("POWER((", dot, "),2)")
  }, character(1))
  norm = c("1.0", paste0("POWER(", x, ",2)"), paste0("POWER(", effects[["means"]], ",2)"))
  paste0("(", paste(score, collapse = "+"), ")<=1e-14*(", paste(norm, collapse = "+"), ")")
}

database_fixed_effects = function(object, state) {
  conn = retained_fit_connection(object[["resource_state"]])
  dictionaries = state[["dictionaries"]]
  result = lapply(dictionaries, function(dictionary) {
    keys = dictionary[["keys"]]
    types = dictionary[["types"]]
    exprs = vapply(seq_along(keys), function(j) {
      if (types[[j]] == "integer") sql_cast_text(conn, keys[[j]]) else keys[[j]]
    }, character(1))
    projection = paste(c(paste(exprs, "AS", keys), "gamma"), collapse = ",")
    rows = DBI::dbGetQuery(conn, paste("SELECT", projection, "FROM", dictionary[["table"]], "ORDER BY id"))
    parts = Map(function(column, type) category_key_part(rows[[column]], type), keys, types)
    encoded = do.call(paste0, parts)
    labels = learned_effect_labels(list(columns = dictionary[["columns"]], keys = encoded))
    stats::setNames(rows[["gamma"]], labels)
  })
  constant = if (is.null(state[["constant"]])) 0 else state[["constant"]]
  result[[1L]] = result[[1L]] + constant - sum(state[["beta"]] * state[["center"]])
  names(result) = vapply(dictionaries, function(dictionary) paste(dictionary[["columns"]], collapse = "^"), character(1))
  result
}
