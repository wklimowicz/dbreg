sql_prediction_sum = function(terms) {
  if (length(terms)) paste0("(", paste(terms, collapse = "+"), ")") else "0.0"
}

sql_prediction_equation = function(object, context, source, dictionaries, conn, row_id, prediction) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  product = function(coef, value) sprintf("(%s)*(%s)", dbglm_number(coef), value)
  square = function(value) if (length(value)) paste0("POWER((", value, "),2)") else character()
  x = if (length(source[["x"]])) paste0("r.", source[["x"]]) else character()
  joins = vapply(seq_along(dictionaries), function(i) {
    dictionary = dictionaries[[i]]
    keys = source[["key_sources"]][[i]]
    comparisons = vapply(seq_along(keys), function(j) {
      value = paste0("r.", source[["raw_aliases"]][[keys[[j]]]])
      cast = sql_prediction_key_cast(conn, value, dictionary[["types"]][[j]], dictionary[["sql_types"]][[j]], source[["types"]][[keys[[j]]]])
      paste0("g", i, ".", dictionary[["keys"]][[j]], "=", cast)
    }, character(1))
    paste0(" LEFT JOIN ", dictionary[["table"]], " g", i, " ON ", paste(comparisons, collapse = " AND "))
  }, character(1))
  known = if (length(dictionaries)) paste0("g", seq_along(dictionaries), ".id IS NOT NULL") else "TRUE"
  state = context[["state"]]
  null = state[["prediction"]][["null_basis"]]
  contribution = score = character()
  norm = "1.0"
  if (context[["kind"]] == "database") {
    contribution = paste0("g", seq_along(dictionaries), ".gamma")
    if (length(dictionaries) == 2L) known = c(known, "g1.component=g2.component")
    if (!is.null(null) && ncol(null)) {
      scale = state[["prediction"]][["design_scale"]]
      if (is.null(scale)) scale = rep(1, length(x))
      means = vapply(seq_along(x), function(j) sql_prediction_sum(vapply(seq_along(dictionaries), function(i)
        paste0("g", i, ".", dictionaries[[i]][["means"]][[j]]), character(1))), character(1))
      difference = sprintf("((%s-(%s))/(%s))", x, means, dbglm_number(scale))
      score = vapply(seq_len(ncol(null)), function(j) sql_prediction_sum(product(null[, j], difference)), character(1))
      norm = sql_prediction_sum(c("1.0", square(sprintf("(%s)/(%s)", x, dbglm_number(scale))),
        square(sprintf("(%s)/(%s)", means, dbglm_number(scale)))))
    }
  } else if (context[["kind"]] == "portable") {
    operations = state[["specification"]][["operations"]]
    multipliers = vapply(operations, function(op) if (is.null(op[["variable"]])) "1.0" else
      paste0("(r.", source[["varying"]][[op[["variable"]]]], "-(",
        dbglm_number(if (is.null(op[["center"]])) 0 else op[["center"]]), "))"), character(1))
    contribution = vapply(seq_along(operations), function(j) paste0("(", multipliers[[j]], ")*g", operations[[j]][["group"]], ".gamma", j), character(1))
    if (!is.null(null) && ncol(null)) {
      fe_n = state[["coefficients"]][["fe_n"]]
      score = vapply(seq_len(ncol(null)), function(k) {
        effects = vapply(seq_along(operations), function(j) paste0("(", multipliers[[j]], ")*g", operations[[j]][["group"]], ".null", j, "_", k), character(1))
        sql_prediction_sum(c(product(null[fe_n + seq_along(x), k], x), effects))
      }, character(1))
      norm = sql_prediction_sum(c(square(x), square(multipliers)))
    }
  }
  eta = sql_prediction_sum(c(product(context[["beta"]], x), contribution, dbglm_number(context[["constant"]]), paste0("r.", source[["offset"]])))
  estimable = if (length(score)) paste0(sql_prediction_sum(square(score)), "<=",
    dbglm_number(state[["prediction"]][["tolerance"]]^2), "*(", norm, ")") else "TRUE"
  factor_valid = vapply(names(context[["design"]][["factor_levels"]]), function(name) {
    levels = context[["design"]][["factor_levels"]][[name]]
    class = prediction_factor_class(context[["design"]], name)
    literals = context[["design"]][["level_sql"]][[name]]
    if (!length(literals)) return("FALSE")
    paste0("r.", q(name), " IN (", paste(literals, collapse = ","), ")")
  }, character(1))
  error = function(message) paste0("CAST(", if (isTRUE(detect_backend(conn)[["is_databricks"]])) "raise_error" else "error",
    "(", sql_value_literals(conn, message), ") AS DOUBLE)")
  cases = c(paste0("WHEN NOT (", source[["valid"]], ") THEN CAST(NULL AS DOUBLE)"),
    if (!is.null(source[["precision"]])) paste0("WHEN ", source[["precision"]], " THEN ", error("Prediction keys or factor values cannot be represented exactly; supply exact integer columns or integer FE strings.")),
    paste0("WHEN NOT (", paste(known, collapse = " AND "), ") THEN CAST(NULL AS DOUBLE)"),
    if (length(factor_valid)) paste0("WHEN NOT (", paste(factor_valid, collapse = " AND "), ") THEN ", error("newdata contains unseen factor levels.")),
    paste0("WHEN NOT (", estimable, ") THEN CAST(NULL AS DOUBLE)"),
    if (context[["kind"]] != "none") paste0("WHEN NOT isfinite(", eta, ") THEN CAST(NULL AS DOUBLE)"))
  output = c(paste0("r.", source[["key_aliases"]], " AS ", q(row_id)),
    paste0("CASE ", paste(cases, collapse = " "), " ELSE ", eta, " END AS ", q(prediction)))
  paste0("SELECT ", paste(output, collapse = ","), " FROM ", source[["frame"]], paste(joins, collapse = ""))
}

sql_prediction_response = function(conn, table, row_id, prediction, link) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  if (link == "log") {
    # R's Poisson linkinv floors tiny means at epsilon and remains finite above
    # eta=700. The conservative fitting guard must not limit new predictions.
    eta = q(prediction)
    mu = paste0("GREATEST(EXP(", eta, "),", dbglm_number(.Machine[["double.eps"]]), ")")
    return(paste0("SELECT ", paste(q(row_id), collapse = ","), ",CASE WHEN ", eta,
      " IS NULL OR isnan(", eta, ") THEN CAST(NULL AS DOUBLE) ELSE ", mu, " END AS ",
      q(prediction), " FROM ", q(table)))
  }
  keys = paste0("dbreg_prediction_key", seq_along(row_id), "__")
  mean_source = paste0("(SELECT ", paste(q(row_id), "AS", keys, collapse = ","), ",", q(prediction),
    " AS __prediction_eta,CAST(0.0 AS DOUBLE) AS __dbglm_y,CAST(1.0 AS DOUBLE) AS __dbglm_a,",
    "CAST(0.0 AS DOUBLE) AS __dbglm_o FROM ", q(table), ")")
  working = dbglm_working_sql(mean_source, character(), numeric(), link, eta = "__prediction_eta")
  paste0(working, "SELECT ", paste(keys, "AS", q(row_id), collapse = ","),
    ",CASE WHEN __prediction_eta IS NULL OR isnan(__prediction_eta) THEN CAST(NULL AS DOUBLE) ELSE mu END AS ",
    q(prediction), " FROM means")
}
