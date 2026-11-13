# Recover fixed effects in bounded training coordinates and certify that they
# reproduce the fitted predictor after subtracting its structural contribution.
recover_learned_effects = function(conn, working, fe, slopes, trends, sizes, x,
                                    beta, centers, fitted_rows, budget, limits) {
  sql_ops = lapply(fe, function(key) list(f = key, v = "1.0"))
  ops = lapply(seq_along(fe), function(i) list(group = i))
  fields = character()
  for (i in seq_along(fe)) for (trend in slopes[[i]]) {
    raw = trends[[trend]]
    group = fe[i]
    name = paste0("vv", length(ops) + 1L)
    partition = paste0(" OVER (PARTITION BY ", group, ")")
    lo = paste0("MIN(", raw, ")")
    hi = paste0("MAX(", raw, ")")
    center = paste0(hi, partition, "/2+", lo, partition, "/2")
    scale = paste0("COALESCE(NULLIF(", hi, partition, "/2-", lo, partition, "/2,0),1)")
    fields = c(fields, paste0("(", raw, "-(", center, "))/(", scale, ") AS ", name))
    key_info = fe_coordinate_fetch(conn, paste0("SELECT ", group, " AS id,",
      hi, "/2+", lo, "/2 AS center,COALESCE(NULLIF(", hi, "/2-", lo,
      "/2,0),1) AS scale FROM ", working, " GROUP BY ", group, " ORDER BY id"),
      3L, budget, limits)
    if (!identical(as.numeric(key_info[["id"]]), as.numeric(seq_len(sizes[i])))) {
      stop("Varying prediction locations do not match fitted ranks.", call. = FALSE)
    }
    sql_ops[[length(sql_ops) + 1L]] = list(f = group, v = name)
    ops[[length(ops) + 1L]] = list(group = i, trend = trend,
      center = as.numeric(key_info[["center"]]), scale = as.numeric(key_info[["scale"]]))
  }
  structural = if (length(x)) paste(paste0("(", dbglm_number(beta), ")*", x), collapse = "+") else "0.0"
  target = paste0("eta-__dbglm_o-(", structural, ") AS prediction_target")
  source_sql = paste("SELECT", paste(c(fe, x, fields, "__dbglm_a AS w", target),
    collapse = ","), "FROM", working)
  source = paste0("(", source_sql, ") varying_prediction_source")
  recovered = fe_coordinate_system(conn, source, sql_ops, c(x, "prediction_target"),
    budget, limits)
  sys = recovered[["sys"]]
  alpha = as.numeric(recovered[["coefficients"]][, length(x) + 1L])
  info = recovered[["info"]]
  null = learned_design_null(conn, sys, info, source, x, budget, limits)
  joins = terms = character()
  for (i in seq_along(ops)) {
    op = sql_ops[[i]]
    j = match(list(op), recovered[["ops"]])
    index = rep(NA_integer_, sizes[ops[[i]][["group"]]])
    if (!is.na(j)) {
      ids = info[[j]][["id"]]
      at = sys[["starts"]][j] + seq_along(ids)
      index[ids] = at
      lookup = matrix(0, length(index), 1L)
      lookup[ids, 1L] = alpha[at]
      sql = fe_coefficient_sql(conn, lookup, budget, limits)
      alias = paste0("vp", i)
      joins = c(joins, paste0(" JOIN (", sql, ") ", alias, " ON s.", op[["f"]], "=", alias, ".id"))
      axis = if (op[["v"]] == "1.0") "1.0" else paste0("s.", op[["v"]])
      terms = c(terms, paste0(alias, ".c1*", axis))
    }
    ops[[i]][["index"]] = index
  }
  certified = DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n,MAX(ABS((",
    paste(terms, collapse = "+"), ")-s.prediction_target)) AS error,",
    "MAX(ABS(s.prediction_target)) AS scale FROM (", source_sql, ") s",
    paste(joins, collapse = "")))
  if (!isTRUE(as.numeric(certified[["n"]]) == fitted_rows) ||
      !is.finite(certified[["error"]]) || certified[["error"]] > 1e-8 * max(1, certified[["scale"]])) {
    projection_abort("certificate_failed", "Recovered effects do not reproduce the fitted predictor.")
  }
  list(alpha = alpha, beta = beta, center = centers, null = null, ops = ops)
}
