# Propagate the smallest first-family rank over the bipartite training graph.
# Pointer compression accelerates propagation without joining observations in R.
database_effect_components = function(conn, source, fe, manager, max_iter) {
  create = manager[["create"]]
  drop = manager[["drop"]]
  edges = create(paste("SELECT DISTINCT", fe[1L], "AS a,", fe[2L], "AS b FROM", source))
  left = create(paste("SELECT DISTINCT a AS id,a AS component FROM", edges))
  right = create(paste("SELECT b AS id,MIN(a) AS component FROM", edges, "GROUP BY b"))
  for (i in seq_len(max_iter)) {
    next_left = create(paste("SELECT e.a AS id,MIN(r.component) AS component FROM", edges,
      "e JOIN", right, "r ON e.b=r.id GROUP BY e.a"))
    compressed = create(paste("SELECT a.id,b.component FROM", next_left, "a JOIN", next_left,
      "b ON a.component=b.id"))
    next_right = create(paste("SELECT e.b AS id,MIN(l.component) AS component FROM", edges,
      "e JOIN", compressed, "l ON e.a=l.id GROUP BY e.b"))
    drop(left)
    drop(right)
    drop(next_left)
    left = compressed
    right = next_right
    mismatch = DBI::dbGetQuery(conn, paste("SELECT 1 AS mismatch FROM", edges, "e LEFT JOIN", left,
      "l ON e.a=l.id LEFT JOIN", right,
      "r ON e.b=r.id WHERE l.component IS NULL OR r.component IS NULL OR l.component<>r.component LIMIT 1"))
    if (!nrow(mismatch)) {
      drop(edges)
      return(c(left, right))
    }
  }
  projection_abort("nonconverged", "Database prediction components did not converge.")
}

database_two_effects = function(conn, working, key_source, fe, key_columns, groups,
                                 x, beta, center, n, resources) {
  p = length(x)
  if ((p + 1)^2 > getOption("dbreg.prediction_max_entries", 2e6)) {
    projection_abort("system_limit", "The database prediction information exceeds its resource budget.")
  }
  manager = fe_cg_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  slope = if (p) paste(paste0("(", dbglm_number(beta), ")*", x), collapse = "+") else "0.0"
  target = paste0("(eta-__dbglm_o-(", slope, "))")
  columns = c("target", x)
  maximum = paste0("MAX(ABS(", c(target, x), ")) AS ", columns)
  scale = pmax(1, as.numeric(DBI::dbGetQuery(conn,
    paste("SELECT", paste(maximum, collapse = ","), "FROM", working))[1L, ]))
  if (any(!is.finite(scale))) projection_abort("certificate_failed", "Invalid database prediction coordinate scale.")
  source = create(paste("SELECT", paste(c(fe, "__dbglm_a AS w",
    paste0("(", c(target, x), ")/(", dbglm_number(scale), ") AS ", columns)), collapse = ","),
    "FROM", working))
  ops = lapply(fe, function(key) list(f = key, v = "1.0"))
  limits = projection_limits()
  budget = list2env(list(rows = 0, cells = 0, summary_bytes = 0), parent = emptyenv())
  control = fe_cg_control()
  fit = fe_sql_cg(conn, source, columns, fe, ops, rep(1, length(columns)),
    1e-12, control, create, drop, budget, limits, retain_coefficients = TRUE)
  certified = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n, MAX(ABS(target)) AS error FROM", fit[["table"]]))
  if (!isTRUE(as.numeric(certified[["n"]]) == n) || !is.finite(certified[["error"]]) ||
      certified[["error"]] > 1e-8) {
    projection_abort("certificate_failed", "Two database dictionaries do not reproduce the fitted predictor.")
  }
  null = NULL
  if (p) {
    H = fe_information_query(conn, fit[["table"]], x, budget, limits)[["cross"]]
    raw = fe_information_query(conn, source, x, budget, limits)[["cross"]]
    null = database_within_null(H, raw)
  }
  labels = database_effect_components(conn, source, fe, manager, control[["fe_max_iter"]])
  means = if (p) paste0("m", seq_len(p))
  dictionaries = lapply(seq_along(fe), function(i) {
    keys = key_columns[[i]]
    key_type = category_key_metadata(conn, key_source, keys)[["types"]]
    aliases = paste0("k", seq_along(keys))
    domain = paste("SELECT DISTINCT", paste(c(paste(fe[i], "AS id"), paste(q(keys), "AS", aliases)),
      collapse = ","), "FROM", key_source)
    fields = c("a.__cg_group AS id", paste0("a.c1*(", dbglm_number(scale[1L]), ") AS gamma"),
      if (p) paste(paste0("a.c", seq_len(p) + 1L), "AS", means),
      "l.component", paste0("d.", aliases))
    table = create(paste("SELECT", paste(fields, collapse = ","), "FROM", fit[["coefficients"]],
      "a JOIN", labels[i], "l ON a.__cg_group=l.id JOIN (", domain,
      ") d ON a.__cg_group=d.id WHERE a.__cg_op=", i))
    list(table = table, columns = groups[[i]], types = key_type, keys = aliases, means = means)
  })
  for (table in c(source, fit[["table"]], fit[["coefficients"]], labels)) drop(table)
  owned = resources[["tables"]][match(vapply(dictionaries, `[[`, character(1), "table"),
    q(resources[["tables"]]))]
  list(storage = "database", tables = owned, dictionaries = dictionaries, beta = beta, center = center,
    null = null, design_scale = scale[-1L])
}

database_additive_effects = function(conn, working, key_source, fe, key_columns, groups,
                                      x, beta, center, n, resources) {
  if (length(fe) == 1L) {
    return(database_one_effect(conn, working, key_source, fe, key_columns[[1L]],
      groups[[1L]], x, beta, center, n, resources))
  }
  if (length(fe) == 2L) {
    database_two_effects(conn, working, key_source, fe, key_columns, groups,
      x, beta, center, n, resources)
  }
}
