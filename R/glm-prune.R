# Recursion operates on a snapshot and never rediscovers keys.
dbglm_prune = function(
  conn, table, fe, y, mode, create, drop, family = "binomial",
  flags = c("__dbglm_single", "__dbglm_infinite")
) {
  rounds = 0L
  counts = c(singleton = 0, infinite_coef = 0)
  if (family == "gaussian") {
    mode = switch(mode, perfect_fit = "singletons", infinite_coef = "none", mode)
  }
  if (!length(fe) || mode == "none") {
    return(list(table = table, rounds = rounds, counts = counts))
  }
  backend_info = detect_backend(conn)
  repeat {
    # Per-FE group summaries are tiny; joining them avoids N-row window sorts.
    singleton = infinite = joins = character()
    for (k in seq_along(fe)) {
      f = fe[k]
      s = paste0("__dbglm_p", k)
      joins = c(joins, paste0(
        " JOIN (SELECT ", f, " AS __dbglm_g, COUNT(*) AS n, MIN(", y, ") AS lo, MAX(", y,
        ") AS hi FROM ", table, " GROUP BY ", f, ") ", s, " ON __dbglm_t.", f, " = ", s, ".__dbglm_g"
      ))
      singleton = c(singleton, paste0(s, ".n = 1"))
      infinite = c(infinite, if (family == "poisson") paste0(s, ".hi = 0") else paste0(s, ".lo = ", s, ".hi"))
    }
    flagged = paste0(
      "(SELECT __dbglm_t.*, (", paste(singleton, collapse = " OR "), ") AS ", flags[[1L]], ", (",
      paste(infinite, collapse = " OR "), ") AS ", flags[[2L]], " FROM ", table, " __dbglm_t",
      paste(joins, collapse = ""), ")"
    )
    condition = switch(mode,
      perfect_fit = paste0("(", flags[[1L]], " OR ", flags[[2L]], ")"),
      singletons = flags[[1L]],
      infinite_coef = flags[[2L]]
    )
    info = DBI::dbGetQuery(conn, paste0(
      "SELECT COUNT(*) FILTER (WHERE ", condition, ") AS removed, ",
      "COUNT(*) FILTER (WHERE ", condition, " AND ", flags[[1L]], ") AS single FROM ", flagged, " __dbglm_f"
    ))
    if (info[["removed"]] == 0) break
    counts = counts + c(singleton = as.numeric(info[["single"]]), infinite_coef = as.numeric(info[["removed"]] - info[["single"]]))
    # Each round refers to its input once per FE. Inlining successive rounds
    # on Databricks grows the plan exponentially, even on a tiny sample.
    # Materialize only after a removal; the common no-pruning case adds no write.
    next_table = create(paste0("SELECT ", sql_select_except(conn, flags, backend_info), " FROM ", flagged, " __dbglm_f WHERE NOT ", condition), materialize = TRUE)
    drop(table)
    table = next_table
    rounds = rounds + 1L
  }
  list(table = table, rounds = rounds, counts = counts)
}

dbglm_ranked_sample = function(conn, resources, from, columns, groups, fe, y,
                                 mode, family) {
  manager = fe_cg_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  ranks = vapply(seq_along(fe), function(i) {
    paste0("DENSE_RANK() OVER (ORDER BY ", paste(q(groups[[i]]), collapse = ","),
      ") AS ", fe[i])
  }, character(1))
  project = function(from) paste("SELECT", paste(c(q(columns), ranks), collapse = ","), from)
  first = create(project(from))
  pruned = dbglm_prune(conn, first, fe, y, mode,
    function(sql, materialize) create(sql), drop, family)
  table = pruned[["table"]]
  if (family != "gaussian" && mode %in% c("none", "singletons")) {
    test = if (family == "poisson") sprintf("MAX(%s) = 0", y) else sprintf("MAX(%s) = MIN(%s)", y, y)
    queries = paste0("SELECT 1 AS bad FROM ", table, " GROUP BY ", fe, " HAVING ", test)
    sql = paste0("SELECT COUNT(*) AS n FROM (", paste(queries, collapse = " UNION ALL "),
      ") glm_separated")
    if (DBI::dbGetQuery(conn, sql)[["n"]] > 0) {
      stop("GLM fixed-effect separation remains under the selected pruning policy.", call. = FALSE)
    }
  }
  if (pruned[["rounds"]]) {
    dense = create(project(paste("FROM", table)))
    drop(table)
    table = dense
  }
  list(from = paste("FROM", table), pruning = pruned[c("rounds", "counts")])
}

dbglm_parameter_counts = function(conn, sample, p, cluster, ssc) {
  sizes = sample[["sizes"]]
  fe = sample[["fe"]]
  slope_df = sum(sizes * lengths(sample[["slopes"]]))
  model = p + slope_df + fe_intercept_df(sizes)
  K = model
  if (ssc[["K.fixef"]] == "none") K = p + slope_df
  if (length(fe) && cluster && ssc[["K.fixef"]] == "nonnested") {
    families = sample[["public_groups"]]
    groups = sample[["default_cluster"]]
    groups = if (is.list(groups)) groups else if (length(groups)) as.list(groups) else families[1L]
    nested = vapply(families, function(fe_g) any(vapply(groups,
      function(cl) all(cl %in% fe_g), logical(1))), logical(1))
    if (!all(nested) && !all(vapply(groups, function(cl) any(vapply(families,
        function(fe_g) identical(cl, fe_g), logical(1))), logical(1)))) {
      queries = paste0("SELECT ", seq_along(fe), " AS idx, MAX(n) AS n FROM (SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM ",
        sample[["base"]], " GROUP BY ", fe, ") groups")
      rows = DBI::dbGetQuery(conn, paste(queries, collapse = " UNION ALL "))
      nested = nested | (rows[["n"]][match(seq_along(fe), rows[["idx"]])] == 1)
    }
    K = max(p + slope_df + 1, model - sum(sizes[nested] - 1))
  }
  list(model = model, covariance = K)
}
