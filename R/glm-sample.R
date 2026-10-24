# GLM-only formula normalization. Keep OLS parsing and contrasts unchanged.
dbglm_formula = function(fml) {
  fml = Formula::Formula(fml)
  if (length(fml)[1] != 1 || length(fml)[2] > 2) stop("Use one outcome and at most one fixed-effect block.", call. = FALSE)
  state = new.env(parent = emptyenv())
  state[["factors"]] = list()
  state[["slopes"]] = list()
  replace = function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1]], as.name("factor"))) {
      if (length(e) != 2 || !is.symbol(e[[2]])) stop("Use factor(column) without additional arguments.", call. = FALSE)
      source = as.character(e[[2]])
      name = paste0("dbglm_factor_", match(source, unique(c(unlist(state[["factors"]]), source))), "__")
      state[["factors"]][[name]] = source
      return(as.name(name))
    }
    as.call(lapply(as.list(e), replace))
  }
  rhs = stats::formula(fml, lhs = 0, rhs = 1)[[2]]
  rhs = replace(rhs)
  fe = if (length(fml)[2] > 1) stats::formula(fml, lhs = 0, rhs = 2)[[2]] else NULL
  extract_slopes = function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1]], as.name("["))) {
      args = as.list(e)[-1]
      if (length(args) < 2 || length(args) > 3 || !all(vapply(args, is.symbol, logical(1)))) stop("Use group[trend] or group[trend,trend2].", call. = FALSE)
      group = as.character(args[[1]])
      state[["slopes"]][[group]] = unique(c(state[["slopes"]][[group]], vapply(args[-1], as.character, character(1))))
      return(args[[1]])
    }
    as.call(lapply(as.list(e), extract_slopes))
  }
  if (!is.null(fe)) fe = replace(extract_slopes(fe))
  lhs = stats::formula(fml, lhs = 1, rhs = 0)[[2]]
  normalized = stats::as.formula(call("~", lhs, if (is.null(fe)) rhs else call("|", rhs, fe)), env = baseenv())
  list(fml = normalized, factors = state[["factors"]], slopes = state[["slopes"]])
}

# Recursion operates on a snapshot and never rediscovers keys.
dbglm_prune = function(conn, table, fe, y, mode, create, drop, family = "binomial",
  flags = c("__dbglm_single", "__dbglm_infinite")) {
  rounds = 0L
  counts = c(singleton = 0, infinite_coef = 0)
  if (!length(fe) || mode == "none") return(list(table = table, rounds = rounds, counts = counts))
  repeat {
    # Per-FE group summaries are tiny; joining them avoids N-row window sorts.
    singleton = infinite = joins = character()
    for (k in seq_along(fe)) {
      f = fe[k]
      s = paste0("__dbglm_p", k)
      joins = c(joins, paste0(" JOIN (SELECT ", f, " AS __dbglm_g, COUNT(*) AS n, MIN(", y, ") AS lo, MAX(", y,
        ") AS hi FROM ", table, " GROUP BY ", f, ") ", s, " ON __dbglm_t.", f, " = ", s, ".__dbglm_g"))
      singleton = c(singleton, paste0(s, ".n = 1"))
      infinite = c(infinite, if (family == "poisson") paste0(s, ".hi = 0") else paste0(s, ".lo = ", s, ".hi"))
    }
    flagged = paste0("(SELECT __dbglm_t.*, (", paste(singleton, collapse = " OR "), ") AS ", flags[[1L]], ", (",
      paste(infinite, collapse = " OR "), ") AS ", flags[[2L]], " FROM ", table, " __dbglm_t",
      paste(joins, collapse = ""), ")")
    condition = switch(mode, perfect_fit = paste0("(", flags[[1L]], " OR ", flags[[2L]], ")"),
      singletons = flags[[1L]], infinite_coef = flags[[2L]])
    info = dbglm_get_query(conn, paste0("SELECT COUNT(*) FILTER (WHERE ", condition, ") AS removed, ",
      "COUNT(*) FILTER (WHERE ", condition, " AND ", flags[[1L]], ") AS single FROM ", flagged, " __dbglm_f"))
    if (info[["removed"]] == 0) break
    counts = counts + c(singleton = as.numeric(info[["single"]]), infinite_coef = as.numeric(info[["removed"]] - info[["single"]]))
    # Each round refers to its input once per FE. Inlining successive rounds
    # on Databricks grows the plan exponentially, even on a tiny sample.
    # Materialize only after a removal; the common no-pruning case adds no write.
    next_table = create(paste0("SELECT * EXCLUDE (", paste(flags, collapse = ","), ") FROM ", flagged, " __dbglm_f WHERE NOT ", condition), materialize = TRUE)
    drop(table)
    table = next_table
    rounds = rounds + 1L
  }
  list(table = table, rounds = rounds, counts = counts)
}

# fixest's default parameter convention is intentionally not exact dummy rank.
dbglm_fe_counts = function(conn, table, fe, p, kind, ssc, fe_names = NULL, cluster_var = NULL, sizes = NULL) {
  if (!length(fe)) return(list(model = p, covariance = p, sizes = integer(), nested = logical()))
  if (is.null(sizes)) {
    size_row = dbglm_get_query(conn, paste("SELECT",
      paste(sprintf("COUNT(DISTINCT %s) AS n%d", fe, seq_along(fe)), collapse = ", "),
      "FROM", table))
    sizes = vapply(size_row, as.numeric, numeric(1))
  }
  slope_df = sum(sizes * lengths(attr(fe, "slopes")))
  model = p + slope_df + sum(sizes) - length(fe) + 1
  nested = rep(FALSE, length(fe))
  explicit = if (length(fe_names) && length(cluster_var)) vapply(strsplit(fe_names, "^", fixed = TRUE), function(parts) cluster_var %in% parts, logical(1)) else rep(FALSE, length(fe))
  if (kind == "cluster" && any(explicit)) nested = explicit
  if (kind == "cluster" && !any(explicit)) {
    nest_queries = vapply(seq_along(fe), function(i) {
      paste0("SELECT ", i, " AS idx, MAX(n) AS max_n FROM (SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM ",
        table, " GROUP BY ", fe[i], ") t")
    }, character(1))
    nest_df = dbglm_get_query(conn, paste(nest_queries, collapse = " UNION ALL "))
    nested[nest_df[["idx"]]] = (nest_df[["max_n"]] == 1)
  }
  K = model
  if (ssc[["K.fixef"]] == "none") K = p + slope_df
  if (kind == "cluster" && ssc[["K.fixef"]] == "nonnested" && any(nested)) {
    K = if (all(nested)) p + slope_df + 1 else model - sum(sizes[nested] - 1)
    K = max(K, p + 1)
  }
  list(model = model, covariance = K, sizes = sizes, nested = nested)
}
