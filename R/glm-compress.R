dbglm_compress_base = function(conn, base, groups, family, center, create, drop,
                               weight_scale = 1) {
  groups = unique(groups)
  if (family == "binomial") return(dbglm_compress_binary(conn, base, groups, create, weight_scale))
  sum = function(expr) sql_stable_sum(conn, expr)
  centered = paste0("(__dbglm_y-(", dbglm_number(center), "))")
  mean = if (family == "gaussian") paste0(
    "(", dbglm_number(center), ")+",
    sum(paste0("__dbglm_a*", centered)), "/", sum("__dbglm_a")
  ) else
    paste0(sum("__dbglm_a*__dbglm_y"), "/", sum("__dbglm_a"))
  totals = paste0(
    "SELECT ", paste(c(
      groups,
      paste(sum("__dbglm_a"), "AS __dbglm_a"),
      paste(mean, "AS __dbglm_y")
    ), collapse = ","),
    " FROM ", base, if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ","))
  )
  means = create(totals)
  on.exit(drop(means), add = TRUE)
  delta = "(b.__dbglm_y-c.__dbglm_y)"
  within = switch(family,
    gaussian = paste0("b.__dbglm_a*", delta, "*", delta),
    poisson = paste0("2.0*b.__dbglm_a*(", sql_poisson_bd0("b.__dbglm_y", "c.__dbglm_y"), ")")
  )
  fields = c(paste0("c.", groups), "c.__dbglm_a", "c.__dbglm_y")
  aggregate = c(
    "COUNT(*) AS __dbglm_cell_n",
    paste(sum("b.__dbglm_a*b.__dbglm_a"), "AS __dbglm_a2"),
    paste(sum(paste0("b.__dbglm_a*b.__dbglm_a*", delta)), "AS __dbglm_a2delta"),
    paste(sum(paste0("b.__dbglm_a*b.__dbglm_a*", delta, "*", delta)), "AS __dbglm_a2delta2"),
    paste(sum(within), "AS __dbglm_dev_const")
  )
  join = if (length(groups)) paste0(
    " JOIN ", means, " c ON ",
    paste(paste0("b.", groups, "=c.", groups), collapse = " AND ")
  ) else paste0(" CROSS JOIN ", means, " c")
  create(paste0(
    "SELECT ", paste(c(fields, aggregate), collapse = ","), " FROM ", base, " b", join,
    " GROUP BY ", paste(fields, collapse = ",")
  ))
}

# A binary outcome needs only positive sums for its two outcome categories.
# Their squared-weight totals give exact HC1 moments without a second row pass.
dbglm_compress_binary = function(conn, base, groups, create, weight_scale) {
  sum = function(expr) sql_stable_sum(conn, expr)
  rounded = sql_round_even(paste0("__dbglm_a*(", dbglm_number(weight_scale), ")"))
  fields = c(
    groups, "COUNT(*) AS __dbglm_cell_n",
    paste(sum("__dbglm_a"), "AS __dbglm_a"),
    paste0(sum("__dbglm_a*__dbglm_y"), "/", sum("__dbglm_a"), " AS __dbglm_y"),
    paste(sum("CASE WHEN __dbglm_y=1 THEN __dbglm_a*__dbglm_a ELSE 0.0 END"), "AS __dbglm_cell_s1"),
    paste(sum("CASE WHEN __dbglm_y=0 THEN __dbglm_a*__dbglm_a ELSE 0.0 END"), "AS __dbglm_cell_s0"),
    paste(sum(rounded), "AS __dbglm_round_a"),
    paste(sum(paste0("(", rounded, ")*__dbglm_y")), "AS __dbglm_round_ay")
  )
  grouped = paste0(
    "SELECT ", paste(fields, collapse = ","), " FROM ", base,
    if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ","))
  )
  fields = c(
    groups, "__dbglm_cell_n", "__dbglm_a", "__dbglm_y",
    "__dbglm_cell_s1", "__dbglm_cell_s0",
    "0.0 AS __dbglm_dev_const", "__dbglm_round_a", "__dbglm_round_ay"
  )
  create(paste0("SELECT ", paste(fields, collapse = ","), " FROM (", grouped, ") binary_cells"))
}

dbglm_cell_score_square = function(link) {
  variance = "GREATEST(0.0,__dbglm_a2delta2+2.0*(__dbglm_y-mu)*__dbglm_a2delta+POWER(__dbglm_y-mu,2)*__dbglm_a2)"
  if (link %in% c("identity", "log")) return(variance)
  paste0("(__dbglm_cell_s1*POWER(1.0-mu,2)+__dbglm_cell_s0*POWER(mu,2))",
    "*POWER(d/(mu*(1.0-mu)),2)")
}

dbglm_cell_source = function(conn, base, groups, n, strategy, family, center,
                             create, drop, weight_scale, ratio, nmax, fe_sizes = numeric()) {
  if (!strategy %in% c("auto", "compress")) {
    return(list(table = base, enabled = FALSE, cells = n, groups = groups))
  }
  cap = if (strategy == "compress") min(n, nmax) else floor(min(n * ratio, nmax))
  decline = function() {
    if (strategy == "compress") stop("Compressed GLM exceeds the cell limit.", call. = FALSE)
    return(list(table = base, enabled = FALSE, cells = n, groups = groups))
  }
  if (cap < 1 || (length(fe_sizes) > 0L && max(fe_sizes) > cap)) return(decline())
  limit_clause = paste("LIMIT", format(cap + 1, scientific = FALSE, trim = TRUE))
  count = if (length(groups)) paste0(
    "SELECT COUNT(*) AS n FROM (SELECT DISTINCT ",
    paste(groups, collapse = ","), " FROM ", base, " ", limit_clause, ") cell_keys"
  ) else "SELECT 1 AS n"
  size = as.numeric(DBI::dbGetQuery(conn, count)[["n"]])
  if (size > cap) return(decline())
  table = dbglm_compress_base(conn, base, groups, family, center, create, drop, weight_scale)
  collected = DBI::dbGetQuery(conn, paste(
    "SELECT COUNT(*) AS cells, SUM(__dbglm_cell_n) AS observations FROM", table))
  if (!isTRUE(as.numeric(collected[["cells"]]) == size) ||
      !isTRUE(as.numeric(collected[["observations"]]) == n)) {
    stop("Compressed GLM cells do not reproduce the original sample counts.", call. = FALSE)
  }
  list(table = table, enabled = TRUE, cells = size, groups = groups)
}
