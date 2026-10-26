# Cell summaries use centered within-cell outcome moments. They preserve the
# original likelihood and individual-row sandwich meat at a common predictor.
dbglm_compress_base = function(conn, base, groups, family, center, create, drop,
                               weight_scale = 1) {
  groups = unique(groups)
  if (family == "binomial") return(dbglm_compress_binary(conn, base, groups, create, weight_scale))
  sum = function(expr) sql_stable_sum(conn, expr)
  centered = paste0("(__dbglm_y-(", dbglm_number(center), "))")
  mean = if (family == "gaussian") paste0("(", dbglm_number(center), ")+",
    sum(paste0("__dbglm_a*", centered)), "/", sum("__dbglm_a")) else
      paste0(sum("__dbglm_a*__dbglm_y"), "/", sum("__dbglm_a"))
  totals = paste0("SELECT ", paste(c(groups, "MIN(__dbglm_key) AS __dbglm_key",
    paste(sum("__dbglm_a"), "AS __dbglm_a"),
    paste(mean, "AS __dbglm_y")), collapse = ","),
    " FROM ", base, if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ",")))
  means = create(totals)
  on.exit(drop(means), add = TRUE)
  delta = "(b.__dbglm_y-c.__dbglm_y)"
  within = switch(family,
    gaussian = paste0("b.__dbglm_a*", delta, "*", delta),
    poisson = paste0("2.0*b.__dbglm_a*(CASE WHEN b.__dbglm_y=0 THEN c.__dbglm_y ELSE ",
      "b.__dbglm_y*LN(b.__dbglm_y/c.__dbglm_y)-b.__dbglm_y+c.__dbglm_y END)"))
  fields = c(paste0("c.", groups), "c.__dbglm_key", "c.__dbglm_a", "c.__dbglm_y")
  aggregate = c("COUNT(*) AS __dbglm_cell_n",
    paste(sum("b.__dbglm_a*b.__dbglm_a"), "AS __dbglm_a2"),
    paste(sum(paste0("b.__dbglm_a*b.__dbglm_a*", delta)), "AS __dbglm_a2delta"),
    paste(sum(paste0("b.__dbglm_a*b.__dbglm_a*", delta, "*", delta)), "AS __dbglm_a2delta2"),
    paste(sum(within), "AS __dbglm_dev_const"))
  join = if (length(groups)) paste0(" JOIN ", means, " c ON ",
    paste(paste0("b.", groups, "=c.", groups), collapse = " AND ")) else paste0(" CROSS JOIN ", means, " c")
  create(paste0("SELECT ", paste(c(fields, aggregate), collapse = ","), " FROM ", base, " b", join,
    " GROUP BY ", paste(fields, collapse = ",")))
}

# A binary outcome needs only positive sums for its two outcome categories.
# Their squared-weight totals give exact HC1 moments without a second row pass.
dbglm_compress_binary = function(conn, base, groups, create, weight_scale) {
  sum = function(expr) sql_stable_sum(conn, expr)
  rounded = sql_round_even(paste0("__dbglm_a*(", dbglm_number(weight_scale), ")"))
  fields = c(groups, "MIN(__dbglm_key) AS __dbglm_key", "COUNT(*) AS __dbglm_cell_n",
    paste(sum("__dbglm_a"), "AS __dbglm_a"),
    paste0(sum("__dbglm_a*__dbglm_y"), "/", sum("__dbglm_a"), " AS __dbglm_y"),
    paste(sum("CASE WHEN __dbglm_y=1 THEN __dbglm_a*__dbglm_a ELSE 0.0 END"), "AS __dbglm_cell_s1"),
    paste(sum("CASE WHEN __dbglm_y=0 THEN __dbglm_a*__dbglm_a ELSE 0.0 END"), "AS __dbglm_cell_s0"),
    paste(sum(rounded), "AS __dbglm_round_a"),
    paste(sum(paste0("(", rounded, ")*__dbglm_y")), "AS __dbglm_round_ay"))
  grouped = paste0("SELECT ", paste(fields, collapse = ","), " FROM ", base,
    if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ",")))
  fields = c(groups, "__dbglm_key", "__dbglm_cell_n", "__dbglm_a", "__dbglm_y",
    "__dbglm_cell_s1+__dbglm_cell_s0 AS __dbglm_a2",
    "__dbglm_cell_s1*(1.0-__dbglm_y)-__dbglm_cell_s0*__dbglm_y AS __dbglm_a2delta",
    "__dbglm_cell_s1*POWER(1.0-__dbglm_y,2)+__dbglm_cell_s0*POWER(__dbglm_y,2) AS __dbglm_a2delta2",
    "0.0 AS __dbglm_dev_const", "__dbglm_round_a", "__dbglm_round_ay")
  create(paste0("SELECT ", paste(fields, collapse = ","), " FROM (", grouped, ") binary_cells"))
}

sql_round_even = function(value) {
  paste0("(CASE WHEN ", value, "-FLOOR(", value, ")<0.5 THEN FLOOR(", value, ") ",
    "WHEN ", value, "-FLOOR(", value, ")>0.5 THEN CEIL(", value, ") ",
    "WHEN MOD(FLOOR(", value, "),2.0)=0 THEN FLOOR(", value, ") ELSE CEIL(", value, ") END)")
}

dbglm_cell_score_square = function(link) {
  variance = "GREATEST(0.0,__dbglm_a2delta2+2.0*(__dbglm_y-mu)*__dbglm_a2delta+POWER(__dbglm_y-mu,2)*__dbglm_a2)"
  if (link %in% c("identity", "log")) return(variance)
  paste0("(", variance, ")*POWER(d/(mu*(1.0-mu)),2)")
}

dbglm_cell_source = function(conn, base, groups, n, strategy, family, center,
                             create, drop, weight_scale, ratio, nmax) {
  if (!strategy %in% c("auto", "compress")) {
    return(list(table = base, enabled = FALSE, cells = n, groups = groups))
  }
  if (strategy == "compress") {
    table = dbglm_compress_base(conn, base, groups, family, center, create, drop, weight_scale)
    size = as.numeric(backend_get_query(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]])
    return(list(table = table, enabled = TRUE, cells = size, groups = groups))
  }
  count = if (length(groups)) paste0("SELECT COUNT(*) AS n FROM (SELECT DISTINCT ",
    paste(groups, collapse = ","), " FROM ", base, ") cell_keys") else "SELECT 1 AS n"
  size = as.numeric(backend_get_query(conn, count)[["n"]])
  if (strategy == "auto" && (size / n > ratio || size > nmax)) {
    return(list(table = base, enabled = FALSE, cells = n, groups = groups))
  }
  table = dbglm_compress_base(conn, base, groups, family, center, create, drop, weight_scale)
  list(table = table, enabled = TRUE, cells = size, groups = groups)
}

# Only requested original-row inference needs this relation. Join the compact
# final projector to original rows; never refit or repeat weighted projection.
dbglm_expand_cell_fit = function(conn, original, fit, x, fe, groups, link, create) {
  if (!length(fe)) return(fit)
  rows = fit[["raw_projected"]]
  if (is.null(rows)) {
    fields = c("b.*", if (length(x)) paste0("p.", x, " AS ", x, "_tilde"))
    rows = paste0("(SELECT ", paste(fields, collapse = ","), " FROM ", fit[["working"]], " b",
      if (length(x)) paste(" JOIN", fit[["projected"]], "p USING (__dbglm_key)"), ")")
  }
  fields = c(groups, "eta AS dbglm_cell_eta__", "mu AS dbglm_cell_mu__",
    "d AS dbglm_cell_derivative__", if (length(x)) paste0(x, "_tilde"))
  lookup = create(paste("SELECT", paste(fields, collapse = ","), "FROM", rows))
  multiplier = if (link %in% c("identity", "log")) "1.0" else
    "c.dbglm_cell_derivative__/(c.dbglm_cell_mu__*(1.0-c.dbglm_cell_mu__))"
  score = paste0("b.__dbglm_a*(b.__dbglm_y-c.dbglm_cell_mu__)*(", multiplier, ") AS score")
  fields = c("b.*", "c.dbglm_cell_eta__ AS eta", "c.dbglm_cell_mu__ AS mu", score,
    if (length(x)) paste0("c.", x, "_tilde AS ", x, "_tilde"))
  fit[["raw_projected"]] = paste0("(SELECT ", paste(fields, collapse = ","), " FROM ", original,
    " b JOIN ", lookup, " c ON ", paste(paste0("b.", groups, "=c.", groups), collapse = " AND "), ")")
  fit
}
