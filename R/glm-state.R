# Capture only immutable fit context. Strategy flags are passed on each call,
# so a moment-cap fallback cannot leave a stale strategy in this callback.
dbglm_fe_state_builder = function(conn, base, x, link, create) {
  backend_info = detect_backend(conn)
  inline_sql = function(sql) {
    if (isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
      create(sql, materialize = TRUE)
    } else {
      paste0("(", dbglm_sql(conn, sql), ")")
    }
  }
  function(beta_val, c_eta_val, gamma_val = NULL, ops = NULL, info = NULL,
           use_moments = FALSE, collect_fe_moments = FALSE) {
    lookup_sub = NULL
    array_terms = if (collect_fe_moments && isTRUE(backend_info[["is_databricks"]]) &&
      !is.null(gamma_val) && !all(gamma_val == 0)) dbglm_databricks_array_terms(ops, info, gamma_val) else NULL
    if (is.null(gamma_val) || all(gamma_val == 0)) {
      eta = paste(c("__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val), sprintf("(%s)*%s", dbglm_number(beta_val), x)), collapse = " + ")
      st = if (use_moments) inline_sql(paste("SELECT *,", eta, "AS __dbglm_eta FROM", base)) else
        create(paste("SELECT *,", eta, "AS __dbglm_eta FROM", base), materialize = FALSE)
    } else if (!is.null(array_terms)) {
      eta = paste(c("s.__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val),
        sprintf("(%s)*s.%s", dbglm_number(beta_val), x), array_terms), collapse = " + ")
      st = inline_sql(paste0("SELECT s.*,CAST((", eta, ") AS DOUBLE) AS __dbglm_eta FROM ", base, " s"))
    } else {
      frame = dbglm_coefficient_frame(info, gamma_val)
      dictionary = dbglm_coefficient_lookup(conn, frame, create,
        inline = !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE)))
      lookup = dictionary[["source"]]
      lookup_sub = if (length(dictionary[["tables"]])) dictionary[["tables"]] else NULL
      ctes = dbglm_coefficient_ctes(lookup, length(ops), "eta")
      fe_terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
        if (v == "1.0") sprintf("COALESCE(eta_g%d.c1,0.0)", i) else sprintf("(%s)*COALESCE(eta_g%d.c1,0.0)", v, i)
      }, character(1))
      base_terms = c("s.__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val), sprintf("(%s)*s.%s", dbglm_number(beta_val), x), fe_terms)
      joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN eta_g", i, " ON s.", ops[[i]][["f"]], "=eta_g", i, ".__glm_group"), character(1))
      broadcast = if (isTRUE(backend_info[["is_databricks"]])) paste0("/*+ BROADCAST(", paste0("eta_g", seq_along(ops), collapse = ","), ") */ ") else ""
      st = inline_sql(paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast, "s.*, CAST((", paste(base_terms, collapse = " + "), ") AS DOUBLE) AS __dbglm_eta FROM ", base, " s", paste(joins, collapse = "")))
    }
    wk = if (use_moments) inline_sql(paste0(dbglm_working_sql(st, character(), numeric(), link, eta = "__dbglm_eta"), "SELECT * FROM working")) else
      dbglm_fe_working(conn, st, link, create)
    list(state = st, working = wk, lookup = lookup_sub)
  }
}

# Working-state construction and diagnostics for absorbed GLMs.
dbglm_fe_working = function(conn, table, link, create, target = NULL) {
  create(paste0(dbglm_working_sql(table, character(), numeric(), link, eta = "__dbglm_eta"), "SELECT * FROM working"), target = target)
}

dbglm_fe_evaluate = function(conn, working, x, fe, collect_fe_moments = FALSE, moment_rows = NULL) {
  cols = c(x, "z")
  terms = c(dbglm_deviance_term(conn), "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
    "SUM(CASE WHEN NOT isfinite(w) OR w<=0 OR NOT isfinite(eta) THEN 1 ELSE 0 END) AS bad",
    dbglm_moment_terms(cols, "w"))
  if (length(x)) {
    x_score = vapply(paste0(x, "*score"), function(e) dbglm_exact_sum(conn, e), character(1))
    terms = c(terms, sprintf("%s AS raw_s%d", x_score, seq_along(x)),
      sprintf("ABS(%s)/SQRT(SUM(w*%s*%s)) AS s%d", x_score, x, x, seq_along(x)))
  }
  if (collect_fe_moments && length(fe) >= 2L && !length(unlist(attr(fe, "slopes")))) {
    max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6)
    max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6)
    dbglm_validate_moment_limits(max_edges, max_entries)
    K = length(fe)
    # The initial projector bounds support on the owned immutable base table.
    # A working predictor preserves every key, so reuse that verified bound.
    # Without it, fetch only the cap plus the global row and an overflow row.
    known_bound = length(moment_rows) == 1L && is.finite(moment_rows) &&
      moment_rows >= 0 && moment_rows == floor(moment_rows) && moment_rows <= max_edges
    res = dbglm_read_fe_moments(conn, dbglm_fe_moments_sql(conn, working, fe,
      dbglm_fe_value_terms(conn, cols, "score"), terms, limit = if (known_bound) NULL else max_edges + 2),
    max_rows = if (known_bound) moment_rows + 1 else max_edges + 1)
    if (is.null(res)) {
      if (known_bound) stop("The fixed-effect moment support changed within fitting.", call. = FALSE)
      return(dbglm_fe_evaluate(conn, working, x, fe, collect_fe_moments = FALSE))
    }
    if (known_bound && nrow(res) != moment_rows + 1) {
      stop("The fixed-effect moment support changed within fitting.", call. = FALSE)
    }
    sets = as.integer(res[["__set"]])
    row = attr(res, "dbglm_globals")
    if (is.null(row)) row = res[sets == 0L, , drop = FALSE]
    fe_moments = NULL
    if (!is.na(row[["bad"]]) && row[["bad"]] == 0 && sum(sets != 0L) <= max_edges) {
      unpacked = dbglm_unpack_fe_moments(res, K, length(cols), 1L)
      if (!is.null(unpacked) && unpacked[["n"]] * max(1L, length(cols)) <= max_entries) fe_moments = unpacked
    }
    singles = sets >= 1L & sets <= K
    h = as.numeric(res[["__h"]][singles])
    # Also certify the Newton-scale group equation. A square-root information
    # score alone vanishes with a group's weight even when its predictor is
    # far from its finite solution. This uses the already collected moments.
    s_val = ifelse(is.finite(h) & h > 0 & is.finite(res[["__ex1"]][singles]),
      abs(res[["__ex1"]][singles]) / sqrt(h) / pmin(1, sqrt(h)), 0)
    result = dbglm_fe_diagnostics(row, x, if (length(s_val)) max(s_val) else 0)
    # Retain the explicit NULL when a bounded summary cannot be used.
    result["fe_moments"] = list(fe_moments)
    return(result)
  }
  fe_queries = character()
  for (i in seq_along(fe)) for (v in c("1.0", attr(fe, "slopes")[[i]])) {
    fe_queries = c(fe_queries, paste0("SELECT COALESCE(ABS(", dbglm_exact_sum(conn, paste0("score*(", v, ")")),
      ")/NULLIF(SQRT(SUM(w*(", v, ")*(", v, ")))*LEAST(1.0,SQRT(SUM(w))),0),0) AS s FROM ", working, " GROUP BY ", fe[i]))
  }
  if (length(fe_queries)) {
    terms = c(terms, paste0("(SELECT MAX(s) FROM (", paste(fe_queries, collapse = " UNION ALL "), ") t) AS fe_score_max"))
  }
  row = dbglm_get_query(conn, paste("SELECT", paste(terms, collapse = ","), "FROM", working))
  fe_score = if (length(fe_queries) && !is.na(row[["fe_score_max"]])) as.numeric(row[["fe_score_max"]]) else 0
  dbglm_fe_diagnostics(row, x, fe_score)
}

# Decode the same scalar diagnostics from fused moments and SQL projections.
dbglm_fe_diagnostics = function(row, x, fe_score) {
  score = if (length(x)) max(as.numeric(row[1, paste0("s", seq_along(x)), drop = FALSE])) else 0
  raw_score = if (length(x)) as.numeric(row[1, paste0("raw_s", seq_along(x)), drop = FALSE]) else numeric()
  raw_moments = dbglm_read_information(row, length(x) + 1L)
  list(deviance = row[["deviance"]], mu_min = row[["mu_min"]], mu_max = row[["mu_max"]], bad = row[["bad"]],
    score_max = max(score, fe_score), raw_score = raw_score, raw_moments = raw_moments,
    scales = sqrt(pmax(diag(raw_moments), .Machine[["double.xmin"]])))
}
