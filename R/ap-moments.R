# Exact weighted within transformation using the bounded FE solver shared with
# GLMs. Raw rows and key dictionaries stay in SQL, including string/BIGINT keys.
# The result obeys the existing AP table contract used by OLS and IV frames.
dbreg_moment_projections = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights", "cluster_var", "verbose"))
  backend = detect_backend(conn)
  if (!(backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) ||
      length(fe) < 2L || length(fe) > 30L) return(NULL)
  if (isTRUE(backend[["is_databricks"]])) validate_databricks_ap_backend(backend)

  if (identical(backend[["name"]], "duckdb") && !is.null(inputs[["retain"]])) {
    fast = dbreg_fwl_moments(inputs, x_sql, x_names, max_iter, tol)
    if (!is.null(fast)) return(fast)
  }
  limits = inputs[["projection_limits"]]
  if (is.null(limits)) limits = projection_limits()
  max_edges = limits[["summary_rows"]]
  max_entries = limits[["design_entries"]]
  max_cells = limits[["transfer_cells"]]
  validate_moment_limits(max_edges, max_entries)
  validate_moment_limits(max_cells, max_cells)
  vars = c(yvar, x_names)
  columns = paste0("dbreg_m_v", seq_along(vars))
  ids = unique(c(fe, cluster_var, inputs[["retained_key_safe"]]))
  id_columns = paste0("dbreg_m_id", seq_along(ids))
  fe_columns = paste0("dbreg_m_f", seq_along(fe))

  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources))
  manager = fit_table_manager(conn, resources, backend)
  create = manager[["create"]]
  drop = manager[["drop"]]
  w = if (is.null(weights)) "1.0" else sql_cast_double(weights, conn)
  source_sql = paste("SELECT", paste(c(
    sprintf("%s AS %s", ids, id_columns),
    sprintf("%s AS %s", sql_cast_double(c(yvar, x_sql), conn), columns),
    paste(w, "AS w")), collapse = ","), from_statement)
  fe_ids = id_columns[match(fe, ids)]
  stats = backend_get_query(conn, paste("SELECT", paste(c(
    "COUNT(*) AS n", "SUM(w) AS sw",
    sprintf("MAX(ABS(%s)) AS original%d", columns, seq_along(columns)),
    sprintf("COUNT(DISTINCT %s) AS n%d", fe_ids, seq_along(fe)),
    sprintf("SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS null%d", fe_ids, seq_along(fe)),
    sprintf("SUM(w*%s)/SUM(w) AS c%d", columns, seq_along(columns))), collapse = ","), "FROM (", source_sql, ") projection_preflight"))
  sizes = vapply(seq_along(fe), function(i) as.numeric(stats[[paste0("n", i)]]), numeric(1))
  centers = vapply(seq_along(columns), function(i) as.numeric(stats[[paste0("c", i)]]), numeric(1))
  nobs = as.numeric(stats[["n"]])
  if (nobs == 0 || any(!is.finite(centers)) || any(sizes < 1) ||
      any(vapply(seq_along(fe), function(i) as.numeric(stats[[paste0("null", i)]]) > 0, logical(1))) ||
      sum(sizes) * length(columns) > max_entries || !fe_dense_feasible(sizes, length(columns), limits)) {
    projection_abort("system_limit", "Fixed-effect system exceeds its projection budget.",
      list(fe_sizes = sizes, columns = length(columns), estimated_bytes = projection_memory(sizes, length(columns)), limits = limits))
  }
  source = create(source_sql)

  # Rank distinct keys, never the N-row source once per FE. Comparisons happen
  # in their SQL types, so IDs above 2^53 and arbitrary text remain distinct.
  dictionaries = vapply(seq_along(fe), function(i) {
    f = fe_ids[i]
    paste0(" JOIN (SELECT ", f, ",ROW_NUMBER() OVER (ORDER BY ", f, ") AS ", fe_columns[i],
      " FROM (SELECT DISTINCT ", f, " FROM ", source, ") keys) d", i, " ON s.", f, "=d", i, ".", f)
  }, character(1))
  fields = c(paste0("s.", id_columns), "s.w",
    sprintf("(s.%s-(%s)) AS %s", columns, dbglm_number(centers), columns),
    sprintf("d%d.%s", seq_along(fe), fe_columns))
  encoded = create(paste0("SELECT ", paste(fields, collapse = ","), " FROM ", source, " s", paste(dictionaries, collapse = "")))
  moments = query_fe_moments(conn, encoded, fe_columns, columns,
    max_edges = max_edges, max_entries = max_entries, max_cells = max_cells, fe_sizes = sizes)
  if (is.null(moments)) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.", list(limits = limits))
  # A stricter solve than the requested absolute group-mean tolerance leaves
  # room for SQL reconstruction roundoff. The SQL certificate below is final.
  target = min(1e-10, tol * sqrt(min(moments[["sys"]][["diag"]])) / 4)
  projected = project_fe_moments(conn, encoded, columns,
    c(id_columns, fe_columns), lapply(fe_columns, function(f) list(f = f, v = "1.0")),
    rep(1, length(columns)), target, list(fe_max_iter = max_iter), create, drop,
    new.env(parent = emptyenv()), precomputed = moments, key = NULL)
  if (is.null(projected)) return(NULL)

  certificate = dbreg_projection_certificate(conn, projected[["table"]], fe_columns, columns, nobs)
  if (!is.finite(certificate) || certificate > tol) {
    if (verbose) message("[dbreg] Bounded FE projection failed row certification; using iterative AP")
    projection_abort("certificate_failed", "Bounded FE projection failed row certification.", list(residual = certificate, target = tol))
  }
  result = create(paste("SELECT", paste(c(
    sprintf("%s AS %s", id_columns, ids), "w AS __w",
    sprintf("%s AS %s_tilde", columns, vars)), collapse = ","), "FROM", projected[["table"]]))
  # The returned physical table is now owned by the frame. Release every
  # intermediate here, even on Databricks, before the frame begins its passes.
  name = gsub('["`]', "", result)
  resources[["tables"]] = setdiff(resources[["tables"]], name)
  if (verbose) message("[dbreg] Bounded FE moments: ", projected[["cycles"]],
    " solve/refinement steps; certified max abs mean = ", sprintf("%.4e", certificate))
  list(table = result, base_table = NULL, all_tables = result,
    iterations = projected[["cycles"]], max_abs = certificate, method = "bounded_moments",
    transfer_cells = moments[["transfer_cells"]], estimated_bytes = projection_memory(sizes, length(columns)),
    nobs = nobs, n_fe_levels = sizes,
    original_scale = setNames(vapply(seq_along(columns), function(i) as.numeric(stats[[paste0("original", i)]]), numeric(1)), vars))
}

# Certify the actual transformed SQL rows, rather than only the small system.
# One grouping-sets pass checks group equations, finite values and row count.
dbreg_projection_certificate = function(conn, table, fe, columns, nobs) {
  means = sprintf("ABS(SUM(w*%s)/SUM(w)) AS r%d", columns, seq_along(columns))
  values = c("w", columns)
  invalid = paste(c("w <= 0", sprintf("%s IS NULL OR NOT isfinite(%s)", values, values)), collapse = " OR ")
  sql = paste0("SELECT ", paste(c(sprintf("MAX(r%d) AS r%d", seq_along(columns), seq_along(columns)),
    "MAX(bad) AS bad", "SUM(n) AS count_sum"), collapse = ","),
    " FROM (SELECT ", paste(c(means, "COUNT(*) AS n", paste0("SUM(CASE WHEN ", invalid, " THEN 1 ELSE 0 END) AS bad")), collapse = ","),
    " FROM ", table, " GROUP BY GROUPING SETS (", paste0("(", fe, ")", collapse = ","), ")) certified")
  row = backend_get_query(conn, sql)
  if (is.na(row[["bad"]]) || row[["bad"]] > 0 ||
      !isTRUE(as.numeric(row[["count_sum"]]) == nobs * length(fe))) return(Inf)
  max(vapply(row[paste0("r", seq_along(columns))], as.numeric, numeric(1)))
}

# General linear nuisance projection, including group-specific slopes. The
# shared projector only solves weighted group equations; no IRLS/family is used.
dbreg_varying_projections = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights", "cluster_var"))
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources), add = TRUE)
  manager = fit_table_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  ids = unique(c(fe, cluster_var, inputs[["retained_key_safe"]]))
  id_names = paste0("__dbreg_vary_id", seq_along(ids))
  columns = paste0("__dbreg_vary_v", seq_along(c(yvar, x_names)))
  trends = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  trend_names = stats::setNames(paste0("__dbreg_vary_t", seq_along(trends)), trends)
  fields = c("ROW_NUMBER() OVER () AS __dbglm_key", paste(ids, "AS", id_names),
    paste(sql_cast_double(c(yvar, x_sql), conn), "AS", columns),
    paste(sql_cast_double(trends, conn), "AS", unname(trend_names)),
    paste(if (is.null(weights)) "1.0" else sql_cast_double(weights, conn), "AS w"))
  source = create(paste("SELECT", paste(fields, collapse = ","), from_statement))
  weight_max = as.numeric(backend_get_query(conn, paste("SELECT MAX(w) AS m FROM", source))[["m"]])
  if (!is.finite(weight_max) || weight_max <= 0) stop("No positive-weight observations remain.", call. = FALSE)
  weight_scale = if (weight_max < 1 || weight_max > 2) 2^floor(log2(weight_max)) else 1
  w = paste0("w/(", dbglm_number(weight_scale), ")")
  centers = as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf(
    "SUM((%s)*%s)/SUM(%s) AS c%d", w, columns, w, seq_along(columns)), collapse = ","), "FROM", source))[1, ])
  coded_fe = paste0("__dbreg_vary_f", seq_along(fe))
  attr(coded_fe, "slopes") = lapply(inputs[["fe_slopes"]], function(vars) unname(trend_names[vars]))
  ranked = paste0("(", sql_rank_fixed_effects(conn, source, as.list(id_names[match(fe, ids)]), coded_fe), ")")
  encoded = create(paste("SELECT", paste(c("__dbglm_key", coded_fe, unname(trend_names),
    sprintf("(%s-(%s)) AS %s", columns, dbglm_number(centers), columns), paste(w, "AS w")), collapse = ","), "FROM", ranked))
  control = dbglm_control(fe_tol = min(tol, 1e-12), fe_max_iter = max_iter)
  projected = solve_partial_wls(conn, encoded, columns, coded_fe, control, create, drop, response = FALSE)
  result = create(paste0("SELECT ", paste(c(paste0("s.", id_names, " AS ", ids),
    paste0("p.", columns, " AS ", c(yvar, x_names), "_tilde"), "s.w AS __w"), collapse = ","),
    " FROM ", projected[["table"]], " p JOIN ", source, " s ON p.__dbglm_key=s.__dbglm_key"))
  name = gsub('["`]', "", result)
  resources[["tables"]] = setdiff(resources[["tables"]], name)
  list(table = result, base_table = NULL, all_tables = result,
    iterations = projected[["cycles"]], max_abs = projected[["residual"]], method = projected[["method"]])
}

# Read-only FWL: aggregate raw centered cross-products and FE information,
# solve the compact group equations, and subtract their projection in R.
dbreg_fwl_moments = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights"))
  limits = inputs[["projection_limits"]]
  if (is.null(limits)) limits = projection_limits()
  vars = c(yvar, x_names)
  taken = unique(c(fe, vars, weights, inputs[["cluster_var"]], inputs[["retained_key_safe"]], "w", "__w"))
  columns = safe_identifiers(paste0("dbreg_fwl_v", seq_along(vars)), taken = taken)
  keys = safe_identifiers(paste0("dbreg_fwl_f", seq_along(fe)), taken = c(taken, columns, paste0(columns, "_tilde")))
  w = if (is.null(weights)) "1.0" else sql_cast_double(weights, conn)
  passthrough = unique(c(fe, inputs[["cluster_var"]], inputs[["retained_key_safe"]]))
  id_columns = safe_identifiers(paste0("dbreg_fwl_id", seq_along(passthrough)),
    taken = c(taken, columns, paste0(columns, "_tilde"), keys))
  fe_ids = id_columns[match(fe, passthrough)]
  output_columns = safe_identifiers(paste0(vars, "_tilde"), taken = passthrough)
  weight_column = safe_identifiers("dbreg_fwl_weight__", taken = c(passthrough, output_columns))
  raw = paste("SELECT", paste(c(paste(passthrough, "AS", id_columns), paste(sql_cast_double(c(yvar, x_sql), conn), "AS", columns),
    paste(w, "AS w")), collapse = ","), from_statement)
  preflight = backend_get_query(conn, paste("SELECT", paste(c("COUNT(*) AS n", "MIN(w) AS wmin", "MAX(w) AS wmax",
    sprintf("COUNT(DISTINCT %s) AS size%d", fe_ids, seq_along(fe)),
    sprintf("MIN(typeof(%s)) AS type%d", fe_ids, seq_along(fe)),
    sprintf(paste0("MAX(CASE WHEN typeof(%1$s) IN ('TINYINT','SMALLINT','INTEGER','BIGINT',",
      "'UTINYINT','USMALLINT','UINTEGER','UBIGINT','DOUBLE','FLOAT','DATE','BOOLEAN') ",
      "THEN 64 ELSE OCTET_LENGTH(encode(CAST(%1$s AS VARCHAR))) END) AS bytes%2$d"), fe_ids, seq_along(fe)),
    sprintf("FAVG(%s) AS center%d", columns, seq_along(columns)),
    sprintf("MAX(ABS(%s)) AS maximum%d", columns, seq_along(columns))), collapse = ","), "FROM (", raw, ") raw"))
  nobs = as.numeric(preflight[["n"]])
  sizes = vapply(seq_along(fe), function(i) as.numeric(preflight[[paste0("size", i)]]), numeric(1))
  centers = vapply(seq_along(columns), function(i) as.numeric(preflight[[paste0("center", i)]]), numeric(1))
  maximum = vapply(seq_along(columns), function(i) as.numeric(preflight[[paste0("maximum", i)]]), numeric(1))
  if (!nobs || any(!is.finite(c(centers, maximum))) || any(sizes < 1) ||
      sum(sizes) * length(columns) > limits[["design_entries"]] ||
      !fe_dense_feasible(sizes, length(columns), limits)) return(NULL)
  wmax = as.numeric(preflight[["wmax"]])
  weight_scale = 2^floor(log2(wmax))
  if (!is.finite(weight_scale) || weight_scale <= 0 || as.numeric(preflight[["wmin"]]) / weight_scale <= 0) return(NULL)
  key_bytes = vapply(seq_along(fe), function(i) as.numeric(preflight[[paste0("bytes", i)]]), numeric(1))
  row_bytes = 4 * (sum(sort(key_bytes, decreasing = TRUE)[1:2]) + 128 + 8 * (length(columns) + 3L))
  system_bytes = projection_memory(sizes, length(columns))
  max_rows = min(limits[["summary_rows"]], floor((limits[["memory_bytes"]] - system_bytes) / row_bytes))
  if (!is.finite(max_rows) || max_rows < 1) return(NULL)
  key_sql = vapply(seq_along(fe), function(i) {
    type = tolower(as.character(preflight[[paste0("type", i)]]))
    if (type %in% c("double", "float", "real")) {
      paste0("CASE WHEN ", fe_ids[[i]], "=0 THEN CAST(0.0 AS DOUBLE) ELSE ", fe_ids[[i]], " END")
    } else fe_ids[[i]]
  }, character(1))
  source = paste0("(SELECT ", paste(c(id_columns, paste(key_sql, "AS", keys),
    sprintf("(%s-(%s)) AS %s", columns, dbglm_number(centers), columns),
    paste0("w/(", dbglm_number(weight_scale), ") AS w")), collapse = ","), " FROM (", raw, ") raw)")
  global_terms = fwl_information_terms(conn, columns)
  moments = query_fe_moments(conn, source, keys, columns, max_edges = max_rows,
    max_entries = limits[["design_entries"]], max_cells = limits[["transfer_cells"]], fe_sizes = sizes,
    global_terms = global_terms, text_keys = TRUE)
  if (is.null(moments)) return(NULL)
  global = lapply(moments[["globals"]], as.numeric)
  if (!identical(as.numeric(global[["fwl_n"]]), nobs) ||
      !identical(as.numeric(vapply(moments[["info"]], function(info) length(info[["id"]]), integer(1))), sizes)) return(NULL)
  target = min(1e-10, tol * sqrt(min(moments[["sys"]][["diag"]])) / 4)
  cache = new.env(parent = emptyenv())
  projected = project_fe_moments(conn, source, columns, keys,
    lapply(keys, function(key) list(f = key, v = "1.0")), rep(1, length(columns)), target,
    list(fe_max_iter = max_iter), create = NULL, drop = NULL, cache = cache,
    precomputed = moments, need_table = FALSE)
  raw_moments = dbglm_read_information(global, length(columns))
  M = raw_moments - projected[["cross_reduction"]]
  M = (M + t(M)) / 2
  positive = diag(raw_moments) > 0
  recompute = any(diag(M)[positive] < 1e-8 * diag(raw_moments)[positive])
  if (any(!is.finite(M))) return(NULL)
  need_rows = retain_sample(inputs[["retain"]]) || recompute || !identical(inputs[["vcov_type_req"]], "iid")
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources))
  table = NULL
  if (need_rows) {
    manager = fit_table_manager(conn, resources)
    rows = project_fe_moments(conn, source, columns, id_columns,
      lapply(keys, function(key) list(f = key, v = "1.0")), rep(1, length(columns)), target,
      list(fe_max_iter = max_iter), manager[["create"]], manager[["drop"]], cache,
      precomputed = moments, materialize = FALSE)
    projected_source = rows[["raw_projected"]]
    if (recompute) {
      actual = backend_get_query(conn, paste("SELECT",
        paste(fwl_information_terms(conn, paste0(columns, "_tilde")), collapse = ","),
        "FROM", projected_source))
      M = dbglm_read_information(lapply(actual, as.numeric), length(columns))
      if (any(!is.finite(M))) return(NULL)
    }
    table = paste0("(SELECT ", paste(c(paste(id_columns, "AS", passthrough),
      paste0("w*(", dbglm_number(weight_scale), ") AS ", weight_column),
      paste(paste0(columns, "_tilde"), "AS", output_columns)), collapse = ","),
      " FROM ", projected_source, ")")
  }
  row = list(n_total = nobs, sum_w = global[["fwl_sw"]], sum_wy = 0, sum_wy_sq = M[1L, 1L])
  for (i in seq_along(fe)) row[[paste0("n_fe", i)]] = sizes[[i]]
  for (i in seq_along(x_names)) {
    row[[sprintf("sum_w%d_y", i)]] = M[i + 1L, 1L]
    for (j in seq_len(i)) row[[sprintf("sum_w%d_%d", j, i)]] = M[j + 1L, i + 1L]
  }
  for (i in seq_along(c(x_names, yvar))) {
    row[[paste0("dbreg_moment_center_", i)]] = 0
    row[[paste0("dbreg_moment_original_", i)]] = maximum[match(c(x_names, yvar)[[i]], vars)]
  }
  row[["dbreg_moment_weight_scale"]] = weight_scale
  tables = resources[["tables"]]
  resources[["tables"]] = character()
  list(table = table, base_table = NULL, all_tables = tables, inline = TRUE, moments = as.data.frame(row),
    cols = setNames(output_columns, vars), weight_column = weight_column,
    iterations = projected[["cycles"]], max_abs = projected[["residual"]] / sqrt(min(moments[["sys"]][["diag"]])),
    method = "bounded_moments", observation_tables = 0L, transfer_cells = moments[["transfer_cells"]],
    estimated_bytes = system_bytes + moments[["rows"]] * row_bytes, nobs = nobs, n_fe_levels = sizes)
}

# Raw and projected moment queries use the same stable scalar aggregations.
fwl_information_terms = function(conn, columns) {
  terms = c("COUNT(*) AS fwl_n", paste(sql_stable_sum(conn, "w"), "AS fwl_sw"))
  for (i in seq_along(columns)) for (j in seq_len(i)) {
    terms = c(terms, paste(sql_stable_sum(conn,
      paste0("w*(", columns[[i]], ")*(", columns[[j]], ")")), "AS", sprintf("h_%d_%d", i, j)))
  }
  terms
}
