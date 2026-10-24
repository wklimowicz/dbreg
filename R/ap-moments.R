# Exact weighted within transformation using the bounded FE solver shared with
# GLMs. Raw rows and key dictionaries stay in SQL, including string/BIGINT keys.
# The result obeys the existing AP table contract used by OLS and IV frames.
dbreg_moment_projections = function(inputs, x_sql, x_names, max_iter, tol) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe", "yvar", "weights", "cluster_var", "verbose"))
  backend = detect_backend(conn)
  if (!(backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) ||
      length(fe) < 2L || length(fe) > 30L) return(NULL)
  if (isTRUE(backend[["is_databricks"]])) validate_databricks_ap_backend(backend)

  max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6)
  max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6)
  max_cells = getOption("dbreg.ap_moment_max_cells", 8e6)
  dbglm_validate_moment_limits(max_edges, max_entries)
  dbglm_validate_moment_limits(max_cells, max_cells)
  vars = c(yvar, x_names)
  columns = paste0("dbreg_m_v", seq_along(vars))
  ids = unique(c(fe, cluster_var, inputs[["retained_key_safe"]]))
  id_columns = paste0("dbreg_m_id", seq_along(ids))
  fe_columns = paste0("dbreg_m_f", seq_along(fe))
  # Pair rows carry NULL RHS fields too. Bound the entire rectangular transfer
  # before fetching; the shared entries cap also bounds FE-by-design matrices.
  row_cap = min(max_edges, floor(max_cells / (length(columns) + 4L)))
  if (row_cap < 1L) return(NULL)

  resources = dbglm_resources()
  resources[["conn"]] = conn
  on.exit(dbglm_release_resources(resources))
  manager = dbglm_table_manager(conn, resources, backend)
  create = manager[["create"]]
  drop = manager[["drop"]]
  w = if (is.null(weights)) "1.0" else sql_cast_double(weights, conn)
  source = create(paste("SELECT", paste(c(
    sprintf("%s AS %s", ids, id_columns),
    sprintf("%s AS %s", sql_cast_double(c(yvar, x_sql), conn), columns),
    paste(w, "AS w")), collapse = ","), from_statement))
  fe_ids = id_columns[match(fe, ids)]
  stats = dbglm_get_query(conn, paste("SELECT", paste(c(
    "COUNT(*) AS n", "SUM(w) AS sw",
    sprintf("COUNT(DISTINCT %s) AS n%d", fe_ids, seq_along(fe)),
    sprintf("SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS null%d", fe_ids, seq_along(fe)),
    sprintf("SUM(w*%s)/SUM(w) AS c%d", columns, seq_along(columns))), collapse = ","), "FROM", source))
  sizes = vapply(seq_along(fe), function(i) as.numeric(stats[[paste0("n", i)]]), numeric(1))
  centers = vapply(seq_along(columns), function(i) as.numeric(stats[[paste0("c", i)]]), numeric(1))
  nobs = as.numeric(stats[["n"]])
  if (nobs == 0 || any(!is.finite(centers)) || any(sizes < 1) ||
      any(vapply(seq_along(fe), function(i) as.numeric(stats[[paste0("null", i)]]) > 0, logical(1))) ||
      sum(sizes) * length(columns) > max_entries || !dbglm_fe_dense_feasible(sizes)) return(NULL)

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
  moments = dbglm_query_fe_moments(conn, encoded, fe_columns, columns,
    max_edges = row_cap, max_entries = max_entries)
  if (is.null(moments)) return(NULL)
  # A stricter solve than the requested absolute group-mean tolerance leaves
  # room for SQL reconstruction roundoff. The SQL certificate below is final.
  target = min(1e-10, tol * sqrt(min(moments[["sys"]][["diag"]])) / 4)
  projected = tryCatch(dbglm_project_moments(conn, encoded, columns,
    c(id_columns, fe_columns), lapply(fe_columns, function(f) list(f = f, v = "1.0")),
    rep(1, length(columns)), target, list(fe_max_iter = max_iter), create, drop,
    new.env(parent = emptyenv()), precomputed = moments, key = NULL),
    error = function(e) {
      msg = conditionMessage(e)
      if (grepl("^GLM fixed-effect (information is singular|projection did not converge)", msg) ||
          identical(msg, "FE coefficient lookup exceeds the SQL transfer limit.")) return(NULL)
      stop(e)
    })
  if (is.null(projected)) return(NULL)

  certificate = dbreg_projection_certificate(conn, projected[["table"]], fe_columns, columns, nobs)
  if (!is.finite(certificate) || certificate > tol) {
    if (verbose) message("[dbreg] Bounded FE projection failed row certification; using iterative AP")
    return(NULL)
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
    iterations = projected[["cycles"]], max_abs = certificate, method = "bounded_moments")
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
  row = dbglm_get_query(conn, sql)
  if (is.na(row[["bad"]]) || row[["bad"]] > 0 ||
      !isTRUE(as.numeric(row[["count_sum"]]) == nobs * length(fe))) return(Inf)
  max(vapply(row[paste0("r", seq_along(columns))], as.numeric, numeric(1)))
}
