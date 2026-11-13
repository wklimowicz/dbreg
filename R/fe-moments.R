# A known non-native Databricks connection cannot execute physical temp DDL.
validate_native_temp_query = function(conn, backend = detect_backend(conn)) {
  if (!isTRUE(backend[["is_databricks"]])) return(invisible(NULL))
  native = backend[["native_query"]]
  fields = unlist(attr(conn, "connectionString", exact = TRUE), use.names = TRUE)
  selected = fields[tolower(sub("^.*\\.", "", names(fields))) == "usenativequery"]
  if (length(selected)) native = !tolower(as.character(selected[[1L]])) %in% c("0", "false", "no")
  if (identical(native, FALSE)) stop("Databricks temporary tables require native query mode; reconnect with useNativeQuery = TRUE.", call. = FALSE)
  invisible(NULL)
}

# AP alternates fresh owned names; no replacement clause is needed or allowed.
create_projection_table_as = function(conn, table_name, select_sql, backend) {
  info = detect_backend(conn)
  if (isTRUE(info[["is_databricks"]])) {
    validate_native_temp_query(conn, info)
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", DBI::dbQuoteIdentifier(conn, table_name), "AS", select_sql))
  } else create_temp_table_as(conn, table_name, select_sql, backend)
}

# SQL owns observations and typed key dictionaries; R receives bounded summaries.
fe_projection_tables = function(conn, resources) {
  function(sql, columns) {
    force(sql)
    # Owned names are freshly randomized. A plain CREATE TEMP TABLE fails on
    # a collision instead of replacing, so no catalog lookup precedes it
    # (on Databricks every SHOW TABLES round trip costs a statement).
    name = basename(tempfile("dbreg_fe_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, sql_sample_table(conn, name, sql, columns))
    as.character(DBI::dbQuoteIdentifier(conn, name))
  }
}

# Store a coefficient lookup in an owned DuckDB table straight from R doubles,
# without parsing a literal for every value.
fe_duckdb_lookup_table = function(conn, coefficients, manager) {
  if (any(!is.finite(coefficients))) projection_abort("singular", "Nonfinite fixed-effect coefficients.")
  data = as.data.frame(cbind(as.numeric(seq_len(nrow(coefficients))), coefficients))
  names(data) = c("id", paste0("c", seq_len(ncol(coefficients))))
  view = basename(tempfile("dbreg_lookup_"))
  duckdb_register(conn, view, data)
  on.exit(duckdb_unregister(conn, view))
  manager[["create"]](paste("SELECT * FROM", DBI::dbQuoteIdentifier(conn, view)))
}

# DuckDB re-parses large inline VALUES subqueries on every downstream query;
# register large lookups into an owned temp table while preserving upload_bytes.
fe_lookup_source = function(conn, coefficients, budget, limits, create, snapshot = create, columns = NULL) {
  if (identical(detect_backend(conn)[["name"]], "duckdb") &&
      nrow(coefficients) * (ncol(coefficients) + 1L) > 1024L) {
    bytes = 8 * nrow(coefficients) * (ncol(coefficients) + 1L) + 1024
    if (budget[["upload_bytes"]] + bytes > min(limits[["upload_bytes"]], limits[["memory_bytes"]] / 8)) {
      projection_abort("upload_limit", "Fixed-effect coefficients exceed the SQL upload budget.", list(bytes = bytes, limits = limits))
    }
    budget[["upload_bytes"]] = budget[["upload_bytes"]] + bytes
    creator = if (is.null(columns)) snapshot else function(sql) snapshot(sql, columns)
    return(fe_duckdb_lookup_table(conn, coefficients, list(create = creator)))
  }
  sql = fe_coefficient_sql(conn, coefficients, budget, limits)
  if (is.null(columns)) create(sql) else create(sql, columns)
}

# Coefficients use bounded dense ordinal keys. Original typed keys stay in SQL.
fe_coefficient_sql = function(conn, coefficients, budget, limits) {
  if (any(!is.finite(coefficients))) projection_abort("singular", "Nonfinite fixed-effect coefficients.")
  values = cbind(seq_len(nrow(coefficients)), coefficients)
  columns = c("id", paste0("c", seq_len(ncol(coefficients))))
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    rows = apply(values, 1L, function(row) paste(sprintf("%.17g", row), collapse = ","))
    payload = paste0("[", paste0("[", rows, "]", collapse = ","), "]")
    fields = c("CAST(item[0] AS BIGINT) AS id", paste0("item[", seq_len(ncol(coefficients)), "] AS ", columns[-1L]))
    sql = paste0("SELECT ", paste(fields, collapse = ","), " FROM explode(from_json('", payload, "','ARRAY<ARRAY<DOUBLE>>')) AS lookup(item)")
  } else {
    rows = apply(values, 1L, function(row) paste(sql_cast_double(sprintf("%.17e", row), conn), collapse = ","))
    sql = paste0("SELECT * FROM (VALUES ", paste0("(", rows, ")", collapse = ","), ") AS lookup(", paste(columns, collapse = ","), ")")
  }
  bytes = nchar(sql, type = "bytes") + 1024
  if (budget[["upload_bytes"]] + bytes > min(limits[["upload_bytes"]], limits[["memory_bytes"]] / 8)) {
    projection_abort("upload_limit", "Fixed-effect coefficients exceed the SQL upload budget.", list(bytes = bytes, limits = limits))
  }
  budget[["upload_bytes"]] = budget[["upload_bytes"]] + bytes
  sql
}

# Certificate shortfall: group-mean residuals in units of max(1, column scale)
# (absolute for unit-scale columns, relative for e.g. dollar outcomes) and
# relative to each column's residual RMS.
fe_certificate_shortfall = function(certificate, tol) {
  value = max(certificate[["max_level"]], certificate[["max_relative"]]) / tol
  if (is.finite(value)) value else Inf
}

fe_projection_certificate = function(conn, table, fe, columns, nobs, original_scale) {
  sum_sql = function(x) if (detect_backend(conn)[["name"]] == "duckdb") paste0("FSUM(", x, ")") else paste0("SUM(", x, ")")
  means = paste0("ABS(", sum_sql(paste0("w*", columns)), ")/", sum_sql("w"), " AS r", seq_along(columns))
  invalid = paste(c("w<=0", paste0(c("w", columns), " IS NULL OR isnan(", c("w", columns), ") OR ABS(", c("w", columns), ")>1.7976931348623157e308")), collapse = " OR ")
  sql = paste0("SELECT ", paste(c(paste0("MAX(r", seq_along(columns), ") AS r", seq_along(columns)), "MAX(bad) AS bad", "SUM(n) AS n"), collapse = ","),
    " FROM (SELECT ", paste(c(means, "COUNT(*) AS n", paste0("SUM(CASE WHEN ", invalid, " THEN 1 ELSE 0 END) AS bad")), collapse = ","),
    " FROM ", table, " GROUP BY GROUPING SETS (", paste0("(", fe, ")", collapse = ","), ")) certified")
  row = DBI::dbGetQuery(conn, sql)
  if (!isTRUE(as.numeric(row[["bad"]]) == 0) || !isTRUE(as.numeric(row[["n"]]) == nobs * length(fe))) return(list(max_abs = Inf, max_level = Inf, max_relative = Inf))
  residual = vapply(row[paste0("r", seq_along(columns))], as.numeric, numeric(1))
  units = sprintf("%.17e", pmax(original_scale, .Machine[["double.xmin"]]))
  rms_terms = paste0("(", units, ")*SQRT(", sum_sql(paste0("w*POWER(", columns, "/(", units, "),2)")), "/", sum_sql("w"), ") AS r", seq_along(columns))
  rms = DBI::dbGetQuery(conn, paste("SELECT", paste(rms_terms, collapse = ","), "FROM", table))
  rms = vapply(rms, as.numeric, numeric(1))
  if (any(!is.finite(rms))) return(list(max_abs = Inf, max_level = Inf, max_relative = Inf))
  identified = rms > 128 * .Machine[["double.eps"]] * original_scale
  list(max_abs = max(residual), max_level = max(residual / pmax(1, original_scale)), max_relative = if (any(identified)) max(residual[identified] / rms[identified]) else 0)
}

# Prepare one native-key source shared by bounded moments and SQL CG.
prepare_fe_projection = function(inputs, x_sql, x_names, materialize = TRUE, coefficient_cg = FALSE) {
  conn = inputs[["conn"]]
  fe = inputs[["fe"]]
  backend = detect_backend(conn)
  if (!(backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) || length(fe) < (if (length(unlist(inputs[["fe_slopes"]]))) 1L else 2L) || (!coefficient_cg && length(fe) > 30L)) return(NULL)
  validate_native_temp_query(conn, backend)
  limits = inputs[["projection_limits"]]
  if (is.null(limits)) limits = projection_limits()
  resources = fit_resources()
  resources[["conn"]] = conn
  success = FALSE
  on.exit(if (!success) release_fit_resources(resources))
  snapshot = fe_projection_tables(conn, resources)
  # Every expression reads the same owned source. FWL keeps transformed rows
  # inline; the ordinary projection consumer can still request physical rows.
  create = if (materialize) snapshot else function(sql, columns) paste0("(", sql, ")")
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  sum_sql = function(x) if (backend[["name"]] == "duckdb") paste0("FSUM(", x, ")") else paste0("SUM(", x, ")")
  vars = c(inputs[["yvar"]], x_names)
  columns = paste0("v", seq_along(vars))
  trends = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  trend_columns = setNames(if (length(trends)) paste0("t", seq_along(trends)) else character(), trends)
  ids = projection_id_columns(inputs)
  id_columns = paste0("id", seq_along(ids))
  fe_ids = id_columns[match(fe, ids)]
  fe_columns = paste0("f", seq_along(fe))
  w = if (is.null(inputs[["weights"]])) "1.0" else sql_cast_double(q(inputs[["weights"]]), conn)
  source = snapshot(paste("SELECT", paste(c(paste(q(ids), "AS", id_columns),
    paste(sql_cast_double(c(q(inputs[["yvar"]]), x_sql), conn), "AS", columns), paste(w, "AS raw_w"), if (length(trends)) paste(sql_cast_double(q(trends), conn), "AS", unname(trend_columns))), collapse = ","), inputs[["from_statement"]]), c(id_columns, columns, "raw_w", unname(trend_columns)))
  stats = DBI::dbGetQuery(conn, paste("SELECT", paste(c("COUNT(*) AS n", "MAX(raw_w) AS sw",
    paste0("MIN(", columns, ") AS lo", seq_along(columns)),
    paste0("MAX(", columns, ") AS hi", seq_along(columns)),
    if (length(trends)) paste0("MAX(ABS(", trend_columns, ")) AS ts", seq_along(trends)),
    paste0("COUNT(DISTINCT ", fe_ids, ") AS n", seq_along(fe)),
    paste0("SUM(CASE WHEN ", fe_ids, " IS NULL THEN 1 ELSE 0 END) AS bad", seq_along(fe))), collapse = ","), "FROM", source))
  sizes = vapply(stats[paste0("n", seq_along(fe))], as.numeric, numeric(1))
  lower = vapply(stats[paste0("lo", seq_along(columns))], as.numeric, numeric(1))
  upper = vapply(stats[paste0("hi", seq_along(columns))], as.numeric, numeric(1))
  anchors = lower/2 + upper/2
  original_scale = upper/2 - lower/2
  scale = pmax(original_scale, .Machine[["double.xmin"]])
  nobs = as.numeric(stats[["n"]]); weight_scale = as.numeric(stats[["sw"]])
  if (!is.finite(nobs) || nobs < 1 || any(!is.finite(scale)) || !is.finite(weight_scale) || weight_scale <= 0 ||
      any(sizes < 1) || any(vapply(stats[paste0("bad", seq_along(fe))], as.numeric, numeric(1)) > 0)) {
    projection_abort("singular", "Fixed-effect snapshot has invalid values or keys.")
  }
  # CG budgets sixteen double-sized vector/buffer slots per coefficient cell;
  # observation-sized working tables stay in SQL. Moments also plans its Schur.
  feasible = if (coefficient_cg) 128 * sum(sizes) * length(columns) <= limits[["memory_bytes"]] else fe_dense_feasible(sizes, length(columns), limits)
  if (sum(sizes) > .Machine[["integer.max"]] - 1 || sum(sizes) * length(columns) > limits[["design_entries"]] || !feasible) {
    projection_abort("system_limit", "Fixed-effect system exceeds its projection budget.", list(fe_sizes = sizes, columns = length(columns), limits = limits))
  }
  trend_scale = if (length(trends)) pmax(vapply(stats[paste0("ts", seq_along(trends))], as.numeric, numeric(1)), .Machine[["double.xmin"]]) else numeric()
  names(trend_scale) = trends
  if (any(!is.finite(trend_scale))) projection_abort("singular", "Varying-slope snapshot has invalid numeric scales.")
  dictionaries = vapply(seq_along(fe), function(i) paste0(" JOIN (SELECT ", fe_ids[i], ",ROW_NUMBER() OVER (ORDER BY ", fe_ids[i], ") AS ", fe_columns[i],
    " FROM (SELECT DISTINCT ", fe_ids[i], " FROM ", source, ") keys) d", i, " ON s.", fe_ids[i], "=d", i, ".", fe_ids[i]), character(1))
  # Center before scaling so a large location cannot quantize a small,
  # representable within-coordinate variation. Every model includes FE
  # intercepts, so this affine change preserves the exact projected span.
  w_expr = paste0("s.raw_w/(", sprintf("%.17e", weight_scale), ")")
  units = paste0("(s.", columns, "-(", sprintf("%.17e", anchors), "))/(", sprintf("%.17e", scale), ")")
  # Compute weighted centers directly from source without FE dictionary joins.
  location = DBI::dbGetQuery(conn, paste("SELECT", paste(paste0(sum_sql(paste0(w_expr, "*", units)), "/", sum_sql(w_expr), " AS ", columns), collapse = ","), "FROM", source, "s"))
  centers = vapply(location, as.numeric, numeric(1))
  if (any(!is.finite(centers))) projection_abort("singular", "Nonfinite fixed-effect projection centers.")
  centered = create(paste0("SELECT ", paste(c(paste0("s.", id_columns), "s.raw_w", paste0(w_expr, " AS w"),
    paste0("(", units, "-(", sprintf("%.17e", centers), ")) AS ", columns),
    if (length(trends)) paste0("s.", trend_columns, " AS ", trend_columns),
    paste0("d", seq_along(fe), ".", fe_columns)), collapse = ","), " FROM ", source, " s", paste(dictionaries, collapse = "")),
    c(id_columns, "raw_w", "w", columns, fe_columns, unname(trend_columns)))
  success = TRUE
  mget(c("conn","fe","limits","resources","snapshot","create","q","vars","columns","ids","id_columns","fe_columns","sizes","original_scale","scale","nobs","weight_scale","source","centered","trend_columns","trend_scale"), envir = environment(), inherits = FALSE)
}

# This result implements the AP frame contract, including owned-table cleanup.
dbreg_moment_projections = function(inputs, x_sql, x_names, max_iter, tol, materialize = TRUE, two_family_pcg = FALSE) {
  if (two_family_pcg) {
    if (length(inputs[["fe"]]) != 2L || detect_backend(inputs[["conn"]])[["name"]] != "duckdb") return(NULL)
    # A per-call copy selects the shared sparse Schur operator. User transfer,
    # cross-block, design, upload and memory caps stay authoritative.
    inputs = list2env(as.list(inputs), parent = emptyenv())
    if (is.null(inputs[["projection_limits"]])) inputs[["projection_limits"]] = projection_limits()
    inputs[["projection_limits"]][["dense_levels"]] = 0
  }
  prepared = prepare_fe_projection(inputs, x_sql, x_names, materialize)
  if (is.null(prepared)) return(NULL)
  env2env(prepared, environment(), names(prepared))
  success = FALSE
  on.exit(if (!success) release_fit_resources(resources))
  budget = list2env(list(rows = 0, cells = 0, upload_bytes = 0), parent = emptyenv())
  summaries = query_fe_summaries(conn, centered, fe_columns, columns, sizes, budget, limits)
  info = summaries[["info"]]
  pairs = summaries[["pairs"]]
  sys = fe_system(info, pairs, limits)
  if (is.null(sys)) projection_abort("system_limit", "Fixed-effect matrix exceeds its projection budget.")
  if (any(!is.finite(sys[["diag"]])) || any(sys[["diag"]] <= 0)) projection_abort("singular", "Nonpositive fixed-effect group weights.")
  B = do.call(rbind, lapply(info, `[[`, "b"))
  cache = new.env(parent = emptyenv()); cache[["projection_limits"]] = limits
  if (two_family_pcg) cache[["pcg_max_iter"]] = max_iter
  solver = fe_solver(sys, cache)
  B = solver[["project"]](B)
  coefficients = solver[["solve"]](B)
  # The final row certificate is authoritative; refinement first removes
  # numerical group-equation error without changing statistical tolerances.
  target = min(1e-12, tol / 8) * max(1, abs(B / sqrt(sys[["diag"]])))
  for (cycle in seq_len(min(max_iter, 8L))) {
    residual = B - fe_multiply(sys, coefficients)
    error = max(abs(residual / sqrt(sys[["diag"]])))
    if (is.finite(error) && error <= target) break
    if (cycle == min(max_iter, 8L)) projection_abort("nonconverged", "Fixed-effect moment refinement did not converge.")
    coefficients = coefficients + solver[["solve"]](residual)
  }
  information = if (!materialize) fe_fwl_information(conn, centered, columns, B, coefficients, sys, budget, limits) else NULL
  lookup_cols = c("id", paste0("c", seq_along(columns)))
  lookups = vapply(seq_along(fe), function(i) {
    positions = seq.int(sys[["starts"]][i] + 1L, sys[["starts"]][i + 1L])
    fe_lookup_source(conn, coefficients[positions, , drop = FALSE], budget, limits, create, snapshot, lookup_cols)
  }, character(1))
  joins = paste0(" JOIN ", lookups, " p", seq_along(fe), " ON s.", fe_columns, "=p", seq_along(fe), ".id")
  subtract = vapply(seq_along(columns), function(j) paste0("(s.", columns[j], "-(", paste0("p", seq_along(fe), ".c", j, collapse = "+"), "))*(", sprintf("%.17e", scale[j]), ") AS ", columns[j]), character(1))
  projected = create(paste0("SELECT ", paste(c(paste0("s.", c(id_columns, "raw_w", "w", fe_columns)), subtract), collapse = ","), " FROM ", centered, " s", paste(joins, collapse = "")), c(id_columns, "raw_w", "w", fe_columns, columns))
  certificate = fe_projection_certificate(conn, projected, fe_columns, columns, nobs, original_scale)
  if (fe_certificate_shortfall(certificate, tol) > 1) {
    projection_abort("certificate_failed", "Fixed-effect moment projection failed its SQL row certificate.", c(certificate, list(target = tol)))
  }
  result = create(paste("SELECT", paste(c(paste(id_columns, "AS", q(ids)), "raw_w AS __w", paste(columns, "AS", q(paste0(vars, "_tilde")))) , collapse = ","), "FROM", projected), c(ids, "__w", paste0(vars, "_tilde")))
  moments = if (!materialize) fe_fwl_moments(conn, projected, columns, information, scale, original_scale, weight_scale, nobs, sizes, budget, limits) else NULL
  success = TRUE
  list(table = result, moments = moments, base_table = NULL, resources = resources,
    observation_tables = if (materialize) 5L else 1L,
    transformed_observation_tables = if (materialize) 4L else 0L,
    moment_method = if (materialize) "rows" else "fwl",
    fwl_recomputed = if (materialize) FALSE else attr(moments, "fwl_recomputed"),
    original_scale = setNames(original_scale, vars),
    iterations = cycle, accelerated_steps = 0L, max_abs = certificate[["max_abs"]], max_relative = certificate[["max_relative"]],
    method = if (two_family_pcg) "bounded_pcg" else "bounded_moments", transfer_cells = budget[["cells"]], summary_rows = budget[["rows"]], upload_bytes = budget[["upload_bytes"]],
    transport = summaries[["transport"]], summary_statements = summaries[["statements"]], wire_rows = summaries[["wire_rows"]], wire_bytes = summaries[["wire_bytes"]],
    estimated_bytes = projection_memory(sizes, length(columns), limits = limits), solver = solver[["method"]],
    solver_iterations = solver[["stats"]][["iterations"]], preconditioner = solver[["preconditioner"]])
}

# Bind independently named source roles once before either projection engine.
# The frame maps its SQL coordinates back; source keys and generated residuals
# never share a table target, even for names such as __w or x_tilde.
projection_frame_inputs = function(inputs, x_sql, x_names) {
  conn = inputs[["conn"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  ids = projection_id_columns(inputs)
  id_map = setNames(paste0("dbreg_p_id", seq_along(ids)), ids)
  x = paste0("dbreg_p_x", seq_along(x_names))
  y = "dbreg_p_y"
  w = if (is.null(inputs[["weights"]])) NULL else "dbreg_p_w"
  trends = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  trend_map = setNames(if (length(trends)) paste0("dbreg_p_t", seq_along(trends)) else character(), trends)
  fields = c(if (length(trends)) paste(q(trends), "AS", unname(trend_map)), paste(q(ids), "AS", id_map), paste(q(inputs[["yvar"]]), "AS", y),
    paste(x_sql, "AS", x), if (!is.null(w)) paste(q(inputs[["weights"]]), "AS", w))
  bound = list2env(as.list(inputs), parent = emptyenv())
  bound[["from_statement"]] = paste0("FROM (SELECT ", paste(fields, collapse = ","), " ", inputs[["from_statement"]], ") dbreg_projection_source")
  bound[["fe_slopes"]] = lapply(inputs[["fe_slopes"]], function(x) unname(trend_map[x]))
  bound[["fe"]] = unname(id_map[inputs[["fe"]]])
  bound[["cluster_var"]] = if (length(inputs[["cluster_var"]])) unname(id_map[inputs[["cluster_var"]]]) else NULL
  bound[["row_id_safe"]] = unname(id_map[inputs[["row_id_safe"]]])
  bound[["yvar"]] = y
  bound[["weights"]] = w
  list(inputs = bound, x = x, y = y, ids = id_map)
}

projection_frame_columns = function(conn, table, bound, variables) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  targets = safe_identifiers(c(paste0("dbreg_projection_v", seq_along(variables)), "dbreg_projection_weight"),
    taken = names(bound[["ids"]]), fold_case = folds_identifier_case(conn))
  projected = setNames(targets[seq_along(variables)], variables)
  weight = targets[length(targets)]
  fields = c(paste(q(unname(bound[["ids"]])), "AS", q(names(bound[["ids"]]))),
    paste(q(paste0(c(bound[["y"]], bound[["x"]]), "_tilde")), "AS", q(projected)), paste("__w AS", q(weight)))
  list(sql = paste("SELECT", paste(fields, collapse = ","), "FROM", table), columns = projected, weights = weight)
}


projection_id_columns = function(inputs) {
  unique(c(inputs[["fe"]], inputs[["cluster_var"]],
    if (!is.null(inputs[["retained_original"]])) inputs[["row_id_safe"]]))
}
