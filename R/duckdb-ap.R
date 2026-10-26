# Work in D2-scaled coordinates so component deflation remains stable for
# unequal group weights. Bounded Schur fill permits symmetric Gauss-Seidel;
# otherwise use an additive block SGS/Jacobi preconditioner without that fill.
two_way_cg_preconditioner = function(C, d1, d2, limits, available_bytes) {
  n1 = nrow(C)
  n2 = ncol(C)
  system = list(n = n1 + n2, big = seq_len(n1), rest = n1 + seq_len(n2),
    Cm = C, Dg = d1, diag = c(d1, d2), sizes = c(n1, n2))
  components = fe_two_components(system)
  if (is.null(components)) return(NULL)
  components = components[n1 + seq_len(n2)]
  total = as.numeric(rowsum(d2, components, reorder = TRUE))
  root = sqrt(d2)
  null = root / sqrt(total[components])
  project = function(value) {
    sums = rowsum(null * value, components, reorder = TRUE)
    value - null * sums[components, , drop = FALSE]
  }
  row_scaled = C
  row_scaled@x = row_scaled@x / d1[row_scaled@i + 1L]
  diagonal = pmax(1 - as.numeric(Matrix::colSums(C * row_scaled)) / d2, .Machine[["double.eps"]])
  walk = function(value) {
    as.matrix(Matrix::crossprod(C, as.matrix(C %*% (value / root)) / d1)) / root
  }
  degree = as.numeric(Matrix::rowSums(C != 0))
  fill = n2 + sum(degree^2)
  if (fill <= limits[["cross_entries"]] && 64 * fill <= available_bytes) {
    scaled = Matrix::Diagonal(x = 1 / sqrt(d1)) %*% C %*% Matrix::Diagonal(x = 1 / root)
    S = Matrix::Diagonal(n2) - Matrix::crossprod(scaled)
    S = S + Matrix::Diagonal(x = diagonal - as.numeric(Matrix::diag(S)))
    lower = Matrix::tril(S)
    upper = Matrix::t(lower)
    precondition = function(value) as.matrix(Matrix::solve(upper,
      diagonal * as.matrix(Matrix::solve(lower, value))))
    method = "schur_sgs"
  } else {
    precondition = function(value) value / diagonal + value + walk(value)
    method = "block_sgs"
  }
  list(project = project, precondition = precondition, walk = walk,
    root = root, method = method, components = max(components))
}

#' Bounded two-way fixed-effect projection for DuckDB and Databricks
#'
#' The observation relation remains in the database. Weighted FE-pair counts and
#' FE-level sums define a sparse Schur complement over the smaller FE family.
#' Only these bounded sufficient statistics are transferred to R.
#'
#' Returns NULL when the bounds or numerical checks require the join/window fallback.
#' @keywords internal
#' @noRd
dbreg_duckdb_cg_projections = function(
  conn, from_statement, fe, yvar, xvars_sql, xvar_names, weights,
  cluster_var = NULL, verbose = FALSE, max_iter = 100L, tol = 1e-10,
  max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6),
  max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6), scales = NULL,
  materialize_result = TRUE, inline_result = FALSE, extra_sums = character()
) {
  valid_cap = function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  }
  if (!valid_cap(max_edges) || !valid_cap(max_entries)) {
    stop("DuckDB AP transfer limits must be positive whole numbers.", call. = FALSE)
  }
  stopifnot(length(fe) == 2L)
  backend_info = detect_backend(conn)
  vars = c(yvar, xvar_names)
  limits = projection_limits()
  width = length(vars) + length(extra_sums)
  key_bytes = function(f) {
    value = paste0("CAST(", f, " AS VARCHAR)")
    if (backend_info[["name"]] == "duckdb") paste0("OCTET_LENGTH(encode(", value, "))") else
      paste0("OCTET_LENGTH(", value, ")")
  }
  w = if (is.null(weights)) "1.0" else weights
  raw_vars = sql_cast_double(c(yvar, xvars_sql), conn)
  w_double = sql_cast_double(w, conn)
  center_terms = paste0(vapply(paste0("(", w_double, ")*(", raw_vars, ")"),
    function(value) sql_stable_sum(conn, value), character(1)), "/(", sql_stable_sum(conn, w_double),
    ") AS center", seq_along(raw_vars))
  preflight = backend_get_query(conn, paste("SELECT COUNT(*) AS source_n,",
    paste(c(sprintf("COUNT(DISTINCT %s) AS n%d", fe, 1:2),
      sprintf("MAX(%s) AS bytes%d", vapply(fe, key_bytes, character(1)), 1:2), center_terms), collapse = ","), from_statement))
  sizes = as.numeric(preflight[c("n1", "n2")][1L, ])
  key_width = as.numeric(preflight[c("bytes1", "bytes2")][1L, ])
  if (any(!is.finite(c(sizes, key_width))) || any(sizes < 1L)) return(NULL)
  centers = vapply(seq_along(vars), function(i) as.numeric(preflight[[paste0("center", i)]]), numeric(1))
  if (any(!is.finite(centers))) return(NULL)
  vars_sql = paste0("(", raw_vars, "-(", dbglm_number(centers), "))")
  if (sum(sizes)*length(vars) > max_entries) projection_abort("system_limit", "Sparse FE system exceeds its design budget.")
  row_bytes = 4*(sum(key_width) + 128 + 8*(width + 3L))
  remaining = limits[["memory_bytes"]] - 8*sum(sizes)*max(1L, width)*12
  max_edges = min(max_edges, floor(limits[["transfer_cells"]]/(width + 3L))-1L,
    floor(remaining/row_bytes)-1L)
  if (max_edges < 1L) projection_abort("system_limit", "Sparse FE system exceeds its memory or transfer budget.")
  equation_metric = !is.null(scales)
  if (equation_metric && (length(scales) != length(vars) ||
      any(!is.finite(scales)) || any(scales <= 0))) {
    stop("Projection scales must be positive and match the projected columns.", call. = FALSE)
  }
  aggregate = if (equation_metric) "FSUM" else "SUM"
  ids = setdiff(unique(c(fe, cluster_var)), vars)
  if ("__w" %in% c(ids, vars)) {
    stop("Model column conflicts with the internal AP column __w.", call. = FALSE)
  }
  seed = paste0(Sys.getpid(), "_", sprintf("%08d", sample.int(1e8, 1L)))
  tables = setNames(
    paste0("dbreg_ap_cg_", seed, "_", c("source", "edges", "a", "b", "result")),
    c("source", "edges", "a", "b", "result")
  )
  views = if(isTRUE(backend_info$is_databricks) &&
    !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates",FALSE))) {
    unname(tables[c("source","edges","a","b", if (!materialize_result) "result")])
  } else character()
  state = new.env(parent = emptyenv())
  state[["created"]] = character()
  drop_state = function(table) {
    if(table %in% views) try(DBI::dbExecute(conn,paste("DROP VIEW IF EXISTS",table)),silent=TRUE) else
      drop_ap_table_if_exists(conn,table,backend_info)
  }
  success = FALSE
  on.exit({
    if (!success) {
      for (table in rev(state[["created"]])) drop_state(table)
    }
  }, add = TRUE)
  create = function(name, sql) {
    state[["created"]] = unique(c(state[["created"]], tables[[name]]))
    if(isTRUE(backend_info$is_databricks)) {
      DBI::dbExecute(conn,backend_sql(conn,paste(if(tables[[name]] %in% views) "CREATE TEMP VIEW" else "CREATE TEMP TABLE",tables[[name]],"AS",sql)))
    } else create_temp_table_as(conn, tables[[name]], sql, "duckdb")
  }
  w = if (is.null(weights)) "1.0" else weights
  # String keys preserve BIGINT identifiers beyond R's exact double range.
  keys = c("dbreg_cg_key_1", "dbreg_cg_key_2")
  key_sql = sprintf("CAST(%s AS VARCHAR)", fe)
  sum_terms = c(
    sprintf("%s((%s) * (%s)) AS s%d", aggregate, w, vars_sql, seq_along(vars)),
    if (length(extra_sums)) sprintf("%s(%s) AS ex%d", aggregate, extra_sums, seq_along(extra_sums))
  )
  edges_sql = paste0(
    "SELECT ", paste(sprintf("%s AS %s", key_sql, keys), collapse = ", "),
    ", ", aggregate, "(", w, ") AS ew, ",
    paste(sum_terms, collapse = ", "), " ", from_statement,
    " GROUP BY ", paste(fe, collapse = ", "),
    " LIMIT ", sprintf("%.0f", max_edges + 1)
  )
  edges = backend_get_query(conn, edges_sql)
  if (!nrow(edges) || nrow(edges) > max_edges) {
    if (verbose) message("[dbreg] FE summaries exceed DuckDB transfer limits; using join AP")
    projection_abort("transfer_limit", "FE pair summaries exceed sparse projection limits.", list(summary_rows = max_edges))
  }
  u1 = unique(edges[[keys[1]]])
  u2 = unique(edges[[keys[2]]])
  entries = (as.numeric(length(u1)) + as.numeric(length(u2))) * length(vars)
  if (!is.finite(entries) || entries > max_entries) {
    if (verbose) message("[dbreg] FE summaries exceed DuckDB transfer limits; using join AP")
    return(NULL)
  }
  scols = paste0("s", seq_along(vars))
  s_mat = as.matrix(edges[scols])
  ew = as.numeric(edges[["ew"]])
  if (any(!is.finite(ew)) || any(ew <= 0) || any(!is.finite(s_mat))) return(NULL)
  idx1 = match(edges[[keys[1]]], u1)
  idx2 = match(edges[[keys[2]]], u2)
  d1 = as.numeric(rowsum(ew, idx1, reorder = FALSE))
  d2 = as.numeric(rowsum(ew, idx2, reorder = FALSE))
  s1 = rowsum(s_mat, idx1, reorder = FALSE)
  s2 = rowsum(s_mat, idx2, reorder = FALSE)
  ex1 = ex2 = NULL
  if (length(extra_sums)) {
    ex_mat = as.matrix(edges[paste0("ex", seq_along(extra_sums))])
    ex1 = rowsum(ex_mat, idx1, reorder = FALSE)
    ex2 = rowsum(ex_mat, idx2, reorder = FALSE)
  }
  if (any(!is.finite(c(d1, d2))) || any(c(d1, d2) <= 0)) return(NULL)
  y1 = s1 / d1
  y2 = s2 / d2
  means = list(
    local({ df = data.frame(k = u1, gw = d1, stringsAsFactors = FALSE); names(df)[1] = keys[1]; df }),
    local({ df = data.frame(k = u2, gw = d2, stringsAsFactors = FALSE); names(df)[1] = keys[2]; df })
  )
  swapped = FALSE
  if (nrow(means[[1]]) < nrow(means[[2]])) {
    fe = rev(fe)
    keys = rev(keys)
    key_sql = rev(key_sql)
    means = rev(means)
    tmp_idx = idx1; idx1 = idx2; idx2 = tmp_idx
    tmp_d = d1; d1 = d2; d2 = tmp_d
    tmp_s = s1; s1 = s2; s2 = tmp_s
    tmp_y = y1; y1 = y2; y2 = tmp_y
    if (length(extra_sums)) { tmp_ex = ex1; ex1 = ex2; ex2 = tmp_ex }
    swapped = TRUE
  }
  C = Matrix::sparseMatrix(
    i = idx1,
    j = idx2, x = ew,
    dims = c(nrow(means[[1]]), nrow(means[[2]]))
  )

  # S = D2 - C' D1^-1 C. The iterative operator never constructs a
  # dense Schur matrix; every vector stays orthogonal to component nulls.
  preconditioner = two_way_cg_preconditioner(C, d1, d2, limits,
    max(0, remaining - nrow(edges) * row_bytes))
  if (is.null(preconditioner)) return(NULL)
  project = preconditioner[["project"]]
  precondition = preconditioner[["precondition"]]
  root = preconditioner[["root"]]
  rhs = s2 - as.matrix(Matrix::crossprod(C, y1))
  rhs_scaled = project(rhs / root)
  apply_s = function(value) value - preconditioner[["walk"]](value)
  normalize = function(value) {
    if (equation_metric) sweep(value, 2L, scales, "/") else value / root
  }
  scaled_b = matrix(0, nrow(y2), ncol(y2))
  r = rhs_scaled
  z = project(precondition(r))
  direction = z
  rz = colSums(r * z)
  residual = max(abs(normalize(r)))
  iterations = 0L
  for (iter in seq_len(max_iter)) {
    if (is.finite(residual) && residual <= tol) break
    active = apply(abs(normalize(r)), 2L, max) > tol
    product = project(apply_s(direction))
    denominator = colSums(direction * product)
    if (any(!is.finite(denominator[active])) || any(denominator[active] <= 0)) return(NULL)
    alpha = numeric(ncol(y2))
    alpha[active] = rz[active] / denominator[active]
    if (any(!is.finite(alpha))) return(NULL)
    scaled_b = project(scaled_b + sweep(direction, 2L, alpha, "*"))
    # Recompute the true residual rather than trusting the CG recurrence.
    r = project(rhs_scaled - apply_s(scaled_b))
    residual = max(abs(normalize(r)))
    iterations = iter
    if (verbose) message("[AP] iter ", iter, ": max abs mean = ", sprintf("%.4e", residual))
    if (is.finite(residual) && residual <= tol) break
    z = project(precondition(r))
    new_rz = colSums(r * z)
    beta = numeric(ncol(y2))
    beta[active] = new_rz[active] / rz[active]
    if (any(!is.finite(beta))) return(NULL)
    direction = project(z + sweep(direction, 2L, beta, "*"))
    rz = new_rz
  }
  if (!is.finite(residual) || residual > tol) return(NULL)
  b = scaled_b / root
  a = y1 - as.matrix(C %*% b) / d1
  # Exact cell-level residual verification in R: r1 = s1 - d1*a - C*b, r2 = s2 - C'*a - d2*b
  r1 = s1 - a * d1 - as.matrix(C %*% b)
  r2 = s2 - as.matrix(Matrix::crossprod(C, a)) - b * d2
  norm1 = if (equation_metric) sweep(abs(r1) / sqrt(d1), 2L, scales, "/") else abs(r1) / d1
  norm2 = if (equation_metric) sweep(abs(r2) / sqrt(d2), 2L, scales, "/") else abs(r2) / d2
  actual = max(norm1, norm2)
  if (!is.finite(actual) || actual > tol) return(NULL)
  # Cross-product reduction from FE projection: V' W V - V_tilde' W V_tilde
  # Since V_tilde = V - D1*a - D2*b, V' W V - V_tilde' W V_tilde = 2*(S1'a + S2'b) - (a' D1' D1 a + b' D2' D2 b + 2 a' C b)
  cross_reduction = crossprod(s1, a) + crossprod(s2, b) + crossprod(r1, a) + crossprod(r2, b)
  cross_reduction = 0.5 * (cross_reduction + t(cross_reduction))
  extra_reduction = if (length(extra_sums)) crossprod(a, ex1) + crossprod(b, ex2) else NULL
  # Public reductions refer to the original columns, while the streamed
  # within rows use centered coordinates. Restore the constant contribution.
  total_centered = colSums(s_mat)
  cross_reduction = cross_reduction + outer(centers, total_centered) +
    outer(total_centered, centers) + sum(ew) * tcrossprod(centers)
  if (length(extra_sums)) extra_reduction = extra_reduction + outer(centers, colSums(ex1))
  raw_a = sweep(a, 2L, centers, "+")
  raw_s1 = s1 + outer(d1, centers)
  raw_s2 = s2 + outer(d2, centers)

  if (inline_result && isTRUE(backend_info$is_databricks) &&
      !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
    lookup_sqls = character(2)
    for (k in 1:2) {
      adjustment = means[[k]][keys[k]]
      effect = if (k == 1L) a else b
      for (j in seq_along(vars)) adjustment[[paste0("a", j)]] = effect[, j]
      lsql = sql_values_frame(conn, adjustment)
      if (is.null(lsql)) return(NULL)
      lookup_sqls[k] = lsql
    }
    w_select = if (is.null(weights)) "1.0 AS __w" else paste0("s.", weights, " AS __w")
    result_sql = paste0(
      "(SELECT /*+ BROADCAST(a, b) */ s.*, ", paste(c(w_select,
        sprintf("(%s) - a.a%d - b.a%d AS %s_tilde", vars_sql, seq_along(vars), seq_along(vars), vars)), collapse = ", "),
      " FROM (SELECT * ", from_statement, ") s JOIN (", lookup_sqls[1],
      ") a ON CAST(s.", fe[1], " AS VARCHAR) = a.", keys[1],
      " JOIN (", lookup_sqls[2], ") b ON CAST(s.", fe[2], " AS VARCHAR) = b.", keys[2], ")"
    )
    success = TRUE
    return(list(table = result_sql, base_table = NULL, all_tables = character(), view_tables = character(),
      iterations = iterations, max_abs = actual, method = "bounded_pcg", backend_info = backend_info, preconditioner = preconditioner[["method"]],
      components = preconditioner[["components"]],
      cross_reduction = cross_reduction, extra_reduction = extra_reduction, s1 = raw_s1, s2 = raw_s2, a = raw_a, b = b,
      transfer_cells = nrow(edges)*(width + 3L), estimated_bytes = nrow(edges)*row_bytes + 8*sum(sizes)*max(1L, width)*12))
  }
  for (k in 1:2) {
    adjustment = means[[k]][keys[k]]
    effect = if (k == 1L) a else b
    for (j in seq_along(vars)) adjustment[[paste0("a", j)]] = effect[, j]
    if(isTRUE(backend_info$is_databricks)) {
      lookup_sql = sql_values_frame(conn,adjustment)
      if(is.null(lookup_sql)) return(NULL)
      create(if(k==1L) "a" else "b",lookup_sql)
    } else dbWriteTable(conn, tables[[if (k == 1L) "a" else "b"]], adjustment,
      temporary = TRUE, overwrite = TRUE)
    state[["created"]] = unique(c(state[["created"]], tables[[if (k == 1L) "a" else "b"]]))
  }
  broadcast_hint = if (isTRUE(backend_info$is_databricks)) "/*+ BROADCAST(a, b) */ " else ""
  w_select = if (is.null(weights)) "1.0 AS __w" else paste0("s.", weights, " AS __w")
  create("result", paste0(
    "SELECT ", broadcast_hint, paste(c(paste0("s.", ids), w_select,
      sprintf("(%s) - a.a%d - b.a%d AS %s_tilde", vars_sql, seq_along(vars), seq_along(vars), vars)), collapse = ", "),
    " FROM (SELECT * ", from_statement, ") s JOIN ", tables[["a"]],
    " a ON CAST(s.", fe[1], " AS VARCHAR) = a.", keys[1],
    " JOIN ", tables[["b"]], " b ON CAST(s.", fe[2], " AS VARCHAR) = b.", keys[2]
  ))
  if (!isTRUE(backend_info$is_databricks) || isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
    checks = vapply(fe, function(f) {
      terms = if (equation_metric) sprintf(
        "FSUM(__w * %s_tilde) / SQRT(FSUM(__w)) / %.17g AS m%d", vars, scales, seq_along(vars)
      ) else sprintf("SUM(__w * %s_tilde) / SUM(__w) AS m%d", vars, seq_along(vars))
      paste0("SELECT MAX(GREATEST(", paste(sprintf("ABS(m%d)", seq_along(vars)), collapse = ", "),
        ")) AS max_abs FROM (SELECT ",
        paste(terms, collapse = ", "),
        " FROM ", tables[["result"]], " GROUP BY ", f, ") q")
    }, character(1))
    actual = backend_get_query(conn, paste0("SELECT MAX(max_abs) AS max_abs FROM (",
      paste(checks, collapse = " UNION ALL "), ") q"))[["max_abs"]]
    n_result = backend_get_query(conn, paste("SELECT COUNT(*) AS n FROM", tables[["result"]]))[["n"]]
    n_source = as.numeric(preflight[["source_n"]])
    if (!is.finite(actual) || actual > tol || n_result != n_source) return(NULL)
  }
  success = TRUE
  list(table = tables[["result"]], base_table = NULL, all_tables = state[["created"]], view_tables = intersect(state[["created"]], views),
    iterations = iterations, max_abs = actual, method = "bounded_pcg", backend_info = backend_info, preconditioner = preconditioner[["method"]],
      components = preconditioner[["components"]],
    cross_reduction = cross_reduction, extra_reduction = extra_reduction, s1 = raw_s1, s2 = raw_s2, a = raw_a, b = b,
    transfer_cells = nrow(edges)*(width + 3L), estimated_bytes = nrow(edges)*row_bytes + 8*sum(sizes)*max(1L, width)*12)
}
