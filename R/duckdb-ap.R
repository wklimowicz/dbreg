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
  equation_metric = !is.null(scales)
  if (equation_metric && (length(scales) != length(vars) ||
      any(!is.finite(scales)) || any(scales <= 0))) {
    stop("Projection scales must be positive and match the projected columns.", call. = FALSE)
  }
  aggregate = if (equation_metric) "FSUM" else "SUM"
  vars_sql = c(yvar, xvars_sql)
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
  created = character()
  drop_state = function(table) {
    if(table %in% views) try(DBI::dbExecute(conn,paste("DROP VIEW IF EXISTS",table)),silent=TRUE) else
      drop_ap_table_if_exists(conn,table,backend_info)
  }
  success = FALSE
  on.exit({
    if (!success) {
      for (table in rev(created)) drop_state(table)
    }
  }, add = TRUE)
  create = function(name, sql) {
    created <<- unique(c(created, tables[[name]]))
    if(isTRUE(backend_info$is_databricks)) {
      DBI::dbExecute(conn,dbglm_sql(conn,paste(if(tables[[name]] %in% views) "CREATE TEMP VIEW" else "CREATE TEMP TABLE",tables[[name]],"AS",sql)))
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
  edges = dbglm_get_query(conn, edges_sql)
  if (!nrow(edges) || nrow(edges) > max_edges) {
    if (verbose) message("[dbreg] FE summaries exceed DuckDB transfer limits; using join AP")
    return(NULL)
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

  # S = D2 - C' D1^-1 C. Applying S never constructs its dense matrix.
  rhs = s2 - as.matrix(Matrix::crossprod(C, y1))
  apply_s = function(b) {
    b * d2 - as.matrix(Matrix::crossprod(C, as.matrix(C %*% b) / d1))
  }
  normalize = function(r) {
    if (equation_metric) sweep(r / sqrt(d2), 2L, scales, "/") else r / d2
  }
  diagonal = d2 - as.numeric(Matrix::colSums(C * (Matrix::Diagonal(x = 1 / d1) %*% C)))
  diagonal = pmax(diagonal, .Machine$double.eps * pmax(d2, 1))
  b = matrix(0, nrow(y2), ncol(y2))
  r = rhs
  z = r / diagonal
  direction = z
  rz = colSums(r * z)
  residual = max(abs(normalize(r)))
  iterations = 0L
  for (iter in seq_len(max_iter)) {
    if (is.finite(residual) && residual <= tol) break
    active = apply(abs(normalize(r)), 2L, max) > tol
    product = apply_s(direction)
    denominator = colSums(direction * product)
    if (any(!is.finite(denominator[active])) || any(denominator[active] <= 0)) return(NULL)
    alpha = numeric(ncol(y2))
    alpha[active] = rz[active] / denominator[active]
    if (any(!is.finite(alpha))) return(NULL)
    b = b + sweep(direction, 2L, alpha, "*")
    # Recompute the true residual rather than trusting the CG recurrence.
    r = rhs - apply_s(b)
    residual = max(abs(normalize(r)))
    iterations = iter
    if (verbose) message("[AP] iter ", iter, ": max abs mean = ", sprintf("%.4e", residual))
    if (is.finite(residual) && residual <= tol) break
    z = r / diagonal
    new_rz = colSums(r * z)
    beta = numeric(ncol(y2))
    beta[active] = new_rz[active] / rz[active]
    if (any(!is.finite(beta))) return(NULL)
    direction = z + sweep(direction, 2L, beta, "*")
    rz = new_rz
  }
  if (!is.finite(residual) || residual > tol) return(NULL)
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

  if (inline_result && isTRUE(backend_info$is_databricks) &&
      !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
    lookup_sqls = character(2)
    for (k in 1:2) {
      adjustment = means[[k]][keys[k]]
      effect = if (k == 1L) a else b
      for (j in seq_along(vars)) adjustment[[paste0("a", j)]] = effect[, j]
      lsql = dbglm_values_sql(conn, adjustment)
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
      iterations = iterations, max_abs = actual, method = "bounded_pcg", backend_info = backend_info,
      cross_reduction = cross_reduction, extra_reduction = extra_reduction, s1 = s1, s2 = s2, a = a, b = b))
  }
  for (k in 1:2) {
    adjustment = means[[k]][keys[k]]
    effect = if (k == 1L) a else b
    for (j in seq_along(vars)) adjustment[[paste0("a", j)]] = effect[, j]
    if(isTRUE(backend_info$is_databricks)) {
      lookup_sql = dbglm_values_sql(conn,adjustment)
      if(is.null(lookup_sql)) return(NULL)
      create(if(k==1L) "a" else "b",lookup_sql)
    } else dbWriteTable(conn, tables[[if (k == 1L) "a" else "b"]], adjustment,
      temporary = TRUE, overwrite = TRUE)
    created = unique(c(created, tables[[if (k == 1L) "a" else "b"]]))
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
    actual = dbglm_get_query(conn, paste0("SELECT MAX(max_abs) AS max_abs FROM (",
      paste(checks, collapse = " UNION ALL "), ") q"))[["max_abs"]]
    n_result = dbglm_get_query(conn, paste("SELECT COUNT(*) AS n FROM", tables[["result"]]))[["n"]]
    n_source = dbglm_get_query(conn, paste("SELECT COUNT(*) AS n", from_statement))[["n"]]
    if (!is.finite(actual) || actual > tol || n_result != n_source) return(NULL)
  }
  success = TRUE
  list(table = tables[["result"]], base_table = NULL, all_tables = created, view_tables = intersect(created, views),
    iterations = iterations, max_abs = actual, method = "bounded_pcg", backend_info = backend_info,
    cross_reduction = cross_reduction, extra_reduction = extra_reduction, s1 = s1, s2 = s2, a = a, b = b)
}
