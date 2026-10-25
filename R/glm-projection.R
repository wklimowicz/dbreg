# SQL owns row state; bounded FE moments or FE-level CG determine projections.
solve_partial_wls = function(conn, table, x, fe, control, create, drop, response = TRUE, cache = NULL, materialize = TRUE, scales = NULL, extra_sums = character(), need_table = TRUE, precomputed = NULL) {
  columns = c(x, if (response) "z")
  backend_info = detect_backend(conn)
  slopes = attr(fe, "slopes")
  if (is.null(cache)) cache = new.env(parent = emptyenv())
  if (length(fe) > 1) {
    dbglm_moments_enabled(conn, table, fe, control, cache, length(columns))
    if (!cache[["moment_enabled"]] && length(fe) == 2L && !length(unlist(slopes)) &&
      identical(control[["projection_method"]], "auto") && length(columns)) {
      identifiers = unique(c(fe, unlist(slopes)))
      if (is.null(scales)) {
        norm = as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf(
          "FSUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)), collapse = ","), "FROM", table))[1, ])
        scales = sqrt(pmax(norm, .Machine[["double.xmin"]]))
      }
      target = control[["fe_tol"]]
      if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
      first = utils::tail(columns, 1)
      rest = utils::head(columns, -1)
      order = match(c(first, rest), columns)
      attempt = projection_attempt("bounded_pcg", dbreg_duckdb_cg_projections(conn, paste("FROM", table), fe, first, rest, rest, "w",
        cluster_var = "__dbglm_key", max_iter = control[["fe_max_iter"]], tol = target, scales = scales[order],
        materialize_result = FALSE, inline_result = !materialize, extra_sums = extra_sums))
      projected = record_projection_attempt(cache, attempt)
      if (!is.null(projected)) {
        on.exit(for (t in rev(projected[["all_tables"]])) {
          if (!materialize && isTRUE(backend_info[["is_databricks"]]) && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) next
          if (t %in% projected[["view_tables"]]) try(DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", t)), silent = TRUE) else
            drop_ap_table_if_exists(conn, t, backend_info)
        }, add = TRUE)
        if (!materialize && isTRUE(backend_info[["is_databricks"]]) && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
          cache[["views"]] = unique(c(cache[["views"]], projected[["view_tables"]]))
          cache[["tables"]] = unique(c(cache[["tables"]], setdiff(projected[["all_tables"]], projected[["view_tables"]])))
        }
        result = create(paste("SELECT __dbglm_key,", paste(c(identifiers, "__w AS w",
          sprintf("%s_tilde AS %s", columns, columns)), collapse = ","), "FROM", projected[["table"]]), materialize = materialize)
        inv_order = match(columns, c(first, rest))
        reduction = if (!is.null(projected[["cross_reduction"]])) projected[["cross_reduction"]][inv_order, inv_order, drop = FALSE] else NULL
        extra_red = if (!is.null(projected[["extra_reduction"]])) projected[["extra_reduction"]][inv_order, , drop = FALSE] else NULL
        raw_proj = if (!materialize && isTRUE(backend_info[["is_databricks"]]) &&
          !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) projected[["table"]] else NULL
        return(list(table = result, raw_projected = raw_proj, cycles = projected[["iterations"]], residual = projected[["max_abs"]], method = "bounded_pcg",
          cross_reduction = reduction, extra_reduction = extra_red, sub_tables = if (!materialize) projected[["all_tables"]] else character(),
          s1 = projected[["s1"]][, inv_order, drop = FALSE], s2 = projected[["s2"]][, inv_order, drop = FALSE],
          a = projected[["a"]][, inv_order, drop = FALSE], b = projected[["b"]][, inv_order, drop = FALSE]))
      }
    }
    if (cache[["moment_enabled"]]) {
      ops = lapply(fe, function(f) list(f = f, v = "1.0"))
      for (i in seq_along(fe)) for (v in slopes[[i]]) ops[[length(ops) + 1]] = list(f = fe[i], v = v)
      identifiers = unique(c(fe, unlist(slopes)))
      if (is.null(scales)) {
        norm = if (length(columns)) as.numeric(backend_get_query(conn, paste("SELECT", paste(sprintf("FSUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)), collapse = ","), "FROM", table))[1, ]) else numeric()
        scales = sqrt(pmax(norm, .Machine[["double.xmin"]]))
      }
      target = control[["fe_tol"]]
      if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
      attempt = projection_attempt("bounded_moments", project_fe_moments(conn, table, columns, identifiers, ops, scales, target, control, create, drop, cache,
        materialize = materialize, extra_sums = extra_sums, need_table = need_table, precomputed = precomputed))
      ans = record_projection_attempt(cache, attempt)
      if (!is.null(ans)) return(ans)
      cache[["moment_enabled"]] = FALSE
    }
    if (isTRUE(backend_info[["is_databricks"]]) && !length(unlist(slopes))) {
      return(dbglm_project_window(conn, table, columns, fe, control, create, drop, response))
    }
  }
  projection = create(paste("SELECT __dbglm_key,", paste(unique(c(fe, unlist(slopes), "w", columns)), collapse = ","), "FROM", table))
  ops = lapply(fe, function(f) list(f = f, v = "1.0"))
  basis = character()
  for (i in seq_along(fe)) {
    local_basis = character()
    for (j in seq_along(slopes[[i]])) {
      v = slopes[[i]][j]
      name = paste0("__dbglm_v", i, "_", j)
      partition = paste0(" OVER (PARTITION BY ", fe[i], ")")
      centered = paste0(v, " - SUM(w*", v, ")", partition, "/SUM(w)", partition)
      next_projection = create(paste("SELECT *,", centered, "AS", name, "FROM", projection))
      drop(projection)
      projection = next_projection
      for (prev in local_basis) {
        adjustment = paste0(name, " - ", prev, "*COALESCE(SUM(w*", prev, "*", name, ")", partition,
          "/NULLIF(SUM(w*", prev, "*", prev, ")", partition, ",0),0)")
        next_projection = create(paste("SELECT * EXCLUDE(", name, "),", adjustment, "AS", name, "FROM", projection))
        drop(projection)
        projection = next_projection
      }
      local_basis = c(local_basis, name)
      basis = c(basis, name)
      ops[[length(ops) + 1]] = list(f = fe[i], v = name)
    }
  }
  identifiers = unique(c(fe, unlist(slopes), basis))
  if (!length(columns)) return(list(table = projection, cycles = 0L, residual = 0, method = "group"))
  norm = backend_get_query(conn, paste("SELECT", paste(sprintf("SUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)), collapse = ","), "FROM", projection))
  scales = sqrt(pmax(as.numeric(norm[1, ]), .Machine[["double.xmin"]]))
  target = control[["fe_tol"]]
  if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
  if (length(fe) > 1) return(dbglm_project_cg(conn, projection, columns, identifiers, ops, scales, target, control, create, drop))
  residual = 0
  for (cycle in seq_len(control[["fe_max_iter"]])) {
    for (op in ops) {
      f = op[["f"]]
      v = op[["v"]]
      values = sprintf("(%1$s - (%3$s)*COALESCE(SUM(w*(%3$s)*%1$s) OVER (PARTITION BY %2$s)/NULLIF(SUM(w*(%3$s)*(%3$s)) OVER (PARTITION BY %2$s),0),0)) AS %1$s", columns, f, v)
      next_projection = create(paste("SELECT __dbglm_key,", paste(c(identifiers, "w", values), collapse = ","), "FROM", projection))
      drop(projection)
      projection = next_projection
    }
    residual = 0
    for (op in ops) {
      f = op[["f"]]
      v = op[["v"]]
      terms = sprintf("COALESCE(ABS(SUM(w*(%3$s)*%1$s))/NULLIF(SQRT(SUM(w*(%3$s)*(%3$s))),0),0)/%2$s AS r%4$d", columns, dbglm_number(scales), v, seq_along(columns))
      row = backend_get_query(conn, paste0("SELECT ", paste(sprintf("MAX(r%d) AS r%d", seq_along(columns), seq_along(columns)), collapse = ","),
        " FROM (SELECT ", paste(terms, collapse = ","), " FROM ", projection, " GROUP BY ", f, ") t"))
      residual = max(residual, as.numeric(row[1, ]))
    }
    if (residual <= target) return(list(table = projection, cycles = cycle, residual = residual, method = "group"))
  }
  stop("GLM fixed-effect projection did not converge.", call. = FALSE)
}

# Direct FE moments need the non-absorbed levels to fit a dense solve in R;
# see fe_dense_feasible(). The decision is cached for the whole fit.
dbglm_moments_enabled = function(conn, table, fe, control, cache, columns = 1L) {
  if (!is.null(cache[["moment_enabled"]])) return(cache[["moment_enabled"]])
  if (identical(control[["projection_method"]], "sql")) {
    cache[["moment_enabled"]] = FALSE
    return(FALSE)
  }
  if (is.null(cache[["projection_limits"]])) cache[["projection_limits"]] = projection_limits()
  sizes = cache[["fe_sizes"]]
  if (is.null(sizes)) {
    count = backend_get_query(conn, paste("SELECT", paste(sprintf("COUNT(DISTINCT %s) AS n%d", fe, seq_along(fe)), collapse = ","), "FROM", table))
    sizes = vapply(count, as.numeric, numeric(1))
  }
  cache[["fe_sizes"]] = sizes
  op_sizes = rep(sizes, 1 + lengths(attr(fe, "slopes")))
  cache[["moment_enabled"]] = sum(op_sizes)*max(1L,columns) <= cache[["projection_limits"]][["design_entries"]] &&
    fe_dense_feasible(op_sizes, columns, limits = cache[["projection_limits"]])
  if (!cache[["moment_enabled"]]) record_projection_attempt(cache,
    list(method = "bounded_moments", status = "declined", reason = "system_limit",
      details = list(fe_sizes = op_sizes, estimated_bytes = projection_memory(op_sizes, columns), limits = cache[["projection_limits"]])))
  cache[["moment_enabled"]]
}
