# Coefficient-space CG: SQL owns all FE vectors and observation transforms.
# R receives one bounded row of scalar diagnostics per query.
fe_cg_fetch = function(conn, sql, width, budget, limits) {
  if (budget[["rows"]] + 1 > limits[["summary_rows"]] || budget[["cells"]] + width > limits[["transfer_cells"]]) {
    projection_abort("transfer_limit", "SQL CG diagnostics exceed the transfer budget.", list(rows = budget[["rows"]] + 1, cells = budget[["cells"]] + width, limits = limits))
  }
  value = DBI::dbGetQuery(conn, sql)
  if (nrow(value) != 1L || ncol(value) != width) stop("Incomplete SQL CG diagnostics.", call. = FALSE)
  budget[["rows"]] = budget[["rows"]] + 1
  budget[["cells"]] = budget[["cells"]] + width
  value
}

fe_cg_manager = function(conn, resources) {
  factory = fe_projection_tables(conn, resources)
  resources[["peak_tables"]] = length(resources[["tables"]])
  list(create = function(sql) {
    table = factory(sql, character())
    resources[["peak_tables"]] = max(resources[["peak_tables"]], length(resources[["tables"]]))
    table
  }, drop = function(table) {
    active = as.character(DBI::dbQuoteIdentifier(conn, resources[["tables"]]))
    at = match(table, active)
    if (is.na(at)) stop("SQL CG tried to drop a borrowed table.", call. = FALSE)
    DBI::dbExecute(conn, paste(if (isTRUE(detect_backend(conn)[["is_databricks"]])) "DROP TEMP TABLE" else "DROP TABLE", table))
    resources[["tables"]] = resources[["tables"]][-at]
    invisible(NULL)
  })
}

fe_cg_control = function() {
  values = list(fe_max_iter = getOption("dbreg.cg_max_iter", 10000L), restart = getOption("dbreg.cg_restart", 500L))
  for (value in values) if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value > .Machine[["integer.max"]] || value != floor(value)) {
    stop("SQL CG iteration/restart controls must be positive whole numbers.", call. = FALSE)
  }
  values
}

fe_sql_cg = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop, budget, limits, allow_zero = FALSE) {
  backend_info = detect_backend(conn)
  sum_sql = function(expr) sql_sum(conn, expr, backend_info)
  number = function(x) sprintf("%.17e", x)
  fetch = function(sql, width) fe_cg_fetch(conn, sql, width, budget, limits)
  index = seq_along(columns)
  r = paste0("__cg_r", index)
  d = paste0("__cg_d", index)
  g = paste0("__cg_g", index)
  a = paste0("__cg_a", index)
  direction_names = paste0("__cg_value", index)
  key = "__cg_op,__cg_group"
  moments = function(source, inputs, outputs, info = NULL) {
    pieces = vapply(seq_along(ops), function(i) {
      op = ops[[i]]
      terms = c(
        if (is.null(info)) paste0("SQRT(", sum_sql(sprintf("w*(%1$s)*(%1$s)", op[["v"]])), ") AS __cg_scale"),
        paste(sum_sql(sprintf("w*(%s)*(%s)", op[["v"]], inputs)), "AS", paste0("u", index))
      )
      paste0(
        "SELECT ", i, " AS __cg_op,", op[["f"]], " AS __cg_group,", paste(terms, collapse = ","),
        " FROM ", source, " GROUP BY ", op[["f"]]
      )
    }, character(1))
    denominator = if (is.null(info)) "u.__cg_scale" else "i.__cg_scale"
    fields = sprintf("COALESCE(u.u%d/NULLIF(%s,0),0) AS %s", index, denominator, outputs)
    create(paste0(
      "SELECT u.__cg_op,u.__cg_group,", denominator, " AS __cg_scale,", paste(fields, collapse = ","),
      " FROM (", paste(pieces, collapse = " UNION ALL "), ") u",
      if (!is.null(info)) paste0(" JOIN ", info, " i USING (", key, ")")
    ))
  }
  rows = function(state, coefficients, outputs, subtract = FALSE) {
    ctes = vapply(seq_along(ops), function(i) paste0(
      "cg", i, " AS (SELECT __cg_group,__cg_scale,",
      paste(coefficients, collapse = ","), " FROM ", state, " WHERE __cg_op=", i, ")"
    ), character(1))
    fields = vapply(index, function(j) {
      terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
        sprintf("(%s)*COALESCE(cg%d.%s/NULLIF(cg%d.__cg_scale,0),0)", v, i, coefficients[j], i)
      }, character(1))
      value = paste(terms, collapse = " + ")
      if (subtract) value = paste0("s.", columns[j], "-(", value, ")")
      paste0("(", value, ") AS ", outputs[j])
    }, character(1))
    joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN cg", i, " ON s.", ops[[i]][["f"]], "=cg", i, ".__cg_group"), character(1))
    create(paste0(
      "WITH ", paste(ctes, collapse = ","), " SELECT ",
      paste(paste0("s.", c(identifiers, "w")), collapse = ","), ",", paste(fields, collapse = ","),
      " FROM ", table, " s", paste(joins, collapse = "")
    ))
  }
  info = moments(table, columns, r)
  checked = fetch(paste0("SELECT COUNT(*) AS n,SUM(CASE WHEN __cg_scale IS NULL OR NOT ", sql_is_finite(conn, "__cg_scale"), if (allow_zero) " OR __cg_scale<0 THEN 1 ELSE 0 END) AS bad FROM " else " OR __cg_scale<=0 THEN 1 ELSE 0 END) AS bad FROM ", info), 2L)
  if (!isTRUE(as.numeric(checked[["bad"]]) == 0) || !isTRUE(as.numeric(checked[["n"]]) > 0)) projection_abort("singular", "SQL CG has invalid group scales.")
  if (allow_zero) {
    active_info = create(paste("SELECT * FROM", info, "WHERE __cg_scale>0"))
    drop(info); info = active_info
  }
  state_rows = as.numeric(checked[["n"]])
  if (state_rows * length(columns) > limits[["design_entries"]] || 128 * state_rows * length(columns) > limits[["memory_bytes"]]) {
    projection_abort("system_limit", "SQL CG coefficient state exceeds its resource budget.")
  }
  restarts = 0L
  state = create(paste("SELECT *,", paste(c(paste(r, "AS", d), paste("0.0 AS", g)), collapse = ","), "FROM", info))
  active = rep(TRUE, length(columns))
  floor_rr = (32 * .Machine[["double.eps"]] * scales)^2
  # Restart interval. Restarting every few cycles discards the Krylov space and
  # stalls on ill-conditioned FE (W5: 5e-7 plateau at 20, exact at 500); each
  # restart also replaces the recursive residual by the certified one.
  restart_every = control[["restart"]]
  zero_steps = 0L
  for (cycle in seq_len(control[["fe_max_iter"]])) {
    direction = rows(state, d, direction_names)
    applied = moments(direction, direction_names, a, info)
    drop(direction)
    joined = create(paste0("SELECT s.*,", paste(paste0("a.", a), collapse = ","), " FROM ", state, " s JOIN ", applied, " a USING (", key, ")"))
    moment = fetch(paste("SELECT", paste(c(
      paste(sum_sql(sprintf("%1$s*%1$s", r)), "AS", paste0("rr", index)),
      paste(sum_sql(sprintf("%s*%s", d, a)), "AS", paste0("da", index))
    ), collapse = ","), "FROM", joined), 2L * length(index))
    rr = vapply(index, function(j) as.numeric(moment[[paste0("rr", j)]]), numeric(1))
    da = vapply(index, function(j) as.numeric(moment[[paste0("da", j)]]), numeric(1))
    if (any(!is.finite(rr)) || any(!is.finite(da)) || any(active & da <= 0 & rr > pmax((target * scales)^2, floor_rr))) {
      projection_abort("singular", "SQL CG lost numerical information.")
    }
    alpha = ifelse(active & da > 0 & rr > floor_rr, rr / da, 0)
    if (any(!is.finite(alpha))) projection_abort("singular", "SQL CG produced a nonfinite step.")
    updated = create(paste("SELECT", paste(c(
      key, "__cg_scale",
      sprintf("(%s+(%s)*%s) AS %s", g, number(alpha), d, g),
      sprintf("(%s-(%s)*%s) AS %s", r, number(alpha), a, r), d
    ), collapse = ","), "FROM", joined))
    drop(joined)
    drop(applied)
    drop(state)
    state = updated
    zero_steps = if (all(alpha == 0)) zero_steps + 1L else 0L
    restart = cycle %% restart_every == 0L || all(alpha == 0)
    if (cycle %% 4L == 0L || cycle == 1L || cycle == control[["fe_max_iter"]] || restart) {
      projected = rows(state, g, columns, subtract = TRUE)
      certificate = moments(projected, columns, r, info)
      row = fetch(paste("SELECT", paste(sprintf("MAX(ABS(%s))/%s AS r%d", r, number(scales), index), collapse = ","), "FROM", certificate), length(index))
      residual_columns = as.numeric(row[1, ])
      residual = max(residual_columns)
      active = residual_columns > target
      if (!is.finite(residual)) projection_abort("singular", "Nonfinite SQL CG certificate.")
      if (residual <= target) {
        drop(certificate)
        drop(state)
        drop(info)
        return(list(table = projected, cycles = cycle, residual = residual, restarts = restarts, method = "cg"))
      }
      drop(projected)
      if (zero_steps >= 2L) projection_abort("nonconverged", "SQL CG stalled above its equation tolerance.", list(residual = residual, target = target))
      if (restart) {
        restarts = restarts + 1L
        restarted = create(paste0(
          "SELECT s.__cg_op,s.__cg_group,s.__cg_scale,", paste(paste0("s.", g), collapse = ","), ",",
          paste(c(paste0("c.", r, " AS ", r), paste0("c.", r, " AS ", d)), collapse = ","),
          " FROM ", state, " s JOIN ", certificate, " c USING (", key, ")"
        ))
        drop(state)
        state = restarted
      }
      drop(certificate)
    }
    if (!restart) {
      new_rr = as.numeric(fetch(paste("SELECT", paste(paste(sum_sql(sprintf("%1$s*%1$s", r)), "AS", paste0("r", index)), collapse = ","), "FROM", state), length(index))[1, ])
      if (any(!is.finite(new_rr)) || any(new_rr < 0)) projection_abort("singular", "SQL CG lost its residual norm.")
      beta = ifelse(active & rr > floor_rr, new_rr / rr, 0)
      if (any(!is.finite(beta))) projection_abort("singular", "SQL CG produced a nonfinite direction.")
      next_state = create(paste(
        "SELECT ", paste(c("__cg_op", "__cg_group", "__cg_scale", g, r), collapse = ","), ",",
        paste(sprintf("((%s)*%s+(%s)*%s) AS %s", number(as.numeric(active)), r, number(beta), d, d), collapse = ","), "FROM", state
      ))
      drop(state)
      state = next_state
    }
  }
  projection_abort("nonconverged", "SQL CG did not converge within its iteration cap.", list(iterations = control[["fe_max_iter"]], target = target))
}

# Auto-dispatch uses this for DuckDB 3+ FE after bounded moments decline.
# The generic SQL kernel is shared with later operator/GLM consumers.
dbreg_cg_projections = function(inputs, x_sql, x_names, tol) {
  control = fe_cg_control()
  prepared = prepare_fe_projection(inputs, x_sql, x_names, materialize = TRUE, coefficient_cg = TRUE)
  if (is.null(prepared)) return(NULL)
  env2env(prepared, environment(), names(prepared))
  success = FALSE
  on.exit(if (!success) release_fit_resources(resources))
  manager = fe_cg_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  budget = list2env(list(rows = 0, cells = 0, upload_bytes = 0), parent = emptyenv())
  weight_groups = vapply(fe_columns, function(key) paste0("SELECT SUM(w) AS sw FROM ", centered, " GROUP BY ", key), character(1))
  total = fe_cg_fetch(conn, paste0("SELECT MIN(sw) AS smallest,SUM(sw)/", length(fe_columns), " AS total FROM (", paste(weight_groups, collapse = " UNION ALL "), ") cg_weights"), 2L, budget, limits)
  share = as.numeric(total[["smallest"]]) / as.numeric(total[["total"]])
  if (!is.finite(share) || share <= 0) projection_abort("singular", "SQL CG has invalid group weight share.")
  target = max(4e-15, min(tol, 1e-12) * min(1, sqrt(share)))
  identifiers = c(id_columns, "raw_w", fe_columns)
  ops = lapply(fe_columns, function(key) list(f = key, v = "1.0"))
  iterations = restarts = 0L
  for (attempt in 1:2) {
    projected = fe_sql_cg(conn, centered, columns, identifiers, ops, rep(1, length(columns)), target, control, create, drop, budget, limits)
    iterations = iterations + projected[["cycles"]]
    restarts = restarts + projected[["restarts"]]
    scaled = create(paste("SELECT", paste(c(identifiers, "w", paste0(columns, "*(", sprintf("%.17e", scale), ") AS ", columns)), collapse = ","), "FROM", projected[["table"]]))
    drop(projected[["table"]])
    certificate = fe_projection_certificate(conn, scaled, fe_columns, columns, nobs, original_scale)
    shortfall = fe_certificate_shortfall(certificate, tol)
    if (is.finite(shortfall) && shortfall <= 1) break
    drop(scaled)
    target = target / (4 * shortfall)
    if (attempt == 2L || !is.finite(shortfall) || !is.finite(target) || target < 5e-16) {
      projection_abort("certificate_failed", "SQL CG failed its independent row certificate.", c(certificate, list(target = tol)))
    }
  }
  result = create(paste("SELECT", paste(c(paste(id_columns, "AS", q(ids)), "raw_w AS __w", paste(columns, "AS", q(paste0(vars, "_tilde")))), collapse = ","), "FROM", scaled))
  drop(scaled)
  success = TRUE
  list(table = result, base_table = NULL, resources = resources, original_scale = setNames(original_scale, vars),
    iterations = iterations, accelerated_steps = 0L, max_abs = certificate[["max_abs"]], max_relative = certificate[["max_relative"]],
    method = "sql_cg", solver = "coefficient_cg", preconditioner = "diagonal", solver_iterations = iterations,
    transfer_cells = budget[["cells"]], summary_rows = budget[["rows"]], upload_bytes = 0,
    estimated_bytes = 128 * sum(sizes) * length(columns), moment_method = "rows",
    observation_tables = length(resources[["tables"]]), transformed_observation_tables = length(resources[["tables"]]) - 1L,
    peak_tables = resources[["peak_tables"]], cg_restarts = restarts, cg_retries = attempt - 1L)
}

