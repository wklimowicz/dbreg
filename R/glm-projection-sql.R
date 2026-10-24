# Preconditioned CG on FE coefficients. Every row correction is reconstructed
# from group lookups, so roundoff cannot introduce arbitrary row-space effects.
# State size is O(number of FE coefficients * number of projected columns).
dbglm_project_cg = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop) {
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
      terms = c(if (is.null(info)) sprintf("SQRT(FSUM(w*(%1$s)*(%1$s))) AS __cg_scale", op[["v"]]),
        sprintf("FSUM(w*(%s)*(%s)) AS u%d", op[["v"]], inputs, index))
      paste0("SELECT ", i, " AS __cg_op,", op[["f"]], " AS __cg_group,", paste(terms, collapse = ","),
        " FROM ", source, " GROUP BY ", op[["f"]])
    }, character(1))
    denominator = if (is.null(info)) "u.__cg_scale" else "i.__cg_scale"
    fields = sprintf("COALESCE(u.u%d/NULLIF(%s,0),0) AS %s", index, denominator, outputs)
    create(paste0("SELECT u.__cg_op,u.__cg_group,", denominator, " AS __cg_scale,", paste(fields, collapse = ","),
      " FROM (", paste(pieces, collapse = " UNION ALL "), ") u",
      if (!is.null(info)) paste0(" JOIN ", info, " i USING (", key, ")")))
  }
  rows = function(state, coefficients, outputs, subtract = FALSE) {
    ctes = vapply(seq_along(ops), function(i) paste0("cg", i, " AS (SELECT __cg_group,__cg_scale,",
      paste(coefficients, collapse = ","), " FROM ", state, " WHERE __cg_op=", i, ")"), character(1))
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
    create(paste0("WITH ", paste(ctes, collapse = ","), " SELECT ",
      paste(paste0("s.", c("__dbglm_key", identifiers, "w")), collapse = ","), ",", paste(fields, collapse = ","),
      " FROM ", table, " s", paste(joins, collapse = "")))
  }
  info = moments(table, columns, r)
  state = create(paste("SELECT *,", paste(c(paste(r, "AS", d), paste("0.0 AS", g)), collapse = ","), "FROM", info))
  active = rep(TRUE, length(columns))
  floor_rr = (32 * .Machine[["double.eps"]] * scales)^2
  for (cycle in seq_len(control[["fe_max_iter"]])) {
    direction = rows(state, d, direction_names)
    applied = moments(direction, direction_names, a, info)
    drop(direction)
    joined = create(paste0("SELECT s.*,", paste(paste0("a.", a), collapse = ","), " FROM ", state, " s JOIN ", applied, " a USING (", key, ")"))
    moment = dbglm_get_query(conn, paste("SELECT", paste(c(sprintf("FSUM(%1$s*%1$s) AS rr%2$d", r, index),
      sprintf("FSUM(%s*%s) AS da%d", d, a, index)), collapse = ","), "FROM", joined))
    rr = vapply(index, function(j) as.numeric(moment[[paste0("rr", j)]]), numeric(1))
    da = vapply(index, function(j) as.numeric(moment[[paste0("da", j)]]), numeric(1))
    if (any(!is.finite(rr)) || any(!is.finite(da)) || any(active & da <= 0 & rr > pmax((target * scales)^2, floor_rr))) {
      stop("GLM fixed-effect projection lost numerical information.", call. = FALSE)
    }
    alpha = ifelse(active & da > 0 & rr > floor_rr, rr / da, 0)
    updated = create(paste("SELECT", paste(c(key, "__cg_scale",
      sprintf("(%s+(%s)*%s) AS %s", g, dbglm_number(alpha), d, g),
      sprintf("(%s-(%s)*%s) AS %s", r, dbglm_number(alpha), a, r), d), collapse = ","), "FROM", joined))
    drop(joined)
    drop(applied)
    drop(state)
    state = updated
    restart = cycle %% 20L == 0L || all(alpha == 0)
    if (cycle %% 4L == 0L || cycle == 1L || cycle == control[["fe_max_iter"]] || restart) {
      projected = rows(state, g, columns, subtract = TRUE)
      certificate = moments(projected, columns, r, info)
      row = dbglm_get_query(conn, paste("SELECT", paste(sprintf("MAX(ABS(%s))/%s AS r%d", r, dbglm_number(scales), index), collapse = ","), "FROM", certificate))
      residual_columns = as.numeric(row[1, ])
      residual = max(residual_columns)
      active = residual_columns > target
      if (isTRUE(getOption("dbreg.glm_trace")) && cycle %% 100L == 0L) message("[projection] cycle ", cycle, " residual ", residual, " target ", target)
      if (residual <= target) {
        drop(certificate)
        drop(state)
        drop(info)
        drop(table)
        return(list(table = projected, cycles = cycle, residual = residual, method = "cg"))
      }
      drop(projected)
      if (restart) {
        restarted = create(paste0("SELECT s.__cg_op,s.__cg_group,s.__cg_scale,", paste(paste0("s.", g), collapse = ","), ",",
          paste(c(paste0("c.", r, " AS ", r), paste0("c.", r, " AS ", d)), collapse = ","),
          " FROM ", state, " s JOIN ", certificate, " c USING (", key, ")"))
        drop(state)
        state = restarted
      }
      drop(certificate)
    }
    if (!restart) {
      new_rr = as.numeric(dbglm_get_query(conn, paste("SELECT", paste(sprintf("FSUM(%1$s*%1$s) AS r%2$d", r, index), collapse = ","), "FROM", state))[1, ])
      beta = ifelse(active & rr > floor_rr, new_rr / rr, 0)
      next_state = create(paste("SELECT * EXCLUDE(", paste(d, collapse = ","), "),",
        paste(sprintf("((%s)*%s+(%s)*%s) AS %s", dbglm_number(as.numeric(active)), r, dbglm_number(beta), d, d), collapse = ","), "FROM", state))
      drop(state)
      state = next_state
    }
  }
  stop("GLM fixed-effect projection did not converge.", call. = FALSE)
}

# Reuse Spark's bounded-plan window AP chunks with the GLM's scaled equation
# certificate and its tighter response-dependent inner tolerance.
dbglm_project_window = function(conn, table, columns, fe, control, create, drop, response) {
  if (!length(columns)) return(list(table = create(paste("SELECT __dbglm_key,",
    paste(c(fe, "w"), collapse = ","), "FROM", table)), cycles = 0L, residual = 0, method = "window"))
  norm = dbglm_get_query(conn, paste("SELECT", paste(sprintf("SUM(w*POWER(%s,2)) AS n%d",
    columns, seq_along(columns)), collapse = ","), "FROM", table))
  scales = sqrt(pmax(as.numeric(norm[1, ]), .Machine[["double.xmin"]]))
  target = control[["fe_tol"]]
  if (response) target = min(target, control[["score_tol"]] / (100 * max(1, utils::tail(scales, 1))))
  first = utils::tail(columns, 1)
  rest = utils::head(columns, -1)
  vars_all = c(first, rest)
  tilde_names = paste0(vars_all, "_tilde")
  id_cols = unique(c(fe, "__dbglm_key"))
  tilde_cols = sprintf("CAST((%s) AS DOUBLE) AS %s", vars_all, tilde_names)
  cur_table = create(paste0("SELECT ", paste(c(id_cols, "w AS __w", tilde_cols), collapse = ", "), " FROM ", table))
  chunk_cycles = getOption("dbreg.spark_ap_chunk_cycles", 3L)
  iter = 0L
  max_abs = Inf
  while (iter < control[["fe_max_iter"]] && max_abs > target) {
    cycles = min(chunk_cycles, control[["fe_max_iter"]] - iter)
    iter = iter + cycles
    window_sql = sql_ap_window_sweeps(cur_table, c(id_cols, "__w"), tilde_names, fe, cycles)
    alt_table = create(window_sql)
    drop(cur_table)
    cur_table = alt_table
    diagnostics = vapply(fe, function(fe_k) {
      aliases = sprintf("dbreg_mean_%04d", seq_along(tilde_names))
      means = sprintf("ABS(SUM(__w * %s))/SQRT(SUM(__w))/%.17g AS %s", tilde_names, scales[match(vars_all, columns)], aliases)
      paste0("SELECT MAX(GREATEST(", paste(aliases, collapse = ", "),
        ")) AS max_abs FROM (SELECT ", fe_k, ", ", paste(means, collapse = ", "),
        " FROM ", cur_table, " GROUP BY ", fe_k, ") dbreg_ap_conv")
    }, character(1))
    convergence_query = paste0("SELECT GREATEST(", paste(sprintf("(%s)", diagnostics), collapse = ", "), ") AS max_abs")
    max_abs = as.numeric(dbglm_get_query(conn, convergence_query)[1, 1L])
  }
  result = create(paste("SELECT __dbglm_key,", paste(c(fe, "__w AS w",
    sprintf("%s_tilde AS %s", columns, columns)), collapse = ","), "FROM", cur_table))
  drop(cur_table)
  list(table = result, cycles = iter, residual = max_abs, method = "window")
}
