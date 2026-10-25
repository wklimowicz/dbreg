# Shared by iterative execution and the inspectable SQL plan.
sql_ap_group_means = function(table, vars, fe) {
  terms = sprintf("SUM(__w * %s_tilde) / SUM(__w) AS %s_mean", vars, vars)
  paste0("SELECT ", paste(c(fe, terms), collapse = ", "),
    " FROM ", table, " GROUP BY ", fe)
}

sql_ap_join_projection = function(table, means, ids, vars, fe) {
  cols = c(paste0("t.", ids), "t.__w",
    sprintf("t.%s_tilde - m.%s_mean AS %s_tilde", vars, vars, vars))
  paste0("SELECT ", paste(cols, collapse = ", "), " FROM ", table,
    " t JOIN ", means, " m ON t.", fe, " = m.", fe)
}

sql_ap_convergence = function(table, vars, fe) {
  means = paste(sprintf("SUM(__w * %s_tilde) / SUM(__w) AS %s_mean", vars, vars), collapse = ", ")
  inner = paste(sprintf("SELECT %s FROM %s GROUP BY %s", means, table, fe), collapse = " UNION ALL ")
  outer = paste(sprintf("MAX(ABS(%s_mean)) AS max_%s", vars, vars), collapse = ", ")
  paste0("SELECT ", outer, " FROM (", inner, ") dbreg_convergence")
}

sql_ap_original_scale = function(table, vars, conn) {
  paste("SELECT", paste(sprintf("MAX(ABS(%s)) AS original_%d", sql_cast_double(vars, conn), seq_along(vars)), collapse = ","), "FROM", table)
}

sql_ap_rms = function(table, vars, conn) {
  values = sql_cast_double(paste0(vars, "_tilde"), conn)
  paste("SELECT", paste(sprintf("SQRT(SUM(__w*(%s)*(%s))/SUM(__w)) AS rms_%d", values, values, seq_along(vars)), collapse = ","), "FROM", table)
}

dbreg_ap_query_plan = function(inputs, x_sql, x_names) {
  conn = inputs[["conn"]]
  fe = inputs[["fe"]]
  y = inputs[["yvar"]]
  vars = c(y, x_names)
  ids = unique(c(fe, inputs[["cluster_var"]]))
  w = sql_weight_expr(inputs[["weights"]])
  if (is.null(w)) w = "1.0"
  backend = detect_backend(conn)[["name"]]
  method = ap_method(backend)
  base = paste0("SELECT ", paste(c(frame_select(y, x_sql, x_names,
    id_cols = setdiff(ids, vars)), paste0(w, " AS __w")), collapse = ", "),
    " ", inputs[["from_statement"]])
  tilde = paste0(vars, "_tilde")
  initial = paste0("SELECT ", paste(c(ids, "__w", sprintf("%s AS %s",
    if (method == "window") sql_cast_double(vars, conn) else vars, tilde)), collapse = ", "),
    " FROM (<model_frame_sql>) dbreg_ap_base")
  projection = setNames(vapply(fe, function(f) {
    if (method == "window") return(sql_ap_window_sweeps("<ap_current>", c(ids, "__w"), tilde, f, 1L))
    paste0("WITH dbreg_means AS (", sql_ap_group_means("<ap_current>", vars, f), ") ",
      sql_ap_join_projection("<ap_current>", "dbreg_means", ids, vars, f))
  }, character(1)), fe)
  frame = list(cte_sql = "WITH demeaned AS (SELECT * FROM <ap_converged>)", cte_name = "demeaned",
    cols = setNames(tilde, vars), weights_expr = "__w", intercept = FALSE, cast = TRUE,
    rank_source = paste0("(", base, ") dbreg_ap_rank_source"),
    prefix_terms = c(sql_count(conn, "n_total"), vapply(seq_along(fe), function(k)
      sql_count(conn, sprintf("n_fe%d", k), fe[[k]], distinct = TRUE), character(1))))
  structure(list(strategy = "demean", backend = backend, method = method,
    tolerance = inputs[["fe_tol"]], max_iterations = inputs[["fe_max_iter"]],
    model_frame_sql = base, initialization_sql = initial,
    design_matrix = data.frame(term = restore_names(standardize_coef_names(x_names), inputs[["name_map"]]), sql = x_sql),
    projection_sql = projection, convergence_sql = sql_ap_convergence("<ap_current>", vars, fe),
    original_scale_sql = sql_ap_original_scale("(<model_frame_sql>) dbreg_ap_original", vars, conn),
    residual_scale_sql = sql_ap_rms("<ap_current>", vars, conn),
    final_moments_sql = frame_moments_sql(frame, y, x_names, conn),
    note = paste("Templates for iterative fixed-effect projection; replace the angle-bracket placeholders.",
      "Repeat the projections until convergence. The automatic solver may instead use bounded FE moments.",
      "Also require group means divided by projected RMS to meet tolerance for columns above their original-value roundoff bound.",
      "No projection or coefficient estimation is executed by sql_only.")), class = "dbreg_query_plan")
}

#' Inspect a fixed-effect projection SQL plan
#' @param x A plan returned by `dbreg()` or `dbivreg()` with `sql_only = TRUE`.
#' @param ... Unused.
#' @export
print.dbreg_query_plan = function(x, ...) {
  cat("Fixed-effect SQL plan (", x[["backend"]], ", ", x[["method"]], ")\n", sep = "")
  for (name in c("model_frame_sql", "initialization_sql", "projection_sql", "convergence_sql", "original_scale_sql", "residual_scale_sql", "final_moments_sql")) {
    cat("\n", name, ":\n", paste(x[[name]], collapse = "\n\n"), "\n", sep = "")
  }
  cat("\n", x[["note"]], "\n", sep = "")
  invisible(x)
}
