# SQL certificates for iterative fixed-effect projection. Absolute group
# means alone can accept a decaying FE-only regressor before its rank is clear.
sql_ap_original_scale = function(table, vars, conn) {
  values = sql_cast_double(vars, conn)
  paste("SELECT", paste(sprintf("(MAX(%s)/2-MIN(%s)/2) AS original_%d", values, values, seq_along(vars)), collapse = ","), "FROM", table)
}

# FE intercepts absorb a global location. Remove it before the first sweep,
# and measure projection roundoff in those same centered coordinates.
sql_ap_initialization = function(source, vars, ids, conn) {
  centers = safe_identifiers(paste0("dbreg_ap_center_", seq_along(vars)),
    taken = c(vars, ids, "__w", paste0(vars, "_tilde")))
  values = sql_cast_double(vars, conn)
  location = paste0("(SELECT ", paste(sprintf("(MIN(%s)/2+MAX(%s)/2) AS %s",
    values, values, centers), collapse = ","), " FROM ", source, " dbreg_ap_location_source) l")
  projected = sprintf("(%s-l.%s) AS %s_tilde", sql_cast_double(paste0("b.", vars), conn), centers, vars)
  paste("SELECT", paste(c(paste0("b.", ids), "b.__w", projected), collapse = ","),
    "FROM", source, "b CROSS JOIN", location)
}

sql_ap_rms = function(table, vars, conn) {
  values = sql_cast_double(paste0(vars, "_tilde"), conn)
  paste("SELECT", paste(sprintf("SQRT(SUM(__w*(%s)*(%s))/SUM(__w)) AS rms_%d", values, values, seq_along(vars)), collapse = ","), "FROM", table)
}

sql_ap_group_means = function(table, vars, fe) {
  terms = sprintf("SUM(__w * %s_tilde) / SUM(__w) AS %s_mean", vars, vars)
  paste0(
    "SELECT ", paste(c(fe, terms), collapse = ", "),
    " FROM ", table, " GROUP BY ", fe
  )
}

sql_ap_join_projection = function(table, means, ids, vars, fe, snapshot = NULL) {
  cols = c(
    paste0("t.", ids), "t.__w",
    if (!is.null(snapshot)) sprintf("t.%s_tilde AS %s_%s", vars, vars, snapshot),
    sprintf("t.%s_tilde - m.%s_mean AS %s_tilde", vars, vars, vars)
  )
  paste0(
    "SELECT ", paste(cols, collapse = ", "), " FROM ", table,
    " t JOIN ", means, " m ON t.", fe, " = m.", fe
  )
}

# Irons-Tuck acceleration of SQL alternating projections, as in fixest's
# demeaning. With X (`_ita`), F(X) (`_itb`) and F(F(X)) (`_tilde`) kept as
# columns, alpha = <dF, d2>_w / <d2, d2>_w per column, where dF = F(F(X)) - F(X)
# and d2 = F(F(X)) - 2 F(X) + X, and X <- F(F(X)) - alpha dF. Every sweep
# subtracts FE-span vectors, so the extrapolation stays in X0 + span(FE) and
# the usual group-mean test still certifies the projection.
# options(dbreg.ap_accelerate = FALSE) restores plain AP.
ap_accelerate = function() {
  value = getOption("dbreg.ap_accelerate", TRUE)
  if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("`dbreg.ap_accelerate` must be TRUE or FALSE.", call. = FALSE)
  value
}

sql_ap_irons_tuck_moments = function(table, tilde_names) {
  base = sub("_tilde$", "", tilde_names)
  d2 = sprintf("(%1$s_tilde - 2 * %1$s_itb + %1$s_ita)", base)
  df = sprintf("(%1$s_tilde - %1$s_itb)", base)
  index = seq_along(base)
  paste0("SELECT ", paste(c(
    sprintf("SUM(__w * %s * %s) AS it_num_%d", df, d2, index),
    sprintf("SUM(__w * %s * %s) AS it_den_%d", d2, d2, index)
  ), collapse = ", "), " FROM ", table)
}

ap_irons_tuck_alpha = function(row, n) {
  vapply(seq_len(n), function(i) {
    num = as.numeric(row[[paste0("it_num_", i)]])
    den = as.numeric(row[[paste0("it_den_", i)]])
    alpha = if (is.finite(num) && is.finite(den) && den > 0) num / den else 0
    if (is.finite(alpha)) alpha else 0
  }, numeric(1))
}

sql_ap_irons_tuck_step = function(table, keep_cols, tilde_names, alpha) {
  base = sub("_tilde$", "", tilde_names)
  paste0("SELECT ", paste(c(keep_cols, sprintf(
    "%1$s_tilde - (%2$.17e) * (%1$s_tilde - %1$s_itb) AS %1$s_tilde", base, alpha
  )), collapse = ", "), " FROM ", table)
}

# Two window AP cycles in one statement, keeping X and F(X) for Irons-Tuck.
sql_ap_window_irons_tuck = function(table, keep_cols, tilde_names, fe) {
  base = sub("_tilde$", "", tilde_names)
  first = paste0("(SELECT *, ", paste(sprintf("%s_tilde AS %s_ita", base, base), collapse = ", "), " FROM ", table, ") it_a")
  one = sql_ap_window_sweeps(first, c(keep_cols, paste0(base, "_ita")), tilde_names, fe, 1L)
  second = paste0("(SELECT *, ", paste(sprintf("%s_tilde AS %s_itb", base, base), collapse = ", "), " FROM (", one, ") it_one) it_b")
  sql_ap_window_sweeps(second, c(keep_cols, paste0(base, "_ita"), paste0(base, "_itb")), tilde_names, fe, 1L)
}

# Validate before creating connections/tables; every estimator shares the same
# per-fit controls and leaves the caller's global options alone.
validate_fe_controls = function(inputs) {
  tol = inputs[["fe_tol"]]
  iterations = inputs[["fe_max_iter"]]
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`fe_tol` must be one finite, positive number.", call. = FALSE)
  }
  if (!is.numeric(iterations) || length(iterations) != 1L || !is.finite(iterations) ||
      iterations < 1 || iterations > .Machine[["integer.max"]] || iterations != floor(iterations)) {
    stop("`fe_max_iter` must be a positive integer.", call. = FALSE)
  }
  inputs[["fe_max_iter"]] = as.integer(iterations)
}

sql_ap_convergence = function(table, vars, fe) {
  means = paste(sprintf("SUM(__w * %s_tilde) / SUM(__w) AS %s_mean", vars, vars), collapse = ", ")
  inner = paste(sprintf("SELECT %s FROM %s GROUP BY %s", means, table, fe), collapse = " UNION ALL ")
  outer = paste(sprintf("MAX(ABS(%s_mean)) AS max_%s", vars, vars), collapse = ", ")
  paste0("SELECT ", outer, " FROM (", inner, ") dbreg_convergence")
}

# An inspectable ordinary-AP program. The execution path shares its projection,
# certificate and moments builders, rather than maintaining a second dialect.
dbreg_ap_query_plan = function(inputs, x_sql, x_names) {
  conn = inputs[["conn"]]
  fe = inputs[["fe"]]
  y = inputs[["yvar"]]
  vars = c(y, x_names)
  ids = projection_id_columns(inputs)
  w = sql_weight_expr(inputs[["weights"]])
  if (is.null(w)) w = "1.0"
  backend = detect_backend(conn)[["name"]]
  method = ap_method(backend)
  base = paste0("SELECT ", paste(c(
    frame_select(y, x_sql, x_names, id_cols = setdiff(ids, vars)),
    paste0(w, " AS __w")
  ), collapse = ", "), " ", inputs[["from_statement"]])
  tilde = paste0(vars, "_tilde")
  initial = sql_ap_initialization("(<model_frame_sql>)", vars, ids, conn)
  projection = setNames(vapply(fe, function(f) {
    if (method == "window") return(sql_ap_window_sweeps("<ap_current>", c(ids, "__w"), tilde, f, 1L))
    paste0("WITH dbreg_means AS (", sql_ap_group_means("<ap_current>", vars, f), ") ",
      sql_ap_join_projection("<ap_current>", "dbreg_means", ids, vars, f))
  }, character(1)), fe)
  frame = list(cte_sql = "WITH demeaned AS (SELECT * FROM <ap_converged>)", cte_name = "demeaned",
    cols = setNames(tilde, vars), weights_expr = "__w", intercept = FALSE, cast = TRUE,
    rank_source = paste0("(", base, ") dbreg_rank_source"),
    prefix_terms = c(sql_count(conn, "n_total"), vapply(seq_along(fe), function(k)
      sql_count(conn, sprintf("n_fe%d", k), fe[[k]], distinct = TRUE), character(1))))
  structure(list(strategy = "demean", backend = backend, method = method,
    tolerance = inputs[["fe_tol"]], max_iterations = inputs[["fe_max_iter"]],
    model_frame_sql = base, initialization_sql = initial,
    design_matrix = data.frame(alias = x_names, sql = x_sql),
    projection_sql = projection, convergence_sql = sql_ap_convergence("<ap_current>", vars, fe),
    original_scale_sql = sql_ap_original_scale("(<model_frame_sql>) dbreg_original", vars, conn),
    residual_scale_sql = sql_ap_rms("<ap_current>", vars, conn),
    final_moments_sql = frame_moments_sql(frame, y, x_names, conn),
    note = paste("Ordinary iterative AP templates; replace angle-bracket placeholders.",
      "Repeat all projections until absolute and RMS-relative group means meet tolerance,",
      "except columns below 128 * double epsilon * their original centered half-range.",
      "Final moments carry the same original-scale rank guard.",
      "Automatic fitting may use acceleration; these templates describe the ordinary equivalent.",
      "Use a caller-owned database table when replaying the plan; package-owned sources may expire.")),
    class = "dbreg_query_plan")
}

#' Print an iterative fixed-effect SQL plan
#' @aliases dbreg_query_plan
#' @details
#' `dbreg()` and `dbivreg()` return this named list for iterative fixed-effect
#' models with `sql_only = TRUE`:
#' - `model_frame_sql` selects the fitted sample and expanded design.
#' - `initialization_sql` centers working columns in double precision.
#' - `projection_sql` contains one ordinary sweep per fixed effect, in order.
#' - `convergence_sql` returns absolute group means; `original_scale_sql` and
#'   `residual_scale_sql` supply the roundoff and relative-convergence checks.
#' - `final_moments_sql` aggregates a converged projection, preserving rank guards.
#' - `design_matrix` maps SQL aliases to expressions; `tolerance` and
#'   `max_iterations` record the per-fit controls.
#'
#' Replace the angle-bracket placeholders with SQL or table names. Reject
#' nonfinite certificates; otherwise use the absolute and RMS-relative checks
#' described in [dbreg-options]. These are ordinary-AP templates, while automatic
#' fitting may accelerate. A caller-owned source table is needed for durable
#' replay because package-owned registrations can expire. No AP tables are
#' created while inspecting a plan.
#' @md
#' @param x A `dbreg_query_plan` object.
#' @param ... Unused.
#' @export
print.dbreg_query_plan = function(x, ...) {
  cat("Iterative fixed-effect SQL plan\n")
  cat("Backend:", x[["backend"]], "  Method:", x[["method"]], "\n")
  cat("Tolerance:", x[["tolerance"]], "  Cycle budget:", x[["max_iterations"]], "\n")
  cat(x[["note"]], "\n")
  invisible(x)
}
