# A row-level warm guess can lie outside the fitted design. Project eta minus
# the known offset in prior-weight coordinates before entering the likelihood
# iteration, so damping always starts from a predictor described by the fit.
dbglm_warm_structural = function(conn, base, x) {
  p = length(x)
  if (!p) return(numeric())
  target = "(__dbglm_etastart-__dbglm_o)"
  rhs = paste(sql_stable_sum(conn, paste0("__dbglm_a*", x, "*", target)),
    "AS", paste0("b", seq_len(p)))
  row = DBI::dbGetQuery(conn, paste("SELECT", paste(c(dbglm_moment_terms(x, "__dbglm_a"), rhs),
    collapse = ","), "FROM", base))
  dbglm_solve(dbglm_read_information(row, p),
    vapply(row[paste0("b", seq_len(p))], as.numeric, numeric(1)))[["beta"]]
}

dbglm_warm_absorbed = function(conn, base, x, fe, sizes, n, control, manager, ops, lookup_manager = NULL) {
  working = paste0("(SELECT *,__dbglm_a AS w,__dbglm_a*(__dbglm_etastart-__dbglm_o) AS score",
    " FROM ", base, ") glm_warm_coordinates")
  tol = dbglm_projection_tolerance(control, strict = TRUE)
  project = if (is.null(manager)) {
    function(columns) dbglm_fe_step(conn, working, columns, fe, sizes, n, tol,
      max_iter = control[["fe_max_iter"]], lookup_manager = lookup_manager)
  } else {
    cg = dbglm_projection_cg_control(control)
    function(columns) dbglm_fe_sql_step(conn, working, columns, fe, n, tol, manager,
      FALSE, FALSE, cg, ops)
  }
  drop = function(state) {
    if (!is.null(manager)) dbglm_fe_sql_release(manager, state)
    for (table in rev(state[["lookup_tables"]])) lookup_manager[["drop"]](table)
  }
  projected = project(x)
  keep = integer()
  if (length(x)) {
    row = DBI::dbGetQuery(conn, paste("SELECT", paste(dbglm_moment_terms(x, "__dbglm_a"), collapse = ","),
      "FROM", base))
    keep = dbglm_identified_slopes(projected[["H"]], dbglm_read_information(row, length(x)),
      control[["collin_tol"]])
  }
  if (length(keep) != length(x)) {
    drop(projected)
    projected = project(x[keep])
  }
  beta = if (length(keep)) dbglm_solve(projected[["H"]], projected[["score"]])[["beta"]] else numeric()
  effects = dbglm_effect_delta(projected, beta, manager, column = "c1")
  drop(projected)
  list(keep = keep, beta = beta, effects = effects)
}
