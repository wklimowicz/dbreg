# Permit only roundoff-sized deviance increases in either IRLS implementation.
dbglm_deviance_allowance = function(deviance) {
  64 * .Machine[["double.eps"]] * max(1, abs(deviance))
}

dbglm_moment_terms = function(x, weight) {
  terms = character()
  for (i in seq_along(x)) for (j in seq_len(i)) {
    terms = c(terms, sprintf("SUM((%s)*(%s)*(%s)) AS h_%d_%d", weight, x[i], x[j], i, j))
  }
  terms
}

dbglm_read_information = function(row, p) {
  H = matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(i)) {
    H[i, j] = H[j, i] = row[[sprintf("h_%d_%d", i, j)]]
  }
  H
}

dbglm_basis = function(H, tol) {
  scale = sqrt(pmax(diag(H), 0))
  keep = integer()
  for (j in seq_len(ncol(H))) {
    if (!is.finite(scale[j]) || scale[j] == 0) next
    indices = c(keep, j)
    C = H[indices, indices, drop = FALSE] / outer(scale[indices], scale[indices])
    R = tryCatch(chol(C), error = function(e) NULL)
    if (!is.null(R) && utils::tail(diag(R), 1)^2 > tol) keep = indices
  }
  if (!length(keep)) stop("No identified GLM design columns remain.", call. = FALSE)
  keep
}

dbglm_solve = function(H, b = NULL) {
  scale = sqrt(diag(H))
  if (any(!is.finite(H)) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("GLM working information is invalid; check identification/separation.", call. = FALSE)
  }
  C = H / outer(scale, scale)
  R = tryCatch(chol(C), error = function(e) NULL)
  if (is.null(R)) stop("GLM working information is singular; check separation.", call. = FALSE)
  inverse = chol2inv(R) / outer(scale, scale)
  beta = if (is.null(b)) NULL else as.numeric(backsolve(R, forwardsolve(t(R), b / scale))) / scale
  list(beta = beta, inverse = inverse)
}

dbglm_evaluate = function(conn, table, x, beta, link) {
  cte = dbglm_working_sql(table, x, beta, link)
  terms = dbglm_moment_terms(x, "w")
  terms = c(terms, sprintf("SUM(w*(%s)*z) AS b_%d", x, seq_along(x)),
    sprintf("%s AS s_%d", vapply(sprintf("(%s)*score", x), function(e) sql_stable_sum(conn, e), character(1)), seq_along(x)),
    dbglm_deviance_term(conn), "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
    "SUM(CASE WHEN NOT isfinite(w) OR w <= 0 OR NOT isfinite(eta) THEN 1 ELSE 0 END) AS bad")
  row = backend_get_query(conn, paste0(cte, "SELECT ", paste(terms, collapse = ", "), " FROM working"))
  H = dbglm_read_information(row, length(x))
  b = vapply(seq_along(x), function(i) as.numeric(row[[paste0("b_", i)]]), numeric(1))
  score = vapply(seq_along(x), function(i) as.numeric(row[[paste0("s_", i)]]), numeric(1))
  list(H = H, b = b, score = score, score_max = max(abs(score) / sqrt(diag(H))),
    deviance = row[["deviance"]], bad = row[["bad"]], mu_min = row[["mu_min"]], mu_max = row[["mu_max"]])
}

dbglm_irls = function(conn, table, x, beta, family, control, verbose, separation = "error") {
  current = dbglm_evaluate(conn, table, x, beta, family[["link"]])
  initial_information = sqrt(diag(current[["H"]]))
  history = current[["deviance"]]
  halvings_total = 0L
  for (iter in seq_len(control[["glm_max_iter"]])) {
    candidate = dbglm_solve(current[["H"]], current[["b"]])[["beta"]]
    step = 1
    accepted = FALSE
    for (halvings in 0:control[["step_max_halvings"]]) {
      next_beta = beta + step * (candidate - beta)
      next_state = dbglm_evaluate(conn, table, x, next_beta, family[["link"]])
      # The no-FE path allows boundary trials for separation = "prediction".
      # Its acceptance policy deliberately differs from the absorbed FE fit.
      allowance = dbglm_deviance_allowance(current[["deviance"]])
      if (next_state[["bad"]] == 0 && is.finite(next_state[["deviance"]]) &&
        next_state[["deviance"]] <= current[["deviance"]] + allowance) {
        accepted = TRUE
        break
      }
      step = step / 2
    }
    if (!accepted) stop("GLM step halving failed to find a valid improving step.", call. = FALSE)
    relative = abs(current[["deviance"]] - next_state[["deviance"]]) / (0.1 + abs(next_state[["deviance"]]))
    beta = next_beta
    current = next_state
    history = c(history, current[["deviance"]])
    halvings_total = halvings_total + halvings
    if (isTRUE(verbose)) message("[dbglm] iteration ", iter, ": deviance=",
      signif(current[["deviance"]], 8), ", score=", signif(current[["score_max"]], 3))
    boundary = dbglm_boundary(current, family[["link"]])
    if (step == 1 && relative <= control[["glm_tol"]] && boundary &&
      max(abs(current[["score"]]) / initial_information) <= control[["score_tol"]]) {
      if (separation == "error") stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      warning("Returning boundary probabilities only: coefficients are not a finite MLE; inference is unavailable.", call. = FALSE)
      return(c(current, list(beta = beta, iterations = iter, deviance_history = history,
        step_halvings = halvings_total, finite_mle = FALSE, converged = FALSE)))
    }
    if (step == 1 && (family[["link"]] == "identity" || relative <= control[["glm_tol"]]) &&
      is.finite(current[["score_max"]]) && current[["score_max"]] <= control[["score_tol"]]) {
      if (dbglm_boundary(current, family[["link"]])) {
        stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      }
      return(c(current, list(beta = beta, iterations = iter,
        deviance_history = history, step_halvings = halvings_total)))
    }
  }
  stop("GLM did not converge within ", control[["glm_max_iter"]], " iterations.", call. = FALSE)
}
