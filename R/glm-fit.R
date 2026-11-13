# Permit only roundoff-sized deviance increases in either IRLS implementation.
dbglm_deviance_allowance = function(deviance) {
  64 * .Machine[["double.eps"]] * max(1, abs(deviance))
}

# A score certificate cannot be made small by invalid or singular information.
dbglm_score_norm = function(score, H) {
  if (!isTRUE(nrow(H) == length(score)) || !isTRUE(ncol(H) == length(score)) ||
      any(!is.finite(H)) || any(!is.finite(score))) return(Inf)
  if (!length(score)) return(0)
  scale = sqrt(pmax(diag(H), 0))
  if (any(scale <= 0)) return(Inf)
  if (is.null(tryCatch(chol(H / outer(scale, scale)), error = function(e) NULL))) return(Inf)
  max(abs(score) / scale)
}

# Gaussian scores carry the outcome's units; measure them against the residual
# scale so large-unit outcomes (prices, revenue) are certified like small ones.
dbglm_score_unit = function(link, deviance, sum_w) {
  if (link != "identity" || !is.finite(deviance) || !is.finite(sum_w) || sum_w <= 0) return(1)
  max(1, sqrt(max(deviance, 0) / sum_w))
}

dbglm_state_unit = function(state) {
  unit = state[["score_unit"]]
  if (is.numeric(unit) && length(unit) == 1L && is.finite(unit) && unit > 0) unit else 1
}

# Full steps need a current certificate. A damped stationary step additionally
# needs the previous independent certificate, never merely a small objective change.
dbglm_converged_step = function(step, relative, link, score, previous_score, control) {
  certified = function(value) is.numeric(value) && length(value) == 1L &&
    is.finite(value) && value >= 0 && value <= control[["score_tol"]]
  (isTRUE(step == 1) || certified(previous_score)) && certified(score) &&
    (link == "identity" || isTRUE(is.finite(relative) && relative <= control[["glm_tol"]]))
}

# Express coefficient corrections in fixed design coordinates. Neither raw
# coefficient units nor weak starting information may determine convergence.
dbglm_step_norm = function(correction, scale) {
  if (is.null(correction) || length(correction) != length(scale) ||
      any(!is.finite(correction))) return(Inf)
  if (length(correction)) max(abs(correction * scale)) else 0
}

# Certify a recession direction on the observed outcome support. Equality on
# positive Poisson counts (or mixed binomial cells) must be literal: a tiny
# nonzero design coordinate can still restore a finite optimum arbitrarily far
# along the proposed ray. Failure to certify leaves the step guard active.
dbglm_direction_boundary = function(conn, table, link, direction) {
  conditions = if (link == "log") {
    c(paste0("__dbglm_y>0 AND (", direction, ")<>0"),
      paste0("__dbglm_y=0 AND (", direction, ")>0"))
  } else if (link %in% c("logit", "probit", "cloglog")) {
    c(paste0("__dbglm_y>0 AND (", direction, ")<0"),
      paste0("__dbglm_y<1 AND (", direction, ")>0"))
  } else return(FALSE)
  row = DBI::dbGetQuery(conn, paste0("SELECT MAX(ABS(", direction, ")) AS amplitude,",
    "SUM(CASE WHEN ", paste(conditions, collapse = " OR "),
    " THEN 1 ELSE 0 END) AS violations FROM ", table))
  isTRUE(is.finite(row[["amplitude"]]) && row[["amplitude"]] > 0 && row[["violations"]] == 0)
}

dbglm_poisson_direction = function(conn, table, x, correction) {
  if (!length(x) || is.null(correction) || any(!is.finite(correction))) return(FALSE)
  value = paste(sprintf("(%s)*(%s)", dbglm_number(correction), x), collapse = "+")
  dbglm_direction_boundary(conn, table, "log", value)
}

dbglm_moment_terms = function(x, weight) {
  terms = character()
  for (i in seq_along(x)) for (j in seq_len(i)) {
    terms = c(terms, sprintf("SUM((%s)*(%s)*(%s)) AS h_%d_%d", weight, x[i], x[j], i, j))
  }
  terms
}

dbglm_read_information = function(row, p, prefix = "h") {
  H = matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(i)) {
    H[i, j] = H[j, i] = row[[sprintf("%s_%d_%d", prefix, i, j)]]
  }
  H
}

dbglm_basis = function(H, tol) {
  if (any(!is.finite(H))) stop("GLM initial information is nonfinite; check predictor magnitudes and weights.", call. = FALSE)
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

dbglm_evaluate = function(conn, table, x, beta, link, eta = NULL, cell = FALSE,
                          materialize = FALSE) {
  backend_info = detect_backend(conn)
  finite = sql_is_finite(conn, c("w", "eta"), backend_info)
  cte = dbglm_working_sql(table, x, beta, link, eta = eta, conn = conn,
    backend_info = backend_info, cell = cell, materialize = materialize)
  terms = c(
    dbglm_moment_terms(x, "w"),
    sprintf("%s AS s_%d", sql_stable_sum(conn, sprintf("(%s)*score", x),
      backend_info), seq_along(x)),
    dbglm_deviance_term(conn), "SUM(w) AS sum_w", "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
    "MIN(eta) AS eta_min", "MAX(eta) AS eta_max", dbglm_misfit_sql(),
    paste0("SUM(CASE WHEN NOT ", finite[[1L]], " OR w <= 0 OR NOT ",
      finite[[2L]], " THEN 1 ELSE 0 END) AS bad")
  )
  row = DBI::dbGetQuery(conn, paste0(cte, "SELECT ", paste(terms, collapse = ", "), " FROM working"))
  H = dbglm_read_information(row, length(x))
  score = vapply(seq_along(x), function(i) as.numeric(row[[paste0("s_", i)]]),
    numeric(1))
  unit = dbglm_score_unit(link, as.numeric(row[["deviance"]]), as.numeric(row[["sum_w"]]))
  list(
    H = H, score = score, score_max = dbglm_score_norm(score, H) / unit, score_unit = unit,
    deviance = row[["deviance"]], bad = row[["bad"]],
    mu_min = row[["mu_min"]], mu_max = row[["mu_max"]],
    eta_min = row[["eta_min"]], eta_max = row[["eta_max"]],
    misfit = row[["misfit"]]
  )
}

dbglm_irls = function(conn, table, x, beta, family, control,
                       verbose = FALSE, separation = "error", cell = FALSE, anchor = numeric(length(beta)), step_scale = rep(1, length(x))) {
  link = family[["link"]]
  initial_halvings = 0L
  repeat {
    current = dbglm_evaluate(conn, table, x, beta, link, cell = cell)
    if (dbglm_valid_start(current, link)) break
    if (initial_halvings >= control[["step_max_halvings"]]) {
      stop("Initial GLM predictor reaches a numerical boundary; choose smaller starts or offsets.", call. = FALSE)
    }
    beta = beta / 2 + anchor / 2
    initial_halvings = initial_halvings + 1L
  }
  initial_information = sqrt(diag(current[["H"]]))
  initial_factor = dbglm_info_factor(current[["H"]])
  history = current[["deviance"]]
  halvings_total = initial_halvings
  needs_step_guard = link %in% c("log", "logit", "probit", "cloglog")
  prev_step_inf = Inf
  next_correction = NULL
  for (iter in seq_len(control[["glm_max_iter"]])) {
    previous_score = current[["score_max"]]
    # Solve for the coefficient correction from the score to avoid cancellation
    # between absolute working-response solutions near convergence.
    correction = if (!is.null(next_correction)) next_correction else dbglm_solve(current[["H"]], current[["score"]])[["beta"]]
    step = 1
    accepted = FALSE
    for (halvings in 0:control[["step_max_halvings"]]) {
      next_beta = beta + step * correction
      next_state = dbglm_evaluate(conn, table, x, next_beta, link, cell = cell)
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
    if (!accepted) {
      if (isTRUE(control[["tolerate_max_iter"]])) return(c(current, list(
        beta = beta, iterations = iter, deviance_history = history,
        step_halvings = halvings_total + halvings, converged = FALSE, stop_reason = "step_halving_failed"
      )))
      stop("GLM step halving failed to find a valid improving step.", call. = FALSE)
    }
    relative = abs(current[["deviance"]] - next_state[["deviance"]]) / (0.1 + abs(next_state[["deviance"]]))
    beta = next_beta
    current = next_state
    history = c(history, current[["deviance"]])
    halvings_total = halvings_total + halvings
    if (isTRUE(verbose)) message(
      "[dbglm] iteration ", iter, ": deviance=",
      signif(current[["deviance"]], 8), ", score=", signif(current[["score_max"]], 3)
    )
    next_correction = tryCatch(dbglm_solve(current[["H"]], current[["score"]])[["beta"]], error = function(e) NULL)
    next_step_inf = dbglm_step_norm(next_correction, step_scale)
    eig_ratio = dbglm_min_info_ratio(current[["H"]], initial_factor[["R"]], initial_factor[["scale"]])
    collapsed = isTRUE(eig_ratio <= 1e-8)
    deep_collapsed = isTRUE(eig_ratio <= 1e-11)
    stalled = needs_step_guard && step == 1 &&
      (is.null(next_correction) || (next_step_inf > 1e-4 && next_step_inf >= 0.5 * prev_step_inf))
    boundary = dbglm_boundary(current, link, collapsed = collapsed, stalled = stalled || is.null(next_correction))
    if (link == "log" && stalled && relative <= control[["glm_tol"]] &&
        current[["score_max"]] <= control[["score_tol"]]) {
      boundary = boundary || dbglm_poisson_direction(conn, table, x, next_correction)
    }
    score_h0 = if (length(initial_information)) max(abs(current[["score"]]) / initial_information) else 0
    if (step == 1 && relative <= control[["glm_tol"]] && boundary &&
      (current[["score_max"]] <= control[["score_tol"]] || (deep_collapsed && score_h0 <= control[["score_tol"]]))) {
      if (separation == "error") stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      warning("Returning boundary probabilities only: coefficients are not a finite MLE; inference is unavailable.", call. = FALSE)
      return(c(current, list(
        beta = beta, iterations = iter, deviance_history = history,
        step_halvings = halvings_total, finite_mle = FALSE, converged = FALSE
      )))
    }
    step_small = !needs_step_guard || (!is.null(next_correction) && next_step_inf <= 1e-4)
    if (step_small && dbglm_converged_step(step, relative, link, current[["score_max"]], previous_score, control)) {
      if (boundary) {
        stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      }
      return(c(current, list(
        beta = beta, iterations = iter, converged = TRUE, stop_reason = "converged",
        deviance_history = history, step_halvings = halvings_total
      )))
    }
    prev_step_inf = if (step == 1) next_step_inf else Inf
  }
  if (isTRUE(control[["tolerate_max_iter"]])) return(c(current, list(
    beta = beta, iterations = control[["glm_max_iter"]], deviance_history = history,
    step_halvings = halvings_total, converged = FALSE, stop_reason = "max_iter"
  )))
  dbglm_nonconvergence(current, link, control[["glm_max_iter"]])
}

dbglm_covariance = function(conn, table, x, fit, family, kind, ssc, n, K = length(x), cell = FALSE) {
  inverse = if (length(x)) dbglm_solve(fit[["H"]])[["inverse"]] else matrix(numeric(), 0L, 0L)
  p = length(x)
  if (ssc[["K.adj"]] && n <= K) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  if (kind == "iid") {
    V = inverse
    if (ssc[["K.adj"]]) V = V * (n - 1) / (n - K)
    return(list(matrix = V, clusters = NULL))
  }
  source = fit[["projected"]]
  if (is.null(source)) {
    cte = dbglm_working_sql(table, x, fit[["beta"]], family[["link"]], conn = conn, cell = cell)
    source = paste0("(", cte, "SELECT * FROM working) glm_covariance_rows")
  } else {
    x = if (p) paste0("v", seq_len(p) + 1L) else character()
  }
  if (!p) {
    G = if (kind == "cluster") as.numeric(DBI::dbGetQuery(conn,
      paste("SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM", source))[["n"]]) else NULL
    if (length(G) && G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
    return(list(matrix = inverse, clusters = G))
  }
  if (kind == "hc1") {
    row = DBI::dbGetQuery(conn, paste0(
      "SELECT ",
      paste(dbglm_moment_terms(x, if (cell) "score_sq" else "score*score"), collapse = ", "), " FROM ", source
    ))
    G = NULL
  } else {
    scores = sprintf("%s AS s_%d", sql_stable_sum(conn, sprintf("(%s)*score", x)),
      seq_along(x))
    terms = dbglm_moment_terms(paste0("s_", seq_along(x)), "1.0")
    row = DBI::dbGetQuery(conn, paste0(
      "WITH cluster_scores AS (SELECT __dbglm_cluster, ", paste(scores, collapse = ", "),
      " FROM ", source, " GROUP BY __dbglm_cluster) SELECT COUNT(*) AS n_clusters, ",
      paste(terms, collapse = ", "), " FROM cluster_scores"
    ))
    G = as.numeric(row[["n_clusters"]])
    if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
  }
  meat = dbglm_read_information(row, p)
  V = inverse %*% meat %*% inverse
  if (ssc[["K.adj"]]) V = V * (if (kind == "hc1") n else n - 1) / (n - K)
  if (kind == "cluster" && ssc[["G.adj"]]) V = V * G / (G - 1)
  list(matrix = V, clusters = G)
}
