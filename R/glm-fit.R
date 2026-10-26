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

dbglm_evaluate = function(conn, table, x, beta, link, cell = FALSE) {
  backend_info = detect_backend(conn)
  finite = sql_is_finite(conn, c("w", "eta"), backend_info)
  cte = dbglm_working_sql(table, x, beta, link, cell = cell)
  terms = dbglm_moment_terms(x, "w")
  terms = c(
    terms, sprintf("SUM(w*(%s)*z) AS b_%d", x, seq_along(x)),
    sprintf("%s AS s_%d", vapply(sprintf("(%s)*score", x), function(e) sql_stable_sum(conn, e, backend_info), character(1)), seq_along(x)),
    dbglm_deviance_term(conn), "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
    paste0("SUM(CASE WHEN NOT ", finite[[1L]], " OR w <= 0 OR NOT ", finite[[2L]], " THEN 1 ELSE 0 END) AS bad")
  )
  row = backend_get_query(conn, paste0(cte, "SELECT ", paste(terms, collapse = ", "), " FROM working"))
  H = dbglm_read_information(row, length(x))
  b = vapply(seq_along(x), function(i) as.numeric(row[[paste0("b_", i)]]), numeric(1))
  score = vapply(seq_along(x), function(i) as.numeric(row[[paste0("s_", i)]]), numeric(1))
  list(
    H = H, b = b, score = score, score_max = max(abs(score) / sqrt(diag(H))),
    deviance = row[["deviance"]], bad = row[["bad"]], mu_min = row[["mu_min"]], mu_max = row[["mu_max"]]
  )
}

dbglm_irls = function(conn, table, x, beta, family, control, verbose, separation = "error", cell = FALSE) {
  current = dbglm_evaluate(conn, table, x, beta, family[["link"]], cell = cell)
  initial_information = sqrt(diag(current[["H"]]))
  history = current[["deviance"]]
  halvings_total = 0L
  for (iter in seq_len(control[["glm_max_iter"]])) {
    candidate = dbglm_solve(current[["H"]], current[["b"]])[["beta"]]
    step = 1
    accepted = FALSE
    for (halvings in 0:control[["step_max_halvings"]]) {
      next_beta = beta + step * (candidate - beta)
      next_state = dbglm_evaluate(conn, table, x, next_beta, family[["link"]], cell = cell)
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
    if (isTRUE(verbose)) message(
      "[dbglm] iteration ", iter, ": deviance=",
      signif(current[["deviance"]], 8), ", score=", signif(current[["score_max"]], 3)
    )
    boundary = dbglm_boundary(current, family[["link"]])
    if (step == 1 && relative <= control[["glm_tol"]] && boundary &&
      max(abs(current[["score"]]) / initial_information) <= control[["score_tol"]]) {
      if (separation == "error") stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      warning("Returning boundary probabilities only: coefficients are not a finite MLE; inference is unavailable.", call. = FALSE)
      return(c(current, list(
        beta = beta, iterations = iter, deviance_history = history,
        step_halvings = halvings_total, finite_mle = FALSE, converged = FALSE
      )))
    }
    if (step == 1 && (family[["link"]] == "identity" || relative <= control[["glm_tol"]]) &&
      is.finite(current[["score_max"]]) && current[["score_max"]] <= control[["score_tol"]]) {
      if (dbglm_boundary(current, family[["link"]])) {
        stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      }
      return(c(current, list(
        beta = beta, iterations = iter,
        deviance_history = history, step_halvings = halvings_total
      )))
    }
  }
  stop("GLM did not converge within ", control[["glm_max_iter"]], " iterations.", call. = FALSE)
}

# Fisher-scoring IRLS with absorbed fixed effects.
dbglm_fe_fit = function(conn, base, x, fe, beta, family, control, create, drop, verbose, cache = NULL, ybar = NULL, vcov_type = "iid", cell = FALSE) {
  if (is.null(ybar)) ybar = backend_get_query(conn, paste("SELECT SUM(__dbglm_a*__dbglm_y)/SUM(__dbglm_a) AS y FROM", base))[["y"]]
  intercept = family[["linkfun"]](ybar)
  slopes = attr(fe, "slopes")
  if (is.null(cache)) cache = new.env(parent = emptyenv())
  if (length(fe) > 1L) dbglm_moments_enabled(conn, base, fe, control, cache, length(x) + 1L)
  use_moments = length(fe) > 1L && isTRUE(cache[["moment_enabled"]])
  collect_fe_moments = use_moments && !length(unlist(slopes))
  evaluate_working = function(table, collect = FALSE) {
    bound = if (identical(cache[["moment_fe"]], as.character(fe))) cache[["moment_rows"]] else NULL
    dbglm_fe_evaluate(conn, table, x, fe,
      collect_fe_moments = collect, moment_rows = bound,
      fe_sizes = cache[["fe_sizes"]], limits = cache[["projection_limits"]]
    )
  }
  build_state = dbglm_fe_state_builder(conn, base, x, family[["link"]], create, cell = cell)
  make_working = function(beta_val, c_eta_val, gamma_val = NULL, ops = NULL, info = NULL) {
    build_state(beta_val, c_eta_val, gamma_val, ops, info,
      use_moments = use_moments, collect_fe_moments = collect_fe_moments
    )
  }
  c_eta = intercept
  gamma = NULL
  wk_obj = make_working(beta, c_eta)
  state = wk_obj[["state"]]
  working = wk_obj[["working"]]
  cur_lookup = wk_obj[["lookup"]]
  current = evaluate_working(working, collect_fe_moments)
  initial_halvings = 0L
  while (current[["bad"]] > 0 || dbglm_boundary(current, family[["link"]])) {
    if (initial_halvings >= control[["step_max_halvings"]]) stop("Initial GLM predictor reaches a numerical boundary; choose smaller starts or offsets.", call. = FALSE)
    beta = beta / 2
    drop(working)
    drop(state)
    drop(cur_lookup)
    wk_obj = make_working(beta, c_eta)
    state = wk_obj[["state"]]
    working = wk_obj[["working"]]
    cur_lookup = wk_obj[["lookup"]]
    current = evaluate_working(working, collect_fe_moments)
    initial_halvings = initial_halvings + 1L
  }
  if (collect_fe_moments && is.null(current[["fe_moments"]])) {
    record_projection_attempt(cache, list(
      method = "bounded_moments", status = "declined", reason = "transfer_limit",
      details = list(limits = cache[["projection_limits"]])
    ))
    cache[["moment_enabled"]] = FALSE
    use_moments = FALSE
    collect_fe_moments = FALSE
    drop(working)
    drop(state)
    drop(cur_lookup)
    wk_obj = make_working(beta, c_eta)
    state = wk_obj[["state"]]
    working = wk_obj[["working"]]
    cur_lookup = wk_obj[["lookup"]]
  }
  history = current[["deviance"]]
  total_halvings = initial_halvings
  total_cycles = 0L
  alt_working = NULL
  for (iter in seq_len(control[["glm_max_iter"]])) {
    projected = solve_partial_wls(conn, working, x, fe, control, create, drop,
      cache = cache,
      materialize = FALSE, scales = current[["scales"]], need_table = !use_moments, precomputed = current[["fe_moments"]]
    )
    total_cycles = total_cycles + projected[["cycles"]]
    p = length(x)
    if (p) {
      if (!is.null(projected[["cross_reduction"]])) {
        M_tilde = current[["raw_moments"]] - projected[["cross_reduction"]]
        H_step = M_tilde[seq_len(p), seq_len(p), drop = FALSE]
        b_step = as.numeric(M_tilde[seq_len(p), p + 1L])
        candidate = dbglm_solve(H_step, b_step)[["beta"]]
      } else {
        row = backend_get_query(conn, paste("SELECT", paste(c(dbglm_moment_terms(x, "w"), sprintf("SUM(w*%s*z) AS b_%d", x, seq_along(x))), collapse = ","), "FROM", projected[["table"]]))
        candidate = dbglm_solve(dbglm_read_information(row, p), vapply(seq_len(p), function(i) as.numeric(row[[paste0("b_", i)]]), numeric(1)))[["beta"]]
      }
    } else candidate = numeric()
    if (!is.null(projected[["coefficients"]])) {
      alpha_z = projected[["coefficients"]][, p + 1L]
      gamma_cand = if (p) as.numeric(alpha_z - projected[["coefficients"]][, seq_len(p), drop = FALSE] %*% candidate) else as.numeric(alpha_z)
      if (is.null(gamma)) gamma = numeric(length(gamma_cand))
      step = 1
      accepted = FALSE
      for (halvings in 0:control[["step_max_halvings"]]) {
        next_beta = beta + step * (candidate - beta)
        next_gamma = gamma + step * (gamma_cand - gamma)
        next_c_eta = (1 - step) * c_eta
        next_wk = make_working(next_beta, next_c_eta, next_gamma, projected[["ops"]], projected[["info"]])
        next_state = next_wk[["state"]]
        next_working = next_wk[["working"]]
        next_lookup = next_wk[["lookup"]]
        next_eval = evaluate_working(next_working, halvings == 0L && collect_fe_moments)
        if (dbglm_fe_step_accepted(current, next_eval, family[["link"]], control)) {
          if (halvings > 0L && collect_fe_moments && is.null(next_eval[["fe_moments"]])) {
            next_eval = evaluate_working(next_working, TRUE)
          }
          accepted = TRUE
          break
        }
        drop(next_working)
        drop(next_state)
        drop(next_lookup)
        step = step / 2
      }
      drop(projected[["table"]])
      for (st in projected[["sub_tables"]]) drop(st)
      if (!accepted) dbglm_fe_step_failed(current, next_eval, verbose)
      relative = abs(current[["deviance"]] - next_eval[["deviance"]]) / (0.1 + abs(next_eval[["deviance"]]))
      beta = next_beta
      gamma = next_gamma
      c_eta = next_c_eta
      drop(working)
      drop(state)
      drop(cur_lookup)
      state = next_state
      working = next_working
      cur_lookup = next_lookup
      current = next_eval
    } else {
      if (!is.null(projected[["raw_projected"]])) {
        residual = paste(c("z_tilde", sprintf("-(%s)*%s_tilde", dbglm_number(candidate), x)), collapse = " + ")
        tilde_cols = paste(c("__w", paste0(c(x, "z"), "_tilde")), collapse = ",")
        candidate_state = create(paste0("SELECT * EXCLUDE(", tilde_cols, "), (z - (", residual, ") + __dbglm_o) AS __dbglm_candidate FROM ", projected[["raw_projected"]]), materialize = FALSE)
      } else {
        residual = paste(c("p.z", sprintf("-(%s)*p.%s", dbglm_number(candidate), x)), collapse = " + ")
        candidate_state = create(paste0("SELECT b.*, (b.z - (", residual, ") + b.__dbglm_o) AS __dbglm_candidate FROM ", working, " b JOIN ", projected[["table"]], " p USING (__dbglm_key)"), materialize = FALSE)
      }
      step = 1
      accepted = FALSE
      for (halvings in 0:control[["step_max_halvings"]]) {
        next_state = create(paste0(
          "SELECT * EXCLUDE(eta,mu,d,w,z,score,dev,__dbglm_eta,__dbglm_candidate",
          if (cell) ",score_sq", "), ",
          "__dbglm_eta + ", dbglm_number(step), "*(__dbglm_candidate-__dbglm_eta) AS __dbglm_eta FROM ", candidate_state
        ), materialize = FALSE)
        next_working = dbglm_fe_working(conn, next_state, family[["link"]], create, target = alt_working, cell = cell)
        next_eval = evaluate_working(next_working)
        # Reject Fisher steps that collapse information before the next WLS.
        if (dbglm_fe_step_accepted(current, next_eval, family[["link"]], control)) {
          accepted = TRUE
          break
        }
        if (is.null(alt_working)) drop(next_working)
        drop(next_state)
        step = step / 2
      }
      # working is materialized before releasing its view dependencies.
      if (accepted) {
        drop(next_state)
        next_state = NULL
      }
      drop(candidate_state)
      drop(projected[["table"]])
      for (st in projected[["sub_tables"]]) drop(st)
      if (!accepted) dbglm_fe_step_failed(current, next_eval, verbose)
      relative = abs(current[["deviance"]] - next_eval[["deviance"]]) / (0.1 + abs(next_eval[["deviance"]]))
      beta = beta + step * (candidate - beta)
      if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
        alt_working = working
      } else {
        drop(working)
      }
      drop(state)
      state = next_state
      working = next_working
      current = next_eval
    }
    history = c(history, current[["deviance"]])
    total_halvings = total_halvings + halvings
    if (isTRUE(verbose)) message("[dbglm] iteration ", iter, ": deviance=", signif(current[["deviance"]], 9), ", score=", signif(current[["score_max"]], 3))
    if (step == 1 && (family[["link"]] == "identity" || relative <= control[["glm_tol"]]) && current[["score_max"]] <= control[["score_tol"]]) {
      if (dbglm_boundary(current, family[["link"]])) stop("GLM reached the numerical mean boundary.", call. = FALSE)
      final = solve_partial_wls(conn, working, x, fe, control, create, drop,
        response = FALSE, cache = cache,
        materialize = FALSE, scales = if (p) current[["scales"]][seq_len(p)] else numeric(),
        extra_sums = if (p && !is.null(current[["raw_score"]])) "score" else character(),
        need_table = !use_moments || (p > 0L && (vcov_type != "iid" || isTRUE(cache[["retain_scores"]]))), precomputed = current[["fe_moments"]]
      )
      H = if (p) {
        if (!is.null(final[["cross_reduction"]])) current[["raw_moments"]][seq_len(p), seq_len(p), drop = FALSE] - final[["cross_reduction"]] else
          dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "w"), collapse = ","), "FROM", final[["table"]])), p)
      } else matrix(numeric(), 0, 0)
      total_cycles = total_cycles + final[["cycles"]]
      profiled_score = 0
      if (p) {
        if (any(!is.finite(diag(H))) || any(diag(H) <= 0)) stop("GLM slope information lost identification.", call. = FALSE)
        if (!is.null(final[["extra_reduction"]])) {
          s_tilde = current[["raw_score"]] - final[["extra_reduction"]][, 1L]
        } else {
          score_row = backend_get_query(conn, paste0(
            "SELECT ", paste(sprintf("%s AS s%d", vapply(sprintf("p.%s*b.score", x), function(e) sql_stable_sum(conn, e), character(1)), seq_along(x)), collapse = ","),
            " FROM ", working, " b JOIN ", final[["table"]], " p USING (__dbglm_key)"
          ))
          s_tilde = as.numeric(score_row[1, ])
        }
        profiled_score = max(abs(s_tilde) / sqrt(diag(H)))
      }
      if (profiled_score > control[["score_tol"]]) {
        if (isTRUE(verbose)) message(
          "[dbglm] profiled score=", signif(profiled_score, 6),
          " exceeds tolerance=", control[["score_tol"]], "; continuing IRLS"
        )
        drop(final[["table"]])
        for (st in final[["sub_tables"]]) drop(st)
        next
      }
      current[["score_max"]] = max(current[["score_max"]], profiled_score)
      current[["score_profile"]] = profiled_score
      return(c(current, list(
        beta = beta, H = H, working = working, projected = final[["table"]], raw_projected = final[["raw_projected"]], iterations = iter, deviance_history = history,
        step_halvings = total_halvings, ap_cycles = total_cycles, ap_residual = final[["residual"]], projection_method = if (is.null(final[["method"]])) "cg" else final[["method"]]
      )))
    }
  }
  stop("GLM did not converge within ", control[["glm_max_iter"]], " iterations.", call. = FALSE)
}

# Both coefficient-space and row-space trials enforce the same FE policy.
# Reject steps that collapse information before the next weighted projection.
dbglm_fe_step_accepted = function(current, candidate, link, control) {
  candidate[["score_max"]] <= max(2 * current[["score_max"]], control[["score_tol"]]) &&
    candidate[["bad"]] == 0 && !dbglm_boundary(candidate, link) &&
    is.finite(candidate[["deviance"]]) &&
    candidate[["deviance"]] <= current[["deviance"]] + dbglm_deviance_allowance(current[["deviance"]])
}

dbglm_fe_step_failed = function(current, candidate, verbose) {
  if (isTRUE(verbose)) message(
    "[dbglm] line search exhausted: deviance change=", dbglm_number(candidate[["deviance"]] - current[["deviance"]]),
    "; allowed=", dbglm_number(dbglm_deviance_allowance(current[["deviance"]])),
    "; current score=", dbglm_number(current[["score_max"]]), "; trial score=", dbglm_number(candidate[["score_max"]]),
    "; mu range=", dbglm_number(candidate[["mu_min"]]), ",", dbglm_number(candidate[["mu_max"]])
  )
  stop("GLM step halving failed to find a valid improving step.", call. = FALSE)
}

# Capture only immutable fit context. Strategy flags are passed on each call,
# so a moment-cap fallback cannot leave a stale strategy in this callback.
dbglm_fe_state_builder = function(conn, base, x, link, create, cell = FALSE) {
  backend_info = detect_backend(conn)
  inline_sql = function(sql) {
    if (isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
      create(sql, materialize = TRUE)
    } else {
      paste0("(", backend_sql(conn, sql), ")")
    }
  }
  function(beta_val, c_eta_val, gamma_val = NULL, ops = NULL, info = NULL,
           use_moments = FALSE, collect_fe_moments = FALSE) {
    lookup_sub = NULL
    array_terms = if (collect_fe_moments && isTRUE(backend_info[["is_databricks"]]) &&
      !is.null(gamma_val) && !all(gamma_val == 0)) dbglm_databricks_array_terms(ops, info, gamma_val) else NULL
    if (is.null(gamma_val) || all(gamma_val == 0)) {
      eta = paste(c("__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val), sprintf("(%s)*%s", dbglm_number(beta_val), x)), collapse = " + ")
      st = if (use_moments) inline_sql(paste("SELECT *,", eta, "AS __dbglm_eta FROM", base)) else
        create(paste("SELECT *,", eta, "AS __dbglm_eta FROM", base), materialize = FALSE)
    } else if (!is.null(array_terms)) {
      eta = paste(c(
        "s.__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val),
        sprintf("(%s)*s.%s", dbglm_number(beta_val), x), array_terms
      ), collapse = " + ")
      st = inline_sql(paste0("SELECT s.*,CAST((", eta, ") AS DOUBLE) AS __dbglm_eta FROM ", base, " s"))
    } else {
      frame = dbglm_coefficient_frame(info, gamma_val)
      dictionary = dbglm_coefficient_lookup(conn, frame, create,
        inline = !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))
      )
      lookup = dictionary[["source"]]
      lookup_sub = if (length(dictionary[["tables"]])) dictionary[["tables"]] else NULL
      ctes = dbglm_coefficient_ctes(lookup, length(ops), "eta")
      fe_terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]])
        if (v == "1.0") sprintf("COALESCE(eta_g%d.c1,0.0)", i) else sprintf("(%s)*COALESCE(eta_g%d.c1,0.0)", v, i)
      }, character(1))
      base_terms = c("s.__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val), sprintf("(%s)*s.%s", dbglm_number(beta_val), x), fe_terms)
      joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN eta_g", i, " ON s.", ops[[i]][["f"]], "=eta_g", i, ".__glm_group"), character(1))
      broadcast = if (isTRUE(backend_info[["is_databricks"]])) paste0("/*+ BROADCAST(", paste0("eta_g", seq_along(ops), collapse = ","), ") */ ") else ""
      st = inline_sql(paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast, "s.*, CAST((", paste(base_terms, collapse = " + "), ") AS DOUBLE) AS __dbglm_eta FROM ", base, " s", paste(joins, collapse = "")))
    }
    wk = if (use_moments) inline_sql(paste0(dbglm_working_sql(st, character(), numeric(), link, eta = "__dbglm_eta", cell = cell), "SELECT * FROM working")) else
      dbglm_fe_working(conn, st, link, create, cell = cell)
    list(state = st, working = wk, lookup = lookup_sub)
  }
}

# Working-state construction and diagnostics for absorbed GLMs.
dbglm_fe_working = function(conn, table, link, create, target = NULL, cell = FALSE) {
  create(paste0(dbglm_working_sql(table, character(), numeric(), link, eta = "__dbglm_eta", cell = cell), "SELECT * FROM working"), target = target)
}

dbglm_fe_evaluate = function(conn, working, x, fe, collect_fe_moments = FALSE, moment_rows = NULL, fe_sizes = NULL, limits = projection_limits()) {
  backend_info = detect_backend(conn)
  finite = sql_is_finite(conn, c("w", "eta"), backend_info)
  cols = c(x, "z")
  terms = c(
    dbglm_deviance_term(conn), "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
    paste0("SUM(CASE WHEN NOT ", finite[[1L]], " OR w<=0 OR NOT ", finite[[2L]], " THEN 1 ELSE 0 END) AS bad"),
    dbglm_moment_terms(cols, "w")
  )
  if (length(x)) {
    x_score = vapply(paste0(x, "*score"), function(e) sql_stable_sum(conn, e, backend_info), character(1))
    terms = c(
      terms, sprintf("%s AS raw_s%d", x_score, seq_along(x)),
      sprintf("ABS(%s)/SQRT(SUM(w*%s*%s)) AS s%d", x_score, x, x, seq_along(x))
    )
  }
  if (collect_fe_moments && length(fe) >= 2L && !length(unlist(attr(fe, "slopes")))) {
    max_edges = limits[["summary_rows"]]
    max_entries = limits[["design_entries"]]
    # Six wire fields include the marginal/global payload slots. Reserve their
    # decoded numbers, the global row and an overflow probe before fetching.
    payload = if (is.null(fe_sizes)) max_entries * (length(cols) + 1L) / max(1L, length(cols)) else sum(fe_sizes) * (length(cols) + 1L)
    row_cap = min(max_edges, floor((limits[["transfer_cells"]] - payload - length(terms)) / 6L) - 2L)
    scalar = function() {
      dbglm_fe_evaluate(conn, working, x, fe, collect_fe_moments = FALSE)
    }
    if (row_cap < 1L || (length(moment_rows) == 1L && is.finite(moment_rows) && moment_rows > row_cap)) return(scalar())
    validate_moment_limits(max_edges, max_entries)
    K = length(fe)
    # The initial projector bounds support on the owned immutable base table.
    # A working predictor preserves every key, so reuse that verified bound.
    # Without it, fetch only the cap plus the global row and an overflow row.
    known_bound = length(moment_rows) == 1L && is.finite(moment_rows) &&
      moment_rows >= 0 && moment_rows == floor(moment_rows) && moment_rows <= row_cap
    res = read_fe_moments(conn, fe_moments_sql(conn, working, fe,
      fe_value_terms(conn, cols, "score"), terms,
      limit = if (known_bound) NULL else row_cap + 2, compact = TRUE
    ),
    max_rows = if (known_bound) moment_rows + 1 else row_cap + 1
    )
    if (is.null(res)) {
      if (known_bound) stop("The fixed-effect moment support changed within fitting.", call. = FALSE)
      return(scalar())
    }
    if (known_bound && nrow(res) != moment_rows + 1) {
      stop("The fixed-effect moment support changed within fitting.", call. = FALSE)
    }
    sets = as.integer(res[["__set"]])
    row = attr(res, "dbglm_globals")
    if (is.null(row)) row = res[sets == 0L, , drop = FALSE]
    fe_moments = NULL
    if (!is.na(row[["bad"]]) && row[["bad"]] == 0 && sum(sets != 0L) <= max_edges) {
      unpacked = unpack_fe_moments(res, K, length(cols), 1L)
      if (!is.null(unpacked) && unpacked[["n"]] * max(1L, length(cols)) <= max_entries) fe_moments = unpacked
    }
    singles = sets >= 1L & sets <= K
    h = as.numeric(res[["__h"]][singles])
    compact = attr(res, "fe_values", exact = TRUE)
    score = if (is.null(compact)) res[["__ex1"]][singles] else
      compact[["values"]][match(which(singles), compact[["rows"]]), length(cols) + 1L]
    # Also certify the Newton-scale group equation. A square-root information
    # score alone vanishes with a group's weight even when its predictor is
    # far from its finite solution. This uses the already collected moments.
    s_val = ifelse(is.finite(h) & h > 0 & is.finite(score),
      abs(score) / sqrt(h) / pmin(1, sqrt(h)), 0
    )
    result = dbglm_fe_diagnostics(row, x, if (length(s_val)) max(s_val) else 0)
    # Retain the explicit NULL when a bounded summary cannot be used.
    result["fe_moments"] = list(fe_moments)
    return(result)
  }
  fe_queries = character()
  for (i in seq_along(fe)) for (v in c("1.0", attr(fe, "slopes")[[i]])) {
    fe_queries = c(fe_queries, paste0(
      "SELECT COALESCE(ABS(", sql_stable_sum(conn, paste0("score*(", v, ")")),
      ")/NULLIF(SQRT(SUM(w*(", v, ")*(", v, ")))*LEAST(1.0,SQRT(SUM(w))),0),0) AS s FROM ", working, " GROUP BY ", fe[i]
    ))
  }
  if (length(fe_queries)) {
    terms = c(terms, paste0("(SELECT MAX(s) FROM (", paste(fe_queries, collapse = " UNION ALL "), ") t) AS fe_score_max"))
  }
  row = backend_get_query(conn, paste("SELECT", paste(terms, collapse = ","), "FROM", working))
  fe_score = if (length(fe_queries) && !is.na(row[["fe_score_max"]])) as.numeric(row[["fe_score_max"]]) else 0
  dbglm_fe_diagnostics(row, x, fe_score)
}

# Decode the same scalar diagnostics from fused moments and SQL projections.
dbglm_fe_diagnostics = function(row, x, fe_score) {
  score = if (length(x)) max(as.numeric(row[1, paste0("s", seq_along(x)), drop = FALSE])) else 0
  raw_score = if (length(x)) as.numeric(row[1, paste0("raw_s", seq_along(x)), drop = FALSE]) else numeric()
  raw_moments = dbglm_read_information(row, length(x) + 1L)
  list(
    deviance = row[["deviance"]], mu_min = row[["mu_min"]], mu_max = row[["mu_max"]], bad = row[["bad"]],
    score_max = max(score, fe_score), raw_score = raw_score, raw_moments = raw_moments,
    scales = sqrt(pmax(diag(raw_moments), .Machine[["double.xmin"]]))
  )
}

# Databricks inlines CTEs, so a UNION branch per FE pair repeatedly plans the
# predictor. Keep three aggregation branches: globals, FE marginals and pairs.
# Pairs need only their weight sum, not every diagnostic/score aggregate.
dbglm_databricks_moments_sql = function(source, fe, value_terms, global_terms = character(), limit = NULL, text_keys = FALSE) {
  K = length(fe)
  # GROUPING_ID uses a signed BIGINT. Keep masks exactly representable in R.
  if (K < 2L || K > 30L) return(NULL)
  aliases = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  id_type = if (text_keys) "STRING" else "BIGINT"
  value_names = aliases(value_terms)
  global_names = aliases(global_terms)
  nulls = function(names) if (length(names)) paste0("CAST(NULL AS DOUBLE) AS ", names) else character()
  grouped = function(groups, sets, terms) {
    masks = vapply(groups, function(g) sum(2^(K - setdiff(seq_len(K), g))), numeric(1))
    key = function(j) {
      cases = vapply(seq_along(groups), function(i) {
        g = groups[[i]]
        if (length(g) < j) return("")
        sprintf("WHEN %.0f THEN CAST(%s AS %s)", masks[[i]], fe[g[[j]]], id_type)
      }, character(1))
      if (!any(nzchar(cases))) return(paste0("CAST(NULL AS ", id_type, ")"))
      paste("CASE __gid", paste(cases, collapse = " "), "END")
    }
    fields = c(
      paste("CASE __gid", paste(sprintf("WHEN %.0f THEN %d", masks, sets), collapse = " "), "END AS __set"),
      paste(key(1L), "AS __i"), paste(key(2L), "AS __j"),
      nulls(global_names), aliases(terms), nulls(setdiff(value_names, aliases(terms)))
    )
    # Keep the UNION schema in the same order even when pairs omit values.
    paste0(
      "SELECT ", paste(fields, collapse = ","), " FROM (SELECT GROUPING_ID(", paste(fe, collapse = ","),
      ") AS __gid,", paste(c(fe, terms), collapse = ","), " FROM ", source, " GROUP BY GROUPING SETS (",
      paste(vapply(groups, function(g) paste0("(", paste(fe[g], collapse = ","), ")"), character(1)), collapse = ","), ")) grouped"
    )
  }
  pairs = utils::combn(K, 2L, simplify = FALSE)
  branches = c(
    if (length(global_terms)) paste0(
      "SELECT 0 AS __set,CAST(NULL AS ", id_type, ") AS __i,CAST(NULL AS ", id_type, ") AS __j,",
      paste(c(global_terms, nulls(value_names)), collapse = ","), " FROM ", source
    ),
    grouped(as.list(seq_len(K)), seq_len(K), value_terms),
    grouped(pairs, K + seq_along(pairs), value_terms[1L])
  )
  sql = paste(branches, collapse = " UNION ALL ")
  # Keep stable set ordering for bounded preflight queries. The immutable FE
  # support lets IRLS reuse this bound without sorting every evaluation.
  if (!is.null(limit)) sql = paste0("SELECT * FROM (", sql, ") limited ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  if (!length(global_terms)) return(sql)
  # Only the global row has these scalars. Avoid sending NULL copies on every
  # pair row; the reader keeps globals separate from the numeric FE moments.
  fields = paste(sprintf("'%s',%s", global_names, global_names), collapse = ",")
  paste0(
    "SELECT __set,__i,__j,", paste(value_names, collapse = ","),
    ",CASE WHEN __set=0 THEN to_json(named_struct(", fields, ")) END AS __globals FROM (", sql, ") packed"
  )
}

# Small, dense FE rank dictionaries can use constant array indexing instead
# of coefficient joins. Large/sparse dictionaries keep the bounded join path.
# These are FE-level values only; observation data never enter SQL literals.
dbglm_databricks_array_terms = function(
  ops, info, coefficients,
  max_entries = getOption("dbreg.glm_databricks_array_max_entries", 65536L), max_bytes = 512 * 1024
) {
  if (!is.numeric(max_entries) || length(max_entries) != 1L || !is.finite(max_entries) ||
    max_entries < 0 || max_entries != floor(max_entries)) {
    stop("GLM coefficient array limit must be a nonnegative whole number.", call. = FALSE)
  }
  if (max_entries == 0 || !length(ops)) return(NULL)
  ids = lapply(info, function(inf) as.numeric(inf[["id"]]))
  sizes = lengths(ids)
  if (length(ids) != length(ops) || any(!sizes) || length(coefficients) != sum(sizes)) {
    stop("Invalid FE coefficient dictionary.", call. = FALSE)
  }
  if (any(!is.finite(coefficients))) stop("Nonfinite projection coefficients.", call. = FALSE)
  if (any(vapply(ids, function(id) any(!is.finite(id) | id < 1 | id != floor(id)) || anyDuplicated(id) > 0L, logical(1)))) return(NULL)
  maxima = vapply(ids, max, numeric(1))
  if (any(maxima > 4 * sizes) || sum(maxima) > max_entries || any(maxima > .Machine[["integer.max"]])) return(NULL)
  starts = cumsum(c(0L, sizes))
  terms = character(length(ops))
  bytes = 0
  for (i in seq_along(ops)) {
    values = numeric(maxima[[i]])
    values[ids[[i]]] = coefficients[starts[[i]] + seq_len(sizes[[i]])]
    payload = paste0("[", paste(dbglm_number(values), collapse = ","), "]")
    encoded = raw_to_base64(charToRaw(enc2utf8(payload)))
    effect = paste0(
      "COALESCE(element_at(from_json(CAST(unbase64('", encoded,
      "') AS STRING),'ARRAY<DOUBLE>'),CAST(s.", ops[[i]][["f"]], " AS INT)),0.0)"
    )
    v = ops[[i]][["v"]]
    terms[[i]] = if (v == "1.0") effect else paste0("s.", v, "*(", effect, ")")
    bytes = bytes + nchar(terms[[i]], type = "bytes")
    if (bytes > max_bytes) return(NULL)
  }
  terms
}
