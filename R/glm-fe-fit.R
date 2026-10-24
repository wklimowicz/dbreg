# Fisher-scoring IRLS with absorbed fixed effects.
dbglm_fe_fit = function(conn, base, x, fe, beta, family, control, create, drop, verbose, cache = NULL, ybar = NULL, vcov_type = "iid") {
  if (is.null(ybar)) ybar = dbglm_get_query(conn, paste("SELECT SUM(__dbglm_a*__dbglm_y)/SUM(__dbglm_a) AS y FROM", base))[["y"]]
  intercept = family[["linkfun"]](ybar)
  slopes = attr(fe, "slopes")
  if (is.null(cache)) cache = new.env(parent = emptyenv())
  if (length(fe) > 1L) dbglm_moments_enabled(conn, base, fe, control, cache)
  use_moments = length(fe) > 1L && isTRUE(cache[["moment_enabled"]])
  collect_fe_moments = use_moments && !length(unlist(slopes))
  evaluate_working = function(table, collect = FALSE) {
    bound = if (identical(cache[["moment_fe"]], as.character(fe))) cache[["moment_rows"]] else NULL
    dbglm_fe_evaluate(conn, table, x, fe, collect_fe_moments = collect, moment_rows = bound)
  }
  build_state = dbglm_fe_state_builder(conn, base, x, family[["link"]], create)
  make_working = function(beta_val, c_eta_val, gamma_val = NULL, ops = NULL, info = NULL) {
    build_state(beta_val, c_eta_val, gamma_val, ops, info,
      use_moments = use_moments, collect_fe_moments = collect_fe_moments)
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
    projected = solve_partial_wls(conn, working, x, fe, control, create, drop, cache = cache,
      materialize = FALSE, scales = current[["scales"]], need_table = !use_moments, precomputed = current[["fe_moments"]])
    total_cycles = total_cycles + projected[["cycles"]]
    p = length(x)
    if (p) {
      if (!is.null(projected[["cross_reduction"]])) {
        M_tilde = current[["raw_moments"]] - projected[["cross_reduction"]]
        H_step = M_tilde[seq_len(p), seq_len(p), drop = FALSE]
        b_step = as.numeric(M_tilde[seq_len(p), p + 1L])
        candidate = dbglm_solve(H_step, b_step)[["beta"]]
      } else {
        row = dbglm_get_query(conn, paste("SELECT", paste(c(dbglm_moment_terms(x, "w"), sprintf("SUM(w*%s*z) AS b_%d", x, seq_along(x))), collapse = ","), "FROM", projected[["table"]]))
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
        next_state = create(paste0("SELECT * EXCLUDE(eta,mu,d,w,z,score,dev,__dbglm_eta,__dbglm_candidate), ",
          "__dbglm_eta + ", dbglm_number(step), "*(__dbglm_candidate-__dbglm_eta) AS __dbglm_eta FROM ", candidate_state), materialize = FALSE)
        next_working = dbglm_fe_working(conn, next_state, family[["link"]], create, target = alt_working)
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
    if (step == 1 && relative <= control[["glm_tol"]] && current[["score_max"]] <= control[["score_tol"]]) {
      if (dbglm_boundary(current, family[["link"]])) stop("GLM reached the numerical mean boundary.", call. = FALSE)
      final = solve_partial_wls(conn, working, x, fe, control, create, drop, response = FALSE, cache = cache,
        materialize = FALSE, scales = if (p) current[["scales"]][seq_len(p)] else numeric(),
        extra_sums = if (p && !is.null(current[["raw_score"]])) "score" else character(),
        need_table = !use_moments || (p > 0L && vcov_type != "iid"), precomputed = current[["fe_moments"]])
      H = if (p) {
        if (!is.null(final[["cross_reduction"]])) current[["raw_moments"]][seq_len(p), seq_len(p), drop = FALSE] - final[["cross_reduction"]] else
          dbglm_read_information(dbglm_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "w"), collapse = ","), "FROM", final[["table"]])), p)
      } else matrix(numeric(), 0, 0)
      total_cycles = total_cycles + final[["cycles"]]
      profiled_score = 0
      if (p) {
        if (any(!is.finite(diag(H))) || any(diag(H) <= 0)) stop("GLM slope information lost identification.", call. = FALSE)
        if (!is.null(final[["extra_reduction"]])) {
          s_tilde = current[["raw_score"]] - final[["extra_reduction"]][, 1L]
        } else {
          score_row = dbglm_get_query(conn, paste0("SELECT ", paste(sprintf("%s AS s%d", vapply(sprintf("p.%s*b.score", x), function(e) dbglm_exact_sum(conn, e), character(1)), seq_along(x)), collapse = ","),
            " FROM ", working, " b JOIN ", final[["table"]], " p USING (__dbglm_key)"))
          s_tilde = as.numeric(score_row[1, ])
        }
        profiled_score = max(abs(s_tilde) / sqrt(diag(H)))
      }
      if (profiled_score > control[["score_tol"]]) {
        if (isTRUE(verbose)) message("[dbglm] profiled score=", signif(profiled_score, 6),
          " exceeds tolerance=", control[["score_tol"]], "; continuing IRLS")
        drop(final[["table"]])
        for (st in final[["sub_tables"]]) drop(st)
        next
      }
      current[["score_max"]] = max(current[["score_max"]], profiled_score)
      current[["score_profile"]] = profiled_score
      return(c(current, list(beta = beta, H = H, working = working, projected = final[["table"]], raw_projected = final[["raw_projected"]], iterations = iter, deviance_history = history,
        step_halvings = total_halvings, ap_cycles = total_cycles, ap_residual = final[["residual"]], projection_method = if (is.null(final[["method"]])) "cg" else final[["method"]])))
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
    "; mu range=", dbglm_number(candidate[["mu_min"]]), ",", dbglm_number(candidate[["mu_max"]]))
  stop("GLM step halving failed to find a valid improving step.", call. = FALSE)
}
