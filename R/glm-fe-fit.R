# Encode a bounded coefficient state against dense FE ranks. The base is
# immutable across all trials and keeps the native key grouping in SQL.
dbglm_fe_predictor_source = function(conn, base, x, fe, beta, effects, ops = NULL,
                                       offset = TRUE, column = "c1", output = "__dbglm_eta", manager = NULL) {
  limits = projection_limits()
  budget = list2env(list(upload_bytes = 0), parent = emptyenv())
  # Every query of an IRLS step re-reads this source; on local engines keep
  # large effect literals in owned tables (released with the fit) instead.
  store = !is.null(manager) && identical(detect_backend(conn)[["name"]], "duckdb")
  keys = if (is.null(ops)) fe else vapply(ops, `[[`, character(1), "f")
  axes = if (is.null(ops)) rep("1.0", length(fe)) else vapply(ops, `[[`, character(1), "v")
  aliases = paste0("a", seq_along(keys))
  joins = vapply(seq_along(keys), function(i) {
    lookup = if (is.character(effects)) {
      paste("SELECT __cg_group AS id,", column, "AS c1 FROM", effects, "WHERE __cg_op =", i)
    } else {
      if (store && length(effects[[i]]) > 1024L) {
        paste("SELECT * FROM", fe_duckdb_lookup_table(conn, matrix(effects[[i]], ncol = 1L), manager))
      } else fe_coefficient_sql(conn, matrix(effects[[i]], ncol = 1L), budget, limits)
    }
    paste0(if (offset) " JOIN (" else " LEFT JOIN (", lookup, ") ", aliases[i], " ON b.", keys[i],
      "=", aliases[i], ".id")
  }, character(1))
  terms = c(if (offset) "b.__dbglm_o",
    paste0(if (offset) paste0(aliases, ".c1") else paste0("COALESCE(", aliases, ".c1,0.0)"),
      "*", ifelse(axes == "1.0", axes, paste0("b.", axes))),
    if (length(x)) paste0("(", dbglm_number(beta), ")*b.", x))
  eta = paste(terms, collapse = "+")
  paste0("(SELECT b.*, (", eta, ") AS ", output, " FROM (SELECT * FROM ", base,
    ") b", paste(joins, collapse = ""), ") glm_trial")
}

dbglm_fe_trial = function(conn, base, x, fe, beta, effects, link, newton = FALSE,
                          cell = FALSE, ops = NULL, manager = NULL, materialize = FALSE) {
  source = dbglm_fe_predictor_source(conn, base, x, fe, beta, effects, ops, manager = manager)
  state = dbglm_evaluate(conn, source, x, beta, link, eta = "__dbglm_eta", cell = cell,
    materialize = materialize)
  sql = dbglm_working_sql(source, x, beta, link, eta = "__dbglm_eta", conn = conn,
    cell = cell, materialize = materialize)
  state[["working"]] = paste0("(", sql, "SELECT * FROM working) glm_fe_working")
  state[["working_columns"]] = x
  state[["information_type"]] = "fisher"
  if (newton) {
    observed = dbglm_working_sql(source, x, beta, link, eta = "__dbglm_eta",
      conn = conn, newton = TRUE, cell = cell, materialize = materialize)
    state[["observed"]] = paste0("(", observed, "SELECT * FROM working) glm_fe_observed")
  }
  state
}

# One group's score equation, as an aggregate over that group's working rows.
dbglm_fe_group_score_sql = function(conn) {
  sum_weight = sql_stable_sum(conn, "w")
  guard = paste0("(CASE WHEN ", sum_weight, "<1.0 THEN SQRT(", sum_weight, ") ELSE 1.0 END)")
  paste0("ABS(", sql_stable_sum(conn, "score"), ")/(SQRT(", sum_weight, ")*", guard, ")")
}

dbglm_fe_group_score = function(conn, working, fe, nobs, ops = NULL) {
  sum_weight = sql_stable_sum(conn, "w")
  # A small group's equation must also pass on the Newton scale. Use base
  # group information, independently of a varying coordinate's measurement units.
  guard = paste0("(CASE WHEN ", sum_weight, "<1.0 THEN SQRT(", sum_weight, ") ELSE 1.0 END)")
  if (!is.null(ops)) {
    scores = vapply(ops, function(op) {
      v = op[["v"]]
      numerator = sql_stable_sum(conn, paste0("score*", v))
      denominator = sql_stable_sum(conn, paste0("w*", v, "*", v))
      paste0("SELECT COALESCE(ABS(", numerator, ")/NULLIF(SQRT(", denominator,
        ")*", guard, ",0),0) AS group_score,COUNT(*) AS n FROM ", working, " GROUP BY ", op[["f"]])
    }, character(1))
    row = DBI::dbGetQuery(conn, paste0("SELECT MAX(group_score) AS score,SUM(n) AS n FROM (",
      paste(scores, collapse = " UNION ALL "), ") varying_scores"))
    if (!isTRUE(as.numeric(row[["n"]]) == nobs * length(ops)) || !is.finite(row[["score"]])) {
      projection_abort("certificate_failed", "GLM varying score certificate is invalid.")
    }
    return(as.numeric(row[["score"]]))
  }
  groups = paste0("(", fe, ")", collapse = ",")
  row = DBI::dbGetQuery(conn, paste0(
    "SELECT MAX(group_score) AS score, SUM(n) AS n FROM (SELECT ", dbglm_fe_group_score_sql(conn),
    " AS group_score, COUNT(*) AS n FROM ", working,
    " GROUP BY GROUPING SETS (", groups, ")) glm_fe_scores"))
  if (!isTRUE(as.numeric(row[["n"]]) == nobs * length(fe)) ||
      !is.finite(row[["score"]])) {
    projection_abort("certificate_failed", "GLM fixed-effect score certificate is invalid.")
  }
  as.numeric(row[["score"]])
}

dbglm_identified_slopes = function(H, raw, tol) {
  keep = which(diag(H) > tol * pmax(diag(raw), .Machine[["double.xmin"]]))
  if (length(keep)) keep = keep[dbglm_basis(H[keep, keep, drop = FALSE], tol)]
  keep
}

dbglm_effect_delta = function(direction, correction, manager, column = "delta", inline = FALSE) {
  if (!is.null(manager)) return(dbglm_fe_sql_delta(direction, correction, manager, column, inline))
  lapply(direction[["projection"]], function(piece) {
    as.numeric(piece[, 1L] - piece[, -1L, drop = FALSE] %*% correction)
  })
}

dbglm_fe_irls = function(conn, base, x, fe, sizes, nobs, beta, effects,
                         family, control, verbose = FALSE, cluster = FALSE, identify_only = FALSE, cell = FALSE, manager = NULL, ops = NULL, score_keys = NULL, initial_intercept = 0, lookup_manager = NULL, step_scale = rep(1, length(x)), initial = NULL) {
  link = family[["link"]]
  tol = dbglm_projection_tolerance(control)
  newton = dbglm_newton_link(link) && identical(control[["irls"]], "newton")
  sql_state = !is.null(manager)
  duckdb_bounded = control[["projection_method"]] == "duckdb_bounded" &&
    inherits(conn, "duckdb_connection") && link == "cloglog" && length(fe) > 1L &&
    is.null(ops)
  project = if (sql_state) {
    cg_control = dbglm_projection_cg_control(control)
    function(working, carry_cluster, keys = NULL, raw_state = NULL, warm = NULL) {
      dbglm_fe_sql_step(conn, working, x, fe, nobs, tol, manager, carry_cluster, cell, cg_control, ops, keys, warm)
    }
  } else {
    function(working, carry_cluster, keys = NULL, raw_state = NULL, warm = NULL) {
      dbglm_fe_step(conn, working, x, fe, sizes, nobs, tol, carry_cluster, cell, keys,
        max_iter = control[["fe_max_iter"]], lookup_manager = lookup_manager,
        raw_state = raw_state, certify_weak = duckdb_bounded)
    }
  }
  drop_projection = function(projection) {
    if (sql_state) dbglm_fe_sql_release(manager, projection)
    for (table in rev(projection[["lookup_tables"]])) lookup_manager[["drop"]](table)
  }
  initial_halvings = 0L
  # An identification pass that kept every column already evaluated and
  # projected this exact start; reuse it instead of repeating the solve.
  reuse = !is.null(initial) && identical(initial[["beta"]], beta) && identical(initial[["x"]], x)
  if (!is.null(initial) && !reuse) drop_projection(initial[["projected"]])
  if (reuse) {
    current = initial[["current"]]
    projected = initial[["projected"]]
  }
  if (!reuse) repeat {
    current = dbglm_fe_trial(conn, base, x, fe, beta, effects, link, newton, cell,
      ops, lookup_manager, materialize = duckdb_bounded)
    projected = NULL
    if (dbglm_valid_start(current, link)) {
      projected = tryCatch(project(current[["working"]], cluster, score_keys, current),
        dbreg_projection_numerical = function(e) {
          if (!e[["reason"]] %in% c("singular", "certificate_failed")) stop(e)
          e
        })
      if (!inherits(projected, "error")) break
    }
    if (isTRUE(verbose)) message("[dbglm] initial predictor recovery ", initial_halvings + 1L,
      if (inherits(projected, "error")) paste0(": ", projected[["reason"]]) else ": numerical boundary")
    if (initial_halvings >= control[["step_max_halvings"]]) {
      if (inherits(projected, "error")) stop(projected)
      stop("Initial GLM predictor reaches a numerical boundary; choose smaller starts or offsets.", call. = FALSE)
    }
    next_beta = beta / 2
    if (sql_state) {
      damped = dbglm_fe_sql_damp(effects, initial_intercept, manager)
      manager[["drop"]](effects)
      effects = damped
    } else {
      next_effects = lapply(seq_along(effects), function(i) effects[[i]] / 2 +
        if (i == 1L) initial_intercept / 2 else 0)
      # Damping cannot repair a failed projection at an unchanged anchor.
      if (identical(next_beta, beta) && identical(next_effects, effects)) {
        if (inherits(projected, "error")) stop(projected)
        stop("Initial GLM predictor reaches a numerical boundary; choose smaller starts or offsets.", call. = FALSE)
      }
      effects = next_effects
    }
    beta = next_beta
    initial_halvings = initial_halvings + 1L
  }
  if (isTRUE(verbose) && initial_halvings) message("[dbglm] initial predictor recovered after ",
    initial_halvings, " halvings")
  # The initial projected information selects a deterministic structural basis.
  keep = dbglm_identified_slopes(projected[["H"]], current[["H"]], control[["collin_tol"]])
  if (identify_only && length(keep) == length(x)) {
    return(list(restart = keep, effects = effects,
      initial = list(current = current, projected = projected, beta = beta, x = x)))
  }
  if (identify_only || length(keep) != length(x)) {
    drop_projection(projected)
    return(list(restart = keep, effects = effects))
  }
  history = current[["deviance"]]
  halvings_total = initial_halvings
  score_max = dbglm_score_norm(current[["score"]], projected[["H"]]) / dbglm_state_unit(current)
  # A SQL projection of these working rows may already carry the group scores.
  group_score = function(projection) {
    value = projection[["fe_group_score"]]
    if (is.null(value)) value = dbglm_fe_group_score(conn, current[["working"]], fe, nobs, ops)
    value / dbglm_state_unit(current)
  }
  fe_score = group_score(projected)
  finish = function(converged, reason) {
    list(beta = beta, effects = effects, H = projected[["H"]],
      working = current[["working"]], projected = projected[["residual_source"]],
      deviance = current[["deviance"]],
      deviance_history = history, score_max = score_max, fe_score_max = fe_score,
      converged = converged, stop_reason = reason, iterations = iter,
      step_halvings = halvings_total, irls = if (newton) "newton" else "fisher",
      projection_tolerance = tol, projection_max_iter = control[["fe_max_iter"]],
      projection_method = if (sql_state) projected[["method"]] else "bounded_moments")
  }
  initial_factor = if (length(x)) dbglm_info_factor(projected[["H"]]) else list(R = NULL, scale = numeric())
  is_binomial_link = link %in% c("logit", "probit", "cloglog")
  needs_step_guard = is_binomial_link || link == "log"
  prev_step_inf = Inf
  for (iter in seq_len(control[["glm_max_iter"]])) {
    previous_score = max(score_max, fe_score)
    direction = if (newton) {
      project(current[["observed"]], FALSE)
    } else {
      projected
    }
    correction = if (length(x)) {
      dbglm_solve(direction[["H"]], direction[["score"]])[["beta"]]
    } else {
      numeric()
    }
    # SQL deltas stay inline subqueries of the projection coefficients, which
    # remain owned until the accepted step's next projection replaces them.
    delta = dbglm_effect_delta(direction, correction, manager, inline = TRUE)
    # Only the fitted predictor is invariant to FE normalizations. Include its
    # nuisance correction even when no structural columns remain; small group
    # scores alone cannot distinguish an escaping FE coefficient from an MLE.
    # Use Fisher coordinates for this guard even when Newton supplies the step:
    # clamped noncanonical tails can make the observed-information step tiny.
    guard_correction = if (newton && length(x)) dbglm_solve(projected[["H"]],
      projected[["score"]])[["beta"]] else correction
    guard_delta = if (newton) dbglm_effect_delta(projected, guard_correction, manager, inline = TRUE) else delta
    direction_source = if (needs_step_guard) dbglm_fe_predictor_source(conn, base,
      x, fe, guard_correction, guard_delta, ops, offset = FALSE, column = "delta", manager = lookup_manager) else NULL
    step_inf = if (needs_step_guard) as.numeric(DBI::dbGetQuery(conn,
      paste("SELECT MAX(ABS(__dbglm_eta)) AS maximum FROM", direction_source))[["maximum"]]) else 0
    direction_boundary = needs_step_guard && is.finite(step_inf) && step_inf > 1e-4 &&
      step_inf >= 0.5 * prev_step_inf && max(score_max, fe_score) <= control[["score_tol"]] &&
      dbglm_direction_boundary(conn, direction_source, link, "__dbglm_eta")
    prev_collapsed = dbglm_info_collapsed(projected[["H"]], initial_factor[["R"]], initial_factor[["scale"]])
    prev_stalled = needs_step_guard && step_inf > 1e-4 && step_inf >= 0.5 * prev_step_inf
    step = 1
    accepted = FALSE
    for (halvings in 0:control[["step_max_halvings"]]) {
      next_beta = beta + step * correction
      next_effects = if (sql_state) dbglm_fe_sql_update(effects, delta, step, manager) else
        Map(function(a, b) a + step * b, effects, delta)
      next_state = dbglm_fe_trial(conn, base, x, fe, next_beta, next_effects, link, newton, cell,
        ops, lookup_manager, materialize = duckdb_bounded)
      if (isTRUE(next_state[["bad"]] == 0) && is.finite(next_state[["deviance"]]) &&
          next_state[["deviance"]] <= current[["deviance"]] +
            dbglm_deviance_allowance(current[["deviance"]])) {
        accepted = TRUE
        break
      }
      if (sql_state) manager[["drop"]](next_effects)
      step = step / 2
    }
    # The inline delta reads the Newton direction's coefficients until here.
    if (newton) drop_projection(direction)
    halvings_total = halvings_total + halvings
    if (!accepted && !control[["tolerate_max_iter"]]) {
      varying = any(vapply(ops, function(op) !identical(op[["v"]], "1.0"), logical(1)))
      stop("GLM step halving failed to find a valid improving step.", if (varying) paste(
        " A varying-slope group may be separated (its outcomes fit exactly along its slope);",
        "fixef_rm does not detect separation in slope groups, so remove such groups from the source."),
        call. = FALSE)
    }
    if (accepted) {
      relative = abs(current[["deviance"]] - next_state[["deviance"]]) /
        (0.1 + abs(next_state[["deviance"]]))
      beta = next_beta
      if (sql_state) {
        manager[["drop"]](effects)
      }
      effects = next_effects
      current = next_state
      history = c(history, current[["deviance"]])
      # The previous projection's coefficients warm-start the next solve;
      # drop its observation-sized residual table before creating the next one.
      previous = projected
      if (sql_state && !is.null(previous[["residual_source"]])) {
        manager[["drop"]](previous[["residual_source"]])
        previous[["residual_source"]] = NULL
      }
      projected = tryCatch(project(current[["working"]], cluster, score_keys, current, warm = previous),
        dbreg_projection_numerical = function(e) {
          if (!e[["reason"]] %in% c("singular", "certificate_failed")) stop(e)
          e
        })
      drop_projection(previous)
      if (inherits(projected, "error")) {
        eps = .Machine[["double.eps"]]
        if (needs_step_guard && (dbglm_boundary(current, link, collapsed = prev_collapsed, stalled = prev_stalled) ||
            (is_binomial_link && isTRUE(current[["mu_min"]] <= eps || current[["mu_max"]] >= 1 - eps)))) {
          stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
        }
        stop(projected)
      }
    }
    score_max = dbglm_score_norm(current[["score"]], projected[["H"]]) / dbglm_state_unit(current)
    fe_score = group_score(projected)
    if (verbose) message("[dbglm] iteration ", iter, ", score ",
      signif(score_max, 4), ", FE score ", signif(fe_score, 4))
    if (!accepted) return(finish(FALSE, "step_halving_failed"))
    eig_ratio = dbglm_min_info_ratio(projected[["H"]], initial_factor[["R"]], initial_factor[["scale"]])
    collapsed = isTRUE(eig_ratio <= 1e-8)
    deep_collapsed = isTRUE(eig_ratio <= 1e-11)
    stalled = needs_step_guard && step == 1 && step_inf > 1e-4 && step_inf >= 0.5 * prev_step_inf
    boundary = dbglm_boundary(current, link, collapsed = collapsed, stalled = stalled) || direction_boundary
    if (step == 1 && relative <= control[["glm_tol"]] && boundary &&
      (max(score_max, fe_score) <= control[["score_tol"]] || deep_collapsed)) {
      stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
    }
    step_small = !needs_step_guard || (is.finite(step_inf) && step_inf <= 1e-4)
    if (step_small && dbglm_converged_step(step, relative, link, max(score_max, fe_score), previous_score, control)) {
      if (boundary) {
        stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      }
      return(finish(TRUE, "converged"))
    }
    prev_step_inf = if (step == 1) step_inf else Inf
  }
  if (control[["tolerate_max_iter"]]) return(finish(FALSE, "max_iter"))
  dbglm_nonconvergence(current, link, control[["glm_max_iter"]])
}

# Intermediate fits may relax projection accuracy. Final convergence still
# depends on independently recomputed structural and fixed-effect scores.
dbglm_projection_tolerance = function(control, strict = FALSE) {
  tolerance = control[["fe_tol"]]
  if (strict || control[["tighten_fe_tol"]]) {
    tolerance = min(tolerance, control[["score_tol"]] / 100)
  }
  tolerance
}

dbglm_projection_cg_control = function(control) {
  cg = fe_cg_control()
  cg[["fe_max_iter"]] = control[["fe_max_iter"]]
  cg
}
