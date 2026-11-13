# SQL kernels for the supported built-in GLM families.
dbglm_family = function(family) {
  if (is.character(family) && length(family) == 1L && !is.na(family)) {
    family = switch(family,
      gaussian = stats::gaussian,
      binomial = stats::binomial,
      poisson = stats::poisson,
      stop("Supported family names are gaussian, binomial and poisson.", call. = FALSE)
    )
  }
  if (is.function(family)) family = family()
  if (!inherits(family, "family") ||
    !is.character(family[["link"]]) || length(family[["link"]]) != 1 ||
    is.na(family[["link"]]) ||
    !(identical(family[["family"]], "binomial") && family[["link"]] %in% c("cloglog", "logit", "probit") ||
      identical(family[["family"]], "poisson") && family[["link"]] == "log" ||
      identical(family[["family"]], "gaussian") && family[["link"]] == "identity")) {
    stop("`family` must be binomial('cloglog'), binomial('logit'), binomial('probit'), poisson('log') or gaussian('identity').", call. = FALSE)
  }
  # Literal link arguments keep canonical family closures independent of this
  # validator's evaluation frame and the caller's family object.
  canonical = switch(paste(family[["family"]], family[["link"]], sep = "/"),
    "binomial/logit" = stats::binomial("logit"),
    "binomial/probit" = stats::binomial("probit"),
    "binomial/cloglog" = stats::binomial("cloglog"),
    "poisson/log" = stats::poisson("log"),
    "gaussian/identity" = stats::gaussian("identity")
  )
  functions = c(
    "linkfun", "linkinv", "mu.eta", "variance", "dev.resids", "aic",
    "validmu", "valideta"
  )
  for (name in functions) {
    if (!is.function(family[[name]]) ||
      !identical(body(family[[name]]), body(canonical[[name]])) ||
      !identical(formals(family[[name]]), formals(canonical[[name]]))) {
      stop("Modified GLM family functions are not supported.", call. = FALSE)
    }
  }
  canonical
}

dbglm_number = function(x) sprintf("%.17e", x)

# R's binomial family$aic uses dbinom(round(w*y), round(w), mu). With a
# binary outcome this reduces to round(w) times the Bernoulli log density.
# SQL ROUND commonly rounds halves away from zero, so implement R's ties to
# even explicitly. Weights here are finite and strictly positive. This is a
# reporting convention only: never use the rounded values in IRLS or scores.
dbglm_binomial_loglik = function(conn, working, weight_scale = 1, cell = FALSE) {
  a = paste0("(__dbglm_a*(", dbglm_number(weight_scale), "))")
  rounded = sql_round_even(a)
  term = if (cell) "__dbglm_round_ay*LN(mu)+(__dbglm_round_a-__dbglm_round_ay)*LN(1.0-mu)" else paste0(rounded, "*(CASE WHEN __dbglm_y=1 THEN LN(mu) ELSE LN(1.0-mu) END)")
  as.numeric(DBI::dbGetQuery(conn, paste0(
    "SELECT ", sql_stable_sum(conn, term),
    " AS ll FROM ", working
  ))[["ll"]])
}

# Standard normal CDF using the Hart rational coefficients, evaluated in SQL
# with Horner polynomials and a continued-fraction tail. The scalar result is
# checked against R's pnorm in the numerical tests, including the joins/tails.
dbglm_normal_sql = function(eta, a = paste0("ABS(", eta, ")"), exp_half = paste0("EXP(-0.5*POWER(", eta, ",2))")) {
  horner = function(coef) {
    value = dbglm_number(utils::tail(coef, 1))
    for (k in rev(utils::head(seq_along(coef), -1))) {
      value = paste0("(", dbglm_number(coef[k]), "+", a, "*", value, ")")
    }
    value
  }
  p = c(
    220.2068679123761, 221.2135961699311, 112.0792914978709,
    33.912866078383, 6.37396220353165, .7003830644436881, .03526249659989109
  )
  q = c(
    440.4137358247522, 793.8265125199484, 637.3336333788311,
    296.5642487796737, 86.78073220294608, 16.06417757920695,
    1.755667163182642, .08838834764831844
  )
  density = paste0(exp_half, "/", dbglm_number(sqrt(2 * pi)))
  fraction = paste0(
    "(", a, "+1.0/(", a, "+2.0/(", a, "+3.0/(", a,
    "+4.0/(", a, "+0.65)))))"
  )
  tail = paste0(
    "(CASE WHEN ", a, ">38 THEN 0.0 WHEN ", a, "<",
    dbglm_number(sqrt(50)), " THEN ", exp_half, "*", horner(p),
    "/", horner(q), " ELSE (", density, ")/", fraction, " END)"
  )
  paste0("CASE WHEN ", eta, "<0 THEN ", tail, " ELSE 1.0-", tail, " END")
}

# Reject unusable initial information before any fixed-effect projection.
# Offsets stay part of the immutable data; only free coefficients are damped.
dbglm_valid_start = function(state, link) {
  isTRUE(state[["bad"]] == 0) && is.finite(state[["deviance"]]) &&
    !dbglm_boundary(state, link)
}

# Mean-offset subtraction is a useful cheap default, but heterogeneous offsets
# can make it invalid. Solve the one-dimensional offset-only score in SQL for
# an independent, usable recovery anchor without collecting observation rows.
dbglm_offset_intercept = function(conn, base, link, target, alpha) {
  eta = function(value) paste0("(__dbglm_o+(", dbglm_number(value), "))")
  initial = dbglm_evaluate(conn, base, character(), numeric(), link, eta = eta(alpha))
  if (dbglm_valid_start(initial, link)) return(alpha)
  limits = DBI::dbGetQuery(conn, paste("SELECT MIN(__dbglm_o) AS lo,MAX(__dbglm_o) AS hi FROM", base))
  center = switch(link, logit = stats::qlogis(target), probit = stats::qnorm(target),
    cloglog = log(-log1p(-target)))
  interval = c(center - limits[["hi"]], center - limits[["lo"]])
  if (any(!is.finite(interval)) || interval[1L] >= interval[2L]) return(alpha)
  score = function(value) {
    mu = dbglm_family_sql(link, eta(value))[["mu"]]
    as.numeric(DBI::dbGetQuery(conn, paste0("SELECT SUM(__dbglm_a*(", mu,
      "))/SUM(__dbglm_a) AS mean FROM ", base))[["mean"]]) - target
  }
  stats::uniroot(score, interval, tol = 1e-11, maxiter = 100L)[["root"]]
}

# Main and optional null fits must use the same offset-aware anchor. For
# Poisson, a shifted exponential sum solves the constant score directly without
# overflowing on large offsets or spending one IRLS step per offset unit.
dbglm_initial_intercept = function(conn, base, family, moments, offset_present = FALSE) {
  alpha = family[["linkfun"]](moments[["ybar"]]) - moments[["offset_mean"]]
  if (!offset_present || !is.finite(alpha)) return(alpha)
  if (family[["family"]] == "binomial") {
    return(dbglm_offset_intercept(conn, base, family[["link"]], moments[["ybar"]], alpha))
  }
  if (family[["family"]] == "poisson") {
    maximum = as.numeric(DBI::dbGetQuery(conn,
      paste("SELECT MAX(__dbglm_o) AS maximum FROM", base))[["maximum"]])
    shifted = paste0("EXP(__dbglm_o-(", dbglm_number(maximum), "))")
    mean = as.numeric(DBI::dbGetQuery(conn, paste0("SELECT ",
      sql_stable_sum(conn, paste0("__dbglm_a*", shifted)), "/",
      sql_stable_sum(conn, "__dbglm_a"), " AS mean FROM ", base))[["mean"]])
    return(log(moments[["ybar"]]) - log(mean) - maximum)
  }
  alpha
}

dbglm_info_factor = function(H) {
  if (!length(H) || any(!is.finite(H))) return(list(R = NULL, scale = numeric()))
  scale = sqrt(pmax(diag(H), 0))
  if (!length(scale) || any(!is.finite(scale)) || any(scale <= 0)) return(list(R = NULL, scale = numeric()))
  R = tryCatch(chol(H / outer(scale, scale)), error = function(e) NULL)
  list(R = R, scale = scale)
}

dbglm_min_info_ratio = function(H, R0, scale0) {
  if (is.null(R0) || !length(scale0)) return(NA_real_)
  if (!isTRUE(nrow(H) == length(scale0)) || any(!is.finite(H))) return(0)
  C = H / outer(scale0, scale0)
  Z = forwardsolve(t(R0), C)
  M = forwardsolve(t(R0), t(Z))
  M = (M + t(M)) / 2
  eigs = eigen(M, symmetric = TRUE, only.values = TRUE)[["values"]]
  if (any(!is.finite(eigs))) return(0)
  min(eigs)
}

dbglm_info_collapsed = function(H, R0, scale0, tol = 1e-8) {
  ratio = dbglm_min_info_ratio(H, R0, scale0)
  is.finite(ratio) && ratio <= tol
}

dbglm_boundary = function(state, link, collapsed = FALSE, stalled = FALSE) {
  if (link == "identity") {
    !is.finite(state[["mu_min"]]) || !is.finite(state[["mu_max"]])
  } else if (link == "log") {
    !is.finite(state[["mu_min"]]) || !is.finite(state[["mu_max"]]) || state[["mu_min"]] <= 0
  } else if (all(vapply(c("eta_min", "eta_max", "misfit"), function(k) is.numeric(state[[k]]) && is.finite(state[[k]]), logical(1)))) {
    # A finite fit can saturate the clamped mean where the outcome agrees
    # (cloglog reaches mu = 1 - eps at eta ~ 3.6, probit at |eta| ~ 8); those
    # rows add ~0 deviance, correctly. A row saturated against its outcome has
    # a clamped, wrong deviance (misfit > 0). True separation collapses Fisher
    # information along a receding direction while Newton steps stay stalled > 0.
    state[["misfit"]] > 0 || (isTRUE(collapsed) && isTRUE(stalled))
  } else state[["mu_min"]] <= .Machine[["double.eps"]] || state[["mu_max"]] >= 1 - .Machine[["double.eps"]]
}

# Rows whose clamped mean saturates against their outcome, so their deviance
# term is clamped and wrong. Such a row carries deviance 2a*ln(1/eps) ~ 72a; a
# row saturated in agreement with its outcome carries ~4e-16a. Reading `dev`
# keeps this to the columns every working relation has.
dbglm_misfit_sql = function(alias = "misfit") paste0(
  "SUM(CASE WHEN (mu >= ", dbglm_number(1 - .Machine[["double.eps"]]), " OR mu <= ",
  dbglm_number(.Machine[["double.eps"]]), ") AND dev > 1e-8 * __dbglm_a THEN 1 ELSE 0 END) AS ", alias)

# Finite fits can saturate the clamped binomial mean (cloglog at eta ~ 3.6), so
# saturation is not separation by itself. It does explain an IRLS that runs out
# of iterations: report that as the boundary rather than as plain nonconvergence.
dbglm_nonconvergence = function(state, link, max_iter) {
  eps = .Machine[["double.eps"]]
  saturated = link %in% c("logit", "probit", "cloglog") &&
    isTRUE(state[["mu_min"]] <= eps || state[["mu_max"]] >= 1 - eps)
  if (saturated || link == "log") stop(
    "GLM did not converge within ", max_iter, " iterations; a finite solution is not established ",
    "(possible separation or numerical probability boundary).", call. = FALSE)
  stop("GLM did not converge within ", max_iter, " iterations.", call. = FALSE)
}

dbglm_family_sql = function(link, eta = "eta", a = NULL, exp_half = NULL) {
  eps = dbglm_number(.Machine[["double.eps"]])
  if (link == "identity") return(list(mu = eta, derivative = "1.0"))
  if (link == "logit") {
    mu = sprintf(
      paste0(
        "CASE WHEN %1$s < -30 THEN %2$s WHEN %1$s > 30 THEN %3$s ",
        "ELSE EXP(%1$s) / (1.0 + EXP(%1$s)) END"
      ),
      eta, dbglm_number(.Machine[["double.eps"]] / (1 + .Machine[["double.eps"]])),
      dbglm_number(1 / (1 + .Machine[["double.eps"]]))
    )
    derivative = sprintf(paste0(
      "CASE WHEN ABS(%1$s) > 30 THEN %2$s ",
      "ELSE EXP(%1$s) / POWER(1.0 + EXP(%1$s), 2) END"
    ), eta, eps)
  } else if (link == "probit") {
    a_expr = if (is.null(a)) paste0("ABS(", eta, ")") else a
    eh_expr = if (is.null(exp_half)) paste0("EXP(-0.5*POWER(", eta, ",2))") else exp_half
    mu = paste0(
      "GREATEST(LEAST((", dbglm_normal_sql(eta, a_expr, eh_expr), "),",
      dbglm_number(1 - .Machine[["double.eps"]]), "),", eps, ")"
    )
    derivative = paste0(
      "CASE WHEN ", a_expr, ">38 THEN ", eps,
      " ELSE GREATEST(", eh_expr, "/", dbglm_number(sqrt(2 * pi)), ",", eps, ") END"
    )
  } else if (link == "log") {
    mu = paste0("CASE WHEN ", eta, ">700 THEN 'Infinity'::DOUBLE ELSE EXP(", eta, ") END")
    derivative = mu
  } else if (link == "cloglog") {
    h = if (is.null(exp_half)) sprintf("EXP(LEAST(%s, 700.0))", eta) else exp_half
    # expm1 is unavailable on some DuckDB builds. This series avoids the
    # cancellation in 1-exp(-h) for tiny h, to double precision at this cutoff.
    probability = sprintf(paste0(
      "CASE WHEN %1$s < -18 THEN ",
      "(%2$s - POWER(%2$s, 2)/2.0 + POWER(%2$s, 3)/6.0) ",
      "ELSE 1.0 - EXP(-(%2$s)) END"
    ), eta, h)
    mu = sprintf(
      "GREATEST(LEAST((%s), %s), %s)", probability,
      dbglm_number(1 - .Machine[["double.eps"]]), eps
    )
    derivative = sprintf("GREATEST((%s) * EXP(-(%s)), %s)", h, h, eps)
  } else stop("Unsupported GLM link.", call. = FALSE)
  list(mu = mu, derivative = derivative)
}

# newton = TRUE replaces the Fisher (expected-information) working weight by the
# observed information where dbglm_newton_link() supports it. Scores and
# deviances are unchanged; covariance consumers keep the default Fisher rows.
dbglm_working_sql = function(table, x, beta, link, eta = NULL, conn,
                             backend_info = detect_backend(conn), newton = FALSE, cell = FALSE,
                             materialize = FALSE) {
  if (!cell) {
    sql = dbglm_working_rows_sql(table, x, beta, link, eta, conn = conn,
      backend_info = backend_info, newton = newton)
  } else {
    rows = dbglm_working_rows_sql(table, x, beta, link, eta, working_name = "cell_working",
      conn = conn, backend_info = backend_info, newton = newton)
    deviance = if (link %in% c("identity", "log")) "dev+__dbglm_dev_const" else
      "-2.0*__dbglm_a*(__dbglm_y*LN(mu)+(1.0-__dbglm_y)*LN(1.0-mu))"
    sql = paste0(rows, ", working AS (SELECT ", sql_select_except(conn, "dev", backend_info),
      ",", deviance, " AS dev,", dbglm_cell_score_square(link),
      " AS score_sq FROM cell_working) ")
  }
  if (materialize) {
    sql = sub(", working AS (", ", working AS MATERIALIZED (", sql, fixed = TRUE)
  }
  sql
}

dbglm_working_rows_sql = function(table, x, beta, link, eta = NULL, working_name = "working", conn, backend_info = detect_backend(conn), newton = FALSE) {
  if (is.null(eta)) eta = paste(c("__dbglm_o", sprintf("(%s)*(%s)", dbglm_number(beta), x)), collapse = " + ")
  if (link == "identity") return(paste0(
    "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
    "means AS (SELECT *,eta AS mu,1.0 AS d FROM predictor), ",
    working_name, " AS (SELECT *,__dbglm_a AS w,__dbglm_y-__dbglm_o AS z,",
    "__dbglm_a*(__dbglm_y-mu) AS score,__dbglm_a*POWER(__dbglm_y-mu,2) AS dev FROM means) "
  ))
  if (link == "log") {
    kernel = dbglm_family_sql(link)
    return(paste0(
      "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
      "means AS (SELECT *, ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM predictor), ",
      working_name, " AS (SELECT *, __dbglm_a*mu AS w, eta-__dbglm_o+(__dbglm_y-mu)/mu AS z, ",
      "__dbglm_a*(__dbglm_y-mu) AS score, ",
      "2.0*__dbglm_a*(", sql_poisson_bd0("__dbglm_y", "mu"), ") AS dev FROM means) "
    ))
  }
  if (link == "probit") {
    kernel = dbglm_family_sql(link, eta = "eta", a = "__dbglm_abs_eta", exp_half = "__dbglm_eh")
    return(paste0(
      "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
      "prep AS (SELECT *, ABS(eta) AS __dbglm_abs_eta, EXP(-0.5*POWER(eta,2)) AS __dbglm_eh FROM predictor), ",
      "means AS (SELECT ", sql_select_except(conn, c("__dbglm_abs_eta", "__dbglm_eh"), backend_info), ", ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM prep), ",
      dbglm_binomial_working_cte(working_name)
    ))
  }
  if (link == "cloglog") {
    kernel = dbglm_family_sql(link, eta = "eta", exp_half = "__dbglm_eh")
    if (isTRUE(newton)) return(paste0(
      "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
      "prep AS (SELECT *, EXP(LEAST(eta, 700.0)) AS __dbglm_eh FROM predictor), ",
      "kernel AS (SELECT *, ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM prep), ",
      "means AS (SELECT ", sql_select_except(conn, "__dbglm_eh", backend_info), ", ", dbglm_cloglog_observed_sql(), " AS __dbglm_wn FROM kernel), ",
      dbglm_binomial_working_cte(working_name, newton = TRUE, conn = conn, backend_info = backend_info)
    ))
    return(paste0(
      "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
      "prep AS (SELECT *, EXP(LEAST(eta, 700.0)) AS __dbglm_eh FROM predictor), ",
      "means AS (SELECT ", sql_select_except(conn, "__dbglm_eh", backend_info), ", ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM prep), ",
      dbglm_binomial_working_cte(working_name)
    ))
  }
  kernel = dbglm_family_sql(link)
  paste0(
    "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
    "means AS (SELECT *, ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM predictor), ",
    dbglm_binomial_working_cte(working_name)
  )
}

# All binomial links share working weights, response, score and deviance.
# Link-specific code only constructs the mean and its derivative. Newton rows
# take the observed information from means.__dbglm_wn instead of d^2/V.
dbglm_binomial_working_cte = function(working_name = "working", newton = FALSE, conn = NULL, backend_info = NULL) {
  score = "__dbglm_a*d*(__dbglm_y-mu)/(mu*(1.0-mu))"
  dev = "-2.0*__dbglm_a*(CASE WHEN __dbglm_y=1 THEN LN(mu) ELSE LN(1.0-mu) END) AS dev "
  if (isTRUE(newton)) return(paste0(
    working_name, " AS (SELECT ", sql_select_except(conn, "__dbglm_wn", backend_info), ", __dbglm_a*__dbglm_wn AS w, ",
    "eta-__dbglm_o+(", score, ")/(__dbglm_a*__dbglm_wn) AS z, ",
    score, " AS score, ", dev, "FROM means) "
  ))
  paste0(
    working_name, " AS (SELECT *, __dbglm_a*d*d/(mu*(1.0-mu)) AS w, ",
    "eta-__dbglm_o+(__dbglm_y-mu)/d AS z, ",
    score, " AS score, ", dev, "FROM means) "
  )
}

# Links whose observed information is non-negative for responses in [0, 1],
# so Newton steps need no safeguard beyond the shared step halving.
dbglm_newton_link = function(link) identical(link, "cloglog")

# Observed information per unit weight for the binomial cloglog likelihood,
# with t = exp(eta) and mu = 1 - exp(-t):
#   -d2l/deta2 = t * ((1 - y) + y * (t - mu) * (1 - mu) / mu^2)  >= 0.
# For eta < -18, (t - mu)(1 - mu)/mu^2 = 1/2 - t/6 + O(t^2) avoids cancellation.
# Reads the cloglog kernel columns __dbglm_eh (t), mu and __dbglm_y.
dbglm_cloglog_observed_sql = function() {
  eps = dbglm_number(.Machine[["double.eps"]])
  paste0(
    "GREATEST(__dbglm_eh*((1.0-__dbglm_y)+__dbglm_y*(CASE WHEN eta < -18 THEN 0.5-__dbglm_eh/6.0 ",
    "ELSE (__dbglm_eh-mu)*(1.0-mu)/(mu*mu) END)), ", eps, ")"
  )
}

dbglm_deviance_term = function(conn) paste(sql_stable_sum(conn, "dev"), "AS deviance")

dbglm_saturated_loglik = function(conn, base) {
  backend_info = detect_backend(conn)
  # Databricks has no lgamma builtin; the portable small-y series replaces it.
  term = dbglm_saturated_loglik_term(lgamma = !isTRUE(backend_info[["is_databricks"]]))
  DBI::dbGetQuery(conn, paste0("SELECT ", sql_sum(conn, term, backend_info), " AS ll FROM ", base))[["ll"]]
}

# Saturated Poisson term y*ln(y) - y - lgamma(y + 1), valid for non-integer
# outcomes. Large y uses Stirling's expansion of the difference (avoiding
# O(y log y) cancellation). Small y uses lgamma, or without a builtin
# lgamma(y + 1) = lgamma(y + 21) - sum(ln(y + k), k = 1..20) with Stirling at y + 21.
dbglm_saturated_loglik_term = function(lgamma = TRUE, y = "__dbglm_y", a = "__dbglm_a") {
  stirling = paste0(
    "-0.5*LN(", dbglm_number(2 * pi), "*", y, ")-1.0/(12.0*", y, ")",
    "+1.0/(360.0*POWER(", y, ",3))-1.0/(1260.0*POWER(", y, ",5))",
    "+1.0/(1680.0*POWER(", y, ",7))"
  )
  small = if (lgamma) paste0("lgamma(", y, "+1)") else {
    z = paste0("(", y, "+21.0)")
    paste0("((", z, "-0.5)*LN(", z, ")-", z, "+", dbglm_number(0.5 * log(2 * pi)),
      "+1.0/(12.0*", z, ")-1.0/(360.0*POWER(", z, ",3))+1.0/(1260.0*POWER(", z, ",5))",
      "-1.0/(1680.0*POWER(", z, ",7))-LN(", paste0("(", y, "+", 1:20, ".0)", collapse = "*"), "))")
  }
  paste0(a, "*(CASE WHEN ", y, "=0 THEN 0.0 WHEN ", y, ">=20 THEN (", stirling, ") ELSE ",
    y, "*LN(", y, ")-", y, "-", small, " END)")
}

# Shared profile Gaussian reporting primitive. Keep normalized RSS and its
# binary weight scale separate; multiplying them first can overflow.
dbglm_gaussian_loglik = function(n, deviance, weight_scale, log_weights) {
  -n / 2 * (log(2 * pi) + log(deviance / n) + log(weight_scale) + 1) + log_weights / 2
}

sql_round_even = function(value) {
  paste0(
    "(CASE WHEN ", value, "-FLOOR(", value, ")<0.5 THEN FLOOR(", value, ") ",
    "WHEN ", value, "-FLOOR(", value, ")>0.5 THEN CEIL(", value, ") ",
    "WHEN MOD(FLOOR(", value, "),2.0)=0 THEN FLOOR(", value, ") ELSE CEIL(", value, ") END)"
  )
}
