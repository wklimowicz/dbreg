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

dbglm_number = function(x) sprintf("%.17g", x)

# R's binomial family$aic uses dbinom(round(w*y), round(w), mu). With a
# binary outcome this reduces to round(w) times the Bernoulli log density.
# SQL ROUND commonly rounds halves away from zero, so implement R's ties to
# even explicitly. Weights here are finite and strictly positive. This is a
# reporting convention only: never use the rounded values in IRLS or scores.
dbglm_binomial_loglik = function(conn, working, weight_scale = 1, cell = FALSE) {
  a = paste0("(__dbglm_a*(", dbglm_number(weight_scale), "))")
  rounded = sql_round_even(a)
  term = if (cell) "__dbglm_round_ay*LN(mu)+(__dbglm_round_a-__dbglm_round_ay)*LN(1.0-mu)" else
    paste0(rounded, "*(CASE WHEN __dbglm_y=1 THEN LN(mu) ELSE LN(1.0-mu) END)")
  as.numeric(backend_get_query(conn, paste0(
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

dbglm_boundary = function(state, link) {
  if (link == "identity") {
    !is.finite(state[["mu_min"]]) || !is.finite(state[["mu_max"]])
  } else if (link == "log") {
    !is.finite(state[["mu_min"]]) || !is.finite(state[["mu_max"]]) || state[["mu_min"]] <= 0
  } else state[["mu_min"]] <= .Machine[["double.eps"]] || state[["mu_max"]] >= 1 - .Machine[["double.eps"]]
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

dbglm_working_sql = function(table, x, beta, link, eta = NULL, cell = FALSE) {
  if (!cell) return(dbglm_working_rows_sql(table, x, beta, link, eta))
  rows = dbglm_working_rows_sql(table, x, beta, link, eta, working_name = "dbglm_cell_working")
  deviance = if (link %in% c("identity", "log")) "dev+__dbglm_dev_const" else
    "-2.0*__dbglm_a*(__dbglm_y*LN(mu)+(1.0-__dbglm_y)*LN(1.0-mu))"
  paste0(
    rows, ", working AS (SELECT * EXCLUDE(dev),", deviance, " AS dev,",
    dbglm_cell_score_square(link), " AS score_sq FROM dbglm_cell_working) "
  )
}

dbglm_working_rows_sql = function(table, x, beta, link, eta = NULL, working_name = "working") {
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
      "2.0*__dbglm_a*(CASE WHEN __dbglm_y=0 THEN mu ELSE ",
      "__dbglm_y*(LN(__dbglm_y)-eta)-__dbglm_y+mu END) AS dev FROM means) "
    ))
  }
  if (link == "probit") {
    kernel = dbglm_family_sql(link, eta = "eta", a = "__dbglm_abs_eta", exp_half = "__dbglm_eh")
    return(paste0(
      "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
      "prep AS (SELECT *, ABS(eta) AS __dbglm_abs_eta, EXP(-0.5*POWER(eta,2)) AS __dbglm_eh FROM predictor), ",
      "means AS (SELECT * EXCLUDE (__dbglm_abs_eta, __dbglm_eh), ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM prep), ",
      dbglm_binomial_working_cte(working_name)
    ))
  }
  if (link == "cloglog") {
    kernel = dbglm_family_sql(link, eta = "eta", exp_half = "__dbglm_eh")
    return(paste0(
      "WITH predictor AS (SELECT *, CAST((", eta, ") AS DOUBLE) AS eta FROM ", table, "), ",
      "prep AS (SELECT *, EXP(LEAST(eta, 700.0)) AS __dbglm_eh FROM predictor), ",
      "means AS (SELECT * EXCLUDE (__dbglm_eh), ", kernel[["mu"]], " AS mu, ", kernel[["derivative"]], " AS d FROM prep), ",
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
# Link-specific code only constructs the mean and its derivative.
dbglm_binomial_working_cte = function(working_name = "working") {
  paste0(
    working_name, " AS (SELECT *, __dbglm_a*d*d/(mu*(1.0-mu)) AS w, ",
    "eta-__dbglm_o+(__dbglm_y-mu)/d AS z, ",
    "__dbglm_a*d*(__dbglm_y-mu)/(mu*(1.0-mu)) AS score, ",
    "-2.0*__dbglm_a*(CASE WHEN __dbglm_y=1 THEN LN(mu) ELSE LN(1.0-mu) END) AS dev ",
    "FROM means) "
  )
}

dbglm_deviance_term = function(conn) paste(sql_stable_sum(conn, "dev"), "AS deviance")

dbglm_saturated_loglik = function(conn, base) {
  backend_info = detect_backend(conn)
  if (!isTRUE(backend_info[["is_databricks"]])) return(backend_get_query(conn, paste0(
    "SELECT ", sql_sum(conn, paste0("__dbglm_a*(CASE WHEN __dbglm_y=0 THEN 0.0 ELSE ",
      "__dbglm_y*LN(__dbglm_y)-__dbglm_y-lgamma(__dbglm_y+1) END)"), backend_info), " AS ll FROM ", base
  ))[["ll"]])
  # Databricks has no lgamma builtin. A bounded count histogram is sufficient
  # for the likelihood constant, with R's canonical lgamma evaluation.
  n = backend_get_query(conn, paste("SELECT COUNT(DISTINCT __dbglm_y) AS n FROM", base))[["n"]]
  if (n > 1e6) stop("Poisson count histogram exceeds the likelihood transfer limit.", call. = FALSE)
  h = backend_get_query(conn, paste("SELECT __dbglm_y AS y,SUM(__dbglm_a) AS a FROM", base, "GROUP BY __dbglm_y"))
  value = numeric(nrow(h))
  positive = h[["y"]] > 0
  value[positive] = h[["y"]][positive] * log(h[["y"]][positive]) - h[["y"]][positive] - lgamma(h[["y"]][positive] + 1)
  sum(h[["a"]] * value)
}
