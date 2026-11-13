# dbglm gives the stats::glm answer: coefficients and deviance from a tightly
# converged glm (fixed effects as dummies), and iid / HC1 / cluster standard
# errors from the glm sandwich under fixest's "full" small-sample convention.
local({
  set.seed(20261117)
  n = 600
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(12, n, TRUE),
    cl = sample(20, n, TRUE), w = runif(n, .5, 2), o = runif(n, -.3, .3))
  eta = .3*d$x - .2*d$z + rnorm(12, 0, .3)[d$g]
  families = list(poisson = poisson(), logit = binomial(), probit = binomial("probit"),
    cloglog = binomial("cloglog"))
  glm_se = function(g, data, v) {
    X = model.matrix(g); n = nrow(X); K = ncol(X); pw = weights(g, "prior")
    fam = family(g); mu = fitted(g); eta = g$linear.predictors
    bread = solve(crossprod(X * sqrt(pw * fam$mu.eta(eta)^2 / fam$variance(mu))))
    S = X * (pw * (g$y - mu) * fam$mu.eta(eta) / fam$variance(mu))
    V = switch(v,
      iid = bread * (n - 1)/(n - K),
      hc1 = bread %*% crossprod(S) %*% bread * n/(n - K),
      cl = { C = rowsum(S, data$cl); G = nrow(C)
        bread %*% crossprod(C) %*% bread * G/(G - 1) * (n - 1)/(n - K) })
    sqrt(diag(V))[c("x", "z")]
  }
  for (fn in names(families)) {
    fam = families[[fn]]
    d$y = if (fn == "poisson") rpois(n, exp(eta + d$o)) else
      rbinom(n, 1, fam$linkinv(eta + d$o - if (fn == "cloglog") 1 else 0))
    for (fe in c(FALSE, TRUE)) for (wt in c(FALSE, TRUE)) {
      gd = transform(d, ww = if (wt) w else 1)
      ref = suppressWarnings(glm(if (fe) y ~ x + z + factor(g) else y ~ x + z, fam, gd,
        weights = ww, offset = o, control = glm.control(epsilon = 1e-15, maxit = 200)))
      for (v in c("iid", "hc1", "cl")) {
        info = paste(fn, if (fe) "fe" else "nofe", if (wt) "weighted" else "", v)
        fit = dbglm(if (fe) y ~ x + z | g else y ~ x + z, data = d, family = fam,
          weights = if (wt) "w", offset = "o", vcov = if (v == "cl") ~cl else v)
        expect_equal(coef(fit)[c("x", "z")], coef(ref)[c("x", "z")], tolerance = 1e-8, info = info)
        expect_equal(sqrt(diag(as.matrix(vcov(fit))))[c("x", "z")], glm_se(ref, gd, v),
          tolerance = 1e-7, info = info)
        if (v == "iid") expect_equal(fit$deviance, deviance(ref), tolerance = 1e-9, info = info)
      }
    }
  }
})

# Learned GLM fixed effects: new-data predictions and fixef() match fixest.
if (requireNamespace("fixest", quietly = TRUE)) local({
  set.seed(20261117)
  n = 800
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(30, n, TRUE), h = sample(10, n, TRUE), w = runif(n, .5, 2))
  d$yp = rpois(n, exp(.3*d$x - .2*d$z + rnorm(30, 0, .3)[d$g] + rnorm(10, 0, .3)[d$h]))
  nd = d[sample(n, 200), ]
  for (rhs in c("x + z | g", "x + z | g + h")) {
    f = as.formula(paste("yp ~", rhs))
    ref = fixest::feglm(f, d, family = poisson(), weights = ~w, glm.tol = 1e-12, fixef.tol = 1e-11, fixef.rm = "none")
    fit = dbglm(f, data = d, family = poisson(), weights = "w", retain = "prediction")
    expect_equal(unname(predict(fit, nd, type = "response")), unname(predict(ref, nd, type = "response")),
      tolerance = 1e-7, info = rhs)
    if (rhs == "x + z | g") expect_equal(fixest::fixef(fit)[["g"]][names(fixest::fixef(ref)$g)],
      fixest::fixef(ref)$g, tolerance = 1e-7)
    close(fit)
  }
})

## Gaussian dbglm certifies convergence on the outcome's scale: prices in the millions fit like OLS.
local({
  set.seed(20261117)
  n = 2000
  d = data.frame(x = rnorm(n), county = sample(40, n, TRUE), year = sample(10, n, TRUE))
  d$y = 1e7 + 4e6 * d$x + 2e6 * rnorm(40)[d$county] + 1e6 * rnorm(10)[d$year] + 3e6 * rnorm(n)
  for (f in list(y ~ x, y ~ x | county, y ~ x | county + year)) {
    ref = fixest::feols(f, d, vcov = "hc1", ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none")
    fit = dbglm(f, data = d, family = gaussian(), vcov = "hc1")
    expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-10, info = deparse(f))
    expect_equal(sqrt(vcov(fit)[["x", "x"]]), fixest::se(ref)[["x"]], tolerance = 1e-8, info = deparse(f))
  }
})

## Keyed fitted values and residuals collect to numeric vectors in row_id order.
local({
  set.seed(3)
  n = 500
  d = data.frame(rid = sample(n), x = rnorm(n), g = sample(20, n, TRUE))
  d$y = rpois(n, exp(.3 * d$x))
  con = DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(con, "d", d)
  ref = stats::glm(y ~ x + factor(g), poisson(), d)
  fit = dbglm(y ~ x | g, conn = con, table = "d", family = poisson(), retain = "inference", row_id = "rid")
  rows = fitted(fit)
  expect_equal(as.numeric(rows), unname(fitted(ref))[order(d$rid)], tolerance = 1e-8)
  expect_true(cor(as.numeric(rows), d$y[order(d$rid)]) > 0)
  expect_error(as.numeric(case.names(fit)), "no values")
  close(rows)
  close(fit)
  DBI::dbDisconnect(con, shutdown = TRUE)
})

## A tolerated non-converged fit warns, like glm(); the separation refusal is a clean error.
local({
  set.seed(4)
  n = 400
  d = data.frame(x = rnorm(n), g = sample(10, n, TRUE))
  d$y = rpois(n, exp(.4 * d$x))
  expect_warning(dbglm(y ~ x | g, data = d, family = poisson(),
    control = dbglm_control(glm_max_iter = 1L, tolerate_max_iter = TRUE)), "did not converge")
  d$y[d$g == 1] = 0
  seen = character()
  err = withCallingHandlers(tryCatch(dbglm(y ~ x | g, data = d, family = poisson(), fixef_rm = "singletons"),
    error = function(e) e), warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_inherits(err, "error")
  expect_equal(length(seen), 0L)
})

## One plain fixed effect uses the exact window projection on every backend.
local({
  set.seed(6)
  n = 3000
  d = data.frame(x = rnorm(n), g = sample(300, n, TRUE), w = runif(n, .5, 2))
  d$y = rpois(n, exp(.3 * d$x + rnorm(300, 0, .3)[d$g]))
  ref = suppressWarnings(stats::glm(y ~ x + factor(g), poisson(), d, weights = w,
    control = stats::glm.control(epsilon = 1e-15, maxit = 100)))
  fit = dbglm(y ~ x | g, data = d, family = poisson(), weights = ~w)
  expect_identical(fit[["projection_method"]], "sql_window")
  expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-9)
  expect_equal(fit[["deviance"]], deviance(ref), tolerance = 1e-9)
})

## Many-level fixed effects keep their coefficient lookups in owned DuckDB tables.
local({
  set.seed(8)
  n = 8000
  d = data.frame(x = rnorm(n), g = sample(1500, n, TRUE), h = sample(8, n, TRUE))
  d$y = rpois(n, exp(.3 * d$x + rnorm(1500, 0, .3)[d$g] + rnorm(8, 0, .3)[d$h]))
  con = DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(con, "d", d)
  ref = fixest::fepois(y ~ x | g + h, d, glm.tol = 1e-12, fixef.tol = 1e-11)
  fit = dbglm(y ~ x | g + h, conn = con, table = "d", family = poisson())
  expect_identical(fit[["projection_method"]], "bounded_moments")
  expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-8)
  expect_identical(DBI::dbListTables(con), "d")
  DBI::dbDisconnect(con, shutdown = TRUE)
})
