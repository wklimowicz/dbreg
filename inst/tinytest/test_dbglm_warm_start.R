if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(611)
  d = expand.grid(g = 1:3, h = 1:2, x = -1:1, trend = -1:1, rep = 1:4)
  d$w = runif(nrow(d), .5, 3)
  d$o = d$h / 10
  d$y = rpois(nrow(d), exp(.6 + d$x / 4 + d$g * d$trend / 10 + d$o))
  cases = list(
    list(y ~ x, y ~ x, "compress"),
    list(y ~ x | g + h, y ~ x + factor(g) + factor(h), "moments"),
    list(y ~ x | (g^h)[trend], y ~ x + interaction(g,h) * trend, "compress")
  )
  for (case in cases) {
    reference = glm(case[[2L]], data = d, family = poisson(), weights = w, offset = o,
      control = glm.control(epsilon = 1e-12))
    d$warm = unname(predict(reference, type = "link"))
    model = dbglm(case[[1L]], data = d, family = poisson(), weights = ~w, offset = ~o,
      etastart = ~warm, strategy = case[[3L]], retain = "prediction")
    expect_true(model$converged)
    expect_equal(model$iterations, 1L)
    expect_equal(model$step_halvings, 0L)
    expect_equal(predict(model, d, type = "link"), d$warm, tolerance = 1e-8)
    close(model)
  }

  d$warm = d$warm + sin(seq_len(nrow(d))) / 4
  d$warm[1:2] = c(NA, Inf)
  expected = glm(y ~ x + factor(g), data = d[-(1:2), ], family = poisson(), weights = w, offset = o)
  model = dbglm(y ~ x | g, data = d, family = poisson(), weights = ~w, offset = ~o,
    etastart = "warm", strategy = "moments", retain = "all")
  expect_equal(nobs(model), nrow(d) - 2)
  expect_equal(case.names(model), seq.int(3L, nrow(d)))
  expect_equal(predict(model, d[-(1:2), ], type = "link"), unname(predict(expected, type = "link")),
    tolerance = 1e-8)
  close(model)

  d$large = d$x + 1e12
  d$warm = d$o + .4 * d$x
  model = dbglm(y ~ large | g, data = d, family = poisson(), weights = ~w, offset = ~o,
    etastart = ~warm)
  expect_true("large" %in% names(coef(model)))
  expect_error(dbglm(y ~ x, data = d, etastart = ~I(warm)), "name one column")
  expect_error(dbglm(y ~ x, data = d, etastart = "text"), "text")
}
