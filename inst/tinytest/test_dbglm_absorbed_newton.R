local({
  set.seed(11144)
  d = expand.grid(g = 1:5, h = 1:4, rep = 1:12)
  d[["x"]] = rnorm(nrow(d))
  d[["w"]] = runif(nrow(d), .3, 3)
  d[["o"]] = runif(nrow(d), -.2, .2)
  d[["cl"]] = sample(1:15, nrow(d), TRUE)
  family = binomial("cloglog")
  eta = -.3 + .7 * d[["x"]] + cos(d[["g"]]) / 4 + sin(d[["h"]]) / 3 + d[["o"]]
  d[["y"]] = rbinom(nrow(d), 1, family[["linkinv"]](eta))
  fits = lapply(c("fisher", "newton"), function(method) {
    dbglm(y ~ x | g + h, data = d, family = family, weights = ~w, offset = ~o,
      vcov = ~cl, control = dbglm_control(irls = method, score_tol = 1e-9))
  })
  expect_identical(vapply(fits, `[[`, character(1), "irls"), c("fisher", "newton"))
  expect_true(all(vapply(fits, `[[`, logical(1), "converged")))
  expect_equal(coef(fits[[1L]]), coef(fits[[2L]]), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fits[[1L]])), as.numeric(vcov(fits[[2L]])), tolerance = 1e-8)
  expect_equal(fits[[1L]][["deviance"]], fits[[2L]][["deviance"]], tolerance = 1e-10)
  expect_true(all(vapply(fits, function(fit) max(fit[["score_max"]], fit[["fe_score_max"]]) <= 1e-9,
    logical(1))))
})
