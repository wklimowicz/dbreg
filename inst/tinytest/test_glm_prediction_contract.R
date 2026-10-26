# Public GLM newdata must be evaluated from source variables, ignoring extras.
local({
  set.seed(1201)
  d = data.frame(x = stats::rnorm(180), f = rep(letters[1:3], 60),
    g = rep(1:9, each = 20), h = rep(1:4, 45))
  eta = .2 + .2 * d[["x"]] + .1 * (d[["f"]] == "b")
  families = list(stats::gaussian(), stats::binomial(), stats::binomial("cloglog"),
    stats::binomial("probit"), stats::poisson())
  for (family in families) for (absorbed in c(FALSE, TRUE)) {
    d[["y"]] = switch(family[["family"]], gaussian = eta + stats::rnorm(nrow(d)),
      poisson = stats::rpois(nrow(d), exp(eta)), binomial = stats::rbinom(nrow(d), 1, family[["linkinv"]](eta)))
    formula = if (absorbed) y ~ x + factor(f) | g + h else y ~ x + factor(f)
    model = dbglm(formula, data = d, family = family, retain = "prediction", fixef_rm = "none")
    reference = predict(model, d, type = "response")
    extra = d
    for (name in model[["model_spec"]][["columns"]][["alias"]]) extra[[name]] = 100
    expect_equal(predict(model, extra, type = "response"), reference, tolerance = 1e-10)
    close(model)
  }
})
