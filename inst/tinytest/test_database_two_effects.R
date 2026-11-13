if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("bit64", quietly = TRUE)) {
  set.seed(942)
  d = expand.grid(a = 1:8, b = 1:6, x = -1:1, rep = 1:4)
  d = d[d$a %% 2 == d$b %% 2, , drop = FALSE]
  d$id = seq_len(nrow(d))
  d$g = bit64::as.integer64("9007199254740992") + d$a
  d$h = bit64::as.integer64("9007199254740992") + d$b
  d$tag = d$a - d$b
  d$o = d$rep / 4
  d$w = runif(nrow(d), 0.5, 2)
  d$y = rpois(nrow(d), exp(0.3 + d$x / 5 + d$tag / 10 + d$o))
  previous = options(dbreg.prediction_max_entries = 9)
  models = list(
    dbreg(y ~ x + tag | g + h, data = d, weights = ~w, offset = ~o, strategy = "demean", retain = "prediction"))
  if (exists("dbglm", mode = "function")) {
    models = c(models, list(dbglm(y ~ x + tag | g + h, data = d, weights = ~w, offset = ~o, family = poisson(), retain = "prediction")))
  }
  options(previous)
  rows = c(16L, 1L, 12L, 40L)
  data = d[rows, c("id", "x", "g", "h", "tag", "o")]
  data$h[2L] = bit64::as.integer64("9007199254740994")
  data$tag[3L] = data$tag[3L] + 1
  for (fit in models) {
    expect_identical(fit$prediction_diagnostics$method, "database_prediction")
    ref = if (inherits(fit, "dbglm")) {
      glm(y ~ x + factor(a) + factor(b), data = d, weights = w, offset = o, family = poisson())
    } else lm(y ~ x + factor(a) + factor(b), data = d, weights = w, offset = o)
    expected = unname(fitted(ref)[rows])
    expected[2:3] = NA_real_
    expect_equal(predict(fit, data), expected, tolerance = 1e-8)
    if (requireNamespace("fixest", quietly = TRUE)) {
      effects = fixest::fixef(fit)
      expect_identical(names(effects), c("g", "h"))
      eta = unname(coef(fit)[["x"]] * d$x + d$o +
        effects$g[as.character(d$g)] + effects$h[as.character(d$h)])
      want = if (inherits(fit, "dbglm")) predict(fit, d, type = "link") else predict(fit, d)
      expect_equal(eta, want, tolerance = 1e-9)
    }
    output = predict(fit, data, row_id = "id")
    close(fit)
    expect_equal(as.data.frame(output)[[output$prediction_column]], expected[order(rows)], tolerance = 1e-8)
    close(output)
  }
}
