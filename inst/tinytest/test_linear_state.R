# Mathematical fit state must not inherit a caller's data/connection lifetime.
local({
  set.seed(614)
  d = data.frame(x = rep(c(-1, 1), 200), g = rep(1:10, each = 40),
    unused = paste0(seq_len(400), strrep("x", 2000)))
  d[["y"]] = .5 * d[["x"]] + d[["g"]] + rnorm(400)
  formula_env = new.env(parent = baseenv())
  formula_env[["captured_panel"]] = d
  fml = stats::as.formula("y ~ x + I(x^2) | g", env = formula_env)
  for (strategy in c("compress", "demean", "mundlak", "moments")) {
    current = if (strategy == "moments") stats::as.formula("y ~ x", env = formula_env) else fml
    fit = dbreg(current, data = d, strategy = strategy)
    expect_identical(environment(fit[["fml"]]), baseenv())
    if (strategy == "compress") {
      expect_identical(environment(fit[["compressed_design"]][["formula"]]), baseenv())
    }
    expect_identical(environment(current), formula_env)
    saved = serialize(fit, NULL)
    expect_true(length(saved) < 100000)
    expected = predict(fit, d)
    # Mutating a captured panel must not mutate the serialized fitted object.
    formula_env[["captured_panel"]][["y"]] = -1e6
    expect_identical(serialize(fit, NULL), saved)
    expect_equal(predict(unserialize(saved), d), expected, tolerance = 1e-10)
  }
})
