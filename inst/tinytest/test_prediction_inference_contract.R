# Prediction uncertainty uses testing df; noise variance still uses residual df.
local({
  set.seed(1204)
  d = data.frame(x = stats::rnorm(200), g = rep(1:5, 40))
  d[["y"]] = 1 + d[["x"]] + stats::rnorm(nrow(d))
  for (strategy in c("moments", "compress")) for (request in list("iid", "hc1", ~g)) {
    for (policy in list(NULL, list(t.df = 3))) {
      model = dbreg(y ~ x, data = d, vcov = request, ssc = policy, strategy = strategy)
      data = d[1:3, , drop = FALSE]
      X = stats::model.matrix(~x, data)
      V = stats::vcov(model)
      variance = as.numeric(Matrix::diag(X %*% V %*% t(X)))
      critical = stats::qt(.975, attr(V, "df.t"))
      noise = sum(stats::residuals(stats::lm(y ~ x, d))^2) / stats::df.residual(stats::lm(y ~ x, d))
      for (interval in c("confidence", "prediction")) {
        actual = predict(model, data, interval = interval)
        expected = critical * sqrt(variance + if (interval == "prediction") noise else 0)
        expect_equal(actual[["upr"]] - actual[["fit"]], expected, tolerance = 1e-10)
        expect_equal(actual[["fit"]] - actual[["lwr"]], expected, tolerance = 1e-10)
      }
    }
  }
})
