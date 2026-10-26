local({
  set.seed(11013)
  d = data.frame(x = stats::rnorm(90), z = stats::rnorm(90))
  for (term in c("main", "power", "product")) {
    column = switch(term, main = d[["x"]], power = d[["x"]]^2, product = d[["x"]] * d[["z"]])
    d[["y"]] = 1.3 - 0.7 * column + stats::rnorm(90, sd = 0.1)
    beta = as.numeric(qr.solve(cbind(1, column), d[["y"]]))
    formula = switch(term, main = y ~ x, power = y ~ I(x^2), product = y ~ x:z)
    for (strategy in c("moments", "compress")) {
      fit = dbreg::dbreg(formula, data = d, strategy = strategy)
      accepted = list(c(100, 200), c(100L, 200L))
      if (requireNamespace("bit64", quietly = TRUE)) accepted[[length(accepted) + 1L]] = bit64::as.integer64(c("100", "200"))
      if (term != "power") {
        accepted[[length(accepted) + 1L]] = as.Date("1970-01-01") + c(100, 200)
        accepted[[length(accepted) + 1L]] = as.POSIXct("1970-01-01", tz = "UTC") + c(100, 200)
      } else accepted[[length(accepted) + 1L]] = c(TRUE, FALSE)
      for (value in accepted) {
        rows = data.frame(z = c(0.3, -0.2))
        rows[["x"]] = value
        numeric = as.numeric(value)
        expected_column = switch(term, main = numeric, power = numeric^2, product = numeric * rows[["z"]])
        expected = as.numeric(cbind(1, expected_column) %*% beta)
        # Logical powers retain the existing model.matrix contrast warning;
        # numerical acceptance is the contract checked here.
        expect_equal(suppressWarnings(stats::predict(fit, rows)), expected, tolerance = 1e-10)
      }
      rejected = list(factor(c("100", "200")), ordered(c("100", "200")),
        c("100", "200"), c(100 + 1i, 200 + 2i))
      if (term != "power") rejected[[length(rejected) + 1L]] = c(TRUE, FALSE) else {
        rejected[[length(rejected) + 1L]] = as.Date("1970-01-01") + c(100, 200)
        rejected[[length(rejected) + 1L]] = as.POSIXct("1970-01-01", tz = "UTC") + c(100, 200)
      }
      for (value in rejected) {
        rows = data.frame(z = c(0.3, -0.2))
        rows[["x"]] = value
        expect_error(suppressWarnings(stats::predict(fit, rows)))
      }
      close(fit)
    }
  }
})
