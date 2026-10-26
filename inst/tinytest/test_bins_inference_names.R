# Public-name restoration must rebind completed covariance coordinates.
local({
  set.seed(40404)
  d = data.frame(x = stats::runif(240))
  d[["raw control"]] = stats::rnorm(240)
  d[["group id"]] = rep(seq_len(8), each = 30)
  d[["cluster"]] = rep(seq_len(12), 20)
  d[["y"]] = 2 + d[["x"]] + 0.7 * d[["raw control"]] + sin(d[["group id"]]) + stats::rnorm(240)
  breaks = seq(0, 1, length.out = 5)
  d[["b"]] = cut(d[["x"]], breaks, include.lowest = TRUE, right = FALSE)
  for (with_fe in c(FALSE, TRUE)) {
    formula = if (with_fe) y ~ x + `raw control` | `group id` else y ~ x + `raw control`
    reference_formula = if (with_fe) y ~ b + `raw control` + factor(`group id`) else y ~ b + `raw control`
    reference = stats::lm(reference_formula, data = d)
    X = stats::model.matrix(reference)
    inverse = solve(crossprod(X))
    for (kind in list("iid", ~`cluster`)) {
      result = dbreg::dbbinsreg(formula, data = d, nbins = 4, binspos = breaks,
        points = c(0, 0), strategy = "compress", vcov = kind, plot = FALSE)
      model = result[["model"]]
      full = stats::vcov(model, fe = TRUE)
      expect_identical(rownames(full), names(stats::coef(model, fe = TRUE)))
      expect_identical(colnames(full), rownames(model[["coeftable"]]))
      expect_identical(model[["inference_result"]][["coefficient_names"]], rownames(full))
      expect_identical(rownames(stats::vcov(model)), names(stats::coef(model)))
      expect_identical(model[["inference_result"]][["status"]], "available")
      expect_identical(model[["inference_result"]][["reason"]], NULL)
      expect_identical(model[["fitted_state"]][["inference"]], NULL)
      expected = stats::vcov(reference)
      if (inherits(kind, "formula")) {
        scores = rowsum(X * as.numeric(stats::residuals(reference)), d[["cluster"]])
        groups = nrow(scores)
        expected = inverse %*% crossprod(scores) %*% inverse * groups / (groups - 1) * (nrow(X) - 1) / (nrow(X) - ncol(X))
        expect_identical(model[["inference_result"]][["clusters"]], list("cluster"))
        expect_identical(model[["cluster_var"]], "cluster")
        expect_equal(attr(full, "n_clusters"), groups)
      }
      # Known generated basis order is intercept, bins 2:4, control, FE contrasts.
      expect_equal(as.numeric(stats::coef(model, fe = TRUE)), as.numeric(stats::coef(reference)), tolerance = 1e-8)
      expect_equal(as.numeric(full), as.numeric(expected), tolerance = 1e-8)
      expect_equal(as.numeric(model[["coeftable"]][, "std.error"]), as.numeric(sqrt(diag(expected))), tolerance = 1e-8)
    }
  }
})

# Pre-finalization producer payloads retain their original path and shape.
local({
  V = matrix(4, 1, 1, dimnames = list("safe", "safe"))
  attr(V, "df.t") = 7
  legacy = list(coeftable = matrix(c(2, 2, 1, 0.4), 1, 4,
    dimnames = list("safe", c("estimate", "std.error", "statistic", "p.value"))), vcov = V)
  map = c(safe = "raw name")
  restored = dbreg:::restore_result_names(legacy, list(name_map = map, source_map = map))
  expect_identical(rownames(restored[["vcov"]]), "`raw name`")
  expect_identical(rownames(restored[["coeftable"]]), "`raw name`")
  expect_equal(as.numeric(restored[["vcov"]]), 4)
  expect_equal(attr(restored[["vcov"]], "df.t"), 7)
  expect_identical(restored[["inference_result"]], NULL)
})
