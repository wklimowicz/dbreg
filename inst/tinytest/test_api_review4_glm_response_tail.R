# The response mean is EXP(eta), including representable means far below
# machine epsilon. Compare ratios so absolute tolerances cannot hide a floor.
local({
  d = data.frame(id = 1:40, y = rep(c(1, 2), 20), o = rep(c(-600, 0), each = 20))
  alpha = log(sum(d$y) / sum(exp(d$o)))
  expected = exp(alpha + d$o)
  for (strategy in c("moments", "compress")) {
    fit = dbglm(y ~ 1, data = d, offset = "o", family = poisson(),
      strategy = strategy, retain = "inference")
    expect_true(fit$converged && fit$finite_mle)
    expect_equal(unname(predict(fit, d, type = "response")) / expected,
      rep(1, nrow(d)), tolerance = 1e-10)
    expect_equal(fitted(fit) / expected, rep(1, nrow(d)), tolerance = 1e-10)
    result = predict(fit, d, type = "response", row_id = "id")
    expect_equal(as.data.frame(result)$prediction / expected,
      rep(1, nrow(d)), tolerance = 1e-10)
    close(result)
    # Keep the log-link overflow handling and ordinary underflow behavior.
    extreme = data.frame(id = 1:2, o = c(710, -1000))
    result = predict(fit, extreme, type = "response", row_id = "id")
    expect_identical(as.data.frame(result)$prediction, c(Inf, 0))
    close(result)
    close(fit)
  }
})
