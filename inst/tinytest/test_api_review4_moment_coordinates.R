local({
  d = expand.grid(x0 = -2:2, z = -1:1, rep = 1:4)
  d$x = 1e12 + d$x0/8
  d$y = .5 + .4*d$x0/8 + .3*d$z + ifelse(d$rep %% 2, -.1, .1)
  nd = data.frame(id = 1:2, x = 1e12 + c(.375, -.125), z = c(.25, -.75))
  expected = .5 + .4*c(.375, -.125) + .3*nd$z
  for (strategy in c("moments", "compress")) {
    fit = dbreg(y ~ x + z, data = d, strategy = strategy)
    expect_equal(unname(predict(fit, nd)), expected, tolerance = 1e-10)
    result = predict(fit, nd, row_id = "id")
    expect_equal(as.data.frame(result)$prediction, expected, tolerance = 1e-10)
    close(result)
  }
})
