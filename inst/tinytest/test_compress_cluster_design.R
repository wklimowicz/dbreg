set.seed(840)
dat = data.frame(x = rnorm(400), z = rnorm(400), f = sample(1:8, 400, TRUE),
  category = sample(c("control", "treated group", "O'Brien"), 400, TRUE),
  cl = rep(1:20, each = 20), w = runif(400, 0.2, 2))
dat[["y"]] = dat[["x"]] * (dat[["category"]] == "treated group") + dat[["z"]] + rnorm(400)
dat[["category"]] = factor(dat[["category"]], levels = c("O'Brien", "control", "treated group"))
for (fml in list(y ~ x * category + z, y ~ x * category + z | f)) {
  for (weighted in c(FALSE, TRUE)) {
    fit = dbreg(fml, data = dat, strategy = "compress", vcov = ~cl,
      weights = if (weighted) "w" else NULL)
    ref = fixest::feols(fml, dat, vcov = ~cl, weights = if (weighted) ~w else NULL,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
    expect_equal(sort(names(coef(fit))), sort(names(coef(ref))))
    expect_equal(coef(fit), coef(ref)[names(coef(fit))], tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[names(coef(fit)), names(coef(fit))]), tolerance = 1e-8)
  }
}
