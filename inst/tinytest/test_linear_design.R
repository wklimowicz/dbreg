set.seed(860)
n = 360L
dat = data.frame(x = rnorm(n), z = rnorm(n), a = sample(letters[1:3], n, TRUE),
  b = sample(1:4, n, TRUE), f = sample(1:7, n, TRUE), cl = sample(1:18, n, TRUE),
  grp = sample(c("control", "treatment", "other"), n, TRUE), w = runif(n, 0.3, 3))
dat[["ab"]] = interaction(dat[["a"]], dat[["b"]], drop = TRUE)
dat[["y"]] = 1 + dat[["x"]] - 0.4 * dat[["x"]]^2 + 0.15 * dat[["x"]]^3 +
  0.3 * dat[["z"]] + (dat[["grp"]] == "treatment") + as.numeric(dat[["ab"]]) + rnorm(n)
reference_data = dat
reference_data[["grp"]] = factor(dat[["grp"]], levels = c("treatment", "control", "other"))
reference_names = function(x) sub("^I\\(I\\((.*)\\)\\)$", "I(\\1)", x)
for (strategy in c("moments", "compress", "demean")) {
  fml = if (strategy == "moments") y ~ x + I(x^2) + I(x^3) + grp + z else y ~ x + I(x^2) + I(x^3) + grp + z | a ^ b + f
  for (weighted in c(FALSE, TRUE)) for (kind in list("iid", "hc1", ~cl)) {
    fit = dbreg(fml, data = dat, weights = if (weighted) "w" else NULL, vcov = kind,
      strategy = strategy, factor_refs = c(grp = "treatment"))
    ref = fixest::feols(fml, reference_data, weights = if (weighted) ~w else NULL, vcov = kind,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
    ref_coef = setNames(coef(ref), reference_names(names(coef(ref))))
    index = match(names(coef(fit)), names(ref_coef))
    expect_equal(sort(names(coef(fit))), sort(names(ref_coef)))
    expect_equal(coef(fit), ref_coef[index], tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[index, index]), tolerance = 1e-7)
    expect_equal(predict(fit, dat), as.numeric(predict(ref, dat)), tolerance = 1e-6)
    if (strategy == "compress") expect_equal(predict(fit), predict(fit, fit[["data"]]), tolerance = 1e-9)
  }
}

# Mundlak forms means from expanded polynomial/factor columns, including weights.
for (weighted in c(FALSE, TRUE)) {
  fit = dbreg(y ~ x + I(x^2) + grp | a ^ b, data = dat, strategy = "mundlak",
    weights = if (weighted) "w" else NULL, factor_refs = c(grp = "treatment"))
  ref = fixest::feols(y ~ x + I(x^2) + grp | a ^ b, reference_data,
    weights = if (weighted) ~w else NULL, fixef.rm = "none")
  ref_coef = setNames(coef(ref), reference_names(names(coef(ref))))
  expect_equal(coef(fit), ref_coef[names(coef(fit))], tolerance = 1e-7)
  X = model.matrix(~ x + I(x^2) + grp, reference_data)[, -1, drop = FALSE]
  weights = if (weighted) dat[["w"]] else rep(1, nrow(dat))
  means = apply(X, 2, function(x) ave(seq_along(x), dat[["ab"]], FUN = function(i) weighted.mean(x[i], weights[i])))
  reference_cre = lm.wfit(cbind(1, X, means), dat[["y"]], weights)
  expect_equal(predict(fit, dat), as.numeric(reference_cre[["fitted.values"]]), tolerance = 1e-6)
}

# Names resembling implementation aliases and non-syntactic source names remain distinct.
dat[["dbreg_power_1__"]] = dat[["z"]]
names(dat)[names(dat) == "x"] = "first slope"
fit = dbreg(y ~ `first slope` + I(`first slope`^2) + dbreg_power_1__, data = dat, strategy = "compress")
ref = lm(y ~ `first slope` + I(`first slope`^2) + dbreg_power_1__, dat)
expect_equal(coef(fit), coef(ref), tolerance = 1e-8)
expect_equal(predict(fit, dat), as.numeric(predict(ref, dat)), tolerance = 1e-8)
expect_equal(predict(fit), predict(fit, fit[["data"]]), tolerance = 1e-8)

# Newdata subsets keep fitted levels and coding even if R's defaults change.
fit = dbreg(y ~ z + grp, data = dat, strategy = "compress", weights = "w", factor_refs = c(grp = "treatment"))
nd = dat[dat[["grp"]] == "other", c("z", "grp"), drop = FALSE]
baseline = predict(fit, nd)
old = options(contrasts = c("contr.sum", "contr.poly"))
expect_equal(predict(fit, nd), baseline)
options(old)
nd[["z"]][2] = NA_real_
prediction = predict(fit, nd)
expect_equal(length(prediction), nrow(nd))
expect_true(is.na(prediction[[2]]))
expect_equal(prediction[-2], baseline[-2])
expect_error(predict(fit, data.frame(z = 1, grp = "unseen")), "unseen")
expect_error(dbreg(y ~ z + grp, data = dat, factor_refs = c(grp = "absent")), "not found")
expect_error(dbreg(y ~ z, data = dat, factor_refs = c(grp = "control")), "factor_refs")
expect_error(dbreg(y ~ I(z^4), data = dat), "Only.*supported")

# No-intercept powers and factors retain the complete first factor coding.
fit = dbreg(y ~ 0 + grp + I(z^2), data = dat, strategy = "moments")
ref = lm(y ~ 0 + grp + I(z^2), dat)
expect_equal(coef(fit), coef(ref)[names(coef(fit))], tolerance = 1e-8)

# Tuple strings cannot collide through a delimiter occurring in a component.
collision = data.frame(a = rep(c("a|b", "a"), 20), b = rep(c("c", "b|c"), 20), x = rnorm(40), y = rnorm(40))
fit = dbreg(y ~ x | a ^ b, data = collision, strategy = "compress", vcov = ~a)
ref = fixest::feols(y ~ x | a ^ b, collision, vcov = ~a, fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
expect_equal(coef(fit), coef(ref), tolerance = 1e-8)
expect_equal(predict(fit, collision), as.numeric(predict(ref, collision)), tolerance = 1e-8)

# Numeric-coded categories, a polynomial-only design and stored-cell prediction.
dat[["category"]] = rep(1:3, length.out = nrow(dat))
ref_data = dat
ref_data[["category"]] = factor(dat[["category"]], levels = c(3, 1, 2))
for (strategy in c("moments", "compress")) {
  fit = dbreg(y ~ I(z^2) * category, data = dat, strategy = strategy,
    factor_refs = c(category = 3), weights = "w", vcov = "hc1")
  ref = lm(y ~ I(z^2) * category, ref_data, weights = w)
  expect_equal(coef(fit)[names(coef(ref))], coef(ref), tolerance = 1e-8)
  expect_equal(predict(fit, dat), as.numeric(predict(ref, ref_data)), tolerance = 1e-8)
  if (strategy == "compress") expect_equal(predict(fit), predict(fit, fit[["data"]]))
}

# FE names sharing a slope prefix must not appear in coef()/vcov() by default.
dat[["z_group"]] = dat[["f"]]
dat[["z_copy"]] = dat[["z"]]
fit = dbreg(y ~ z + z_copy | z_group, data = dat, strategy = "compress")
expect_equal(names(coef(fit)), "z")
expect_equal(rownames(vcov(fit)), "z")
expect_equal(names(coef(fit, collin = TRUE)), c("z", "z_copy"))

# Composite labels retain the user's spelling even when SQL aliases are needed.
dat[["first FE"]] = dat[["a"]]
fit = dbreg(y ~ z | `first FE` ^ b, data = dat, strategy = "demean")
expect_true(grepl("first FE", fit[["fe_labels"]], fixed = TRUE))
