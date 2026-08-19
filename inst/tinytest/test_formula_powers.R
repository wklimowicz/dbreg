set.seed(91)
n = 300
dat_power = data.frame(
  x = runif(n, -2, 2),
  grp = sample(c("control", "treatment", "other"), n, TRUE),
  fe = sample(letters[1:12], n, TRUE)
)
dat_power$y = 1 + 0.8 * dat_power$x - 0.5 * dat_power$x^2 +
  0.2 * dat_power$x^3 + ifelse(dat_power$grp == "treatment", 1, 0) + rnorm(n)

lm_power = lm(y ~ x + I(x^2) + I(x^3), data = dat_power)
fit_power_moments = dbreg(
  y ~ x + I(x^2) + I(x^3),
  data = dat_power,
  strategy = "moments"
)
expect_equal(coef(fit_power_moments), coef(lm_power), tolerance = 1e-8)
expect_true(all(c("I(x^2)", "I(x^3)") %in% tidy(fit_power_moments)$term))

fit_power_demean = dbreg(
  y ~ x + I(x^2) + I(x^3) | fe,
  data = dat_power,
  strategy = "demean"
)
fixest_power = fixest::feols(y ~ x + I(x^2) + I(x^3) | fe, data = dat_power, fixef.rm = "none")
expect_equal(unname(coef(fit_power_demean)), unname(coef(fixest_power)), tolerance = 1e-7)
expect_equal(rownames(vcov(fit_power_demean)), names(coef(fit_power_demean)))
expect_equal(colnames(vcov(fit_power_demean)), names(coef(fit_power_demean)))

fit_power_compress = dbreg(
  y ~ x + I(x^2),
  data = dat_power,
  strategy = "compress"
)
expect_equal(
  unname(coef(fit_power_compress)),
  unname(coef(lm(y ~ x + I(x^2), dat_power))),
  tolerance = 1e-8
)

fit_power_interaction = dbreg(
  y ~ grp:I(x^2) + x,
  data = dat_power,
  strategy = "moments"
)
lm_power_interaction = lm(y ~ grp:I(x^2) + x, data = dat_power)
expect_equal(
  unname(coef(fit_power_interaction)),
  unname(coef(lm_power_interaction)),
  tolerance = 1e-8
)

# Explicit references work for character and numeric-coded categories.
fit_ref = dbreg(
  y ~ grp + x,
  data = dat_power,
  factor_refs = c(grp = "treatment"),
  strategy = "moments"
)
lm_ref = lm(y ~ relevel(factor(grp), "treatment") + x, data = dat_power)
expect_equal(unname(coef(fit_ref)), unname(coef(lm_ref)), tolerance = 1e-8)
expect_false(any(grepl("treatment", names(coef(fit_ref)))))

dat_power$group_num = match(dat_power$grp, c("control", "treatment", "other"))
fit_numeric_ref = dbreg(
  y ~ group_num + x,
  data = dat_power,
  factor_refs = c(group_num = 2),
  strategy = "moments"
)
expect_equal(length(coef(fit_numeric_ref)), 4L)
expect_false(any(grepl("group_num2", names(coef(fit_numeric_ref)))))

# Transformed names remain aligned in covariance, collinearity, and AP plans.
dat_power$fe2 = sample(LETTERS[1:7], nrow(dat_power), TRUE)
dat_power$x_duplicate = dat_power$x
fit_power_collin = dbreg(
  y ~ x + I(x^2) + x_duplicate | fe,
  data = dat_power,
  strategy = "compress"
)
expect_equal(
  names(coef(fit_power_collin, collin = TRUE)),
  c("x", "I(x^2)", "x_duplicate")
)
expect_equal(sum(is.na(coef(fit_power_collin, collin = TRUE))), 1L)
expect_equal(rownames(vcov(fit_power_collin)), names(coef(fit_power_collin)))
expect_equal(colnames(vcov(fit_power_collin)), names(coef(fit_power_collin)))

power_ap_plan = dbreg(
  y ~ x + I(x^2) | fe + fe2,
  data = dat_power[-1, ],
  strategy = "demean",
  sql_only = TRUE
)
expect_true(inherits(power_ap_plan, "dbreg_query_plan"))
expect_equal(power_ap_plan$design_matrix$term, c("x", "I(x^2)"))
expect_true(grepl("POWER(x, 2)", power_ap_plan$model_frame_sql, fixed = TRUE))

expect_error(dbreg(y ~ I(x^4), data = dat_power), "Only .* are supported")
expect_error(dbreg(y ~ x^2, data = dat_power), "Unsupported regressor expression")
expect_error(dbreg(y ~ log(x), data = dat_power), "Unsupported regressor expression")
expect_error(
  dbreg(y ~ grp, data = dat_power, factor_refs = c(grp = "missing")),
  "was not found"
)

# Prediction reuses the exact fitted expansion for every strategy.
for (strategy in c("moments", "compress")) {
  fit = dbreg(y ~ x + I(x^2) + I(x^3), data = dat_power, strategy = strategy)
  expect_equal(
    predict(fit, newdata = dat_power),
    as.numeric(predict(lm_power, newdata = dat_power)),
    tolerance = 1e-6
  )
}

fit_power_demean = dbreg(
  y ~ x + I(x^2) + I(x^3) | fe,
  data = dat_power,
  strategy = "demean"
)
expect_equal(
  predict(fit_power_demean, newdata = dat_power),
  as.numeric(predict(fixest_power, newdata = dat_power)),
  tolerance = 1e-5
)

# A source name that resembles the old power alias remains an ordinary term.
dat_collision = transform(dat_power, x__pow2 = rev(x))
fit_collision = dbreg(
  y ~ x__pow2 + I(x^2),
  data = dat_collision,
  strategy = "moments"
)
lm_collision = lm(y ~ x__pow2 + I(x^2), data = dat_collision)
expect_equal(names(coef(fit_collision)), names(coef(lm_collision)))
expect_equal(
  predict(fit_collision, newdata = dat_collision),
  as.numeric(predict(lm_collision, newdata = dat_collision)),
  tolerance = 1e-6
)

# Opaque design aliases avoid user regressors and weight columns.
dat_internal_name = transform(
  dat_power,
  dbreg_internal_x_0001__ = rev(x),
  dbreg_internal_x_0002__ = seq_len(nrow(dat_power)) / 10
)
fit_internal_name = dbreg(
  y ~ x + dbreg_internal_x_0001__,
  data = dat_internal_name,
  weights = "dbreg_internal_x_0002__",
  strategy = "moments"
)
expect_equal(
  names(coef(fit_internal_name)),
  c("(Intercept)", "x", "dbreg_internal_x_0001__")
)
