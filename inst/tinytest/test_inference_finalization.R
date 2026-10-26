# Completed inference records preserve established numerical/reference policies.
local({
  set.seed(12012)
  d = data.frame(x = stats::rnorm(160), z = stats::rnorm(160), g = rep(seq_len(8), each = 20))
  d[["a"]] = 0.8 * d[["z"]] + stats::rnorm(160)
  d[["y"]] = 1 + 0.4 * d[["x"]] + 1.2 * d[["a"]] + stats::rnorm(160)
  d[["b"]] = stats::rbinom(160, 1, stats::plogis(0.2 + 0.3 * d[["x"]]))
  d[["p"]] = stats::rpois(160, exp(0.1 + 0.2 * d[["x"]]))
  d[["g^raw"]] = d[["g"]]
  for (kind in list("iid", "hc1", ~`g^raw`)) {
    fit = dbreg::dbreg(y ~ x, data = d, vcov = kind)
    reference = fixest::feols(y ~ x, data = d, vcov = kind, notes = FALSE)
    expect_equal(c(stats::vcov(fit)), c(stats::vcov(reference)), tolerance = 1e-8)
    expect_equal(fit[["coeftable"]][, "estimate"], stats::coef(reference), tolerance = 1e-8)
    expect_identical(fit[["inference_result"]][["status"]], "available")
    expect_identical(fit[["fitted_state"]][["inference"]], NULL)
    if (inherits(kind, "formula")) expect_identical(fit[["inference_result"]][["clusters"]], list("g^raw"))
    iv = dbreg::dbivreg(y ~ x | a ~ z, data = d, vcov = kind)
    iv_reference = fixest::feols(y ~ x | a ~ z, data = d, vcov = kind, notes = FALSE)
    reference_coefficients = stats::coef(iv_reference)
    names(reference_coefficients)[names(reference_coefficients) == "fit_a"] = "a"
    reference_covariance = stats::vcov(iv_reference)
    dimnames(reference_covariance) = rep(list(names(reference_coefficients)), 2)
    coefficient_order = names(stats::coef(iv))
    expect_equal(c(stats::vcov(iv)), c(reference_covariance[coefficient_order, coefficient_order]), tolerance = 1e-7)
    expect_equal(stats::coef(iv), reference_coefficients[coefficient_order], tolerance = 1e-8)
    expect_identical(iv[["inference_result"]][["status"]], "available")
  }
  for (family in c("gaussian", "binomial", "poisson")) {
    fml = if (family == "gaussian") y ~ x else if (family == "binomial") b ~ x else p ~ x
    fit = dbreg::dbglm(fml, data = d, family = family)
    reference = fixest::feglm(fml, data = d, family = family, glm.tol = 1e-12, notes = FALSE)
    expect_equal(c(stats::vcov(fit)), c(stats::vcov(reference)), tolerance = 1e-5)
    expect_identical(fit[["inference_result"]][["status"]], "available")
    expect_equal(fit[["df_t"]], attr(stats::vcov(fit), "df.t"))
    expect_equal(fit[["n_params_vcov"]], attr(stats::vcov(fit), "df.K"))
    expect_identical(fit[["fitted_state"]], NULL)
    if (family == "gaussian") {
      expect_equal(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(reference)), tolerance = 1e-7)
      expect_equal(fit[["dispersion_scaled"]], fit[["deviance"]] / stats::df.residual(fit), tolerance = 1e-10)
    } else {
      expect_equal(fit[["coeftable"]][, "p.value"], 2 * stats::pnorm(-abs(fit[["coeftable"]][, "statistic"])))
    }
  }

  joint = dbreg::dbreg(y ~ x | g, data = d, strategy = "compress")
  reference = stats::lm(y ~ x + factor(g), data = d)
  expect_equal(as.numeric(stats::coef(joint, fe = TRUE)), as.numeric(stats::coef(reference)), tolerance = 1e-8)
  expect_equal(c(stats::vcov(joint, fe = TRUE)), c(stats::vcov(reference)), tolerance = 1e-8)
  expect_equal(nrow(joint[["coeftable"]]), nrow(stats::vcov(joint, fe = TRUE)))
  expect_equal(c(stats::vcov(joint)), as.numeric(stats::vcov(reference)["x", "x", drop = FALSE]), tolerance = 1e-8)
  for (family in c("ols", "gaussian", "binomial")) {
    empty = if (family == "ols") dbreg::dbreg(y ~ 1 | g, data = d, strategy = "demean") else
      dbreg::dbglm(if (family == "gaussian") y ~ 1 | g else b ~ 1 | g, data = d, family = family)
    expect_identical(dim(stats::vcov(empty)), c(0L, 0L))
    expect_identical(empty[["inference_result"]][["status"]], "available")
  }

  retained = dbreg::dbreg(y ~ x, data = d, retain = "inference", vcov = ~g)
  V = stats::vcov(retained)
  ci = stats::confint(retained)
  close(retained)
  expect_equal(stats::vcov(retained), V)
  expect_equal(stats::confint(retained), ci)
  expect_error(summary(retained, vcov = "hc1"), "released")

  portable = dbreg::dbglm(y ~ x | g, data = d, family = "gaussian", retain = "prediction")
  expect_identical(portable[["fitted_state"]][["inference"]], NULL)
  expect_true(!is.null(portable[["fitted_state"]][["prediction"]]))
  expect_error(summary(portable, vcov = "hc1"), "retain")
  close(portable)

  # Independent constructor/installer seam: no fit, Matrix class or fake recipe.
  V = Matrix::Diagonal(x = c(4, 9))
  dimnames(V) = rep(list(c("a", "b")), 2)
  coefficients = stats::setNames(c(2, -3), c("a", "b"))
  record = dbreg:::new_inference_result(V, "iid", c("a", "b"), list(), NULL, "full", 2, 7, "t")
  object = structure(list(), class = "dbreg")
  object = dbreg:::install_inference_result(object, record, coefficients)
  expect_identical(object[["fitted_state"]], NULL)
  expect_identical(class(stats::vcov(object)), class(V))
  expect_equal(object[["coeftable"]][, "std.error"], stats::setNames(c(2, 3), c("a", "b")))
  expect_equal(object[["coeftable"]][, "p.value"], 2 * stats::pt(-abs(c(1, -1)), 7), check.attributes = FALSE)
})

local({
  d = data.frame(x = rep(c(-1, 1), each = 12), y = rep(c(0, 1), each = 12), g = rep(seq_len(4), 6))
  for (kind in list("iid", ~g)) {
    fit = suppressWarnings(dbreg::dbglm(y ~ x, data = d, family = stats::binomial(), vcov = kind,
      separation = "prediction", ssc = list(t.df = 5)))
    expect_identical(fit[["inference_result"]][["status"]], "unavailable")
    expect_identical(fit[["inference_result"]][["reason"]], "nonfinite_mle")
    expect_true(all(is.na(stats::vcov(fit))))
    expect_true(all(is.na(fit[["coeftable"]][, c("std.error", "statistic", "p.value")])))
    expect_equal(fit[["df_t"]], 5)
    expect_equal(attr(stats::vcov(fit), "df.t"), 5)
    expect_identical(fit[["n_clusters"]], NULL)
    if (inherits(kind, "formula")) expect_identical(fit[["inference_result"]][["clusters"]], list("g"))
    expect_true(max(abs(stats::predict(fit, d, type = "response") - d[["y"]])) < 1e-7)
    expect_true(is.finite(as.numeric(stats::logLik(fit))))
    expect_true(is.na(attr(stats::logLik(fit), "df")))
    expect_true(is.na(stats::AIC(fit)))
    expect_true(is.na(stats::BIC(fit)))
    expect_identical(fit[["fitted_state"]], NULL)
  }
  expect_warning(dbreg::dbglm(y ~ x, data = d, family = stats::binomial(),
    vcov = ~g, separation = "prediction"), "boundary probabilities")
  default_df = suppressWarnings(dbreg::dbglm(y ~ x, data = d, family = stats::binomial(),
    vcov = ~g, separation = "prediction"))
  expect_equal(default_df[["df_t"]], nrow(d) - default_df[["n_params_vcov"]])
  expect_error(suppressWarnings(dbreg::dbglm(y ~ x, data = d, family = stats::binomial(),
    separation = "prediction", retain = "inference")), "finite GLM fit")
})
