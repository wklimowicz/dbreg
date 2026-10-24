# The same wide, rare-outcome synthetic design runs locally and on native SQL.
dbglm_test_wide = function(conn, verbose = FALSE) {
  set.seed(2118)
  n = 44000L
  level_c = rep(0:43, each = n / 44L)
  fe_b = sample.int(30, n, TRUE)
  fe_a = pmin(12L, 1L + level_c %/% 4L + rbinom(n, 1, .2))
  u = rnorm(n)
  v = .6 * u + .8 * rnorm(n)
  d = data.frame(id = seq_len(n), fe_c = level_c,
    x1 = 10.8 + .2 * u + .05 * sin(fe_b), x2 = 10.5 + .2 * v + fe_a / 30,
    fe_a = fe_a, fe_b = fe_b, fe_d = pmin(64L, 21L + level_c + sample(0:9, n, TRUE)),
    fe_e = rbinom(n, 1, .72), fe_f = rbinom(n, 1, .25))
  eta = -4.2 + (.3 + .05 * sin(level_c / 10)) * u - .2 * v +
    .1 * sin(fe_a) + .15 * cos(fe_b) + .004 * d[["fe_d"]]
  d[["y"]] = rbinom(n, 1, binomial("cloglog")[["linkinv"]](eta))
  formula = y ~ x1:factor(fe_c) + x2:factor(fe_c) |
    fe_a + fe_b + fe_c + fe_d + fe_e + fe_f
  before = dbglm_validation_tables(conn)
  table = dbglm_validation_upload(conn, d)
  retained = basename(tempfile("glm_wide_sample_"))
  on.exit({dbglm_validation_drop(conn, retained); dbglm_validation_drop(conn, table)})
  if (verbose) message("[GLM wide] ", n, " rows; exits=", sum(d[["y"]]), "; 88 slopes, 6 FEs")
  # Constant shifts of either regressor interaction lie in the absorbed fe_c
  # span. Center the reference independently, keeping its slopes and likelihood
  # unchanged while avoiding almost-collinear FE/regressor columns in fixest.
  reference_data = d
  for (name in c("x1", "x2")) reference_data[[name]] = d[[name]] - mean(d[[name]])
  # At this size, explicit nuisance indicators give an independent fixest
  # oracle without iterative demeaning error in the clustered covariance.
  # Count all nuisance coefficients in both small-sample corrections.
  reference_formula = y ~ x1:factor(fe_c) + x2:factor(fe_c) +
    factor(fe_a) + factor(fe_b) + factor(fe_c) + factor(fe_d) + factor(fe_e) + factor(fe_f)
  reference = dbglm_reference_feglm(reference_formula, reference_data,
    family = binomial("cloglog"), vcov = ~fe_b, nthreads = 2)
  slopes = grep(":", names(coef(reference)), fixed = TRUE, value = TRUE)
  expect_equal(length(slopes), 88L)
  if (verbose) message("[GLM wide] independent reference complete; fitting database model")
  fit = dbglm(formula, conn = conn, table = table, vcov = ~fe_b, ssc = list(K.fixef = "full"),
    row_id = "id", sample_table = retained, verbose = verbose)
  expect_equal(length(coef(fit)), 88L)
  dbglm_validation_compare(fit, reference, "wide rare outcome", terms = slopes)
  keys = sort(as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained))[["id"]]))
  expect_identical(keys, sort(as.numeric(d[["id"]][fixest::obs(reference)])))
  expect_identical(dbglm_validation_tables(conn), sort(c(before, table, retained)))
}
