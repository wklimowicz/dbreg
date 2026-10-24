# The same wide, rare-outcome teacher design runs locally and on native SQL.
dbglm_test_wide = function(conn, verbose = FALSE) {
  set.seed(2118)
  n = 44000L
  experience = rep(0:43, each = n / 44L)
  district = sample.int(30, n, TRUE)
  spine = pmin(12L, 1L + experience %/% 4L + rbinom(n, 1, .2))
  u = rnorm(n)
  v = .6 * u + .8 * rnorm(n)
  d = data.frame(id = seq_len(n), Experience = experience,
    outside = 10.8 + .2 * u + .05 * sin(district), teacher = 10.5 + .2 * v + spine / 30,
    spine = spine, district = district, AGE = pmin(64L, 21L + experience + sample(0:9, n, TRUE)),
    SEX = rbinom(n, 1, .72), FTPTWK = rbinom(n, 1, .25))
  eta = -4.2 + (.3 + .05 * sin(experience / 10)) * u - .2 * v +
    .1 * sin(spine) + .15 * cos(district) + .004 * d[["AGE"]]
  d[["y"]] = rbinom(n, 1, binomial("cloglog")[["linkinv"]](eta))
  formula = y ~ outside:factor(Experience) + teacher:factor(Experience) |
    spine + district + Experience + AGE + SEX + FTPTWK
  before = dbglm_validation_tables(conn)
  table = dbglm_validation_upload(conn, d)
  retained = basename(tempfile("glm_wide_sample_"))
  on.exit({dbglm_validation_drop(conn, retained); dbglm_validation_drop(conn, table)})
  if (verbose) message("[GLM wide] ", n, " rows; exits=", sum(d[["y"]]), "; 88 slopes, 6 FEs")
  # Constant shifts of either wage interaction lie in the absorbed Experience
  # span. Center the reference independently, keeping its slopes and likelihood
  # unchanged while avoiding almost-collinear FE/wage columns in fixest.
  reference_data = d
  for (name in c("outside", "teacher")) reference_data[[name]] = d[[name]] - mean(d[[name]])
  # At this size, explicit nuisance indicators give an independent fixest
  # oracle without iterative demeaning error in the clustered covariance.
  # Count all nuisance coefficients in both small-sample corrections.
  reference_formula = y ~ outside:factor(Experience) + teacher:factor(Experience) +
    factor(spine) + factor(district) + factor(Experience) + factor(AGE) + factor(SEX) + factor(FTPTWK)
  reference = dbglm_reference_feglm(reference_formula, reference_data,
    family = binomial("cloglog"), vcov = ~district, nthreads = 2)
  slopes = grep(":", names(coef(reference)), fixed = TRUE, value = TRUE)
  expect_equal(length(slopes), 88L)
  if (verbose) message("[GLM wide] independent reference complete; fitting database model")
  fit = dbglm(formula, conn = conn, table = table, vcov = ~district, ssc = list(K.fixef = "full"),
    row_id = "id", sample_table = retained, verbose = verbose)
  expect_equal(length(coef(fit)), 88L)
  dbglm_validation_compare(fit, reference, "wide rare teacher", terms = slopes)
  keys = sort(as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained))[["id"]]))
  expect_identical(keys, sort(as.numeric(d[["id"]][fixest::obs(reference)])))
  expect_identical(dbglm_validation_tables(conn), sort(c(before, table, retained)))
}
