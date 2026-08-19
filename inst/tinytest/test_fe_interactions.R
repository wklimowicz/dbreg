set.seed(2026)
n = 240
dat_fe = data.frame(
  subject = sample(c("math", "english", "science"), n, TRUE),
  cohort = sample(2019:2022, n, TRUE),
  fsm = sample(c("no", "yes"), n, TRUE),
  school = sample(letters[1:8], n, TRUE),
  x = rnorm(n)
)
dat_fe$y = 0.7 * dat_fe$x +
  as.numeric(interaction(dat_fe$subject, dat_fe$cohort)) / 4 + rnorm(n)
dat_fe$subject_cohort = interaction(dat_fe$subject, dat_fe$cohort, drop = TRUE)
dat_fe$cohort_fsm = interaction(dat_fe$cohort, dat_fe$fsm, drop = TRUE)

# Composite syntax is equivalent to a manually created collision-safe ID.
fit_composite = dbreg(
  y ~ x | subject ^ cohort,
  data = dat_fe,
  strategy = "demean"
)
fit_manual = dbreg(
  y ~ x | subject_cohort,
  data = dat_fe,
  strategy = "demean"
)
expect_equal(coef(fit_composite), coef(fit_manual), tolerance = 1e-8)
expect_equal(fit_composite$fe_labels, "subject^cohort")
expect_equal(fit_composite$n_fe_levels, fit_manual$n_fe_levels)

# Composite dimensions work in AP and compression paths.
fit_ap_composite = dbreg(
  y ~ x | subject ^ cohort + school,
  data = dat_fe[-1, ],
  strategy = "demean",
  fe_tol = 1e-7
)
fit_ap_manual = dbreg(
  y ~ x | subject_cohort + school,
  data = dat_fe[-1, ],
  strategy = "demean",
  fe_tol = 1e-7
)
expect_equal(coef(fit_ap_composite), coef(fit_ap_manual), tolerance = 1e-6)

fit_compress_composite = dbreg(
  y ~ x | cohort ^ fsm,
  data = dat_fe,
  strategy = "compress"
)
fit_compress_manual = dbreg(
  y ~ x | cohort_fsm,
  data = dat_fe,
  strategy = "compress"
)
expect_equal(coef(fit_compress_composite), coef(fit_compress_manual), tolerance = 1e-8)

# Tuple ranking distinguishes values that delimiter concatenation would collide.
collision_dat = data.frame(
  a = c("a_b", "a", "a_b", "a"),
  b = c("c", "b_c", "c", "b_c"),
  x = c(1, 2, 3, 4),
  y = c(2, 5, 5, 9)
)
collision_fit = dbreg(y ~ x | a ^ b, data = collision_dat, strategy = "demean")
expect_equal(unname(collision_fit$n_fe_levels), 2L)

sql_composite = dbreg(
  y ~ x | subject ^ cohort,
  data = dat_fe,
  strategy = "compress",
  sql_only = TRUE
)
expect_true(grepl("DENSE_RANK\\(\\) OVER \\(ORDER BY subject, cohort\\)", sql_composite))

# Compression retains tuple components for stored-data and newdata prediction.
expect_equal(
  predict(fit_compress_composite, newdata = dat_fe),
  predict(fit_compress_manual, newdata = dat_fe),
  tolerance = 1e-8
)
expect_silent(predict(fit_compress_composite))

# Generated identifiers never shadow similarly named source columns.
name_collision = transform(
  dat_fe,
  dbreg_fe_1__ = 1L,
  dbreg_internal_fe_1__ = 1L
)
fit_name_collision = dbreg(
  y ~ x | subject ^ cohort,
  data = name_collision,
  strategy = "demean"
)
expect_equal(coef(fit_name_collision), coef(fit_manual), tolerance = 1e-8)
