dat_diag = data.frame(
  y = c(1, 2, NA, 4, 5, 6, 7),
  x = c(1, NA, 3, 4, 5, 6, 7),
  fe = c("a", "a", "b", NA, "c", "c", "d"),
  cl = c("a", "a", "b", "b", NA, "c", "d"),
  w = c(1, 1, 1, 1, 0, 1, 1)
)

fit_diag = dbreg(
  y ~ x | fe,
  data = dat_diag,
  weights = "w",
  cluster = ~cl,
  strategy = "demean",
  ssc = "nested",
  diagnostics = TRUE
)

g_diag = glance(fit_diag)
expect_equal(g_diag$nobs.input, 7)
expect_equal(g_diag$nobs, 3)
expect_equal(g_diag$nobs.dropped, 4)
expect_equal(g_diag$vcov.type, "cluster")
expect_equal(g_diag$cluster, "cl")
expect_equal(g_diag$ssc, "nested")
expect_equal(g_diag$n.fe, 1)
expect_false(is.na(g_diag$n.params))

audit = fit_diag$sample_diagnostics$detail
expect_equal(audit$n_missing[audit$variable == "y"], 1)
expect_equal(audit$n_missing[audit$variable == "x"], 1)
expect_equal(audit$n_missing[audit$variable == "fe"], 1)
expect_equal(audit$n_missing[audit$variable == "cl"], 1)
expect_equal(fit_diag$sample_diagnostics$n_nonpositive_weight, 1)

set.seed(42)
ap_dat = expand.grid(fe1 = 1:4, fe2 = 1:3)
ap_dat = ap_dat[-1, ]
ap_dat$x = rnorm(nrow(ap_dat))
ap_dat$y = ap_dat$x + rnorm(nrow(ap_dat))
ap_fit = dbreg(
  y ~ x | fe1 + fe2,
  data = ap_dat,
  strategy = "demean",
  fe_tol = 1e-7,
  fe_max_iter = 200
)
expect_true(ap_fit$ap_converged)
expect_equal(ap_fit$ap_tolerance, 1e-7)
expect_equal(ap_fit$ap_max_iter, 200L)
expect_true(glance(ap_fit)$ap.converged)

# AP models return an inspectable, structured SQL plan without executing AP.
ap_plan = dbreg(
  y ~ x | fe1 + fe2,
  data = ap_dat,
  strategy = "demean",
  sql_only = TRUE,
  fe_tol = 1e-7,
  fe_max_iter = 200
)
expect_true(inherits(ap_plan, "dbreg_query_plan"))
expect_true(grepl("SELECT", ap_plan$model_frame_sql, fixed = TRUE))
expect_equal(ap_plan$design_matrix$term, "x")
expect_equal(length(ap_plan$projection_sql), 2L)
expect_equal(length(ap_plan$convergence_sql), 2L)
expect_true(grepl("<ap_converged>", ap_plan$final_moments_sql, fixed = TRUE))

# Public coefficient and covariance methods expose exactly the same names.
cov_fit = dbreg(y ~ x | fe1 + fe2, data = ap_dat, strategy = "compress")
expect_equal(rownames(vcov(cov_fit)), names(coef(cov_fit)))
expect_equal(colnames(vcov(cov_fit)), names(coef(cov_fit)))
expect_equal(rownames(vcov(cov_fit, fe = TRUE)), names(coef(cov_fit, fe = TRUE)))
expect_equal(colnames(vcov(cov_fit, fe = TRUE)), names(coef(cov_fit, fe = TRUE)))

# Collinear terms can be restored as NA without changing the default methods.
collin_dat = transform(ap_dat, x_duplicate = x)
collin_fit = dbreg(
  y ~ x + x_duplicate | fe1,
  data = collin_dat,
  strategy = "demean"
)
expect_equal(length(collin_fit$collin.var), 1L)
expect_equal(length(coef(collin_fit)), 1L)
expect_equal(length(coef(collin_fit, collin = TRUE)), 2L)
expect_equal(sum(is.na(coef(collin_fit, collin = TRUE))), 1L)
expect_false(anyNA(tidy(collin_fit)$estimate))
expect_equal(sum(is.na(tidy(collin_fit, collin = TRUE)$estimate)), 1L)

# Mundlak means participate in the same collinearity handling as base terms.
mundlak_collin = transform(ap_dat, group_constant = fe1)
mundlak_collin_fit = dbreg(
  y ~ x + group_constant | fe1,
  data = mundlak_collin,
  strategy = "mundlak"
)
expect_true(length(mundlak_collin_fit$collin.var) >= 1L)

expect_error(
  dbreg(
    y ~ x + x_bar_fe1 | fe1,
    data = transform(ap_dat, x_bar_fe1 = rnorm(nrow(ap_dat))),
    strategy = "mundlak"
  ),
  "duplicate coefficient names"
)
