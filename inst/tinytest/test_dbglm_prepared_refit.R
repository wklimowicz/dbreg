if (requireNamespace("duckdb", quietly = TRUE)) local({
  set.seed(61)
  d = expand.grid(g = 1:5, h = 1:4, x = -1:1, rep = 1:5)
  d$id = seq_len(nrow(d))
  if (requireNamespace("bit64", quietly = TRUE)) {
    d$id = bit64::as.integer64("9007199254740992") + d$id
  }
  d$cl = d$rep
  d$o = d$h / 15
  d$y = rpois(nrow(d), exp(.4 + d$x / 5 + d$g / 12 + d$o))
  d$x[5L] = NA_real_
  conn = DBI::dbConnect(duckdb::duckdb(bigint = if (requireNamespace("bit64", quietly = TRUE)) "integer64" else "numeric"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "design_input", d)
  baseline = sort(DBI::dbListTables(conn))
  p = dbglm_prepare(y ~ x | g + h, conn, "design_input", family = poisson(),
    row_id = "id", offset = ~o, vcov = ~cl)
  s = data.frame(id = d$id, w = runif(nrow(d), .1, 5), eta = 0)
  d$w = s$w
  d$eta = s$eta
  DBI::dbWriteTable(conn, "state", s[rev(seq_len(nrow(s))), ])
  fresh = function(control = dbglm_control(), warm = FALSE) {
    dbglm(y ~ x | g + h, data = d, family = poisson(), offset = ~o, weights = ~w, vcov = ~cl,
      etastart = if (warm) ~eta, strategy = "moments", control = control, retain = "inference")
  }
  ref = fresh()
  a = dbglm_refit(p, "state", etastart = FALSE)
  expect_true(a$prepared_reused)
  expect_equal(coef(a), coef(ref), tolerance = 1e-10)
  expect_equal(vcov(a), vcov(ref), tolerance = 1e-10)
  expect_equal(logLik(a, type = "weighted"), logLik(ref, type = "weighted"), tolerance = 1e-10)
  saved = predict(a, type = "link")
  expect_equal(as.data.frame(saved)$.link, predict(ref, type = "link"), tolerance = 1e-10)
  close(saved)
  close(a)
  close(ref)

  DBI::dbExecute(conn, "UPDATE design_input SET x=999,y=1")
  d$w = s$w = runif(nrow(d), 1, 12)
  d$eta = s$eta = d$o + 0.15 * ifelse(is.na(d$x), 0, d$x)
  DBI::dbWriteTable(conn, "state", s, overwrite = TRUE)
  control = dbglm_control(compute_inference = FALSE, tolerate_max_iter = TRUE, glm_max_iter = 1L)
  ref = fresh(control, TRUE)
  a = dbglm_refit(p, "state", control = control, compact = TRUE)
  expect_true(a$prepared_reused)
  expect_true(a$prepared_compact)
  expect_equal(coef(a), coef(ref), tolerance = 1e-10)
  expect_equal(a$deviance, ref$deviance, tolerance = 1e-10)
  augmented = generics::augment(a)
  expect_identical(names(augmented), c("id", "w", ".fitted", ".resid"))
  expect_equal(augmented$w, s$w[-5L], tolerance = 1e-14)
  export = predict(a, type = "link", prediction_table = "em_eta")
  exported = as.data.frame(export)
  expect_identical(as.character(exported$id), as.character(d$id[-5L]))
  expect_equal(exported$prediction, predict(ref, type = "link"), tolerance = 1e-10)
  close(export)
  close(a)
  expect_true(DBI::dbExistsTable(conn, "em_eta"))
  DBI::dbRemoveTable(conn, "em_eta")
  close(ref)
  expect_error(dbglm_refit(p, "state", compact = TRUE), "compute_inference")
  expect_error(dbglm_refit(p, "state", vcov = ~h), "Cluster")

  DBI::dbWriteTable(conn, "bad", s[-1L, ], overwrite = TRUE)
  expect_error(dbglm_refit(p, "bad"), "keys")
  DBI::dbWriteTable(conn, "bad", s[c(1L, seq_len(nrow(s) - 1L)), ], overwrite = TRUE)
  expect_error(dbglm_refit(p, "bad"), "unique")
  invalid = s
  invalid$w[1L] = -1
  DBI::dbWriteTable(conn, "bad", invalid, overwrite = TRUE)
  expect_error(dbglm_refit(p, "bad"), "nonnegative")
  s$w[seq_len(15L)] = d$w[seq_len(15L)] = 0
  s$eta[40L] = d$eta[40L] = Inf
  DBI::dbWriteTable(conn, "state", s, overwrite = TRUE)
  a = dbglm_refit(p, "state")
  ref = fresh(warm = TRUE)
  expect_false(a$prepared_reused)
  expect_equal(coef(a), coef(ref), tolerance = 1e-10)
  expect_equal(nobs(a), nobs(ref))
  close(p)
  expect_error(dbglm_refit(p, "state"), "no live database session")
  DBI::dbRemoveTable(conn, "state")
  expect_equal(vcov(a, vcov = "hc1"), vcov(ref, vcov = "hc1"), tolerance = 1e-10)
  close(a)
  close(ref)
  DBI::dbRemoveTable(conn, "bad")
  expect_identical(sort(DBI::dbListTables(conn)), baseline)
})

if (requireNamespace("duckdb", quietly = TRUE)) local({
  set.seed(155)
  d = data.frame(id = 1:96, x = rnorm(96), g = rep(1:6, length.out = 96))
  d$w = runif(nrow(d), .5, 3)
  d$y = .3 * d$x + d$g / 4 + rnorm(nrow(d))
  d$cl = rep(1:8, length.out = nrow(d))
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "design_input", d)
  DBI::dbWriteTable(conn, "state", d[c("id", "w")])
  for (formula in list(y ~ x, y ~ x | g)) {
    p = dbglm_prepare(formula, conn, "design_input", row_id = "id")
    fit = dbglm_refit(p, "state", etastart = FALSE)
    close(p)
    reference = lm(if (length(fit$fe)) y ~ x + factor(g) else y ~ x, data = d, weights = w)
    expect_identical(fit$prepared_reused, length(fit$fe) > 0L)
    expect_equal(coef(fit), coef(reference)[names(coef(fit))], tolerance = 1e-10)
    reference_vcov = vcov(reference)[names(coef(fit)), names(coef(fit))]
    # The explicit IID SSC policy also applies to absorbed Gaussian fits.
    reference_vcov = reference_vcov * (nrow(d) - 1) / df.residual(reference)
    expect_equal(as.numeric(vcov(fit)), as.numeric(reference_vcov), tolerance = 1e-10)
    actual = fitted(fit)
    expect_equal(as.data.frame(actual)$.fitted, unname(fitted(reference)), tolerance = 1e-10)
    close(actual)
    close(fit)
  }
  DBI::dbExecute(conn, "UPDATE design_input SET cl=NULL WHERE id=1")
  expect_error(dbglm_prepare(y ~ x | g, conn, "design_input", row_id = "id", vcov = ~cl), "cluster keys")
  expect_error(dbglm_prepare(y ~ x | g, conn, "design_input", row_id = "id", weights = "x"), "immutable model")
  expect_identical(sort(DBI::dbListTables(conn)), c("design_input", "state"))
})

# Centered Gaussian refits must transform the supplied predictor too.
local({
  set.seed(4)
  d_warm = data.frame(id = 1:120, x = rnorm(120), g = rep(1:6, each = 20))
  d_warm$y = 1e8 + 2 + d_warm$x + d_warm$g / 10 + rnorm(120)
  d_warm$eta = fitted(lm(y ~ x + factor(g), d_warm))
  d_warm$w = 1
conn_prep = DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(conn_prep, "review_design", d_warm)
DBI::dbWriteTable(conn_prep, "review_state", d_warm[c("id", "w", "eta")])
prep = dbglm_prepare(y ~ x | g, conn_prep, "review_design", row_id = "id")
refit_cold = dbglm_refit(prep, "review_state", etastart = FALSE)
refit_warm = dbglm_refit(prep, "review_state")
expect_true(refit_warm$converged)
expect_equal(unname(coef(refit_warm)), unname(coef(refit_cold)), tolerance = 1e-8)
close(refit_cold)
close(refit_warm)
close(prep)
DBI::dbDisconnect(conn_prep, shutdown = TRUE)

})
