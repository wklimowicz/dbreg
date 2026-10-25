library(dbreg)
library(fixest)
library(DBI)

set.seed(1927)
n_glm <- 480L
glm_dat <- data.frame(
  x = rnorm(n_glm), z = rnorm(n_glm),
  grp = factor(sample(c("a", "b", "c"), n_glm, TRUE)),
  sex = factor(sample(c("f", "m"), n_glm, TRUE)),
  history = factor(sample(c("no", "yes"), n_glm, TRUE)),
  cluster = rep(1:24, each = 20),
  a = sample(1:3, n_glm, TRUE), o = runif(n_glm, -0.4, 0.4)
)
close_glm <- function(actual, reference, atol, rtol) {
  expect_equal(dim(actual), dim(reference))
  expect_true(all(is.finite(actual)))
  expect_true(all(abs(actual - reference) <= atol + rtol * abs(reference)),
    info = paste(glm_case, paste(deparse(substitute(actual)), collapse = ""),
                 "max absolute error", signif(max(abs(actual - reference)), 8)))
}
glm_case <- "initial"
source(system.file('tinytest/helper_glm_reference.R',package='dbreg'),local=TRUE)
reference_feglm <- dbglm_reference_feglm
for (link in c("cloglog", "logit")) {
  family_glm <- binomial(link)
  eta_glm <- -1 + 0.35 * glm_dat$x - 0.15 * glm_dat$z + glm_dat$o
  glm_dat$y <- rbinom(n_glm, 1, family_glm$linkinv(eta_glm))
  for (weighted in c(FALSE, TRUE)) {
    for (kind in c("iid", "hc1", "cluster")) {
      glm_case <- paste(link, "weighted", weighted, kind)
      covariance <- if (kind == "cluster") ~cluster else kind
      covariance_ref <- if (kind == "hc1") "hetero" else covariance
      fit_glm <- dbglm(y ~ x + z, data = glm_dat, family = family_glm,
                      weights = if (weighted) "a" else NULL,
                      offset = "o", vcov = covariance)
      ref_glm <- reference_feglm(y ~ x + z, data = glm_dat, family = family_glm,
                      weights = if (weighted) ~a else NULL,
                      offset = ~o, vcov = covariance_ref)
      expect_identical(names(coef(fit_glm)), names(coef(ref_glm)))
      close_glm(coef(fit_glm), coef(ref_glm), 1e-7, 1e-6)
      close_glm(sqrt(diag(vcov(fit_glm))), sqrt(diag(vcov(ref_glm))), 1e-8, 1e-6)
      close_glm(vcov(fit_glm), vcov(ref_glm), 1e-8, 1e-5)
      expect_equal(as.numeric(logLik(fit_glm)), as.numeric(logLik(ref_glm)), tolerance = 1e-8)
      expect_equal(nobs(fit_glm), nobs(ref_glm))
      expect_true(fit_glm$converged)
      expect_true(fit_glm$score_max <= fit_glm$control$score_tol)
      expect_true(all(diff(fit_glm$deviance_history) <= 1e-10))
      close_glm(predict(fit_glm, glm_dat, type = "response"),
                predict(ref_glm, glm_dat, type = "response"), 1e-7, 1e-6)
      expect_equal(fit_glm$coeftable[, "p.value"],
                   2 * pnorm(-abs(coef(fit_glm) / sqrt(diag(vcov(fit_glm))))))
    }
  }
  # Promotion-like categorical/polynomial/triple interaction design.
  glm_case <- paste(link, "promotion design")
  fml_glm <- y ~ x + I(x^2) + z * sex * history + grp
  fit_glm <- dbglm(fml_glm, data = glm_dat, family = family_glm,
                  ssc = list(K.adj = FALSE))
  ref_glm <- glm(fml_glm, data = glm_dat, family = family_glm,
                 control = glm.control(epsilon = 1e-14, maxit = 100))
  expect_identical(names(coef(fit_glm)), names(coef(ref_glm)))
  close_glm(coef(fit_glm), coef(ref_glm), 1e-7, 1e-6)
  close_glm(vcov(fit_glm), vcov(ref_glm), 1e-8, 1e-5)
  close_glm(predict(fit_glm, glm_dat[1, ], type = "response"),
            predict(ref_glm, glm_dat[1, ], type = "response"), 1e-7, 1e-6)
  expect_true(all(is.finite(confint(fit_glm))))
  expect_equal(tidy(fit_glm, conf.int = TRUE)$estimate, unname(coef(fit_glm)))
  expect_equal(glance(fit_glm)$nobs, n_glm)
  expect_equal(unname(gof(fit_glm)["deviance"]), fit_glm$deviance)
  expect_equal(attr(logLik(fit_glm), "df"), length(coef(fit_glm)))
  expect_true(length(capture.output(print(fit_glm))) > 0)

  for (fml_glm in list(y ~ 1, y ~ 0 + x + z, y ~ x + z - 1)) {
    glm_case <- paste(link, paste(deparse(fml_glm), collapse = ""))
    fit_glm <- dbglm(fml_glm, data = glm_dat, family = family_glm,
                    ssc = list(K.adj = FALSE))
    ref_glm <- glm(fml_glm, data = glm_dat, family = family_glm,
                   control = glm.control(epsilon = 1e-14, maxit = 100))
    close_glm(coef(fit_glm), coef(ref_glm), 1e-7, 1e-6)
    close_glm(vcov(fit_glm), vcov(ref_glm), 1e-8, 1e-5)
    close_glm(predict(fit_glm, glm_dat[1:3, ]),
              predict(ref_glm, glm_dat[1:3, ]), 1e-7, 1e-6)
  }
}

# Caller-owned connection and all fit-owned registrations/tables are cleaned.
glm_con <- dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
dbWriteTable(glm_con, "glm_input", glm_dat)
before_glm <- dbListTables(glm_con)
fit_glm <- dbglm(family = stats::binomial("cloglog"), y ~ x + z, conn = glm_con, table = "glm_input", offset = "o")
expect_identical(dbListTables(glm_con), before_glm)
expect_true(dbIsValid(glm_con))
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = glm_con, data = glm_dat,
                   control = dbglm_control(glm_max_iter = 1)), "converge")
expect_identical(dbListTables(glm_con), before_glm)
expect_true(dbIsValid(glm_con))

# Exclusions are aligned across regressors, offset, weights and clustering.
glm_missing <- glm_dat
glm_missing$x[1] <- NA
glm_missing$o[2] <- Inf
glm_missing$a[3] <- 0
glm_missing$cluster[4] <- NA
fit_glm <- dbglm(family = stats::binomial("cloglog"), y ~ x + z, conn = glm_con, data = glm_missing,
                 weights = "a", offset = "o", vcov = ~cluster)
ref_glm <- reference_feglm(y ~ x + z, glm_missing[-(1:4), ], family = binomial("cloglog"),
                 weights = ~a, offset = ~o, vcov = ~cluster)
expect_equal(nobs(fit_glm), n_glm - 4)
close_glm(coef(fit_glm), coef(ref_glm), 1e-7, 1e-6)
close_glm(vcov(fit_glm), vcov(ref_glm), 1e-8, 1e-5)
expect_identical(dbListTables(glm_con), before_glm)
expect_error(predict(dbglm(family = stats::binomial("cloglog"), y ~ x | grp, data = glm_dat), glm_dat), "absorbed")
expect_error(dbglm(y ~ x, data = glm_dat, family = gaussian("log")), "gaussian")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, data = transform(glm_dat, y = 2)), "zero or one")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, data = transform(glm_dat, a = -1), weights = "a"), "weights")
expect_equal(coef(dbglm(family = stats::binomial("cloglog"), y ~ x, data = transform(glm_dat, a = 0.5), weights = "a")),
             coef(dbglm(family = stats::binomial("cloglog"), y ~ x, data = glm_dat)), tolerance = 1e-7)
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, data = glm_dat, control = dbglm_control(glm_tol = 0)), "positive")
expect_error(dbreg(y ~ x, data = glm_dat, family = binomial()), "dbglm")
expect_error(predict(fit_glm, transform(glm_dat, z = NULL)), "missing")
dbDisconnect(glm_con, shutdown = TRUE)
