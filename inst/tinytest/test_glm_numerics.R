library(dbreg)
library(DBI)

con_glm <- dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
eta_glm <- c(-1000, -100, -36, -31, -30, -18, -3, 0, 3, 4, 30, 31, 100, 700, 1000)
for (link_glm in c("cloglog", "logit")) {
  kernel_glm <- dbreg:::dbglm_family_sql(link_glm)
  sql_glm <- paste0("WITH values_glm(eta) AS (VALUES ",
    paste0("(", eta_glm, ")", collapse = ","), ") SELECT ",
    kernel_glm$mu, " AS mu, ", kernel_glm$derivative, " AS d FROM values_glm")
  got_glm <- dbGetQuery(con_glm, sql_glm)
  fam_glm <- binomial(link_glm)
  expect_true(all(is.finite(as.matrix(got_glm))))
  expect_true(all(abs(got_glm$mu - fam_glm$linkinv(eta_glm)) < 5e-16))
  expect_true(all(abs(got_glm$d - fam_glm$mu.eta(eta_glm)) < 5e-16))
}

set.seed(111)
bad_start_dat <- data.frame(x = rnorm(180), z = rnorm(180))
bad_start_dat$y <- rbinom(180, 1, plogis(-0.4 + 0.5 * bad_start_dat$x - 0.3 * bad_start_dat$z))
for (link_glm in c("cloglog", "logit")) {
  fit_good <- dbglm(y ~ x + z, conn = con_glm, data = bad_start_dat,
                    family = binomial(link_glm))
  fit_bad <- dbglm(y ~ x + z, conn = con_glm, data = bad_start_dat,
                   family = binomial(link_glm), start = c(2, -1, 1))
  expect_true(fit_bad$step_halvings > 0)
  expect_true(all(diff(fit_bad$deviance_history) <= 1e-10))
  expect_equal(coef(fit_bad), coef(fit_good), tolerance = 1e-7)
  expect_equal(vcov(fit_bad), vcov(fit_good), tolerance = 1e-7)
}

# Lock the basis, preserve original coefficient order, and solve in scaled units.
alias_glm <- transform(bad_start_dat, duplicate = x, large = 1e6 * z)
fit_alias <- dbglm(y ~ x + duplicate + large, conn = con_glm, data = alias_glm,
                   family = binomial("logit"), ssc = list(K.adj = FALSE))
ref_alias <- glm(y ~ x + large, alias_glm, family = binomial("logit"),
                 control = glm.control(epsilon = 1e-14))
expect_identical(names(coef(fit_alias)), c("(Intercept)", "x", "large"))
expect_identical(names(coef(fit_alias, collin = TRUE)), c("(Intercept)", "x", "duplicate", "large"))
expect_true(is.na(coef(fit_alias, collin = TRUE)["duplicate"]))
expect_equal(coef(fit_alias), coef(ref_alias), tolerance = 1e-7)
# Completed dbreg covariance carries inference policy metadata beyond glm's V.
expect_identical(dim(vcov(fit_alias)), dim(vcov(ref_alias)))
expect_identical(dimnames(vcov(fit_alias)), dimnames(vcov(ref_alias)))
expect_equal(as.numeric(vcov(fit_alias)), as.numeric(vcov(ref_alias)), tolerance = 1e-7)
expect_equal(predict(fit_alias, alias_glm), unname(predict(ref_alias, alias_glm)), tolerance = 1e-7)

ref_dat <- transform(bad_start_dat, code = rep(1:3, length.out = 180))
fit_ref <- dbglm(y ~ x + code, conn = con_glm, data = ref_dat, factor_refs = c(code = 2),
                 family = binomial("logit"), ssc = list(K.adj = FALSE))
glm_ref <- glm(y ~ x + code, transform(ref_dat, code = relevel(factor(code), "2")),
                family = binomial("logit"), control = glm.control(epsilon = 1e-14))
expect_identical(names(coef(fit_ref)), names(coef(glm_ref)))
expect_equal(coef(fit_ref), coef(glm_ref), tolerance = 1e-7)
expect_equal(predict(fit_ref, ref_dat[1, ], type = "response"),
              unname(predict(glm_ref, transform(ref_dat[1, ],
                code = factor(code, levels = glm_ref$xlevels$code)), type = "response")), tolerance = 1e-7)
expect_error(predict(fit_ref, transform(ref_dat[1, ], code = 99)), "unseen")
missing_new <- ref_dat[1:3, ]
missing_new$x[2] <- NA
expect_true(is.na(predict(fit_ref, missing_new)[2]))

modified_family <- binomial("logit")
modified_family$linkinv <- function(eta) rep(0.5, length(eta))
expect_error(dbglm(y ~ x, conn = con_glm, data = ref_dat, family = modified_family), "Modified")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = con_glm, data = transform(ref_dat, y = 0)), "finite")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = con_glm, data = ref_dat, start = c(NA, 0)), "finite")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = con_glm, data = ref_dat, vcov = ~x + z), "single-variable")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = con_glm, data = transform(ref_dat, c = 1), vcov = ~c), "two clusters")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = con_glm, data = ref_dat, ssc = list(K.exact = TRUE)), "Unsupported")
expect_error(dbglm(family = stats::binomial("cloglog"), y ~ x, conn = con_glm, data = ref_dat,
                   control = dbglm_control(glm_max_iter = 0.5)), "integer")
expect_error(confint(fit_ref, level = 1), "between")
expect_error(predict(fit_ref), "newdata")

separated_glm <- data.frame(x = rep(c(-1, 1), each = 10), y = rep(c(0, 1), each = 10))
expect_error(dbglm(y ~ x, conn = con_glm, data = separated_glm, family = binomial("logit")),
              "boundary|converge|singular")
expect_equal(length(dbListTables(con_glm)), 0L)
dbDisconnect(con_glm, shutdown = TRUE)
