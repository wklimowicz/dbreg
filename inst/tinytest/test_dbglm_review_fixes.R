library(tinytest)

# 1. Structural-only shifted predictors must not be dropped or diverge
d1 = data.frame(x = 1e8 + seq_len(20))
d1$z = d1$x - 1e8
d1$y = 2 * d1$z + cos(d1$z)
for (strategy in c("moments", "compress")) {
  a = dbglm(y ~ x, data = d1, strategy = strategy, retain = "inference")
  b = dbglm(y ~ z, data = d1, strategy = strategy, retain = "inference")
  expect_true(a$converged)
  expect_equal(unname(coef(a)[["x"]]), unname(coef(b)[["z"]]), tolerance = 1e-8)
  expect_equal(a$deviance, b$deviance, tolerance = 1e-8)
  expect_equal(unname(vcov(a)[2L, 2L]), unname(vcov(b)[2L, 2L]), tolerance = 1e-10)
  expect_equal(unname(vcov(a, vcov = "hc1")[2L, 2L]),
    unname(vcov(b, vcov = "hc1")[2L, 2L]), tolerance = 1e-10)
  close(a)
  close(b)
}

# 2. Saturated probabilities agreeing with outcome are not separation
d2 = data.frame(
  x = c(rep(-1, 100), rep(1, 100), 40),
  y = c(rep(c(0, 0, 0, 1), 25), rep(c(1, 1, 1, 0), 25), 1)
)
g2 = glm(y ~ x, data = d2, family = binomial())
f2 = dbglm(y ~ x, data = d2, family = binomial())
expect_true(f2$finite_mle)
expect_true(f2$converged)
expect_false(anyNA(vcov(f2)))
expect_equal(unname(coef(f2)), unname(coef(g2)), tolerance = 1e-6)

# Large known offset agreeing with outcome must still initialize and converge
d2_off = data.frame(
  x = c(rep(-1, 50), rep(1, 50), 0),
  off = c(rep(0, 100), 40),
  y = c(rep(c(0, 0, 0, 1), 12), 0, 1, rep(c(1, 1, 1, 0), 12), 1, 0, 1)
)
g2_off = glm(y ~ x, data = d2_off, family = binomial(), offset = off)
f2_off = dbglm(y ~ x, data = d2_off, family = binomial(), offset = "off")
expect_true(f2_off$finite_mle)
expect_true(f2_off$converged)
expect_equal(unname(coef(f2_off)), unname(coef(g2_off)), tolerance = 1e-6)

# 3. Fitting-time and retained clustered covariance use identical K.fixef policy
set.seed(12)
d3 = data.frame(
  id = 1:200, x = rnorm(200),
  g = rep(1:5, each = 40), h = rep(1:10, each = 20)
)
d3$y = rpois(200, exp(0.2 + 0.3 * d3$x + d3$g / 10))
f3 = dbglm(y ~ x | g + h, data = d3, family = poisson(),
  vcov = ~g, ssc = "nested", retain = "inference")
v3_fit = vcov(f3)
v3_ret = vcov(f3, vcov = ~g)
expect_equal(attr(v3_fit, "df.K"), attr(v3_ret, "df.K"))
expect_equal(as.numeric(v3_fit), as.numeric(v3_ret), tolerance = 1e-12)
close(f3)

# 4. Large Gaussian response location converges cleanly
set.seed(801)
d4 = data.frame(x = rnorm(200), g = rep(1:10, each = 20))
d4$y = 3 + 2 * d4$x + rnorm(200)
f4_0 = dbglm(y ~ x, data = d4)
d4$y = d4$y + 1e8
f4_8 = dbglm(y ~ x, data = d4, retain = "all")
f4_fe = dbglm(y ~ x | g, data = d4, retain = "all")
expect_true(f4_8$converged)
expect_true(f4_fe$converged)
expect_equal(unname(coef(f4_8)[["x"]]), unname(coef(f4_0)[["x"]]), tolerance = 1e-8)
expect_equal(fitted(f4_fe), predict(f4_fe, d4), tolerance = 1e-6)
close(f4_8)
close(f4_fe)

# 5. Poisson deviance avoids cancellation for large y near mu
k = 1e9
d5 = data.frame(
  x = rep(c(0, 1), each = 10),
  y = c(rep(k, 5), rep(k + 1, 5), rep(k + 2, 5), rep(k + 3, 5))
)
for (strategy in c("moments", "compress")) {
  f5 = dbglm(y ~ x, data = d5, family = poisson(), strategy = strategy)
  expect_true(f5$converged)
  expect_true(f5$deviance > 0)
  expect_equal(as.numeric(f5$deviance), 5e-9, tolerance = 1e-7)
}

# 6. Gaussian logLik(type = "weighted") includes dispersion df (+1) and matches BIC(glm)
d6 = data.frame(x = rnorm(200))
d6$y = 3 + 2 * d6$x + rnorm(200)
g6 = glm(y ~ x, data = d6)
f6 = dbglm(y ~ x, data = d6)
ll_w = logLik(f6, type = "weighted")
expect_equal(attr(ll_w, "df"), attr(logLik(g6), "df"))
expect_equal(BIC(ll_w), BIC(g6), tolerance = 1e-10)

# 7. Complete separation must be rejected even at relaxed score_tol, and
# separation = "prediction" must return boundary fits across binomial links
d_sep = data.frame(x = rep(c(-1, 1), each = 20), y = rep(c(0, 1), each = 20))
for (lnk in c("logit", "probit", "cloglog")) {
  expect_error(
    dbglm(y ~ x, data = d_sep, family = binomial(link = lnk),
      control = dbglm_control(score_tol = 1e-7)),
    "probability boundary"
  )
  pred_fit = suppressWarnings(dbglm(y ~ x, data = d_sep, family = binomial(link = lnk),
    separation = "prediction"))
  expect_false(pred_fit$finite_mle)
  expect_false(pred_fit$converged)
  expect_true(all(is.na(vcov(pred_fit))))
}

# 8. Gaussian warm starts in centered response coordinates
set.seed(4)
d_warm = data.frame(id = 1:120, x = rnorm(120), g = rep(1:6, each = 20))
d_warm$y = 1e8 + 2 + d_warm$x + d_warm$g / 10 + rnorm(120)
d_warm$eta = fitted(lm(y ~ x + factor(g), d_warm))
d_warm$w = 1
warm_cold = dbglm(y ~ x | g, data = d_warm)
warm_direct = dbglm(y ~ x | g, data = d_warm, etastart = "eta")
expect_true(warm_direct$converged)
expect_equal(unname(coef(warm_direct)), unname(coef(warm_cold)), tolerance = 1e-8)

# 9. Large known Gaussian offsets and partition-of-unity categorical models (y ~ 0 + g)
set.seed(801)
d_off = data.frame(x = rnorm(200), g = factor(rep(1:10, each = 20)))
d_off$o = 1e8
d_off$y = d_off$o + 3 + 2 * d_off$x + rnorm(200)
d_off$adjusted = d_off$y - d_off$o
f_off = dbglm(y ~ x, data = d_off, offset = "o")
f_adj = dbglm(adjusted ~ x, data = d_off)
expect_true(f_off$converged)
expect_equal(unname(coef(f_off)), unname(coef(f_adj)), tolerance = 1e-8)

f_off_fe = dbglm(y ~ x | g, data = d_off, offset = "o")
f_adj_fe = dbglm(adjusted ~ x | g, data = d_off)
expect_true(f_off_fe$converged)
expect_equal(unname(coef(f_off_fe)), unname(coef(f_adj_fe)), tolerance = 1e-8)

d_part = data.frame(g = factor(rep(1:4, each = 25)))
d_part$y = 1e8 + as.numeric(d_part$g) + rnorm(100, sd = 0.1)
f_part = dbglm(y ~ 0 + g, data = d_part)
expect_true(f_part$converged)
expect_equal(unname(coef(f_part)), unname(coef(lm(y ~ 0 + g, data = d_part))), tolerance = 1e-8)

# 10. R and keyed SQL prediction preserve stable fitting coordinates on shifted predictors
d_pred = data.frame(id = 1:60, x = 1e15 + seq_len(60))
d_pred$z = d_pred$x - 1e15
d_pred$y = 2 * d_pred$z + cos(d_pred$z)
conn_pred = DBI::dbConnect(duckdb::duckdb())
f_pred = dbglm(y ~ x, data = d_pred, conn = conn_pred, retain = "inference")
saved_pred = fitted(f_pred)
r_pred = predict(f_pred, d_pred)
sql_pred = predict(f_pred, conn = conn_pred, newdata = d_pred, row_id = "id")
expect_equal(r_pred, saved_pred, tolerance = 1e-12)
expect_equal(as.data.frame(sql_pred)$prediction, saved_pred, tolerance = 1e-12)
close(sql_pred)
close(f_pred)
DBI::dbDisconnect(conn_pred, shutdown = TRUE)

# 11. Retained Gaussian robust inference and residuals preserve precision at large y
set.seed(5)
d_rob = data.frame(x = rnorm(200), g = rep(1:10, each = 20))
d_rob$y = 1e12 + 0.002 * d_rob$x + rnorm(200, sd = 0.0005)
f_rob = dbglm(y ~ x | g, data = d_rob, vcov = "hc1", retain = "inference")
expect_equal(as.numeric(vcov(f_rob)), as.numeric(vcov(f_rob, vcov = "hc1")), tolerance = 1e-12)
expect_equal(as.numeric(residuals(f_rob, type = "response")),
  unname(residuals(lm(I(y - 1e12) ~ x + factor(g), data = d_rob))), tolerance = 1e-10)
close(f_rob)

# 12. Poisson saturated log-likelihood remains accurate at large counts (k = 1e15)
k15 = 1e15
d_p15 = data.frame(
  x = rep(c(0, 1), each = 10),
  y = c(rep(k15, 5), rep(k15 + 1, 5), rep(k15 + 2, 5), rep(k15 + 3, 5))
)
f_p15 = dbglm(y ~ x, data = d_p15, family = poisson(), control = dbglm_control(score_tol = 1e-6))
exact_ll15 = sum(dpois(d_p15$y, exp(predict(f_p15, d_p15)), log = TRUE))
expect_equal(as.numeric(logLik(f_p15)), exact_ll15, tolerance = 1e-10)

# 13. Quasi-separation and relaxed score_tol separation vs finite extreme-leverage saturation
d_quasi = data.frame(
  x1 = c(0, 0, 0, 0, 1, 1, 1, 1),
  x2 = c(0, 1, 0, 1, 0, 1, 0, 1),
  y  = c(0, 0, 1, 1, 1, 1, 1, 1)
)
for (st in c("moments", "compress")) {
  for (lnk in c("logit", "probit", "cloglog")) {
    expect_error(
      dbglm(y ~ x1 + x2, data = d_quasi, family = binomial(link = lnk), strategy = st),
      "probability boundary"
    )
    pf_q = suppressWarnings(dbglm(y ~ x1 + x2, data = d_quasi,
      family = binomial(link = lnk), strategy = st, separation = "prediction"))
    expect_false(pf_q$finite_mle)
    expect_false(pf_q$converged)
  }
}
d_small = data.frame(x = c(-1, -1, 1, 1), y = c(0, 0, 1, 1))
for (lnk in c("logit", "probit", "cloglog")) {
  expect_error(
    dbglm(y ~ x, data = d_sep, family = binomial(link = lnk), control = dbglm_control(score_tol = 1e-6)),
    "probability boundary"
  )
  expect_error(
    dbglm(y ~ x, data = d_small, family = binomial(link = lnk), control = dbglm_control(score_tol = 1e-6)),
    "probability boundary"
  )
}
d_fin_sat = data.frame(x = c(-1e9, -1, 1, 1e9), y = c(0, 1, 0, 1))
f_fin_sat = dbglm(y ~ x, data = d_fin_sat, family = binomial(link = "cloglog"))
expect_true(f_fin_sat$converged)
expect_true(f_fin_sat$finite_mle)

# 14. Partition-of-unity centering with reordered terms, 1-level factors, shifted x, and null_model
set.seed(801)
d_p_extra = data.frame(x = rnorm(200), g = factor(rep(1:10, each = 20)))
d_p_extra$y_big = 1e12 + as.numeric(d_p_extra$g) + 2 * d_p_extra$x + rnorm(200, sd = 0.001)
f_xg = dbglm(y_big ~ 0 + x + g, data = d_p_extra)
expect_true(f_xg$converged)
expect_equal(unname(coef(f_xg)[["x"]]),
  unname(coef(lm(I(y_big - 1e12) ~ 0 + x + g, data = d_p_extra))[["x"]]), tolerance = 1e-10)

d_1lvl = data.frame(g = factor(rep("a", 50)), y = 1e12 + rnorm(50, sd = 0.001))
f_1lvl = dbglm(y ~ 0 + g, data = d_1lvl)
expect_true(f_1lvl$converged)
expect_equal(unname(coef(f_1lvl)[1L]), mean(d_1lvl$y), tolerance = 1e-4)

d_p_extra$y_n = 1e12 + as.numeric(d_p_extra$g) + rnorm(200, sd = 0.001)
f_0g_null = dbglm(y_n ~ 0 + g, data = d_p_extra, null_model = TRUE)
expect_equal(f_0g_null$null_diagnostics$weighted$logLik,
  as.numeric(logLik(lm(y_n ~ 0, data = d_p_extra))), tolerance = 1e-6)

d_p_extra$x_big = 1e15 + seq_len(200)
d_p_extra$y_xb = 2 * (d_p_extra$x_big - 1e15) + as.numeric(d_p_extra$g) + rnorm(200, sd = 0.01)
f_0g_xb = dbglm(y_xb ~ 0 + g + x_big, data = d_p_extra)
expect_true(f_0g_xb$converged)
expect_equal(unname(coef(f_0g_xb)[["x_big"]]),
  unname(utils::tail(coef(lm(y_xb ~ 0 + g + I(x_big - 1e15), data = d_p_extra)), 1L)),
  tolerance = 1e-10)

# 15. Absorbed FE quasi-separation raises numerical probability boundary error
d_fe_q = data.frame(
  g = factor(rep(1:2, each = 20)),
  x = c(rep(0, 10), rep(1, 10), rep(0, 10), rep(1, 10)),
  y = c(rep(c(0, 1), 5), rep(1, 10), rep(c(0, 1), 5), rep(1, 10))
)
for (lnk in c("logit", "probit", "cloglog")) {
  expect_error(
    dbglm(y ~ x | g, data = d_fe_q, family = binomial(link = lnk)),
    "probability boundary"
  )
}
