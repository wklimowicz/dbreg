library(DBI)

# Counts and integer-valued sums can both arrive as BIGINT through DBI.
# Check fitted values against lm with the same data and both DuckDB modes.
(function() {
  set.seed(20261023)
  d = data.frame(x = rnorm(200), z = sample.int(7L, 200, TRUE))
  d[["y"]] = .5 + .7 * d[["x"]] - .3 * d[["z"]] + rnorm(200)
  d[["count"]] = as.integer(round(10 + 2 * d[["z"]] + rnorm(200)))
  d[["g"]] = rep(seq_len(10), 20)
  d[["binary"]] = rbinom(200, 1, binomial("cloglog")[["linkinv"]](
    -.5 + .4 * d[["x"]] + .3 * sin(d[["g"]])))
  d[["instrument"]] = d[["x"]] + .3 * rnorm(200)
  bigint_modes = "numeric"
  if (requireNamespace("bit64", quietly = TRUE)) bigint_modes = c(bigint_modes, "integer64")
  for (bigint in bigint_modes) {
    conn = dbConnect(duckdb::duckdb(bigint = bigint))
    tryCatch({
      dbWriteTable(conn, "input", d)
      for (fml in list(y ~ x + z, count ~ z)) {
        ref = lm(fml, d)
        for (kind in c("iid", "hc1")) {
          got = dbreg(fml, conn = conn, table = "input", strategy = "moments", vcov = kind)
          expect_equal(coef(got), coef(ref), tolerance = 1e-8)
          if (kind == "iid") expect_equal(as.numeric(vcov(got)), as.numeric(vcov(ref)), tolerance = 1e-8)
          else {
            X = model.matrix(ref)
            bread = solve(crossprod(X))
            scores = X * residuals(ref)
            robust = bread %*% crossprod(scores) %*% bread * nrow(d) / df.residual(ref)
            expect_equal(as.numeric(vcov(got)), as.numeric(robust), tolerance = 1e-8)
          }
        }
      }
      # FE ranks arrive as BIGINT as well. Convert each key vector before
      # unlisting: unlist() by itself strips the integer64 representation.
      ref_glm = fixest::feglm(binary ~ x | g + z, d,
        family = binomial("cloglog"), fixef.rm = "none", glm.tol = 1e-12)
      for (kind in c("iid", "hc1")) {
        got_glm = dbglm(binary ~ x | g + z, conn = conn, table = "input",
          family = binomial("cloglog"), vcov = kind)
        ref_vcov = vcov(ref_glm, vcov = if (kind == "hc1") "hetero" else "iid")
        expect_equal(coef(got_glm), coef(ref_glm), tolerance = 1e-6)
        expect_equal(as.numeric(vcov(got_glm)), as.numeric(ref_vcov), tolerance = 1e-6)
      }
      for (kind in c("iid", "hc1")) {
        got_iv = dbivreg(y ~ z | x ~ instrument, conn = conn, table = "input",
          strategy = "moments", vcov = kind)
        ref_iv = fixest::feols(y ~ z | x ~ instrument, d,
          vcov = if (kind == "hc1") "hetero" else "iid", fixef.rm = "none")
        ref_names = sub("^fit_", "", names(coef(ref_iv)))
        ref_order = match(names(coef(got_iv)), ref_names)
        expect_equal(unname(coef(got_iv)), unname(coef(ref_iv)[ref_order]), tolerance = 1e-8)
        expect_equal(as.numeric(vcov(got_iv)), as.numeric(vcov(ref_iv)[ref_order, ref_order]), tolerance = 1e-8)
      }
    }, finally = dbDisconnect(conn, shutdown = TRUE))
  }
})()
