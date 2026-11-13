# Explicitly centered row designs provide an independent oracle when raw
# covariance algorithms cannot represent their moving means accurately.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(76592)
  n = 240L
  d = data.frame(x = 1e14 + runif(n, -1, 1), z = rnorm(n),
    instrument = rnorm(n), extra = rnorm(n), cluster = rep(1:20, each = 12))
  d[["dx"]] = d[["x"]] - 1e14
  d[["endog"]] = .7*d[["instrument"]] + .3*d[["extra"]] + rnorm(n)
  d[["y"]] = 2*d[["dx"]] + .4*d[["z"]] + .8*d[["endog"]] + rnorm(n, sd = .05)
  for (weighted in c(FALSE, TRUE)) {
    d[["w"]] = if (weighted) seq(.3, 1.7, length.out = n) else 1
    DBI::dbWriteTable(conn, "review4_moments", d, overwrite = TRUE)
    centered = d
    centered[["x"]] = d[["dx"]]
    for (kind in list("iid", "hc1", ~cluster)) {
      for (iv in c(FALSE, TRUE)) {
        fml = if (iv) y ~ x + z | endog ~ instrument + extra else y ~ x + z + endog
        fun = if (iv) dbivreg else dbreg
        fit = fun(fml, conn = conn, table = "review4_moments",
          weights = if (weighted) "w" else NULL, strategy = "moments", vcov = kind)
        ref = fixest::feols(fml, centered, weights = ~w,
          vcov = if (identical(kind, "hc1")) "hetero" else kind, notes = FALSE)
        positions = match(names(coef(fit))[-1L], sub("^fit_", "", names(coef(ref))))
        expect_equal(unname(coef(fit)[-1L]), unname(coef(ref)[positions]), tolerance = 1e-9)
        expect_equal(as.numeric(vcov(fit)[-1L, -1L]),
          as.numeric(vcov(ref)[positions, positions]), tolerance = 1e-9)
        if (iv) expect_equal(fit[["diagnostics"]][["first_stage_f"]][["endog"]][["stat"]],
          fixest::fitstat(ref, "ivf1", simplify = TRUE)[["stat"]], tolerance = 1e-9)
      }
    }
  }
  exact = data.frame(x = 1e14 + c(0, 1, 1, 0, 2, 2)/64)
  exact[["y"]] = 2*(exact[["x"]]-1e14)
  expect_equal(unname(coef(dbreg(y ~ x, data = exact, strategy = "moments"))[["x"]]), 2,
    tolerance = 1e-12)
})
