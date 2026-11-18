# A fixed effect absorbs a factor's intercepts, not its slopes: x:f | f keeps every level.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(9); n = 1500
  d = data.frame(id = sample(60, n, TRUE), edu = sample(c("hs", "ba", "ma"), n, TRUE),
    reg = sample(c("n", "s", "e"), n, TRUE), hours = sample(20:60, n, TRUE))
  d$wage = 40 + c(hs = 6, ba = 9, ma = 11)[d$edu] * d$hours + rnorm(60, 0, 15)[d$id] + rnorm(n, 0, 10)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  for (f in list(wage ~ hours:edu | id + edu, wage ~ hours + hours:edu | edu)) {
    ref = fixest::feols(f, d, ssc = fixest::ssc(K.fixef = "full"), fixef.tol = 1e-11, notes = FALSE)
    for (s in c("auto", "demean", "compress")) {
      fit = dbreg(f, conn = conn, table = "d", strategy = s)
      k = names(coef(ref))
      expect_equal(sort(names(coef(fit))), sort(k), info = paste(deparse(f), s))
      expect_equal(unname(coef(fit)[k]), unname(coef(ref)), tolerance = 1e-8, info = paste(deparse(f), s))
      expect_equal(as.numeric(sqrt(diag(vcov(fit)))[k]), as.numeric(fixest::se(ref)),
        tolerance = 1e-8, info = paste(deparse(f), s))
    }
  }
  # Factor-by-factor cells beside a factor FE span the same space as fixest's.
  f = wage ~ hours + edu:reg | id + reg
  ref = fixest::feols(f, d, ssc = fixest::ssc(K.fixef = "full"), fixef.tol = 1e-11, notes = FALSE)
  for (s in c("auto", "demean", "compress")) {
    fit = dbreg(f, conn = conn, table = "d", strategy = s)
    expect_equal(length(coef(fit)), length(coef(ref)), info = s)
    expect_equal(unname(coef(fit)["hours"]), unname(coef(ref)["hours"]), tolerance = 1e-8, info = s)
    expect_equal(gof(fit)[["r2"]], unname(fixest::r2(ref, "r2")), tolerance = 1e-8, info = s)
  }
})
