# Per-fit projection controls override options without mutating them. Both
# public estimators use the same projected frame and validate before DB setup.
local({
  saved = options(dbreg.ap_max_iter = 1L,dbreg.ap_tol = 1e-10,
    dbreg.ap_method = getOption("dbreg.ap_method"))
  on.exit(options(saved),add = TRUE)
  set.seed(723)
  d = data.frame(x = rnorm(180),z = rnorm(180),g = sample(1:8,180,TRUE),h = sample(1:6,180,TRUE),w = runif(180,.4,2))
  d[["endog"]] = d[["z"]]+rnorm(180)
  d[["y"]] = .5*d[["x"]]+.8*d[["endog"]]+sin(d[["g"]])+d[["h"]]+rnorm(180)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn,shutdown = TRUE),add = TRUE)
  DBI::dbWriteTable(conn,"control_fixture",d)
  for (method in c("join","window")) {
    options(dbreg.ap_method = method)
    for (iv in c(FALSE,TRUE)) {
      fml = if (iv)y~x|g+h|endog~z else y~x|g+h
      estimator = if (iv)dbivreg else dbreg
      fit = estimator(fml,conn = conn,table = "control_fixture",weights = "w",strategy = "demean",fe_max_iter = 100L,fe_tol = 1e-9)
      ref = fixest::feols(fml,d,weights = ~w,fixef.rm = "none")
      b = coef(ref);names(b) = sub("^fit_","",names(b))
      expect_equal(as.numeric(coef(fit)),as.numeric(b[names(coef(fit))]),tolerance = 1e-8)
      expect_identical(getOption("dbreg.ap_max_iter"),1L)
      expect_identical(getOption("dbreg.ap_tol"),1e-10)
      expect_true(fit[["ap_iterations"]]<=100L)
      expect_identical(fit[["ap_tolerance"]],1e-9)
      expect_identical(fit[["ap_max_iter"]],100L)
      expect_true(fit[["ap_max_relative"]]<=1e-9)
      expect_error(estimator(fml,conn = conn,table = "control_fixture",weights = "w",strategy = "demean",fe_max_iter = 1L),"within 1 iterations")
      expect_identical(DBI::dbListTables(conn),"control_fixture")
      for (bad in list(0,-1,NA_real_,Inf,"small",c(1,2),NULL)) {
        expect_error(estimator(fml,conn = conn,table = "control_fixture",fe_tol = bad),"fe_tol")
        expect_error(estimator(fml,conn = conn,table = "control_fixture",fe_max_iter = bad),"fe_max_iter")
      }
      expect_error(estimator(fml,conn = conn,table = "control_fixture",fe_max_iter = 1.5),"fe_max_iter")
      expect_error(estimator(fml,conn = conn,table = "control_fixture",fe_max_iter = 2^32),"fe_max_iter")
      expect_true(DBI::dbIsValid(conn))
      expect_identical(DBI::dbListTables(conn),"control_fixture")
    }
  }
})
