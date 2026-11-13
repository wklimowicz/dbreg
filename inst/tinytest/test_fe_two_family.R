local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  saved = options(dbreg.ap_solver = "auto", dbreg.ap_method = NULL, dbreg.duckdb_ap_method = "cg")
  on.exit(options(saved), add = TRUE)
  set.seed(11125)
  n = 360L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = sample(1:12,n,TRUE),
    h = sample(1:9,n,TRUE), w = runif(n,.2,2), o = rnorm(n))
  d[["endog"]] = .8*d[["z"]]+rnorm(n)
  d[["y"]] = d[["o"]]+.7*d[["x"]]+.4*d[["endog"]]+sin(d[["g"]])+cos(d[["h"]])+rnorm(n)
  DBI::dbWriteTable(conn, "two_family", d)
  initial = sort(DBI::dbListTables(conn))
  adjusted = d; adjusted[["y"]] = d[["y"]]-d[["o"]]
  for (iv in c(FALSE,TRUE)) for (kind in list("iid","hc1",~g)) {
    fml = if (iv) y~x|g+h|endog~z else y~x+z|g+h
    fun = if (iv) dbivreg else dbreg
    fit = fun(fml, conn=conn, table="two_family", weights=~w, offset=~o,
      vcov=kind, strategy="demean", fe_max_iter=100L, fe_tol=1e-9)
    ref = fixest::feols(fml, adjusted, weights=~w, vcov=if (identical(kind,"hc1")) "hetero" else kind,
      ssc=fixest::ssc(K.fixef="full"), fixef.rm="none", fixef.tol=1e-11, notes=FALSE)
    ix = match(names(coef(fit)),sub("^fit_","",names(coef(ref))))
    expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)[ix]),tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)[ix,ix,drop=FALSE]),tolerance=1e-8)
    expect_identical(fit[["ap_method"]],"bounded_pcg")
    expect_identical(fit[["ap_solver"]],"schur_pcg")
    expect_true(fit[["ap_solver_iterations"]]>0L)
    expect_true(fit[["ap_max_abs"]]<=1e-9 && fit[["ap_max_relative"]]<=1e-9)
    expect_identical(fit[["ap_observation_tables"]],1L)
    expect_identical(sort(DBI::dbListTables(conn)),initial)
  }
  inputs = list2env(list(conn=conn, from_statement="FROM two_family", fe=c("g","h"), yvar="y",
    weights="w", cluster_var="id", verbose=FALSE, projection_limits=dbreg:::projection_limits()),parent=emptyenv())
  fit = dbreg:::dbreg_moment_projections(inputs,c("x","z"),c("x","z"),100L,1e-9,two_family_pcg=TRUE)
  actual = DBI::dbGetQuery(conn,paste("SELECT y_tilde,x_tilde,z_tilde FROM",fit[["table"]],"ORDER BY id"))
  D = model.matrix(~factor(g)+factor(h),d); sw=sqrt(d[["w"]])
  ref_rows = qr.resid(qr(D*sw),as.matrix(d[c("y","x","z")])*sw)/sw
  expect_equal(unname(as.matrix(actual)),unname(ref_rows),tolerance=1e-8)
  dbreg:::release_fit_resources(fit[["resources"]])
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  attempt = dbreg:::projection_attempt("pcg",dbreg:::dbreg_moment_projections(inputs,c("x","z"),c("x","z"),1L,1e-9,two_family_pcg=TRUE))
  expect_identical(attempt[["reason"]],"nonconverged")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  expect_error(dbreg(y~x+z|g+h,conn=conn,table="two_family",weights=~w,strategy="demean",fe_max_iter=1L),"within 1 iterations")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  # Explicit moments takes precedence; an explicit join bypasses bounded paths.
  prior = options(dbreg.ap_solver="moments")
  fit = dbreg(y~x+z|g+h,conn=conn,table="two_family",weights=~w,strategy="demean")
  expect_identical(fit[["ap_method"]],"bounded_moments")
  options(prior)
  prior = options(dbreg.duckdb_ap_method="join")
  fit = dbreg(y~x+z|g+h,conn=conn,table="two_family",weights=~w,strategy="demean")
  expect_identical(fit[["ap_method"]],"row_tables")
  options(prior)
  # All canonical user budgets survive the per-call sparse representation choice.
  for (name in c("summary_rows","transfer_cells","design_entries","memory_bytes","upload_bytes","cross_entries")) {
    limited = list2env(as.list(inputs),parent=emptyenv())
    limited[["projection_limits"]][[name]] = if (name=="cross_entries") 0 else 1
    attempt = dbreg:::projection_attempt("pcg",dbreg:::dbreg_moment_projections(limited,c("x","z"),c("x","z"),100L,1e-9,two_family_pcg=TRUE))
    expect_identical(attempt[["status"]],"declined")
    expect_identical(sort(DBI::dbListTables(conn)),initial)
  }
  # Disconnected components and native keys do not invent identification axes.
  d[["g"]] = d[["g"]]+12*(d[["id"]]%%2L)
  d[["h"]] = d[["h"]]+9*(d[["id"]]%%2L)
  fit = dbreg(y~x+z|g+h,data=d,weights=~w,offset=~o,strategy="demean",fe_tol=1e-9)
  d[["y"]] = d[["y"]]-d[["o"]]
  ref = fixest::feols(y~x+z|g+h,d,weights=~w,fixef.rm="none",ssc=fixest::ssc(K.fixef="full"),notes=FALSE)
  expect_identical(fit[["ap_method"]],"bounded_pcg")
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-8)
  expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)),tolerance=1e-8)
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  DBI::dbExecute(conn,"CREATE TEMP TABLE typed_two AS SELECT * EXCLUDE(g,h),CAST('9007199254740992' AS BIGINT)+CAST(g AS BIGINT) AS g,'a|/'||CAST(h AS VARCHAR) AS h FROM two_family")
  fit = dbreg(y~x+z|g+h,conn=conn,table="typed_two",weights=~w,offset=~o,strategy="demean",
    row_id=c("id","g","h"),sample_table="two_keys",fe_tol=1e-9)
  keys = DBI::dbGetQuery(conn,"SELECT id,CAST(g AS VARCHAR) AS g,h FROM two_keys ORDER BY id")
  expect_identical(fit[["ap_method"]],"bounded_pcg")
  expect_identical(keys[["g"]],paste0("900719925474",sprintf("%04d",992+(d[["g"]]-12*(d[["id"]]%%2L)))))
  DBI::dbRemoveTable(conn,"two_keys");DBI::dbRemoveTable(conn,"typed_two")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  # The final independent row certificate cannot be replaced by PCG residuals.
  project = dbreg:::dbreg_moment_projections
  environment(project) = new.env(parent=environment(project))
  environment(project)[["fe_projection_certificate"]] = function(...) list(max_abs=Inf,max_relative=Inf)
  attempt = dbreg:::projection_attempt("pcg",project(inputs,c("x","z"),c("x","z"),100L,1e-9,two_family_pcg=TRUE))
  expect_identical(attempt[["reason"]],"certificate_failed")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  environment(project)[["fe_projection_certificate"]] = function(...) stop("unexpected two-family driver error")
  expect_error(dbreg:::projection_attempt("pcg",project(inputs,c("x","z"),c("x","z"),100L,1e-9,two_family_pcg=TRUE)),"unexpected two-family driver error")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  for (bad in list("unknown",NA_character_,c("join","cg"),1L)) {
    prior=options(dbreg.duckdb_ap_method=bad)
    expect_error(dbreg(y~x|g+h,conn=conn,table="two_family",strategy="demean"),"duckdb_ap_method")
    options(prior)
  }
})
