local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11169)
  core = expand.grid(g = 1:3, h = 1:3, x = -1:1, rep = 1:3)
  tail = data.frame(g = c(4,4,5,5,6), h = c(3,4,4,5,5), x = c(1,-1,0,1,-1), rep = 1)
  d = rbind(core, tail)
  n = nrow(d)
  d$w = runif(n, .3, 2)
  d$o = sin(d$rep) / 10
  d$cluster = rep(1:7, length.out = n)
  d$y = 1 + .4*d$x + cos(d$g) + sin(d$h) + d$o + rnorm(n)
  d$x[1L] = NA
  d$w[2L] = 0
  d$id = seq_len(n)
  used = seq.int(3L, nrow(core))
  dense = lm(y ~ x + factor(g) + factor(h), data = d[used, ], weights = w, offset = o)
  X = model.matrix(dense)
  B = solve(crossprod(X, X*d$w[used]))
  scores = rowsum(X*(d$w[used]*residuals(dense)),d$cluster[used])
  clustered = B %*% crossprod(scores) %*% B
  for (strategy in c("demean", "compress")) {
    fit = dbreg(y ~ x | g + h, conn = conn, data = d, strategy = strategy,
      weights = ~w, offset = ~o, fixef_rm = "singletons", retain = "all")
    expect_equal(nobs(fit), length(used))
    expect_equal(fit$nobs_pruned, 5)
    expect_equal(fit$pruning_rounds, 5L)
    expect_identical(fit$fixef_rm, "singletons")
    expect_equal(coef(fit)[["x"]],coef(dense)[["x"]],tolerance=1e-9)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(dense)["x","x"]),tolerance=1e-9)
    expect_identical(case.names(fit), used)
    expect_equal(fitted(fit), unname(fitted(dense)), tolerance=1e-9)
    original = predict(fit, sample = "original")
    expect_identical(which(is.na(original)), setdiff(seq_len(n),used))
    expect_equal(original[used], unname(fitted(dense)), tolerance=1e-9)
    expect_equal(as.numeric(vcov(fit,cluster=~cluster,ssc=list(K.adj=FALSE,G.adj=FALSE))),
      as.numeric(clustered["x","x"]),tolerance=1e-9)
    close(fit)
    expect_identical(DBI::dbListTables(conn),character())
  }
  d$t = sin(d$id)/2
  varying_dense = lm(y ~ x + factor(g)*t + factor(h), data=d[used,], weights=w,offset=o)
  for (strategy in c("demean","compress")) {
    varying = dbreg(y ~ x | g[t] + h,conn=conn,data=d,strategy=strategy,
      weights=~w,offset=~o,fixef_rm="singletons",retain="all")
    expect_equal(nobs(varying),length(used))
    expect_equal(coef(varying)[["x"]],coef(varying_dense)[["x"]],tolerance=1e-8)
    expect_equal(fitted(varying),unname(fitted(varying_dense)),tolerance=1e-8)
    close(varying)
    expect_identical(DBI::dbListTables(conn),character())
  }
  default = dbreg(y ~ x | g + h, conn=conn,data=d,strategy="compress",weights=~w,offset=~o)
  expect_equal(nobs(default),n-2L)
  expect_identical(default$fixef_rm,"none")
  DBI::dbWriteTable(conn,"source",d,temporary=TRUE)
  fit = dbreg(y ~ x | g + h, conn=conn,table="source",strategy="compress",weights=~w,offset=~o,
    row_id="id",sample_table="kept",fixef_rm="singletons")
  expect_equal(DBI::dbGetQuery(conn,"SELECT id FROM kept ORDER BY id")$id,used)
  close(fit)
  DBI::dbExecute(conn,"DROP TABLE kept")
  expect_identical(DBI::dbListTables(conn),"source")
  expect_error(dbreg(y ~ x | g,conn=conn,data=d,fixef_rm="all"),"arg")
  expect_error(dbreg(y ~ x | g,conn=conn,data=d,fixef_rm="singletons",sql_only=TRUE),"requires execution")
  expect_error(dbreg(y ~ x | g,conn=conn,data=d,fixef_rm="singletons",drop_missings=FALSE),"requires execution")
  expect_error(dbreg(y ~ x | id,conn=conn,data=d,fixef_rm="singletons"),"No observations remain")
  expect_identical(DBI::dbListTables(conn),"source")
})
