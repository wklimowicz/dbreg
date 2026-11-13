local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  saved = options(
    dbreg.ap_solver = "auto",
    dbreg.ap_method = NULL,
    dbreg.duckdb_ap_method = NULL
  )
  on.exit(options(saved), add = TRUE)
  set.seed(11126)
  n = 360L
  d = data.frame(
    id = seq_len(n),
    x = rnorm(n),
    z = rnorm(n),
    q = rnorm(n),
    g = sample(1:12, n, TRUE),
    h = sample(1:7, n, TRUE),
    w = runif(n, .2, 2),
    o = rnorm(n)
  )
  d[["y"]] = d[["o"]] +
    .7 * d[["x"]] +
    sin(d[["g"]]) +
    cos(d[["h"]]) +
    (.2 + d[["g"]] / 20) * d[["z"]] +
    d[["g"]] / 30 * d[["q"]] +
    rnorm(n)
  DBI::dbWriteTable(conn, "varying", d)
  initial = sort(DBI::dbListTables(conn))
  adjusted = d
  adjusted[["y"]] = d[["y"]] - d[["o"]]
  for (policy in c("auto", "iterative")) {
    for (kind in list("iid", "hc1", ~g)) {
      prior = options(dbreg.ap_solver = policy)
      fit = dbreg(
        y ~ x | g[z, q] + h,
        conn = conn,
        table = "varying",
        weights = ~w,
        offset = ~o,
        vcov = kind,
        strategy = "demean",
        fe_max_iter = 100L,
        fe_tol = 1e-9
      )
      ref = fixest::feols(
        y ~ x | g[z, q] + h,
        adjusted,
        weights = ~w,
        vcov = if (identical(kind, "hc1")) "hetero" else kind,
        ssc = fixest::ssc(K.fixef = "full"),
        fixef.rm = "none",
        fixef.tol = 1e-11,
        notes = FALSE
      )
      expect_equal(
        as.numeric(coef(fit)),
        as.numeric(coef(ref)),
        tolerance = 1e-8
      )
      expect_equal(
        as.numeric(vcov(fit)),
        as.numeric(vcov(ref)),
        tolerance = 1e-8
      )
      expect_identical(
        fit[["ap_method"]],
        if (policy == "auto") "coordinate_moments" else "coordinate_cg"
      )
      expect_true(
        fit[["ap_max_abs"]] <= 1e-9 && fit[["ap_max_relative"]] <= 1e-9
      )
      expect_identical(fit[["varying_slopes"]], list(c("z", "q"), character()))
      expect_identical(sort(DBI::dbListTables(conn)), initial)
      D = model.matrix(~ factor(g) + factor(g):z + factor(g):q + factor(h), d)
      sw = sqrt(d[["w"]])
      within = qr.resid(qr(D * sw), as.matrix(adjusted[c("x", "y")]) * sw) / sw
      beta = drop(
        crossprod(within[, 1] * sw, within[, 2] * sw) /
          crossprod(within[, 1] * sw)
      )
      expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
      expected = d[["y"]] - within[, 2] + within[, 1] * beta
      expect_equal(predict(fit, d), expected, tolerance = 1e-8)
      no_y = d[setdiff(names(d), "y")]
      expect_message(
        {
          pred = predict(fit, no_y)
        },
        "Outcome"
      )
      expect_equal(pred, within[, 1] * beta + d[["o"]], tolerance = 1e-8)
      options(prior)
    }
  }
  # A raw trend can also be a structural regressor or a native FE source.
  for (fml in list(
    y ~ x + z | g[z] + h,
    y ~ x | g[z] + g[q] + h,
    y ~ x | g[z, q] + h[z]
  )) {
    fit = dbreg(fml, data = d, weights = ~w, offset = ~o, fe_tol = 1e-9)
    ref = fixest::feols(
      fml,
      adjusted,
      weights = ~w,
      fixef.rm = "none",
      ssc = fixest::ssc(K.fixef = "full"),
      fixef.tol = 1e-11,
      notes = FALSE
    )
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  }
  # Exact null axes and duplicate trends are harmless to the actual projection.
  d[["copy"]] = 2 * d[["z"]]
  d[["constant"]] = 1
  ref = dbreg(
    y ~ x | g[z] + h,
    data = d,
    weights = ~w,
    offset = ~o,
    fe_tol = 1e-9
  )
  for (policy in c("auto", "iterative")) {
    prior = options(dbreg.ap_solver = policy)
    fit = dbreg(
      y ~ x | g[z, copy, constant] + h,
      data = d,
      weights = ~w,
      offset = ~o,
      fe_tol = 1e-9
    )
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
    expect_true(fit[["ap_max_relative"]] <= 1e-9)
    expect_identical(
      fit[["ap_method"]],
      if (policy == "auto") "coordinate_moments" else "coordinate_cg"
    )
    options(prior)
  }
  # Missing trend rows enter audit/filtering. Zero weights retain prediction
  # roles but are not used to estimate local coordinate coefficients.
  d[["z"]][1L] = NA_real_
  fit = dbreg(
    y ~ x | g[z] + h,
    data = d,
    weights = ~w,
    offset = ~o,
    fe_tol = 1e-9
  )
  expect_equal(nobs(fit), n - 1L)
  pred = predict(fit, d)
  expect_true(is.na(pred[1L]))
  expect_true(all(is.finite(pred[-1L])))
  zero = d[-1L, ]
  zero[["w"]][1L] = 0
  pred = predict(fit, zero)
  expect_true(all(is.finite(pred)))
  for (factor in c(1e-100, 1e100)) {
    scaled = adjusted
    scaled[["z"]] = scaled[["z"]] * factor
    scaled[["q"]] = scaled[["q"]] * factor
    scaled[["w"]] = scaled[["w"]] * factor
    fit = dbreg(
      y ~ x | g[z, q] + h,
      data = scaled,
      weights = ~w,
      fe_tol = 1e-9
    )
    ref = fixest::feols(
      y ~ x | g[z, q] + h,
      adjusted,
      weights = ~w,
      ssc = fixest::ssc(K.fixef = "full"),
      fixef.rm = "none",
      fixef.tol = 1e-11,
      notes = FALSE
    )
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  }
  # Frozen native tuple IDs, numeric producers and exported sample keys stay
  # separate, including source names used by internal state.
  DBI::dbExecute(
    conn,
    "CREATE TEMP TABLE typed_vary AS SELECT * EXCLUDE(g,h),CAST('9007199254740992' AS BIGINT)+CAST(g AS BIGINT) AS g,'a,b;|/'||CAST(h AS VARCHAR) AS h,z AS __w,q AS x_tilde FROM varying"
  )
  fit = dbreg(
    y ~ x | (g^h)[`__w`] + h[x_tilde],
    conn = conn,
    table = "typed_vary",
    weights = ~w,
    offset = ~o,
    fe_tol = 1e-9,
    row_id = c("id", "g", "h"),
    sample_table = "vary_keys"
  )
  ref = fixest::feols(
    y ~ x | (g^h)[z] + h[q],
    adjusted,
    weights = ~w,
    fixef.rm = "none",
    ssc = fixest::ssc(K.fixef = "full"),
    fixef.tol = 1e-11,
    notes = FALSE
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  keys = DBI::dbGetQuery(
    conn,
    "SELECT id,CAST(g AS VARCHAR) AS g,h FROM vary_keys ORDER BY id"
  )
  expect_identical(
    keys[["g"]],
    paste0("900719925474", sprintf("%04d", 992 + adjusted[["g"]]))
  )
  expect_identical(keys[["h"]], paste0("a,b;|/", adjusted[["h"]]))
  DBI::dbRemoveTable(conn, "vary_keys")
  DBI::dbRemoveTable(conn, "typed_vary")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # A small real within-group trend must survive a very large group offset.
  shifted = adjusted
  shifted[["z"]] = 1e12 + shifted[["z"]]
  fit = dbreg(y ~ x | g[z] + h, data = shifted, weights = ~w, fe_tol = 1e-9)
  D = model.matrix(~ factor(g) + factor(h) + factor(g):I(z - 1e12), shifted)
  sw = sqrt(shifted[["w"]])
  within = qr.resid(qr(D * sw), as.matrix(shifted[c("x", "y")]) * sw) / sw
  beta = drop(
    crossprod(within[, 1] * sw, within[, 2] * sw) / crossprod(within[, 1] * sw)
  )
  expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
  expected = shifted[["y"]]-within[,2]+within[,1]*beta
  expect_equal(predict(fit,shifted),expected,tolerance=1e-8)
  inputs=list2env(list(conn=conn,from_statement="FROM varying",fe=c("g","h"),fe_slopes=list(c("z","q"),character()),
    yvar="y",weights="w",cluster_var=NULL,verbose=FALSE,projection_limits=dbreg:::projection_limits()),parent=emptyenv())
  for(name in c("summary_rows","transfer_cells","design_entries","memory_bytes")) {
    limited=list2env(as.list(inputs),parent=emptyenv());limited[["projection_limits"]][[name]]=1
    attempt=dbreg:::projection_attempt("coordinates",dbreg:::dbreg_varying_projections(limited,"x","x",100L,1e-9))
    expect_identical(attempt[["status"]],"declined")
    expect_identical(sort(DBI::dbListTables(conn)),initial)
  }
  project=dbreg:::dbreg_varying_projections
  environment(project)=new.env(parent=environment(project))
  environment(project)[["fe_coordinate_certificate"]]=function(...)list(max_abs=Inf,max_relative=Inf)
  attempt=dbreg:::projection_attempt("coordinates",project(inputs,"x","x",100L,1e-9))
  expect_identical(attempt[["reason"]],"certificate_failed")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  environment(project)[["fe_coordinate_certificate"]]=function(...)stop("unexpected coordinate driver error")
  expect_error(dbreg:::projection_attempt("coordinates",project(inputs,"x","x",100L,1e-9)),"unexpected coordinate driver error")
  expect_identical(sort(DBI::dbListTables(conn)),initial)
  shifted_x=adjusted;shifted_x[["x"]]=1e12+shifted_x[["x"]]
  D=model.matrix(~factor(g)+factor(h)+factor(g):z+factor(g):q,shifted_x)
  sw=sqrt(shifted_x[["w"]]);within=qr.resid(qr(D*sw),cbind(shifted_x[["x"]]-1e12,shifted_x[["y"]])*sw)/sw
  beta=drop(crossprod(within[,1]*sw,within[,2]*sw)/crossprod(within[,1]*sw))
  fit=dbreg(y~x|g[z,q]+h,data=shifted_x,weights=~w,fe_tol=1e-9)
  expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
  # Canonical limits can decline moments and route to the same certified CG.
  prior = options(dbreg.fe_dense_max = 0, dbreg.fe_cross_max = 0)
  fit = dbreg(
    y ~ x | g[z, q] + h,
    conn = conn,
    table = "varying",
    weights = ~w,
    offset = ~o,
    fe_tol = 1e-9
  )
  expect_identical(fit[["ap_method"]], "coordinate_cg")
  expect_identical(
    vapply(fit[["projection_diagnostics"]], `[[`, character(1), "status"),
    c("declined", "completed")
  )
  options(prior)
  prior = options(dbreg.ap_solver = "iterative")
  expect_error(
    dbreg(
      y ~ x | g[z, q] + h,
      conn = conn,
      table = "varying",
      weights = ~w,
      fe_max_iter = 1L
    ),
    "iteration cap"
  )
  options(prior)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  expect_error(dbreg(y ~ x | g[z], data = d, sql_only = TRUE), "SQL-only")
  expect_error(
    dbreg(y ~ x | g[z], data = d, strategy = "mundlak"),
    "Varying slopes"
  )
  expect_error(
    dbivreg(y ~ x | g[z] | q ~ z, conn = conn, table = "varying"),
    "Varying slopes are not supported"
  )
  d[["text"]] = "bad"
  expect_error(dbreg(y ~ x | g[text], data = d), "must be numeric")
  expect_error(dbreg(y ~ x | g[log(z)], data = d), "Use group")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})


# A perfect fit cannot acquire a negative IID covariance from cancellation.
local({
  S=matrix(c(1,1,1,1-2*.Machine[["double.eps"]]),2L,dimnames=list(c("x","y"),c("x","y")))
  expect_equal(dbreg:::moment_rss(S,"y","x",1),0)
  S["y","y"]=-1
  expect_error(dbreg:::moment_rss(S,"y","x",1),"invalid negative")
})

# Declining coordinate_moments on transfer_limit must not starve coordinate_cg.
local({
  set.seed(20261119)
  n = 200L
  d = data.frame(
    g1 = rep(1:20, each = 10),
    g2 = sample(1:15, n, replace = TRUE),
    t = rnorm(n),
    x = rnorm(n)
  )
  d[["y"]] = 1 + 0.5 * d[["x"]] + 0.2 * d[["g1"]] + 0.1 * d[["g1"]] * d[["t"]] - 0.3 * d[["g2"]] + rnorm(n, sd = 0.1)
  old = options(dbreg.projection_max_rows = 80L)
  on.exit(options(old), add = TRUE)
  fit_auto = dbreg(y ~ x | g1[t] + g2, data = d, strategy = "demean")
  old2 = options(dbreg.ap_solver = "iterative")
  on.exit(options(old2), add = TRUE)
  fit_iter = dbreg(y ~ x | g1[t] + g2, data = d, strategy = "demean")
  expect_equal(as.numeric(coef(fit_auto)), as.numeric(coef(fit_iter)), tolerance = 1e-8)
})
