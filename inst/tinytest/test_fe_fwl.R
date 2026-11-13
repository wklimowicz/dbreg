local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  saved = options(dbreg.ap_solver = "moments", dbreg.ap_method = NULL)
  on.exit(options(saved), add = TRUE)
  set.seed(11121)
  n = 500L
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(1:13, n, TRUE),
    h = sample(1:9, n, TRUE), k = sample(1:5, n, TRUE), w = runif(n, .2, 2), o = rnorm(n))
  d[["endog"]] = .8 * d[["z"]] + rnorm(n)
  d[["y"]] = d[["o"]] + .6 * d[["x"]] + .4 * d[["endog"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n)
  DBI::dbWriteTable(conn, "fwl_fixture", d)
  initial = sort(DBI::dbListTables(conn))
  # Both consumers share summaries and the same solver/certificate. Only the
  # sufficient-statistic consumer avoids transformed observation tables.
  for (K in 2:3) {
    inputs = list2env(list(conn = conn, from_statement = "FROM fwl_fixture",
      fe = c("g", "h", "k")[seq_len(K)], yvar = "y", weights = "w", cluster_var = NULL,
      verbose = FALSE, projection_limits = dbreg:::projection_limits()), parent = emptyenv())
    rows = dbreg:::dbreg_moment_projections(inputs, c("x", "z"), c("x", "z"), 100L, 1e-9)
    row_cells = rows[["transfer_cells"]]
    dbreg:::release_fit_resources(rows[["resources"]])
    fit = dbreg:::dbreg_moment_projections(inputs, c("x", "z"), c("x", "z"), 100L, 1e-9, materialize = FALSE)
    expect_identical(length(fit[["resources"]][["tables"]]), 1L)
    expect_identical(fit[["moment_method"]], "fwl")
    expect_identical(fit[["transformed_observation_tables"]], 0L)
    expect_equal(fit[["summary_rows"]], rows[["summary_rows"]] + 1)
    expect_equal(fit[["transfer_cells"]], row_cells + 8)
    # Independent weighted dummy incidence projection, no fixest dependency
    # for the compact cross-product oracle.
    D = model.matrix(reformulate(c("factor(g)", "factor(h)", if (K == 3) "factor(k)")), d)
    V = as.matrix(d[c("y", "x", "z")])
    sw = sqrt(d[["w"]] / max(d[["w"]]))
    residual = qr.resid(qr(D * sw), V * sw)
    ref = crossprod(residual)
    S = dbreg:::cross_moments(fit[["moments"]], c("x", "z"), "y", FALSE)
    expect_equal(unname(S[c("y", "x", "z"), c("y", "x", "z")]), unname(ref), tolerance = 1e-10)
    dbreg:::release_fit_resources(fit[["resources"]])
    expect_identical(sort(DBI::dbListTables(conn)), initial)
    low = inputs[["projection_limits"]]
    low[["transfer_cells"]] = row_cells + 1
    inputs[["projection_limits"]] = low
    attempt = dbreg:::projection_attempt("fwl", dbreg:::dbreg_moment_projections(inputs, c("x", "z"), c("x", "z"), 100L, 1e-9, materialize = FALSE))
    expect_identical(attempt[["reason"]], "transfer_limit")
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  # A changed caller source between summary and global reads cannot mix
  # samples: both stages must continue reading the owned snapshot.
  DBI::dbWriteTable(conn, "mutable_fwl", d)
  inputs = list2env(list(conn = conn, from_statement = "FROM mutable_fwl", fe = c("g", "h"),
    yvar = "y", weights = "w", cluster_var = NULL, verbose = FALSE,
    projection_limits = dbreg:::projection_limits()), parent = emptyenv())
  project = dbreg:::dbreg_moment_projections
  environment(project) = new.env(parent = environment(project))
  environment(project)[["query_fe_summaries"]] = function(...) {
    summaries = dbreg:::query_fe_summaries(...)
    DBI::dbExecute(conn, "UPDATE mutable_fwl SET y=y+100")
    summaries
  }
  fit = project(inputs, c("x", "z"), c("x", "z"), 100L, 1e-9, materialize = FALSE)
  D = model.matrix(~factor(g) + factor(h), d)
  sw = sqrt(d[["w"]] / max(d[["w"]]))
  ref = crossprod(qr.resid(qr(D * sw), as.matrix(d[c("y", "x", "z")]) * sw))
  S = dbreg:::cross_moments(fit[["moments"]], c("x", "z"), "y", FALSE)
  expect_equal(unname(S[c("y", "x", "z"), c("y", "x", "z")]), unname(ref), tolerance = 1e-10)
  dbreg:::release_fit_resources(fit[["resources"]])
  DBI::dbRemoveTable(conn, "mutable_fwl")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  rm("query_fe_summaries", envir = environment(project))
  environment(project)[["fe_fwl_information"]] = function(...) stop("unexpected FWL driver error")
  inputs[["from_statement"]] = "FROM fwl_fixture"
  expect_error(dbreg:::projection_attempt("fwl", project(inputs, c("x", "z"), c("x", "z"), 100L, 1e-9, materialize = FALSE)), "unexpected FWL driver error")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  adjusted = d
  adjusted[["y"]] = adjusted[["y"]] - adjusted[["o"]]
  for (K in 2:3) for (kind in list("iid", "hc1", ~g)) for (iv in c(FALSE, TRUE)) {
    fe = paste(c("g", "h", "k")[seq_len(K)], collapse = "+")
    fml = as.formula(paste0("y~x", if (!iv) "+z", "|", fe, if (iv) "|endog~z"))
    fun = if (iv) dbivreg else dbreg
    fit = fun(fml, conn = conn, table = "fwl_fixture", weights = ~w, offset = ~o,
      strategy = "demean", vcov = kind, fe_max_iter = 2L, fe_tol = 1e-9)
    ref = fixest::feols(fml, adjusted, weights = ~w,
      vcov = if (identical(kind, "hc1")) "hetero" else kind,
      ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
    ix = match(names(coef(fit)), sub("^fit_", "", names(coef(ref))))
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)[ix]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[ix, ix, drop = FALSE]), tolerance = 1e-8)
    expect_equal(nobs(fit), n)
    expect_identical(fit[["ap_moment_method"]], "fwl")
    expect_identical(fit[["ap_transformed_observation_tables"]], 0L)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  # Near absorption forces actual projected cross-products, rather than
  # trusting cancellation of two large sufficient-statistic matrices.
  for (amplitude in c(1e-4, 1e-6)) {
    near = d
    near[["x"]] = near[["g"]] + amplitude * near[["x"]]
    near[["y"]] = 2 * near[["x"]] + cos(near[["h"]]) + .01 * rnorm(n)
    fit = dbreg(y ~ x | g + h, data = near, weights = ~w, strategy = "demean", fe_tol = 1e-9)
    D = model.matrix(~factor(g) + factor(h), near)
    sw = sqrt(near[["w"]])
    design = qr(D * sw)
    xr = qr.resid(design, near[["x"]] * sw)
    yr = qr.resid(design, near[["y"]] * sw)
    beta = sum(xr * yr) / sum(xr * xr)
    variance = sum((yr - beta * xr)^2) / (n - design[["rank"]] - 1) / sum(xr * xr)
    if (amplitude == 1e-4 || identical(fit[["ap_method"]], "bounded_moments")) {
      expect_true(fit[["ap_fwl_recomputed"]])
    } else {
      # Extremely small residual variation may honestly decline the bounded
      # solve's unchanged relative certificate and use ordinary AP.
      expect_identical(fit[["projection_diagnostics"]][[1]][["reason"]], "certificate_failed")
    }
    expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-5)
    expect_equal(as.numeric(vcov(fit)), variance, tolerance = 1e-5)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }

})
