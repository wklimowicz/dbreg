local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  saved = options(dbreg.ap_solver = "auto", dbreg.ap_method = NULL,
    dbreg.fe_dense_max = 0, dbreg.fe_cross_max = 0,
    dbreg.cg_max_iter = 1000L, dbreg.cg_restart = 500L)
  on.exit(options(saved), add = TRUE)
  set.seed(11124)
  n = 360L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = sample(1:12, n, TRUE),
    h = sample(1:9, n, TRUE), k = sample(1:7, n, TRUE), w = runif(n, .2, 2), o = rnorm(n))
  d[["endog"]] = .8 * d[["z"]] + rnorm(n)
  d[["y"]] = d[["o"]] + .7 * d[["x"]] + .4 * d[["endog"]] + sin(d[["g"]]) + cos(d[["h"]]) + d[["k"]] / 8 + rnorm(n)
  DBI::dbWriteTable(conn, "cg_fixture", d)
  initial = sort(DBI::dbListTables(conn))
  adjusted = d;adjusted[["y"]] = d[["y"]] - d[["o"]]
  for (iv in c(FALSE, TRUE)) for (kind in list("iid", "hc1", ~g)) {
    fml = if (iv) y ~ x | g + h + k | endog ~ z else y ~ x + z | g + h + k
    fun = if (iv) dbivreg else dbreg
    fit = fun(fml, conn = conn, table = "cg_fixture", weights = ~w, offset = ~o,
      vcov = kind, strategy = "demean", fe_max_iter = 1L, fe_tol = 1e-9)
    ref = fixest::feols(fml, adjusted, weights = ~w, vcov = if (identical(kind, "hc1")) "hetero" else kind,
      ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
    ix = match(names(coef(fit)), sub("^fit_", "", names(coef(ref))))
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)[ix]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[ix, ix, drop = FALSE]), tolerance = 1e-8)
    expect_identical(fit[["ap_method"]], "sql_cg")
    expect_identical(fit[["ap_preconditioner"]], "diagonal")
    expect_true(fit[["ap_solver_iterations"]] > 1L)
    expect_true(fit[["ap_peak_tables"]] <= 12L)
    expect_true(fit[["ap_max_abs"]] <= 1e-9 && fit[["ap_max_relative"]] <= 1e-9)
    expect_identical(vapply(fit[["projection_diagnostics"]], `[[`, character(1), "status"), c("declined", "completed"))
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  # Compare actual projected rows to an independent weighted dummy QR.
  inputs = list2env(list(conn = conn, from_statement = "FROM cg_fixture", fe = c("g", "h", "k"),
    yvar = "y", weights = "w", cluster_var = "id", verbose = FALSE,
    projection_limits = dbreg:::projection_limits()), parent = emptyenv())
  fit = dbreg:::dbreg_cg_projections(inputs, c("x", "z"), c("x", "z"), 1e-9)
  actual = DBI::dbGetQuery(conn, paste("SELECT y_tilde,x_tilde,z_tilde FROM", fit[["table"]], "ORDER BY id"))
  D = model.matrix(~factor(g) + factor(h) + factor(k), d)
  sw = sqrt(d[["w"]])
  ref = qr.resid(qr(D * sw), as.matrix(d[c("y", "x", "z")]) * sw) / sw
  # Explicit native passthrough IDs establish row identity for this oracle.
  expect_equal(unname(as.matrix(actual)), unname(ref), tolerance = 1e-8)
  expect_equal(unname(crossprod(as.matrix(actual), d[["w"]] * as.matrix(actual))), unname(crossprod(ref, d[["w"]] * ref)), tolerance = 1e-8)
  # All source/encoded/final table ownership is explicit and releasable.
  expect_true(length(fit[["resources"]][["tables"]]) <= 5L)
  dbreg:::release_fit_resources(fit[["resources"]])
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  prior = options(dbreg.cg_restart = 4L)
  restarted = dbreg(y ~ x + z | g + h + k, conn = conn, table = "cg_fixture", weights = ~w, offset = ~o, strategy = "demean", fe_max_iter = 1L, fe_tol = 1e-9)
  options(prior)
  ref = fixest::feols(y ~ x + z | g + h + k, adjusted, weights = ~w, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
  expect_identical(restarted[["ap_method"]], "sql_cg")
  expect_equal(as.numeric(coef(restarted)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_true(restarted[["ap_peak_tables"]] <= 12L)
  expect_true(restarted[["ap_cg_restarts"]] > 0L)
  # Distinct native keys and exported rows survive coefficient-space state.
  DBI::dbExecute(conn, "CREATE TEMP TABLE typed_cg AS SELECT * EXCLUDE(g,h),CAST('9007199254740992' AS BIGINT)+CAST(g AS BIGINT) AS g,'a|/'||CAST(h AS VARCHAR) AS h FROM cg_fixture")
  typed = dbreg(y ~ x + z | g + h + k, conn = conn, table = "typed_cg", weights = ~w, offset = ~o,
    strategy = "demean", fe_max_iter = 1L, fe_tol = 1e-9, row_id = c("id", "g", "h"), sample_table = "cg_keys")
  expect_identical(typed[["ap_method"]], "sql_cg")
  expect_equal(as.numeric(coef(typed)), as.numeric(coef(ref)), tolerance = 1e-8)
  keys = DBI::dbGetQuery(conn, "SELECT id,CAST(g AS VARCHAR) AS g,h FROM cg_keys ORDER BY id")
  expect_identical(keys[["g"]], paste0("900719925474", sprintf("%04d", 992 + d[["g"]])))
  expect_identical(keys[["h"]], paste0("a|/", d[["h"]]))
  DBI::dbRemoveTable(conn, "cg_keys");DBI::dbRemoveTable(conn, "typed_cg")
  for (factor in c(1e-100, 1e100)) {
    DBI::dbExecute(conn, paste0("CREATE TEMP TABLE cg_scaled AS SELECT * EXCLUDE(w),w*", sprintf("%.17g", factor), " AS w FROM cg_fixture"))
    scaled = dbreg(y ~ x + z | g + h + k, conn = conn, table = "cg_scaled", weights = ~w, offset = ~o, strategy = "demean", fe_max_iter = 1L, fe_tol = 1e-9)
    expect_identical(scaled[["ap_method"]], "sql_cg")
    expect_equal(as.numeric(coef(scaled)), as.numeric(coef(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(scaled)), as.numeric(vcov(ref)), tolerance = 1e-8)
    DBI::dbRemoveTable(conn, "cg_scaled")
  }
  separate = d
  block = seq_len(n) %% 2L
  separate[["g"]] = separate[["g"]] + 12 * block
  separate[["h"]] = separate[["h"]] + 9 * block
  separate[["k"]] = separate[["k"]] + 7 * block
  disconnected = dbreg(y ~ x + z | g + h + k, data = separate, weights = ~w, offset = ~o, strategy = "demean", fe_max_iter = 1L, fe_tol = 1e-9)
  separate[["y"]] = separate[["y"]] - separate[["o"]]
  separate_ref = fixest::feols(y ~ x + z | g + h + k, separate, weights = ~w, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
  expect_identical(disconnected[["ap_method"]], "sql_cg")
  expect_equal(as.numeric(coef(disconnected)), as.numeric(coef(separate_ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(disconnected)), as.numeric(vcov(separate_ref)), tolerance = 1e-8)
  # CG does not use grouping-ID bit masks and must cover >30 FE families.
  count = 256L
  many = data.frame(x = rnorm(count), w = runif(count, .2, 2))
  families = paste0("family", seq_len(31L))
  for (name in families) many[[name]] = sample(1:2, count, TRUE)
  many[["y"]] = 1.3 * many[["x"]] + as.numeric(as.matrix(many[families]) %*% seq(.01, .31, .01)) + rnorm(count)
  fml = as.formula(paste("y~x|", paste(families, collapse = "+")))
  fit = dbreg(fml, data = many, weights = ~w, strategy = "demean", fe_max_iter = 1L, fe_tol = 1e-9)
  ref_many = fixest::feols(fml, many, weights = ~w, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
  expect_identical(fit[["ap_method"]], "sql_cg")
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref_many)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref_many)), tolerance = 1e-8)
  # CG caps and typed capacity declines cannot leak state or claim convergence.
  prior = options(dbreg.cg_max_iter = 1L)
  attempt = dbreg:::projection_attempt("cg", dbreg:::dbreg_cg_projections(inputs, c("x", "z"), c("x", "z"), 1e-9))
  expect_identical(attempt[["reason"]], "nonconverged")
  options(prior)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  limited = list2env(as.list(inputs), parent = emptyenv())
  limited[["projection_limits"]][["transfer_cells"]] = 1
  attempt = dbreg:::projection_attempt("cg", dbreg:::dbreg_cg_projections(limited, c("x", "z"), c("x", "z"), 1e-9))
  expect_identical(attempt[["reason"]], "transfer_limit")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # Unexpected database failures propagate; final row certificates are independent.
  project = dbreg:::dbreg_cg_projections
  environment(project) = new.env(parent = environment(project))
  environment(project)[["fe_projection_certificate"]] = function(...) list(max_abs = Inf, max_relative = Inf)
  attempt = dbreg:::projection_attempt("cg", project(inputs, c("x", "z"), c("x", "z"), 1e-9))
  expect_identical(attempt[["reason"]], "certificate_failed")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  environment(project)[["fe_projection_certificate"]] = function(...) stop("unexpected CG driver error")
  expect_error(dbreg:::projection_attempt("cg", project(inputs, c("x", "z"), c("x", "z"), 1e-9)), "unexpected CG driver error")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # A finite but vanishing curvature must produce a typed numerical decline,
  # never an Inf SQL literal or an unexpected parser failure.
  kernel = dbreg:::fe_sql_cg
  environment(kernel) = new.env(parent = environment(kernel))
  environment(kernel)[["fe_cg_fetch"]] = function(conn, sql, width, budget, limits) {
    value = dbreg:::fe_cg_fetch(conn, sql, width, budget, limits)
    if (all(c("rr1", "da1") %in% names(value))) {
      value[["rr1"]] = 1
      value[["da1"]] = .Machine[["double.xmin"]] * .Machine[["double.eps"]]
    }
    value
  }
  project = dbreg:::dbreg_cg_projections
  environment(project) = new.env(parent = environment(project))
  environment(project)[["fe_sql_cg"]] = kernel
  attempt = dbreg:::projection_attempt("cg", project(inputs, c("x", "z"), c("x", "z"), 1e-9))
  expect_identical(attempt[["reason"]], "singular")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  manager_resources = dbreg:::fit_resources();manager_resources[["conn"]] = conn
  manager = dbreg:::fe_cg_manager(conn, manager_resources)
  expect_error(manager[["drop"]](as.character(DBI::dbQuoteIdentifier(conn, "cg_fixture"))), "borrowed table")
  for (bad in list(0, -1, NA_real_, Inf, "many", c(1, 2))) {
    prior = options(dbreg.cg_restart = bad)
    expect_error(dbreg:::fe_cg_control(), "positive whole")
    options(prior)
  }
  # Explicit iterative/moments choices keep their existing backend policy.
  for (choice in list(list(dbreg.ap_solver = "iterative"), list(dbreg.ap_solver = "moments"), list(dbreg.ap_method = "window"))) {
    prior = options(choice)
    expect_error(dbreg(y ~ x | g + h + k, conn = conn, table = "cg_fixture", weights = ~w, strategy = "demean", fe_max_iter = 1L), "within 1 iterations")
    options(prior)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})
