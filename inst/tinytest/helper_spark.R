dbreg_test_spark_connection = function() {
  out = Sys.getenv("DBREG_SPARK_TEST_OUT", tempdir())
  config = sparklyr::spark_config()
  config[["sparklyr.shell.driver-memory"]] = "1G"
  config[["spark.sql.shuffle.partitions"]] = 2L
  config[["spark.default.parallelism"]] = 2L
  config[["spark.ui.enabled"]] = "false"
  config[["sparklyr.connect.enablehivesupport"]] = FALSE
  config[["spark.sql.warehouse.dir"]] = file.path(out, "warehouse")
  config[["spark.local.dir"]] = file.path(out, "local")
  # Select the actual installation: older sparklyr catalogs list 4.2 previews.
  conn = sparklyr::spark_connect(master = "local[2]", spark_home = Sys.getenv("SPARK_HOME", "/opt/spark"), config = config)
  version = as.character(DBI::dbGetQuery(conn, "SELECT version() AS version")[["version"]])
  if (!startsWith(version, "4.2.0 ")) {
    sparklyr::spark_disconnect(conn)
    stop("Connected engine is not the final Spark 4.2.0 release.", call. = FALSE)
  }
  writeLines(version, file.path(out, "spark-version.txt"))
  conn
}

dbreg_test_spark_algorithms = function(conn) {
  tables = function() sort(DBI::dbGetQuery(conn, "SHOW TABLES")[["tableName"]])
  before = tables()
  set.seed(812)
  n = 48L
  d = expand.grid(g = 1:3, h = 1:2, replicate = 1:8)
  d[["id"]] = seq_len(n)
  d[["w"]] = (d[["g"]]+1)*(d[["h"]]+1)
  d[["x"]] = rnorm(n)
  d[["z"]] = rnorm(n)
  d[["y"]] = .5*d[["x"]] -.3*d[["z"]] + sin(d[["g"]]) + rnorm(n)
  d[["x value"]] = d[["x"]]
  d[["instrument"]] = rnorm(n)
  d[["endog"]] = d[["instrument"]] + .2*d[["z"]] + rnorm(n)
  name = basename(tempfile("dbreg_spark_fixture_"))
  raw = paste0(name, "_raw")
  sparklyr::copy_to(conn, d, raw, overwrite = FALSE)
  # sparklyr normalizes R column names during upload. Create the nonstandard
  # SQL identifier explicitly so this tests the actual source schema.
  DBI::dbExecute(conn, paste("CREATE TEMP VIEW", name,
    "AS SELECT g,h,replicate,id,w,x,z,y,instrument,endog,x AS `x value`,CAST(50000+id AS INT) AS int_x,CAST(60000+id AS INT) AS int_z FROM", raw))
  on.exit({
    DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", name))
    DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", raw))
  })
  inventory = tables()
  expect_identical(dbreg:::detect_backend(conn)[["name"]], "spark")
  lookup = data.frame(id = 1:2, value = 1+c(0, .Machine[["double.eps"]]), label = c("quote'\\slash", "Łódź"))
  values = DBI::dbGetQuery(conn, paste("SELECT * FROM (", dbreg:::sql_values_frame(conn, lookup), ") lookup_values ORDER BY id"))
  expect_identical(values[["value"]], lookup[["value"]])
  expect_identical(values[["label"]], lookup[["label"]])
  exact_keys = DBI::dbGetQuery(conn, paste("SELECT COUNT(DISTINCT CAST(9007199254740992 AS BIGINT)+id) AS n FROM", name))
  expect_equal(as.numeric(exact_keys[["n"]]), n)
  for (link in c("logit", "probit", "cloglog")) {
    kernel = dbreg:::dbglm_family_sql(link, eta = "x")
    values = DBI::dbGetQuery(conn, paste("SELECT id,", kernel[["mu"]], "AS mu,", kernel[["derivative"]], "AS derivative FROM", name, "ORDER BY id"))
    family = stats::binomial(link)
    expect_true(max(abs(values[["mu"]]-family[["linkinv"]](d[["x"]]))) < 1e-12)
    expect_true(max(abs(values[["derivative"]]-family[["mu.eta"]](d[["x"]]))) < 1e-12)
  }
  # Actual Spark GROUPING SETS, JSON arrays and decimal SUM; no SQL emulator.
  moments = dbreg:::query_fe_moments(conn, name, c("g", "h"), c("x", "y"))
  expect_true(!is.null(moments))
  solver = dbreg:::fe_solver(moments[["sys"]], new.env(parent = emptyenv()))
  coefficients = solver[["solve"]](solver[["project"]](moments[["B"]]))
  indicators = cbind(stats::model.matrix(~factor(g)-1, d), stats::model.matrix(~factor(h)-1, d))
  projected = as.matrix(d[c("x", "y")]) - indicators %*% coefficients
  dummy = stats::model.matrix(~factor(g)+factor(h), d)
  reference = stats::lm.fit(dummy*sqrt(d[["w"]]), as.matrix(d[c("x", "y")])*sqrt(d[["w"]]))[["residuals"]]/sqrt(d[["w"]])
  expect_true(max(abs(projected-reference)) < 1e-9)
  expect_true(max(abs(crossprod(indicators, projected*d[["w"]]))) < 1e-8)

  for (formula in list(y ~ x + z, y ~ `x value` + z, y ~ x + z | g)) {
    actual = dbreg(formula, conn = conn, table = name, strategy = if (length(Formula::Formula(formula))[[2L]] > 1L) "demean" else "moments", weights = "w", vcov = "hc1")
    expected = fixest::feols(formula, d, weights = ~w, vcov = "hetero", fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
    expect_true(max(abs(coef(actual)-coef(expected))) < 1e-7)
    expect_true(max(abs(vcov(actual)-vcov(expected))) < 1e-7)
  }
  integer_fit = dbreg(y ~ int_x:int_z, conn = conn, table = name, strategy = "moments")
  integer_reference = stats::lm(y ~ I((50000+id)*(60000+id)), d)
  expect_equal(unname(coef(integer_fit)), unname(coef(integer_reference)), tolerance = 1e-6)
  expect_equal(as.numeric(vcov(integer_fit)), as.numeric(vcov(integer_reference)), tolerance = 1e-6)
  iv = dbivreg(y ~ z | endog ~ instrument, conn = conn, table = name, strategy = "moments")
  reference_iv = fixest::feols(y ~ z | endog ~ instrument, d, notes = FALSE)
  names = sub("^fit_", "", names(coef(reference_iv)))
  order = match(names(coef(iv)), names)
  expect_true(max(abs(coef(iv)-coef(reference_iv)[order])) < 1e-7)
  expect_true(max(abs(vcov(iv)-vcov(reference_iv)[order, order])) < 1e-7)

  # Tiny separably weighted panel converges quickly while exercising real
  # physical AP tables, window batching and owned-table cleanup.
  saved = options(dbreg.ap_solver = "iterative", dbreg.ap_method = "window", dbreg.ap_batch_size = 1L)
  on.exit(options(saved), add = TRUE)
  within = dbreg(y ~ x + z | g + h, conn = conn, table = name, strategy = "demean", weights = "w", fe_tol = 1e-8, fe_max_iter = 3L)
  reference_within = fixest::feols(y ~ x + z | g + h, d, weights = ~w,
    fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
  expect_true(max(abs(coef(within)-coef(reference_within))) < 1e-7)
  expect_true(max(abs(vcov(within)-vcov(reference_within))) < 1e-7)
  expect_identical(tables(), inventory)
  expect_error(dbreg(y ~ x + z | g + h, conn = conn, table = name, strategy = "demean", weights = "w", fe_tol = 1e-30, fe_max_iter = 1L), "did not converge")
  expect_identical(tables(), inventory)
  DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", name))
  DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", raw))
  expect_identical(tables(), before)
}
