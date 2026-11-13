local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11097)
  n = 160L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n),
    g = rep(1:8, 20), cl = rep(1:10, 16), w = runif(n, .2, 2),
    keep = rep(c(TRUE, FALSE, NA, TRUE), 40), u = rep(c(NA, NaN, -2, 1, 3), 32),
    f = rep(c("NaN", "a'b", "c", "d"), 40))
  d[["endog"]] = .8*d[["z"]] + rnorm(n)
  d[["y"]] = .5*d[["x"]] + d[["endog"]] + sin(d[["g"]]) + rnorm(n)
  DBI::dbWriteTable(conn, "sample_filter_source", d)
  # R uploads may replace NaN with SQL NULL, so introduce real NaNs in SQL.
  DBI::dbExecute(conn, "UPDATE sample_filter_source SET u=CAST('NaN' AS DOUBLE) WHERE id%5=2")
  schema = DBI::dbGetQuery(conn, "SELECT * FROM sample_filter_source LIMIT 0")
  predicates = list(~u>0, ~u %in% c(NA, 1), ~u %in% c(NaN, 1), ~is.na(u),
    ~f %in% c("a'b", "NaN"), ~!(f=="c"), ~keep==TRUE,
    ~keep & (x + 2*z > -1), ~(keep | is.na(u)) & !(x^2 <= .2),
    ~u %in% c(), ~u %in% NULL, ~u %in% c(-2, +3), ~-x > 0, ~+x < 0,
    ~x/2 - z*3 >= 0, ~keep != FALSE, ~x < Inf)
  for (predicate in predicates) {
    sql = dbreg:::sql_subset_predicate(conn, predicate, schema)
    actual = DBI::dbGetQuery(conn, paste("SELECT id FROM sample_filter_source WHERE", sql))[["id"]]
    selected = eval(predicate[[2L]], d, baseenv())
    selected[is.na(selected)] = FALSE
    expect_equal(sort(actual), d[["id"]][selected], info = paste(deparse(predicate), collapse = ""))
  }
  selected = !is.na(d[["keep"]]) & d[["keep"]] & d[["x"]] > -.7
  for (strategy in c("moments", "compress", "demean", "mundlak")) {
    formula = if (strategy %in% c("demean", "mundlak")) y~x|g else y~x
    for (sql in c(FALSE, TRUE)) {
      args = if (sql) list(conn = conn, table = "sample_filter_source") else list(data = d)
      fit = do.call(dbreg, c(list(fml = formula, subset = ~keep & x > -.7,
        weights = "w", vcov = ~cl, strategy = strategy), args))
      # The same strategy on R-selected rows also covers Mundlak's CRE model.
      reference = dbreg(formula, data = d[selected, ], weights = "w", vcov = ~cl, strategy = strategy)
      expect_equal(nobs(fit), sum(selected))
      expect_equal(coef(fit, fe = TRUE), coef(reference, fe = TRUE), tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit, fe = TRUE)), as.numeric(vcov(reference, fe = TRUE)), tolerance = 1e-9)
      if (strategy == "moments") {
        independent = lm(y~x, d[selected, ], weights = w)
        expect_equal(as.numeric(coef(fit)), as.numeric(coef(independent)), tolerance = 1e-9)
      }
    }
  }
  for (strategy in c("moments", "compress", "demean")) {
    formula = if (strategy == "demean") y~x|g|endog~z else y~x|endog~z
    fit = dbivreg(formula, conn = conn, table = "sample_filter_source", subset = ~keep & x > -.7,
      weights = "w", vcov = ~cl, strategy = strategy)
    reference = fixest::feols(formula, d[selected, ], weights = ~w, vcov = ~cl,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
    names_reference = sub("^fit_", "", names(coef(reference)))
    expect_equal(nobs(fit), sum(selected))
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(reference)[match(names(coef(fit)), names_reference)]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)[match(names(coef(fit)), names_reference), match(names(coef(fit)), names_reference)]), tolerance = 1e-7)
  }
  # Bare logical column selectors, non-syntactic source names and environments.
  unsafe = d; names(unsafe)[names(unsafe)=="keep"] = "keep row"
  fit = dbreg(y~x, data = unsafe, subset = "keep row", strategy = "moments")
  reference = lm(y~x, d[!is.na(d[["keep"]]) & d[["keep"]], ])
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(reference)), tolerance = 1e-9)
  expect_equal(coef(dbreg(y~x, data = unsafe, subset = ~`keep row`, strategy = "moments")), coef(fit), tolerance = 1e-9)
  for (subset in list(1:3, ~log(x)>0, ~x, ~unknown>0, ~keep & x,
      ~f>2, ~x %in% c("a"), ~x %in% seq(1,3), ~system("false"))) {
    expect_error(dbreg(y~x, data = d, subset = subset, strategy = "moments"), info = paste(deparse(subset), collapse = ""))
    expect_error(dbivreg(y~x|endog~z, data = d, subset = subset, strategy = "moments"))
  }
  # Source weights remain validated even if a predicate excludes the bad row.
  bad = d; bad[["w"]][1] = -1
  expect_error(dbreg(y~x, data = bad, weights = "w", subset = ~id>1), "non-negative")
  tables = DBI::dbListTables(conn)
  expect_error(dbreg(y~x, conn = conn, data = d, subset = ~unknown>0))
  expect_equal(DBI::dbListTables(conn), tables)
  expect_true(DBI::dbIsValid(conn))
  # Existing source clauses and inspection modes use the same selected sample.
  table = "(SELECT * FROM sample_filter_source WHERE id>10 ORDER BY id) s"
  fit = dbreg(y~x, conn = conn, table = table, subset = ~keep, strategy = "moments")
  expect_equal(nobs(fit), sum(d[["id"]]>10 & !is.na(d[["keep"]]) & d[["keep"]]))
  query = dbreg(y~x, conn = conn, table = table, subset = ~keep, strategy = "moments", sql_only = TRUE)
  moments = dbreg(y~x, conn = conn, table = table, subset = ~keep, strategy = "moments", data_only = TRUE)
  expect_equal(DBI::dbGetQuery(conn, query), moments, tolerance = 1e-12)
})

local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11098); n = 120L
  d = data.frame(id=seq_len(n), x=rnorm(n), z=rnorm(n), g=as.numeric(rep(1:6,20)),
    cl=as.numeric(rep(1:12,10)), w=runif(n,.2,2), f=rep(c("NaN","a","b"),40))
  d[["endog"]] = d[["z"]] + rnorm(n); d[["y"]] = d[["x"]] + d[["endog"]] + rnorm(n)
  for (column in c("x", "y", "cl", "w", "z", "endog")) {
    DBI::dbWriteTable(conn, "nan_source", d, overwrite = TRUE)
    DBI::dbExecute(conn, paste0("UPDATE nan_source SET ", column, "=CAST('NaN' AS DOUBLE) WHERE id=1"))
    formula = y~x+f
    fit = dbreg(formula, conn=conn, table="nan_source", weights="w", vcov=~cl, strategy="moments")
    required = column %in% c("x","y","cl","w")
    reference = dbreg(formula, data=d[if(required) -1 else seq_len(n), ], weights="w", vcov=~cl, strategy="moments")
    expect_equal(nobs(fit), nobs(reference))
    expect_equal(coef(fit), coef(reference), tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance=1e-8)
    iv = dbivreg(y~x|endog~z, conn=conn, table="nan_source", weights="w", vcov=~cl, strategy="moments")
    reference_iv = dbivreg(y~x|endog~z, data=d[-1, ], weights="w", vcov=~cl, strategy="moments")
    expect_equal(nobs(iv), n-1)
    expect_equal(coef(iv), coef(reference_iv), tolerance=1e-8)
  }
  DBI::dbWriteTable(conn, "nan_source", d, overwrite=TRUE)
  DBI::dbExecute(conn, "UPDATE nan_source SET g=CAST('NaN' AS DOUBLE) WHERE id=1")
  for (strategy in c("demean", "compress", "mundlak")) {
    fit = dbreg(y~x+f|g, conn=conn, table="nan_source", strategy=strategy)
    reference = dbreg(y~x+f|g, data=d[-1, ], strategy=strategy)
    expect_equal(nobs(fit), n-1)
    expect_equal(coef(fit), coef(reference), tolerance=1e-8)
  }
})


local({
  expect_true(dbreg:::databricks_metadata(list(dbms.name="Databricks")))
  expect_true(dbreg:::databricks_metadata(list(driver.name="Databricks ODBC")))
  expect_true(dbreg:::databricks_metadata(connection_string=list(HOST="example.azuredatabricks.net")))
  expect_true(dbreg:::databricks_metadata(connection_string=list(driver="Databricks")))
  expect_false(dbreg:::databricks_metadata(list(dbms.name="Spark SQL"), list(host="spark.example.org")))
  expect_false(dbreg:::databricks_metadata(list(dbms.name="Spark SQL"), list(host="spark.example.org", pwd="databricks")))
})


local({
  class = "dbreg_sample_backend_connection"
  methods::setClass(class, contains="DBIConnection", slots=c(metadata="list"))
  methods::setMethod(DBI::dbGetInfo, class, function(dbObj, ...) dbObj@metadata)
  on.exit({methods::removeMethod(DBI::dbGetInfo, class); methods::removeClass(class)})
  for (case in list(list(dbms.name="Spark SQL", driver.name="Databricks ODBC"),
      list(dbms.name="Spark SQL"), list(dbms.name="Microsoft SQL Server"))) {
    connection = methods::new(class, metadata=case)
    backend = dbreg:::detect_backend(connection)
    expect_identical(backend[["name"]], if (case[["dbms.name"]]=="Microsoft SQL Server") "sqlserver" else "spark")
    expect_identical(backend[["is_databricks"]], identical(case[["driver.name"]], "Databricks ODBC"))
  }
})

## Integer ranges are literal %in% values.
local({
  set.seed(5)
  d = data.frame(y = rnorm(300), x = rnorm(300), year = sample(2005:2020, 300, TRUE))
  ref = stats::lm(y ~ x, d, subset = year %in% 2010:2015)
  fit = dbreg(y ~ x, data = d, subset = ~ year %in% 2010:2015)
  expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-10)
  expect_equal(nobs(fit), nobs(ref))
  fit = dbreg(y ~ x, data = d, subset = ~ year %in% c(2010:2012, 2015))
  expect_equal(nobs(fit), sum(d$year %in% c(2010:2012, 2015)))
  expect_error(dbreg(y ~ x, data = d, subset = ~ year %in% 1:1e6), "literal")
  # Exact double literals in [1e-4, 1) must not round via DuckDB DECIMAL parsing.
  set.seed(2026)
  k = min(runif(120))
  d_sub = data.frame(y = c(1, 2, 3, 4, 5, 6), x = c(0, 1, 2, 3, 0, 1), v = c(k, k, k, k, 1.0, 1.0))
  fit = dbreg(y ~ x, data = d_sub, subset = stats::as.formula(substitute(~ v == val, list(val = k))))
  expect_equal(nobs(fit), 4)
})
