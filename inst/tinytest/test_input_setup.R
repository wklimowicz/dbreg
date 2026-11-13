local({
  # Setup owns resources until a successful return transfers them to its caller.
  # Inject failure after registration, exercising the actual native relation.
  setup_env = environment(dbreg:::setup_db_connection)
  original = get("validate_integer64_import", setup_env)
  state = new.env(parent = emptyenv())
  replacement = function(conn, table, columns) {
    state[["conn"]] = conn
    state[["table"]] = table
    state[["count"]] = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", table))[[1L]]
    stop("injected import verification failure")
  }
  assignInNamespace("validate_integer64_import", replacement, "dbreg")
  on.exit(assignInNamespace("validate_integer64_import", original, "dbreg"), add = TRUE)
  if (requireNamespace("bit64", quietly = TRUE)) {
    data = data.frame(x = bit64::as.integer64(1:4))
    expect_error(dbreg:::setup_db_connection(NULL, NULL, data, NULL), "injected import")
    expect_equal(as.numeric(state[["count"]]), 4)
    expect_false(DBI::dbIsValid(state[["conn"]]))
    conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    DBI::dbWriteTable(conn, "caller_table", data.frame(z = 1:3))
    before = DBI::dbListTables(conn)
    expect_error(dbreg:::setup_db_connection(conn, NULL, data, NULL), "injected import")
    expect_true(DBI::dbIsValid(conn))
    expect_identical(DBI::dbListTables(conn), before)
    expect_equal(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM caller_table")[[1L]], 3)
  }
  assignInNamespace("validate_integer64_import", original, "dbreg")
})

if (requireNamespace("bit64", quietly = TRUE)) local({
  set.seed(11051)
  n = 120L
  labels = bit64::as.integer64(c("9007199254740992", "9007199254740993", "9007199254740994", "9007199254740995"))
  g = rep(1:4, length.out = n)
  d = data.frame(g = labels[g], q = bit64::as.integer64(rep(-3:3, length.out = n)),
    x = rnorm(n), z = rnorm(n), w = bit64::as.integer64(rep(1:3, length.out = n)), cl = labels[g])
  d[["endog"]] = .7*d[["z"]] + rnorm(n)
  d[["y"]] = .4+.6*as.numeric(d[["q"]])+.3*d[["x"]]+sin(g)+rnorm(n)
  setup = dbreg:::setup_db_connection(NULL, NULL, d, NULL)
  on.exit(dbreg:::release_db_setup(setup[["conn"]], setup[["own_conn"]], setup[["registered_table"]]), add = TRUE)
  imported = DBI::dbGetQuery(setup[["conn"]], paste("SELECT g FROM", setup[["table_name"]]))[["g"]]
  expect_identical(as.character(imported), as.character(d[["g"]]))
  expect_true(inherits(imported, "integer64"))
  # The verification rejects incompatible and missing native schema entries.
  DBI::dbWriteTable(setup[["conn"]], "wrong_type", data.frame(g = as.numeric(1:3)))
  expect_error(dbreg:::validate_integer64_import(setup[["conn"]], "wrong_type", "g"), "did not preserve integer64")
  expect_error(dbreg:::validate_integer64_import(setup[["conn"]], "wrong_type", "missing"), "did not preserve integer64")
  reference = d
  for (name in c("q", "w")) reference[[name]] = as.numeric(reference[[name]])
  reference[["g"]] = as.character(d[["g"]])
  reference[["cl"]] = as.character(d[["cl"]])
  for (strategy in c("moments", "compress", "demean")) {
    formula = if (strategy == "demean") y~q+x|g else y~q+x
    oracle = fixest::feols(formula, reference, weights = ~w, vcov = ~cl,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
    fit = dbreg(formula, data = d, weights = "w", vcov = ~cl, strategy = strategy)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(oracle)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(oracle)), tolerance = 1e-8)
    if (strategy != "demean") expect_equal(predict(fit, d), as.numeric(predict(oracle, reference)), tolerance = 1e-8)
    if (strategy == "compress") expect_equal(sort(predict(fit)), sort(as.numeric(predict(oracle, reference))), tolerance = 1e-8)
  }
  for (mode in c("numeric", "integer64")) {
    supplied = DBI::dbConnect(duckdb::duckdb(bigint = mode, config = list(threads = "2")))
    tryCatch({
      # DuckDB >= 1.5.2 registers integer64 as DOUBLE when bigint = "numeric",
      # which validate_integer64_import() must reject; native SQL BIGINT tables
      # still work under both read modes.
      preserves_i64 = tryCatch({
        duckdb::duckdb_register(supplied, "__check_i64", data.frame(g = labels[1:2]))
        on.exit(duckdb::duckdb_unregister(supplied, "__check_i64"), add = TRUE)
        identical(DBI::dbGetQuery(supplied, "DESCRIBE __check_i64")[["column_type"]], "BIGINT")
      }, finally = try(duckdb::duckdb_unregister(supplied, "__check_i64"), silent = TRUE))
      if (preserves_i64) {
        for (strategy in c("demean", "compress")) {
          fit = dbreg(y~q+x|g, conn = supplied, data = d, weights = "w", vcov = ~cl, strategy = strategy)
          oracle = fixest::feols(y~q+x|g, reference, weights = ~w, vcov = ~cl,
            fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
          expect_equal(as.numeric(coef(fit)), as.numeric(coef(oracle)), tolerance = 1e-8)
          expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(oracle)), tolerance = 1e-8)
        }
        fit = dbreg(y~q+x+i(g), conn = supplied, data = d, weights = "w", strategy = "moments")
        oracle = fixest::feols(y~q+x|g, reference, weights = ~w, fixef.rm = "none")
        expect_equal(unname(coef(fit)[c("q", "x")]), unname(coef(oracle)), tolerance = 1e-8)
      } else {
        expect_error(dbreg(y~q+x|g, conn = supplied, data = d, weights = "w", vcov = ~cl, strategy = "demean"), "did not preserve integer64")
        values = sprintf("(CAST('%s' AS BIGINT),CAST('%s' AS BIGINT),%.17e,CAST('%s' AS BIGINT),CAST('%s' AS BIGINT),%.17e)",
          as.character(d[["g"]]), as.character(d[["q"]]), d[["x"]], as.character(d[["w"]]), as.character(d[["cl"]]), d[["y"]])
        DBI::dbExecute(supplied, paste0("CREATE TEMP TABLE native_i64 AS SELECT * FROM (VALUES ",
          paste(values, collapse = ","), ") v(g, q, x, w, cl, y)"))
        for (strategy in c("demean", "compress")) {
          fit = dbreg(y~q+x|g, conn = supplied, table = "native_i64", weights = "w", vcov = ~cl, strategy = strategy)
          oracle = fixest::feols(y~q+x|g, reference, weights = ~w, vcov = ~cl,
            fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
          expect_equal(as.numeric(coef(fit)), as.numeric(coef(oracle)), tolerance = 1e-8)
          expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(oracle)), tolerance = 1e-8)
        }
        fit = dbreg(y~q+x+i(g), conn = supplied, table = "native_i64", weights = "w", strategy = "moments")
        oracle = fixest::feols(y~q+x|g, reference, weights = ~w, fixef.rm = "none")
        expect_equal(unname(coef(fit)[c("q", "x")]), unname(coef(oracle)), tolerance = 1e-8)
        DBI::dbExecute(supplied, "DROP TABLE native_i64")
      }
      raw = DBI::dbGetQuery(supplied, "SELECT CAST(9007199254740993 AS BIGINT) AS k")[["k"]]
      expect_identical(inherits(raw, "integer64"), mode == "integer64")
      exact = dbreg:::read_keyed_query(supplied,
        "SELECT * FROM (VALUES (9007199254740992::BIGINT), (9007199254740993::BIGINT), (NULL::BIGINT), (-9007199254740993::BIGINT)) v(k)", "k")[["k"]]
      expect_identical(sort(as.character(exact), na.last = TRUE), sort(c("9007199254740992", "9007199254740993", NA_character_, "-9007199254740993"), na.last = TRUE))
      expect_identical(DBI::dbListTables(supplied), character())
      expect_true(DBI::dbIsValid(supplied))
    }, finally = DBI::dbDisconnect(supplied, shutdown = TRUE))
  }
  # The default BIGINT read mode also reaches SQL i() literal generation.
  X = model.matrix(~q+x+factor(g), reference)
  beta = as.numeric(qr.solve(X*sqrt(reference[["w"]]),d[["y"]]*sqrt(reference[["w"]])))
  for (strategy in c("moments", "compress")) {
    fit = dbreg(y~q+x+i(g), data = d, weights = "w", strategy = strategy)
    expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
    expect_equal(predict(fit,d), as.numeric(X%*%beta), tolerance = 1e-8)
    iv = dbivreg(y~q|endog~z, data = d, weights = "w", strategy = strategy)
    oracle = fixest::feols(y~q|endog~z, reference, weights = ~w)
    expected = coef(oracle)
    names(expected) = sub("^fit_", "", names(expected))
    expect_equal(coef(iv), expected[names(coef(iv))], tolerance = 1e-8)
    expect_equal(predict(iv,d), as.numeric(predict(oracle,reference)), tolerance = 1e-8)
    if(strategy == "compress") expect_equal(sort(predict(iv)),sort(predict(iv,d)),tolerance=1e-8)
  }
})
