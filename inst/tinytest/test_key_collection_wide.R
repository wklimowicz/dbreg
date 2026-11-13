# Native exact SQL keys must survive drivers that collect them as doubles.
# This direct collector test does not depend on retained or prediction APIs.
local({
  collect_exact = function(conn, sql, keys) {
    dbreg:::read_keyed_query(conn, sql, keys, exact_keys = TRUE)
  }
  run = function(mode) {
    conn = DBI::dbConnect(duckdb::duckdb(bigint = mode))
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
    query = function(type, keys) {
      values = paste0("(", seq_along(keys), ",", ifelse(is.na(keys), "NULL",
        as.character(DBI::dbQuoteString(conn, keys))), ")")
      paste0("SELECT id,CAST(k AS ", type, ") AS k,CAST(1.25 AS DOUBLE) AS statistic ",
        "FROM (VALUES ", paste(values, collapse = ","), ") source(id,k) ORDER BY id")
    }
    for (type in c("BIGINT", "DECIMAL(20,0)", "HUGEINT")) {
      keys = c("9007199254740992", "9007199254740993", NA, "-9007199254740993")
      values = collect_exact(conn, query(type, keys), "k")
      expect_identical(as.character(values$k), keys, info = paste(mode, type))
      expect_true(inherits(values$k, "integer64"))
      expect_identical(values$statistic, rep(1.25, length(keys)))
    }
    for (type in c("DECIMAL(38,0)", "HUGEINT", "UHUGEINT", "UBIGINT")) {
      keys = c("9223372036854775808", "9223372036854775809", "18446744073709551615")
      values = collect_exact(conn, query(type, keys), "k")
      expect_identical(values$k, keys, info = paste(mode, type))
      expect_identical(length(unique(values$k)), length(keys))
    }
    # bit64 reserves the minimum signed bit pattern for NA. Exact SQL keys
    # containing that value therefore remain text even on integer64 drivers.
    keys = c("-9223372036854775808", "-9223372036854775807", "9223372036854775807", NA)
    values = collect_exact(conn, query("BIGINT", keys), "k")
    expect_identical(values$k, keys)
    keys = c("0.10000000000000000001", "0.10000000000000000002", "1.00000000000000000000")
    values = collect_exact(conn, query("DECIMAL(38,20)", keys), "k")
    expect_identical(values$k, keys)
    expect_identical(length(unique(values$k)), length(keys))
    # A large native DOUBLE is an intentional floating-point key, not an exact
    # decimal declaration. Preserve its numeric contract.
    values = collect_exact(conn,
      "SELECT * FROM (VALUES (1.0::DOUBLE),(10000000000000000.0::DOUBLE)) source(k)", "k")
    expect_identical(values$k, c(1, 1e16))

    # Exercise Spark's typeof/STRING branch with actual local SQL. This is a
    # controlled adapter check, not a native Databricks certificate.
    collect = dbreg:::read_keyed_query
    env = new.env(parent = environment(collect))
    env$detect_backend = function(...) list(name = "spark", is_databricks = TRUE)
    environment(collect) = env
    keys = c("9223372036854775808", "9223372036854775809")
    values = collect(conn, query("DECIMAL(38,0)", keys), "k", exact_keys = TRUE)
    expect_identical(values$k, keys)
    keys = c("0.10000000000000000001", "0.10000000000000000002")
    values = collect(conn, query("DECIMAL(38,20)", keys), "k", exact_keys = TRUE)
    expect_identical(values$k, keys)
    # Numeric GROUP BY columns are also read through this helper during
    # compressed estimation. Their arithmetic type must remain numeric.
    basis = dbreg:::read_keyed_query(conn,
      "SELECT * FROM (VALUES (0.0::DECIMAL(2,1)),(1.0::DECIMAL(2,1))) source(k)", "k")
    expect_identical(basis$k, c(0, 1))
    decimal_sql = query("DECIMAL(38,20)", keys)
    grouped = dbreg:::read_keyed_query(conn, decimal_sql, "k")
    expect_identical(grouped$k, DBI::dbGetQuery(conn, decimal_sql)$k)
    bigint_sql = query("BIGINT", c("9007199254740992", "9007199254740993"))
    grouped = dbreg:::read_keyed_query(conn, bigint_sql, "k")
    expect_identical(as.character(grouped$k), c("9007199254740992", "9007199254740993"))
    expect_true(inherits(grouped$k, "integer64"))
    expect_identical(DBI::dbListTables(conn), character())
  }
  if (requireNamespace("bit64", quietly = TRUE)) {
    for (mode in c("numeric", "integer64")) run(mode)
  }
})
