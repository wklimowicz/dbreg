# Unsupported connections must fail before falling through to generic SQL.
for (name in c("SQLite", "PostgreSQL")) {
  expect_error(dbreg:::backend_info_from_metadata(list(dbms.name = name)),
    "not supported by dbreg")
}
for (class in c("SQLiteConnection", "PostgreSQLConnection", "PqConnection")) {
  expect_error(dbreg:::backend_info_from_metadata(connection_class = class),
    "not supported by dbreg")
}
expect_identical(dbreg:::backend_info_from_metadata(list(dbms.name = "DuckDB"))[["name"]], "duckdb")
expect_identical(dbreg:::backend_info_from_metadata(list(dbms.name = "MySQL"))[["name"]], "mysql")
expect_identical(dbreg:::backend_info_from_metadata(list(dbms.name = "Microsoft SQL Server"))[["name"]], "sqlserver")
