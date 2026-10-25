if (identical(tolower(Sys.getenv("DBREG_TEST_SPARK")), "true")) {
  source(system.file("tinytest/helper_spark.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbreg_test_spark_connection()
    on.exit(sparklyr::spark_disconnect(conn))
    dbreg_test_spark_algorithms(conn)
  })
}
