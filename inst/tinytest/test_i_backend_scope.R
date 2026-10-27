library(DBI)
source(system.file("tinytest/helper_i_backend_scope.R", package = "dbreg"), local = TRUE)
local({
  actual = dbreg_test_existing_i_backends()
  data = actual[["data"]]
  for (name in names(actual[["captures"]])) {
    capture = actual[["captures"]][[name]]
    index = as.integer(sub("^.*-", "", name))
    if (startsWith(name, "BigQuery")) {
      expect_identical(capture[["backend"]], "bigquery")
      queries = capture[["queries"]]
      expect_true(grepl(" LIMIT 0$", queries[[1L]]))
      expect_true(any(startsWith(queries, "SELECT typeof(`")))
      expect_false(any(grepl("SQL_VARIANT_PROPERTY|native_decimal_text| TOP ", queries)))
      if (index >= 3L) {
        expected_type = if (index == 3L) "int64" else "float64"
        expect_true(grepl(paste0("key SQL type '", expected_type, "'"), capture[["error"]], fixed = TRUE))
        expect_null(capture[["compiled"]])
        expect_false(any(grepl(" GROUP BY ", queries, fixed = TRUE)))
      } else {
        expect_null(capture[["error"]])
        spec = capture[["compiled"]][["i_specs"]][[1L]]
        expect_identical(spec[["type"]], "character")
        expect_identical(spec[["values"]], sort(unique(data[["text"]])))
        keys = sort(unique(data[["text"]]))
        if (index == 1L) keys = keys[-1L]
        expected = outer(data[["text"]], keys, `==`) * 1
        expect_equal(unname(as.matrix(capture[["values"]])), unname(expected))
        expect_equal(unname(dbreg:::evaluate_i_spec(spec, data[["text"]])), unname(expected))
      }
      next
    }
    expect_identical(capture[["backend"]], if (startsWith(name, "Unclassified")) "other" else "athena")
    queries = capture[["queries"]]
    # The existing zero-row and typeof queries survive unchanged. This fixture
    # rejects dbSendQuery(), so a driver-metadata adapter cannot pass unnoticed.
    expect_true(grepl(" LIMIT 0$", queries[[1L]]))
    expect_true(any(startsWith(queries, "SELECT typeof(")))
    expect_false(any(grepl("SQL_VARIANT_PROPERTY|native_decimal_text| TOP ", queries)))
    compiled = capture[["compiled"]]
    spec = compiled[["i_specs"]][[1L]]
    all = index == 6L
    source = switch(index, data[["text"]], data[["raw_key"]], data[["adjacent"]], data[["text"]], data[["day"]], data[["text"]])
    keys = switch(index, sort(unique(data[["text"]])), unique(data[["raw_key"]]),
      unique(data[["adjacent"]]), sort(unique(data[["text"]])), sort(unique(data[["day"]])), sort(unique(data[["text"]])))
    if (!all && index != 4L) keys = keys[-1L]
    expected = outer(source, keys, `==`) * 1
    if (index == 4L) expected = expected * data[["flag"]]
    expect_equal(unname(as.matrix(capture[["values"]])), unname(expected))
    expect_equal(unname(dbreg:::evaluate_i_spec(spec, source, if (index == 4L) data[["flag"]])), unname(expected))
    if (index %in% c(1L, 4L, 6L)) expect_identical(spec[["values"]], sort(unique(data[["text"]])))
    if (index == 2L) expect_identical(spec[["values"]], unique(data[["raw_key"]]))
  }
  for (removed in actual[["removed"]]) {
    expect_true(grepl("not supported by dbreg", removed[["error"]], fixed = TRUE))
    expect_identical(removed[["queries"]], character())
  }
})
