local({
  fixtures = stats::setNames(c("", "Zg==", "Zm8=", "Zm9v", "Zm9vYg==", "Zm9vYmE=", "Zm9vYmFy"),
    c("", "f", "fo", "foo", "foob", "fooba", "foobar"))
  for (i in seq_along(fixtures)) {
    expect_identical(dbreg:::raw_to_base64(charToRaw(names(fixtures)[[i]])), unname(fixtures[[i]]))
  }
  expect_identical(dbreg:::read_numeric_json("[1e-200,3.1415926535897931,null]"), c(1e-200, pi, NA_real_))
  expect_identical(dbreg:::read_numeric_json('{"large":1e200,"small":1e-200}'),
    list(large = 1e200, small = 1e-200))
  for (value in c("[1,,2]", "[1,]", "{\"x\":1,}", "[01]", "[1+2]", "[[1]]", '["1"]',
    '{"x":1,"x":2}', '{"x":TRUE}')) {
    expect_error(dbreg:::read_numeric_json(value), "Invalid numeric SQL payload")
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    set.seed(1010)
    for (n in c(0:12, 255, 65536)) {
      bytes = as.raw(sample.int(256, n, TRUE) - 1L)
      expect_identical(dbreg:::raw_to_base64(bytes),
        gsub("[\r\n]", "", jsonlite::base64_enc(bytes)))
    }
    strings = c("café", "雪", "quote\"slash\\", intToUtf8(1:31), NA_character_)
    payload = paste0("[", paste(dbreg:::json_string(strings), collapse = ","), "]")
    expect_identical(jsonlite::fromJSON(payload), strings)
    values = c(0, pi, 1e-200, 1e200, -1e200, .Machine[["double.xmin"]])
    payload = dbreg:::json_lookup_rows(data.frame(x = values))
    expect_identical(as.numeric(jsonlite::fromJSON(payload)), values)
  }
})
