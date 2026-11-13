# A retained child on a package-owned connection participates in the same
# session lifetime, regardless of which object is closed first.
local({
  d = data.frame(id = 1:40, x = sin(1:40), y = 1:40/10 + cos(1:40))
  base = dbreg(y ~ x, data = d, strategy = "moments")
  for (parent_first in c(FALSE, TRUE)) {
    p = predict(base, d, row_id = "id")
    conn = p[["conn"]]
    DBI::dbWriteTable(conn, "child_source", d)
    expect_error(dbreg(y ~ x, conn = conn, table = "child_source", row_id = "id",
      sample_table = "forbidden"), "caller-owned")
    expect_false(DBI::dbExistsTable(conn, "forbidden"))
    child = dbreg(y ~ x, conn = conn, table = "child_source", row_id = "id",
      retain = "inference", strategy = "moments")
    expected = as.numeric(fitted(lm(y ~ x, d)))
    if (parent_first) {
      close(p)
      expect_true(DBI::dbIsValid(conn))
      rows = fitted(child)
      expect_equal(as.data.frame(rows)[[".fitted"]], expected, tolerance = 1e-10)
      close(rows)
      expect_equal(as.numeric(vcov(child, vcov = "hc1")),
        as.numeric(vcov(dbreg(y ~ x, data = d, strategy = "moments", vcov = "hc1"))), tolerance = 1e-10)
      close(child)
    } else {
      close(child)
      expect_true(DBI::dbIsValid(conn))
      expect_equal(as.data.frame(p)[[p[["prediction_column"]]]], expected, tolerance = 1e-10)
      close(p)
    }
    expect_false(DBI::dbIsValid(conn))
  }
})
