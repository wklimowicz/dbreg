if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(817)
  d = data.frame(id = 1:60, sub = rep(c("a", "b"), 30), x = rep(-1:1, 20), f = rep(letters[1:3], each = 20))
  d$o = d$id / 20
  d$y = rpois(nrow(d), exp(d$x / 4 + d$o))
  conn = DBI::dbConnect(duckdb::duckdb(bigint = "numeric"))
  DBI::dbWriteTable(conn, "points", d)
  DBI::dbExecute(conn, paste("CREATE TEMP TABLE keyed_points AS SELECT *,",
    "CAST('9007199254740992' AS BIGINT) + CAST(id AS BIGINT) AS \"row key\" FROM points"))
  before = sort(DBI::dbListTables(conn))
  newdata = d[c(4L, 2L, 9L, 6L), ]
  newdata$x[2L] = NA_real_
  for (model in list(dbreg(y ~ x * factor(f), data = d, offset = ~o),
                     dbglm(y ~ x * factor(f), data = d, offset = ~o, family = poisson()))) {
    res = predict(model, newdata, row_id = c("id", "sub"))
    expect_true(inherits(res, "dbreg_prediction"))
    out = as.data.frame(res)
    order = order(newdata$id)
    expect_identical(out$id, newdata$id[order])
    expect_equal(out[[res$prediction_column]], unname(predict(model, newdata))[order], tolerance = 1e-10)
    expect_identical(head(res, 2L), utils::head(out, 2L))
    second = predict(model, conn = res$conn, newdata = newdata, row_id = "id")
    owner = res$resource_state$conn
    close(res)
    expect_true(DBI::dbIsValid(owner))
    expect_equal(as.data.frame(second)[[second$prediction_column]], out[[res$prediction_column]])
    close(second)
    expect_false(DBI::dbIsValid(owner))
    expect_error(as.data.frame(res), "no live database session")
    if (requireNamespace("bit64", quietly = TRUE)) {
      pred = predict(model, conn = conn, table = "keyed_points", row_id = c("row key", "sub"))
      got = as.data.frame(pred)
      expect_equal(got[[pred$prediction_column]], unname(predict(model, d)), tolerance = 1e-10)
      expect_identical(as.character(got[["row key"]]),
          as.character(bit64::as.integer64("9007199254740992") + d$id))
      close(pred)
    }
    expect_identical(sort(DBI::dbListTables(conn)), before)
  }
  model = dbreg(y ~ x * factor(f), data = d, offset = ~o)
  unknown = newdata
  unknown$f[1L] = "unknown"
  expect_error(predict(model, unknown, row_id = "id"), "unseen factor levels")
  duplicate = rbind(newdata, newdata[1L, ])
  expect_error(predict(model, duplicate, row_id = "id"), "unique")
  expect_error(predict(model, conn = conn, table = "points", row_id = "id", interval = "confidence"), "interval")
  exported = predict(model, conn = conn, table = "points", row_id = "id", prediction_table = "exported")
  expect_equal(as.data.frame(exported)[[exported$prediction_column]], unname(predict(model, d)))
  close(exported)
  expect_true(DBI::dbExistsTable(conn, "exported"))
  expect_error(predict(model, conn = conn, table = "points", row_id = "id", prediction_table = "exported"), "already exists")
  if (requireNamespace("bit64", quietly = TRUE)) {
    exact = bit64::as.integer64("9007199254740992") + rep(0:2, each = 20)
    train = transform(d, phase = exact)
    typed = dbreg(y ~ i(phase) + x, data = train)
    points = data.frame(id = 1:4, x = 1,
      phase = c("9007199254740992", "9007199254740993", "9007199254740994", "9007199254740995"))
    keyed = predict(typed, points, row_id = "id")
    expect_equal(as.data.frame(keyed)[[keyed$prediction_column]], unname(predict(typed, points)), tolerance = 1e-10)
    close(keyed)
    points$phase = 9007199254740994
    expect_error(predict(typed, points, row_id = "id"), "integer keys require exact")
  }
  DBI::dbDisconnect(conn)
}
