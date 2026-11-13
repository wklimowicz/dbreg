if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("bit64", quietly = TRUE)) {
  set.seed(1568)
  d = expand.grid(g = 1:16, x = -1:1, rep = 1:8)
  d$g = bit64::as.integer64("9007199254740992") + d$g
  d$id = seq_len(nrow(d))
  d$absorbed = rep(seq_len(16), length.out = nrow(d))
  d$y = rpois(nrow(d), exp(0.5 + d$x / 5 + d$absorbed / 30))
  d$w = rep(c(0.5, 1, 2, 4), length.out = nrow(d))
  d$o = d$rep / 40
  option = options(dbreg.prediction_max_entries = 9)
  models = list(dbreg(y ~ x + absorbed | g, data = d, weights = ~w, offset = ~o,
    strategy = "demean", retain = "prediction"),
    dbglm(y ~ x + absorbed | g, data = d, weights = ~w, offset = ~o,
      family = poisson(), retain = "prediction"))
  options(option)
  pred = d[c(19L, 2L, 108L, 99L), c("id", "x", "absorbed", "g", "o")]
  pred$g[2L] = bit64::as.integer64("9007199254741992")
  pred$absorbed[3L] = pred$absorbed[3L] + 1
  for (model in models) {
    expect_identical(model$prediction_diagnostics$method, "database_prediction")
    expect_identical(model$prediction_diagnostics$portable_reason, "system_limit")
    ref = if (inherits(model, "dbglm")) {
      glm(y ~ x + factor(absorbed), data = d, weights = w, offset = o, family = poisson())
    } else lm(y ~ x + factor(absorbed), data = d, weights = w, offset = o)
    actual = predict(model, pred)
    expected = if (inherits(model, "dbglm")) predict(ref, type = "link") else fitted(ref)
    expected = expected[c(19L, 2L, 108L, 99L)]
    expected[2:3] = NA_real_
    expect_equal(actual, unname(expected), tolerance = 1e-8)
    if (requireNamespace("fixest", quietly = TRUE)) {
      values = fixest::fixef(model)
      expect_identical(names(values$g), unique(as.character(d$g)))
    }
    keyed = predict(model, pred, row_id = "id")
    conn = model$resource_state$conn
    other = DBI::dbConnect(duckdb::duckdb())
    expect_error(predict(model, pred, conn = other, row_id = "id"), "fitted owning session")
    DBI::dbDisconnect(other)
    close(model)
    expect_true(DBI::dbIsValid(conn))
    expect_error(predict(model, pred), "no live database session")
    expect_equal(as.data.frame(keyed)[[keyed$prediction_column]], actual[order(pred$id)])
    close(keyed)
    expect_false(DBI::dbIsValid(conn))
  }
}

if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("bit64", quietly = TRUE)) {
  conn = DBI::dbConnect(duckdb::duckdb(bigint = "integer64"))
  resources = dbreg:::fit_resources()
  resources$conn = conn
  d = data.frame(id = bit64::as.integer64("9007199254740992") + seq_len(530L),
    text = rep(c("a'b", "path\\file", NA), length.out = 530L),
    date = as.Date("2025-01-01") + seq_len(530L),
    flag = rep(c(TRUE, FALSE), 265L))
  table = dbreg:::sql_native_prediction_input(conn, d, names(d), resources)
  actual = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "ORDER BY id"))
  expect_identical(as.character(actual$id), as.character(d$id))
  expect_identical(actual$text, d$text)
  expect_identical(actual$date, d$date)
  expect_identical(actual$flag, d$flag)
  dbreg:::release_fit_resources(resources)
  expect_identical(DBI::dbListTables(conn), character())
  DBI::dbDisconnect(conn)
}

if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("fixest", quietly = TRUE)) local({
  set.seed(20261119)
  n = 80L
  d = data.frame(
    x = rnorm(n),
    fl = rep(c(FALSE, TRUE), length.out = n),
    fn = rep(c(1.0, 2.0, 3.5, 4.0), length.out = n)
  )
  d[["y"]] = 1 + 0.5 * d[["x"]] + 0.4 * d[["fl"]] + 0.2 * d[["fn"]] + rnorm(n, sd = 0.1)
  for (fml in list(y ~ x | fn, y ~ x | fl ^ fn)) {
    portable = dbreg(fml, data = d, strategy = "demean", retain = "all")
    old = options(dbreg.prediction_max_entries = 4L)
    db_fit = tryCatch(dbreg(fml, data = d, strategy = "demean", retain = "all"), finally = options(old))
    on.exit({ close(portable); close(db_fit) }, add = TRUE)
    expect_identical(names(fixest::fixef(db_fit)[[1L]]), names(fixest::fixef(portable)[[1L]]))
    expect_equal(unname(fixest::fixef(db_fit)[[1L]]), unname(fixest::fixef(portable)[[1L]]), tolerance = 1e-10)
    close(portable)
    close(db_fit)
  }
})
