if (requireNamespace("duckdb", quietly = TRUE)) local({
  set.seed(648)
  d = data.frame(id = 1:45, x = rnorm(45), y = rpois(45, 3), o = rnorm(45, 0, .1),
    g = rep(1:3, length.out = 45), w = 1)
  d$x[3L] = NA_real_
  d$w[9L] = 0
  d$o[10L] = NA_real_
  d$keep = d$id != 11L
  omitted = c(3L, 9L, 10L, 11L)
  lm_fit = dbreg(y ~ x | g, data = d, weights = ~w, offset = ~o, subset = ~keep,
    retain = "all", strategy = "demean")
  glm_fit = dbglm(y ~ x | g, data = d, weights = ~w, offset = ~o, subset = ~keep,
    retain = "all", family = poisson())
  for (fit in list(lm_fit, glm_fit)) {
    observed = predict(fit, sample = "original")
    expect_identical(which(is.na(observed)), omitted)
    expected = if (inherits(fit, "dbglm")) fit$family$linkfun(fitted(fit)) else fitted(fit)
    expect_equal(observed[-omitted], expected)
    expect_equal(predict(fit), expected)
    expect_error(predict(fit, d, sample = "original"), "retained fitting source")
    close(fit)
  }
  expect_error(predict(dbreg(y ~ x, data = d), sample = "original"), "retain")
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "source", d)
  model = dbglm(y ~ x | g, conn = conn, table = "source", row_id = "id", weights = ~w, offset = ~o,
    subset = ~keep, retain = "inference", family = poisson())
  output = predict(model, type = "link", sample = "original")
  exported = predict(model, type = "link", prediction_table = "em eta")
  original = predict(model, type = "response", sample = "original", prediction_table = "original_mu")
  expect_inherits(exported, "dbreg_prediction")
  expect_error(predict(model, prediction_table = "em eta"), "already exists")
  expect_error(predict(model, row_id = "g", prediction_table = "other"), "explicit row_id")
  other = DBI::dbConnect(duckdb::duckdb())
  expect_error(predict(model, conn = other, prediction_table = "other"), "fitted owning session")
  DBI::dbDisconnect(other, shutdown = TRUE)
  close(model)
  expect_error(predict(model, prediction_table = "other"), "no live database session")
  observed = as.data.frame(output)
  expect_identical(observed$id, d$id)
  expect_identical(which(is.na(observed$.link)), omitted)
  expect_identical(head(output, 2L), utils::head(observed, 2L))
  keyed = DBI::dbGetQuery(conn,
    'SELECT s.id, e.prediction FROM source s JOIN "em eta" e USING (id) ORDER BY s.id')
  expect_identical(keyed$id, d$id[-omitted])
  expect_equal(keyed$prediction, observed$.link[-omitted])
  expect_identical(as.data.frame(exported), keyed)
  response = as.data.frame(original)
  expect_identical(response$id, d$id)
  expect_equal(response$prediction, exp(observed$.link))
  close(output)
  close(exported)
  close(original)
  expect_identical(sort(DBI::dbListTables(conn)), c("em eta", "original_mu", "source"))
  DBI::dbRemoveTable(conn, "em eta")
  DBI::dbRemoveTable(conn, "original_mu")
  expect_identical(DBI::dbListTables(conn), "source")
})
