# Bound actual DBI result sizes, rather than testing the spelling of SQL.
local({
  state = new.env(parent = emptyenv())
  state[["rows"]] = numeric()
  state[["queries"]] = 0L
  state[["observe"]] = function(value) {
    state[["queries"]] = state[["queries"]]+1L
    if (is.data.frame(value)) state[["rows"]] = c(state[["rows"]], nrow(value))
  }
  key = ".dbreg_review_transfer__"
  stopifnot(!exists(key, envir = .GlobalEnv, inherits = FALSE))
  assign(key, state, envir = .GlobalEnv)
  invisible(capture.output(trace("dbGetQuery", where = asNamespace("DBI"),
    exit = quote(get(".dbreg_review_transfer__", envir = .GlobalEnv)[["observe"]](returnValue())), print = FALSE)))
  on.exit({
    invisible(capture.output(untrace("dbGetQuery", where = asNamespace("DBI"))))
    rm(list = key, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(451)
  n = 600L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), y = rnorm(n))
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  clustered = dbreg(y ~ x+z, data = d, conn = conn, vcov = ~id, strategy = "moments")
  expect_equal(nobs(clustered), n)
  expect_true(max(state[["rows"]]) < n, info = paste("cluster score maximum returned rows:", max(state[["rows"]])))
  state[["rows"]] = numeric()
  bins = dbbinsreg(y ~ x+z, data = d, conn = conn, nbins = 5, points = c(1,0), plot = FALSE)
  expect_equal(bins[["opt"]][["N"]], n)
  expect_true(max(state[["rows"]]) < n, info = paste("binscatter maximum returned rows:", max(state[["rows"]]), "strategy:", bins[["model"]][["strategy"]]))

  state[["queries"]] = 0L
  automatic = dbreg(y ~ x+z, data = d, conn = conn)
  auto_queries = state[["queries"]]
  state[["queries"]] = 0L
  forced = dbreg(y ~ x+z, data = d, conn = conn, strategy = "moments")
  expect_equal(coef(automatic), coef(forced))
  expect_equal(auto_queries, state[["queries"]])

  # Decimal keys require exact text, but collecting N results should fetch N once.
  DBI::dbWriteTable(conn, "review_prediction", d)
  DBI::dbExecute(conn, "CREATE TEMP TABLE review_decimal AS SELECT CAST(id AS DECIMAL(30,0)) AS id,x,z,y FROM review_prediction")
  fit = dbreg(y ~ x+z, conn = conn, table = "review_decimal", row_id = "id", retain = "prediction", strategy = "moments")
  on.exit(close(fit), add = TRUE)
  predictions = predict(fit, conn = conn, table = "review_decimal", row_id = "id")
  on.exit(close(predictions), add = TRUE)
  state[["rows"]] = numeric()
  values = as.data.frame(predictions)
  expect_equal(nrow(values), n)
  expect_identical(values[["id"]], as.character(seq_len(n)))
  expect_equal(sum(state[["rows"]] >= n), 1L)
})
