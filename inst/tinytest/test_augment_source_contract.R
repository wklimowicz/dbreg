local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(9218)
  d = data.frame(id = sample(1:96), g = rep(1:6, each = 16),
    x = rep(c(-1, 0, 1), 32), w = rep(c(.5, 1, 2), 32),
    o = rep(c(-.1, .1), 48), use = rep(c(TRUE, TRUE, TRUE, FALSE), 24))
  d$y = .3 + .4 * d$x + d$g / 5 + d$o + rnorm(nrow(d), sd = .3)
  d$yp = rpois(nrow(d), exp(.2 + .2 * d$x + d$g / 20 + d$o))
  d$category = factor(rep(letters[1:3], 32))
  d$date = as.Date("2026-01-01") + seq_len(nrow(d))
  d$x[3L] = NA_real_
  d$w[5L] = 0
  keep = which(d$use & !is.na(d$x) & d$w > 0)
  forms = list(
    list(FALSE, y ~ x, y ~ x, "moments"),
    list(FALSE, y ~ x | g, y ~ x + factor(g), "demean"),
    list(FALSE, y ~ x | g, y ~ x + factor(g), "compress")
  )
  if (exists("dbglm", mode = "function")) {
    forms = c(forms, list(
      list(TRUE, yp ~ x, yp ~ x, "compress"),
      list(TRUE, yp ~ x | g, yp ~ x + factor(g), "compress")
    ))
  }
  for (source in c("R", "database")) for (explicit in c(FALSE, TRUE)) {
    if (source == "database" && !explicit) next
    for (case in forms) {
      initial = sort(DBI::dbListTables(conn))
      if (source == "database") DBI::dbWriteTable(conn, "original", d)
      maker = if (case[[1L]]) dbglm else dbreg
      args = list(fml = case[[2L]], conn = conn, weights = ~w, offset = ~o,
        subset = ~use, retain = "all", strategy = case[[4L]])
      if (case[[1L]]) args$family = poisson()
      if (explicit) args$row_id = "id"
      if (source == "R") args$data = d else args$table = "original"
      fit = do.call(maker, args)
      dense = if (case[[1L]]) glm(case[[3L]], data = d[keep, ], weights = w, offset = o,
        family = poisson(), control = glm.control(epsilon = 1e-14, maxit = 1000)) else
          lm(case[[3L]], data = d[keep, ], weights = w, offset = o)
      order = if (explicit) order(d$id[keep]) else seq_along(keep)
      expected = d[keep[order], , drop = FALSE]
      mu = unname(if (case[[1L]]) predict(dense, type = "response") else fitted(dense))[order]
      if (source == "database") DBI::dbExecute(conn, "UPDATE original SET y=999,yp=999,x=99")
      actual = generics::augment(fit)
      expect_identical(names(actual), c(names(d), ".fitted", ".resid"))
      expect_equal(actual$id, expected$id)
      expect_equal(actual$x, expected$x)
      expect_equal(as.character(actual$category), as.character(expected$category))
      expect_equal(as.character(actual$date), as.character(expected$date))
      expect_equal(actual$.fitted, mu, tolerance = 1e-7)
      response = if (case[[1L]]) expected$yp else expected$y
      expect_equal(actual$.resid, response - mu, tolerance = 1e-7)
      expect_equal(nobs(fit), length(keep))
      if (!explicit && source == "R") {
        supplied = generics::augment(fit, data = d)
        expect_identical(supplied[names(d)], expected)
      } else expect_error(generics::augment(fit, data = d), "original R input")
      result = predict(fit, d[c(1, 11, 18), ], conn = conn, row_id = "id")
      saved = as.data.frame(result)
      close(fit)
      expect_identical(as.data.frame(result), saved)
      expect_error(generics::augment(fit), "no live database session")
      close(result)
      if (source == "database") DBI::dbRemoveTable(conn, "original")
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
  }
  # These are caller columns, not grounds for a name-prefix exclusion rule.
  d$v_dbreg_input_position = seq_len(nrow(d)) / 100
  d$v_dbreg_input_position_1 = letters[rep(1:3, 32)]
  d$dbreg_input_position = rep(c(FALSE, TRUE), 48)
  glm_flags = if (exists("dbglm", mode = "function")) c(FALSE, TRUE) else FALSE
  for (glm in glm_flags) {
    f = if (glm) dbglm(yp ~ x, data = d, conn = conn, family = poisson(),
      retain = "inference") else dbreg(y ~ x, data = d, conn = conn, retain = "inference")
    actual = generics::augment(f)
    expect_identical(names(actual), c(names(d), ".fitted", ".resid"))
    expect_equal(actual$v_dbreg_input_position, d$v_dbreg_input_position[!is.na(d$x)])
    expect_equal(actual$v_dbreg_input_position_1, d$v_dbreg_input_position_1[!is.na(d$x)])
    expect_equal(actual$dbreg_input_position, d$dbreg_input_position[!is.na(d$x)])
    close(f)
  }
})

local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(x = rep(-1:1, 20), y = rep(c(0, 1, 1), 20))
  makers = if (exists("dbglm", mode = "function")) list(dbreg, dbglm) else list(dbreg)
  for (maker in makers) {
    for (args in list(list(fml = y ~ v_dbreg_input_position),
        list(fml = y ~ x, weights = "v_dbreg_input_position"),
        list(fml = y ~ x, offset = "v_dbreg_input_position"),
        list(fml = y ~ x, vcov = ~v_dbreg_input_position),
        list(fml = y ~ x, subset = ~v_dbreg_input_position > 2))) {
      result = tryCatch(do.call(maker, c(args, list(data = d, conn = conn, retain = "all"))), error = identity)
      expect_true(inherits(result, "error"))
      if (inherits(result, "error")) expect_true(grepl("missing|Missing|unknown|Unknown|absent|not found", conditionMessage(result)))
      else close(result)
      expect_identical(DBI::dbListTables(conn), character())
    }
    d$.fitted = seq_len(nrow(d))
    f = maker(y ~ x, data = d, conn = conn, retain = "inference")
    expect_error(generics::augment(f), "conflict")
    expect_error(generics::augment(f, data = d), "conflict")
    safe = d
    names(safe)[names(safe) == ".fitted"] = "original_fitted"
    actual = generics::augment(f, data = safe)
    expect_equal(actual$original_fitted, safe$original_fitted)
    expect_true(all(c(".fitted", ".resid") %in% names(actual)))
    close(f)
    d$.fitted = NULL
  }
})

local({
  # Exercise native collection's exact-integer branch against real BIGINT
  # storage. This controlled branch test is not a native Databricks certificate.
  conn = DBI::dbConnect(duckdb::duckdb(bigint = "numeric"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbExecute(conn, paste("CREATE TEMP TABLE original AS SELECT i AS id,",
    "'9007199254740992'::BIGINT+i%3 AS big, i::DOUBLE AS x,",
    "(i%5)::DOUBLE AS y FROM range(30) t(i)"))
  fit = dbreg(y ~ x, conn = conn, table = "original", row_id = "id", retain = "inference")
  on.exit(close(fit), add = TRUE, after = FALSE)
  collect = getS3method("augment", "dbreg", envir = asNamespace("dbreg"))
  e = new.env(parent = environment(collect))
  e$detect_backend = function(...) list(is_databricks = TRUE)
  environment(collect) = e
  actual = collect(fit)
  keys = c("9007199254740992", "9007199254740993", "9007199254740994")
  expect_identical(actual$id, as.character(0:29))
  expect_identical(actual$big, keys[(0:29) %% 3 + 1L])
  expect_identical(names(actual), c("id", "big", "x", "y", ".fitted", ".resid"))
})
