local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  saved = options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(saved), add = TRUE)
  set.seed(11128)
  n = 600L
  d = data.frame(
    id = seq_len(n),
    x = sample(-1:1, n, TRUE),
    z = sample(-2:2, n, TRUE),
    q = sample(-1:2, n, TRUE),
    g = sample(1:8, n, TRUE),
    h = sample(1:5, n, TRUE),
    w = runif(n, .2, 2)
  )
  d[["y"]] = .7 *
    d[["x"]] +
    sin(d[["g"]]) +
    cos(d[["h"]]) +
    (.2 + d[["g"]] / 20) * d[["z"]] +
    d[["g"]] / 15 * d[["q"]] +
    rnorm(n)
  DBI::dbWriteTable(conn, "compressed_vary", d)
  initial = sort(DBI::dbListTables(conn))
  # Independent weighted raw dummy regression; inference must match its actual
  # structural coordinates, not any displayed coefficient prefix.
  for (policy in c("full", "none", "nested")) {
    for (kind in list("iid", "hc1", ~g)) {
      fml = y ~ x | g[z, q] + h
      fit = dbreg(
        fml,
        conn = conn,
        table = "compressed_vary",
        weights = ~w,
        vcov = kind,
        ssc = policy,
        strategy = "compress"
      )
      ref = fixest::feols(
        fml,
        d,
        weights = ~w,
        vcov = if (identical(kind, "hc1")) "hetero" else kind,
        ssc = fixest::ssc(
          K.fixef = if (policy == "nested") "nonnested" else policy
        ),
        fixef.rm = "none",
        fixef.tol = 1e-11,
        notes = FALSE
      )
      expect_equal(
        as.numeric(coef(fit)),
        as.numeric(coef(ref)),
        tolerance = 1e-8
      )
      expect_equal(
        as.numeric(vcov(fit)),
        as.numeric(vcov(ref)),
        tolerance = 1e-8
      )
      D = model.matrix(
        ~ x + factor(g) + factor(h) + factor(g):z + factor(g):q,
        d
      )
      beta = qr.coef(qr(D * sqrt(d[["w"]])), d[["y"]] * sqrt(d[["w"]]))
      pred = unname(drop(D %*% beta))
      expect_equal(predict(fit, d), pred, tolerance = 1e-8)
      expect_equal(predict(fit), predict(fit, fit[["data"]]), tolerance = 1e-8)
      cells = fit[["data"]][rev(seq_len(nrow(fit[["data"]]))), , drop = FALSE]
      expect_equal(predict(fit, cells), rev(predict(fit)), tolerance = 1e-8)
      expect_equal(
        predict(unserialize(serialize(fit, NULL)), d),
        pred,
        tolerance = 1e-8
      )
      bad = d[1:4, ]
      bad[["g"]][1] = 999
      expect_true(is.na(predict(fit, bad)[1L]))
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
  }
  # Expanded structural producers, duplicate public labels and absorbed trends.
  d[["group"]] = rep(c("a", "b"), n / 2)
  d[["groupb"]] = rnorm(n)
  for (fml in list(
    y ~ x + z | g[z, q] + h,
    y ~ x + group + groupb | g[z] + h,
    y ~ x | g[z] + g[q] + h
  )) {
    fit = dbreg(fml, data = d, weights = ~w, strategy = "compress")
    if ("group" %in% all.vars(fml)) {
      D = model.matrix(
        ~ x + group + groupb + factor(g) + factor(h) + factor(g):z,
        d,
        contrasts.arg = list(group = contr.treatment(c("a", "b")))
      )
      sw = sqrt(d[["w"]])
      beta = qr.coef(qr(D * sw), d[["y"]] * sw)
      residual = d[["y"]] - drop(D %*% beta)
      V = solve(crossprod(D, d[["w"]] * D)) *
        sum(d[["w"]] * residual^2) /
        (n - ncol(D))
      expect_equal(
        as.numeric(coef(fit)),
        as.numeric(beta[2:4]),
        tolerance = 1e-8
      )
      expect_equal(
        as.numeric(vcov(fit)),
        as.numeric(V[2:4, 2:4]),
        tolerance = 1e-8
      )
      expect_equal(predict(fit, d), as.numeric(D %*% beta), tolerance = 1e-8)
    } else {
      ref = fixest::feols(
        fml,
        d,
        weights = ~w,
        fixef.rm = "none",
        ssc = fixest::ssc(K.fixef = "full"),
        fixef.tol = 1e-11,
        notes = FALSE
      )
      expect_equal(
        as.numeric(coef(fit)),
        as.numeric(coef(ref)),
        tolerance = 1e-8
      )
      expect_equal(
        as.numeric(vcov(fit)),
        as.numeric(vcov(ref)),
        tolerance = 1e-8
      )
    }
    expect_true(all(is.finite(predict(fit, d))))
  }
  # A single observed FE level is a constant producer, not an R contrasts error.
  constant = d
  constant[["g"]] = 7L
  fit = dbreg(
    y ~ x | g[z],
    data = constant,
    weights = ~w,
    strategy = "compress"
  )
  ref = fixest::feols(
    y ~ x | g[z],
    constant,
    weights = ~w,
    fixef.rm = "none",
    ssc = fixest::ssc(K.fixef = "full"),
    notes = FALSE
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(
    predict(fit, constant),
    as.numeric(predict(ref)),
    tolerance = 1e-8
  )
  bad = constant[1:2, ]
  bad[["g"]][1] = 8L
  expect_true(is.na(predict(fit, bad)[1L]))
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})


local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11129)
  n = 500L
  d = data.frame(
    x = sample(-1:1, n, TRUE),
    g = sample(1:9, n, TRUE),
    h = sample(1:5, n, TRUE),
    z = sample(-2:2, n, TRUE),
    q = sample(-2:2, n, TRUE),
    w = runif(n, .2, 2),
    o = rnorm(n)
  )
  d[["y"]] = d[["o"]] +
    .8 * d[["x"]] +
    sin(d[["g"]]) +
    cos(d[["h"]]) +
    d[["g"]] / 20 * d[["z"]] +
    rnorm(n)
  DBI::dbWriteTable(conn, "raw_cv", d)
  initial = sort(DBI::dbListTables(conn))
  shifted = d
  shifted[["z"]] = 1e12 + shifted[["z"]]
  fit = dbreg(
    y ~ x | g[z] + h,
    data = shifted,
    weights = ~w,
    offset = ~o,
    strategy = "compress"
  )
  adjusted = d
  adjusted[["y"]] = d[["y"]] - d[["o"]]
  ref = fixest::feols(
    y ~ x | g[z] + h,
    adjusted,
    weights = ~w,
    ssc = fixest::ssc(K.fixef = "full"),
    fixef.rm = "none",
    fixef.tol = 1e-11,
    notes = FALSE
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  expect_equal(
    predict(fit, shifted),
    as.numeric(predict(ref)) + d[["o"]],
    tolerance = 1e-8
  )
  # Moment-field spellings can be FE, raw numeric RHS and trend roles.
  collision = d
  collision[["n"]] = collision[["g"]]
  collision[["sum_w"]] = collision[["z"]]
  fit = dbreg(
    y ~ x + n | n[sum_w] + h,
    data = collision,
    weights = ~w,
    offset = ~o,
    strategy = "compress"
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  expect_equal(
    predict(fit, collision),
    as.numeric(predict(ref)) + d[["o"]],
    tolerance = 1e-8
  )
  expect_equal(anyDuplicated(names(fit[["data"]])), 0L)
  cells = dbreg(
    y ~ x | n[sum_w] + h,
    data = collision,
    weights = ~w,
    strategy = "compress",
    data_only = TRUE
  )
  expect_true(all(
    c("n", "n_source", "sum_w", "sum_w_source") %in% names(cells)
  ))
  expect_equal(sum(cells[["n"]]), n)
  expect_equal(sum(cells[["sum_w"]]), sum(collision[["w"]]), tolerance = 1e-10)
  expect_equal(
    sum(cells[["sum_wy"]]),
    sum(collision[["w"]] * collision[["y"]]),
    tolerance = 1e-10
  )
  unadjusted = dbreg(
    y ~ x | n[sum_w] + h,
    data = collision,
    weights = ~w,
    strategy = "compress"
  )
  expect_equal(
    predict(unadjusted, cells),
    predict(unadjusted),
    tolerance = 1e-8
  )
  changed = unadjusted[["data"]]
  changed[["x"]] = as.numeric(changed[["x"]]) + 1
  expect_equal(
    predict(unadjusted, changed) - predict(unadjusted),
    rep(coef(unadjusted)[[1L]], nrow(changed)),
    tolerance = 1e-8
  )
  # Exact typed keys and tuple ranks survive both prediction adapters.
  DBI::dbExecute(
    conn,
    "CREATE TEMP TABLE typed_cv AS SELECT * EXCLUDE(g,h),CAST('9007199254740992' AS BIGINT)+CAST(g AS BIGINT) AS g,'a,b;|/'||CAST(h AS VARCHAR) AS h FROM raw_cv"
  )
  fit = dbreg(
    y ~ x | (g^h)[z],
    conn = conn,
    table = "typed_cv",
    weights = ~w,
    offset = ~o,
    strategy = "compress",
    vcov = ~g
  )
  ref = fixest::feols(
    y ~ x | (g^h)[z],
    adjusted,
    weights = ~w,
    vcov = ~g,
    ssc = fixest::ssc(K.fixef = "full"),
    fixef.rm = "none",
    fixef.tol = 1e-11,
    notes = FALSE
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  raw = d
  raw[["g"]] = paste0("900719925474", sprintf("%04d", 992 + d[["g"]]))
  raw[["h"]] = paste0("a,b;|/", d[["h"]])
  expect_equal(
    predict(fit, raw),
    as.numeric(predict(ref)) + d[["o"]],
    tolerance = 1e-8
  )
  DBI::dbRemoveTable(conn, "typed_cv")
  # Exact null trends preserve coefficient values and declared dimension counts.
  constant = d
  constant[["copy"]] = 2 * constant[["z"]]
  constant[["constant"]] = 1
  fit = dbreg(
    y ~ x | g[z, copy, constant] + h,
    data = constant,
    weights = ~w,
    offset = ~o,
    strategy = "compress"
  )
  ref = fixest::feols(
    y ~ x | g[z, copy, constant] + h,
    transform(constant, y = y - o),
    weights = ~w,
    ssc = fixest::ssc(K.fixef = "full"),
    fixef.rm = "none",
    fixef.tol = 1e-11,
    notes = FALSE
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  for (cap in c(
    "dbreg.projection_max_rows",
    "dbreg.projection_max_cells",
    "dbreg.projection_max_entries",
    "dbreg.projection_max_bytes"
  )) {
    prior = options(setNames(list(1), cap))
    expect_error(
      dbreg(
        y ~ x | g[z] + h,
        conn = conn,
        table = "raw_cv",
        weights = ~w,
        strategy = "compress"
      ),
      "budget"
    )
    options(prior)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})


local({
  set.seed(11130)
  n = 360L
  d = data.frame(
    x = sample(-1:1, n, TRUE),
    z = sample(-2:2, n, TRUE),
    g = sample(1:8, n, TRUE),
    q = 0,
    h = sample(1:6, n, TRUE)
  )
  d[["y"]] = .7 *
    d[["x"]] +
    d[["g"]] / 10 * d[["z"]] +
    sin(d[["g"]]) +
    cos(d[["h"]]) +
    rnorm(n)
  fit = dbreg(y ~ x | g[z, q] + h, data = d, strategy = "compress")
  bad = d[1:4, ]
  bad[["q"]][1L] = 1
  expect_true(is.na(predict(fit, bad)[1L]))
  expect_true(all(is.finite(predict(fit, bad)[-1L])))
  # Producer-looking raw names have no authority over private source bindings.
  d[["cf1"]] = d[["x"]]
  d[["ct1"]] = d[["z"]]
  fit = dbreg(y ~ cf1 | g[ct1] + h, data = d, strategy = "compress")
  ref = fixest::feols(
    y ~ x | g[z] + h,
    d,
    ssc = fixest::ssc(K.fixef = "full"),
    fixef.rm = "none",
    fixef.tol = 1e-11,
    notes = FALSE
  )
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(predict(fit, d), as.numeric(predict(ref)), tolerance = 1e-8)
})

local({
  set.seed(11131)
  d = expand.grid(f = 1:4, g = 1:2, z = c(-1, 1), repeat_id = 1:8)
  d[["y"]] = .8 *
    (d[["f"]] == 2) +
    d[["g"]] / 5 * d[["z"]] +
    d[["g"]] +
    rnorm(nrow(d))
  fit = dbreg(y ~ i(f, keep = 2) | g[z], data = d, strategy = "compress")
  ref = fixest::feols(
    y ~ i(f, keep = 2) | g[z],
    d,
    fixef.rm = "none",
    ssc = fixest::ssc(K.fixef = "full"),
    notes = FALSE
  )
  expect_equal(fit[["nobs"]], 8L)
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
  expect_equal(predict(fit, d), as.numeric(predict(ref)), tolerance = 1e-8)
  expect_equal(predict(fit), predict(fit, fit[["data"]]), tolerance = 1e-8)
})

local({
  set.seed(11132)
  n = 300L
  index = sample(1:3, n, TRUE)
  d = data.frame(
    x = sample(-1:1, n, TRUE),
    g = (1 + c(0, 1, 2) * .Machine[["double.eps"]])[index],
    z = sample(-2:2, n, TRUE),
    w = runif(n, .2, 2)
  )
  d[["y"]] = .7 * d[["x"]] + sin(index) + index / 10 * d[["z"]] + rnorm(n)
  fit = dbreg(y ~ x | g[z], data = d, weights = ~w, strategy = "compress")
  D = model.matrix(~ x + factor(index) + factor(index):z, d)
  sw = sqrt(d[["w"]])
  b = qr.coef(qr(D * sw), d[["y"]] * sw)
  expect_equal(as.numeric(coef(fit)), as.numeric(b[2]), tolerance = 1e-8)
  expect_equal(predict(fit, d), as.numeric(D %*% b), tolerance = 1e-8)
  expect_equal(length(unique(fit[["data"]][["g"]])), 3L)
})

local({
  set.seed(11133)
  n = 400L
  d = data.frame(
    x = sample(-1:1, n, TRUE),
    g = sample(1:4, n, TRUE),
    z = sample(-2:2, n, TRUE)
  )
  d[["h"]] = ifelse(d[["g"]] <= 2, sample(1:2, n, TRUE), sample(3:4, n, TRUE))
  d[["y"]] = .8 *
    d[["x"]] +
    d[["g"]] / 10 * d[["z"]] +
    sin(d[["g"]]) +
    cos(d[["h"]]) +
    rnorm(n)
  fit = dbreg(y ~ x | g[z] + h, data = d, strategy = "compress")
  expect_true(all(is.finite(predict(fit, d))))
  raw = d[1:2, ]
  raw[["g"]][1] = 1
  raw[["h"]][1] = 4
  expect_true(is.na(predict(fit, raw)[1L]))
})

local({
  d = data.frame(
    x = rep(-1:1, 10),
    g = rep(c("short", paste(rep("x", 65535), collapse = "")), 15),
    z = rep(-2:2, 6),
    y = seq_len(30)
  )
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "wide_vary_key", d)
  before = sort(DBI::dbListTables(conn))
  expect_error(
    dbreg(
      y ~ x | g[z],
      conn = conn,
      table = "wide_vary_key",
      strategy = "compress"
    ),
    "byte budget"
  )
  expect_identical(sort(DBI::dbListTables(conn)), before)
})

local({
  set.seed(11134)
  d = expand.grid(x = c(-1, 1), g = 1:2, z = c(-1, 1), repetition = 1:16)
  d[["w"]] = runif(nrow(d), .2, 2)
  d[["y"]] = .8 * d[["x"]] + d[["g"]] * d[["z"]] + rnorm(nrow(d))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "cv_plan", d)
  initial = sort(DBI::dbListTables(conn))
  fit = dbreg(y ~ x | g[z], conn = conn, table = "cv_plan", weights = ~w)
  expect_identical(fit[["strategy"]], "compress")
  prior = options(dbreg.projection_max_rows = 7)
  fallback = dbreg(y ~ x | g[z], conn = conn, table = "cv_plan", weights = ~w)
  options(prior)
  expect_identical(fallback[["strategy"]], "demean")
  expect_equal(
    as.numeric(coef(fallback)),
    as.numeric(coef(fit)),
    tolerance = 1e-8
  )
  expect_identical(
    fallback[["projection_diagnostics"]][[1L]][["reason"]],
    "system_limit"
  )
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # SQL-only cells remain executable after the planner has returned.
  invisible(capture.output({
    sql = dbreg(
      y ~ x | g[z],
      conn = conn,
      table = "cv_plan",
      weights = ~w,
      strategy = "compress",
      sql_only = TRUE
    )
  }))
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  expect_equal(sum(DBI::dbGetQuery(conn, sql)[["n"]]), nrow(d))
  ordinary = dbreg(
    y ~ x | g[z],
    conn = conn,
    table = "cv_plan",
    weights = ~w,
    vcov = ~g,
    strategy = "compress"
  )
  prior = options(dbreg.projection_max_upload_bytes = 1)
  limited = dbreg(
    y ~ x | g[z],
    conn = conn,
    table = "cv_plan",
    weights = ~w,
    vcov = ~g,
    strategy = "compress"
  )
  expect_equal(
    as.numeric(vcov(limited)),
    as.numeric(vcov(ordinary)),
    tolerance = 1e-8
  )
  expect_identical(
    ordinary[["compression_diagnostics"]][["cluster_score"]],
    "sql"
  )
  expect_identical(
    limited[["compression_diagnostics"]][["cluster_score"]],
    "cluster_cells"
  )
  prior_rows = options(dbreg.projection_max_rows = 9)
  expect_error(
    dbreg(
      y ~ x | g[z],
      conn = conn,
      table = "cv_plan",
      weights = ~w,
      vcov = ~g,
      strategy = "compress"
    ),
    "cluster cells exceed"
  )
  options(prior_rows)
  options(prior)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
