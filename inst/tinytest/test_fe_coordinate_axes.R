local({
  set.seed(11021)
  n = 240L
  d = data.frame(grp = rep(c("a", "b"), n / 2L), grpb = stats::rnorm(n),
    g = rep(seq_len(8), each = 30), cl = rep(seq_len(12), 20))
  Z = stats::model.matrix(~factor(g), d)
  X = cbind(Z[, 1], indicator_b = as.numeric(d[["grp"]] == "b"), raw_grpb = d[["grpb"]], Z[, -1, drop = FALSE])
  d[["y"]] = as.numeric(X %*% c(1, 2, -0.8, seq_len(7) / 7)) + stats::rnorm(n, sd = 0.3)
  beta = as.numeric(qr.solve(X, d[["y"]]))
  residual = as.numeric(d[["y"]] - X %*% beta)
  inverse = solve(crossprod(X))
  slopes = 2:3
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  models = list()
  on.exit({for (model in models) close(model); DBI::dbDisconnect(conn, shutdown = TRUE)}, add = TRUE)
  for (strategy in c("compress", "demean")) for (policy in c("full", "none")) for (kind in c("iid", "hc1", "cluster")) {
    request = if (kind == "cluster") ~cl else kind
    K = if (policy == "full") ncol(X) else 2L
    expected = if (kind == "iid") inverse * sum(residual^2) / (n - K) else if (kind == "hc1") {
      inverse %*% crossprod(X * residual) %*% inverse * n / (n - K)
    } else {
      scores = rowsum(X * residual, d[["cl"]])
      inverse %*% crossprod(scores) %*% inverse * 12 / 11 * (n - 1) / (n - K)
    }
    model = dbreg::dbreg(y ~ grp + grpb | g, conn = conn, data = d,
      strategy = strategy, vcov = request, ssc = list(K.fixef = policy), retain = "inference")
    models[[length(models) + 1L]] = model
    expect_equal(as.numeric(stats::coef(model)), beta[slopes], tolerance = 1e-8)
    expect_equal(length(stats::coef(model)), 2L)
    expect_equal(as.numeric(stats::vcov(model)), as.numeric(expected[slopes, slopes]), tolerance = 1e-8)
    expect_equal(stats::df.residual(model), n - ncol(X))
    expect_equal(attr(stats::vcov(model), "df.K"), K)
    if (strategy == "compress") {
      expect_equal(as.numeric(stats::coef(model, fe = TRUE)), beta, tolerance = 1e-8)
      expect_equal(as.numeric(stats::vcov(model, fe = TRUE)), as.numeric(expected), tolerance = 1e-8)
      cells = model[["data"]]
      cell_Z = stats::model.matrix(~factor(g), cells)
      cell_X = cbind(cell_Z[, 1], as.numeric(cells[["grp"]] == "b"), cells[["grpb"]], cell_Z[, -1, drop = FALSE])
      expect_equal(stats::predict(model), as.numeric(cell_X %*% beta), tolerance = 1e-8)
    }
    errors = sqrt(diag(expected)[slopes])
    critical = stats::qt(0.975, if (kind == "cluster") 11 else n - K)
    expected_ci = cbind(beta[slopes] - critical * errors, beta[slopes] + critical * errors)
    expect_equal(as.numeric(stats::confint(model)), as.numeric(expected_ci), tolerance = 1e-8)
    expect_equal(as.numeric(stats::confint(model, parm = 2)), as.numeric(expected_ci[2, ]), tolerance = 1e-8)
    table = generics::tidy(model, conf.int = TRUE)
    expect_equal(table[["estimate"]], beta[slopes], tolerance = 1e-8)
    expect_identical(table[["term"]], c("grpb", "grpb"))
    expect_equal(as.numeric(table[["conf.low"]]), as.numeric(expected_ci[, 1]), tolerance = 1e-8)
    printed = capture.output({ value = print(model) })
    expect_equal(as.numeric(value[, "Estimate"]), beta[slopes], tolerance = 1e-8)
    changed = summary(model, vcov = "hc1")
    expected_hc1 = inverse %*% crossprod(X * residual) %*% inverse * n / (n - K)
    expect_equal(as.numeric(stats::coef(changed)), beta[slopes], tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(changed)), as.numeric(expected_hc1[slopes, slopes]), tolerance = 1e-8)
    expect_equal(as.numeric(stats::coef(changed, fe = TRUE)), as.numeric(stats::coef(model, fe = TRUE)), tolerance = 1e-10)
    expect_equal(nrow(stats::vcov(changed)), 2L)
    expect_error(stats::vcov(changed, fe = TRUE), "Post-fit absorbed-effect covariance")
    if (requireNamespace("fixest", quietly = TRUE)) {
      expect_equal(as.numeric(fixest::coeftable(model)[, 1]), beta[slopes], tolerance = 1e-8)
      expect_equal(as.numeric(fixest::coeftable(model, keep = "grpb")[, 1]), beta[slopes], tolerance = 1e-8)
    }
  }
})

local({
  set.seed(11022)
  d = data.frame(grp = rep(c("a", "b"), 90), grpb = stats::rnorm(180), g = rep(1:6, each = 30))
  d[["copy"]] = d[["grpb"]]
  d[["y"]] = 1 + 2 * (d[["grp"]] == "b") - 0.8 * d[["grpb"]] + sin(d[["g"]]) + stats::rnorm(180)
  for (strategy in c("compress", "demean")) {
    model = dbreg::dbreg(y ~ grp + copy + grpb | g, data = d, strategy = strategy)
    coefficients = stats::coef(model, collin = TRUE)
    expect_equal(length(coefficients), 3L)
    expect_equal(sum(is.na(coefficients)), 1L)
    table = generics::tidy(model, collin = TRUE, conf.int = TRUE)
    expect_equal(nrow(table), 3L)
    expect_equal(sum(is.na(table[["estimate"]])), 1L)
    expect_equal(sum(is.na(table[["conf.low"]])), 1L)
    close(model)
  }
})

local({
  V = Matrix::Diagonal(x = c(4, 9))
  dimnames(V) = rep(list(c("same", "same")), 2L)
  record = dbreg:::new_inference_result(V, "iid", c("same", "same"), list(), NULL,
    dbreg:::inference_ssc("full"), 2, 7, "t", coefficient_ids = c("first", "second"))
  expect_identical(class(dbreg:::inference_covariance(record)), class(V))
  subset = dbreg:::inference_covariance(record, coefficient_ids = "second")
  expect_true(inherits(subset, "Matrix"))
  expect_equal(as.numeric(subset), 9)
  expect_identical(attr(subset, "dbreg_coefficient_ids"), "second")
  empty = dbreg:::inference_covariance(record, coefficient_ids = character())
  expect_identical(dim(empty), c(0L, 0L))
  expect_equal(attr(empty, "df.t"), 7)
})
