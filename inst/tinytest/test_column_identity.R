local({
  set.seed(11041)
  d = data.frame(
    a = rnorm(160),
    b = rnorm(160),
    a_x_b = rnorm(160),
    dbreg_i1_1 = rnorm(160),
    f = rep(c("a", "b"), 80),
    w = runif(160, .4, 2),
    cl = rep(1:16, 10)
  )
  d[["y"]] = .7 +
    .8 * d[["a"]] * d[["b"]] -
    .6 * d[["a_x_b"]] +
    .5 * (d[["f"]] == "b") +
    rnorm(160)
  cases = list(
    y ~ a:b + a_x_b,
    y ~ a_x_b + a:b,
    y ~ i(f, ref = "a") + dbreg_i1_1,
    y ~ dbreg_i1_1 + i(f, ref = "a")
  )
  design = function(formula, data) {
    terms = attr(terms(formula), "term.labels")
    columns = lapply(terms, function(term) {
      if (grepl("i(", term, fixed = TRUE)) {
        return(as.numeric(data[["f"]] == "b"))
      }
      if (term == "a:b") {
        return(data[["a"]] * data[["b"]])
      }
      data[[term]]
    })
    do.call(cbind, c(list(rep(1, nrow(data))), columns))
  }
  for (formula in cases) {
    for (strategy in c("moments", "compress")) {
      for (kind in list("iid", "hc1", ~cl)) {
        X = design(formula, d)
        w = d[["w"]]
        beta = as.numeric(qr.solve(X * sqrt(w), d[["y"]] * sqrt(w)))
        residual = d[["y"]] - as.numeric(X %*% beta)
        bread = solve(crossprod(X, w * X))
        meat = if (identical(kind, "hc1")) {
          crossprod(X * (w * residual))
        } else if (inherits(kind, "formula")) {
          crossprod(rowsum(X * (w * residual), d[["cl"]]))
        } else {
          NULL
        }
        V = if (is.null(meat)) {
          bread * sum(w * residual^2) / (160 - ncol(X))
        } else {
          bread %*%
            meat %*%
            bread *
            if (inherits(kind, "formula")) {
              16 / 15 * 159 / (160 - ncol(X))
            } else {
              160 / (160 - ncol(X))
            }
        }
        fit = tryCatch(
          dbreg(
            formula,
            data = d,
            weights = "w",
            vcov = kind,
            strategy = strategy
          ),
          error = function(e) e
        )
        expect_false(
          inherits(fit, "error"),
          info = paste(deparse(formula), strategy, deparse(kind))
        )
        if (inherits(fit, "error")) {
          next
        }
        expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
        expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-8)
        prediction = tryCatch(predict(fit, d), error = function(e) e)
        expect_false(inherits(prediction, "error"))
        if (!inherits(prediction, "error")) {
          expect_equal(
            as.numeric(prediction),
            as.numeric(X %*% beta),
            tolerance = 1e-8
          )
        }
        if (strategy == "compress") {
          prediction = tryCatch(predict(fit), error = function(e) e)
          expect_false(inherits(prediction, "error"))
          if (!inherits(prediction, "error")) {
            # Compressed rows follow cell order rather than input order.
            stored = fit[["data"]]
            expect_equal(length(prediction), nrow(stored))
            expect_equal(
              sort(as.numeric(prediction)),
              sort(as.numeric(X %*% beta)),
              tolerance = 1e-8
            )
          }
        }
      }
    }
  }
  # Source roles outside the RHS reserve their names beside generated columns.
  for (role in c("response", "weight", "cluster")) {
    for (strategy in c("moments", "compress")) {
      data = d
      data[["a_x_b"]] = switch(
        role,
        response = d[["y"]],
        weight = d[["w"]],
        cluster = d[["cl"]]
      )
      formula = if (role == "response") a_x_b ~ a:b else y ~ a:b
      X = cbind(1, data[["a"]] * data[["b"]])
      w = if (role == "weight") data[["a_x_b"]] else rep(1, 160)
      y = if (role == "response") data[["a_x_b"]] else data[["y"]]
      beta = as.numeric(qr.solve(X * sqrt(w), y * sqrt(w)))
      args = list(formula, data = data, strategy = strategy)
      if (role == "weight") {
        args[["weights"]] = "a_x_b"
      }
      if (role == "cluster") {
        args[["vcov"]] = ~a_x_b
      }
      fit = tryCatch(do.call(dbreg, args), error = function(e) e)
      expect_false(inherits(fit, "error"), info = paste(role, strategy))
      if (!inherits(fit, "error")) {
        expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
      }
    }
  }
})

local({
  set.seed(11042)
  d = data.frame(
    a = rnorm(180),
    b_x_c = rnorm(180),
    a_x_b = rnorm(180),
    c = rnorm(180),
    z = rnorm(180)
  )
  d[["endog"]] = .7 * d[["a_x_b"]] * d[["c"]] + .4 * d[["z"]] + rnorm(180)
  d[["y"]] = 1 + .4 * d[["a"]] * d[["b_x_c"]] + 1.2 * d[["endog"]] + rnorm(180)
  X = cbind(1, d[["a"]] * d[["b_x_c"]], d[["endog"]])
  Z = cbind(1, d[["a"]] * d[["b_x_c"]], d[["a_x_b"]] * d[["c"]], d[["z"]])
  projection = crossprod(X, Z) %*% solve(crossprod(Z))
  bread = solve(projection %*% crossprod(Z, X))
  beta = as.numeric(bread %*% projection %*% crossprod(Z, d[["y"]]))
  V = bread * sum((d[["y"]] - X %*% beta)^2) / 177
  for (strategy in c("moments", "compress")) {
    fit = tryCatch(
      dbivreg(y ~ a:b_x_c | endog ~ a_x_b:c + z, data = d, strategy = strategy),
      error = function(e) e
    )
    expect_false(inherits(fit, "error"), info = strategy)
    if (inherits(fit, "error")) {
      next
    }
    expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-8)
    prediction = tryCatch(predict(fit, d), error = function(e) e)
    expect_false(inherits(prediction, "error"))
    if (!inherits(prediction, "error")) {
      expect_equal(
        as.numeric(prediction),
        as.numeric(X %*% beta),
        tolerance = 1e-8
      )
    }
  }
})

# Column construction is positional, including repeated public labels, trained
# category coding, dropped terms, row positions and interval covariance.
local({
  set.seed(1702)
  d = data.frame(
    x = rnorm(120),
    z = rnorm(120),
    f = rep(c("a", "b", "c"), 40),
    fb = rnorm(120)
  )
  d[["y"]] = .6 +
    .4 * d[["x"]] +
    .8 * (d[["f"]] == "b") -
    .3 * (d[["f"]] == "c") +
    .2 * d[["fb"]] +
    rnorm(120)
  X = cbind(
    1,
    d[["x"]],
    as.numeric(d[["f"]] == "b"),
    as.numeric(d[["f"]] == "c"),
    d[["fb"]]
  )
  beta = as.numeric(qr.solve(X, d[["y"]]))
  V = solve(crossprod(X)) * sum((d[["y"]] - X %*% beta)^2) / 115
  rows = d[c(6, 1, 2, 8, 4), ]
  rows[["f"]] = factor(rows[["f"]], levels = c("c", "b", "a"))
  rows[["x"]][3] = NA_real_
  rows[["irrelevant"]] = NA_real_
  newX = cbind(
    1,
    rows[["x"]],
    as.numeric(rows[["f"]] == "b"),
    as.numeric(rows[["f"]] == "c"),
    rows[["fb"]]
  )
  expected = as.numeric(newX %*% beta)
  se = sqrt(rowSums((newX %*% V) * newX))
  for (strategy in c("moments", "compress")) {
    fit = dbreg(y ~ x + f + fb, data = d, strategy = strategy)
    expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
    expect_identical(names(coef(fit)), c("(Intercept)", "x", "fb", "fc", "fb"))
    expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-8)
    expect_equal(predict(fit, rows), expected, tolerance = 1e-8)
    interval = predict(fit, rows, interval = "confidence")
    expect_equal(
      interval[["lwr"]],
      expected - qt(.975, 115) * se,
      tolerance = 1e-8
    )
    expect_equal(
      interval[["upr"]],
      expected + qt(.975, 115) * se,
      tolerance = 1e-8
    )
    unknown = rows
    unknown[["f"]] = as.character(unknown[["f"]])
    unknown[["f"]][1] = "new category"
    expect_true(is.na(predict(fit, unknown)[1]))
    changed = rows
    changed[["x"]] = factor(changed[["x"]])
    expect_error(predict(fit, changed), "must be numeric")
    expect_equal(length(predict(fit, rows[FALSE, , drop = FALSE])), 0)
    broken = fit
    broken[["linear_design"]][["fitted_ids"]][1] = "missing identity"
    expect_error(predict(broken, rows), "Invalid fitted column")
  }
  # R's logical matrix basis includes both FALSE and TRUE, even when only TRUE
  # is observed. Its no-intercept promotion is a producer decision.
  d[["flag"]] = TRUE
  fit = dbreg(y ~ 0 + flag + x, data = d, strategy = "compress")
  ref = lm(y ~ 0 + flag + x, d)
  expect_equal(
    as.numeric(coef(fit)),
    as.numeric(coef(ref)[!is.na(coef(ref))]),
    tolerance = 1e-8
  )
  expect_equal(predict(fit, d), as.numeric(predict(ref, d)), tolerance = 1e-8)
})

local({
  set.seed(1703)
  d = data.frame(
    a = rnorm(120),
    b = rnorm(120),
    a_b = rnorm(120),
    check.names = FALSE
  )
  d[["a:b"]] = rnorm(120)
  d[["y"]] = .7 +
    .3 * d[["a:b"]] +
    .5 * d[["a_b"]] +
    .8 * d[["a"]] * d[["b"]] +
    rnorm(120)
  X = cbind(1, d[["a:b"]], d[["a_b"]], d[["a"]] * d[["b"]])
  beta = as.numeric(qr.solve(X, d[["y"]]))
  for (strategy in c("moments", "compress")) {
    fit = dbreg(y ~ `a:b` + a_b + a:b, data = d, strategy = strategy)
    expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
    expect_identical(names(coef(fit)), c("(Intercept)", "`a:b`", "a_b", "a:b"))
    expect_equal(predict(fit, d), as.numeric(X %*% beta), tolerance = 1e-8)
  }
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "quoted_columns", d)
  compiled = sql_model_matrix(~ `a:b` + a_b + a:b, conn, "quoted_columns")
  expect_identical(compiled[["aliases"]][1:2], c("a_b_1", "a_b"))
  actual = DBI::dbGetQuery(
    conn,
    paste(
      "SELECT",
      paste(compiled[["select_exprs"]], collapse = ","),
      "FROM quoted_columns"
    )
  )
  expect_equal(unname(as.matrix(actual)), unname(X[, -1]), tolerance = 1e-12)
  lhs = sql_model_matrix(Formula::Formula(a_x_b ~ a:b), conn, "quoted_columns")
  expect_identical(lhs[["aliases"]], "a_x_b_1")
  reserved = sql_model_matrix(
    ~ a:b,
    conn,
    "quoted_columns",
    reserved_names = "a_x_b"
  )
  expect_identical(reserved[["aliases"]], "a_x_b_1")
})

local({
  saved = options("contrasts")
  on.exit(options(saved), add = TRUE)
  set.seed(1704)
  d = data.frame(x = rep(1:4, 30), f = rep(c("a", "b", "c"), 40))
  d[["y"]] = .4 +
    .7 * d[["x"]] +
    .6 * (d[["f"]] == "b") -
    .2 * (d[["f"]] == "c") +
    rnorm(120)
  X = cbind(
    1,
    d[["x"]],
    as.numeric(d[["f"]] == "b"),
    as.numeric(d[["f"]] == "c")
  )
  beta = as.numeric(qr.solve(X, d[["y"]]))
  options(contrasts = c("contr.sum", "contr.poly"))
  fit = dbreg(y ~ x + f, data = d, strategy = "compress")
  expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
  expect_equal(predict(fit, d), as.numeric(X %*% beta), tolerance = 1e-8)
  options(contrasts = c("contr.helmert", "contr.poly"))
  cells = fit[["data"]]
  expected = cbind(
    1,
    cells[["x"]],
    as.numeric(cells[["f"]] == "b"),
    as.numeric(cells[["f"]] == "c")
  ) %*%
    beta
  expect_equal(as.numeric(predict(fit)), as.numeric(expected), tolerance = 1e-8)
  expect_equal(length(fit[["linear_design"]][["dictionaries"]]), 1)
})

local({
  set.seed(11042)
  n = 180L
  d = data.frame(
    a = rnorm(n),
    b_x_c = rnorm(n),
    a_x_b = rnorm(n),
    c = rnorm(n),
    dbreg_i = rep(c("0", "1"), 90),
    g = rep(c("0", "1"), 90),
    dbreg_i1 = rnorm(n),
    f = rep(c("a", "b", "c"), 60),
    z = rnorm(n),
    w = runif(n, .4, 2),
    cl = rep(1:15, 12)
  )
  d[["endog"]] = .7 * d[["a_x_b"]] * d[["c"]] + .4 * d[["z"]] + rnorm(n)
  d[["y"]] = 1 + .4 * d[["a"]] * d[["b_x_c"]] + 1.2 * d[["endog"]] + rnorm(n)
  formulas = list(
    mixed = y ~ dbreg_i | endog ~ i(f, ref = "a") + z,
    factor_i = y ~ g + i(f, ref = "a") | endog ~ z,
    raw = y ~ dbreg_i1 | endog ~ i(f, ref = "a") + z,
    overlap = y ~ i(f, ref = "a") + dbreg_i1 | endog ~ z + i(f, ref = "a")
  )
  for (case in names(formulas)) {
    exog = switch(
      case,
      mixed = as.numeric(d[["dbreg_i"]] == "1"),
      raw = d[["dbreg_i1"]],
      factor_i = cbind(
        as.numeric(d[["g"]] == "1"),
        as.numeric(d[["f"]] == "b"),
        as.numeric(d[["f"]] == "c")
      ),
      overlap = cbind(
        as.numeric(d[["f"]] == "b"),
        as.numeric(d[["f"]] == "c"),
        d[["dbreg_i1"]]
      )
    )
    instruments = if (case %in% c("factor_i", "overlap")) {
      d[["z"]]
    } else {
      cbind(as.numeric(d[["f"]] == "b"), as.numeric(d[["f"]] == "c"), d[["z"]])
    }
    X = cbind(1, exog, d[["endog"]])
    Z = cbind(1, exog, instruments)
    w = d[["w"]]
    P = crossprod(X, w * Z) %*% solve(crossprod(Z, w * Z))
    bread = solve(P %*% crossprod(Z, w * X))
    beta = as.numeric(bread %*% P %*% crossprod(Z, w * d[["y"]]))
    residual = d[["y"]] - as.numeric(X %*% beta)
    for (kind in list("iid", "hc1", ~cl)) {
      meat = if (identical(kind, "hc1")) {
        crossprod(Z * (w * residual))
      } else if (inherits(kind, "formula")) {
        crossprod(rowsum(Z * (w * residual), d[["cl"]]))
      } else {
        NULL
      }
      V = if (is.null(meat)) {
        bread * sum(w * residual^2) / (n - ncol(X))
      } else {
        bread %*%
          P %*%
          meat %*%
          t(P) %*%
          bread *
          if (inherits(kind, "formula")) {
            15 / 14 * 179 / (n - ncol(X))
          } else {
            180 / (n - ncol(X))
          }
      }
      for (strategy in c("moments", "compress")) {
        if (case == "overlap") {
          expect_warning(
            {
              fit = dbivreg(
                formulas[[case]],
                data = d,
                weights = "w",
                vcov = kind,
                strategy = strategy
              )
            },
            "already valid instruments"
          )
        } else {
          fit = dbivreg(
            formulas[[case]],
            data = d,
            weights = "w",
            vcov = kind,
            strategy = strategy
          )
        }
        expect_equal(
          as.numeric(coef(fit)),
          beta,
          tolerance = 1e-8,
          info = paste(case, strategy, deparse(kind))
        )
        expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-8)
        expect_equal(predict(fit, d), as.numeric(X %*% beta), tolerance = 1e-8)
        if (strategy == "compress") {
          expect_equal(
            sort(predict(fit)),
            sort(as.numeric(X %*% beta)),
            tolerance = 1e-8
          )
        }
      }
    }
  }
  expect_error(
    dbivreg(
      y ~ f | endog ~ i(f, ref = "a") + z,
      data = d,
      strategy = "moments"
    ),
    "Instrument matrix is rank deficient"
  )
  expect_error(
    dbivreg(y ~ 1 | endog ~ endog, data = d, strategy = "moments"),
    "Excluded instruments cannot also appear"
  )
})

# A point-prediction matrix with many categories must remain sparse rather than
# allocating a dense N-by-G buffer (including missing/unknown row handling).
local({
  set.seed(1705)
  d = data.frame(x = rnorm(320), f = rep(sprintf("g%03d", 1:80), 4))
  d[["y"]] = .4 *
    d[["x"]] +
    sin(match(d[["f"]], sort(unique(d[["f"]])))) +
    rnorm(320)
  fit = dbreg(y ~ x + f, data = d, strategy = "moments")
  rows = d[rep(1:320, length.out = 10000), , drop = FALSE]
  rows[["f"]][1] = NA_character_
  rows[["f"]][2] = "unknown"
  mm = dbreg:::linear_prediction_matrix(fit, rows)
  expect_true(inherits(mm, "sparseMatrix"))
  expect_true(Matrix::nnzero(mm) <= 3 * nrow(rows))
  prediction = predict(fit, rows)
  expect_true(all(is.na(prediction[1:2])))
  expect_equal(length(prediction), 10000)
  reference = lm(
    y ~ x + f,
    d,
    contrasts = list(f = contr.treatment(sort(unique(d[["f"]]))))
  )
  expect_equal(
    prediction[3:10000],
    as.numeric(predict(reference, rows[3:10000, ])),
    tolerance = 1e-8
  )
})

# A literal category named _constant is distinct from the compiler's appended
# singleton sentinel; its target must not absorb a raw predictor of that name.
local({
  set.seed(1706)
  d = data.frame(
    f = rep("_constant", 60),
    f_constant = rnorm(60),
    y = rnorm(60)
  )
  X = cbind(1, d[["f_constant"]])
  beta = as.numeric(qr.solve(X, d[["y"]]))
  fit = dbreg(y ~ 0 + f + f_constant, data = d, strategy = "moments")
  expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-8)
  expect_equal(predict(fit, d), as.numeric(X %*% beta), tolerance = 1e-8)
  unknown = d
  unknown[["f"]][1] = "other"
  expect_true(is.na(predict(fit, unknown)[1]))
})
