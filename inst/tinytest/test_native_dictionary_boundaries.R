# Four focused independent value witnesses for producer/role boundary amendments.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  index = seq_len(120)
  fit_case = function(formula, data, matrix, slopes, strategy, learned, change, inspect) {
    old = options(dbreg.prediction_max_entries = 100000)
    on.exit(options(old))
    X = matrix(data)
    reference = stats::lm.wfit(X, data[["y"]], data[["w"]])
    expect_identical(reference[["rank"]], ncol(X))
    beta = reference[["coefficients"]]
    V = solve(crossprod(X, data[["w"]] * X)) *
      sum(data[["w"]] * reference[["residuals"]]^2) / (nrow(X) - ncol(X))
    model = dbreg(formula, data = data, conn = conn, weights = "w", strategy = strategy,
      retain = if (learned) "prediction" else "none", ssc = list(K.fixef = "full"))
    on.exit(close(model), add = TRUE, after = FALSE)
    expect_equal(as.numeric(coef(model)), as.numeric(beta[slopes]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(model)), as.numeric(V[slopes, slopes, drop = FALSE]), tolerance = 1e-7)
    rows = change(data[c(17, 2, 55, 9, 37, 79), ])
    rows[["y"]] = -999
    expected = as.numeric(matrix(rows) %*% beta)
    r = predict(model, rows)
    expect_identical(is.na(r), is.na(expected))
    expect_equal(as.numeric(r), expected, tolerance = 1e-7)
    out = predict(model, rows, conn = conn, row_id = "id")
    table = as.data.frame(out)
    actual = table[[out[["prediction_column"]]]][match(rows[["id"]], table[["id"]])]
    close(out)
    expect_identical(is.na(actual), is.na(expected))
    expect_equal(as.numeric(actual), expected, tolerance = 1e-7)
    inspect(model, rows)
  }
  base = data.frame(id = index, w = 1 + (index %% 5) / 4,
    f = rep(c(1, 1, 2, 2), 30), g = rep(seq_len(20), each = 6),
    x = sin(index * .19) + cos(index * .37))
  base[["factor(f)"]] = rep(c("a", "b"), 60)
  base[["trend value"]] = sin(index * .43) + .1 * cos(index * .61)

  # Trend metadata completeness is checked against actual varying_slopes.
  data = base
  trend_X = function(rows) {
    D = outer(rows[["g"]], seq_len(20), `==`) * 1
    cbind(as.numeric(rows[["factor(f)"]] == "b"), as.numeric(rows[["f"]] == 2), D, D * rows[["trend value"]])
  }
  data[["y"]] = as.numeric(trend_X(data) %*% c(.7, -.4, sin(seq_len(20)), .3 * cos(seq_len(20)))) + .03 * sin(index * .71)
  fit_case(y ~ `factor(f)` + factor(f) | g[`trend value`], data, trend_X, 1:2, "demean", TRUE,
    change = function(rows) { rows[["trend value"]] = rows[["trend value"]] + .15; rows },
    inspect = function(model, rows) {
      expect_identical(unname(model[["linear_design"]][["trend_sources"]]), "trend value")
      broken = model
      broken[["linear_design"]][["trend_sources"]] = model[["linear_design"]][["trend_sources"]][-1L]
      expect_error(predict(broken, rows), "descriptor|source")
      expect_error(predict(broken, rows, conn = conn, row_id = "id"), "descriptor|source")
      broken = model
      broken[["linear_design"]][["fe_sources"]] = model[["linear_design"]][["fe_sources"]][-1L]
      expect_error(predict(broken, rows), "descriptor|source")
      expect_error(predict(broken, rows, conn = conn, row_id = "id"), "descriptor|source")
    })

  # CAST(x AS DOUBLE) in prepared cells cannot overwrite raw integer authority.
  data = base
  data[["x"]] = rep(seq_len(4), 30)
  integer_X = function(rows) cbind(1, rows[["x"]], as.numeric(rows[["f"]] == 2))
  data[["y"]] = as.numeric(integer_X(data) %*% c(1, .8, .3)) + .04 * sin(index * .53)
  fit_case(y ~ x + factor(f), data, integer_X, 1:3, "compress", FALSE,
    change = function(rows) { rows[["x"]] = rows[["x"]] + 2L; rows },
    inspect = function(model, rows) {
      layout = model[["linear_design"]]
      expect_identical(layout[["source_classes"]][["x"]], "integer")
      expect_true(all(names(layout[["source_classes"]]) %in% names(layout[["source_map"]])))
      forced = names(layout[["factor_roles"]])[layout[["factor_roles"]] == "forced"]
      expect_identical(model[["factor_classes"]][[forced]], layout[["source_classes"]][[layout[["factor_sources"]][[forced]]]])
    })

  # Native i() classes remain authoritative when that integer source is also FE.
  data = base
  data[["group code"]] = rep(seq_len(4), each = 30)
  integer_fe_X = function(rows) cbind(rows[["x"]], outer(rows[["group code"]], seq_len(4), `==`) * 1)
  data[["y"]] = as.numeric(integer_fe_X(data) %*% c(.9, .3, -.2, .5, .1)) + .03 * sin(index * .79)
  fit_case(y ~ x + i(`group code`, ref = 1) | `group code`, data, integer_fe_X, 1L, "compress", TRUE,
    change = function(rows) { rows[["x"]] = rows[["x"]] + .2; rows },
    inspect = function(model, rows) {
      layout = model[["linear_design"]]
      spec = model[["i_specs"]][[1L]]
      first = spec[["first"]]
      expect_identical(model[["factor_classes"]][[first]], layout[["source_classes"]][[first]])
      expect_identical(spec[["source_type"]], "numeric")
      expect_false(identical(first, "group code"))
      expect_identical(unname(layout[["fe_sources"]]), "group code")
    })

  # A structural categorical g dictionary stays available after FE drops its beta.
  data = base
  data[["g"]] = rep(letters[1:4], each = 30)
  character_fe_X = function(rows) cbind(rows[["x"]], outer(rows[["g"]], letters[1:4], `==`) * 1)
  data[["y"]] = as.numeric(character_fe_X(data) %*% c(.9, .3, -.2, .5, .1)) + .03 * cos(index * .67)
  fit_case(y ~ g + x | g, data, character_fe_X, 1L, "demean", TRUE,
    change = function(rows) { rows[["x"]] = rows[["x"]] + .2; rows },
    inspect = function(model, rows) {
      layout = model[["linear_design"]]
      state = model[["fitted_state"]]
      expect_identical(layout[["factor_roles"]][["g"]], "ordinary")
      expect_true(length(state[["specification"]][["design_ids"]]) == 4L && length(coef(model)) == 1L)
      expect_identical(layout[["source_classes"]][["g"]], "character")
      joint = state[["coefficients"]][["joint"]]
      tail = joint[state[["coefficients"]][["fe_n"]] + seq_len(4)]
      expect_equal(tail[1:3], rep(0, 3), tolerance = 1e-12)
    })
})
