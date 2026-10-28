# Independent full-recovery-coordinate regression fixture.
# All statistical expectations use explicit dummy/trend X, never another dbreg path.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  originals = lapply(c("dbreg", "finalize_dbreg_result", "finalize_linear_prediction",
    "finalize_fwl_prediction", "build_prediction_state"), function(name) get(name, asNamespace("dbreg")))
  names(originals) = c("dbreg", "finalize_dbreg_result", "finalize_linear_prediction",
    "finalize_fwl_prediction", "build_prediction_state")
  original_environments = lapply(originals, environment)
  # Private function copies count actual FWL reuse without namespace mutation or
  # a production hook. Numerical functions/arguments remain the actual owners.
  fwl_runner = function(counter) {
    scope = new.env(parent = asNamespace("dbreg"))
    fwl = dbreg:::finalize_fwl_prediction
    fresh = dbreg:::build_prediction_state
    scope[["finalize_fwl_prediction"]] = function(...) {
      out = fwl(...)
      counter[["reused"]] = counter[["reused"]] + as.integer(!is.null(out))
      out
    }
    scope[["build_prediction_state"]] = function(...) {
      counter[["fresh"]] = counter[["fresh"]] + 1L
      fresh(...)
    }
    for (name in c("finalize_linear_prediction", "finalize_dbreg_result", "dbreg")) {
      fn = get(name, envir = asNamespace("dbreg"))
      environment(fn) = scope
      scope[[name]] = fn
    }
    scope[["dbreg"]]
  }
  check = function(collision = FALSE, storage = "portable", dropped = FALSE,
                   twoway = FALSE, slopes = TRUE, raw = FALSE, trend = FALSE,
                   strategy = "demean", retain = "prediction", forced = FALSE, typed = FALSE) {
    old = options(dbreg.prediction_max_entries = if (storage == "database") 64 else 100000,
      dbreg.ap_solver = "moments")
    on.exit(options(old))
    id = seq_len(240)
    groups = if (twoway) expand.grid(g = seq_len(8), h = seq_len(6), r = seq_len(5)) else
      data.frame(g = rep(seq_len(40), each = 6), h = 1L, r = id)
    d = data.frame(id = id, g = groups[["g"]], h = groups[["h"]],
      f = if (twoway) ifelse((groups[["g"]] + groups[["h"]] + groups[["r"]]) %% 2, "a", "b") else rep(c("a", "b"), 120),
      x = sin(id * .17) + cos(id * .39), z = sin(id * .43) + .1 * cos(id * .61),
      z0 = 0, w = 1 + (id %% 5) / 4)
    d[["fb"]] = d[["x"]]
    d[["y"]] = .6 * (d[["f"]] == "b") + 1.7 * d[["x"]] + 1 + sin(d[["g"]] * .23) +
      (if (twoway) .2 * cos(d[["h"]] * .37) else 0)
    if (trend) d[["y"]] = d[["y"]] + .3 * cos(d[["g"]] * .11) * d[["z"]]
    d[["y"]] = d[["y"]] + .05 * sin(id * .71) + .03 * cos(id * .31)
    f = "f"
    x = if (collision) "fb" else "x"
    g = "g"
    if (raw) {
      names(d)[match(c("f", "x", "g"), names(d))] = c("f field", "raw x", "group id")
      f = "f field"
      x = "raw x"
      g = "group id"
    }
    G = if (twoway) 8L else 40L
    full_X = function(data) {
      D = outer(data[[g]], seq_len(G), `==`) * 1
      nuisance = cbind(D, if (twoway) outer(data[["h"]], 2:6, `==`) * 1,
        if (trend) D * data[["z"]])
      if (slopes) cbind(as.numeric(data[[f]] == "b"), data[[x]], nuisance) else nuisance
    }
    X = full_X(d)
    reference = stats::lm.wfit(X, d[["y"]], d[["w"]])
    expect_identical(reference[["rank"]], ncol(X))
    beta = reference[["coefficients"]]
    V = solve(crossprod(X, d[["w"]] * X)) *
      sum(d[["w"]] * reference[["residuals"]]^2) / (nrow(X) - reference[["rank"]])
    quoted = function(x) paste(deparse(as.name(x), backtick = TRUE), collapse = "")
    first = if (forced) paste0("factor(", quoted(f), ")") else if (typed) paste0("i(", quoted(f), ", ref = 'a')") else quoted(f)
    rhs = if (slopes) paste(first, "+", quoted(x), if (dropped) "+ z0" else "") else "1"
    fe = paste0(quoted(g), if (trend) "[z]" else if (twoway) " + h" else "")
    formula = stats::as.formula(paste("y ~", rhs, "|", fe))
    counter = new.env(parent = emptyenv())
    counter[["reused"]] = counter[["fresh"]] = 0L
    run = if (twoway) fwl_runner(counter) else dbreg::dbreg
    model = run(formula, data = d, conn = conn, weights = "w", strategy = strategy,
      retain = retain, vcov = "iid", ssc = list(K.fixef = "full"))
    on.exit(close(model), add = TRUE, after = FALSE)
    state = model[["fitted_state"]]
    if (twoway) {
      expect_true(all(vapply(names(originals), function(name)
        identical(get(name, asNamespace("dbreg")), originals[[name]]), logical(1))))
      expect_true(all(vapply(names(originals), function(name)
        identical(environment(get(name, asNamespace("dbreg"))), original_environments[[name]]), logical(1))))
    }
    expect_identical(identical(state[["prediction"]][["storage"]], "database"), storage == "database")
    if (twoway && storage == "portable") {
      expect_equal(counter[["reused"]], 1L)
      expect_equal(counter[["fresh"]], 0L)
    }
    selected = if (slopes) 1:2 else integer()
    expect_equal(as.numeric(coef(model)), as.numeric(beta[selected]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(model)), as.numeric(V[selected, selected, drop = FALSE]), tolerance = 1e-8)
    rows = d[c(17, 2, 201, 55, 120, 9, 37, 89, 160, 230), ]
    if (slopes) rows[[x]] = rows[[x]] + .15 * sin(seq_len(nrow(rows)))
    rows[["y"]] = -999
    rows[["fb_1"]] = 999
    rows[["dbreg_internal_x_0001__"]] = -999
    expected = as.numeric(full_X(rows) %*% beta)
    if (dropped) {
      rows[["z0"]][4L] = 1
      expected[4L] = NA_real_
    }
    # Complete NA masks plus every value; duplicate labels are not uniquified.
    r = predict(model, rows)
    expect_identical(is.na(r), is.na(expected))
    expect_equal(as.numeric(r), expected, tolerance = 1e-8)
    out = predict(model, rows, conn = conn, row_id = "id")
    table = as.data.frame(out)
    sql = table[[out[["prediction_column"]]]][match(rows[["id"]], table[["id"]])]
    close(out)
    expect_identical(is.na(sql), is.na(expected))
    expect_equal(as.numeric(sql), expected, tolerance = 1e-8)
    # FULL recovery IDs link to existing recipes.
    ids = state[["specification"]][["design_ids"]]
    columns = model[["linear_design"]][["columns"]]
    positions = match(ids, columns[["id"]])
    expect_true(is.character(ids) && !anyNA(positions) && !anyDuplicated(ids))
    expect_true(all(columns[["role"]][positions] == "slope"))
    expect_identical(length(ids), if (!slopes) 0L else if (dropped) 3L else 2L)
    expect_identical(length(state[["specification"]][["design_center"]]), length(ids))
    tail = if (storage == "database") state[["coefficients"]][["slopes"]] else
      state[["coefficients"]][["joint"]][state[["coefficients"]][["fe_n"]] + seq_along(ids)]
    expect_identical(length(tail), length(ids))
    if (slopes && !dropped) {
      bad = rows
      bad[[f]][1L] = "unseen"
      if (forced) bad[[paste0("factor(", f, ")")]] = "a"
      expect_error(predict(model, bad), "unseen|level|dictionary")
      bad = rows
      bad[[x]] = as.character(bad[[x]])
      expect_error(predict(model, bad), "numeric|incompatible|design")
    }
    if (slopes && !dropped && !twoway && !trend && !raw && !forced && !typed && strategy == "demean") {
      corrupted = model
      corrupted[["fitted_state"]][["specification"]][["design_ids"]][2L] = ids[[1L]]
      expect_error(predict(corrupted, rows), "coordinates")
      expect_error(predict(corrupted, rows, conn = conn, row_id = "id"), "coordinates")
      corrupted[["fitted_state"]][["specification"]][["design_ids"]][2L] = "unknown-coordinate"
      expect_error(predict(corrupted, rows), "coordinates")
      expect_error(predict(corrupted, rows, conn = conn, row_id = "id"), "coordinates|design")
    }
    if (collision) {
      expect_true(anyDuplicated(names(coef(model))) > 0L)
      expect_true(anyDuplicated(state[["specification"]][["design_names"]]) > 0L)
    }
    if (dropped) {
      expect_equal(length(coef(model)), 2L)
      expect_equal(tail[[3L]], 0, tolerance = 1e-12)
    }
    if (retain == "all") expect_equal(as.numeric(predict(model)), as.numeric(X %*% beta), tolerance = 1e-8)
    if (storage == "portable" && retain == "prediction" && !dropped) {
      saved = tempfile(fileext = ".rds")
      on.exit(unlink(saved), add = TRUE)
      saveRDS(model, saved)
      restored = readRDS(saved)
      expect_equal(predict(restored, rows), expected, tolerance = 1e-8)
    }
    if (storage == "database") {
      close(model)
      expect_error(predict(model, rows), "closed|released")
    }
  }
  check() # Unique-label control.
  check(collision = TRUE)
  check(collision = TRUE, storage = "database")
  check(dropped = TRUE)
  check(dropped = TRUE, storage = "database")
  check(twoway = TRUE, collision = TRUE) # Actual FWL reuse return is observed.
  check(twoway = TRUE, collision = TRUE, storage = "database")
  check(collision = TRUE, strategy = "compress")
  check(slopes = FALSE) # Valid empty full-slope axis.
  check(raw = TRUE)
  check(forced = TRUE)
  check(typed = TRUE)
  check(trend = TRUE)
  check(collision = TRUE, retain = "all")
})
