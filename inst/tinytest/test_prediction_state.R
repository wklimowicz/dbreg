source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_prediction_state.R", package = "dbreg"), local = TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbglm_test_prediction_state(conn)
  dbglm_test_prediction_keys(conn)
  dbglm_test_prediction_weight_scale(conn)
})

local({
  # A recovery transfer can exceed the retained-matrix budget while fitting
  # safely inside the independent shared transfer budget. Keep both bounds.
  saved = options(dbreg.prediction_max_entries = 1200L, dbreg.ap_moment_max_cells = 8000L)
  on.exit(options(saved))
  set.seed(8211)
  n = 1200L
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample.int(13, n, TRUE),
    h = sample.int(9, n, TRUE), j = sample.int(5, n, TRUE), k = sample.int(3, n, TRUE),
    sex = sample(0:1, n, TRUE), part = sample(0:1, n, TRUE))
  d[["y"]] = .3*d[["x"]] -.2*d[["z"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n, sd = .1)
  formula = y ~ x + z | g + h + j + k + sex + part
  fit = dbreg(formula, data = d, strategy = "demean", retain = "prediction")
  reference = fixest::feols(formula, d, fixef.rm = "none", notes = FALSE)
  input = d[1:30, setdiff(names(d), "y"), drop = FALSE]
  expected = as.numeric(predict(reference, input))
  expect_equal(predict(fit, input), expected, tolerance = 1e-7)
  restored = unserialize(serialize(fit, NULL))
  close(fit)
  expect_equal(predict(restored, input), expected, tolerance = 1e-7)
})
