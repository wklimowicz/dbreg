dbreg_test_inference_numerics = function(conn) {
  before = dbglm_validation_tables(conn)
  resources = dbreg:::dbglm_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::dbglm_release_resources(resources), add = TRUE)
  set.seed(399)
  n = 240L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:12, each = 20),
    h = sample(1:6, n, TRUE), w = 1)
  d[["y"]] = .4*d[["x"]] -.2*d[["z"]] + sin(d[["g"]]) + rnorm(n)
  reference = fixest::feols(y ~ x + z, d, notes = FALSE)
  for (scale in c(1, 1e-200, 1e200)) {
    message("[Inference numerics] OLS uniform weight scale ", scale)
    d[["w"]] = scale
    table = dbreg:::prediction_upload_frame(conn, d, resources)
    fit = dbreg(y ~ x + z, conn = conn, table = table, weights = "w", row_id = "id",
      strategy = "moments", retain = "all")
    expect_equal(coef(fit), coef(reference), tolerance = 1e-7)
    # Uniform scaling leaves OLS coefficients/sandwich covariance unchanged.
    # fixest's default collinearity gate rejects the uniformly tiny fit itself;
    # the independent unit-weight fit is the scale-invariant numerical oracle.
    for (request in list("iid", "hc1", ~g, ~g+h)) {
      V = vcov(fit, vcov = request, ssc = fixest::ssc())
      expected = suppressWarnings(vcov(reference, vcov = request))
      expect_true(all(is.finite(V)))
      expect_true(all(abs(V-expected) <= 1e-9 + 1e-6*abs(expected)))
    }
    close(fit)
  }
  d[["w"]] = rep(c(.2,.7,1), length.out = n)
  d[["b"]] = stats::rbinom(n, 1L, stats::plogis(.2*d[["x"]]))
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  fit = dbglm(b ~ x + z | g + h, conn = conn, table = table, weights = "w",
    family = stats::binomial("logit"), fixef_rm = "none", retain = "all")
  reference = fixest::feglm(b ~ x + z | g + h, d, weights = ~w, family = stats::binomial("logit"),
    fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
  for (request in list("hc1", ~g)) {
    V = vcov(fit, vcov = request)
    expected = vcov(reference, vcov = request)
    expect_true(all(abs(V-expected) <= 1e-8 + 1e-5*abs(expected)))
  }
  close(fit)

  message("[Inference numerics] SQL NaN identities versus literal text NaN")
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  create = function(query) {
    name = basename(tempfile("inference_numeric_keys_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, name, query, c(names(d), "floating_key", "text_key")))
    q(name)
  }
  bad = create(paste0("SELECT *,CASE WHEN id=1 THEN CAST('NaN' AS DOUBLE) ELSE CAST(id AS DOUBLE) END AS floating_key,",
    "'NaN' AS text_key FROM ", table))
  inventory = dbglm_validation_tables(conn)
  expect_error(dbreg(y ~ x, conn = conn, table = bad, row_id = "floating_key", strategy = "moments", retain = "all"), "unique and nonmissing")
  expect_identical(dbglm_validation_tables(conn), inventory)
  valid = create(paste0("SELECT *,CAST(id AS DOUBLE) AS floating_key,'NaN' AS text_key FROM ", table))
  fit = dbreg(y ~ x, conn = conn, table = valid, row_id = c("floating_key", "text_key"), strategy = "moments", retain = "all")
  output = dbreg::obs(fit)
  actual = as.data.frame(output)
  expect_equal(sort(actual[["floating_key"]]), as.numeric(seq_len(n)))
  expect_true(all(actual[["text_key"]] == "NaN"))
  close(output)
  close(fit)
  dbreg:::dbglm_release_resources(resources)
  resources[["tables"]] = character()
  expect_identical(dbglm_validation_tables(conn), before)
  expect_true(DBI::dbIsValid(conn))
}
