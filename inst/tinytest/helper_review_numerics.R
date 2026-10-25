dbreg_test_review_numerics = function(conn) {
  # Keep ownership in this function's environment without nonlocal assignment.
  owned = new.env(parent = emptyenv())
  owned[["tables"]] = character()
  upload = function(data) {
    table = dbglm_validation_upload(conn, data)
    owned[["tables"]] = c(owned[["tables"]], table)
    table
  }
  on.exit(for (table in rev(owned[["tables"]])) dbglm_validation_drop(conn, table), add = TRUE)
  set.seed(447)
  n = 72L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:3, each = 24), w = 1)
  d[["y"]] = .4*d[["x"]]-.2*d[["z"]]+rnorm(n)
  for (column in c("x", "y")) {
    shifted = d
    shifted[[column]] = shifted[[column]]+if (column == "x") 1e6 else 1e8
    table = upload(shifted)
    actual = dbreg(y ~ x+z, conn = conn, table = table, strategy = "moments")
    ref = lm(y ~ x+z, shifted)
    expect_true(all(abs(coef(actual)-coef(ref)) <= 1e-6+1e-6*abs(coef(ref))))
    expect_true(all(abs(vcov(actual)-vcov(ref)) <= 1e-7+1e-6*abs(vcov(ref))))
    expect_equal(attr(vcov(actual), "rss"), sum(residuals(ref)^2), tolerance = 1e-6)
    robust = dbreg(y ~ x+z, conn = conn, table = table, strategy = "moments", vcov = "hc1")
    baseline = fixest::feols(y ~ x+z, d, vcov = "hc1", notes = FALSE)
    A = diag(3)
    if (column == "x") A[1,2] = -1e6
    expected = A %*% vcov(baseline) %*% t(A)
    expect_true(all(abs(vcov(robust)-expected) <= 1e-7+1e-6*abs(expected)))
  }
  for (scale in c(1e-200, 1e200)) {
    scaled = d
    scaled[["w"]] = scale
    table = upload(scaled)
    for (kind in list("iid", "hc1", ~g)) {
      actual = dbreg(y ~ x+z, conn = conn, table = table, weights = "w", vcov = kind, strategy = "moments")
      ref = fixest::feols(y ~ x+z, d, vcov = kind, notes = FALSE)
      expect_equal(coef(actual), coef(ref), tolerance = 1e-7)
      expect_true(all(abs(vcov(actual)-vcov(ref)) <= 1e-8+1e-6*abs(vcov(ref))))
    }
  }
  table = upload(d)
  integer_source = paste0("(SELECT y,CAST(50000+id AS INT) AS int_x,CAST(60000+id AS INT) AS int_z FROM ", table, ") review_integer")
  integer_fit = dbreg(y ~ int_x:int_z, conn = conn, table = integer_source, strategy = "moments")
  integer_reference = stats::lm(y ~ I((50000+id)*(60000+id)), d)
  expect_equal(unname(coef(integer_fit)), unname(coef(integer_reference)), tolerance = 1e-6)
  expect_equal(as.numeric(vcov(integer_fit)), as.numeric(vcov(integer_reference)), tolerance = 1e-6)
  fit = dbreg(y ~ x+z, conn = conn, table = table, vcov = ~g, retain = "all", row_id = "id", strategy = "moments")
  on.exit(close(fit), add = TRUE)
  recomputed = summary(fit, vcov = ~g)
  expect_equal(as.numeric(fit[["coeftable"]]), as.numeric(recomputed[["coeftable"]]), tolerance = 1e-7)
  expect_equal(as.numeric(confint(fit)), as.numeric(confint(recomputed)), tolerance = 1e-7)
  bins = dbbinsreg(y ~ x, conn = conn, table = table, nbins = 3, plot = FALSE, vcov = "hc1")
  local = dbbinsreg(y ~ x, data = d, nbins = 3, plot = FALSE, vcov = "hc1")
  expect_true(all(is.finite(bins[["points"]][["fit"]])))
  expect_equal(bins[["points"]], local[["points"]], tolerance = 1e-7)
  dbreg_test_review_compression(conn)
}

# Focused native coverage for the sparse compressed-design follow-up.
dbreg_test_review_compression = function(conn) {
  before = dbglm_validation_tables(conn)
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE)
  set.seed(731)
  d = data.frame(id = 1:72, x = rnorm(72), g = rep(1:6, each = 12),
    f = rep(c("red", "blue"), 36), cl = rep(1:3, 24), w = 1e-200)
  d[["y"]] = .4*d[["x"]]+sin(d[["g"]])+(d[["f"]] == "red")+rnorm(72)
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  fit = dbreg(y ~ x+f | g, conn = conn, table = table, weights = "w", vcov = ~cl,
    strategy = "compress", row_id = "id", retain = "all")
  on.exit(close(fit), add = TRUE)
  ref = fixest::feols(y ~ x+f | g, d, vcov = ~cl, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
  expect_equal(coef(fit), coef(ref), tolerance = 1e-7)
  expect_true(all(abs(vcov(fit)-vcov(ref)) <= 1e-8+1e-6*abs(vcov(ref))))
  updated = summary(fit, vcov = ~cl)
  expect_equal(as.numeric(updated[["coeftable"]]), as.numeric(fit[["coeftable"]][rownames(updated[["coeftable"]]), , drop = FALSE]), tolerance = 1e-7)
  expect_equal(predict(fit, d), predict(ref, d), tolerance = 1e-7)
  close(fit)
  dbreg:::release_fit_resources(resources)
  expect_identical(dbglm_validation_tables(conn), before)
  expect_true(DBI::dbIsValid(conn))
}
