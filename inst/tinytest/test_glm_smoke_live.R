# One small full-model native GLM for merge smoke; broad cases run separately.
if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbglm_validation_connection(TRUE)
    on.exit(DBI::dbDisconnect(conn))
    before = dbglm_validation_tables(conn)
    set.seed(721)
    n = 360L
    d = data.frame(x = rnorm(n), z = rnorm(n), g = rep(1:6, length.out = n),
      h = sample.int(4, n, TRUE), cl = rep(1:9, length.out = n), w = runif(n, .5, 2), o = runif(n, -.2, .2))
    d[["y"]] = rbinom(n, 1, stats::plogis(.3*d[["x"]]-.2*d[["z"]]+sin(d[["g"]])+d[["o"]]))
    table = dbglm_validation_upload(conn, d)
    on.exit(dbglm_validation_drop(conn, table), add = TRUE, after = FALSE)
    reference = dbglm_reference_feglm(y ~ x + z | g + h, d, family = stats::binomial("logit"),
      weights = ~w, offset = ~o, vcov = ~cl, fixef.rm = "perfect_fit")
    fit = dbglm(y ~ x + z | g + h, conn = conn, table = table, family = stats::binomial("logit"),
      weights = "w", offset = "o", vcov = ~cl)
    dbglm_validation_compare(fit, reference, "merge smoke")
    close(fit)
    expect_identical(dbglm_validation_tables(conn), sort(c(before, table)))
  })
}
