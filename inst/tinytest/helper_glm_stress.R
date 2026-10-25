# A seeded design grid, not repeated draws from one easy fixture. No outcomes
# are forced to keep FE groups: independent fixest pruning defines the oracle.
dbglm_stress_data = function(seed, profile) {
  settings = list(
    rare_correlated = c(n = 6000, rate = .012, dependence = .9, imbalance = 1.5, noise = 1),
    very_rare_weighted = c(n = 10000, rate = .003, dependence = .75, imbalance = 2, noise = 1),
    near_collinear = c(n = 5000, rate = .05, dependence = .95, imbalance = 1.4, noise = .04),
    unbalanced = c(n = 6000, rate = .02, dependence = .6, imbalance = 2.8, noise = 1))
  cfg = settings[[profile]]
  set.seed(seed)
  n = as.integer(cfg[["n"]])
  f = sample.int(12, n, TRUE, prob = seq_len(12)^cfg[["imbalance"]])
  g = ifelse(runif(n) < cfg[["dependence"]], (f - 1L) %% 9L + 1L, sample.int(9, n, TRUE))
  d = data.frame(id = seq_len(n), f = f, g = g, h = sample.int(6, n, TRUE),
    j = sample.int(3, n, TRUE), k = sample(0:1, n, TRUE), l = sample(0:1, n, TRUE),
    cluster = sample.int(23, n, TRUE), x = rnorm(n) + .3 * sin(f),
    a = sample(c(1, 2, 10), n, TRUE, c(.7, .2, .1)), o = runif(n, -.8, .8))
  d[["z"]] = if (profile == "near_collinear") .998 * d[["x"]] + cfg[["noise"]] * rnorm(n) else rnorm(n)
  if (profile == "rare_correlated") d[["a"]] = 1L
  eta = binomial("cloglog")[["linkfun"]](cfg[["rate"]]) + .35 * d[["x"]] - .2 * d[["z"]] +
    .2 * sin(f) + .1 * cos(g) + .05 * d[["h"]] + d[["o"]]
  d[["y"]] = rbinom(n, 1, binomial("cloglog")[["linkinv"]](eta))
  d
}

dbglm_test_stress = function(conn, seeds = c(104L, 611L, 2401L), verbose = FALSE) {
  formula = y ~ x + z | f + g + h + j + k + l
  profiles = c("rare_correlated", "very_rare_weighted", "near_collinear", "unbalanced")
  for (seed in seeds) for (i in seq_along(profiles)) {
    profile = profiles[[i]]
    label = paste(profile, seed)
    if (verbose) message("[GLM stress] ", label)
    d = dbglm_stress_data(seed, profile)
    covariance = list("iid", "hc1", ~cluster)[[(i + match(seed, seeds) - 2L) %% 3L + 1L]]
    before = dbglm_validation_tables(conn)
    table = dbglm_validation_upload(conn, d)
    retained = basename(tempfile("glm_retained_"))
    tryCatch({
      ref = dbglm_reference_feglm(formula, d, family = binomial("cloglog"), weights = ~a, offset = ~o,
        vcov = if (identical(covariance, "hc1")) "hetero" else covariance, fixef.rm = "perfect_fit")
      fit = dbglm(formula, conn = conn, table = table, family = binomial("cloglog"),
        weights = "a", offset = "o", vcov = covariance, row_id = "id", sample_table = retained)
      dbglm_validation_compare(fit, ref, label)
      keys = as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained, "ORDER BY id"))[["id"]])
      expect_identical(keys, sort(as.numeric(d[["id"]][fixest::obs(ref)])), info = label)
      expect_equal(fit[["nobs_dropped"]], nrow(d) - length(keys), info = label)
      expect_identical(dbglm_validation_tables(conn), sort(c(before, table, retained)), info = label)
    }, finally = {
      dbglm_validation_drop(conn, retained)
      dbglm_validation_drop(conn, table)
    })
    expect_identical(dbglm_validation_tables(conn), before, info = label)
  }
}

dbglm_test_invariance = function(conn, verbose = FALSE) {
  d = dbglm_stress_data(611L, "rare_correlated")
  formula = y ~ x + z | f + g + h + j + k + l
  before = dbglm_validation_tables(conn)
  base = dbglm_validation_upload(conn, d)
  retained = basename(tempfile("glm_invariant_keys_"))
  on.exit({dbglm_validation_drop(conn, retained); dbglm_validation_drop(conn, base)})
  fit = dbglm(family = stats::binomial("cloglog"), formula, conn = conn, table = base, weights = "a", offset = "o", vcov = ~cluster,
    row_id = "id", sample_table = retained)
  keys = sort(as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained))[["id"]]))
  variants = list(permutation = d[sample.int(nrow(d)), ], renamed_fe = d, units = d)
  for (f in c("f", "g", "h", "j", "k", "l")) {
    variants[["renamed_fe"]][[f]] = paste0("FE' / ", f, "_", 1000L - 13L * d[[f]])
  }
  variants[["units"]][["x"]] = 1000 * d[["x"]]
  variants[["units"]][["z"]] = d[["z"]] / 100
  for (name in names(variants)) {
    if (verbose) message("[GLM invariance] ", name)
    table = dbglm_validation_upload(conn, variants[[name]])
    sample = basename(tempfile("glm_variant_keys_"))
    tryCatch({
      got = dbglm(family = stats::binomial("cloglog"), formula, conn = conn, table = table, weights = "a", offset = "o", vcov = ~cluster,
        row_id = "id", sample_table = sample)
      scale = if (name == "units") c(1000, .01) else c(1, 1)
      expect_true(max(abs(coef(got) * scale - coef(fit))) < 1e-7, info = name)
      expect_true(max(abs(vcov(got) * outer(scale, scale) - vcov(fit))) < 1e-7, info = name)
      expect_true(abs(as.numeric(logLik(got) - logLik(fit))) / nobs(fit) < 1e-10, info = name)
      got_keys = sort(as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", sample))[["id"]]))
      expect_identical(got_keys, keys, info = name)
      expect_true(got[["score_max"]] <= got[["control"]][["score_tol"]], info = name)
    }, finally = {dbglm_validation_drop(conn, sample); dbglm_validation_drop(conn, table)})
  }
  expect_identical(dbglm_validation_tables(conn), sort(c(before, base, retained)))
}
