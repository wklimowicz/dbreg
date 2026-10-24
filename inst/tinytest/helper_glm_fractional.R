# The same value comparisons run locally and on the native warehouse.
dbglm_fractional_phase = function(label, expr, verbose = FALSE) {
  if (verbose) message("[Fractional weights] starting ", label)
  started = proc.time()[["elapsed"]]
  on.exit(if (verbose) message("[Fractional weights] ", label, ": ",
    round(proc.time()[["elapsed"]] - started, 2), "s"))
  force(expr)
}

dbglm_test_fractional_weights = function(conn, verbose = FALSE) {
  set.seed(240719)
  n = 720L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n),
    g = rep(seq_len(12), each = 60), h = rep(seq_len(9), 80),
    cl = rep(seq_len(24), each = 30), o = runif(n, -.2, .2))
  d[["a"]] = rep(c(0, .1, .5, .9, 1.5, 2.5, 3.5, 1e-8), length.out = n)
  d[["a"]][2] = NA_real_
  d[["x"]][3] = NA_real_
  d[["b"]] = rep(c(.2, .7, 0), length.out = n)
  # Generate outcomes before inserting missingness into the design predictor.
  eta = -.4 + .25 * replace(d[["x"]], 3, 0) - .2 * d[["z"]] +
    .15 * sin(d[["g"]]) + .1 * cos(d[["h"]]) + d[["o"]]
  families = list(binomial("cloglog"), binomial("logit"), binomial("probit"), poisson())
  before = dbglm_validation_tables(conn)
  close = function(actual, expected, atol = 1e-7, rtol = 1e-6, label) {
    expect_identical(dim(actual), dim(expected), info = label)
    expect_true(all(is.finite(actual)) &&
      all(abs(actual - expected) <= atol + rtol * abs(expected)), info = label)
  }
  for (family in families) {
    d[["y"]] = if (family[["family"]] == "poisson") rpois(n, exp(eta)) else
      rbinom(n, 1, family[["linkinv"]](eta))
    table = dbglm_fractional_upload(conn, d)
    retained = basename(tempfile("fractional_retained_"))
    tryCatch({
      for (fml in list(y ~ x + z, y ~ x + z | g + h)) {
        for (covariance in list("iid", "hc1", ~cl)) {
          label = paste(family[["family"]], family[["link"]], deparse(fml), deparse(covariance))
          reference = dbglm_fractional_phase(paste(label, "reference"), suppressWarnings(dbglm_reference_feglm(fml, d, family = family,
            weights = ~a, offset = ~o,
            vcov = if (identical(covariance, "hc1")) "hetero" else covariance)), verbose)
          fit = dbglm_fractional_phase(paste(label, "dbglm"), dbglm(fml, conn = conn, table = table, family = family,
            weights = "a", offset = "o", vcov = covariance,
            row_id = "id", sample_table = retained), verbose)
          expect_identical(names(coef(fit)), names(coef(reference)), info = label)
          close(coef(fit), coef(reference), label = label)
          close(vcov(fit), vcov(reference), atol = 1e-8, rtol = 1e-5, label = label)
          close(sqrt(diag(vcov(fit))), sqrt(diag(vcov(reference))), atol = 1e-8, label = label)
          rows = fixest::obs(reference)
          keys = as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained, "ORDER BY id"))[["id"]])
          expect_identical(keys, as.numeric(d[["id"]][rows]), info = label)
          expect_equal(nobs(fit), length(rows), info = label)
          expect_equal(fit[["n_params_model"]], reference[["nparams"]], info = label)
          expect_equal(fit[["n_params_vcov"]], attr(vcov(reference, attr = TRUE), "df.K"), info = label)
          close(as.numeric(logLik(fit)), as.numeric(logLik(reference)), label = label)
          close(fit[["deviance"]], reference[["deviance"]], label = label)
          mu = predict(reference, newdata = d[rows, ], type = "response")
          density = if (family[["family"]] == "poisson") dpois(d[["y"]][rows], mu, log = TRUE) else
            dbinom(d[["y"]][rows], 1, mu, log = TRUE)
          close(as.numeric(logLik(fit, type = "weighted")), sum(d[["a"]][rows] * density), label = label)
          expect_equal(attr(logLik(fit, type = "weighted"), "df"), reference[["nparams"]])
          expect_true(fit[["converged"]] && fit[["score_max"]] <= fit[["control"]][["score_tol"]], info = label)
          if (!length(fit[["fe"]])) close(predict(fit, newdata = d[rows, ], type = "response"), mu, label = label)
          dbglm_validation_drop(conn, retained)
        }
      }
      # Repeat fits on one source with different posterior-weight columns.
      first = dbglm_fractional_phase("repeated weights: first fit",
        dbglm(y ~ x + z | g + h, conn = conn, table = table, family = family, weights = "a"), verbose)
      saved = serialize(first, NULL)
      second = dbglm_fractional_phase("repeated weights: changed column",
        dbglm(y ~ x + z | g + h, conn = conn, table = table, family = family, weights = "b"), verbose)
      second_ref = suppressWarnings(dbglm_reference_feglm(y ~ x + z | g + h, d, family = family, weights = ~b))
      close(coef(second), coef(second_ref), label = "changed weights")
      close(vcov(second), vcov(second_ref), label = "changed weights")
      again = dbglm_fractional_phase("repeated weights: restored column",
        dbglm(y ~ x + z | g + h, conn = conn, table = table, family = family, weights = "a"), verbose)
      close(coef(again), coef(first), label = "repeated weights")
      expect_identical(serialize(first, NULL), saved)
      # Force the SQL projector too: fractional weights must not rely on the
      # bounded moments solver retaining all FE summaries in R.
      sql_fit = dbglm_fractional_phase("forced SQL projection",
        dbglm(y ~ x + z | g + h, conn = conn, table = table, family = family,
          weights = "a", control = dbglm_control(projection_method = "sql"), verbose = verbose), verbose)
      close(coef(sql_fit), coef(first), label = "SQL projection")
      close(vcov(sql_fit), vcov(first), label = "SQL projection")
      close(as.numeric(logLik(sql_fit)), as.numeric(logLik(first)), label = "SQL projection")
      expect_error(dbglm(y ~ x + z | g + h, conn = conn, table = table, family = family,
        weights = "a", control = dbglm_control(glm_max_iter = 1), sample_table = retained), "converge")
      expect_identical(dbglm_validation_tables(conn), sort(c(before, table)))
    }, finally = {
      dbglm_validation_drop(conn, retained)
      dbglm_validation_drop(conn, table)
    })
  }
  # R ties-to-even at half-integer weights, including the adjacent doubles,
  # checked against the independent binomial AIC implementation on both engines.
  w = c(.5 - 2^-54, .5, .5 + 2^-53, 1.5, 2.5, 3.5, .1, .9, 1e-12)
  kernel = data.frame("__dbglm_a" = w, "__dbglm_y" = rep(c(0, 1), length.out = length(w)),
    mu = seq(.1, .9, length.out = length(w)), check.names = FALSE)
  table = dbglm_validation_upload(conn, kernel)
  tryCatch({
    expected = -binomial()[["aic"]](kernel[["__dbglm_y"]], rep(1, length(w)), kernel[["mu"]], w, 0) / 2
    expect_equal(dbreg:::dbglm_binomial_loglik(conn, table), expected, tolerance = 1e-13)
  }, finally = dbglm_validation_drop(conn, table))
  expect_identical(dbglm_validation_tables(conn), before)
}

# The bounded VALUES uploader is also used for coefficient dictionaries and
# deliberately rejects NA. Inject the fixture's missing values in SQL on both
# engines, without weakening that production transport check.
dbglm_fractional_upload = function(conn, data) {
  uploaded = data
  uploaded[["a"]][2] = 0
  uploaded[["x"]][3] = 0
  source = dbglm_validation_upload(conn, uploaded)
  on.exit(dbglm_validation_drop(conn, source))
  table = basename(tempfile("fractional_input_"))
  fields = names(data)
  fields[fields == "a"] = "CASE WHEN id=2 THEN CAST(NULL AS DOUBLE) ELSE a END AS a"
  fields[fields == "x"] = "CASE WHEN id=3 THEN CAST(NULL AS DOUBLE) ELSE x END AS x"
  DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS SELECT", paste(fields, collapse = ","), "FROM", source))
  table
}

# Independent saturated two-cell MLE and Fisher information. Tiny uniform
# weights used to satisfy absolute convergence checks after only one step.
dbglm_test_fractional_scale = function(conn, verbose = FALSE) {
  d = data.frame(x = rep(rep(c(-1, 1), each = 100), 3), g = rep(1:3, each = 200),
    y = rep(c(rep(c(0, 0, 0, 1), 25), rep(c(1, 1, 1, 1, 0), 20)), 3),
    a = 1, tiny = 1e-20, cl = rep(1:20, each = 30))
  table = dbglm_validation_upload(conn, d)
  on.exit(dbglm_validation_drop(conn, table))
  X = model.matrix(~x, d)
  mu = ifelse(d[["x"]] == -1, .25, .8)
  for (family in list(binomial("cloglog"), binomial("logit"), binomial("probit"), poisson())) {
    eta = family[["linkfun"]](c(.25, .8))
    expected = c(mean(eta), diff(eta) / 2)
    information_w = family[["mu.eta"]](family[["linkfun"]](mu))^2 / family[["variance"]](mu)
    V = solve(crossprod(X, information_w * X))
    density = if (family[["family"]] == "poisson") dpois(d[["y"]], mu, log = TRUE) else
      dbinom(d[["y"]], 1, mu, log = TRUE)
    for (fml in list(y ~ x, y ~ x | g)) {
      label = paste(family[["family"]], family[["link"]], deparse(fml))
      plain = dbglm_fractional_phase(paste(label, "unit scale"),
        dbglm(fml, conn = conn, table = table, family = family, weights = "a", ssc = list(K.adj = FALSE)), verbose)
      tiny = dbglm_fractional_phase(paste(label, "tiny scale"),
        dbglm(fml, conn = conn, table = table, family = family, weights = "tiny", ssc = list(K.adj = FALSE)), verbose)
      expected_beta = if (length(tiny[["fe"]])) expected[2] else expected
      expected_V = if (length(tiny[["fe"]])) V[2, 2, drop = FALSE] else V
      expect_true(max(abs(unname(coef(tiny)) - expected_beta)) < 1e-7)
      expect_true(max(abs(unname(vcov(tiny)) * 1e-20 - expected_V)) < 1e-8)
      expect_equal(coef(tiny), coef(plain), tolerance = 1e-7)
      expect_true(abs(as.numeric(logLik(tiny, type = "weighted")) / 1e-20 - sum(density)) < 1e-6)
      expect_equal(tiny[["deviance"]] / 1e-20, plain[["deviance"]], tolerance = 1e-8)
      expect_equal(nobs(tiny), nrow(d))
      expect_true(tiny[["weight_scale"]] > 0 && tiny[["weight_scale"]] < 1e-20)
      if (family[["family"]] == "binomial") expect_equal(as.numeric(logLik(tiny)), 0)
    }
  }
  # Sandwich covariance is invariant to multiplying every original weight;
  # the IID inverse-information adjustment above must not be applied to it.
  for (fml in list(y ~ x, y ~ x | g)) {
    for (kind in list("hc1", ~cl)) {
      label = paste("sandwich scale", deparse(fml), deparse(kind))
      plain = dbglm_fractional_phase(paste(label, "unit"),
        dbglm(fml, conn = conn, table = table, weights = "a", vcov = kind), verbose)
      tiny = dbglm_fractional_phase(paste(label, "tiny"),
        dbglm(fml, conn = conn, table = table, weights = "tiny", vcov = kind), verbose)
      expect_equal(vcov(tiny), vcov(plain), tolerance = 1e-7)
    }
  }
}

# Replace values in the SAME source relation/column, including support changes.
# The six extra rows form a mixed-outcome FE cycle. Zeroing its final row
# triggers recursive pruning; zeroing the last core replicates also makes z=x.
dbglm_test_fractional_mutation = function(conn, verbose = FALSE) {
  before = dbglm_validation_tables(conn)
  d = expand.grid(g = 1:4, h = 1:3, x = c(-1, 1), replicate = 1:20)
  d[["z"]] = ifelse(d[["replicate"]] > 16, -d[["x"]], d[["x"]])
  d[["y"]] = as.numeric(ifelse(d[["x"]] > 0, d[["replicate"]] %% 5 != 0, d[["replicate"]] %% 4 == 0))
  d = rbind(d, data.frame(g = c(5, 5, 6, 6, 7, 7), h = c(1, 4, 4, 5, 5, 1),
    x = 0, replicate = 0, z = 0, y = c(0, 1, 0, 1, 0, 1)))
  d[["id"]] = seq_len(nrow(d))
  d[["a"]] = .2 + .1 * (d[["id"]] %% 4)
  d[["changed"]] = as.numeric(d[["replicate"]] > 16 | d[["id"]] == nrow(d))
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("fractional_mutable_"))
  retained = basename(tempfile("fractional_mutable_retained_"))
  on.exit({
    dbglm_validation_drop(conn, retained)
    dbglm_validation_drop(conn, table)
    dbglm_validation_drop(conn, source)
  })
  replace_weights = function(changed) {
    fields = as.character(DBI::dbQuoteIdentifier(conn, setdiff(names(d), "a")))
    weight = if (changed) "CASE WHEN changed=1 THEN 0.0 ELSE a END AS a" else "a"
    DBI::dbExecute(conn, paste("CREATE OR REPLACE TEMP TABLE", table,
      "AS SELECT", paste(c(fields, weight), collapse = ","), "FROM", source))
  }
  for (family in list(binomial("cloglog"), poisson())) {
    original = saved = NULL
    for (changed in c(FALSE, TRUE, FALSE)) {
      replace_weights(changed)
      current = d
      if (changed) current[["a"]][current[["changed"]] == 1] = 0
      label = paste("source mutation", family[["family"]], if (changed) "pruned" else "full")
      reference = suppressWarnings(dbglm_reference_feglm(y ~ x + z | g + h, current,
        family = family, weights = ~a))
      fit = dbglm_fractional_phase(label, suppressWarnings(dbglm(y ~ x + z | g + h,
        conn = conn, table = table, family = family, weights = "a", row_id = "id",
        sample_table = retained)), verbose)
      dbglm_validation_compare(fit, reference, info = label)
      keys = as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained, "ORDER BY id"))[["id"]])
      expect_identical(keys, as.numeric(current[["id"]][fixest::obs(reference)]), info = label)
      expect_identical(names(coef(fit)), if (changed) "x" else c("x", "z"), info = label)
      expect_equal(unname(fit[["n_fe_levels"]]), if (changed) c(4, 3) else c(7, 5), info = label)
      if (changed) {
        expect_true(fit[["pruning_rounds"]] > 1, info = label)
        expect_equal(fit[["nobs_pruned"]], 5, info = label)
      }
      if (is.null(original)) {
        original = fit
        saved = serialize(fit, NULL)
      } else {
        expect_identical(serialize(original, NULL), saved, info = label)
        if (!changed) {
          expect_equal(coef(fit), coef(original), tolerance = 1e-7, info = label)
          expect_equal(vcov(fit), vcov(original), tolerance = 1e-7, info = label)
        }
      }
      dbglm_validation_drop(conn, retained)
      expect_identical(dbglm_validation_tables(conn), sort(c(before, source, table)))
    }
  }
}
