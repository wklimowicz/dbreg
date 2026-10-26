local({
  set.seed(952)
  d = data.frame(x = rnorm(240), z = rnorm(240), g = rep(1:12, each = 20), w = runif(240, 0.2, 1.5))
  d[["y"]] = 0.5 * d[["x"]] + sin(d[["g"]]) + cos(d[["g"]]) * d[["z"]] + rnorm(240)
  policies = list("full", "nested", "none", list(K.adj = FALSE, G.adj = FALSE),
    list(K.fixef = "none", t.df = 7), list(K.fixef = "full", t.df = "conventional"))
  for (fml in list(y ~ x | g, y ~ x | g[z])) for (strategy in c("demean", "compress")) {
    for (kind in list("iid", "hc1", ~g)) for (policy in policies) {
      info = paste(paste(deparse(fml), collapse=""), strategy, paste(deparse(kind),collapse=""), paste(deparse(policy),collapse=""))
      native = dbreg:::dbglm_ssc(policy)
      reference_policy = do.call(fixest::ssc, native)
      fit = dbreg::dbreg(fml, data = d, weights = ~w, strategy = strategy, vcov = kind, ssc = policy)
      reference = fixest::feols(fml, data = d, weights = ~w, vcov = if (identical(kind, "hc1")) "hetero" else kind,
        ssc = reference_policy, fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7, info = info)
      expect_equal(fit[["coeftable"]][names(stats::coef(fit)), "p.value"], as.numeric(fixest::coeftable(reference)[, 4L]), check.attributes = FALSE, tolerance = 1e-7, info = info)
      expect_equal(attr(stats::vcov(fit), "df.K"), fixest::degrees_freedom(reference, "k"), info = info)
      if (identical(kind, "iid")) {
        retained = dbreg::dbreg(fml, data = d, weights = ~w, strategy = strategy, retain = "all")
        expect_equal(as.numeric(stats::vcov(retained, vcov = kind, ssc = policy)), as.numeric(stats::vcov(fit)), tolerance = 1e-7, info = info)
        close(retained)
      }
    }
  }
})

local({
  set.seed(961)
  d = data.frame(x = rnorm(240), z = rnorm(240), g = rep(1:12, each = 20), w = runif(240, 0.2, 1.5))
  first = seq(1, 240, by = 20)
  last = first + 19
  for (family in list(stats::gaussian(), stats::binomial("logit"), stats::poisson())) {
    eta = 0.2 * d[["x"]] + sin(d[["g"]]) / 5
    d[["y"]] = switch(family[["family"]], gaussian = eta + rnorm(240),
      binomial = stats::rbinom(240, 1, stats::plogis(eta)), poisson = stats::rpois(240, exp(eta)))
    if (family[["family"]] == "binomial") { d[["y"]][first] = 0; d[["y"]][last] = 1 }
    if (family[["family"]] == "poisson") d[["y"]][first] = 1
    for (kind in list("iid", "hc1", ~g)) for (policy in list("full", "nested", "none", list(K.adj = FALSE, G.adj = FALSE, t.df = 7))) {
      normalized = dbreg:::dbglm_ssc(policy)
      fit = dbreg::dbglm(y ~ x | g[z], data = d, family = family, weights = ~w, vcov = kind,
        ssc = policy, fixef_rm = "none")
      reference = fixest::feglm(y ~ x | g[z], data = d, family = family, weights = ~w,
        vcov = if (identical(kind, "hc1")) "hetero" else kind, ssc = do.call(fixest::ssc, normalized),
        fixef.rm = "none", fixef.tol = 1e-11, glm.tol = 1e-12, notes = FALSE)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-6)
      expect_equal(as.numeric(fit[["coeftable"]][, "p.value"]), as.numeric(fixest::coeftable(reference)[, 4L]), tolerance = 1e-6)
    }
  }
})

local({
  set.seed(964)
  d = data.frame(x = rnorm(240), z = rnorm(240), g = rep(1:12, each = 20), w = runif(240, 0.2, 1.5))
  d[["d"]] = d[["z"]] + rnorm(240)
  d[["y"]] = 0.5 * d[["x"]] + 0.7 * d[["d"]] + sin(d[["g"]]) + rnorm(240)
  for (kind in list("iid", "hc1", ~g)) for (policy in list("full", "none", list(K.adj = FALSE, G.adj = FALSE, t.df = 9))) {
    normalized = dbreg:::dbglm_ssc(policy)
    request = if (identical(kind, "hc1")) "hetero" else kind
    for (strategy in c("moments", "compress")) {
      fit = dbreg::dbreg(y ~ x, data = d, weights = ~w, strategy = strategy, vcov = kind, ssc = policy)
      reference = fixest::feols(y ~ x, data = d, weights = ~w, vcov = request,
        ssc = do.call(fixest::ssc, normalized), notes = FALSE)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
    }
    for (fml in list(y ~ x | d ~ z, y ~ x | g | d ~ z)) {
      fit = dbreg::dbivreg(fml, data = d, weights = ~w, vcov = kind, ssc = policy)
      reference = fixest::feols(fml, data = d, weights = ~w, vcov = request,
        ssc = do.call(fixest::ssc, normalized), fixef.rm = "none", notes = FALSE)
      V = stats::vcov(reference)
      dimnames(V) = lapply(dimnames(V), function(names) sub("^fit_", "", names))
      names = names(stats::coef(fit))
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(V[names, names]), tolerance = 1e-7)
    }
  }
})
