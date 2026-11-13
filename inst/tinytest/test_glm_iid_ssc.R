# Independent full-dummy information, not another call to the fitted projection.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  set.seed(43118)
  d = expand.grid(g = 1:4, h = 1:3, x = -1:1, trend = -1:1, rep = 1:3)
  d$w = rep(c(.25, 1.5, 4), length.out = nrow(d))
  d$o = d$h / 20
  d$copy = 2 * d$x
  d$absorbed = d$g
  eta = .1 + .25 * d$x - .1 * d$x^2 + d$g / 20 + d$o
  forms = list(
    additive = list(y ~ x + I(x^2) + copy + absorbed | g + h,
      y ~ x + I(x^2) + factor(g) + factor(h), 0L),
    tuple = list(y ~ x + I(x^2) + copy | g^h,
      y ~ x + I(x^2) + interaction(g, h), 0L),
    varying = list(y ~ x + I(x^2) + copy | g[trend] + h,
      y ~ x + I(x^2) + factor(g) / trend + factor(h), 4L)
  )
  bare = function(x) matrix(as.numeric(x), nrow(x), ncol(x))
  for (family in list(gaussian(), poisson(), binomial("logit"),
      binomial("probit"), binomial("cloglog"))) {
    mu = family$linkinv(eta)
    d$y = switch(family$family, gaussian = rnorm(nrow(d), mu, .6),
      poisson = rpois(nrow(d), mu), binomial = rbinom(nrow(d), 1, mu))
    for (name in names(forms)) {
      spec = forms[[name]]
      dense = suppressWarnings(glm(spec[[2L]], data = d, weights = w, offset = o,
        family = family, control = glm.control(epsilon = 1e-14, maxit = 1000L)))
      X = model.matrix(dense)
      N = nrow(X)
      rank = qr(X)$rank
      information = crossprod(X, X * (d$w * family$mu.eta(predict(dense))^2 /
        family$variance(fitted(dense))))
      bread = solve(information)
      positions = match(c("x", "I(x^2)"), colnames(X))
      # Gaussian dispersion estimates a full-model variance before applying the
      # separately requested SSC convention. Other supported families fix phi=1.
      phi = if (family$family == "gaussian") sum(d$w * residuals(dense, type = "response")^2) /
        (N - rank) else 1
      expected = bread[positions, positions, drop = FALSE] * phi
      strategy = if (name == "tuple") "moments" else "compress"
      noK = dbglm(spec[[1L]], data = d, conn = conn, family = family,
        weights = ~w, offset = ~o, retain = "inference", strategy = strategy,
        ssc = list(K.adj = FALSE))
      fit = dbglm(spec[[1L]], data = d, conn = conn, family = family,
        weights = ~w, offset = ~o, retain = "inference", strategy = strategy)
      expect_identical(names(coef(fit)), c("x", "I(x^2)"))
      expect_equal(unname(coef(fit)), unname(coef(dense)[positions]), tolerance = 1e-7)
      expect_equal(fit$n_params_model, rank)
      expect_equal(df.residual(fit), N - rank)
      expect_equal(bare(vcov(noK)), bare(expected), tolerance = 1e-7)
      for (policy in c("full", "nonnested", "none")) {
        K = if (policy == "none") length(positions) + spec[[3L]] else rank
        adjusted = expected * (N - 1) / (N - K)
        V = vcov(fit, vcov = "iid", ssc = list(K.fixef = policy))
        expect_equal(bare(V), bare(adjusted), tolerance = 1e-7)
        expect_equal(attr(V, "df.K"), K)
        expect_equal(attr(V, "df.t"), N - K)
      }
      adjusted = expected * (N - 1) / (N - rank)
      expect_equal(bare(vcov(fit)), bare(adjusted), tolerance = 1e-7)
      expect_equal(bare(vcov(fit, vcov = "iid")), bare(adjusted), tolerance = 1e-7)
      se = sqrt(diag(adjusted))
      selected = summary(noK, vcov = "iid", ssc = list(K.adj = TRUE))
      expect_equal(selected$coeftable[, "std.error"], setNames(se, names(coef(fit))), tolerance = 1e-7)
      critical = if (family$family == "gaussian") qt(.975, N - rank) else qnorm(.975)
      intervals = cbind(coef(dense)[positions] - critical * se, coef(dense)[positions] + critical * se)
      expect_equal(bare(confint(selected)), bare(intervals), tolerance = 1e-7)
      expect_equal(generics::tidy(selected)$std.error, unname(se), tolerance = 1e-7)
      close(noK)
      close(fit)
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
  }
})

local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  set.seed(7162)
  d = expand.grid(g = 1:4, h = 1:3, x = -1:1, rep = 1:5)
  d$w = rep(c(.3, 1, 4), length.out = nrow(d))
  d$o = d$rep / 20
  d$y = rpois(nrow(d), exp(.2 + .2 * d$x + .04 * d$g + d$o))
  reference = glm(y ~ x + factor(g) + factor(h), data = d, family = poisson(),
    weights = w, offset = o, control = glm.control(epsilon = 1e-14))
  X = model.matrix(reference)
  B = solve(crossprod(X, X * (d$w * fitted(reference))))
  scores = rowsum(X * (d$w * (d$y - fitted(reference))), d$g)
  meat = crossprod(scores)
  N = nrow(d)
  G = length(unique(d$g))
  model_K = qr(X)$rank
  # Group g is nested in cluster g; h is not. This expected K is obtained
  # from the independent full design, not the model's returned df metadata.
  nested_K = model_K - (nlevels(factor(d$g)) - 1)
  fit = dbglm(y ~ x | g + h, data = d, conn = conn, family = poisson(),
    weights = ~w, offset = ~o, vcov = ~g, retain = "inference",
    control = dbglm_control(projection_method = "sql"))
  for (policy in c("full", "nonnested", "none")) {
    K = switch(policy, full = model_K, nonnested = nested_K, none = 1L)
    expected = (B %*% meat %*% B)["x", "x"] * G / (G - 1) * (N - 1) / (N - K)
    V = vcov(fit, vcov = ~g, ssc = list(K.fixef = policy))
    expect_equal(as.numeric(V), expected, tolerance = 1e-7)
    expect_equal(attr(V, "df.K"), K)
    expect_equal(attr(V, "df.t"), G - 1)
  }
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(fit, vcov = ~g)), tolerance = 1e-10)
  expected = B["x", "x"] * (N - 1) / (N - model_K)
  expect_equal(as.numeric(vcov(fit, vcov = "iid")), expected, tolerance = 1e-7)
  close(fit)
  # Multiplying all prior weights changes known-dispersion information inversely,
  # even when the solver normalizes by vastly different internal weight scales.
  for (scale in c(1e-6, 1e6)) {
    data = d
    data$w = data$w * scale
    f = dbglm(y ~ x | g + h, data = data, conn = conn, family = poisson(),
      weights = ~w, offset = ~o, retain = "inference")
    expect_equal(as.numeric(vcov(f)) * scale, expected, tolerance = 1e-7)
    expect_equal(as.numeric(vcov(f, vcov = "iid")) * scale, expected, tolerance = 1e-7)
    expect_equal(coef(f)[["x"]], coef(reference)[["x"]], tolerance = 1e-7)
    close(f)
  }
  expect_error(dbglm(y ~ x | g, data = d, conn = conn,
    ssc = list(K.exact = TRUE)), "Unsupported")
  expect_error(dbglm(y ~ x | g, data = d[d$rep == 1L & d$h == 1L & d$x == -1, ],
    conn = conn, fixef_rm = "none"), "degrees of freedom")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
