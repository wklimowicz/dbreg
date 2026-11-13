local({
  set.seed(11149)
  d = expand.grid(g = 1:3, h = 1:2, x = c(-1, 1), z = c(-1, 1), rep = 1:8)
  n = nrow(d)
  d[["shifted"]] = 1e12 + d[["z"]]
  d[["constant"]] = 3
  d[["w"]] = runif(n, .5, 2)
  d[["cl"]] = d[["g"]] + d[["h"]]
  eta = .5 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) * d[["z"]]
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  for (family in list(gaussian(), binomial("logit"))) {
    mu = family[["linkinv"]](eta)
    d[["y"]] = if (family[["family"]] == "gaussian") rnorm(n, mu) else rbinom(n, 1, mu)
    model = dbglm(y ~ x | (g^h)[shifted, constant], conn = conn, data = d,
      family = family, weights = ~w, vcov = ~cl, strategy = "compress",
      ssc = list(K.adj = FALSE, G.adj = TRUE),
      control = dbglm_control(score_tol = 1e-9))
    ref = suppressWarnings(glm(y ~ x + interaction(g, h)*z, data = d, weights = w,
      family = family, control = glm.control(epsilon = 1e-13)))
    X = model.matrix(ref)
    mu = fitted(ref)
    derivative = family[["mu.eta"]](predict(ref, type = "link"))
    var = family[["variance"]](mu)
    bread = solve(crossprod(X, X * (d[["w"]] * derivative^2 / var)))
    scores = rowsum(X * (d[["w"]] * derivative * (d[["y"]] - mu) / var), d[["cl"]])
    G = nrow(scores)
    V = bread %*% crossprod(scores) %*% bread * G / (G - 1)
    expect_equal(coef(model)[["x"]], coef(ref)[["x"]], tolerance = 1e-8)
    expect_equal(vcov(model)["x", "x"], V["x", "x"], tolerance = 1e-8)
    expect_equal(model[["deviance"]], deviance(ref), tolerance = 1e-10)
    expect_identical(model[["projection_method"]], "sql_cg")
    expect_equal(model[["n_cells"]], 24)
    expect_equal(model[["n_params_model"]], 19)
    expect_true(model[["converged"]] && max(model[["score_max"]], model[["fe_score_max"]]) <= 1e-9)
    expect_identical(DBI::dbListTables(conn), character())
  }
})
