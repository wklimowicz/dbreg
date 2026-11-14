# Analytic one-/balanced-two-family projections recover the within model,
# including response translations and FE-only numeric controls.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(19203)
  d = expand.grid(g = 1:7, h = 1:5, repeat_id = 1:8)
  n = nrow(d)
  d[["x"]] = 1e3 + runif(n, -1, 1)
  d[["dx"]] = d[["x"]] - 1e3
  d[["z"]] = rnorm(n)
  d[["w"]] = runif(n, .2, 2)
  d[["absorbed"]] = 1e3 + d[["g"]]
  d[["y"]] = 1e3 + 2*d[["dx"]] + .4*d[["z"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n, sd = .05)
  d[["dy"]] = d[["y"]] - 1e3
  DBI::dbWriteTable(conn, "review4_analytic", d)
  for (two in c(FALSE, TRUE)) for (weighted in c(FALSE, TRUE)) {
    if (two && weighted) next
    fml = if (two) y ~ x + z + absorbed | g + h else y ~ x + z + absorbed | g
    ref = lm(if (two) dy ~ dx + z + factor(g) + factor(h) else dy ~ dx + z + factor(g),
      d, weights = if (weighted) w else NULL)
    for (kind in list("iid", "hc1", ~g)) {
      fit = dbreg(fml, conn = conn, table = "review4_analytic", strategy = "demean",
        weights = if (weighted) "w" else NULL, vcov = kind)
      X = model.matrix(ref)
      w = if (weighted) d[["w"]] else rep(1, n)
      bread = solve(crossprod(X, w*X))
      scores = X*as.numeric(w*residuals(ref))
      K = ncol(X)
      V = if (identical(kind, "iid")) vcov(ref) else if (identical(kind, "hc1")) {
        n/(n-K)*bread %*% crossprod(scores) %*% bread
      } else {
        G = length(unique(d[["g"]]))
        G/(G-1)*(n-1)/(n-K)*bread %*% crossprod(rowsum(scores, d[["g"]])) %*% bread
      }
      expect_identical(names(coef(fit)), c("x", "z"))
      expect_identical(fit[["collin.var"]], "absorbed")
      expect_equal(unname(coef(fit)), unname(coef(ref)[c("dx", "z")]), tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(V[c("dx", "z"), c("dx", "z")]), tolerance = 1e-9)
    }
  }
  # A source column can be both a structural regressor and an FE/cluster key.
  fit = dbreg(y ~ x + g | g, conn = conn, table = "review4_analytic", weights = "w", vcov = ~g)
  ref = lm(dy ~ dx + factor(g), d, weights = w)
  expect_identical(names(coef(fit)), "x")
  expect_equal(unname(coef(fit)), unname(coef(ref)[["dx"]]), tolerance = 1e-9)
})
