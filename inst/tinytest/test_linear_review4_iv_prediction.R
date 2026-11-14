# Structural predictions and confidence intervals must use centered 2SLS
# coefficients and the projected structural meat, preserving raw public coef.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(1441)
  cells = 60L
  cell = data.frame(x = 1e14 + sample(-60:60, cells, TRUE)/64,
    z = rnorm(cells), q = rnorm(cells))
  cell[["endog"]] = round(8*(.7*cell[["z"]] + .3*cell[["q"]] + rnorm(cells)))/8
  d = cell[rep(seq_len(cells), each = 4L), , drop = FALSE]
  n = nrow(d)
  d[["dx"]] = d[["x"]] - 1e14
  d[["cluster"]] = rep(1:20, length.out = n)
  d[["o"]] = rnorm(n, sd = .2)
  d[["y"]] = d[["o"]] + 2*d[["dx"]] + .5*d[["endog"]] + rnorm(n, sd = .1)
  X = cbind(1, d[["dx"]], d[["endog"]])
  Z = cbind(1, d[["dx"]], d[["z"]], d[["q"]])
  y = d[["y"]] - d[["o"]]
  for (weighted in c(FALSE, TRUE)) {
    d[["w"]] = if (weighted) runif(n, .3, 2) else 1
    w = d[["w"]]
    P = crossprod(X, w*Z) %*% solve(crossprod(Z, w*Z))
    bread = solve(P %*% crossprod(Z, w*X))
    beta = bread %*% P %*% crossprod(Z, w*y)
    residual = as.numeric(y-X %*% beta)
    scores = Z*as.numeric(w*residual)
    expected = as.numeric(X %*% beta) + d[["o"]]
    DBI::dbWriteTable(conn, "review4_iv_prediction", d, overwrite = TRUE)
    for (strategy in c("moments", "compress")) for (kind in list("iid", "hc1", ~cluster)) {
      fit = dbivreg(y ~ x | endog ~ z + q, conn = conn, table = "review4_iv_prediction",
        weights = if (weighted) "w" else NULL, offset = "o", strategy = strategy, vcov = kind)
      centered = dbivreg(y ~ dx | endog ~ z + q, conn = conn, table = "review4_iv_prediction",
        weights = if (weighted) "w" else NULL, offset = "o", strategy = strategy, vcov = kind)
      G = length(unique(d[["cluster"]]))
      K = ncol(X)
      V = if (identical(kind, "iid")) sum(w*residual^2)/(n-K)*bread else {
        meat = if (identical(kind, "hc1")) crossprod(scores) else crossprod(rowsum(scores, d[["cluster"]]))
        correction = if (identical(kind, "hc1")) n/(n-K) else G/(G-1)*(n-1)/(n-K)
        correction*bread %*% P %*% meat %*% t(P) %*% bread
      }
      se = sqrt(rowSums((X %*% V)*X))
      df = if (inherits(kind, "formula")) G-1 else n-K
      intervals = predict(fit, d, interval = "confidence")
      expect_equal(predict(fit, d), expected, tolerance = 1e-9)
      expect_equal(predict(fit, d), predict(centered, d), tolerance = 1e-9)
      expect_equal(intervals, predict(centered, d, interval = "confidence"), tolerance = 1e-9)
      expect_equal(intervals[["fit"]], expected, tolerance = 1e-9)
      expect_equal(intervals[["lwr"]], expected-qt(.975, df)*se, tolerance = 1e-9)
      expect_equal(intervals[["upr"]], expected+qt(.975, df)*se, tolerance = 1e-9)
      expect_equal(unname(coef(fit)[-1L]), as.numeric(beta[-1L]), tolerance = 1e-9)
    }
  }
})
