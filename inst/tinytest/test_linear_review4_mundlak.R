# The oracle constructs CRE means on small centered source values. Public
# coefficients retain raw identities; prediction uses the saved coordinates.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(46399)
  n = 240L
  d = data.frame(g = rep(1:8, each = 30), h = rep(1:5, length.out = n),
    z = rnorm(n), w = runif(n, .2, 2))
  d[["x"]] = 1e14 + runif(n, -1, 1)
  d[["dx"]] = d[["x"]] - 1e14
  d[["y"]] = 2*d[["dx"]] + .4*d[["z"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n, sd = .1)
  for (two in c(FALSE, TRUE)) for (weighted in c(FALSE, TRUE)) {
    w = if (weighted) d[["w"]] else rep(1, n)
    means = function(value, group) ave(seq_len(n), group,
      FUN = function(rows) weighted.mean(value[rows], w[rows]))
    X = cbind("(Intercept)" = 1, x = d[["dx"]], z = d[["z"]],
      x_bar_g = means(d[["dx"]], d[["g"]]), z_bar_g = means(d[["z"]], d[["g"]]))
    if (two) X = cbind(X, x_bar_h = means(d[["dx"]], d[["h"]]), z_bar_h = means(d[["z"]], d[["h"]]))
    ref = lm.wfit(X, d[["y"]], w)
    bread = solve(crossprod(X, w*X))
    K = ncol(X)
    for (kind in list("iid", "hc1", ~g)) {
      DBI::dbWriteTable(conn, "review4_mundlak", d, overwrite = TRUE)
      fit = dbreg(if (two) y ~ x + z | g + h else y ~ x + z | g,
        conn = conn, table = "review4_mundlak", weights = if (weighted) "w" else NULL,
        strategy = "mundlak", vcov = kind)
      scores = X*as.numeric(w*ref[["residuals"]])
      V = if (identical(kind, "iid")) sum(w*ref[["residuals"]]^2)/(n-K)*bread else if (identical(kind, "hc1")) {
        n/(n-K)*bread %*% crossprod(scores) %*% bread
      } else {
        G = length(unique(d[["g"]]))
        G/(G-1)*(n-1)/(n-K)*bread %*% crossprod(rowsum(scores, d[["g"]])) %*% bread
      }
      expect_equal(unname(coef(fit, fe = TRUE)[-1L]), unname(ref[["coefficients"]][-1L]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit, fe = TRUE)[-1L, -1L]), as.numeric(V[-1L, -1L]), tolerance = 1e-8)
      source_coordinates = diag(K)
      source_coordinates[1L, grepl("^x($|_bar_)", colnames(X))] = -1e14
      expect_equal(unname(coef(fit, fe = TRUE)), as.numeric(source_coordinates %*% ref[["coefficients"]]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit, fe = TRUE)), as.numeric(source_coordinates %*% V %*% t(source_coordinates)), tolerance = 1e-8)
      expected = drop(X %*% ref[["coefficients"]])
      predicted = predict(fit, d, interval = "confidence")
      se = sqrt(rowSums((X %*% V)*X))
      df = if (inherits(kind, "formula")) 7 else n-K
      expect_equal(predicted[["fit"]], expected, tolerance = 1e-8)
      expect_equal(predicted[["lwr"]], expected-qt(.975, df)*se, tolerance = 1e-8)
      expect_equal(predicted[["upr"]], expected+qt(.975, df)*se, tolerance = 1e-8)
      restored = unserialize(serialize(fit, NULL))
      expect_equal(predict(restored, d), expected, tolerance = 1e-8)
    }
  }
})
