# Post-fit covariance changes must refresh interval uncertainty in fitting
# coordinates. A large raw-origin covariance cannot safely be transformed back.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(6283)
  dx = rep(seq(-1, 1, by = 1/16), 5)
  n = length(dx)
  d = data.frame(dx = dx, w = runif(n, .4, 2), g = rep(1:11, length.out = n),
    h = rep(1:7, length.out = n), o = .1*sin(seq_len(n)))
  d[["y"]] = 1 + 2*dx + d[["o"]] + rnorm(n, sd = exp(dx))
  X = cbind(1, dx)
  bread = solve(crossprod(X, d[["w"]]*X))
  beta = bread %*% crossprod(X, d[["w"]]*(d[["y"]]-d[["o"]]))
  e = as.numeric(d[["y"]]-d[["o"]]-X %*% beta)
  score = X*as.numeric(d[["w"]]*e)
  rss = sum(d[["w"]]*e^2)
  clusters = function(group) crossprod(rowsum(score, group))
  psd = function(V) {
    eig = eigen(V, symmetric = TRUE)
    tcrossprod(sweep(eig[["vectors"]], 2L, sqrt(pmax(eig[["values"]], 0)), "*"))
  }
  cases = list(
    list(vcov = "hc1", V = n/(n-2)*bread %*% crossprod(score) %*% bread, df = n-2),
    list(vcov = ~g, V = 11/10*(n-1)/(n-2)*bread %*% clusters(d[["g"]]) %*% bread, df = 10),
    list(vcov = ~g+h, V = 7/6*(n-1)/(n-2)*bread %*%
      (clusters(d[["g"]])+clusters(d[["h"]])-clusters(interaction(d[["g"]], d[["h"]]))) %*% bread, df = 6),
    list(vcov = "iid", V = rss/(n-2)*bread, df = n-2),
    list(vcov = "hc1", ssc = list(K.adj = FALSE, G.adj = FALSE, t.df = 37),
      V = bread %*% crossprod(score) %*% bread, df = 37)
  )
  for (location in c(0, 1e14)) for (strategy in c("moments", "compress")) {
    d[["x"]] = location + dx
    for (initial in c("iid", "hc1")) {
      fit = dbreg(y ~ x, conn = conn, data = d, weights = "w", offset = "o",
        strategy = strategy, retain = "inference", vcov = initial)
      initial_V = vcov(fit)
      initial_prediction_V = fit[["prediction_coordinates"]][["covariance"]]
      initial_intervals = predict(fit, d, interval = "confidence")
      local_center = fit[["prediction_coordinates"]][["centers"]][2L] - location
      C = matrix(c(1, 0, local_center, 1), 2, 2)
      W = cbind(1, dx-local_center)
      for (case in cases) {
        updated = suppressWarnings(summary(fit, vcov = case[["vcov"]], ssc = case[["ssc"]]))
        predicted = predict(updated, d, interval = "confidence")
        V = C %*% case[["V"]] %*% t(C)
        if (inherits(case[["vcov"]], "formula") && length(all.vars(case[["vcov"]])) == 2L) V = psd(V)
        expected = as.numeric(X %*% beta) + d[["o"]]
        half = qt(.975, case[["df"]])*sqrt(rowSums((W %*% V)*W))
        A = matrix(c(1, 0, -(location+local_center), 1), 2, 2)
        expect_equal(predicted[["fit"]], expected, tolerance = 1e-9)
        expect_equal(predicted[["lwr"]], expected-half, tolerance = 1e-9)
        expect_equal(predicted[["upr"]], expected+half, tolerance = 1e-9)
        expect_equal(as.numeric(vcov(updated)), as.numeric(A %*% V %*% t(A)), tolerance = 1e-9)
      }
      # Disabling two-way repair retains the unmodified fitting-coordinate V.
      uncorrected = suppressWarnings(summary(fit, vcov = ~g+h, vcov_fix = FALSE))
      local_center = fit[["prediction_coordinates"]][["centers"]][2L] - location
      C = matrix(c(1, 0, local_center, 1), 2, 2)
      expect_equal(as.numeric(uncorrected[["prediction_coordinates"]][["covariance"]]),
        as.numeric(C %*% cases[[3L]][["V"]] %*% t(C)), tolerance = 1e-9)
      expect_false(isTRUE(attr(vcov(uncorrected), "repair")[["repaired"]]))
      expect_identical(vcov(fit), initial_V)
      expect_identical(fit[["prediction_coordinates"]][["covariance"]], initial_prediction_V)
      expect_identical(predict(fit, d, interval = "confidence"), initial_intervals)
      close(fit)
    }
  }
})

# CRE means have their own raw coordinate origins; refresh their covariance
# without treating the already-centered mean recipe as a raw source column.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(4773)
  d = data.frame(g = rep(1:8, each = 15), w = runif(120, .3, 2))
  d[["dx"]] = sample(-16:16, 120, TRUE)/16
  d[["x"]] = 1e14 + d[["dx"]]
  d[["y"]] = 1 + 2*d[["dx"]] + sin(d[["g"]]) + rnorm(120, sd = exp(d[["dx"]]))
  means = ave(seq_len(nrow(d)), d[["g"]], FUN = function(i) weighted.mean(d[["dx"]][i], d[["w"]][i]))
  X = cbind(1, d[["dx"]], means)
  B = solve(crossprod(X, d[["w"]]*X))
  beta = B %*% crossprod(X, d[["w"]]*d[["y"]])
  e = as.numeric(d[["y"]]-X %*% beta)
  V = 120/117*B %*% crossprod(X*as.numeric(d[["w"]]*e)) %*% B
  fit = dbreg(y ~ x | g, conn = conn, data = d, weights = "w", strategy = "mundlak", retain = "inference")
  updated = summary(fit, vcov = "hc1")
  actual = predict(updated, d, interval = "confidence")
  expected = as.numeric(X %*% beta)
  half = qt(.975, 117)*sqrt(rowSums((X %*% V)*X))
  expect_equal(actual[["fit"]], expected, tolerance = 1e-9)
  expect_equal(actual[["lwr"]], expected-half, tolerance = 1e-9)
  expect_equal(actual[["upr"]], expected+half, tolerance = 1e-9)
  A = diag(3); A[1L, 2:3] = -1e14
  expect_equal(as.numeric(vcov(updated, fe = TRUE)), as.numeric(A %*% V %*% t(A)), tolerance = 1e-9)
  close(fit)
})

# An actually indefinite two-way covariance is repaired on the reported
# coefficients; disabling repair preserves the same unrepaired score covariance.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(1)
  n = 36L
  dx = rep(seq(-1, 1, length.out = 9), 4)
  d = data.frame(x = 1e14+dx, g = rep(1:3, each = 12), h = rep(1:4, 9))
  X = cbind(1, dx)
  d[["y"]] = as.numeric(X %*% c(1, 2)) + rnorm(n)
  B = solve(crossprod(X))
  beta = B %*% crossprod(X, d[["y"]])
  score = X*as.numeric(d[["y"]]-X %*% beta)
  meat = function(g) crossprod(rowsum(score, g))
  V = 3/2*(n-1)/(n-2)*B %*% (meat(d[["g"]])+meat(d[["h"]])-meat(interaction(d[["g"]], d[["h"]]))) %*% B
  eig = eigen(V, symmetric = TRUE)
  expect_true(min(eig[["values"]]) < -1e-5)
  fit = dbreg(y ~ x, conn = conn, data = d, strategy = "moments", retain = "inference")
  repaired = suppressWarnings(summary(fit, vcov = ~g+h))
  unrepaired = suppressWarnings(summary(fit, vcov = ~g+h, vcov_fix = FALSE))
  A = matrix(c(1, 0, -1e14, 1), 2, 2)
  # Reported covariance is repaired in public coordinates, as fixest does;
  # intervals use the repair in fitting coordinates.
  public = A %*% V %*% t(A)
  pe = eigen(public, symmetric = TRUE)
  expect_equal(as.numeric(vcov(repaired)), as.numeric(tcrossprod(sweep(pe[["vectors"]], 2L, sqrt(pmax(pe[["values"]], 0)), "*"))), tolerance = 1e-9)
  expect_equal(as.numeric(vcov(unrepaired)), as.numeric(A %*% V %*% t(A)), tolerance = 1e-9)
  expect_equal(as.numeric(unrepaired[["prediction_coordinates"]][["covariance"]]), as.numeric(V), tolerance = 1e-9)
  expect_false(isTRUE(attr(vcov(unrepaired), "repair")[["repaired"]]))
  # Interval widths of a repaired covariance at a 1e14 origin are beyond double
  # precision; only the point predictions are checked here.
  actual = predict(repaired, d, interval = "confidence")
  expect_equal(actual[["fit"]], as.numeric(X %*% beta), tolerance = 1e-9)
  close(fit)
})
