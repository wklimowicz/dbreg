# The two-cycle window projection and its extrapolation agree with an
# independent weighted R calculation, including a zero-denominator column.
local({
  d = data.frame(id = 1:12, g = rep(1:3, 4), h = rep(1:2, each = 6),
    z1_tilde = sin(1:12), z2_tilde = cos(1:12), zero_tilde = 1,
    w = rep(c(.5, 1, 2), 4), check.names = FALSE)
  names(d)[names(d) == "w"] = "__w"
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "ap_seed", d)
  columns = c("z1_tilde", "z2_tilde", "zero_tilde")
  project = function(X) {
    for (fe in c("g", "h")) {
      for (j in seq_len(ncol(X))) {
        means = tapply(X[, j]*d[["__w"]], d[[fe]], sum) / tapply(d[["__w"]], d[[fe]], sum)
        X[, j] = X[, j] - means[as.character(d[[fe]])]
      }
    }
    X
  }
  A = as.matrix(d[columns]); B = project(A); C = project(B)
  delta = C - B; delta2 = C - 2*B + A
  denominator = colSums(d[["__w"]] * delta2^2)
  alpha = colSums(d[["__w"]] * delta * delta2) / denominator
  alpha[!is.finite(alpha)] = 0
  expected = C - sweep(delta, 2, alpha, "*")
  two = dbreg:::sql_ap_window_irons_tuck("ap_seed", c("id", "g", "h", "__w"), columns, c("g", "h"))
  source = paste0("(", two, ") ap_two")
  moments = DBI::dbGetQuery(conn, dbreg:::sql_ap_irons_tuck_moments(source, columns))
  actual_alpha = dbreg:::ap_irons_tuck_alpha(moments, length(columns))
  expect_equal(actual_alpha, as.numeric(alpha), tolerance = 1e-12)
  actual = DBI::dbGetQuery(conn, dbreg:::sql_ap_irons_tuck_step(source, c("id", "g", "h", "__w"), columns, actual_alpha))
  actual = actual[match(d[["id"]], actual[["id"]]), ]
  expect_equal(as.numeric(as.matrix(actual[columns])), as.numeric(expected), tolerance = 1e-12)
  expect_identical(actual[["__w"]], d[["__w"]])
  expect_identical(actual[["g"]], d[["g"]])
  expect_identical(DBI::dbListTables(conn), "ap_seed")
  for (num in c(0, Inf, NA_real_, 1e300)) {
    row = list(it_num_1 = num, it_den_1 = if (is.finite(num) && num > 0) 1e-300 else 0)
    expect_identical(dbreg:::ap_irons_tuck_alpha(row, 1L), 0)
  }
})

# Shared weighted OLS/IV consumers retain the same inference with acceleration
# on or off. Two-cycle steps must respect odd and one-cycle iteration caps.
local({
  saved = options(dbreg.ap_method = getOption("dbreg.ap_method"),
    dbreg.ap_accelerate = getOption("dbreg.ap_accelerate"),
    dbreg.ap_max_iter = 101L, dbreg.ap_tol = 1e-10)
  on.exit(options(saved), add = TRUE)
  set.seed(913)
  d = data.frame(x = rnorm(240), z = rnorm(240), g = sample(1:9, 240, TRUE),
    h = sample(1:7, 240, TRUE), w = runif(240, .4, 2))
  d[["endog"]] = .7*d[["z"]] + .3*d[["x"]] + rnorm(240)
  d[["y"]] = .5*d[["x"]] + .8*d[["endog"]] + d[["g"]] + sin(d[["h"]]) + rnorm(240)
  for (method in c("join", "window")) {
    options(dbreg.ap_method = method)
    for (iv in c(FALSE, TRUE)) {
      fml = if (iv) y ~ x | g + h | endog ~ z else y ~ x | g + h
      estimator = if (iv) dbivreg else dbreg
      ref = fixest::feols(fml, d, weights = ~w, vcov = "hc1", fixef.rm = "none",
        ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
      for (accelerate in c(FALSE, TRUE)) {
        options(dbreg.ap_accelerate = accelerate)
        fit = estimator(fml, data = d, weights = "w", vcov = "hc1", strategy = "demean")
        b = coef(ref); names(b) = sub("^fit_", "", names(b))
        V = vcov(ref); dimnames(V) = lapply(dimnames(V), function(x) sub("^fit_", "", x))
        terms = names(coef(fit))
        expect_equal(as.numeric(coef(fit)), as.numeric(b[terms]), tolerance = 1e-8)
        expect_equal(as.numeric(vcov(fit)), as.numeric(V[terms, terms]), tolerance = 1e-8)
        expect_true(fit[["ap_iterations"]] <= 101L)
        expect_identical(fit[["ap_accelerated_steps"]] > 0L, accelerate)
      }
    }
  }
  options(dbreg.ap_accelerate = TRUE, dbreg.ap_max_iter = 1L)
  expect_error(dbreg(y ~ x | g + h, data = d, weights = "w", strategy = "demean"), "within 1 iterations")
})
