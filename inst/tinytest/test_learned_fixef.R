local({
  keys = data.frame(a = c("\u00e9:1", ""), b = c("9007199254740993", "9007199254740994"))
  dictionary = list(columns = names(keys), keys = dbreg:::composite_key(keys, names(keys),
    c("character", "integer")))
  expect_identical(dbreg:::learned_effect_labels(dictionary),
    c('("\u00e9:1", "9007199254740993")', '("", "9007199254740994")'))
})

if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("fixest", quietly = TRUE)) {
  set.seed(61)
  d = expand.grid(g = c('group "one"', "group:two"), h = 1:3, t = 1:3, x = -1:1, rep = 1:3)
  d$o = (d$t + d$h) / 10
  d$y = rpois(nrow(d), exp(d$x / 4 + d$o))
  d$z = log1p(d$y)
  keys = paste0("(", encodeString(as.character(d$g), quote = '"'), ", ", encodeString(as.character(d$h), quote = '"'), ")")
  for (fit in list(
    dbreg(z ~ x | (g^h)[t], data = d, retain = "prediction", offset = ~o, strategy = "compress"),
    dbglm(y ~ x | (g^h)[t], data = d, retain = "prediction", offset = ~o, family = poisson())
  )) {
    close(fit)
    saved = fixest::fixef(fit)
    expect_identical(names(saved), c("g^h", "g^h[[t]]"))
    expect_equal(length(saved[[1L]]), 6L)
    expect_true(all(unique(keys) %in% names(saved[[1L]])))
    eta = unname(coef(fit)[["x"]] * d$x + saved[[1L]][keys] + d$t * saved[[2L]][keys] + d$o)
    expect_equal(eta, if (inherits(fit, "dbglm")) predict(fit, d, type = "link") else predict(fit, d),
      tolerance = 1e-10)
  }
  additive = dbglm(y ~ x | g + h, data = d, retain = "prediction", offset = ~o, family = poisson())
  close(additive)
  saved = fixest::fixef(additive)
  expect_identical(names(saved), c("g", "h"))
  eta = unname(coef(additive)[["x"]] * d$x + saved$g[as.character(d$g)] + saved$h[as.character(d$h)] + d$o)
  expect_equal(eta, predict(additive, d, type = "link"), tolerance = 1e-10)
  expect_identical(fixest::fixef(dbreg(y ~ x, data = d)), list())
  expect_error(fixest::fixef(dbreg(y ~ x | g, data = d)), "portable learned-effect state")
}
