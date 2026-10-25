# A many-level factor must not create a dense cells-by-FE panel in R.
if (isTRUE(capabilities("profmem"))) local({
  set.seed(732)
  n = 3000L
  groups = 200L
  d = data.frame(x = rnorm(n), f = factor(rep(seq_len(groups), length.out = n)))
  d[["y"]] = .4*d[["x"]]+as.numeric(d[["f"]])/groups+rnorm(n)
  # Initialize sparse Matrix methods before measuring the actual fit.
  invisible(dbreg(y ~ x | f, data = d[1:400, ], strategy = "compress"))
  path = tempfile("dbreg_compression_allocations_")
  on.exit(unlink(path), add = TRUE)
  Rprofmem(path)
  fit = tryCatch(dbreg(y ~ x | f, data = d, strategy = "compress"), finally = Rprofmem(NULL))
  lines = readLines(path)
  allocations = as.numeric(sub(" .*", "", lines[grepl("^[0-9]+", lines)]))
  # Dense coefficient matrices are O(G^2); an O(N*G) design allocation is a regression.
  expect_true(max(c(0, allocations)) < 4*n*groups)
  reference = fixest::feols(y ~ x | f, d, fixef.rm = "none", notes = FALSE)
  expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-8)
})
