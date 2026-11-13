# Exact FE dependencies hidden by SQL summation roundoff (Spark) must still be
# dropped from the dense Schur basis, else a later IRLS iteration cannot factor
# the cached basis. Structure mirrors the real experience models: a composite
# FE (spine x band x year) whose band is also determined by district.
library(tinytest)
fe_pivot_rank = dbreg:::fe_pivot_rank
expect_equal(fe_pivot_rank(c(1, 0.4, 7e-11, 5e-11)), 2L)
expect_equal(fe_pivot_rank(c(1, 0.4, 0.2)), 3L)
# Weakly identified but not a jump: kept.
expect_equal(fe_pivot_rank(c(1, 1e-3, 1e-6, 1e-9)), 4L)
expect_equal(fe_pivot_rank(numeric()), 0L)

set.seed(11)
n = 4000
district = sample(40, n, replace = TRUE)
band = (district %% 4) + 1
year = sample(3, n, replace = TRUE)
spine = sample(8, n, replace = TRUE)
f = list(sby = interaction(spine, band, year, drop = TRUE), district = factor(district),
  exper = factor(sample(10, n, replace = TRUE)), sex = factor(sample(2, n, replace = TRUE)))
D = do.call(cbind, lapply(f, function(x) Matrix::sparse.model.matrix(~ x - 1)))
w = runif(n, 0.05, 0.25)
H = as.matrix(Matrix::crossprod(D, D * w))
sizes = vapply(f, nlevels, 1L)
big = seq_len(sizes[[1]]); rest = setdiff(seq_len(ncol(H)), big)
deps = ncol(H) - qr(H, tol = 1e-9)$rank
expect_true(deps > length(f) - 1L)
noisy = function(seed) {
  set.seed(seed)
  E = H[rest, rest]
  N = matrix(rnorm(length(E)), nrow(E)) * 1e-10 * sqrt(outer(diag(E), diag(E))) * (E != 0)
  list(Dg = diag(H)[big], Cm = Matrix::Matrix(H[big, rest], sparse = TRUE), E = Matrix::Matrix(E + (N + t(N)) / 2),
    sizes = unname(sizes), big = big, rest = rest, diag = diag(H), n = ncol(H))
}
cache = new.env()
first = dbreg:::fe_solver(noisy(1), cache)
expect_equal(length(cache[["moment_basis"]]), length(rest) - deps)
second = tryCatch(dbreg:::fe_solver(noisy(2), cache), error = function(e) NULL)
expect_false(is.null(second))
