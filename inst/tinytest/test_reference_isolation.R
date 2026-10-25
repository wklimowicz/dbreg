# Reference calls must stay independent when dbreg masks fixest exports.
local({
  reference_symbols = c("feols", "feglm", "ssc", "coeftable", "se", "pvalue", "tstat", "degrees_freedom", "obs")
  unqualified = function(file) {
    tokens = getParseData(parse(file, keep.source = TRUE))
    selected = tokens[tokens[["terminal"]] &
      tokens[["token"]] %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL") &
      tokens[["text"]] %in% reference_symbols, ]
    selected[vapply(selected[["parent"]], function(parent) {
      !any(tokens[tokens[["parent"]] == parent, "token"] %in%
        c("NS_GET", "NS_GET_INT", "'$'", "'@'"))
    }, logical(1)), ]
  }
  files = list.files(system.file("tinytest", package = "dbreg"), pattern = "\\.R$", full.names = TRUE)
  for (file in files) expect_equal(nrow(unqualified(file)), 0L, info = basename(file))

  data = data.frame(x = seq_len(30), y = sin(seq_len(30)), b = rep(0:1, 15))
  ols = fixest::feols(y ~ x, data, notes = FALSE)
  glm = fixest::feglm(b ~ x, data, family = stats::binomial(), notes = FALSE)
  expect_true(inherits(ols, "fixest") && !inherits(ols, "dbreg"))
  expect_true(inherits(glm, "fixest") && !inherits(glm, "dbreg"))
  expect_true(!identical(fixest::feols, dbreg::feols))
  expect_true(!identical(fixest::feglm, dbreg::feglm))
})
