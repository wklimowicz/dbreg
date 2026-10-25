local({
  d = data.frame(y = rep(0:1, each = 10), x = seq_len(20), g = rep(1:2, each = 10))
  expect_error(dbreg::dbglm(y ~ x | g, data = d, family = stats::binomial()),
    "All observations were removed by fixed-effect singleton/separation pruning.", fixed = TRUE)
  warnings = new.env(parent = emptyenv())
  warnings[["messages"]] = character()
  tryCatch(withCallingHandlers(dbreg::dbglm(y ~ x, data = 12),
    warning = function(w) {
      warnings[["messages"]] = c(warnings[["messages"]], conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  expect_identical(warnings[["messages"]], character())
})
