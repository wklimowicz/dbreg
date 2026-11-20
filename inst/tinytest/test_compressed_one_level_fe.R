d = data.frame(x = rep(1:3, 40), g = "only", h = rep(letters[1:4], each = 30))
d$y = 2 + .7*d$x + match(d$h, letters) + sin(seq_len(nrow(d)))
for (f in list(y ~ x | g, y ~ 0 + x | g, y ~ x | g + h, y ~ 0 + x | g + h)) {
  ref = lm(if ("h" %in% all.vars(f)) y ~ x + h else y ~ x, d)
  fit = dbreg(f, data = d, strategy = "compress")
  expect_equal(unname(coef(fit)["x"]), unname(coef(ref)["x"]))
  expect_equal(predict(fit, d), unname(fitted(ref)))
  other = d
  other$g[1:2] = c("unseen", NA)
  expect_true(all(is.na(predict(fit, other)[1:2])))
}
