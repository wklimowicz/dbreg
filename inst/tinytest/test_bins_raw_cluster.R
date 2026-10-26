# A covariance request names raw source columns; the inner fit reads safe aliases.
local({
  set.seed(41001)
  d = data.frame(x = stats::runif(240), g = rep(seq_len(8), each = 30))
  d[["bin2"]] = stats::rnorm(240) # Forces collision-safe basis/source allocation.
  d[["cluster id"]] = rep(seq_len(12), 20)
  d[["cluster^id"]] = d[["cluster id"]]
  d[["y"]] = 2 + d[["x"]] + 0.7 * d[["bin2"]] + sin(d[["g"]]) + stats::rnorm(240)
  breaks = seq(0, 1, length.out = 5)
  d[["b"]] = cut(d[["x"]], breaks, include.lowest = TRUE, right = FALSE)
  knots = as.numeric(tapply(d[["x"]], d[["b"]], max))[1:3]
  d[["k1"]] = pmax(d[["x"]] - knots[1], 0)
  d[["k2"]] = pmax(d[["x"]] - knots[2], 0)
  d[["k3"]] = pmax(d[["x"]] - knots[3], 0)
  for (with_fe in c(FALSE, TRUE)) {
    formula = if (with_fe) y ~ x + bin2 | g else y ~ x + bin2
    for (smooth in c(FALSE, TRUE)) {
      reference_formula = if (smooth) {
        if (with_fe) y ~ x + k1 + k2 + k3 + bin2 + factor(g) else y ~ x + k1 + k2 + k3 + bin2
      } else {
        if (with_fe) y ~ b + bin2 + factor(g) else y ~ b + bin2
      }
      reference = stats::lm(reference_formula, data = d)
      X = stats::model.matrix(reference)
      inverse = solve(crossprod(X))
      scores = rowsum(X * as.numeric(stats::residuals(reference)), d[["cluster id"]])
      groups = nrow(scores)
      covariance = inverse %*% crossprod(scores) %*% inverse * groups / (groups - 1) * (nrow(X) - 1) / (nrow(X) - ncol(X))
      for (cluster in c("cluster id", "cluster^id")) {
        request = stats::as.formula(paste0("~`", cluster, "`"))
        result = dbreg::dbbinsreg(formula, data = d, nbins = 4, binspos = breaks,
          points = if (smooth) c(1, 1) else c(0, 0), sample_fit = FALSE,
          vcov = request, plot = FALSE)
        model = result[["model"]]
        # Absorbed FE exposes slopes; explicit bin models retain the joint V.
        indices = if (smooth && with_fe) 2:6 else seq_len(ncol(X))
        V = if (!smooth && with_fe) stats::vcov(model, fe = TRUE) else stats::vcov(model)
        beta = if (!smooth && with_fe) stats::coef(model, fe = TRUE) else stats::coef(model)
        expect_equal(as.numeric(V), as.numeric(covariance[indices, indices, drop = FALSE]), tolerance = 1e-7)
        expect_equal(as.numeric(beta), as.numeric(stats::coef(reference)[indices]), tolerance = 1e-7)
        expect_identical(model[["inference_result"]][["clusters"]], list(cluster))
        expect_identical(all.vars(result[["opt"]][["vcov"]]), cluster)
        expect_equal(attr(V, "n_clusters"), groups)
      }
    }
    # Separate points/line calls share the translated request. Line-only CI is
    # disabled by the existing planner; its IID policy must remain unchanged.
    separated = dbreg::dbbinsreg(formula, data = d, nbins = 4, binspos = breaks,
      points = c(0, 0), line = c(1, 1), sample_fit = FALSE, vcov = ~`cluster id`, plot = FALSE)
    expect_identical(separated[["model"]][["points"]][["cluster_var"]], "cluster id")
    expect_identical(separated[["model"]][["line"]][["vcov_type"]], "iid")
    expect_identical(all.vars(separated[["opt"]][["vcov"]]), "cluster id")
  }
  expect_error(dbreg::dbbinsreg(y ~ x, data = d, vcov = ~absent, plot = FALSE), "absent|column|found")
  expect_error(dbreg::dbbinsreg(y ~ x, data = d, vcov = ~g + `cluster id`, plot = FALSE), "single-variable")
  expect_error(dbreg::dbbinsreg(y ~ x, data = d, vcov = 4, plot = FALSE), "character string or a formula")
})
