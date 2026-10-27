# Wire decoding is independent of numerical assembly. Weighted incidence
# products supply an oracle independent of either transport representation.
local({
  groups = list(
    c("9007199254740993", "9007199254740994", "9007199254740993", "9007199254740994", "9007199254740995", "9007199254740995"),
    c("west", "east", "east", "west", "east", "west"),
    c("u", "v", "v", "u", "u", "v")
  )
  w = c(0.2, 2, 0.5, 3, 0.8, 1.1)
  V = cbind(c(1, 3, -2, 0.5, 4, 2), c(0.5, -1, 2, 3, 1, -2))
  extra = c(0.3, -0.2, 0.7, -1, 0.4, 0.8)
  K = length(groups)
  rows = list()
  for (i in seq_len(K)) for (id in sort(unique(groups[[i]]))) {
    take = groups[[i]] == id
    rows[[length(rows) + 1L]] = data.frame("__set" = i, "__i" = id, "__j" = NA_character_,
      "__h" = sum(w[take]), "__b1" = sum(w[take] * V[take, 1]), "__b2" = sum(w[take] * V[take, 2]),
      "__ex1" = sum(extra[take]), check.names = FALSE)
  }
  pairs = utils::combn(K, 2L, simplify = FALSE)
  for (k in seq_along(pairs)) {
    pair = pairs[[k]]
    observed = unique(data.frame(a = groups[[pair[1]]], b = groups[[pair[2]]]))
    for (j in seq_len(nrow(observed))) {
      take = groups[[pair[1]]] == observed[["a"]][j] & groups[[pair[2]]] == observed[["b"]][j]
      rows[[length(rows) + 1L]] = data.frame("__set" = K + k, "__i" = observed[["a"]][j], "__j" = observed[["b"]][j],
        "__h" = sum(w[take]), "__b1" = NA_real_, "__b2" = NA_real_, "__ex1" = NA_real_, check.names = FALSE)
    }
  }
  wide = do.call(rbind, rows)
  set.seed(90001)
  wide = wide[sample.int(nrow(wide)), , drop = FALSE]
  globals = data.frame(example = 2)
  attr(wide, "dbglm_globals") = globals
  marginal = which(wide[["__set"]] <= K)
  packed = wide[, c("__set", "__i", "__j", "__h"), drop = FALSE]
  attr(packed, "fe_values") = list(rows = marginal, values = as.matrix(wide[marginal, c("__b1", "__b2", "__ex1"), drop = FALSE]))
  attr(packed, "dbglm_globals") = globals
  D = do.call(cbind, lapply(groups, function(g) stats::model.matrix(~factor(g) - 1)))
  expected_H = crossprod(D, w * D)
  expected_B = crossprod(D, w * V)
  expected_extra = crossprod(D, extra)
  decoded = dbreg:::decode_fe_moments(wide, K, 2L, 1L)
  compact = dbreg:::decode_fe_moments(packed, K, 2L, 1L)
  expect_equal(compact, decoded)
  expect_identical(names(decoded), c("info", "pairs", "rows", "n_extra"))
  expect_identical(decoded[["info"]][[1L]][["id"]], sort(unique(groups[[1L]])))
  expect_identical(attr(wide, "dbglm_globals"), globals)
  expect_identical(attr(packed, "fe_values")[["rows"]], marginal)
  for (packet in list(decoded, compact)) {
    actual = dbreg:::assemble_fe_moments(packet)
    expect_equal(as.numeric(dbreg:::fe_multiply(actual[["sys"]], diag(ncol(D)))), as.numeric(expected_H), tolerance = 1e-13)
    expect_equal(as.numeric(actual[["B"]]), as.numeric(expected_B), tolerance = 1e-13)
    expect_equal(as.numeric(actual[["B_extra"]]), as.numeric(expected_extra), tolerance = 1e-13)
    expect_identical(dim(actual[["B"]]), dim(expected_B))
    expect_equal(actual[["fe_score_max"]], max(abs(expected_extra[, 1L]) / sqrt(diag(expected_H))), tolerance = 1e-13)
    expect_equal(actual[["rows"]], nrow(wide))
  }
  empty = wide[, c("__set", "__i", "__j", "__h"), drop = FALSE]
  zero = dbreg:::assemble_fe_moments(dbreg:::decode_fe_moments(empty, K, 0L, 0L))
  expect_identical(dim(zero[["B"]]), c(ncol(D), 0L))
  expect_identical(dim(zero[["B_extra"]]), c(ncol(D), 0L))
  expect_identical(zero[["fe_score_max"]], 0)
  expect_identical(dbreg:::assemble_fe_moments(NULL), NULL)
  expect_identical(dbreg:::decode_fe_moments(wide[wide[["__set"]] != 1L, ], K, 2L, 1L), NULL)
  duplicate = rbind(wide, wide[marginal[1L], , drop = FALSE])
  expect_identical(dbreg:::decode_fe_moments(duplicate, K, 2L, 1L), NULL)
  unmatched = wide
  unmatched[["__j"]][which(unmatched[["__set"]] > K)[1L]] = "absent"
  expect_identical(dbreg:::decode_fe_moments(unmatched, K, 2L, 1L), NULL)
  bad_width = packed
  attr(bad_width, "fe_values")[["values"]] = attr(bad_width, "fe_values")[["values"]][, 1:2, drop = FALSE]
  expect_identical(dbreg:::decode_fe_moments(bad_width, K, 2L, 1L), NULL)
  old = options(dbreg.projection_max_bytes = 1)
  on.exit(options(old))
  expect_equal(dbreg:::decode_fe_moments(packed, K, 2L, 1L), compact)
  expect_identical(dbreg:::assemble_fe_moments(compact), NULL)
  expect_error(dbreg:::read_numeric_json("[1,]"), "Invalid numeric SQL payload")
})
