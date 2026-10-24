local({
  X = matrix(c(1, 2), ncol = 1L, dimnames = list(NULL, "b1"))
  cells = data.frame(cl = c("A", "B"), a = c("a\rb", "a"), b = c("c", "b\rc"),
    sum_w_gc = 1, sum_wy_gc = 1)
  compressed = cells[c("a", "b")]
  meat = dbreg:::cluster_meat_from_cells(cells, compressed, c("a", "b"), "cl", X, c(0, 0))
  expect_equal(as.numeric(meat), 5)
  expect_equal(attr(meat, "n_clusters"), 2)
  # Close numeric FE values keep their fitted raw keys behind display labels.
  cells = data.frame(cl = c(1, 1 + 2^-52), a = c(1, 1 + 2^-52), sum_w_gc = 1, sum_wy_gc = 1)
  compressed = data.frame(a = dbreg:::exact_numeric_factor(cells[["a"]]))
  meat = dbreg:::cluster_meat_from_cells(cells, compressed, "a", "cl", X, c(0, 0))
  expect_equal(as.numeric(meat), 5)
  expect_equal(attr(meat, "n_clusters"), 2)
  # Preserve the cell-meat helper's existing explicit missing-cluster group.
  cells[["cl"]] = c(NA_real_, 1)
  meat = dbreg:::cluster_meat_from_cells(cells, compressed, "a", "cl", X, c(0, 0))
  expect_equal(as.numeric(meat), 5)
  expect_equal(attr(meat, "n_clusters"), 2)
  cells[["a"]][2] = 99
  expect_error(dbreg:::cluster_meat_from_cells(cells, compressed, "a", "cl", X, c(0, 0)), "fitted compression sample")
})
