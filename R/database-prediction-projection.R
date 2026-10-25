# SQL projections recover learned effects and regressor relations without
# collecting FE levels. Reusing native table names bounds the number of live
# intermediates; the fit manager owns all allocations until final adoption.
prediction_two_fe_projection = function(conn, source, fe, targets, create, drop, tol = 1e-10, max_iter = 100L) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  names = paste0("r", seq_along(targets))
  scale_row = backend_get_query(conn, paste("SELECT", paste(sprintf("MAX(ABS(%s)) AS s%d", targets, seq_along(targets)), collapse = ","), "FROM", source))
  scales = pmax(1, as.numeric(scale_row[1L, ]))
  residual = create(paste("SELECT", paste(c(fe, "w", sprintf("(%s)/(%s) AS %s", targets, dbglm_number(scales), names)), collapse = ","), "FROM", source))
  dictionaries = lapply(fe, function(f) create(paste("SELECT", f, "AS id,", paste0("CAST(0.0 AS DOUBLE) AS ", names, collapse = ","), "FROM", source, "GROUP BY", f)))
  dictionary_spare = rep(list(NULL), 2L)
  means = rep(list(NULL), 2L)
  spare_residual = NULL
  error = Inf
  for (iteration in seq_len(max_iter)) {
    for (i in seq_along(fe)) {
      f = fe[[i]]
      next_mean = create(paste("SELECT", f, "AS id,", paste(sprintf("FSUM(w*%s)/FSUM(w) AS %s", names, names), collapse = ","),
        "FROM", residual, "GROUP BY", f), target = means[[i]])
      next_dictionary = create(paste("SELECT d.id,", paste(sprintf("d.%s+m.%s AS %s", names, names, names), collapse = ","),
        "FROM", dictionaries[[i]], "d JOIN", next_mean, "m ON d.id=m.id"), target = dictionary_spare[[i]])
      next_residual = create(paste("SELECT", paste(c(paste0("s.", fe), "s.w",
        sprintf("s.%s-m.%s AS %s", names, names, names)), collapse = ","), "FROM", residual,
        "s JOIN", next_mean, "m ON", paste0("s.", f, "=m.id")), target = spare_residual)
      old_residual = residual
      residual = next_residual
      spare_residual = old_residual
      old_dictionary = dictionaries[[i]]
      dictionaries[[i]] = next_dictionary
      dictionary_spare[[i]] = old_dictionary
      if (!native) {
        drop(old_residual)
        drop(old_dictionary)
        if (!is.null(means[[i]])) drop(means[[i]])
      }
      means[[i]] = next_mean
    }
    terms = paste(sprintf("ABS(FSUM(w*%s)/FSUM(w)) AS %s", names, names), collapse = ",")
    branches = vapply(fe, function(f) paste("SELECT", terms, "FROM", residual, "GROUP BY", f), character(1))
    diagnostic = backend_get_query(conn, paste0("SELECT ", paste(sprintf("MAX(%s) AS %s", names, names), collapse = ","),
      " FROM (", paste(branches, collapse = " UNION ALL "), ") errors"))
    error = max(as.numeric(diagnostic[1L, ]))
    if (is.finite(error) && error <= tol) break
  }
  if (!is.finite(error) || error > tol) stop("Two-FE prediction recovery did not converge.")
  list(dictionaries = dictionaries, residual = residual, scales = scales,
    columns = names, iterations = iteration, error = error)
}

# Component labels are always IDs from the first FE. Pointer compression
# shortens propagation while every update stays inside its graph component.
prediction_two_fe_components = function(conn, source, fe, create, drop, max_iter = 100L) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  edges = create(paste("SELECT DISTINCT", fe[[1L]], "AS a,", fe[[2L]], "AS b FROM", source))
  left = create(paste("SELECT a AS id,a AS component FROM", edges, "GROUP BY a"))
  right = create(paste("SELECT b AS id,MIN(a) AS component FROM", edges, "GROUP BY b"))
  spare_left = spare_right = spare_compressed = NULL
  for (iteration in seq_len(max_iter)) {
    next_left = create(paste("SELECT e.a AS id,MIN(r.component) AS component FROM", edges,
      "e JOIN", right, "r ON e.b=r.id GROUP BY e.a"), target = spare_left)
    compressed = create(paste("SELECT a.id,b.component FROM", next_left, "a JOIN", next_left, "b ON a.component=b.id"), target = spare_compressed)
    next_right = create(paste("SELECT e.b AS id,MIN(l.component) AS component FROM", edges,
      "e JOIN", compressed, "l ON e.a=l.id GROUP BY e.b"), target = spare_right)
    old_left = left
    old_right = right
    left = compressed
    right = next_right
    spare_left = old_left
    spare_right = old_right
    spare_compressed = next_left
    if (!native) {
      drop(old_left)
      drop(old_right)
      drop(next_left)
    }
    mismatch = backend_get_query(conn, paste("SELECT 1 AS mismatch FROM", edges, "e JOIN", left,
      "l ON e.a=l.id JOIN", right, "r ON e.b=r.id WHERE l.component<>r.component LIMIT 1"))
    if (!nrow(mismatch)) break
  }
  if (nrow(mismatch)) stop("Two-FE component labels did not converge.")
  list(left = left, right = right, iterations = iteration)
}
