# Exact classification checks shared by offline and native Databricks tests.
# End-to-end recursive pruning and retained-row parity are tested by the caller.
dbglm_test_databricks_pruning = function(conn) {
  stored = basename(tempfile("glm_prune_preflight_"))
  on.exit(try(DBI::dbExecute(conn, paste("DROP TEMP TABLE", stored)), silent = TRUE))
  d = expand.grid(f = c("O'Brien", "é / ^"), g = c(10, 20), y = c(0, 1), stringsAsFactors = FALSE)
  d[["id"]] = seq_len(nrow(d))
  DBI::dbExecute(conn, paste("CREATE TEMP TABLE", stored, "AS", dbreg:::dbglm_values_sql(conn, d)))
  f = as.character(DBI::dbQuoteIdentifier(conn, "f group"))
  # Delta's default physical column mapping forbids spaces; a query alias
  # exercises quoted identifier handling without changing storage properties.
  table = paste0("(SELECT f AS ", f, ",g,y,id FROM ", stored, ")")
  separate = list(f, "g")
  composite = list(c(f, "g"), f)
  check = function(source, groups, family, expected) {
    for (mode in names(expected)) expect_identical(
      dbreg:::dbglm_databricks_pruning_needed(conn, source, groups, "y", mode, family),
      expected[[mode]])
  }
  modes = c(singletons = FALSE, infinite_coef = FALSE, perfect_fit = FALSE)
  check(table, separate, "binomial", modes)
  check(table, composite, "binomial", modes)
  # A pure group is not a singleton. Repeated/shared grouping columns must
  # preserve the classification for compound FE keys.
  pure = paste0("(SELECT ", f, ",g,CASE WHEN g=10 THEN 0 ELSE y END AS y FROM ", table, ")")
  check(pure, separate, "binomial", c(singletons=FALSE, infinite_coef=TRUE, perfect_fit=TRUE))
  check(pure, composite, "binomial", c(singletons=FALSE, infinite_coef=TRUE, perfect_fit=TRUE))
  singleton = paste0("(SELECT * FROM ", table, " WHERE id != 5)")
  check(singleton, separate, "binomial", modes)
  check(singleton, composite, "binomial", c(singletons=TRUE, infinite_coef=TRUE, perfect_fit=TRUE))
  # Constant positive Poisson groups have finite effects; all-zero groups do not.
  positive = paste0("(SELECT ", f, ",g,2 AS y FROM ", table, ")")
  check(positive, composite, "poisson", modes)
  zero = paste0("(SELECT ", f, ",g,CASE WHEN g=10 THEN 0 ELSE 2 END AS y FROM ", table, ")")
  check(zero, composite, "poisson", c(singletons=FALSE, infinite_coef=TRUE, perfect_fit=TRUE))
  check(paste0("(SELECT * FROM ", table, " WHERE 1=0)"), composite, "binomial", modes)
  expect_false(dbreg:::dbglm_databricks_pruning_needed(NULL, table, composite, "y", "none", "binomial"))
  expect_false(dbreg:::dbglm_databricks_pruning_needed(NULL, table, list(), "y", "perfect_fit", "binomial"))
  expect_true(dbreg:::dbglm_databricks_pruning_needed(NULL, table, as.list(paste0("f",1:31)), "y", "perfect_fit", "binomial"))
}
