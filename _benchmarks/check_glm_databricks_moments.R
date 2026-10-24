# Numerical portability checks for the experimental query generators.
# These are not warehouse performance measurements.
pkgload::load_all()
source('_benchmarks/glm_databricks_moments.R')
check = function() {
  original = dbreg:::dbglm_fe_moments_sql
  on.exit(assignInNamespace('dbglm_fe_moments_sql',original,'dbreg'))
  for(generator in list(grouped_moments_sql,hybrid_moments_sql)) {
    assignInNamespace('dbglm_fe_moments_sql',generator,'dbreg')
    result = tinytest::run_test_file('inst/tinytest/test_glm_databricks.R')
    print(result)
    stopifnot(all(as.data.frame(result)[['result']]))
  }
  # The packed global row is the only JSON value in a result; pair rows retain
  # numeric columns. Check double precision and the NULL rows of the adapter.
  d = data.frame(set=c(0L,1L),h=c(NA,2),globals=c(
    '{"dev":1.2345678901234567,"score":1.25e-200,"bad":0}',NA))
  names(d) = c('__set','__h','__globals')
  unpacked = unpack_global_moments(d)
  stopifnot(identical(unpacked[['dev']],c(1.2345678901234567,NA_real_)),
    identical(unpacked[['score']],c(1.25e-200,NA_real_)),
    identical(unpacked[['bad']],c(0,NA_real_)))
}
check()
