# Synthetic teacher-year rows, not confidential data or a calibrated population.
# Uniform innovations use independently salted xxhash64 values on Databricks.
# The salary regressor is a district/year/age/sex cell value, as external wages
# commonly are. The DGP has exactly the fitted additive FE structure and beta=.8.
teacher_cloglog_sql = function(n = 5000000L) {
  uniform = function(seed) sprintf(
    "CAST(pmod(xxhash64(id, CAST(%d AS BIGINT)), 1000000007) AS DOUBLE)/1000000007.0", seed)
  innovations = sprintf("%s AS u%d", vapply(1:8, uniform, character(1)), 1:8)
  paste0(
    "WITH innovations AS (SELECT id,", paste(innovations, collapse=","), " FROM range(", n, ")), ",
    "location AS (SELECT *, CAST(FLOOR(23*u1)+1 AS INT) AS year, ",
    "CAST(FLOOR(296*POWER(u2,1.3))+1 AS INT) AS district FROM innovations), ",
    "career AS (SELECT *, CAST(FLOOR(44*POWER(GREATEST(0.0,LEAST(0.999999999,",
    "u3 + .12*SIN(district*.17))),1.4)) AS INT) AS Experience, ",
    "CAST(u4 < .72+.04*SIN(district*.11) AS INT) AS SEX FROM location), ",
    "demographics AS (SELECT *, CAST(LEAST(64,21+Experience+FLOOR(12*u5)+",
    "CASE WHEN u6<.15 THEN FLOOR(12*u6/.15) ELSE 0 END) AS INT) AS AGE FROM career), ",
    "employment AS (SELECT *, CAST(u7 < 1/(1+EXP(-(-2.1+1.05*SEX+.025*(AGE-35)+",
    ".2*SIN(district*.11)))) AS INT) AS FTPTWK, ",
    "CAST(1+FLOOR(GREATEST(0.0,LEAST(86.999999999,2*Experience+",
    "6*SIN(district*.17)+16*(u6-.5)-2*SEX))) AS INT) AS band FROM demographics), ",
    "covariates AS (SELECT *, CAST((year-1)*87+band AS INT) AS spine_band_year, ",
    "10.8+.08*SIN(district*.09)+.006*(AGE-40)+.012*(year-12)-.04*SEX+",
    ".08*SIN(.23*district+.51*year)+.05*COS(.17*AGE+.13*year)+",
    ".035*SEX*SIN(.37*year+.13*district) AS log_outside_salary FROM employment), ",
    "predictor AS (SELECT *, -4.0+.8*(log_outside_salary-10.8)+",
    ".20*SIN(.13*band)+.10*COS(.3*year)+.18*COS(.1*district)+",
    ".25*EXP(-Experience/4.0)-.005*Experience+",
    ".65*POWER(GREATEST(AGE-50,0)/14.0,2)-.10*SEX+.20*FTPTWK AS eta FROM covariates) ",
    "SELECT id,year,band,spine_band_year,district,Experience,AGE,SEX,FTPTWK,",
    "CAST(log_outside_salary AS DOUBLE) AS log_outside_salary,",
    "CAST(u8 < 1-EXP(-EXP(eta)) AS INT) AS wastage, ",
    "CAST(1-EXP(-EXP(eta)) AS DOUBLE) AS true_probability FROM predictor"
  )
}

# Explicit DOUBLE operands avoid Spark's finite-scale DECIMAL division. The
# binary-fraction wage innovation must be identical in remote/local fixtures.
teacher_cloglog_salary_sql = function() {
  paste0("10.5+CAST(band AS DOUBLE)/128.0+CAST(Experience AS DOUBLE)/512.0+",
    "CAST((id % 97)-48 AS DOUBLE)/64.0 AS log_teacher_salary")
}
