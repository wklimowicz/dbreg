directory = "_benchmarks/database_prediction_20261024/performance"
base = readRDS(file.path(directory, "none-values.rds"))
retained = readRDS(file.path(directory, "prediction-values.rds"))
for (value in c(base, retained)) {
  stopifnot(identical(names(value[["coef"]]), names(base[[1L]][["coef"]])),
    isTRUE(all.equal(value[["coef"]], base[[1L]][["coef"]], tolerance = 1e-7)),
    isTRUE(all.equal(value[["vcov"]], base[[1L]][["vcov"]], tolerance = 1e-7)),
    value[["nobs"]] == 1e6)
}
summary = do.call(rbind, lapply(c("none", "prediction"), function(mode) {
  fits = read.csv(file.path(directory, paste0(mode, ".csv")))
  timing = readLines(file.path(directory, paste0(mode, "-time.txt")))
  rss = as.numeric(sub(".*: *", "", timing[grepl("Maximum resident set size", timing, fixed = TRUE)]))
  data.frame(mode = mode, repetitions = nrow(fits), median_seconds = median(fits[["elapsed_seconds"]]),
    peak_rss_kib = rss, maximum_state_bytes = max(fits[["state_bytes"]]),
    maximum_prediction_fit_query_rows = max(fits[["maximum_prediction_fit_query_rows"]]),
    maximum_prediction_fit_query_bytes = max(fits[["prediction_fit_query_bytes"]]))
}))
write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
cat("All six measured fits match coefficients, full covariance and N.\n")
print(summary, row.names = FALSE)
