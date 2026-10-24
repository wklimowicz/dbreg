# Run from the repository root after all six benchmark batches complete.
directory = "_benchmarks/linear_prediction_20261024/performance"
batches = c("baseline-default", "current-default", "current-prediction",
  "recheck-baseline", "recheck-default", "recheck-prediction")
configurations = rep(c("baseline", "current", "prediction"), 2L)
results = lapply(seq_along(batches), function(i) {
  x = read.csv(file.path(directory, paste0(batches[[i]], ".csv")))
  stopifnot(nrow(x) == 15L, all(is.finite(x[["elapsed_seconds"]])))
  x[["batch"]] = batches[[i]]
  x[["configuration"]] = configurations[[i]]
  x
})
fits = do.call(rbind, results)
decode = function(x) jsonlite::fromJSON(x)
close_values = function(x, y, atol, rtol) {
  length(x) == length(y) && all(is.finite(x)) && all(abs(x - y) <= atol + rtol * abs(y))
}
for (model in unique(fits[["model"]])) {
  reference = fits[fits[["model"]] == model & fits[["configuration"]] == "baseline", ][1L, ]
  candidates = fits[fits[["model"]] == model, ]
  for (i in seq_len(nrow(candidates))) {
    candidate = candidates[i, ]
    stopifnot(identical(decode(candidate[["coef_names"]]), decode(reference[["coef_names"]])),
      candidate[["nobs"]] == reference[["nobs"]],
      candidate[["strategy"]] == reference[["strategy"]],
      close_values(decode(candidate[["coef_values"]]), decode(reference[["coef_values"]]), 1e-7, 1e-6),
      close_values(decode(candidate[["vcov_values"]]), decode(reference[["vcov_values"]]), 1e-10, 1e-5))
  }
}
summary = aggregate(elapsed_seconds ~ model + configuration, fits, median)
summary[["repetitions"]] = 10L
summary[["state_bytes"]] = vapply(seq_len(nrow(summary)), function(i) {
  max(fits[["prediction_state_bytes"]][fits[["model"]] == summary[["model"]][[i]] &
    fits[["configuration"]] == summary[["configuration"]][[i]]])
}, numeric(1))
write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
memory = data.frame(batch = batches, configuration = configurations,
  peak_rss_kib = vapply(batches, function(batch) {
    lines = readLines(file.path(directory, paste0(batch, "-time.txt")))
    value = lines[grepl("Maximum resident set size", lines, fixed = TRUE)]
    stopifnot(length(value) == 1L)
    as.numeric(sub(".*: *", "", value))
  }, numeric(1)))
write.csv(memory, file.path(directory, "memory.csv"), row.names = FALSE)
cat("All 90 measured fits match baseline names, coefficients, covariance, N and strategy.\n")
print(summary, row.names = FALSE)
print(memory, row.names = FALSE)
