# Run from the repository root after both default/prediction batches finish.
out = "_benchmarks/prediction_state_20261024/performance"
cases = c("current-default", "current-prediction", "recheck-default", "recheck-prediction")
data = lapply(cases, function(case) read.csv(file.path(out, paste0(case, ".csv"))))
names(data) = cases
reference = data[["current-default"]]
checks = list()
for (case in cases) {
  frame = data[[case]]
  models = if (grepl("prediction", case)) c("glm", "glm_fe") else c("ols", "ols_fe", "glm", "glm_fe")
  stopifnot(nrow(frame) == 5L * length(models), setequal(unique(frame[["model"]]), models))
  for (i in seq_len(nrow(frame))) {
    row = frame[i, ]
    ref = reference[reference[["model"]] == row[["model"]], ][1L, ]
    b = jsonlite::fromJSON(row[["coef_values"]])
    br = jsonlite::fromJSON(ref[["coef_values"]])
    v = jsonlite::fromJSON(row[["vcov_values"]])
    vr = jsonlite::fromJSON(ref[["vcov_values"]])
    passed = identical(jsonlite::fromJSON(row[["coef_names"]]), jsonlite::fromJSON(ref[["coef_names"]])) &&
      all(abs(b - br) <= 1e-7 + 1e-6 * abs(br)) && all(abs(v - vr) <= 1e-8 + 1e-5 * abs(vr)) &&
      row[["nobs"]] == ref[["nobs"]] && row[["strategy"]] == ref[["strategy"]]
    if (startsWith(row[["model"]], "glm")) {
      passed = passed && abs(row[["weighted_loglik"]] - ref[["weighted_loglik"]]) <= 1e-8 * (1 + abs(ref[["weighted_loglik"]]))
    }
    checks[[length(checks) + 1L]] = data.frame(case = case, model = row[["model"]], repetition = row[["repetition"]], passed = passed)
  }
}
checks = do.call(rbind, checks)
write.csv(checks, file.path(out, "value-checks.csv"), row.names = FALSE)
stopifnot(all(checks[["passed"]]))
default = do.call(rbind, data[c("current-default", "recheck-default")])
prediction = do.call(rbind, data[c("current-prediction", "recheck-prediction")])
comparison = do.call(rbind, lapply(c("glm", "glm_fe"), function(model) {
  ref = default[default[["model"]] == model, ]
  cur = prediction[prediction[["model"]] == model, ]
  data.frame(model = model, default_seconds = median(ref[["elapsed_seconds"]]),
    prediction_seconds = median(cur[["elapsed_seconds"]]), ratio = median(cur[["elapsed_seconds"]]) / median(ref[["elapsed_seconds"]]),
    prediction_min_seconds = min(cur[["elapsed_seconds"]]), prediction_max_seconds = max(cur[["elapsed_seconds"]]),
    state_bytes = median(cur[["prediction_state_bytes"]]),
    max_recovery_error = if (all(is.na(cur[["prediction_recovery_error"]]))) NA_real_ else max(cur[["prediction_recovery_error"]], na.rm = TRUE))
}))
write.csv(comparison, file.path(out, "comparison.csv"), row.names = FALSE)
print(comparison)
cat(nrow(checks), "measured fits passed value checks.\n")
