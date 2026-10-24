# Validate fit equivalence before interpreting warmed performance measurements.
args = commandArgs(TRUE)
if (length(args) != 1L) stop("Expected the warmed-results directory.")
out = args[[1L]]
cases = c("baseline-default", "baseline-integer", "current-default", "current-integer",
  "current-row_id", "current-sample", "current-fraction")
rechecks = c("recheck-baseline-1", "recheck-current-1", "recheck-current-2", "recheck-baseline-2")
if (all(file.exists(file.path(out, paste0(rechecks, ".csv"))))) cases = c(cases, rechecks)
data = lapply(cases, function(case) read.csv(file.path(out, paste0(case, ".csv")), stringsAsFactors = FALSE))
names(data) = cases
decode = function(x) jsonlite::fromJSON(x)
checks = list()
for (case in cases) {
  frame = data[[case]]
  models = if (case == "current-fraction") c("glm", "glm_fe") else c("ols", "ols_fe", "glm", "glm_fe")
  stopifnot(setequal(unique(frame[["model"]]), models), nrow(frame) == 5L * length(models))
  reference_case = switch(case,
    "current-default" = "baseline-default", "current-integer" = "baseline-integer",
    "current-row_id" = "current-default", "current-sample" = "current-default",
    "current-fraction" = "current-integer", case)
  if (startsWith(case, "recheck-")) reference_case = "baseline-integer"
  scale = if (case == "current-fraction") .25 else 1
  for (i in seq_len(nrow(frame))) {
    row = frame[i, ]
    reference = data[[reference_case]]
    reference = reference[reference[["model"]] == row[["model"]], ][1L, ]
    b = decode(row[["coef_values"]])
    b_ref = decode(reference[["coef_values"]])
    V = decode(row[["vcov_values"]]) * scale
    V_ref = decode(reference[["vcov_values"]])
    same_names = identical(decode(row[["coef_names"]]), decode(reference[["coef_names"]]))
    close_b = length(b) == length(b_ref) && all(is.finite(b)) &&
      all(abs(b - b_ref) <= 1e-7 + 1e-6 * abs(b_ref))
    close_V = length(V) == length(V_ref) && all(is.finite(V)) &&
      all(abs(V - V_ref) <= 1e-8 + 1e-5 * abs(V_ref))
    close_ll = TRUE
    ll_delta = NA_real_
    if (startsWith(row[["model"]], "glm")) {
      expected_ll = scale * reference[["weighted_loglik"]]
      ll_delta = abs(row[["weighted_loglik"]] - expected_ll)
      close_ll = is.finite(ll_delta) && ll_delta <= 1e-8 + 1e-8 * abs(expected_ll)
    }
    passed = same_names && close_b && close_V && close_ll &&
      row[["nobs"]] == reference[["nobs"]] && row[["strategy"]] == reference[["strategy"]]
    checks[[length(checks) + 1L]] = data.frame(case = case, model = row[["model"]],
      repetition = row[["repetition"]], coefficient_delta = max(abs(b - b_ref)),
      covariance_delta = max(abs(V - V_ref)), weighted_likelihood_delta = ll_delta, passed = passed)
  }
}
checks = do.call(rbind, checks)
write.csv(checks, file.path(out, "value-checks.csv"), row.names = FALSE)
if (any(!checks[["passed"]])) stop("Fit equivalence failed; do not interpret performance ratios.")
medians = do.call(rbind, lapply(cases, function(case) {
  frame = data[[case]]
  do.call(rbind, lapply(unique(frame[["model"]]), function(model) {
    selected = frame[frame[["model"]] == model, ]
    data.frame(case = case, model = model, seconds = median(selected[["elapsed_seconds"]]),
      min_seconds = min(selected[["elapsed_seconds"]]), max_seconds = max(selected[["elapsed_seconds"]]),
      likelihood_seconds = median(selected[["likelihood_seconds"]]),
      iterations = if (startsWith(model, "glm")) median(selected[["iterations"]]) else NA_real_)
  }))
}))
write.csv(medians, file.path(out, "medians.csv"), row.names = FALSE)
pairs = list(default = c("baseline-default", "current-default"),
  integer = c("baseline-integer", "current-integer"),
  row_id = c("current-default", "current-row_id"),
  sample = c("current-default", "current-sample"),
  fraction = c("current-integer", "current-fraction"))
comparison = do.call(rbind, lapply(names(pairs), function(name) {
  pair = pairs[[name]]
  target = medians[medians[["case"]] == pair[[2L]], ]
  reference = medians[medians[["case"]] == pair[[1L]], ]
  reference = reference[match(target[["model"]], reference[["model"]]), ]
  data.frame(comparison = name, model = target[["model"]],
    reference_seconds = reference[["seconds"]], current_seconds = target[["seconds"]],
    ratio = target[["seconds"]] / reference[["seconds"]],
    likelihood_pass_seconds = target[["likelihood_seconds"]],
    reference_iterations = reference[["iterations"]], current_iterations = target[["iterations"]])
}))
if (all(rechecks %in% cases)) {
  pooled = lapply(c("baseline", "current"), function(label) {
    selected = c(paste0(label, "-integer"), rechecks[startsWith(rechecks, paste0("recheck-", label))])
    do.call(rbind, data[selected])
  })
  extra = do.call(rbind, lapply(unique(pooled[[1L]][["model"]]), function(model) {
    ref = pooled[[1L]][pooled[[1L]][["model"]] == model, ]
    cur = pooled[[2L]][pooled[[2L]][["model"]] == model, ]
    data.frame(comparison = "integer_pooled", model = model,
      reference_seconds = median(ref[["elapsed_seconds"]]), current_seconds = median(cur[["elapsed_seconds"]]),
      ratio = median(cur[["elapsed_seconds"]]) / median(ref[["elapsed_seconds"]]),
      likelihood_pass_seconds = 0,
      reference_iterations = if (startsWith(model, "glm")) median(ref[["iterations"]]) else NA_real_,
      current_iterations = if (startsWith(model, "glm")) median(cur[["iterations"]]) else NA_real_)
  }))
  comparison = rbind(comparison, extra)
}
write.csv(comparison, file.path(out, "comparison.csv"), row.names = FALSE)
print(comparison, digits = 3, row.names = FALSE)
cat("\nValidated", nrow(checks), "completed fits before timing comparisons.\n")
