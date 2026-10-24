# Independent, full-indicator Fisher-scoring check for the synthetic teacher
# benchmark. This is deliberately separate from dbglm's SQL/FE projections.
# Supports this complete, unweighted binary design and one cluster only.
# Rscript _benchmarks/teacher_cloglog_sparse_reference.R INPUT OUT [N]
teacher_sparse_reference = function(data, formula, cluster, covariance_df = NULL,
                                    initial_probability = NULL) {
  started = proc.time()[["elapsed"]]
  y = data[["y"]]
  stopifnot(all(y %in% c(0, 1)), !anyNA(data), length(unique(y)) == 2L)
  X = Matrix::sparse.model.matrix(formula, data)
  slopes = grep(":", colnames(X), fixed = TRUE)
  stopifnot(length(slopes) > 0L, nrow(X) == nrow(data))
  coefficients = setNames(numeric(ncol(X)), colnames(X))
  if (is.null(initial_probability)) {
    coefficients[["(Intercept)"]] = binomial("cloglog")[["linkfun"]](mean(y))
  } else {
    stopifnot(length(initial_probability) == nrow(X),
      all(initial_probability > 0 & initial_probability < 1))
    # The known generating predictor supplies only initialization, never a fit
    # from dbglm. Solve in the full indicator design without FE absorption.
    initial_eta = binomial("cloglog")[["linkfun"]](initial_probability)
    coefficients[] = as.numeric(Matrix::solve(Matrix::crossprod(X), Matrix::crossprod(X, initial_eta)))
  }
  likelihood = function(eta) {
    hazard = exp(eta)
    sum(ifelse(y == 1, log(-expm1(-hazard)), -hazard))
  }
  converged = FALSE
  for (iteration in seq_len(50L)) {
    eta = as.numeric(X %*% coefficients)
    hazard = exp(eta)
    ratio = hazard / expm1(hazard)
    weights = hazard * ratio
    scores = ifelse(y == 1, ratio, -hazard)
    information = Matrix::crossprod(Matrix::Diagonal(x = sqrt(weights)) %*% X)
    gradient = as.numeric(Matrix::crossprod(X, scores))
    increment = as.numeric(Matrix::solve(information, gradient))
    predictor_step = as.numeric(X %*% increment)
    score = max(abs(gradient) / sqrt(Matrix::diag(information)))
    slope_step = max(abs(increment[slopes]))
    message("Sparse reference iteration ", iteration, ": score=", signif(score, 4),
      "; slope step=", signif(slope_step, 4))
    if (score < 1e-8 && slope_step < 1e-10 && max(abs(predictor_step)) < 1e-8) {
      converged = TRUE
      break
    }
    before = likelihood(eta)
    allowance = 64 * .Machine[["double.eps"]] * max(1, abs(before))
    accepted = FALSE
    for (halving in 0:20) {
      fraction = 2^(-halving)
      after = likelihood(eta + fraction * predictor_step)
      if (is.finite(after) && after >= before - allowance) {
        coefficients = coefficients + fraction * increment
        accepted = TRUE
        break
      }
    }
    if (!accepted) stop("Independent sparse reference line search failed.")
  }
  if (!converged) stop("Independent sparse reference failed to converge.")
  selector = Matrix::sparseMatrix(i = slopes, j = seq_along(slopes), x = 1,
    dims = c(ncol(X), length(slopes)))
  bread = Matrix::solve(information, selector)
  groups = as.integer(factor(data[[cluster]]))
  stopifnot(max(groups) > 1L)
  cluster_scores = Matrix::sparseMatrix(i = groups, j = seq_len(nrow(X)), x = scores,
    dims = c(max(groups), nrow(X))) %*% X
  K = if (is.null(covariance_df)) ncol(X) else covariance_df
  stopifnot(K > 0, K < nrow(X))
  covariance = as.matrix(Matrix::crossprod(cluster_scores %*% bread)) *
    (nrow(X) - 1) / (nrow(X) - K) * max(groups) / (max(groups) - 1)
  dimnames(covariance) = list(names(coefficients)[slopes], names(coefficients)[slopes])
  stopifnot(all(is.finite(coefficients)), all(is.finite(covariance)))
  list(coef = coefficients[slopes], vcov = covariance, logLik = likelihood(eta),
    nobs = nrow(X), n_input = nrow(data), keys = as.numeric(data[["id"]]), exits = sum(y),
    model_df = ncol(X), covariance_df = K, score = score, slope_step = slope_step,
    elapsed = proc.time()[["elapsed"]] - started, iterations = iteration,
    reference_engine = "independent full-indicator sparse Fisher scoring")
}

if (sys.nframe() == 0L) {
  source("_benchmarks/teacher_cloglog_data.R")
  args = commandArgs(TRUE)
  input = args[[1]]
  out = args[[2]]
  n = if (length(args) > 2L) as.integer(args[[3]]) else 5000000L
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")))
  data = DBI::dbGetQuery(conn, paste0("SELECT *,wastage AS y,",
    teacher_cloglog_salary_sql(), " FROM read_parquet(",
    DBI::dbQuoteString(conn, file.path(input, "teachers.parquet")), ") WHERE id < ", n, " ORDER BY id"))
  DBI::dbDisconnect(conn, shutdown = TRUE)
  stopifnot(identical(as.numeric(data[["id"]]), as.numeric(seq_len(n) - 1L)))
  for (name in c("log_outside_salary", "log_teacher_salary")) data[[name]] = data[[name]] - mean(data[[name]])
  formula = y ~ log_outside_salary:factor(Experience) + log_teacher_salary:factor(Experience) +
    factor(spine_band_year) + factor(district) + factor(Experience) + factor(AGE) + factor(SEX) + factor(FTPTWK)
  fes = c("spine_band_year", "district", "Experience", "AGE", "SEX", "FTPTWK")
  # Only district is nested in district clustering in this fixture.
  for (fe in setdiff(fes, "district")) {
    stopifnot(any(vapply(split(data[["district"]], data[[fe]]), function(x) length(unique(x)) > 1L, logical(1))))
  }
  full_df = 88L + sum(vapply(data[fes], function(x) length(unique(x)), integer(1))) - 5L
  covariance_df = full_df - length(unique(data[["district"]])) + 1L
  reference = teacher_sparse_reference(data, formula, "district", covariance_df,
    initial_probability = data[["true_probability"]])
  stopifnot(reference[["model_df"]] == full_df, reference[["nobs"]] == n)
  saveRDS(reference, file.path(out, "sparse_reference.rds"))
  writeLines(c(capture.output(sessionInfo()), paste("Rows:", n)), file.path(out, "sparse_reference_environment.txt"))
  print(reference[c("nobs", "exits", "model_df", "covariance_df", "score", "slope_step", "elapsed")])
}
