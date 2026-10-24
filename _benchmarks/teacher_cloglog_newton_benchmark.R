# Interleave production Fisher scoring and an adaptive Newton/Fisher prototype.
# Rscript _benchmarks/teacher_cloglog_newton_benchmark.R /tmp/teacher-cloglog-5m /tmp/teacher-cloglog-newton
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("_benchmarks/teacher_cloglog_profile.R")
source("_benchmarks/glm_databricks_more_helpers.R")
args = commandArgs(TRUE); input = args[[1]]; out = args[[2]]
dir.create(file.path(out, "sql"), recursive = TRUE, showWarnings = FALSE)
reference = readRDS(file.path(input, "fixest_strict_1.rds"))
summary = read.csv(file.path(input, "data_summary.csv"))
log = new.env(parent = emptyenv())
log[["out"]] = out; log[["id"]] = 0L; log[["tag"]] = paste0("teacher_newton_", Sys.getpid())
log[["run"]] = "setup"; log[["queries"]] = list()
run = function() {
  native = DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
  on.exit(DBI::dbDisconnect(native))
  conn = methods::new("teacher_profile_connection", native = native, log = log)
  DBI::dbExecute(native, "SET use_cached_result=false")
  table = paste0("teacher_newton_input_", Sys.getpid())
  sql = teacher_cloglog_sql(summary[["n"]])
  stopifnot(identical(sql, paste(readLines(file.path(input, "data_generation.sql")), collapse = "\n")))
  DBI::dbExecute(native, paste("CREATE TEMP TABLE", table, "AS", sql))
  counts = DBI::dbGetQuery(native, paste("SELECT COUNT(*) AS n,SUM(wastage) AS exits FROM", table))
  stopifnot(as.numeric(counts[["n"]])==summary[["n"]], as.numeric(counts[["exits"]])==summary[["exits"]])
  formula = wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
  environment(formula) = baseenv()
  results = list()
  for (variant in c("control", "adaptive_newton", "adaptive_newton", "control")) {
    experiment = if(variant=="adaptive_newton") glm_experiment_newton("adaptive") else NULL
    rep = sum(vapply(results, function(r) r[["variant"]]==variant, logical(1)))+1L
    log[["run"]] = paste0(variant,"_",rep)
    DBI::dbGetQuery(native,"SELECT 1 AS warmup"); gc()
    cat("START",log[["run"]],format(Sys.time()),"\n"); flush.console()
    elapsed = system.time({fit=tryCatch(dbglm(formula,conn=conn,table=table,
      family=binomial("cloglog"),vcov=~district,verbose=TRUE),
      finally=if(!is.null(experiment)) experiment[["restore"]]())})[["elapsed"]]
    if(!is.null(experiment)) stopifnot(experiment[["state"]][["last_fisher"]])
    row = data.frame(variant=variant,repetition=rep,elapsed=elapsed,
      coefficient=unname(coef(fit)[[1]]),se=sqrt(vcov(fit)[1,1]),
      coef_error=max(abs(coef(fit)-reference[["coefficients"]])),
      se_error=max(abs(sqrt(diag(vcov(fit)))-sqrt(diag(reference[["vcov"]])))),
      nobs=nobs(fit),iterations=fit[["iterations"]],score=fit[["score_max"]],
      halvings=fit[["step_halvings"]],pruning_rounds=fit[["pruning_rounds"]],deviance=deviance(fit))
    results[[length(results)+1L]]=row
    write.csv(do.call(rbind,results),file.path(out,"fits.csv"),row.names=FALSE)
    saveRDS(fit,file.path(out,paste0(log[["run"]],".rds")))
    print(row);flush.console()
    stopifnot(nobs(fit)==summary[["n"]],row[["coef_error"]]<1e-8,row[["se_error"]]<1e-7,
      abs(deviance(fit)-reference[["deviance"]])<1e-5,fit[["score_max"]]<=fit[["control"]][["score_tol"]])
  }
  DBI::dbExecute(native,paste("DROP TEMP TABLE",table))
}
run()
