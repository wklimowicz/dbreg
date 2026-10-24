# Verify local settings at runtime and rerun on the exact exported data.
# Rscript _benchmarks/teacher_cloglog_local_benchmark.R /tmp/teacher-cloglog-5m /tmp/teacher-cloglog-local 2
pkgload::load_all()
args = commandArgs(TRUE)
input = args[[1]]; out = args[[2]]
repetitions = if(length(args)>2L) as.integer(args[[3]]) else 2L
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reference = readRDS(file.path(input,'fixest_strict_1.rds'))
run = function() {
  conn = DBI::dbConnect(duckdb::duckdb(config=list(threads='4',memory_limit='4GB',temp_directory=file.path(out,'duckdb_tmp'))))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  settings = DBI::dbGetQuery(conn,"SELECT current_setting('threads') AS threads,current_setting('memory_limit') AS memory_limit")
  stopifnot(as.integer(settings[['threads']])==4L)
  write.csv(settings,file.path(out,'settings.csv'),row.names=FALSE)
  print(settings)
  load = system.time(DBI::dbExecute(conn,paste0('CREATE TEMP TABLE teachers AS SELECT * FROM read_parquet(',
    DBI::dbQuoteString(conn,file.path(input,'teachers.parquet')),')')))[['elapsed']]
  write.csv(data.frame(load_seconds=load),file.path(out,'load.csv'),row.names=FALSE)
  fml = wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
  environment(fml) = baseenv()
  results = list()
  for(rep in seq_len(repetitions)) {
    gc()
    cat('START DuckDB 4 threads',rep,format(Sys.time()),'\n'); flush.console()
    timing = system.time({fit = dbglm(fml,conn=conn,table='teachers',family=binomial('cloglog'),vcov=~district,verbose=TRUE)})
    results[[rep]] = data.frame(engine='DuckDB',variant='current_4threads',repetition=rep,
      elapsed=timing[['elapsed']],user_cpu=timing[['user.self']],system_cpu=timing[['sys.self']],
      coefficient=unname(coef(fit)[[1]]),se=sqrt(vcov(fit)[1,1]),
      coef_error=max(abs(coef(fit)-reference[['coefficients']])),
      se_error=max(abs(sqrt(diag(vcov(fit)))-sqrt(diag(reference[['vcov']])))),
      nobs=nobs(fit),iterations=fit[['iterations']],score=fit[['score_max']],halvings=fit[['step_halvings']],
      method=fit[['projection_method']],deviance=deviance(fit),error='',warnings='')
    write.csv(do.call(rbind,results),file.path(out,'fits.csv'),row.names=FALSE)
    saveRDS(fit,file.path(out,paste0('DuckDB_4threads_',rep,'.rds')))
    print(results[[rep]]); flush.console()
  }
}
run()
