# Run from the package root. Remote generation is reproducible; only synthetic
# session-scoped tables are created. Identical downloaded doubles feed DuckDB
# and fixest. Fit timers exclude preparation, transfer, startup and diagnostics.
# Rscript _benchmarks/teacher_cloglog_benchmark.R /tmp/teacher-cloglog-5m 5000000 2
pkgload::load_all()
source('_benchmarks/teacher_cloglog_data.R')
source('_benchmarks/teacher_cloglog_profile.R')
source('_benchmarks/glm_databricks_moments.R')
args = commandArgs(TRUE)
out = if(length(args)) args[[1]] else '/tmp/teacher-cloglog-5m'
n = if(length(args)>1L) as.integer(args[[2]]) else 5000000L
repetitions = if(length(args)>2L) as.integer(args[[3]]) else 2L
dir.create(file.path(out,'sql'),recursive=TRUE,showWarnings=FALSE)
log = new.env(parent=emptyenv())
log[['out']] = out; log[['id']] = 0L; log[['tag']] = paste0('teacher_perf_',Sys.getpid())
log[['run']] = 'setup'; log[['queries']] = list(); log[['fits']] = list()

run = function() {
  native = DBI::dbConnect(odbc::databricks(),
    driver=Sys.getenv('DATABRICKS_ODBC_DRIVER','/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so'),
    httpPath=Sys.getenv('DATABRICKS_SQL_PATH'),useNativeQuery=TRUE,DefaultStringColumnLength=65536)
  on.exit(DBI::dbDisconnect(native))
  local = DBI::dbConnect(duckdb::duckdb(config=list(threads='4',memory_limit='4GB',temp_directory=file.path(out,'duckdb_tmp'))))
  on.exit(DBI::dbDisconnect(local,shutdown=TRUE),add=TRUE)
  local_settings = DBI::dbGetQuery(local,"SELECT current_setting('threads') AS threads,current_setting('memory_limit') AS memory_limit")
  stopifnot(as.integer(local_settings[['threads']])==4L)
  write.csv(local_settings,file.path(out,'duckdb_settings.csv'),row.names=FALSE)
  remote = methods::new('teacher_profile_connection',native=native,log=log)
  DBI::dbExecute(native,'SET use_cached_result=false')
  fixest::setFixest_nthreads(4)
  metadata = list(date=Sys.time(),n=n,repetitions=repetitions,session=sessionInfo(),
    duckdb_settings=local_settings,fixest_threads=fixest::getFixest_nthreads(),
    controls=dbglm_control(),remote_version=DBI::dbGetQuery(native,'SELECT current_version() AS version'))
  saveRDS(metadata,file.path(out,'metadata.rds'))
  table = paste0('teacher_bench_',Sys.getpid())
  sql = teacher_cloglog_sql(n)
  writeLines(sql,file.path(out,'data_generation.sql'))
  preparation = system.time(DBI::dbExecute(native,paste('CREATE TEMP TABLE',table,'AS',sql)))[['elapsed']]
  cat('Remote data generation:',preparation,'seconds\n'); flush.console()
  download = system.time({d = DBI::dbGetQuery(native,paste('SELECT * FROM',table))})[['elapsed']]
  cat('Download identical data:',download,'seconds; bytes:',as.numeric(object.size(d)),'\n'); flush.console()
  local_load = system.time(DBI::dbWriteTable(local,'teachers',d,temporary=TRUE))[['elapsed']]
  # Preserve the exact benchmark data locally for follow-up profiling, without
  # keeping a persistent remote table. Parquet export is outside fit timings.
  DBI::dbExecute(local,paste0("COPY teachers TO ",DBI::dbQuoteString(local,file.path(out,'teachers.parquet'))," (FORMAT PARQUET)"))
  write.csv(data.frame(generation=preparation,download=download,duckdb_load=local_load),file.path(out,'preparation.csv'),row.names=FALSE)
  fe = c('spine_band_year','district','Experience','AGE','SEX','FTPTWK')
  counts = vapply(d[fe],function(x) length(unique(x)),integer(1))
  summary = data.frame(n=nrow(d),exits=sum(d[['wastage']]),exit_rate=mean(d[['wastage']]),
    true_probability_min=min(d[['true_probability']]),true_probability_max=max(d[['true_probability']]),
    age_experience_cor=cor(d[['AGE']],d[['Experience']]),band_experience_cor=cor(d[['band']],d[['Experience']]),
    female_share=mean(d[['SEX']]),part_time_share=mean(d[['FTPTWK']]),
    part_time_female=mean(d[['FTPTWK']][d[['SEX']]==1]),part_time_male=mean(d[['FTPTWK']][d[['SEX']]==0]))
  write.csv(summary,file.path(out,'data_summary.csv'),row.names=FALSE)
  write.csv(data.frame(fe=names(counts),levels=counts),file.path(out,'fe_sizes.csv'),row.names=FALSE)
  print(summary); print(counts); flush.console()
  # Pair support and exact compressibility are computed once, outside timings.
  pair_counts = lapply(utils::combn(fe,2,simplify=FALSE),function(pair) {
    q = paste0('SELECT COUNT(*) AS n FROM (SELECT ',paste(pair,collapse=','),' FROM teachers GROUP BY ',paste(pair,collapse=','),') q')
    data.frame(first=pair[[1]],second=pair[[2]],edges=as.numeric(DBI::dbGetQuery(local,q)[['n']]))
  })
  write.csv(do.call(rbind,pair_counts),file.path(out,'fe_pairs.csv'),row.names=FALSE)
  cells = DBI::dbGetQuery(local,paste0('SELECT COUNT(*) AS cells FROM (SELECT ',paste(c(fe,'log_outside_salary'),collapse=','),
    ' FROM teachers GROUP BY ',paste(c(fe,'log_outside_salary'),collapse=','),') q'))
  write.csv(cells,file.path(out,'cells.csv'),row.names=FALSE)
  # Audit pure FE cells before deciding the shared estimation sample.
  pure = lapply(fe,function(f) DBI::dbGetQuery(local,paste0("SELECT '",f,"' AS fe,COUNT(*) AS pure_groups,COALESCE(SUM(n),0) AS rows FROM ",
    '(SELECT COUNT(*) AS n FROM teachers GROUP BY ',f,' HAVING MIN(wastage)=MAX(wastage)) q')))
  write.csv(do.call(rbind,pure),file.path(out,'pure_groups.csv'),row.names=FALSE)
  if(any(vapply(pure,function(x) as.numeric(x[['pure_groups']]),numeric(1))>0)) {
    stop('Synthetic sample contains pure FE groups; align recursive pruning before benchmarking.')
  }
  fml = wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
  environment(fml) = baseenv()
  original_moments = dbreg:::dbglm_fe_moments_sql
  original_query = dbreg:::dbglm_get_query
  on.exit(assignInNamespace('dbglm_fe_moments_sql',original_moments,'dbreg'),add=TRUE)
  on.exit(assignInNamespace('dbglm_get_query',original_query,'dbreg'),add=TRUE)
  assignInNamespace('dbglm_get_query',function(conn,sql) unpack_global_moments(original_query(conn,sql)),'dbreg')
  log[['reference']] = NULL
  record = function(engine,variant,rep,thunk) {
    log[['run']] = paste(engine,variant,rep,sep='_')
    log[['error']] = ''; log[['warnings']] = character()
    gc()
    fit = NULL
    cat('START',log[['run']],format(Sys.time()),'\n'); flush.console()
    timing = system.time(tryCatch(withCallingHandlers({fit = thunk()},warning=function(w) {
      log[['warnings']] = c(log[['warnings']],conditionMessage(w))
    }),error=function(e) {log[['error']] = conditionMessage(e)}))
    reference = log[['reference']]
    if(!is.null(fit) && is.null(reference)) {
      if(isFALSE(fit[['convStatus']])) stop('The independent reference did not converge.')
      reference = list(coef=coef(fit),vcov=vcov(fit))
      log[['reference']] = reference
    }
    scalar = function(name,default=NA_real_) if(!is.null(fit) && length(fit[[name]])==1L) fit[[name]] else default
    log[['fits']][[length(log[['fits']])+1L]] = data.frame(engine=engine,variant=variant,repetition=rep,
      elapsed=timing[['elapsed']],user_cpu=timing[['user.self']],system_cpu=timing[['sys.self']],
      coefficient=if(is.null(fit)) NA else unname(coef(fit)[[1]]),se=if(is.null(fit)) NA else sqrt(vcov(fit)[1,1]),
      coef_error=if(is.null(fit)) NA else max(abs(coef(fit)-reference[['coef']])),
      se_error=if(is.null(fit)) NA else max(abs(sqrt(diag(vcov(fit)))-sqrt(diag(reference[['vcov']])))),
      nobs=if(is.null(fit)) NA else nobs(fit),iterations=scalar('iterations'),score=scalar('score_max'),halvings=scalar('step_halvings'),
      method=scalar('projection_method',''),deviance=if(is.null(fit)) NA else deviance(fit),
      error=log[['error']],warnings=paste(log[['warnings']],collapse='; '))
    write.csv(do.call(rbind,log[['fits']]),file.path(out,'fits.csv'),row.names=FALSE)
    if(!is.null(fit)) saveRDS(if(engine=='fixest') list(coefficients=coef(fit),vcov=vcov(fit),iterations=fit[['iterations']],
      nobs=nobs(fit),deviance=deviance(fit),convStatus=fit[['convStatus']]) else fit,file.path(out,paste0(log[['run']],'.rds')))
    print(utils::tail(do.call(rbind,log[['fits']]),1)); flush.console()
    invisible(fit)
  }
  # Independent strict oracle is timed and reported separately from defaults.
  record('fixest','strict',1L,function() fixest::feglm(fml,d,family=binomial('cloglog'),vcov=~district,
    fixef.rm='none',glm.tol=1e-12,fixef.tol=1e-9,glm.iter=100,fixef.iter=10000,nthreads=4))
  if(is.null(log[['reference']])) stop('Independent reference failed.')
  for(rep in seq_len(repetitions)) {
    record('fixest','default',rep,function() fixest::feglm(fml,d,family=binomial('cloglog'),vcov=~district,fixef.rm='none',nthreads=4))
    assignInNamespace('dbglm_fe_moments_sql',original_moments,'dbreg')
    record('DuckDB','current',rep,function() dbglm(fml,conn=local,table='teachers',family=binomial('cloglog'),vcov=~district,verbose=TRUE))
    # Reverse current/optimized order for the second repetition.
    for(variant in if(rep %% 2L) c('current','packed') else c('packed','current')) {
      assignInNamespace('dbglm_fe_moments_sql',if(variant=='current') original_moments else packed_moments_sql,'dbreg')
      DBI::dbGetQuery(native,'SELECT 1 AS warmup')
      record('Databricks',variant,rep,function() dbglm(fml,conn=remote,table=table,family=binomial('cloglog'),vcov=~district,verbose=TRUE))
    }
  }
  DBI::dbExecute(native,paste('DROP TEMP TABLE',table))
}
run()
