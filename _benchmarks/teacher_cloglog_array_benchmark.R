# Follow-up on the exact deterministic teacher sample from the main benchmark.
# Rscript _benchmarks/teacher_cloglog_array_benchmark.R /tmp/teacher-cloglog-5m /tmp/teacher-cloglog-centered 2
# Optional fourth argument 'include_uncentered' reproduces the numerically
# stalled array-only experiment; the default tests the centered equivalent.
pkgload::load_all()
source('_benchmarks/teacher_cloglog_data.R')
source('_benchmarks/teacher_cloglog_profile.R')
source('_benchmarks/glm_databricks_moments.R')
source('_benchmarks/teacher_cloglog_arrays.R')
args = commandArgs(TRUE)
input = args[[1]]; out = args[[2]]
repetitions = if(length(args)>2L) as.integer(args[[3]]) else 2L
reference = readRDS(file.path(input,'fixest_strict_1.rds'))
summary = read.csv(file.path(input,'data_summary.csv'))
dir.create(file.path(out,'sql'),recursive=TRUE,showWarnings=FALSE)
log = new.env(parent=emptyenv())
log[['out']] = out; log[['id']] = 0L; log[['tag']] = paste0('teacher_arrays_',Sys.getpid())
log[['run']] = 'setup'; log[['queries']] = list(); log[['fits']] = list()
run = function() {
  native = DBI::dbConnect(odbc::databricks(),
    driver=Sys.getenv('DATABRICKS_ODBC_DRIVER','/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so'),
    httpPath=Sys.getenv('DATABRICKS_SQL_PATH'),useNativeQuery=TRUE,DefaultStringColumnLength=65536)
  on.exit(DBI::dbDisconnect(native))
  remote = methods::new('teacher_profile_connection',native=native,log=log)
  DBI::dbExecute(native,'SET use_cached_result=false')
  teacher_array_check(native,out)
  table = paste0('teacher_arrays_input_',Sys.getpid())
  sql = teacher_cloglog_sql(summary[['n']])
  stopifnot(identical(sql,paste(readLines(file.path(input,'data_generation.sql')),collapse='\n')))
  DBI::dbExecute(native,paste0('CREATE TEMP TABLE ',table,' AS SELECT *,log_outside_salary-10.8 AS centered_salary FROM (',sql,') q'))
  check = DBI::dbGetQuery(native,paste('SELECT COUNT(*) AS n,SUM(wastage) AS exits FROM',table))
  stopifnot(as.numeric(check[['n']])==summary[['n']],as.numeric(check[['exits']])==summary[['exits']])
  saved_moments = dbreg:::dbglm_fe_moments_sql
  saved_query = dbreg:::dbglm_get_query
  saved_fit = dbreg:::dbglm_fe_fit
  on.exit(assignInNamespace('dbglm_fe_moments_sql',saved_moments,'dbreg'),add=TRUE)
  on.exit(assignInNamespace('dbglm_get_query',saved_query,'dbreg'),add=TRUE)
  on.exit(assignInNamespace('dbglm_fe_fit',saved_fit,'dbreg'),add=TRUE)
  assignInNamespace('dbglm_fe_moments_sql',packed_moments_sql,'dbreg')
  assignInNamespace('dbglm_get_query',function(conn,sql) unpack_global_moments(saved_query(conn,sql)),'dbreg')
  assignInNamespace('dbglm_fe_fit',teacher_array_fit(saved_fit),'dbreg')
  fml = wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
  environment(fml) = baseenv()
  tasks = c(if(length(args)>3L && args[[4]]=='include_uncentered') 'arrays_packed',
    'packed_centered','arrays_centered',rep('arrays_split_centered',repetitions))
  for(index in seq_along(tasks)) {
    variant = tasks[[index]]
    rep = sum(tasks[seq_len(index)]==variant)
    model = if(grepl('centered',variant,fixed=TRUE)) {
      wastage ~ centered_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
    } else fml
    environment(model) = baseenv()
    assignInNamespace('dbglm_fe_fit',if(startsWith(variant,'arrays')) teacher_array_fit(saved_fit) else saved_fit,'dbreg')
    assignInNamespace('dbglm_fe_moments_sql',if(variant=='arrays_split_centered') teacher_split_packed_sql else packed_moments_sql,'dbreg')
    log[['run']] = paste0('Databricks_',variant,'_',rep)
    log[['error']] = ''; log[['warnings']] = character()
    fit = NULL
    gc(); DBI::dbGetQuery(native,'SELECT 1 AS warmup')
    cat('START',log[['run']],format(Sys.time()),'\n'); flush.console()
    timing = system.time(tryCatch(withCallingHandlers({
      fit = dbglm(model,conn=remote,table=table,family=binomial('cloglog'),vcov=~district,verbose=TRUE)
    },warning=function(w) log[['warnings']] = c(log[['warnings']],conditionMessage(w))),
      error=function(e) log[['error']] = conditionMessage(e)))
    scalar = function(name,default=NA_real_) if(!is.null(fit) && length(fit[[name]])==1L) fit[[name]] else default
    log[['fits']][[index]] = data.frame(engine='Databricks',variant=variant,repetition=rep,
      elapsed=timing[['elapsed']],user_cpu=timing[['user.self']],system_cpu=timing[['sys.self']],
      coefficient=if(is.null(fit)) NA else unname(coef(fit)[[1]]),se=if(is.null(fit)) NA else sqrt(vcov(fit)[1,1]),
      coef_error=if(is.null(fit)) NA else max(abs(coef(fit)-reference[['coefficients']])),
      se_error=if(is.null(fit)) NA else max(abs(sqrt(diag(vcov(fit)))-sqrt(diag(reference[['vcov']])))),
      nobs=if(is.null(fit)) NA else nobs(fit),iterations=scalar('iterations'),score=scalar('score_max'),halvings=scalar('step_halvings'),
      method=scalar('projection_method',''),deviance=if(is.null(fit)) NA else deviance(fit),
      error=log[['error']],warnings=paste(log[['warnings']],collapse='; '))
    write.csv(do.call(rbind,log[['fits']]),file.path(out,'fits.csv'),row.names=FALSE)
    if(!is.null(fit)) saveRDS(fit,file.path(out,paste0(log[['run']],'.rds')))
    print(log[['fits']][[index]]); flush.console()
  }
  DBI::dbExecute(native,paste('DROP TEMP TABLE',table))
}
run()
