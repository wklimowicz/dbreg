# Run from the package root with configured DATABRICKS_* credentials:
# Rscript _benchmarks/glm_databricks_performance.R /tmp/glm-perf 100000 4 union,grouping,hybrid
# Uses synthetic data and session-scoped tables only. Connection/startup, data
# generation/download, and the fixest oracle are excluded from fit timings.
# SQL result reuse is disabled; data and warehouse caches can still be warm.
# Variants replace a private function in this R process, not package files.
pkgload::load_all()
library(DBI)
source('_benchmarks/glm_databricks_moments.R')
args = commandArgs(TRUE)
out = if(length(args)) args[[1]] else '/tmp/dbglm-profile'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
state = new.env(parent=emptyenv())
state[['rows']] = list()
state[['run']] = 'setup'
state[['tag']] = paste0('dbglm_perf_', Sys.getpid())
state[['id']] = 0L
methods::setClass('glm_profile_connection', contains='DBIConnection', slots=c(native='ANY'))
methods::setMethod('dbGetInfo','glm_profile_connection',function(dbObj,...) DBI::dbGetInfo(dbObj@native,...))
methods::setMethod('dbIsValid','glm_profile_connection',function(dbObj,...) DBI::dbIsValid(dbObj@native,...))
methods::setMethod('dbQuoteIdentifier',c('glm_profile_connection','character'),function(conn,x,...) DBI::dbQuoteIdentifier(conn@native,x,...))
methods::setMethod('dbQuoteString',c('glm_profile_connection','character'),function(conn,x,...) DBI::dbQuoteString(conn@native,x,...))
methods::setMethod('dbExistsTable',c('glm_profile_connection','character'),function(conn,name,...) DBI::dbExistsTable(conn@native,name,...))
profile_query = function(conn, statement, execute, ...) {
  state[['id']] = state[['id']] + 1L
  id = state[['id']]
  calls = vapply(sys.calls(), function(x) as.character(x[[1]])[[1]], character(1))
  phase = if('dbglm_fe_evaluate' %in% calls) 'evaluation' else if('dbglm_fe_covariance' %in% calls) 'covariance' else 'setup_cleanup'
  tag = paste(state[['tag']], state[['run']], id, sep='_')
  writeLines(statement, file.path(out, paste0(tag,'.sql')))
  start = proc.time()[['elapsed']]
  result = if(execute) DBI::dbExecute(conn@native,paste0('/* ',tag,' */ ',statement),...) else DBI::dbGetQuery(conn@native,paste0('/* ',tag,' */ ',statement),...)
  elapsed = proc.time()[['elapsed']] - start
  state[['rows']][[length(state[['rows']])+1L]] = data.frame(tag=tag, run=state[['run']], phase=phase,
    elapsed=elapsed, sql_bytes=nchar(statement,type='bytes'), rows=if(execute) NA else nrow(result),
    result_bytes=as.numeric(object.size(result)))
  write.csv(do.call(rbind,state[['rows']]),file.path(out,'queries.csv'),row.names=FALSE)
  cat(tag,phase,round(elapsed,3),'s\n'); flush.console()
  result
}
methods::setMethod('dbGetQuery',c('glm_profile_connection','character'),function(conn,statement,...) profile_query(conn,statement,FALSE,...))
methods::setMethod('dbExecute',c('glm_profile_connection','character'),function(conn,statement,...) profile_query(conn,statement,TRUE,...))
run = function() {
  native = DBI::dbConnect(odbc::databricks(),driver=Sys.getenv('DATABRICKS_ODBC_DRIVER','/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so'),
    httpPath=Sys.getenv('DATABRICKS_SQL_PATH'),useNativeQuery=TRUE,DefaultStringColumnLength=65536)
  on.exit(DBI::dbDisconnect(native))
  conn = methods::new('glm_profile_connection',native=native)
  DBI::dbExecute(native,'SET use_cached_result=false')
  print(DBI::dbGetQuery(native,'SELECT current_version() AS version'))
  saveRDS(list(date=Sys.time(),session=sessionInfo(),controls=dbglm_control(),
    result_cache=FALSE),file.path(out,'metadata.rds'))
  orig = dbreg:::dbglm_fe_moments_sql
  on.exit(assignInNamespace('dbglm_fe_moments_sql',orig,'dbreg'),add=TRUE)
  original_query = dbreg:::dbglm_get_query
  assignInNamespace('dbglm_get_query',function(conn,sql) unpack_global_moments(original_query(conn,sql)),'dbreg')
  on.exit(assignInNamespace('dbglm_get_query',original_query,'dbreg'),add=TRUE)
  results = list()
  sizes = if(length(args)>1L) as.integer(args[[2]]) else 100000L
  K = if(length(args)>2L) as.integer(args[[3]]) else 4L
  stopifnot(K %in% c(4L,6L), length(sizes)==1L, is.finite(sizes), sizes>0L)
  for(n in sizes) {
    table = paste0('dbglm_perf_input_',Sys.getpid())
    fs = if(K==4L) c(100,20,5,3) else c(2001,296,44,44,2,2)
    factors = paste(sprintf('CAST(pmod(id*%d+floor(id/%d),%d)+1 AS INT) AS f%d',c(17,13,7,11,19,23)[seq_len(K)],c(101,37,23,17,31,43)[seq_len(K)],fs,seq_len(K)),collapse=',')
    DBI::dbExecute(native,paste0('CREATE TEMP TABLE ',table,' AS WITH raw AS (SELECT id,',factors,
      ',(sin(id*12.9898)*43758.5453-floor(sin(id*12.9898)*43758.5453)-.5)*4 AS x,',
      '(sin(id*78.233)*19341.167-floor(sin(id*78.233)*19341.167)-.5)*4 AS z,',
      'sin(id*19.913)*7193.31-floor(sin(id*19.913)*7193.31) AS v FROM range(',n,')),',
      'pred AS (SELECT *, -.5+.4*x-.2*z+.3*sin(f1)+.2*cos(f2) AS eta FROM raw) ',
      'SELECT *, CAST(v < 1-exp(-exp(eta)) AS INT) AS y FROM pred'))
    d = DBI::dbGetQuery(native,paste('SELECT * FROM',table))
    fml = as.formula(paste('y~x+z|',paste0('f',seq_len(K),collapse='+')))
    ref = fixest::feglm(fml,d,family=binomial('cloglog'),vcov=~f2,fixef.rm='none',glm.tol=1e-12,glm.iter=100)
    saveRDS(list(coef=coef(ref),vcov=vcov(ref),iterations=ref[['iterations']]),file.path(out,'reference.rds'))
    variants = if(length(args)>3L) strsplit(args[[4]],',',fixed=TRUE)[[1]] else c('union','grouping','grouping','union')
    for(variant in variants) {
      state[['run']] = paste0(variant,'_',length(results)+1L)
      assignInNamespace('dbglm_fe_moments_sql',switch(variant,union=orig,grouping=grouped_moments_sql,hybrid=hybrid_moments_sql,packed=packed_moments_sql),'dbreg')
      fit = NULL
      state[['error']] = ''
      gc()
      elapsed = system.time(tryCatch({fit=dbglm(fml,conn=conn,table=table,family=binomial('cloglog'),vcov=~f2,verbose=TRUE)},error=function(e){state[['error']] = conditionMessage(e)}))[['elapsed']]
      results[[length(results)+1L]] = data.frame(run=state[['run']],n=n,fe=K,elapsed=elapsed,
        iterations=if(is.null(fit)) NA else fit[['iterations']],halvings=if(is.null(fit)) NA else fit[['step_halvings']],
        score=if(is.null(fit)) NA else fit[['score_max']],
        coef_error=if(is.null(fit)) NA else max(abs(coef(fit)-coef(ref))),
        se_error=if(is.null(fit)) NA else max(abs(sqrt(diag(vcov(fit)))-sqrt(diag(vcov(ref))))),error=state[['error']])
      write.csv(do.call(rbind,results),file.path(out,'fits.csv'),row.names=FALSE)
      if(!is.null(fit)) saveRDS(fit,file.path(out,paste0(state[['run']],'.rds')))
      print(utils::tail(do.call(rbind,results),1)); flush.console()
    }
    DBI::dbExecute(native,paste('DROP TEMP TABLE',table))
  }
}
run()
