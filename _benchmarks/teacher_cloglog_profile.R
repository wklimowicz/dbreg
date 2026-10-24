# Lightweight instrumentation for remote SQL only. This keeps native DuckDB
# connections intact, so registration and local materialization are unchanged.
library(DBI)
methods::setClass('teacher_profile_connection', contains='DBIConnection', slots=c(native='ANY', log='environment'))
methods::setMethod('dbGetInfo','teacher_profile_connection',function(dbObj,...) DBI::dbGetInfo(dbObj@native,...))
methods::setMethod('dbIsValid','teacher_profile_connection',function(dbObj,...) DBI::dbIsValid(dbObj@native,...))
methods::setMethod('dbQuoteIdentifier',c('teacher_profile_connection','character'),function(conn,x,...) DBI::dbQuoteIdentifier(conn@native,x,...))
methods::setMethod('dbQuoteString',c('teacher_profile_connection','character'),function(conn,x,...) DBI::dbQuoteString(conn@native,x,...))
methods::setMethod('dbExistsTable',c('teacher_profile_connection','character'),function(conn,name,...) DBI::dbExistsTable(conn@native,name,...))
teacher_profile_query = function(conn, statement, execute, ...) {
  log = conn@log
  log[['id']] = log[['id']]+1L
  calls = vapply(sys.calls(),function(x) as.character(x[[1]])[[1]],character(1))
  phase = if('dbglm_fe_evaluate' %in% calls) 'evaluation' else if('dbglm_fe_covariance' %in% calls) 'covariance' else 'setup_cleanup'
  tag = paste(log[['tag']],log[['run']],log[['id']],sep='_')
  writeLines(statement,file.path(log[['out']],'sql',paste0(tag,'.sql')))
  start = proc.time()[['elapsed']]
  sql = paste0('/* ',tag,' */ ',statement)
  result = if(execute) DBI::dbExecute(conn@native,sql,...) else DBI::dbGetQuery(conn@native,sql,...)
  log[['queries']][[length(log[['queries']])+1L]] = data.frame(tag=tag,run=log[['run']],phase=phase,
    elapsed=proc.time()[['elapsed']]-start,sql_bytes=nchar(statement,type='bytes'),
    rows=if(execute) NA else nrow(result),result_bytes=as.numeric(object.size(result)))
  write.csv(do.call(rbind,log[['queries']]),file.path(log[['out']],'queries.csv'),row.names=FALSE)
  cat(tag,phase,round(utils::tail(log[['queries']],1)[[1]][['elapsed']],3),'s\n'); flush.console()
  result
}
methods::setMethod('dbGetQuery',c('teacher_profile_connection','character'),function(conn,statement,...) teacher_profile_query(conn,statement,FALSE,...))
methods::setMethod('dbExecute',c('teacher_profile_connection','character'),function(conn,statement,...) teacher_profile_query(conn,statement,TRUE,...))
