# Local SQL portability fixture. This does not emulate distributed aggregation
# order, Spark planning, ODBC transport or warehouse cleanup; use live tests too.
dbglm_test_databricks_sql = function(statement, native) {
 # Spark's JSON reader is emulated with R in this SQL fixture; actual Spark
 # decoding, extreme doubles and key escaping are checked remotely as well.
 while(grepl(" FROM explode(from_json(",statement,fixed=TRUE)) {
  m_pos=regexpr("SELECT CAST(item[",statement,fixed=TRUE)[1]
  sub_stmt=substr(statement,m_pos,nchar(statement))
  end_pos=regexpr("AS lookup(item)",sub_stmt,fixed=TRUE)[1]+nchar("AS lookup(item)")-1L
  block=substr(sub_stmt,1,end_pos)
  encoded=regmatches(block,regexec("unbase64\\('([^']*)'\\)",block,perl=TRUE))[[1]][2]
  rows=jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(encoded)),simplifyVector=FALSE)
  fields=regmatches(block,gregexpr('CAST\\(item\\[[0-9]+\\] AS (?:DOUBLE|STRING)\\) AS (?:"[^\"]+"|[A-Za-z_][A-Za-z0-9_]*)',block,perl=TRUE))[[1]]
  frame=as.data.frame(lapply(seq_along(fields),function(j) {
   if(grepl(' AS DOUBLE)',fields[j],fixed=TRUE)) vapply(rows,function(r)as.numeric(r[[j]]),numeric(1)) else
     vapply(rows,function(r)as.character(r[[j]]),character(1))
  }),stringsAsFactors=FALSE)
  names(frame)=gsub('"','',sub('.* AS ','',fields),fixed=TRUE)
  replacement=dbreg:::dbglm_values_sql(native,frame,format='values')
  statement=paste0(substr(statement,1,m_pos-1),replacement,substr(sub_stmt,end_pos+1,nchar(sub_stmt)))
 }
 # The production backend owns these JSON/array functions; the local fixture
 # adapts only their syntax. The live suite tests the actual Spark decoder.
 pattern = "from_json\\(CAST\\(unbase64\\('([A-Za-z0-9+/=]+)'\\) AS STRING\\),'ARRAY<DOUBLE>'\\)"
 repeat {
  match = regexpr(pattern, statement, perl=TRUE)
  if (match[[1]] < 0) break
  block = regmatches(statement, match)
  encoded = regmatches(block, regexec(pattern, block, perl=TRUE))[[1]][[2]]
  values = jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(encoded)))
  regmatches(statement, match) = paste0('[',paste(dbreg:::dbglm_number(values),collapse=','),']::DOUBLE[]')
 }
 statement = gsub('element_at(', 'list_extract(', statement, fixed=TRUE)
 repeat {
  match = regexpr('to_json\\(named_struct\\([^()]*\\)\\)', statement, perl=TRUE)
  if (match[[1]] < 0) break
  block = regmatches(statement, match)
  fields = trimws(strsplit(sub('^to_json\\(named_struct\\((.*)\\)\\)$','\\1',block),',',fixed=TRUE)[[1]])
  keys = gsub("'",'',fields[seq(1,length(fields),2)],fixed=TRUE)
  values = fields[seq(2,length(fields),2)]
  # Avoid requiring DuckDB's optional JSON extension in ordinary offline CI.
  pieces = sprintf("'\"%s\":' || COALESCE(CAST(%s AS VARCHAR),'null')",keys,values)
  regmatches(statement,match) = paste0("('{' || ",paste(pieces,collapse=" || ',' || ")," || '}')")
 }
 statement=gsub('monotonically_increasing_id()','ROW_NUMBER() OVER ()',statement,fixed=TRUE)
 statement=gsub(' EXCEPT',' EXCLUDE',statement,fixed=TRUE)
 sub('^DROP TEMP TABLE','DROP TABLE',statement)
}

dbglm_test_databricks_connection = function(native) {
 if (!methods::isClass('dbglm_sql_fixture')) methods::setClass('dbglm_sql_fixture',contains='DBIConnection',slots=c(native='ANY'))
 methods::setMethod('dbGetInfo','dbglm_sql_fixture',function(dbObj,...)list(dbms.name='Databricks SQL'))
 methods::setMethod('dbIsValid','dbglm_sql_fixture',function(dbObj,...)DBI::dbIsValid(dbObj@native))
 methods::setMethod('dbQuoteIdentifier',c('dbglm_sql_fixture','character'),function(conn,x,...)DBI::dbQuoteIdentifier(conn@native,x,...))
 methods::setMethod('dbQuoteString',c('dbglm_sql_fixture','character'),function(conn,x,...)DBI::dbQuoteString(conn@native,x,...))
 methods::setMethod('dbExistsTable',c('dbglm_sql_fixture','character'),function(conn,name,...)DBI::dbExistsTable(conn@native,name,...))
 methods::setMethod('dbGetQuery',c('dbglm_sql_fixture','character'),function(conn,statement,...)DBI::dbGetQuery(conn@native,dbglm_test_databricks_sql(statement,conn@native),...))
 methods::setMethod('dbExecute',c('dbglm_sql_fixture','character'),function(conn,statement,...)DBI::dbExecute(conn@native,dbglm_test_databricks_sql(statement,conn@native),...))
 methods::new('dbglm_sql_fixture',native=native)
}
