# Rscript THIS existing-paired-noise-duckdb-file
args = commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[[1L]]))
pkgload::load_all()
conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = '2', memory_limit = '2GB')), dbdir = args[[1L]])
root = '_benchmarks/linear_offset_20261024/'
sources = c('DESCRIPTION','NAMESPACE',list.files('R',full.names=TRUE))
write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),paste0(root,'performance-source-md5.csv'),row.names=FALSE)
writeLines(c(paste('UTC:',format(Sys.time(),tz='UTC')),capture.output(sessionInfo())),paste0(root,'performance-environment.txt'))
DBI::dbExecute(conn, 'CREATE TEMP VIEW offset_panel AS SELECT *,y+(.17*x-.09*z) AS shifted_y,.17*x-.09*z AS o FROM panel')
original = sort(DBI::dbListTables(conn))
metrics = new.env(parent=emptyenv());metrics[['active']]=FALSE
trace('dbGetQuery',where=asNamespace('DBI'),print=FALSE,exit=quote({
if(isTRUE(metrics[['active']])) metrics[['rows']]=max(metrics[['rows']],nrow(returnValue()))
}))
set.seed(562)
records = list(); values=list()
for(i in 0:5) for(mode in sample(c('baseline','offset'))) {
metrics[['rows']]=0L;metrics[['active']]=TRUE
elapsed=system.time({fit=if(mode=='baseline') dbreg(y~x+z,conn=conn,table='offset_panel',strategy='moments') else
 dbreg(shifted_y~x+z,conn=conn,table='offset_panel',strategy='moments',offset='o')})[['elapsed']]
metrics[['active']]=FALSE
value=list(coef=coef(fit),V=vcov(fit),n=nobs(fit))
if(is.null(values[['reference']])) values[['reference']]=value
reference=values[['reference']]
stopifnot(identical(names(value[['coef']]),names(reference[['coef']])),max(abs(value[['coef']]-reference[['coef']]))<1e-8,
 max(abs(value[['V']]-reference[['V']]))<1e-10,identical(value[['n']],reference[['n']]),metrics[['rows']]<=1L,
 identical(sort(DBI::dbListTables(conn)),original))
if(i>0L) records[[length(records)+1L]]=data.frame(repetition=i,mode=mode,elapsed=elapsed,max_rows=metrics[['rows']],object_bytes=as.numeric(object.size(fit)))
}
untrace('dbGetQuery',where=asNamespace('DBI'))
write.csv(do.call(rbind,records),paste0(root,'performance.csv'),row.names=FALSE)
saveRDS(values,paste0(root,'performance-values.rds'))
DBI::dbDisconnect(conn,shutdown=TRUE)
