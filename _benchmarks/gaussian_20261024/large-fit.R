# Rscript THIS existing-paired-noise-duckdb-file baseline|gaussian
args = commandArgs(TRUE)
stopifnot(length(args) == 2L, file.exists(args[[1L]]), args[[2L]] %in% c('baseline','gaussian'))
pkgload::load_all()
mode = args[[2L]]
root = paste0('_benchmarks/gaussian_20261024/',mode)
conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = '2', memory_limit = '2GB')), dbdir = args[[1L]])
original = sort(DBI::dbListTables(conn))
sources = c('DESCRIPTION','NAMESPACE',list.files('R',full.names=TRUE))
write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),paste0(root,'-source-md5.csv'),row.names=FALSE)
writeLines(c(paste('UTC:',format(Sys.time(),tz='UTC')),capture.output(sessionInfo())),paste0(root,'-environment.txt'))
metrics = new.env(parent=emptyenv())
metrics[['active']] = FALSE
trace('dbGetQuery',where=asNamespace('DBI'),print=FALSE,exit=quote({
if(isTRUE(metrics[['active']])) metrics[['rows']]=max(metrics[['rows']],nrow(returnValue()))
}))
reference = dbreg(y~x+z,conn=conn,table='panel',strategy='moments')
records = list()
for(i in 0:5) {
metrics[['rows']] = 0L
metrics[['active']] = TRUE
elapsed = system.time({fit=if(mode=='baseline') dbreg(y~x+z,conn=conn,table='panel',strategy='moments') else
 dbglm(y~x+z,conn=conn,table='panel',family=gaussian())})[['elapsed']]
metrics[['active']] = FALSE
adjustment = if(mode=='gaussian') (1e6-1)/(1e6-length(coef(fit))) else 1
stopifnot(identical(names(coef(fit)),names(coef(reference))),max(abs(coef(fit)-coef(reference)))<1e-8,
 max(abs(vcov(fit)-vcov(reference)*adjustment))<1e-10,identical(nobs(fit),nobs(reference)),metrics[['rows']]<=1L,
 identical(sort(DBI::dbListTables(conn)),original))
if(mode=='gaussian') stopifnot(fit[['iterations']]==1L)
if(i>0L) records[[i]]=data.frame(repetition=i,mode=mode,elapsed=elapsed,max_rows=metrics[['rows']],object_bytes=as.numeric(object.size(fit)))
}
untrace('dbGetQuery',where=asNamespace('DBI'))
write.csv(do.call(rbind,records),paste0(root,'.csv'),row.names=FALSE)
saveRDS(list(coef=coef(fit),V=vcov(fit),n=nobs(fit)),paste0(root,'-values.rds'))
DBI::dbDisconnect(conn,shutdown=TRUE)
