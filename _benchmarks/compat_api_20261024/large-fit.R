# Rscript THIS DATABASE_FILE none|two baseline|lean|full
args = commandArgs(TRUE)
stopifnot(length(args)==3L,file.exists(args[[1L]]),args[[2L]] %in% c('none','two'),args[[3L]] %in% c('baseline','lean','full'))
pkgload::load_all()
dimension = args[[2L]]
mode = args[[3L]]
root = paste0('_benchmarks/compat_api_20261024/',dimension,'-',mode)
conn = DBI::dbConnect(duckdb::duckdb(config=list(threads='2',memory_limit='2GB')),dbdir=args[[1L]])
original = sort(DBI::dbListTables(conn))
fml = if(dimension=='two') y~x+z|g+h else y~x+z
strategy = if(dimension=='two') 'demean' else 'moments'
sources = c('DESCRIPTION','NAMESPACE',list.files('R',full.names=TRUE))
write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),paste0(root,'-source-md5.csv'),row.names=FALSE)
writeLines(c(paste('UTC:',format(Sys.time(),tz='UTC')),capture.output(sessionInfo())),paste0(root,'-environment.txt'))
reference = dbreg(fml,conn=conn,table='panel',strategy=strategy,row_id='id')
reference_ll = as.numeric(logLik(reference))
metrics = new.env(parent=emptyenv())
metrics[['active']] = FALSE
trace('dbGetQuery',where=asNamespace('DBI'),print=FALSE,exit=quote({
if(isTRUE(metrics[['active']])) {
 value=returnValue()
 metrics[['rows']]=max(metrics[['rows']],nrow(value))
 metrics[['bytes']]=metrics[['bytes']]+as.numeric(object.size(value))
}
}))
capture = get('retain_score_sample',asNamespace('dbreg'))
measured = function(...) {
 metrics[['active']]=TRUE
 on.exit({metrics[['active']]=FALSE})
 capture(...)
}
assignInNamespace('retain_score_sample',measured,'dbreg')
records=values=list()
for(i in 0:3) {
 gc()
 metrics[['rows']]=metrics[['bytes']]=0
 elapsed=system.time({fit=if(mode=='baseline') dbreg(fml,conn=conn,table='panel',strategy=strategy,row_id='id') else
  dbreg::feols(fml,conn=conn,table='panel',strategy=strategy,row_id='id',lean=mode=='lean',notes=FALSE)})[['elapsed']]
 stopifnot(identical(names(coef(fit)),names(coef(reference))),max(abs(coef(fit)-coef(reference)))<1e-8,
  max(abs(vcov(fit)-vcov(reference)))<1e-10,identical(nobs(fit),nobs(reference)),
  abs(as.numeric(logLik(fit))-reference_ll)<1e-4,metrics[['rows']]<=1L)
 if(mode=='lean') stopifnot(is.null(fit[['resource_state']]),is.null(fit[['fitted_state']][['prediction']]),
  identical(sort(DBI::dbListTables(conn)),original))
 if(mode=='full') {
  state=fit[['fitted_state']][['inference']]
  stopifnot(state[['nobs']]==1e6,state[['nobs_input']]==1e6)
  if(dimension=='two') {
   q=function(x) as.character(DBI::dbQuoteIdentifier(conn,x))
   check=DBI::dbGetQuery(conn,paste0('SELECT MAX(ABS(s.',q(state[['eta']]),'-(2+.3*p.x-.2*p.z+SIN(p.g*.05)+COS(p.h*.03)))) AS error FROM ',
    q(state[['table']]),' s JOIN panel p ON s.id=p.id WHERE s.',q(state[['estimated']])))
   stopifnot(check[['error']]<1e-7,identical(fit[['fitted_state']][['prediction']][['storage']],'database'))
  }
 }
 if(i>0L) {
 records[[i]]=data.frame(repetition=i,dimension=dimension,mode=mode,elapsed=elapsed,
  capture_max_rows=metrics[['rows']],capture_object_bytes=metrics[['bytes']],model_object_bytes=as.numeric(object.size(fit)))
 values[[i]]=list(coef=coef(fit),V=vcov(fit),N=nobs(fit),LL=as.numeric(logLik(fit)))
 }
 close(fit)
 stopifnot(DBI::dbIsValid(conn),identical(sort(DBI::dbListTables(conn)),original))
 message(dimension,' ',mode,' repetition ',i,': ',elapsed,'s')
}
untrace('dbGetQuery',where=asNamespace('DBI'))
DBI::dbDisconnect(conn,shutdown=TRUE)
write.csv(do.call(rbind,records),paste0(root,'.csv'),row.names=FALSE)
saveRDS(values,paste0(root,'-values.rds'))
