root='_benchmarks/compat_api_20261024'
final=c('local-final-boolean/full-source-md5.csv','native-final/source-md5.csv','native-ssc/source-md5.csv')
rows=list()
for(file in final) {
 manifest=read.csv(file.path(root,file))
 differences=unname(tools::md5sum(manifest[['file']]))!=manifest[['md5']]
 stopifnot(!anyNA(differences),!any(differences))
 rows[[length(rows)+1L]]=data.frame(manifest=file,files=nrow(manifest),differences=0,scope='final implementation')
}
performance=paste0(rep(c('none','two'),each=3),'-',rep(c('baseline','lean','full'),2),'-source-md5.csv')
for(file in c('local/full-source-md5.csv','local-final/full-source-md5.csv','native/source-md5.csv',performance)) {
 manifest=read.csv(file.path(root,file))
 differences=unname(tools::md5sum(manifest[['file']]))!=manifest[['md5']]
 allowed=if(file %in% performance) c('R/compat-extractors.R','R/sql-literals.R') else
  if(startsWith(file,'local-final/')) c('R/sql-literals.R','inst/tinytest/helper_compat_api.R') else
   c('R/compat-extractors.R','R/sql-literals.R','inst/tinytest/helper_compat_api.R','inst/tinytest/test_compat_api.R')
 stopifnot(!anyNA(differences),setequal(manifest[['file']][differences],allowed))
 archived=file.path(root,'pre-ssc-source',manifest[['file']][differences])
 stopifnot(identical(unname(tools::md5sum(archived)),manifest[['md5']][differences]))
 rows[[length(rows)+1L]]=data.frame(manifest=file,files=nrow(manifest),differences=sum(differences),
  scope=if(file %in% performance) 'before SSC tag/Boolean literal patches; numeric IID fixture paths unchanged' else 'superseded development gate; exact original files archived')
}
audit=do.call(rbind,rows)
write.csv(audit,file.path(root,'source-audit.csv'),row.names=FALSE)
print(audit,row.names=FALSE)
