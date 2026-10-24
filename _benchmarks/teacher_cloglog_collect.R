# Consolidate the initial run, verified local settings, and SQL experiments.
# Rscript _benchmarks/teacher_cloglog_collect.R MAIN LOCAL CENTERED [INTERRUPTED]
args = commandArgs(TRUE)
folders = args[1:3]
fits = lapply(seq_along(folders),function(i) {
  d = read.csv(file.path(folders[[i]],'fits.csv'),colClasses=c(error='character',warnings='character'))
  d[['batch']] = c('main','verified_local','array_experiments')[[i]]
  if(i==1L) d[['variant']][d[['engine']]=='DuckDB'] = 'initial_default_config'
  d[['local_workers']] = ifelse(d[['engine']]=='fixest',4,
    ifelse(d[['engine']]=='DuckDB',if(i==2L) 4 else 16,NA))
  d
})
fits = do.call(rbind,fits)
write.csv(fits,'_benchmarks/teacher_cloglog_fits.csv',row.names=FALSE)
ok = subset(fits,error=='')
groups = split(ok,interaction(ok[['engine']],ok[['variant']],drop=TRUE))
summaries = lapply(groups,function(d) {
  if(!nrow(d)) return(NULL)
  data.frame(engine=d[['engine']][[1]],variant=d[['variant']][[1]],repetitions=nrow(d),
    median_seconds=median(d[['elapsed']]),min_seconds=min(d[['elapsed']]),max_seconds=max(d[['elapsed']]),
    min_iterations=min(d[['iterations']]),max_iterations=max(d[['iterations']]),
    max_coef_error=max(d[['coef_error']]),max_se_error=max(d[['se_error']]))
})
summary = do.call(rbind,summaries)
rownames(summary) = NULL
write.csv(summary,'_benchmarks/teacher_cloglog_summary.csv',row.names=FALSE)
print(summary,row.names=FALSE)
query_folders = c(folders[c(1,3)],if(length(args)>3L) args[[4]])
queries = do.call(rbind,lapply(query_folders,function(folder) read.csv(file.path(folder,'server_queries.csv'))))
write.csv(queries,'_benchmarks/teacher_cloglog_queries.csv',row.names=FALSE)
for(name in c('data_summary','fe_sizes','fe_pairs','cells','pure_groups','preparation')) {
  file.copy(file.path(folders[[1]],paste0(name,'.csv')),paste0('_benchmarks/teacher_cloglog_',name,'.csv'),overwrite=TRUE)
}
