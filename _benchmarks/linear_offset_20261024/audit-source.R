root = '_benchmarks/linear_offset_20261024'
for(file in c('local/full-source-md5.csv','native/source-md5.csv','performance-source-md5.csv')) {
manifest=read.csv(file.path(root,file))
stopifnot(identical(unname(tools::md5sum(manifest[['file']])),manifest[['md5']]))
cat(file,':',nrow(manifest),'source files match\n')
}
