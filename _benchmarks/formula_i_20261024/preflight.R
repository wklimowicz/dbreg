pkgload::load_all()
root = "_benchmarks/formula_i_20261024"
out = file.path(root,"preflight")
dir.create(out,showWarnings=FALSE)
sources = c("DESCRIPTION","NAMESPACE",list.files("R",full.names=TRUE),
 list.files("inst/tinytest",full.names=TRUE),file.path(root,c("preflight.R","test-preflight.R")))
write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),file.path(out,"source-md5.csv"),row.names=FALSE)
writeLines(c(paste("UTC:",format(Sys.time(),tz="UTC")),capture.output(sessionInfo())),file.path(out,"environment.txt"))
started = Sys.time()
r = tinytest::run_test_file(file.path(root,"test-preflight.R"))
print(r)
results = as.data.frame(r)
record = data.frame(assertions=nrow(results),failures=sum(is.na(results[["result"]]) | !results[["result"]]),
 elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")))
write.csv(record,file.path(out,"summary.csv"),row.names=FALSE)
saveRDS(r,file.path(out,"results.rds"))
stopifnot(record[["assertions"]]==4L,record[["failures"]]==0L)
