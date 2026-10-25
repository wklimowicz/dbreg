pkgload::load_all()
local({
 root = "_benchmarks/prediction_overflow_20261025"
 set.seed(932)
 n = 100000L
 d = as.data.frame(matrix(rnorm(n*5),n,5))
 names(d) = paste0("x",1:5)
 d[["y"]] = d[["x1"]]-.2*d[["x2"]]+rnorm(n)
 fit = dbglm(y~x1+x2+x3+x4+x5,data=d,family=gaussian())
 on.exit(close(fit))
 original = get("prediction_matrix_product",asNamespace("dbreg"))
 on.exit(assignInNamespace("prediction_matrix_product",original,"dbreg"),add=TRUE)
 records = list()
 for (mode in c("matrix_product","guarded")) {
  assignInNamespace("prediction_matrix_product",if(mode=="guarded") original else function(x,beta) as.vector(x%*%beta),"dbreg")
  value = predict(fit,d)
  if(mode=="matrix_product") reference = value else stopifnot(identical(value,reference))
  elapsed = replicate(5,system.time({for(i in 1:5) value=predict(fit,d)})[["elapsed"]]/5)
  records[[mode]] = data.frame(mode=mode,rows=n,columns=5,median_seconds=median(elapsed),min_seconds=min(elapsed),max_seconds=max(elapsed))
 }
 write.csv(do.call(rbind,records),file.path(root,"performance.csv"),row.names=FALSE)
 sources = c("DESCRIPTION","NAMESPACE",list.files("R",full.names=TRUE),file.path(root,"check-performance.R"))
 write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),file.path(root,"performance-source-md5.csv"),row.names=FALSE)
 writeLines(c(paste("UTC:",format(Sys.time(),tz="UTC")),capture.output(sessionInfo())),file.path(root,"performance-environment.txt"))
 print(do.call(rbind,records))
})
