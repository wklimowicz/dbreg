args=commandArgs(TRUE)
stopifnot(length(args)==1L,args[[1L]] %in% c('before','after'))
if(args[[1L]]=='before') loadNamespace('fixest')
pkgload::load_all()
if(args[[1L]]=='after') loadNamespace('fixest')
fit=dbreg::feols(mpg~hp+wt,mtcars,notes=FALSE,lean=TRUE)
reference=fixest::feols(mpg~hp+wt,mtcars,notes=FALSE)
stopifnot(!inherits(fit,'fixest'),max(abs(fixest::se(fit)-as.numeric(fixest::se(reference))))<1e-8,
 max(abs(as.numeric(fixest::coeftable(fit))-as.numeric(as.matrix(fixest::coeftable(reference)))))<1e-8,
 abs(AIC(fit)-AIC(reference))<1e-8)
close(fit)
cat('fixest loaded',args[[1L]],'dbreg: extraction dispatch and values pass\n')
