library(DBI)
library(dbreg)
cw <- dbConnect(duckdb::duckdb(config=list(temp_directory="")))
set.seed(117)
for(nfe in c(1,2,6)) {
 wdat <- data.frame(key=1:160,check.names=FALSE)
 names(wdat) <- '__dbglm_key'
 wdat$w <- runif(160,.1,3)
 wdat$x <- rnorm(160);wdat$large <- 1e6*rnorm(160);wdat$z <- rnorm(160)
 for(j in seq_len(nfe)) wdat[[paste0('f',j)]] <- sample(1:(j+2),160,TRUE)
 fe <- paste0('f',seq_len(nfe));attr(fe,'slopes') <- rep(list(NULL),nfe)
 dbWriteTable(cw,'raw',wdat)
 owned <- character()
 create <- function(sql, materialize = TRUE) {
   name <- basename(tempfile('wls_'));owned <<- c(owned,name)
   dbExecute(cw,paste('CREATE TEMP TABLE',dbQuoteIdentifier(cw,name),'AS',sql))
   as.character(dbQuoteIdentifier(cw,name))
 }
 drop <- function(table) {dbExecute(cw,paste('DROP TABLE',table));owned <<- setdiff(owned,gsub('"','',table,fixed=TRUE))}
 ans <- dbreg:::solve_partial_wls(cw,'raw',c('x','large'),fe,dbglm_control(),create,drop)
 got <- dbGetQuery(cw,paste('SELECT * FROM',ans$table,'ORDER BY __dbglm_key'))
 D <- model.matrix(as.formula(paste('~',paste(sprintf('factor(f%d)',seq_len(nfe)),collapse='+'))),wdat)
 ref <- lm.wfit(D,as.matrix(wdat[,c('x','large','z')]),wdat$w)$residuals
 actual <- as.matrix(got[,c('x','large','z')])
 # Projection tolerances are relative to each column's units. Check the
 # million-scaled regressor in its original units, keeping the span test below.
 normalized_actual <- sweep(actual,2,c(1,1e6,1),"/")
 normalized_ref <- sweep(ref,2,c(1,1e6,1),"/")
 expect_true(all(abs(normalized_actual-normalized_ref)<=1e-10+1e-10*abs(normalized_ref)))
 # The removed component lies in the additive FE span, not just joint cells.
 correction <- as.matrix(wdat[,c('x','large','z')])-actual
 span_error <- lm.wfit(D,correction,wdat$w)$residuals
 expect_true(max(abs(span_error[,c('x','z')]))<1e-8)
 expect_true(max(abs(span_error[,'large']))<1e-7,info=paste(nfe,max(abs(span_error[,'large']))))
 types <- dbGetQuery(cw,paste('DESCRIBE SELECT x,large,z,w FROM',ans$table))$column_type
 expect_true(all(types=='DOUBLE'))
 drop(ans$table);dbRemoveTable(cw,'raw')
 expect_equal(length(owned),0L)
 expect_equal(length(dbListTables(cw)),0L)
}
dbDisconnect(cw,shutdown=TRUE)
