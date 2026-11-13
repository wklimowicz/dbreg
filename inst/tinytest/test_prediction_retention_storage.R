local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(11171)
  d = expand.grid(g=1:8,x=-1:1,rep=1:5)
  d$id=seq_len(nrow(d));d$w=runif(nrow(d),.4,2);d$o=d$rep/20
  d$y=rpois(nrow(d),exp(.3+.3*d$x+cos(d$g)/5+d$o))
  DBI::dbWriteTable(conn,"source",d,temporary=TRUE)
  dense=list(lm(y~x+factor(g),data=d,weights=w,offset=o),
    glm(y~x+factor(g),data=d,family=poisson(),weights=w,offset=o))
  for (database in c(FALSE, TRUE)) {
    prior=options(dbreg.prediction_max_entries=if(database)4 else 2e6)
    models=list(dbreg(y~x|g,conn=conn,table="source",weights=~w,offset=~o,strategy="demean",retain="prediction"))
    if (exists("dbglm", mode = "function")) {
      models=c(models, list(dbglm(y~x|g,conn=conn,table="source",family=poisson(),weights=~w,offset=~o,retain="prediction")))
    }
    options(prior)
    for (i in seq_along(models)) {
      model=models[[i]]
      expect_identical(model$retained_rows,NULL)
      expect_identical(model$retained_inference,NULL)
      expect_error(fitted(model),"inference")
      expected=unname(predict(dense[[i]],type=if(i==2L)"link" else "response"))
      expect_equal(predict(model,d),expected,tolerance=1e-8)
      if (database) {
        tables=model$resource_state$tables
        expect_equal(length(tables),1L)
        expect_equal(as.numeric(DBI::dbGetQuery(conn,paste("SELECT COUNT(*) AS n FROM",DBI::dbQuoteIdentifier(conn,tables)))$n),8)
        output=predict(model,d,row_id="id")
        close(model)
        expect_equal(as.data.frame(output)[[output$prediction_column]],expected,tolerance=1e-8)
        close(output)
      } else {
        expect_identical(model$resource_state,NULL)
        close(model)
        expect_equal(predict(unserialize(serialize(model,NULL)),d),expected,tolerance=1e-8)
      }
    }
    expect_identical(DBI::dbListTables(conn),"source")
  }
})
