# Runnable example: source(system.file("examples/fixest-compatibility.R", package="dbreg"))
# During development first use pkgload::load_all().
local({
  conn = NULL
  models = outputs = list()
  on.exit({
    for (output in outputs) close(output)
    for (model in models) close(model)
    if (!is.null(conn) && DBI::dbIsValid(conn)) DBI::dbDisconnect(conn,shutdown=TRUE)
  })
  set.seed(2410)
  n = 600L
  d = data.frame(id=seq_len(n), x=rnorm(n), z=rnorm(n),
    g=rep(1:20,each=30), period=rep(1:10,60),
    group=rep(c("control","treated","other"),200), w=runif(n,.2,1))
  d[["y"]] = .4*d[["x"]]-.2*d[["z"]]+sin(d[["g"]])+rnorm(n,sd=.3)
  d[["b"]] = rbinom(n,1,1-exp(-exp(.2*d[["x"]]+.1*sin(d[["g"]]))))

  # Retention is explicit; compact fits remain the default.
  fit = dbreg::dbreg(y~x+i(group,z,ref="control")|g+period,data=d,weights=~w,retain="all")
  models[[length(models)+1L]] = fit
  print(fixest::coeftable(fit,cluster=~g+period))
  print(df.residual(fit))
  print(stats::case.names(fit)[1:6])
  newdata = d[c(2,7,5),c("id","x","z","group","g","period")]
  print(predict(fit,newdata)) # Neither outcome nor fitting weights are needed.

  # A small portable learned-effect fit can predict after close/save/restore.
  predicted = predict(fit,newdata)
  portable = unserialize(serialize(fit,NULL))
  close(fit)
  stopifnot(isTRUE(all.equal(predict(portable,newdata),predicted)))

  # The caller computes each weight update and stores it in a source column.
  for (weights in list(d[["w"]],d[["w"]]^2)) {
    d[["current_weight"]] = weights
    refit = dbreg::dbglm(b~x+z|g,data=d,family=binomial("cloglog"),weights=~current_weight)
    print(coef(refit))
    close(refit)
  }

  # SQL sources and output use explicit keys and have no implicit row order.
  conn = DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(conn,"observations",d)
  sql_fit = dbreg::dbreg(y~x+z|g,conn=conn,table="observations",weights=~w,row_id="id",retain="all")
  models[[length(models)+1L]] = sql_fit
  output = predict(sql_fit,conn=conn,table="observations",row_id="id")
  outputs[[length(outputs)+1L]] = output
  print(head(output))
  close(sql_fit) # Output owns a separate lease and remains readable.
  print(head(output))
  close(output)
})
