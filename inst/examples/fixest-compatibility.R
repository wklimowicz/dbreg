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

  # Full adapters retain fitted rows and learned effects in the owning session.
  fit = dbreg::feols(y~x+i(group,z,ref="control")|g+period,d,weights=~w)
  models[[length(models)+1L]] = fit
  print(dbreg::coeftable(fit,cluster=~g+period))
  print(dbreg::degrees_freedom(fit,"k"))
  print(dbreg::obs(fit)[1:6])
  newdata = d[c(2,7,5),c("id","x","z","group","g","period")]
  print(predict(fit,newdata)) # Neither outcome nor fitting weights are needed.

  # A small portable learned-effect fit can predict after close/save/restore.
  predicted = predict(fit,newdata)
  portable = unserialize(serialize(fit,NULL))
  close(fit)
  stopifnot(isTRUE(all.equal(predict(portable,newdata),predicted)))

  # The caller computes and supplies each new weight vector. Lean releases
  # sample-dependent state after fitting; it avoids learned-effect recovery.
  for (weights in list(d[["w"]],d[["w"]]^2)) {
    refit = dbreg::feglm(b~x+z|g,d,family=binomial("cloglog"),weights=weights,lean=TRUE,notes=FALSE)
    print(coef(refit))
    close(refit)
  }

  # SQL sources and output use explicit keys and have no implicit row order.
  conn = DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(conn,"observations",d)
  sql_fit = dbreg::feols(y~x+z|g,conn=conn,table="observations",weights=~w,row_id="id")
  models[[length(models)+1L]] = sql_fit
  output = predict(sql_fit,conn=conn,table="observations",row_id="id")
  outputs[[length(outputs)+1L]] = output
  print(head(output))
  close(sql_fit) # Output owns a separate lease and remains readable.
  print(head(output))
  close(output)
})
