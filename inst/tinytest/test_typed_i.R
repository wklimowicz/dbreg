local({
  set.seed(11061)
  n = 144L
  d = data.frame(f = rep(c("alpha", "Beta", "gamma"), 48), j = rep(1:3, 48),
    flag = rep(c(FALSE, TRUE), 72), x = rnorm(n), z = rnorm(n), w = runif(n, .4, 2), cl = rep(1:12, 12))
  d[["y"]] = .7 + .5*d[["x"]] + .4*(d[["f"]] == "gamma") + rnorm(n)
  cases = list(y~x+i(f), y~i(j,2)+x, y~i(f,ref=TRUE)+x,
    y~i(f,ref="@^[aB]")+x, y~i(f,keep="B")+x,
    y~i(f,x,keep=c("alpha","gamma")),y~i(f,flag),y~i(flag)+x,y~0+i(f),y~x+fixest::i(f),
    y~i(f,ref=NULL)+x,y~i(f,ref=c("alpha","Beta","gamma"))+x)
  for (formula in cases) {
    for (strategy in c("moments", "compress")) {
      fit = dbreg(formula, data = d, weights = "w", vcov = ~cl, strategy = strategy)
      # Independent ordinary dummy bases span explicit NULL and empty-reference cases.
      oracle_formula = if (identical(formula,y~x+fixest::i(f))) y~x+i(f) else
        if (identical(formula,y~i(f,ref=NULL)+x)) y~f+x else
        if (identical(formula,y~i(f,ref=c("alpha","Beta","gamma"))+x)) y~x else
        if (identical(formula,y~i(f,ref="@^[aB]")+x)) y~i(f,ref=c("alpha","Beta"))+x else formula
      oracle = fixest::feols(oracle_formula, d, weights = ~w, vcov = ~cl, notes = FALSE)
      expect_equal(predict(fit,d), as.numeric(predict(oracle,d)), tolerance = 1e-8,
        info=paste(deparse(formula),strategy))
      if (strategy == "compress") expect_equal(sort(predict(fit)), sort(as.numeric(predict(oracle,d))[!duplicated(model.matrix(oracle))]), tolerance = 1e-8)
      if (!identical(formula,y~i(f,ref=NULL)+x)) {
        expect_equal(coef(fit),coef(oracle)[names(coef(fit))], tolerance = 1e-8)
        expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(oracle)[names(coef(fit)),names(coef(fit)),drop=FALSE]),tolerance=1e-8)
      }
    }
  }
  # Stored selector results survive changes to the formula's environment and RDS.
  reference = "Beta"
  formula = y~x+i(f,ref=reference)
  fit = dbreg(formula,data=d,strategy="moments")
  expected = predict(fit,d)
  reference = "gamma"
  expect_equal(predict(fit,d),expected,tolerance=1e-12)
  serialized = unserialize(serialize(fit,NULL))
  expect_equal(predict(serialized,d),expected,tolerance=1e-12)
  unknown=d[1:3,];unknown[["f"]][1]="new";unknown[["x"]][2]=NA_real_
  expect_equal(predict(serialized,unknown),c(NA_real_,NA_real_,expected[3]),tolerance=1e-12)
  # Empty i() expansion reduces to the intercept; no malformed SELECT is emitted.
  for(strategy in c("moments","compress")) {
    empty=dbreg(y~i(f,ref=c("alpha","Beta","gamma")),data=d,strategy=strategy)
    expect_equal(as.numeric(coef(empty)),mean(d[["y"]]),tolerance=1e-8)
    expect_equal(predict(empty,d),rep(mean(d[["y"]]),n),tolerance=1e-8)
  }
  expect_error(dbreg(y~i(f,ref=NA),data=d),"finite, nonmissing")
  expect_error(dbreg(y~i(f,keep=Inf),data=d),"finite, nonmissing")
  expect_error(dbreg(y~i(f,ref="missing"),data=d),"match no level")
  expect_error(dbreg(y~i(f,bin=2),data=d),"not supported yet")
  expect_error(dbreg(y~i(f,2,ref="Beta"),data=d),"positional or named")
  expect_error(dbreg(y~i(f,f),data=d),"second i.*numeric or logical")
  conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  old=options(dbreg.i_max_columns=2)
  on.exit(options(old),add=TRUE)
  expect_error(dbreg(y~i(f,x),conn=conn,data=d),"i.*columns exceed")
  expect_identical(DBI::dbListTables(conn),character())
  old2=options(dbreg.i_max_columns=0)
  expect_error(dbreg(y~i(f),conn=conn,data=d),"positive integer")
  options(old2)
  options(dbreg.i_max_columns=2)
  expect_error(dbivreg(y~i(f)|x~i(j),conn=conn,data=d),"columns exceed.*IV blocks")
  # Non-i fits do not depend on i()'s resource option.
  options(dbreg.i_max_columns=0)
  expect_true(all(is.finite(coef(dbreg(y~x,conn=conn,data=d)))))
  options(old)
})

local({
  # Adjacent representable doubles have separate identities and exact selectors.
  set.seed(11062)
  k = 1+c(0,1,2)*.Machine[["double.eps"]]
  index=rep(1:3,40)
  d=data.frame(k=k[index],x=rnorm(120),w=runif(120,.4,1.8),cl=rep(1:10,12))
  d[["y"]]=d[["x"]]+index+rnorm(120)
  reference=k[2]
  X=cbind(1,d[["x"]],as.numeric(index==1),as.numeric(index==3))
  beta=as.numeric(qr.solve(X*sqrt(d[["w"]]),d[["y"]]*sqrt(d[["w"]])))
  residual=d[["y"]]-as.numeric(X%*%beta);bread=solve(crossprod(X,d[["w"]]*X))
  meat=crossprod(rowsum(X*(d[["w"]]*residual),d[["cl"]]))
  V=bread%*%meat%*%bread*10/9*119/(120-ncol(X))
  for(strategy in c("moments","compress")) {
    fit=dbreg(y~x+i(k,ref=reference),data=d,weights="w",vcov=~cl,strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(V),tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    expect_false(anyDuplicated(names(coef(fit)))>0L)
    keep=k[c(1,3)]
    selected=dbreg(y~x+i(k,keep=keep),data=d,weights="w",strategy=strategy)
    expect_equal(predict(selected,d),as.numeric(X%*%beta),tolerance=1e-8)
  }
})

local({
  # SQL integer keys span the whole signed domain and bypass DBI double reads.
  set.seed(11063)
  conn=DBI::dbConnect(duckdb::duckdb(bigint="numeric",config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  keys=c("-9223372036854775808","9007199254740993","9223372036854775807")
  index=rep(1:3,30)
  d=data.frame(k=keys[index],x=rnorm(90),y=rnorm(90))
  DBI::dbWriteTable(conn,"key_text",d)
  DBI::dbExecute(conn,"CREATE TABLE typed_keys AS SELECT CAST(k AS BIGINT) AS k,x,y FROM key_text")
  X=cbind(1,d[["x"]],as.numeric(index==2),as.numeric(index==3))
  beta=as.numeric(qr.solve(X,d[["y"]]))
  for(strategy in c("moments","compress")) {
    fit=dbreg(y~x+i(k),conn=conn,table="typed_keys",strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    if(strategy=="compress")expect_equal(sort(predict(fit)),sort(as.numeric(X%*%beta)),tolerance=1e-8)
    lossy=d;lossy[["k"]]=as.numeric(lossy[["k"]])
    expect_error(predict(fit,lossy),"must be exact")
  }
  for(type in c("REAL","DECIMAL(10,2)")) {
    table=paste0("(SELECT CAST(x AS ",type,") AS k,y,x FROM typed_keys) unsupported_key")
    expect_error(dbreg(y~x+i(k),conn=conn,table=table),"floating-point|SQL type")
  }
})

local({
  set.seed(11064)
  dates=as.Date("2020-01-01")+rep(0:2,40)
  d=data.frame(day=dates,x=rnorm(120),y=rnorm(120))
  X=model.matrix(~x+factor(day),d)
  beta=as.numeric(qr.solve(X,d[["y"]]))
  for(strategy in c("moments","compress")) {
    fit=dbreg(y~x+i(day),data=d,strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    strings=d;strings[["day"]]=as.character(dates)
    expect_equal(predict(fit,strings),as.numeric(X%*%beta),tolerance=1e-8)
  }
  # Small doubles whose %.17g fixed-point literal exceeds 18 decimal digits must
  # still match exact DOUBLE values in DuckDB rather than rounding via DECIMAL.
  set.seed(2026)
  k = min(runif(120))
  df = data.frame(y = c(1, 2, 3, 4), x = c(1, 2, 1, 2), g = c(k, k, 1.0, 1.0))
  fit = dbreg(y ~ i(g, ref = 1), data = df)
  expect_equal(unname(coef(fit)), c(3.5, -2.0), tolerance = 1e-12)
})
