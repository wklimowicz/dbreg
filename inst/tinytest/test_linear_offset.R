local({
  set.seed(11100);n=180L
  d=data.frame(x=rnorm(n),z=rnorm(n),o=rnorm(n),g=rep(1:9,20),h=rep(1:6,30),w=runif(n,.2,2),cl=rep(1:12,15),keep=rep(c(TRUE,TRUE,FALSE),60))
  d[["endog"]]=.7*d[["z"]]+rnorm(n)
  d[["o"]]=d[["o"]]+.8*d[["x"]]
  d[["y"]]=d[["o"]]+.5*d[["x"]]+.7*d[["endog"]]+sin(d[["g"]])+rnorm(n)
  d[["o"]][1:3]=c(NA,Inf,NaN);d[["w"]][4]=0
  selected=d[["keep"]] & is.finite(d[["o"]]) & d[["w"]]>0
  data=d[selected,]
  for(strategy in c("moments","compress","demean"))for(kind in list("iid","hc1",~cl)){
    formula=if(strategy=="demean")y~x+z|g+h else y~x+z
    fit=dbreg(formula,data=d,weights="w",offset=~o,subset=~keep,strategy=strategy,vcov=kind)
    reference=fixest::feols(formula,data,weights=~w,offset=~o,vcov=kind,fixef.rm="none",ssc=fixest::ssc(K.fixef="full"),notes=FALSE)
    expect_equal(coef(fit),coef(reference),tolerance=1e-7)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(reference)),tolerance=1e-7)
    expect_equal(nobs(fit),nrow(data))
    expect_equal(predict(fit,data),as.numeric(predict(reference,data)),tolerance=1e-7)
    expect_identical(fit[["yvar"]],"y");expect_identical(fit[["offset"]],"o")
    expect_identical(all.vars(stats::formula(fit[["fml"]]))[1],"y")
    expect_error(predict(fit,data[setdiff(names(data),"o")]),"missing required offset")
    expect_error(predict(fit),"newdata")
    if(strategy!="demean"){
      interval=predict(fit,data,interval="confidence")
      expect_equal(interval[["fit"]],predict(fit,data),tolerance=1e-10)
      # Known offsets change location but not conditional standard errors.
      changed=data;changed[["o"]]=changed[["o"]]+2
      shifted=predict(fit,changed,interval="confidence")
      expect_equal(unname(as.matrix(shifted)-as.matrix(interval)),matrix(2,nrow(data),3),tolerance=1e-9)
    }
  }
  for(strategy in c("moments","compress","demean")){
    formula=if(strategy=="demean")y~x|g+h|endog~z else y~x|endog~z
    fit=dbivreg(formula,data=d,offset="o",weights="w",subset=~keep,strategy=strategy,vcov=~cl)
    # Compare to explicitly adjusted outcomes: no-FE IV offset= is incorrect
    # in the installed fixest version on this independent 2SLS witness.
    adjusted=data;adjusted[["y"]]=adjusted[["y"]]-adjusted[["o"]]
    reference=fixest::feols(formula,adjusted,weights=~w,vcov=~cl,fixef.rm="none",ssc=fixest::ssc(K.fixef="full"),notes=FALSE)
    reference_names=sub("^fit_","",names(coef(reference)));i=match(names(coef(fit)),reference_names)
    expect_equal(as.numeric(coef(fit)),as.numeric(coef(reference)[i]),tolerance=1e-7)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(reference)[i,i,drop=FALSE]),tolerance=1e-7)
    expect_equal(nobs(fit),nrow(data))
    expect_equal(predict(fit,data),as.numeric(predict(reference,adjusted))+data[["o"]],tolerance=1e-7)
    expect_identical(fit[["yvar"]],"y")
    if(strategy!="demean"){
      X=model.matrix(~x+endog,data);Z=model.matrix(~x+z,data);w=data[["w"]]
      P=crossprod(X,w*Z)%*%solve(crossprod(Z,w*Z))
      A=solve(P%*%crossprod(Z,w*X))%*%P
      beta=as.numeric(A%*%crossprod(Z,w*adjusted[["y"]]))
      e=adjusted[["y"]]-as.numeric(X%*%beta);G=length(unique(data[["cl"]]))
      V=A%*%crossprod(rowsum(Z*(w*e),data[["cl"]]))%*%t(A)*G/(G-1)*(nrow(data)-1)/(nrow(data)-ncol(X))
      expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
      expect_equal(as.numeric(vcov(fit)),as.numeric(V),tolerance=1e-8)
    }
  }
  # Weighted CRE prediction/covariance uses adjusted Y and raw source recipes.
  meanw=function(x)ave(seq_len(nrow(data)),data[["g"]],FUN=function(i)weighted.mean(x[i],data[["w"]][i]))
  X=cbind(1,data[["x"]],data[["z"]],meanw(data[["x"]]),meanw(data[["z"]]))
  w=data[["w"]];y=data[["y"]]-data[["o"]]
  beta=as.numeric(qr.solve(X*sqrt(w),y*sqrt(w)));e=y-as.numeric(X%*%beta)
  bread=solve(crossprod(X,w*X));G=length(unique(data[["cl"]]))
  V=bread%*%crossprod(rowsum(X*(w*e),data[["cl"]]))%*%bread*G/(G-1)*(nrow(data)-1)/(nrow(data)-ncol(X))
  fit=dbreg(y~x+z|g,data=d,offset="o",weights="w",subset=~keep,strategy="mundlak",vcov=~cl)
  expect_equal(as.numeric(coef(fit,fe=TRUE)),beta,tolerance=1e-8)
  expect_equal(as.numeric(vcov(fit,fe=TRUE)),as.numeric(V),tolerance=1e-8)
  expect_equal(predict(fit,data),as.numeric(X%*%beta)+data[["o"]],tolerance=1e-8)
  # Demean prediction masks bad offsets before jointly projecting adjusted Y.
  fit=dbreg(y~x+z|g+h,data=data,weights="w",offset="o",strategy="demean")
  changed=data;changed[["o"]][c(1,3)]=c(NA,Inf)
  found=predict(fit,changed);expect_true(all(is.na(found[c(1,3)])))
  good=changed[-c(1,3),];F=model.matrix(~factor(g)+factor(h),good);U=as.matrix(good[c("x","z")]);y=good[["y"]]-good[["o"]];w=good[["w"]]
  projected=cbind(U,y)-qr.fitted(qr(F*sqrt(w)),cbind(U,y)*sqrt(w))/sqrt(w)
  expected=as.numeric(projected[,1:2]%*%coef(fit)+y-projected[,3]+good[["o"]])
  expect_equal(found[-c(1,3)],expected,tolerance=1e-7)
  good[["y"]]=NULL
  expect_equal(suppressMessages(predict(fit,good)),as.numeric(projected[,1:2]%*%coef(fit))+good[["o"]],tolerance=1e-7)
})

local({
  set.seed(11101);n=100L
  d=data.frame(x=rnorm(n),o=rnorm(n),f=rep(1:4,25),g=rep(1:10,10),w=runif(n,.2,2))
  d[["y"]]=d[["o"]]+.5*d[["x"]]+.2*d[["x"]]^2+rnorm(n)
  for(strategy in c("moments","compress","demean","mundlak")){
    formula=if(strategy %in% c("demean","mundlak"))y~x+I(x^2)+factor(f)|g else y~x+I(x^2)+factor(f)
    fit=dbreg(formula,data=d,offset="o",weights="w",strategy=strategy)
    adjusted=d;adjusted[["y"]]=adjusted[["y"]]-adjusted[["o"]]
    reference=dbreg(formula,data=adjusted,weights="w",strategy=strategy)
    expect_equal(coef(fit,fe=TRUE),coef(reference,fe=TRUE),tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit,fe=TRUE)),as.numeric(vcov(reference,fe=TRUE)),tolerance=1e-8)
    expect_equal(predict(fit,d),predict(reference,adjusted)+d[["o"]],tolerance=1e-7)
  }
  # The original Y remains an independent regressor even when offset=Y.
  fit=dbreg(y~x+y,data=d,offset="y",strategy="moments")
  expect_equal(as.numeric(coef(fit)),rep(0,3),tolerance=1e-9)
  expect_equal(predict(fit,d),d[["y"]],tolerance=1e-9)
  # A raw predictor uses the usual generated target name, and unsafe names map.
  names(d)[names(d)=="x"]="dbreg_offset_outcome__"
  names(d)[names(d)=="o"]="offset value"
  fit=dbreg(y~dbreg_offset_outcome__,data=d,offset=~`offset value`,strategy="moments")
  ref=lm(y~dbreg_offset_outcome__+offset(`offset value`),d)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-9)
  expect_equal(predict(fit,d),as.numeric(predict(ref,d)),tolerance=1e-9)
  expect_identical(fit[["offset"]],"offset value")
  expect_error(dbreg(y~f,data=d,offset=~f+g),"must name")
  expect_error(dbivreg(y~f|g~w,data=d,offset=c("f","g")),"must name")
  bad=d;bad[["offset value"]]="bad"
  expect_error(dbreg(y~f,data=bad,offset="offset value"),"must be numeric")
  expect_error(predict(fit,bad),"must be numeric")
})

local({
  conn=DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(11102);n=80L
  d=data.frame(id=seq_len(n),x=rnorm(n),o=rnorm(n))
  d[["y"]]=d[["o"]]+.6*d[["x"]]+rnorm(n)
  DBI::dbWriteTable(conn,"offset_source",d)
  DBI::dbExecute(conn,"UPDATE offset_source SET o=CASE WHEN id=1 THEN CAST('NaN' AS DOUBLE) WHEN id=2 THEN CAST('Infinity' AS DOUBLE) WHEN id=3 THEN -1e308 ELSE o END, y=CASE WHEN id=3 THEN 1e308 ELSE y END")
  # Finite operands can still overflow their adjusted difference.
  fit=dbreg(y~x,conn=conn,table="offset_source",offset="o",strategy="moments",drop_missings=FALSE)
  ref=lm(y~x+offset(o),d[-(1:3),])
  expect_equal(nobs(fit),n-3)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-9)
  query=dbreg(y~x,conn=conn,table="offset_source",offset="o",strategy="moments",sql_only=TRUE)
  moments=dbreg(y~x,conn=conn,table="offset_source",offset="o",strategy="moments",data_only=TRUE)
  expect_equal(DBI::dbGetQuery(conn,query),moments,tolerance=1e-12)
  tables=DBI::dbListTables(conn)
  expect_error(dbreg(y~x,conn=conn,data=d,offset="missing"),"Offset column not found")
  expect_equal(DBI::dbListTables(conn),tables)
  expect_true(DBI::dbIsValid(conn))
  # Cast before subtraction, rather than overflowing a 32-bit SQL integer.
  d[["o"]]=as.integer(2000000000L-d[["id"]]);d[["y"]]=as.integer(-2000000000L+d[["id"]])
  DBI::dbWriteTable(conn,"integer_offset",d)
  fit=dbreg(y~x,conn=conn,table="integer_offset",offset="o",strategy="moments")
  ref=lm(I(as.numeric(y)-as.numeric(o))~x,d)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-7)
  expect_equal(unname(coef(fit)[["x"]]),unname(coef(ref)[["x"]]),tolerance=1e-6)
  d[["y"]]=d[["x"]]>0
  fit=dbreg(y~x,data=d,offset="o",strategy="moments")
  ref=lm(y~x+offset(o),d)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-7)
})
