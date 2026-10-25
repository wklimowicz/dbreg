pkgload::load_all()
local({
 set.seed(714)
 n = 90L
 for (kind in c("integer64","enum","date")) {
  key = switch(kind,integer64=bit64::as.integer64(rep(c("9007199254740993","9007199254740994","9007199254740995"),30)),
   enum=factor(rep(c("b","c","a"),30),levels=c("c","b","a")),
   date=rep(as.Date(c("2020-01-01","2020-02-01","2020-03-01")),30))
  d = data.frame(key=key,y=rep(1:3,30)+rnorm(n,sd=.1),id=seq_len(n))
  fit = dbreg::feols(y~i(key),d)
  truth = as.numeric(predict(lm(y~factor(as.character(key)),d),d))
  stopifnot(max(abs(predict(fit,d)-truth))<1e-8)
  conn = DBI::dbConnect(duckdb::duckdb())
  output = predict(fit,d,conn=conn,row_id="id")
  value = as.data.frame(output)
  stopifnot(max(abs(value[[output[["prediction_column"]]]][match(d[["id"]],value[["id"]])]-truth))<1e-8)
  if(kind=="enum") stopifnot(identical(names(coef(fit)),c("(Intercept)","key::b","key::a")))
  close(output)
  close(fit)
  DBI::dbDisconnect(conn,shutdown=TRUE)
  message(kind," predictions match dummy regression")
 }
})
local({
 old = options(dbreg.i_max_columns=2)
 on.exit(options(old))
 set.seed(715)
 d = data.frame(f=rep(letters[1:3],20),y=rnorm(60))
 fit = dbreg::feols(y~i(f),d)
 stopifnot(length(coef(fit))==3L,max(abs(predict(fit,d)-as.numeric(predict(lm(y~f,d),d))))<1e-8)
 close(fit)
 message("Two generated indicator columns fit at the configured two-column limit")
})
