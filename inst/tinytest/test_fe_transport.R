library(DBI)
local({
  conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  d=data.frame(f1=c(1,2,3,1,2,3,1,3),f2=c(1,2,1,2,1,2,1,2),f3=c(1,1,2,2,2,1,2,1),w=c(.2,2,.5,3,.8,1.1,.7,1.3),v1=c(1,3,-2,.5,4,2,-1,.7),v2=c(.5,-1,2,3,1,-2,.3,.9))
  DBI::dbWriteTable(conn,"summary_source",d)
  sizes=c(3,2,2);K=length(sizes);width=2L
  limits=dbreg:::projection_limits()
  budget=function()list2env(list(rows=0,cells=0,upload_bytes=0),parent=emptyenv())
  actual=dbreg:::query_fe_summaries(conn,"summary_source",c("f1","f2","f3"),c("v1","v2"),sizes,budget(),limits)
  D=do.call(cbind,lapply(d[c("f1","f2","f3")],function(group)model.matrix(~factor(group)-1)))
  expected_H=crossprod(D,d[["w"]]*D)
  expected_B=crossprod(D,d[["w"]]*as.matrix(d[c("v1","v2")]))
  system=dbreg:::fe_system(actual[["info"]],actual[["pairs"]])
  expect_equal(as.numeric(dbreg:::fe_multiply(system,diag(sum(sizes)))),as.numeric(expected_H),tolerance=1e-13)
  expect_equal(as.numeric(do.call(rbind,lapply(actual[["info"]],`[[`,"b"))),as.numeric(expected_B),tolerance=1e-13)
  expect_identical(actual[["statements"]],1L);expect_identical(actual[["transport"]],"lists")
  # Build an independent packet oracle from weighted incidence products.
  # Neither packet order nor entry order carries a sorting contract.
  combos=utils::combn(K,2L,simplify=FALSE)
  rows=list()
  for(i in seq_len(K)){
    ids=rev(seq_len(sizes[i]));where=sum(sizes[seq_len(i-1L)])+ids
    entries=vapply(seq_along(ids),function(j)paste(sprintf("%.17g",c(ids[j],diag(expected_H)[where[j]],expected_B[where[j],])),collapse=","),character(1))
    rows[[i]]=data.frame(kind=i,n=length(ids),bytes=nchar(paste(entries,collapse=";"),type="bytes"),overflow=0,payload=paste(entries,collapse=";"))
  }
  for(k in seq_along(combos)){
    pair=combos[[k]];pairs=unique(d[paste0("f",pair)]);pairs=pairs[rev(seq_len(nrow(pairs))),,drop=FALSE]
    entries=vapply(seq_len(nrow(pairs)),function(j){take=d[[paste0("f",pair[1])]]==pairs[[1]][j]&d[[paste0("f",pair[2])]]==pairs[[2]][j];paste(sprintf("%.17g",c(pairs[[1]][j],pairs[[2]][j],sum(d[["w"]][take]))),collapse=",")},character(1))
    payload=paste(entries,collapse=";")
    rows[[K+k]]=data.frame(kind=K+k,n=nrow(pairs),bytes=nchar(payload,type="bytes"),overflow=0,payload=payload)
  }
  packets=do.call(rbind,rev(rows))
  decode=function(result,lim=limits)dbreg:::fe_decode_summaries(result,sizes,width,budget(),lim,packets=TRUE)
  decoded=decode(packets)
  system=dbreg:::fe_system(decoded[["info"]],decoded[["pairs"]])
  expect_equal(as.numeric(dbreg:::fe_multiply(system,diag(sum(sizes)))),as.numeric(expected_H),tolerance=1e-13)
  expect_equal(as.numeric(do.call(rbind,lapply(decoded[["info"]],`[[`,"b"))),as.numeric(expected_B),tolerance=1e-13)
  # Truncating the final numeric field can leave a syntactically valid number:
  # a field-count check alone would miss the lost digits.
  bad=packets;bad[["payload"]][1]=substr(bad[["payload"]][1],1,nchar(bad[["payload"]][1])-1)
  expect_error(decode(bad),"byte count is incomplete")
  bad=packets;bad[["payload"]][1]="1,2";bad[["bytes"]][1]=3
  expect_error(decode(bad),"field count is incomplete")
  bad=packets;bad[["kind"]][1]=999
  expect_error(decode(bad),"summary metadata")
  bad=packets;bad[["n"]][1]=-1
  expect_error(decode(bad),"summary metadata")
  bad=packets;bad[["overflow"]][1]=NA_real_
  expect_error(decode(bad),"overflow metadata")
  expect_error(decode(packets[setdiff(names(packets),"bytes")]),"summary schema")
  bad=packets;bad[["overflow"]][1]=1
  attempt=dbreg:::projection_attempt("wire",decode(bad))
  expect_identical(attempt[["reason"]],"transfer_limit")
  small=limits;small[["summary_rows"]]=1
  expect_identical(dbreg:::projection_attempt("wire",decode(packets,small))[["reason"]],"transfer_limit")
  small=limits;small[["transfer_cells"]]=1
  expect_identical(dbreg:::projection_attempt("wire",decode(packets,small))[["reason"]],"transfer_limit")
  modify_entries=function(result,kind,fun){
    i=which(result[["kind"]]==kind);entries=strsplit(result[["payload"]][i],";",fixed=TRUE)[[1]]
    result[["payload"]][i]=paste(fun(entries),collapse=";")
    result[["bytes"]][i]=nchar(result[["payload"]][i],type="bytes");result
  }
  bad=modify_entries(packets,1,function(entries){entries[2]=entries[1];entries})
  expect_error(decode(bad),"Duplicate fixed-effect marginal")
  bad=modify_entries(packets,K+1,function(entries){entries[2]=entries[1];entries})
  expect_error(decode(bad),"Duplicate fixed-effect pair")
  bad=modify_entries(packets,1,function(entries){entries[1]=sub("^[^,]+,","999,",entries[1]);entries})
  expect_error(decode(bad),"ordinal keys")
  bad=modify_entries(packets,1,function(entries){entries[1]=sub("^([^,]+),[^,]+,","\\1,NaN,",entries[1]);entries})
  expect_identical(dbreg:::projection_attempt("wire",decode(bad))[["reason"]],"singular")
  bad=modify_entries(packets,1,function(entries){entries[1]=sub("^([^,]+),[^,]+,","\\1,abc,",entries[1]);entries})
  expect_error(decode(bad),"numeric fixed-effect packet")
  # The project requires cells strictly below 65,535 bytes. Exercise the
  # last safe entry and the first unsafe entry, independently of packet width.
  edge_width=16382L
  safe=paste(c("1","1.0","-1.0",rep("0.0",edge_width-1L)),collapse=",")
  wire=data.frame(kind=1,n=1,bytes=nchar(safe,type="bytes"),overflow=0,payload=safe)
  edge=dbreg:::fe_decode_summaries(wire,1,edge_width,budget(),limits,packets=TRUE)
  expect_equal(as.numeric(edge[["info"]][[1]][["b"]]),c(-1,numeric(edge_width-1L)))
  wire[["payload"]]=sub("-1.0","-10.0",safe,fixed=TRUE);wire[["bytes"]]=65535;wire[["overflow"]]=1
  expect_identical(dbreg:::projection_attempt("wire",dbreg:::fe_decode_summaries(wire,1,edge_width,budget(),limits,packets=TRUE))[["reason"]],"transfer_limit")
  wire[["overflow"]]=0
  expect_error(dbreg:::fe_decode_summaries(wire,1,edge_width,budget(),limits,packets=TRUE),"byte count is incomplete")
  # A real DBI delegate counts aggregate statements, without changing results.
  methods::setClass("dbreg_transport_test",contains="DBIConnection",slots=c(native="ANY",state="environment"))
  on.exit(methods::removeClass("dbreg_transport_test"),add=TRUE)
  methods::setMethod("dbGetInfo","dbreg_transport_test",function(dbObj,...)list(dbms.name=dbObj@state[["dbms"]]))
  methods::setMethod("dbGetQuery",c("dbreg_transport_test","character"),function(conn,statement,...){conn@state[["queries"]]=c(conn@state[["queries"]],statement);DBI::dbGetQuery(conn@native,statement,...)})
  state=new.env(parent=emptyenv());state[["dbms"]]="DuckDB";state[["queries"]]=character()
  spy=methods::new("dbreg_transport_test",native=conn,state=state)
  observed=dbreg:::query_fe_summaries(spy,"summary_source",c("f1","f2","f3"),c("v1","v2"),sizes,budget(),limits)
  expect_equal(length(state[["queries"]]),1)
  expect_equal(observed[["info"]],actual[["info"]])
  # Runtime overflow probes remain, even when marginal cardinalities are known.
  state[["queries"]]=character();small=limits;small[["summary_rows"]]=sum(sizes)+1
  attempt=dbreg:::projection_attempt("wire",dbreg:::query_fe_summaries(spy,"summary_source",c("f1","f2","f3"),c("v1","v2"),sizes,budget(),small))
  expect_identical(attempt[["reason"]],"transfer_limit")
  expect_equal(length(state[["queries"]]),1)
  expect_true(grepl(paste("LIMIT",sum(sizes)+2),state[["queries"]][1],fixed=TRUE))
  expect_error(dbreg:::query_fe_summaries(conn,"summary_source",c("f1","f2","f3"),c("v1","v2"),c(2,2,2),budget(),limits),"summary support changed")
  # Numeric JSON lookups round-trip doubles without escaping source text or
  # widening ordinal keys. They retain the same cumulative upload limit.
  state[["dbms"]]="Databricks"
  coefficients=matrix(c(.17,1e-200,-3.25,1e100,-.3,0),3,2)
  sql=dbreg:::fe_coefficient_sql(spy,coefficients,budget(),limits)
  expect_true(grepl("explode(from_json",sql,fixed=TRUE))
  payload=strsplit(strsplit(sql,"from_json('",fixed=TRUE)[[1]][2],"','ARRAY<ARRAY<DOUBLE>>'",fixed=TRUE)[[1]][1]
  expect_true(startsWith(payload,"[[")&&endsWith(payload,"]]"))
  encoded_rows=strsplit(substr(payload,3,nchar(payload)-2),"],[",fixed=TRUE)[[1]]
  parsed=lapply(encoded_rows,function(row){fields=strsplit(row,",",fixed=TRUE)[[1]];expect_true(all(grepl("^-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?$",fields)));as.numeric(fields)})
  expect_equal(do.call(rbind,parsed),cbind(seq_len(3),coefficients),tolerance=1e-15,check.attributes=FALSE)
  low=limits;low[["upload_bytes"]]=1
  expect_identical(dbreg:::projection_attempt("wire",dbreg:::fe_coefficient_sql(spy,coefficients,budget(),low))[["reason"]],"upload_limit")
  for(bad in list(0,-1,32001,NA_real_,Inf,"small",c(1,2))){
    expect_error(dbreg:::fe_packet_sql(conn,"SELECT 1",2,2,packet_bytes=bad),"packet_bytes")
  }
})
