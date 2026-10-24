# Experimental SQL only. No modifications to package source or normal defaults.
glm_experiment_source = function(sql) {
  start = regexpr("FROM (WITH predictor AS",sql,fixed=TRUE)[[1]] + 5L
  stopifnot(start > 5L)
  chars = strsplit(substr(sql,start,nchar(sql)),"",fixed=TRUE)[[1]]
  depth = 0L; quote = ""
  for (i in seq_along(chars)) {
    ch = chars[[i]]
    if (nzchar(quote)) { if (ch==quote) quote=""; next }
    if (ch %in% c("'",'"',"`")) {quote=ch;next}
    if (ch=="(") depth=depth+1L
    if (ch==")") depth=depth-1L
    if (depth==0L) return(substr(sql,start,start+i-1L))
  }
  stop("Could not find the working source.")
}

glm_experiment_terms = function(conn, x) {
  c(dbreg:::dbglm_deviance_term(conn),"MIN(mu) AS mu_min","MAX(mu) AS mu_max",
    "SUM(CASE WHEN NOT isfinite(w) OR w<=0 OR NOT isfinite(eta) THEN 1 ELSE 0 END) AS bad",
    dbreg:::dbglm_moment_terms(c(x,"z"),"w"),
    sprintf("%s AS raw_s%d",vapply(paste0(x,"*score"),function(e) dbreg:::dbglm_exact_sum(conn,e),character(1)),seq_along(x)),
    sprintf("ABS(%s)/SQRT(SUM(w*%s*%s)) AS s%d",
      vapply(paste0(x,"*score"),function(e) dbreg:::dbglm_exact_sum(conn,e),character(1)),x,x,seq_along(x)))
}

glm_experiment_two_branches = function(source,fe,value_terms,global_terms) {
  K=length(fe);pairs=utils::combn(K,2L,simplify=FALSE)
  alias=function(x) sub(".* AS ([A-Za-z_0-9]+)$","\\1",x)
  values=alias(value_terms);globals=alias(global_terms)
  grouped=function(groups,sets,terms) {
    masks=vapply(groups,function(g) sum(2^(K-setdiff(seq_len(K),g))),numeric(1))
    key=function(j) {
     cases=vapply(seq_along(groups),function(i) {
      g=groups[[i]];if(length(g)<j) return("")
      sprintf("WHEN %.0f THEN CAST(%s AS BIGINT)",masks[[i]],fe[g[[j]]])
     },character(1))
     if(!any(nzchar(cases))) return("CAST(NULL AS BIGINT)")
     paste("CASE __gid",paste(cases,collapse=" "),"ELSE CAST(NULL AS BIGINT) END")
    }
    fields=c(paste("CASE __gid",paste(sprintf("WHEN %.0f THEN %d",masks,sets),collapse=" "),"END AS __set"),
      paste(key(1L),"AS __i"),paste(key(2L),"AS __j"),
      vapply(c(globals,values),function(n) if(n %in% alias(terms)) n else paste("CAST(NULL AS DOUBLE) AS",n),character(1)))
    paste0("SELECT ",paste(fields,collapse=",")," FROM (SELECT GROUPING_ID(",paste(fe,collapse=","),") AS __gid,",
      paste(c(fe,terms),collapse=",")," FROM ",source," GROUP BY GROUPING SETS (",
      paste(vapply(groups,function(g) paste0("(",paste(fe[g],collapse=","),")"),character(1)),collapse=","),")) grouped")
  }
  sql=paste(grouped(c(list(integer()),as.list(seq_len(K))),0:K,c(global_terms,value_terms)),
    grouped(pairs,K+seq_along(pairs),value_terms[1]),sep=" UNION ALL ")
  fields=vapply(values,function(v) paste0("CASE WHEN __set=0 THEN CAST(NULL AS DOUBLE) ELSE ",v," END AS ",v),character(1))
  paste0("SELECT __set,__i,__j,",paste(fields,collapse=","),",CASE WHEN __set=0 THEN to_json(named_struct(",
    paste(sprintf("'%s',%s",globals,globals),collapse=","),")) END AS __globals FROM (",sql,") packed")
}

# Observed-curvature directions for a bounded number of initial cloglog steps.
# The final steps, convergence certificate and inference remain Fisher-based.
# This is deliberately a benchmark-only prototype, not a public solver option.
glm_experiment_newton = function(switch_at = "adaptive", score_tol = 1e-8) {
  old_working=dbreg:::dbglm_working_sql
  old_evaluate=dbreg:::dbglm_fe_evaluate
  old_fit=dbreg:::dbglm_fe_fit
  state=new.env(parent=emptyenv());state[["calls"]]=0L;state[["last_fisher"]]=FALSE;state[["switch"]]=FALSE
  working=function(table,x,beta,link,eta=NULL) {
    sql=old_working(table,x,beta,link,eta)
    if(link!="cloglog") return(sql)
    state[["calls"]]=state[["calls"]]+1L
    if(state[["switch"]] || is.numeric(switch_at) && state[["calls"]]>=switch_at) return(sql)
    sql=sub("working AS (","__experiment_fisher AS (",sql,fixed=TRUE)
    curvature=paste0("__dbglm_a*(CASE WHEN __dbglm_y=0 THEN __t WHEN __t<0.0001 THEN ",
      "__t*(0.5-__t/6.0+POWER(__t,3)/180.0) ELSE ",
      "(__t*EXP(-__t)/mu)*(__t/mu-1.0) END)")
    paste0(sql,", __experiment_rate AS (SELECT *,EXP(LEAST(eta,700.0)) AS __t FROM __experiment_fisher),",
      "__experiment_newton AS (SELECT *,",curvature," AS __observed_w FROM __experiment_rate),",
      "working AS (SELECT * EXCLUDE(w,z,__t,__observed_w),__observed_w AS w,",
      "eta-__dbglm_o+score/__observed_w AS z FROM __experiment_newton) ")
  }
  evaluate=function(conn,working,x,fe,collect_fe_moments=FALSE,moment_rows=NULL) {
    value=old_evaluate(conn,working,x,fe,collect_fe_moments,moment_rows)
    state[["last_fisher"]]=!grepl("__experiment_newton",working,fixed=TRUE)
    if(identical(switch_at,"adaptive") && is.finite(value[["score_max"]]) && value[["score_max"]]<=sqrt(score_tol)) state[["switch"]]=TRUE
    value
  }
  # Add a final-state eligibility condition without touching the line-search
  # score guard or any actual score/curvature/deviance tolerance.
  needle=quote(step==1 && relative<=control$glm_tol && current$score_max<=control$score_tol)
  replace=function(e) {
    if(identical(e,needle)) return(call("&&",e,quote(!grepl("__experiment_newton",working,fixed=TRUE))))
    if(is.call(e)) as.call(lapply(as.list(e),replace)) else e
  }
  fit=old_fit;body(fit)=replace(body(old_fit))
  stopifnot(!identical(body(fit),body(old_fit)))
  assignInNamespace("dbglm_working_sql",working,"dbreg")
  assignInNamespace("dbglm_fe_evaluate",evaluate,"dbreg")
  assignInNamespace("dbglm_fe_fit",fit,"dbreg")
  list(state=state,restore=function() {
    assignInNamespace("dbglm_working_sql",old_working,"dbreg")
    assignInNamespace("dbglm_fe_evaluate",old_evaluate,"dbreg")
    assignInNamespace("dbglm_fe_fit",old_fit,"dbreg")
  })
}
