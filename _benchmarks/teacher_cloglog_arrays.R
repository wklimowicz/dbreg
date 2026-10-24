# Experimental Databricks predictor for dense integer FE ranks. The existing
# fit uses broadcast joins for FE coefficients. Here each FE's coefficients
# form a bounded constant array indexed by its rank, eliminating those joins.
# This changes only the predictor construction, not IRLS, its certificates,
# the moment solver, the sample or covariance formulas.
teacher_array_working = function(conn,base,x,beta,c_eta,gamma,ops,info,link,inline_sql) {
  ids = lapply(info,function(inf) as.numeric(inf[['id']]))
  maxima = vapply(ids,max,numeric(1))
  sizes = lengths(ids)
  if(any(!is.finite(maxima)) || any(maxima>1000000) || sum(maxima)>1000000 ||
     any(maxima>4*pmax(sizes,1))) return(NULL)
  if(any(vapply(ids,function(id) any(id<1 | id!=floor(id)),logical(1)))) return(NULL)
  starts = cumsum(c(0L,sizes))
  terms = vapply(seq_along(ops),function(i) {
    values = numeric(maxima[[i]])
    values[ids[[i]]] = gamma[starts[[i]]+seq_len(sizes[[i]])]
    if(any(!is.finite(values))) stop('Nonfinite FE coefficients in array predictor.')
    payload = jsonlite::toJSON(values,digits=17,auto_unbox=FALSE)
    encoded = gsub('[\r\n]','',jsonlite::base64_enc(charToRaw(enc2utf8(payload))))
    a = paste0("from_json(CAST(unbase64('",encoded,"') AS STRING),'ARRAY<DOUBLE>')")
    effect = paste0('COALESCE(element_at(',a,',CAST(s.',ops[[i]][['f']],' AS INT)),0.0)')
    v = ops[[i]][['v']]
    if(v=='1.0') effect else paste0('s.',v,'*(',effect,')')
  },character(1))
  eta = paste(c('s.__dbglm_o',if(c_eta!=0) dbreg:::dbglm_number(c_eta),
    sprintf('(%s)*s.%s',dbreg:::dbglm_number(beta),x),terms),collapse=' + ')
  st = inline_sql(paste0('SELECT s.*,CAST((',eta,') AS DOUBLE) AS __dbglm_eta FROM ',base,' s'))
  wk = inline_sql(paste0(dbreg:::dbglm_working_sql(st,character(),numeric(),link,eta='__dbglm_eta'),'SELECT * FROM working'))
  list(state=st,working=wk,lookup=NULL)
}

teacher_array_fit = function(original = dbreg:::dbglm_fe_fit) {
  expressions = as.list(body(original))
  positions = which(vapply(expressions,function(e) is.call(e) &&
    identical(e[[1]],as.name('=')) && identical(e[[2]],as.name('make_working')),logical(1)))
  stopifnot(length(positions)==1L)
  i = positions[[1]]
  old = expressions[[i]][[3]]
  expressions[[i]] = substitute({
    original_make_working = ORIGINAL
    make_working = function(beta_val,c_eta_val,gamma_val=NULL,ops=NULL,info=NULL) {
      if(use_moments && isTRUE(backend_info[['is_databricks']]) &&
         !is.null(gamma_val) && !all(gamma_val==0)) {
        result = teacher_array_working(conn,base,x,beta_val,c_eta_val,gamma_val,ops,info,family[['link']],inline_sql)
        if(!is.null(result)) return(result)
      }
      original_make_working(beta_val,c_eta_val,gamma_val,ops,info)
    }
  },list(ORIGINAL=old))
  result = original
  body(result) = as.call(expressions)
  environment(result) = list2env(list(teacher_array_working=teacher_array_working),parent=environment(original))
  result
}

teacher_array_check = function(conn,out) {
  base = '(SELECT * FROM VALUES (1,1,1,0.0,1.0,0.1,-1.0),(2,3,2,1.0,2.0,-0.1,2.0),(3,2,1,0.0,3.0,0.0,0.0) AS t(__dbglm_key,f1,f2,__dbglm_y,__dbglm_a,__dbglm_o,x))'
  inline = function(sql) paste0('(',dbreg:::dbglm_sql(conn,sql),')')
  args = list(conn=conn,base=base,x='x',beta=.3,c_eta=-3,gamma=c(.2,-.4,.1,-.2),
    ops=list(list(f='f1',v='1.0'),list(f='f2',v='1.0')),info=list(list(id=c(1,3)),list(id=c(1,2))),
    link='cloglog',inline_sql=inline)
  w = do.call(teacher_array_working,args)
  query = paste('SELECT eta,mu,score FROM',w[['working']],'ORDER BY __dbglm_key')
  actual = DBI::dbGetQuery(conn,query)
  eta = -3+c(.1,-.1,0)+.3*c(-1,2,0)+c(.2,-.4,0)+c(.1,-.2,.1)
  family = binomial('cloglog')
  mu = family[['linkinv']](eta)
  score = c(1,2,3)*family[['mu.eta']](eta)*(c(0,1,0)-mu)/(mu*(1-mu))
  stopifnot(isTRUE(all.equal(actual,data.frame(eta=eta,mu=mu,score=score),tolerance=1e-13)))
  # Very sparse ranks retain the old join path instead of allocating an array.
  args[['info']][[1]][['id']] = c(1,1000001)
  stopifnot(is.null(do.call(teacher_array_working,args)))
  plan = DBI::dbGetQuery(conn,paste('EXPLAIN FORMATTED',query))
  writeLines(as.character(plan[[1]]),file.path(out,'array_predictor_explain.txt'))
  cat('Array predictor: eta/mu/score parity and sparse-rank fallback passed\n'); flush.console()
}

# On larger samples, global diagnostics need not be evaluated for every FE
# marginal. Separate global, marginal and pair aggregates into three branches.
teacher_split_packed_sql = function(conn,source,fe,value_terms,global_terms=character(),limit=NULL) {
  if(!length(global_terms)) return(hybrid_moments_sql(conn,source,fe,value_terms,global_terms,limit))
  alias = function(terms) sub('.* AS ([A-Za-z_0-9]+)$','\\1',terms)
  null = function(terms) paste0('CAST(NULL AS DOUBLE) AS ',alias(terms))
  global = paste0('SELECT 0 AS __set,CAST(NULL AS BIGINT) AS __i,CAST(NULL AS BIGINT) AS __j,',
    paste(c(global_terms,null(value_terms)),collapse=','),' FROM ',source)
  singles = paste0('SELECT __set,__i,__j,',paste(c(null(global_terms),alias(value_terms)),collapse=','),' FROM (',
    grouped_moments_sql(conn,source,fe,value_terms,kind='main'),') singles')
  pairs = grouped_moments_sql(conn,source,fe,value_terms,global_terms,kind='pairs')
  sql = paste(global,singles,pairs,sep=' UNION ALL ')
  if(!is.null(limit)) sql = paste0('SELECT * FROM (',sql,') limited ORDER BY __set LIMIT ',sprintf('%.0f',limit))
  names = alias(global_terms)
  json_fields = paste(sprintf("'%s',%s",names,names),collapse=',')
  paste0('SELECT __set,__i,__j,',paste(alias(value_terms),collapse=','),
    ',CASE WHEN __set=0 THEN to_json(named_struct(',json_fields,')) END AS __globals FROM (',sql,') packed')
}
