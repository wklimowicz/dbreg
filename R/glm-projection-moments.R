# Bounded direct nuisance solves for moderate FE dimensions. The database sends
# only FE cross-products and FE-by-design moments; no observation vectors enter R.
dbglm_project_moments = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop, cache) {
  info = lapply(ops,function(op) dbglm_get_query(conn,paste0("SELECT ",op$f," AS id,FSUM(w*(",op$v,")*(",op$v,")) AS h",
    if(length(columns)) paste0(",",paste(sprintf("FSUM(w*(%s)*%s) AS b%d",op$v,columns,seq_along(columns)),collapse=",")),
    " FROM ",table," GROUP BY ",op$f," ORDER BY ",op$f)))
  ids = unlist(lapply(seq_along(ops),function(i) paste(i,info[[i]]$id,sep="/")),use.names=FALSE)
  if(!is.null(cache$moment_keys) && !identical(ids,cache$moment_keys)) stop("The fixed-effect sample changed within fitting.",call.=FALSE)
  cache$moment_keys = ids
  n = length(ids)
  H = matrix(0,n,n)
  B = matrix(0,n,length(columns))
  starts = cumsum(c(0L,vapply(info,nrow,integer(1))))
  for(i in seq_along(ops)) {
    ii = seq.int(starts[i]+1L,starts[i+1L])
    H[cbind(ii,ii)] = info[[i]]$h
    if(length(columns)) B[ii,] = as.matrix(info[[i]][,paste0("b",seq_along(columns)),drop=FALSE])
    if(i>1) for(j in seq_len(i-1L)) {
      oi = ops[[i]];oj = ops[[j]]
      edge = dbglm_get_query(conn,paste0("SELECT ",oi$f," AS i,",oj$f," AS j,FSUM(w*(",oi$v,")*(",oj$v,")) AS h FROM ",
        table," GROUP BY ",oi$f,",",oj$f))
      ir = starts[i]+match(edge$i,info[[i]]$id)
      jr = starts[j]+match(edge$j,info[[j]]$id)
      H[cbind(ir,jr)] = edge$h
      H[cbind(jr,ir)] = edge$h
    }
  }
  diagonal = sqrt(pmax(diag(H),0))
  if(is.null(cache$moment_basis)) {
    positive = which(diagonal>0 & is.finite(diagonal))
    C = H[positive,positive,drop=FALSE]/outer(diagonal[positive],diagonal[positive])
    factor = suppressWarnings(chol(C,pivot=TRUE,tol=1e-12))
    rank = attr(factor,"rank")
    if(rank<1) stop("No identified fixed-effect span remains.",call.=FALSE)
    cache$moment_basis = positive[attr(factor,"pivot")[seq_len(rank)]]
  }
  keep = cache$moment_basis
  scale = diagonal[keep]
  factor = tryCatch(chol(H[keep,keep,drop=FALSE]/outer(scale,scale)),error=function(e) NULL)
  if(is.null(factor)) stop("GLM fixed-effect information is singular.",call.=FALSE)
  solve_rhs = function(rhs) backsolve(factor,forwardsolve(t(factor),rhs[keep,,drop=FALSE]/scale))/scale
  coefficients = matrix(0,n,length(columns))
  if(length(columns)) coefficients[keep,] = solve_rhs(B)
  lookup = NULL
  on.exit(if(!is.null(lookup)) try(drop(lookup),silent=TRUE),add=TRUE)
  for(cycle in seq_len(min(control$fe_max_iter,8L))) {
    if(!length(columns)) {
      projected = create(paste("SELECT __dbglm_key,",paste(c(identifiers,"w"),collapse=","),"FROM",table))
      drop(table)
      return(list(table=projected,cycles=cycle,residual=0,method="moments"))
    }
    frame = data.frame(op=rep(seq_along(ops),vapply(info,nrow,integer(1))),check.names=FALSE)
    names(frame) = "__glm_op"
    frame[["__glm_group"]] = as.numeric(unlist(lapply(info,`[[`,"id"),use.names=FALSE))
    frame[paste0("c",seq_along(columns))] = as.data.frame(coefficients)
    if(isTRUE(detect_backend(conn)$is_databricks)) {
      lookup_sql = dbglm_values_sql(conn,frame)
      if(is.null(lookup_sql)) stop("FE coefficient lookup exceeds the SQL transfer limit.",call.=FALSE)
      lookup = create(lookup_sql)
    } else {
      registered = basename(tempfile("dbglm_coef_"))
      tryCatch({
        duckdb::duckdb_register(conn,registered,frame)
        lookup = create(paste("SELECT * FROM",DBI::dbQuoteIdentifier(conn,registered)))
      },finally=try(duckdb::duckdb_unregister(conn,registered),silent=TRUE))
    }
    ctes = vapply(seq_along(ops),function(i) paste0("g",i," AS (SELECT * FROM ",lookup," WHERE __glm_op=",i,")"),character(1))
    expressions = vapply(seq_along(columns),function(j) {
      terms = vapply(seq_along(ops),function(i) {
        v = if(ops[[i]]$v=="1.0") "1.0" else paste0("s.",ops[[i]]$v)
        sprintf("(%s)*COALESCE(g%d.c%d,0)",v,i,j)
      },character(1))
      paste0("(s.",columns[j],"-(",paste(terms,collapse=" + "),")) AS ",columns[j])
    },character(1))
    joins = vapply(seq_along(ops),function(i) paste0(" LEFT JOIN g",i," ON s.",ops[[i]]$f,"=g",i,".__glm_group"),character(1))
    projected = create(paste0("WITH ",paste(ctes,collapse=",")," SELECT ",paste(paste0("s.",c("__dbglm_key",identifiers,"w")),collapse=","),",",
      paste(expressions,collapse=",")," FROM ",table," s",paste(joins,collapse="")))
    drop(lookup);lookup = NULL
    residual_rhs = matrix(0,n,length(columns))
    for(i in seq_along(ops)) {
      op = ops[[i]]
      terms = sprintf("FSUM(w*(%s)*%s) AS r%d",op$v,columns,seq_along(columns))
      residuals = dbglm_get_query(conn,paste0("SELECT ",op$f," AS id,",paste(terms,collapse=",")," FROM ",projected," GROUP BY ",op$f," ORDER BY ",op$f))
      ii = seq.int(starts[i]+1L,starts[i+1L])
      residual_rhs[ii,] = as.matrix(residuals[,paste0("r",seq_along(columns)),drop=FALSE])
    }
    normalized = residual_rhs/pmax(diagonal,.Machine$double.xmin)
    residual = max(abs(sweep(normalized,2,scales,"/")))
    if(residual<=target) {
      drop(table)
      return(list(table=projected,cycles=cycle,residual=residual,method="moments"))
    }
    drop(projected)
    coefficients[keep,] = coefficients[keep,,drop=FALSE]+solve_rhs(residual_rhs)
  }
  stop("GLM fixed-effect projection did not converge to its group-equation tolerance.",call.=FALSE)
}
