# Translate only generated SQL tokens, preserving quoted source identifiers and
# literals. DuckDB retains its compensated aggregates; Spark uses double SUM.
dbglm_sql = function(conn, sql) {
  if (!isTRUE(detect_backend(conn)$is_databricks)) return(sql)
  tokens = regmatches(sql, gregexpr("'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|`(?:``|[^`])*`|[A-Za-z_][A-Za-z_0-9]*|[[:space:]]+|.", sql, perl=TRUE))[[1]]
  tokens[tokens == "FSUM"] = "SUM"
  tokens[tokens == "EXCLUDE"] = "EXCEPT"
  tokens[tokens == "VARCHAR"] = "STRING"
  repeat {
    positions = which(tokens == "isfinite")
    if (!length(positions)) break
    start = positions[1]
    open = start + 1L
    while (grepl("^[[:space:]]+$", tokens[open])) open = open + 1L
    stopifnot(tokens[open] == "(")
    depth = 1L; end = open
    while (depth > 0L) {
      end = end + 1L
      if (tokens[end] == "(") depth = depth + 1L
      if (tokens[end] == ")") depth = depth - 1L
    }
    arg = paste0(tokens[seq.int(open+1L,end-1L)],collapse="")
    value = paste0("(NOT isnan(CAST((",arg,") AS DOUBLE)) AND ABS(CAST((",arg,") AS DOUBLE)) <= 1.7976931348623157e308)")
    tokens = c(if(start>1L) tokens[seq_len(start-1L)],value,if(end<length(tokens)) tokens[seq.int(end+1L,length(tokens))])
  }
  paste0(tokens,collapse="")
}

dbglm_get_query = function(conn, sql) DBI::dbGetQuery(conn, dbglm_sql(conn, sql))

# Bounded lookup literals are FE-level coefficients, never observation data.
# A size limit prevents constructing an unbounded SQL statement. The caller
# owns materialization and cleanup just as for every other projector state.
dbglm_values_sql = function(conn, frame, max_bytes = 8*1024^2,
                              format = getOption("dbreg.glm_databricks_lookup", "json")) {
  if (!nrow(frame)) stop("A coefficient lookup cannot be empty.",call.=FALSE)
  if(!is.character(format) || length(format)!=1 || !format %in% c("json","values")) {
    stop("GLM FE lookup format must be json or values.",call.=FALSE)
  }
  if(any(vapply(frame,function(x) is.numeric(x) && any(!is.finite(x)),logical(1)))) {
    stop("Nonfinite projection coefficients.",call.=FALSE)
  }
  if(format=="json") {
    # Arrays avoid repeating field names for every FE. Seventeen digits preserve
    # double round trips (digits=NA rounds some extremes to infinity). Base64
    # keeps SQL string escaping independent of arbitrary key text/Unicode.
    payload = jsonlite::toJSON(frame,dataframe="values",digits=17,na="null")
    encoded = gsub("[\r\n]","",jsonlite::base64_enc(charToRaw(enc2utf8(payload))))
    if(nchar(encoded,type="bytes")+1024>max_bytes) return(NULL)
    fields = sprintf("CAST(item[%d] AS %s) AS %s",seq_along(frame)-1L,
      ifelse(vapply(frame,is.numeric,logical(1)),"DOUBLE","STRING"),
      DBI::dbQuoteIdentifier(conn,names(frame)))
    return(paste0("SELECT ",paste(fields,collapse=","),
      " FROM explode(from_json(CAST(unbase64('",encoded,"') AS STRING), ",
      "'ARRAY<ARRAY<STRING>>')) AS lookup(item)"))
  }
  fields = lapply(frame,function(x) {
    if(is.numeric(x)) {
      if(any(!is.finite(x))) stop("Nonfinite projection coefficients.",call.=FALSE)
      ifelse(x==0,"0.0",sprintf("%.17ge0",x))
    } else as.character(DBI::dbQuoteString(conn,as.character(x)))
  })
  # sprintf can already emit an exponent; append e0 only to decimal literals.
  for(j in which(vapply(frame,is.numeric,logical(1)))) {
    fields[[j]] = sub("(e[+-]?[0-9]+)e0$","\\1",fields[[j]])
  }
  rows = do.call(paste,c(fields,sep=","))
  if(sum(nchar(rows,type="bytes"))+3*nrow(frame)>max_bytes) return(NULL)
  names = paste(DBI::dbQuoteIdentifier(conn,names(frame)),collapse=",")
  paste0("SELECT * FROM VALUES ",paste0("(",rows,")",collapse=",")," AS lookup(",names,")")
}

dbglm_saturated_loglik = function(conn, base) {
  if (!isTRUE(detect_backend(conn)$is_databricks)) return(dbglm_get_query(conn,paste0(
    "SELECT FSUM(__dbglm_a*(CASE WHEN __dbglm_y=0 THEN 0.0 ELSE ",
    "__dbglm_y*LN(__dbglm_y)-__dbglm_y-lgamma(__dbglm_y+1) END)) AS ll FROM ",base))$ll)
  # Databricks has no lgamma builtin. A bounded count histogram is sufficient
  # for the likelihood constant, with R's canonical lgamma evaluation.
  n = dbglm_get_query(conn,paste("SELECT COUNT(DISTINCT __dbglm_y) AS n FROM",base))$n
  if(n>1e6) stop("Poisson count histogram exceeds the likelihood transfer limit.",call.=FALSE)
  h = dbglm_get_query(conn,paste("SELECT __dbglm_y AS y,SUM(__dbglm_a) AS a FROM",base,"GROUP BY __dbglm_y"))
  value = numeric(nrow(h)); positive = h$y>0
  value[positive] = h$y[positive]*log(h$y[positive])-h$y[positive]-lgamma(h$y[positive]+1)
  sum(h$a*value)
}

# Reuse Spark's bounded-plan window AP chunks with the GLM's scaled equation
# certificate and its tighter response-dependent inner tolerance.
dbglm_project_window = function(conn, table, columns, fe, control, create, drop, response) {
  if(!length(columns)) return(list(table=create(paste("SELECT __dbglm_key,",
    paste(c(fe,"w"),collapse=","),"FROM",table)),cycles=0L,residual=0,method="window"))
  norm = dbglm_get_query(conn,paste("SELECT",paste(sprintf("SUM(w*POWER(%s,2)) AS n%d",
    columns,seq_along(columns)),collapse=","),"FROM",table))
  scales = sqrt(pmax(as.numeric(norm[1,]),.Machine$double.xmin))
  target = control$fe_tol
  if(response) target = min(target,control$score_tol/(100*max(1,utils::tail(scales,1))))
  first = utils::tail(columns,1); rest = utils::head(columns,-1)
  vars_all = c(first, rest)
  tilde_names = paste0(vars_all, "_tilde")
  id_cols = unique(c(fe, "__dbglm_key"))
  tilde_cols = sprintf("CAST((%s) AS DOUBLE) AS %s", vars_all, tilde_names)
  cur_table = create(paste0("SELECT ", paste(c(id_cols, "w AS __w", tilde_cols), collapse = ", "), " FROM ", table))
  chunk_cycles = getOption("dbreg.spark_ap_chunk_cycles", 3L)
  iter = 0L
  max_abs = Inf
  while (iter < control$fe_max_iter && max_abs > target) {
    cycles = min(chunk_cycles, control$fe_max_iter - iter)
    iter = iter + cycles
    window_sql = sql_ap_window_sweeps(cur_table, c(id_cols, "__w"), tilde_names, fe, cycles)
    alt_table = create(window_sql)
    drop(cur_table)
    cur_table = alt_table
    diagnostics = vapply(fe, function(fe_k) {
      aliases = sprintf("dbreg_mean_%04d", seq_along(tilde_names))
      means = sprintf("ABS(SUM(__w * %s))/SQRT(SUM(__w))/%.17g AS %s", tilde_names, scales[match(vars_all, columns)], aliases)
      paste0("SELECT MAX(GREATEST(", paste(aliases, collapse = ", "),
             ")) AS max_abs FROM (SELECT ", fe_k, ", ", paste(means, collapse = ", "),
             " FROM ", cur_table, " GROUP BY ", fe_k, ") dbreg_ap_conv")
    }, character(1))
    convergence_query = paste0("SELECT GREATEST(", paste(sprintf("(%s)", diagnostics), collapse = ", "), ") AS max_abs")
    max_abs = as.numeric(dbglm_get_query(conn, convergence_query)[1, 1L])
  }
  result = create(paste("SELECT __dbglm_key,", paste(c(fe, "__w AS w",
    sprintf("%s_tilde AS %s", columns, columns)), collapse=","), "FROM", cur_table))
  drop(cur_table)
  list(table=result, cycles=iter, residual=max_abs, method="window")
}
