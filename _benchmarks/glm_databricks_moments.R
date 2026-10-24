# Experimental alternative with the same result schema as the production UNION.
grouped_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL, kind = 'all') {
  K = length(fe)
  pairs = utils::combn(K, 2L, simplify = FALSE)
  groups = c(if (length(global_terms)) list(integer()), as.list(seq_len(K)), pairs)
  sets = c(if (length(global_terms)) 0L, seq_len(K + length(pairs)))
  keep = if(kind=='pairs') sets>K else if(kind=='main') sets<=K else rep(TRUE,length(sets))
  groups = groups[keep]
  sets = sets[keep]
  masks = vapply(groups, function(g) sum(2^(K - setdiff(seq_len(K), g))), numeric(1))
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  mapping = paste(sprintf("WHEN %.0f THEN %d", masks, sets), collapse = " ")
  key_case = function(j) {
    cases = vapply(seq_along(groups), function(k) {
      g = groups[[k]]
      if (length(g) < j) return("")
      sprintf("WHEN %.0f THEN CAST(%s AS BIGINT)", masks[[k]], fe[g[[j]]])
    }, character(1))
    if(!any(nzchar(cases))) return('CAST(NULL AS BIGINT)')
    paste("CASE __gid", paste(cases, collapse = " "), "ELSE CAST(NULL AS BIGINT) END")
  }
  global_names = alias(global_terms)
  value_names = alias(value_terms)
  cols = c(paste0("CASE __gid ", mapping, " END AS __set"),
    paste(key_case(1L), "AS __i"), paste(key_case(2L), "AS __j"),
    if (length(global_names)) if(kind=='pairs') paste0('CAST(NULL AS DOUBLE) AS ',global_names) else sprintf("CASE WHEN __gid=%.0f THEN %s END AS %s", 2^K-1, global_names, global_names),
    sprintf("CASE WHEN __gid<>%.0f THEN %s END AS %s", 2^K-1, value_names[[1]], value_names[[1]]),
    if (length(value_names)>1L) if(kind=='pairs') paste0('CAST(NULL AS DOUBLE) AS ',value_names[-1]) else sprintf("CASE WHEN __gid IN (%s) THEN %s END AS %s",
      paste(2^K-1-2^(K-seq_len(K)), collapse=","), value_names[-1], value_names[-1]))
  sql = paste0("SELECT ", paste(cols, collapse=","), " FROM (SELECT GROUPING_ID(", paste(fe, collapse=","), ") AS __gid,",
    paste(c(fe, if(kind=='pairs') value_terms[1] else c(global_terms, value_terms)), collapse=","), " FROM ", source,
    " GROUP BY GROUPING SETS (", paste(vapply(groups, function(g) paste0("(", paste(fe[g], collapse=","), ")"), character(1)), collapse=","), ")) __dbglm_grouped")
  if (!is.null(limit)) sql = paste0(sql, " ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  sql
}

hybrid_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL) {
  sql = paste(grouped_moments_sql(conn,source,fe,value_terms,global_terms,kind='main'),
    grouped_moments_sql(conn,source,fe,value_terms,global_terms,kind='pairs'),sep=' UNION ALL ')
  if(!is.null(limit)) sql = paste0('SELECT * FROM (',sql,') __dbglm_hybrid ORDER BY __set LIMIT ',sprintf('%.0f',limit))
  sql
}

# Pack the single global summary into JSON instead of transferring its mostly
# NULL scalar columns on every FE-pair row. The benchmark unpacks it before the
# production evaluator sees it, preserving that evaluator's interface.
packed_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL) {
  sql = hybrid_moments_sql(conn,source,fe,value_terms,global_terms,limit)
  if(!length(global_terms)) return(sql)
  alias = function(terms) sub('.* AS ([A-Za-z_0-9]+)$','\\1',terms)
  fields = alias(global_terms)
  json_fields = paste(sprintf("'%s',%s",fields,fields),collapse=',')
  paste0('SELECT __set,__i,__j,',paste(alias(value_terms),collapse=','),
    ',CASE WHEN __set=0 THEN to_json(named_struct(',json_fields,')) END AS __globals FROM (',sql,') __dbglm_packed')
}

unpack_global_moments = function(result) {
  if(!'__globals' %in% names(result)) return(result)
  global_row = which(as.integer(result[['__set']])==0L)
  values = jsonlite::fromJSON(result[['__globals']][global_row])
  for(name in names(values)) {
    result[[name]] = rep(NA_real_,nrow(result))
    result[[name]][global_row] = as.numeric(values[[name]])
  }
  result[['__globals']] = NULL
  result
}
