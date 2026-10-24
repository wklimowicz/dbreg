# Fetch sanitized server timings for the tagged queries from one benchmark run:
# Rscript _benchmarks/glm_databricks_history.R /tmp/glm-perf
# This only reads query history; tokens, hostnames and unrelated query text are
# never written to the result files. Run within three hours of the benchmark.
args = commandArgs(TRUE)
out = args[[1]]
q = read.csv(file.path(out,'queries.csv'))
host = sub('/+$','',Sys.getenv('DATABRICKS_HOST'))
if(!startsWith(host,'http')) host = paste0('https://',host)
req = httr2::request(paste0(host,'/api/2.0/sql/history/queries'))
req = httr2::req_auth_bearer_token(req,Sys.getenv('DATABRICKS_TOKEN'))
req = httr2::req_url_query(req,max_results=999,include_metrics='true',
  'filter_by.query_start_time_range.start_time_ms'=format((as.numeric(Sys.time())-10800)*1000,scientific=FALSE,trim=TRUE))
body = httr2::resp_body_json(httr2::req_perform(req),simplifyVector=FALSE)
rows = lapply(body[['res']],function(r) {
  tag = sub('^/\\* ([A-Za-z0-9_]+) \\*/.*','\\1',r[['query_text']])
  if(!tag %in% q[['tag']]) return(NULL)
  fields = c('compilation_time_ms','execution_time_ms','result_fetch_time_ms','total_time_ms','read_bytes','rows_read_count','read_cache_bytes','read_remote_bytes','write_remote_bytes','network_sent_bytes','result_from_cache','photon_total_time_ms','task_total_time_ms','spill_to_disk_bytes')
  data = lapply(fields,function(n) if(is.null(r[['metrics']][[n]])) NA else r[['metrics']][[n]])
  names(data) = fields
  cbind(data.frame(tag=tag,status=r[['status']]),as.data.frame(data))
})
rows = do.call(rbind,rows)
if(!is.null(rows)) {
  merged = merge(q,rows,by='tag',all.x=TRUE)
  write.csv(merged,file.path(out,'server_queries.csv'),row.names=FALSE)
  print(aggregate(cbind(elapsed,compilation_time_ms,execution_time_ms,result_fetch_time_ms,rows_read_count)~run+phase,merged,sum))
} else cat('No matching query metrics yet\n')
