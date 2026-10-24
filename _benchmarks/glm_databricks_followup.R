# Targeted experiments on the exact saved production queries and synthetic data.
# No package mutations. Compare transport and a raw-key pruning preflight.
# Rscript _benchmarks/glm_databricks_followup.R /tmp/teacher-cloglog-production /tmp/glm-followup
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("_benchmarks/teacher_cloglog_profile.R")
args = commandArgs(TRUE)
out = args[[2]]
dir.create(file.path(out, "sql"), recursive = TRUE, showWarnings = FALSE)
log = new.env(parent = emptyenv())
log[["out"]] = out; log[["id"]] = 0L
log[["tag"]] = paste0("glm_followup_", Sys.getpid())
log[["run"]] = "setup"; log[["queries"]] = list()
saved = function(i) paste(readLines(file.path(args[[1]], "sql",
  paste0("teacher_integrated_174818_Databricks_integrated_1_", i, ".sql"))), collapse = "\n")
run = function() {
  native = DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE,
    DefaultStringColumnLength = 10000000)
  on.exit(DBI::dbDisconnect(native))
  conn = methods::new("teacher_profile_connection", native = native, log = log)
  DBI::dbExecute(native, "SET use_cached_result=false")
  DBI::dbExecute(native, paste("CREATE TEMP TABLE teacher_integrated_input_174818 AS", teacher_cloglog_sql(5e6)))
  DBI::dbExecute(conn, saved(5))
  prune = paste0("SELECT COUNT(*) AS groups_to_remove FROM (SELECT COUNT(*) AS n, MIN(wastage) AS lo,MAX(wastage) AS hi ",
    "FROM teacher_integrated_input_174818 GROUP BY GROUPING SETS ",
    "((spine_band_year),(district),(Experience),(AGE),(SEX),(FTPTWK))) g WHERE n=1 OR lo=hi")
  # Transport is isolated from mathematical aggregation. Each packed row carries
  # its keys with its values; collect_list order is deliberately irrelevant.
  eval = saved(10)
  fields = c("__i", "__j", "__h", "__b1", "__b2", "__ex1")
  packed = paste0("SELECT __set,COUNT(*) AS __n,MAX(__globals) AS __globals,",
    "to_json(collect_list(array(", paste(fields, collapse=","), "))) AS __payload FROM (",
    eval, ") p GROUP BY __set,pmod(COALESCE(__i,0),256)")
  rows = list(); reference = NULL
  for (rep in 1:2) for (variant in c("prune_original", "prune_preflight", "moments_original", "moments_packed")) {
    log[["run"]] = paste0(variant, "_", rep)
    sql = switch(variant, prune_original=saved(3), prune_preflight=prune,
      moments_original=eval, moments_packed=packed)
    gc()
    elapsed = system.time({res = DBI::dbGetQuery(conn, sql)})[["elapsed"]]
    decode = 0; status = "passed"
    if (variant == "moments_original") reference = res
    if (variant == "moments_packed") {
      decode = system.time({error = tryCatch({
        frames = lapply(seq_len(nrow(res)), function(i) {
          values = jsonlite::fromJSON(res[["__payload"]][[i]])
          stopifnot(nrow(values) == as.numeric(res[["__n"]][[i]]))
          frame = as.data.frame(values); names(frame) = fields
          frame[["__set"]] = res[["__set"]][[i]]
          frame
        })
        decoded = do.call(rbind, frames)
        ordered = function(d) d[order(d[["__set"]],d[["__i"]],d[["__j"]],na.last=TRUE), ]
        a = ordered(reference); b = ordered(decoded)
        stopifnot(nrow(a)==nrow(b))
        for (field in c("__set",fields)) {
          x = as.numeric(a[[field]]); y = as.numeric(b[[field]])
          stopifnot(identical(is.na(x),is.na(y)), isTRUE(all.equal(x,y,tolerance=1e-12)))
        }
        stopifnot(isTRUE(all.equal(jsonlite::fromJSON(reference[["__globals"]][reference[["__set"]]==0]),
          jsonlite::fromJSON(res[["__globals"]][res[["__set"]]==0]), tolerance=1e-12)))
        NULL
      }, error=function(e) conditionMessage(e))})[["elapsed"]]
      if (!is.null(error)) status = paste("rejected:", error)
    }
    rows[[length(rows)+1L]] = data.frame(variant=variant,repetition=rep,query_seconds=elapsed,
      decode_and_check_seconds=decode,rows=nrow(res),result_bytes=as.numeric(object.size(res)),status=status)
    print(utils::tail(rows,1)[[1]]); flush.console()
    write.csv(do.call(rbind,rows),file.path(out,"experiments.csv"),row.names=FALSE)
  }
}
run()
