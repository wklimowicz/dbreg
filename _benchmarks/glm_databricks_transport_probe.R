# Diagnose truncation without changing the user's installed packages.
pkgload::load_all()
out = commandArgs(TRUE)[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows = list()
for (variant in c("numeric_large", "string_large", "string_bounded")) {
  setting = switch(variant, numeric_large=1000000, string_large="1000000", string_bounded="65536")
  conn = DBI::dbConnect(odbc::databricks(),
    driver=Sys.getenv("DATABRICKS_ODBC_DRIVER","/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath=Sys.getenv("DATABRICKS_SQL_PATH"),useNativeQuery=TRUE,DefaultStringColumnLength=setting)
  DBI::dbExecute(conn,"SET use_cached_result=false")
  for (size in c(512L,4096L,60000L,100000L)) {
    sql = paste0("SELECT REPEAT('x',",size,") AS txt,CAST(REPEAT('x',",size,") AS BINARY) AS bin")
    data = DBI::dbGetQuery(conn,sql)
    row = data.frame(variant=variant,setting=as.character(setting),expected=size,
      text_bytes=nchar(data[["txt"]][[1]],type="bytes"),binary_bytes=length(data[["bin"]][[1]]))
    print(row);flush.console()
    rows[[length(rows)+1L]]=row
    write.csv(do.call(rbind,rows),file.path(out,"probe.csv"),row.names=FALSE)
  }
  DBI::dbDisconnect(conn)
}
