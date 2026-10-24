# Own only fit-created resources; caller connections and sample exports survive.
dbglm_resources = function() {
  resources = new.env(parent = emptyenv())
  resources[["conn"]] = NULL
  resources[["own_conn"]] = FALSE
  resources[["registered_table"]] = NULL
  resources[["tables"]] = character()
  resources[["views"]] = character()
  resources
}

dbglm_release_resources = function(resources) {
  if (!is.null(resources[["conn"]])) {
    for (name in rev(resources[["views"]])) {
      try(DBI::dbExecute(resources[["conn"]], paste("DROP VIEW IF EXISTS", DBI::dbQuoteIdentifier(resources[["conn"]], name))), silent = TRUE)
    }
    for (name in rev(resources[["tables"]])) {
      try(drop_ap_table_if_exists(resources[["conn"]], as.character(DBI::dbQuoteIdentifier(resources[["conn"]], name)),
        detect_backend(resources[["conn"]])), silent = TRUE)
    }
    if (!is.null(resources[["registered_table"]])) {
      try(duckdb_unregister(resources[["conn"]], resources[["registered_table"]]), silent = TRUE)
    }
    if (isTRUE(resources[["own_conn"]])) try(DBI::dbDisconnect(resources[["conn"]], shutdown = TRUE), silent = TRUE)
  }
  invisible(NULL)
}

# Databricks inline expressions can still reference earlier physical tables.
# Defer their removal until fit cleanup; DuckDB can release them immediately.
dbglm_table_manager = function(conn, resources, backend_info = detect_backend(conn)) {
  quote_id = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  create = function(sql, materialize = TRUE, target = NULL) {
    view = isTRUE(backend_info[["is_databricks"]]) && !materialize &&
      !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))
    if (view) return(paste0("(", dbglm_sql(conn, sql), ")"))
    if (!is.null(target) && isTRUE(backend_info[["is_databricks"]])) {
      DBI::dbExecute(conn, dbglm_sql(conn, paste("CREATE OR REPLACE TEMP TABLE", target, "AS", sql)))
      return(target)
    }
    name = basename(tempfile("dbglm_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, dbglm_sql(conn, paste("CREATE TEMP TABLE", quote_id(name), "AS", sql)))
    quote_id(name)
  }
  drop = function(table) {
    if (is.null(table) || startsWith(table, "(")) return(invisible(NULL))
    name = gsub('["`]', "", table)
    if (isTRUE(backend_info[["is_databricks"]])) return(invisible(NULL))
    if (name %in% resources[["views"]]) {
      DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", table))
      resources[["views"]] = setdiff(resources[["views"]], name)
    } else {
      drop_ap_table_if_exists(conn, table, backend_info)
      resources[["tables"]] = setdiff(resources[["tables"]], name)
    }
    invisible(NULL)
  }
  list(create = create, drop = drop)
}

# One coefficient-frame/lookup contract serves prediction and projection.
# Decode integer64 keys before unlist() can strip their class.
dbglm_coefficient_frame = function(info, coefficients) {
  coefficients = as.matrix(coefficients)
  frame = data.frame(op = rep(seq_along(info), vapply(info, function(group) length(group[["id"]]), integer(1))))
  names(frame) = "__glm_op"
  frame[["__glm_group"]] = unlist(lapply(info, function(group) as.numeric(group[["id"]])), use.names = FALSE)
  frame[paste0("c", seq_len(ncol(coefficients)))] = as.data.frame(coefficients)
  frame
}

# `tables` lists physical objects owned by the fit; inline SQL owns nothing.
dbglm_coefficient_lookup = function(conn, frame, create, inline = FALSE) {
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    sql = dbglm_values_sql(conn, frame)
    if (is.null(sql)) stop("FE coefficient lookup exceeds the SQL transfer limit.", call. = FALSE)
    if (inline) return(list(source = paste0("(", sql, ")"), tables = character()))
    lookup = create(sql, materialize = TRUE)
  } else {
    registered = basename(tempfile("dbglm_coef_"))
    tryCatch(
      {
        duckdb::duckdb_register(conn, registered, frame)
        lookup = create(paste("SELECT * FROM", DBI::dbQuoteIdentifier(conn, registered)), materialize = TRUE)
      },
      finally = try(duckdb::duckdb_unregister(conn, registered), silent = TRUE))
  }
  list(source = lookup, tables = lookup)
}

# Define a coefficient lookup once, then filter it per FE operation. Repeating
# an inline JSON lookup in every operation can exceed the warehouse query-text
# limit even when the lookup itself is within its transfer bound.
dbglm_coefficient_ctes = function(source, n_ops, prefix) {
  inline = startsWith(source, "(")
  lookup = if (inline) paste0(prefix, "_lookup") else source
  c(if (inline) paste0(lookup, " AS (SELECT * FROM ", source, ")"),
    vapply(seq_len(n_ops), function(i) paste0(prefix, "_g", i,
      " AS (SELECT * FROM ", lookup, " WHERE __glm_op=", i, ")"), character(1)))
}
