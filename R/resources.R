# Own only fit-created resources; caller connections and sample exports survive.
fit_resources = function() {
  resources = new.env(parent = emptyenv())
  resources[["conn"]] = NULL
  resources[["own_conn"]] = FALSE
  resources[["registered_table"]] = NULL
  resources[["tables"]] = character()
  resources[["views"]] = character()
  resources
}

release_fit_resources = function(resources, quiet = TRUE) {
  errors = character()
  if (!is.null(resources[["conn"]])) {
    for (kind in c("views", "tables")) for (name in rev(resources[[kind]])) {
      error = tryCatch({ drop_owned_relation(resources[["conn"]], name, kind); NULL }, error = identity)
      if (is.null(error)) resources[[kind]] = setdiff(resources[[kind]], name) else
        errors = c(errors, conditionMessage(error))
    }
    if (!is.null(resources[["registered_table"]])) {
      error = tryCatch({ duckdb_unregister(resources[["conn"]], resources[["registered_table"]]); NULL }, error = identity)
      if (is.null(error)) resources[["registered_table"]] = NULL else errors = c(errors, conditionMessage(error))
    }
    if (isTRUE(resources[["own_conn"]])) {
      error = tryCatch({ DBI::dbDisconnect(resources[["conn"]], shutdown = TRUE); NULL }, error = identity)
      if (is.null(error)) {
        resources[["own_conn"]] = FALSE
        resources[["conn"]] = NULL
        resources[["tables"]] = resources[["views"]] = character()
        resources[["registered_table"]] = NULL
        errors = character()
      } else errors = c(errors, conditionMessage(error))
    }
  }
  resources[["cleanup_errors"]] = errors
  if (length(errors) && !quiet) stop("Could not release fit resources: ", paste(errors, collapse = "; "), call. = FALSE)
  invisible(NULL)
}

# Databricks inline expressions can still reference earlier physical tables.
# Defer their removal until fit cleanup; DuckDB can release them immediately.
fit_table_manager = function(conn, resources, backend_info = detect_backend(conn)) {
  quote_id = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  create = function(sql, materialize = TRUE, target = NULL) {
    view = isTRUE(backend_info[["is_databricks"]]) && !materialize &&
      !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))
    if (view) return(paste0("(", backend_sql(conn, sql), ")"))
    if (!is.null(target) && isTRUE(backend_info[["is_databricks"]])) {
      DBI::dbExecute(conn, backend_sql(conn, paste("CREATE OR REPLACE TEMP TABLE", target, "AS", sql)))
      return(target)
    }
    name = basename(tempfile("dbglm_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, backend_sql(conn, paste("CREATE TEMP TABLE", quote_id(name), "AS", sql)))
    quote_id(name)
  }
  drop = function(table) {
    if (is.null(table) || startsWith(table, "(")) return(invisible(NULL))
    name = gsub('["`]', "", table)
    if (isTRUE(backend_info[["is_databricks"]])) return(invisible(NULL))
    if (name %in% resources[["views"]]) {
      drop_owned_relation(conn, name, "views")
      resources[["views"]] = setdiff(resources[["views"]], name)
    } else {
      drop_owned_relation(conn, name)
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
    sql = sql_values_frame(conn, frame)
    if (is.null(sql)) projection_abort("upload_limit", "FE coefficient lookup exceeds the SQL transfer limit.")
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
