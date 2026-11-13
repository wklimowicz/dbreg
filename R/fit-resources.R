# This tracker contains only owned database resources, never a fitting frame.
fit_resources = function() {
  resources = new.env(parent = emptyenv())
  resources[["conn"]] = NULL
  resources[["own_conn"]] = FALSE
  resources[["registered_table"]] = NULL
  resources[["tables"]] = resources[["views"]] = character()
  resources[["cleanup_errors"]] = character()
  resources
}

.dbreg_retired_resources = new.env(parent = emptyenv())

resource_connection_valid = function(conn) {
  !is.null(conn) && isTRUE(tryCatch(DBI::dbIsValid(conn), error = function(e) FALSE))
}

resource_same_connection = function(a, b) {
  for (slot in c("conn_ref", "ptr")) {
    pointer = attr(a, slot, exact = TRUE)
    if (typeof(pointer) == "externalptr" && identical(pointer, attr(b, slot, exact = TRUE))) return(TRUE)
  }
  attr(a, "dbreg_resource_lease") = attr(b, "dbreg_resource_lease") = NULL
  identical(a, b)
}

resource_pending = function(resources) {
  isTRUE(resources[["own_conn"]]) || !is.null(resources[["registered_table"]]) ||
    length(resources[["tables"]]) > 0L || length(resources[["views"]]) > 0L
}

resource_clear = function(resources) {
  resources[["conn"]] = resources[["registered_table"]] = NULL
  resources[["own_conn"]] = FALSE
  resources[["tables"]] = resources[["views"]] = resources[["cleanup_errors"]] = character()
}

drop_owned_relation = function(conn, name, kind = "tables") {
  quoted = as.character(DBI::dbQuoteIdentifier(conn, name))
  command = if (kind == "views") "DROP VIEW IF EXISTS" else
    if (isTRUE(detect_backend(conn)[["is_databricks"]])) "DROP TEMP TABLE IF EXISTS" else "DROP TABLE IF EXISTS"
  DBI::dbExecute(conn, paste(command, quoted))
}

release_fit_resources = function(resources, quiet = TRUE, retire = TRUE) {
  if (is.null(resources)) return(invisible(NULL))
  errors = character()
  if (resource_pending(resources)) {
    conn = resources[["conn"]]
    if (!resource_connection_valid(conn)) {
      resource_clear(resources)
    } else {
      for (kind in c("views", "tables")) for (name in rev(resources[[kind]])) {
        error = tryCatch({drop_owned_relation(conn, name, kind); NULL}, error = identity)
        if (is.null(error)) resources[[kind]] = setdiff(resources[[kind]], name) else errors = c(errors, conditionMessage(error))
      }
      if (!is.null(resources[["registered_table"]])) {
        error = tryCatch({duckdb::duckdb_unregister(conn, resources[["registered_table"]]); NULL}, error = identity)
        if (is.null(error)) resources[["registered_table"]] = NULL else errors = c(errors, conditionMessage(error))
      }
      if (isTRUE(resources[["own_conn"]])) {
        error = tryCatch({DBI::dbDisconnect(conn, shutdown = TRUE); NULL}, error = identity)
        if (is.null(error)) {
          resource_clear(resources)
          errors = character()
        } else errors = c(errors, conditionMessage(error))
      }
    }
  }
  resources[["cleanup_errors"]] = errors
  key = resources[["retired_key"]]
  if (isTRUE(retire)) {
    if (resource_pending(resources)) {
      if (is.null(key)) key = basename(tempfile("retired_"))
      resources[["retired_key"]] = key
      .dbreg_retired_resources[[key]] = resources
    } else if (!is.null(key)) {
      if (exists(key, envir = .dbreg_retired_resources, inherits = FALSE)) rm(list = key, envir = .dbreg_retired_resources)
      resources[["retired_key"]] = NULL
    }
  }
  if (length(errors) && !quiet) stop("Could not release fit resources: ", paste(errors, collapse = "; "), call. = FALSE)
  invisible(NULL)
}

# Cleanup retries happen synchronously, before another source is allocated.
# No driver calls or connection closes are performed by garbage collection.
flush_retired_resources = function(conn = NULL, quiet = FALSE) {
  errors = character()
  for (key in ls(.dbreg_retired_resources, all.names = TRUE)) {
    resources = .dbreg_retired_resources[[key]]
    if (!isTRUE(resources[["own_conn"]]) && resource_connection_valid(resources[["conn"]]) &&
        (is.null(conn) || !resource_same_connection(conn, resources[["conn"]]))) next
    rm(list = key, envir = .dbreg_retired_resources)
    release_fit_resources(resources, quiet = TRUE, retire = FALSE)
    if (resource_pending(resources)) .dbreg_retired_resources[[key]] = resources else resources[["retired_key"]] = NULL
    errors = c(errors, resources[["cleanup_errors"]])
  }
  if (length(errors) && !quiet) stop("Could not release retired fit resources: ", paste(errors, collapse = "; "), call. = FALSE)
  invisible(NULL)
}
