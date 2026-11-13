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
    length(resources[["tables"]]) > 0L || length(resources[["views"]]) > 0L ||
    !is.null(resources[["session"]])
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
  session = resources[["session"]]
  if (!is.null(session) && !length(resources[["tables"]]) &&
      !length(resources[["views"]]) && is.null(resources[["registered_table"]])) {
    resources[["session"]] = NULL
    session[["leases"]] = session[["leases"]] - 1L
    if (session[["leases"]] == 0L && session[["own_conn"]] &&
        resource_connection_valid(session[["conn"]])) {
      error = tryCatch({DBI::dbDisconnect(session[["conn"]], shutdown = TRUE); NULL},
        error = identity)
      if (!is.null(error)) {
        resources[["conn"]] = session[["conn"]]
        resources[["own_conn"]] = TRUE
        errors = c(errors, conditionMessage(error))
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
    owns_session = isTRUE(resources[["session"]][["own_conn"]])
    if (!isTRUE(resources[["own_conn"]]) && !owns_session && resource_connection_valid(resources[["conn"]]) &&
        (is.null(conn) || !resource_same_connection(conn, resources[["conn"]]))) next
    rm(list = key, envir = .dbreg_retired_resources)
    release_fit_resources(resources, quiet = TRUE, retire = FALSE)
    if (resource_pending(resources)) .dbreg_retired_resources[[key]] = resources else resources[["retired_key"]] = NULL
    errors = c(errors, resources[["cleanup_errors"]])
  }
  if (length(errors) && !quiet) stop("Could not release retired fit resources: ", paste(errors, collapse = "; "), call. = FALSE)
  invisible(NULL)
}

# A lease takes ownership of selected physical tables and, when applicable,
# the fit-created connection. The temporary fit tracker still owns all other
# tables and the R registration.
retain_fit_resources = function(resources, tables, session = NULL) {
  if (!all(tables %in% resources[["tables"]])) {
    stop("Cannot retain tables the fit does not own.", call. = FALSE)
  }
  lease = fit_resources()
  lease[["conn"]] = resources[["conn"]]
  if (is.null(session)) {
    session = list2env(list(conn = resources[["conn"]], own_conn = resources[["own_conn"]],
      leases = 0L), parent = emptyenv())
  }
  session[["leases"]] = session[["leases"]] + 1L
  lease[["session"]] = session
  lease[["tables"]] = tables
  lease[["closed"]] = FALSE
  resources[["tables"]] = setdiff(resources[["tables"]], tables)
  resources[["own_conn"]] = FALSE
  reg.finalizer(lease, function(state) {
    if (resource_pending(state)) {
      key = state[["retired_key"]]
      if (is.null(key)) key = basename(tempfile("retired_"))
      state[["retired_key"]] = key
      .dbreg_retired_resources[[key]] = state
    }
  }, onexit = FALSE)
  lease
}

retained_fit_connection = function(lease) {
  if (is.null(lease) || isTRUE(lease[["closed"]]) ||
      !resource_connection_valid(lease[["conn"]])) {
    stop("The retained fit has no live database session.", call. = FALSE)
  }
  lease[["conn"]]
}

#' Release tables owned by a fitted database regression
#'
#' @param con A fitted database regression object.
#' @param ... Unused.
#' @return Invisibly, NULL.
#' @export
close.dbreg = function(con, ...) close_retained_lease(con[["resource_state"]])

close_retained_lease = function(lease) {
  if (!is.null(lease)) {
    lease[["closed"]] = TRUE
    release_fit_resources(lease, quiet = FALSE)
  }
  invisible(NULL)
}
