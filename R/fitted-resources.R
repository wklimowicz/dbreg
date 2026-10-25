# Mathematical fitted state stays in plain lists. Only this separate lease
# owns database handles; its finalizer must never capture a fitting frame.
drop_owned_relation = function(conn, name, kind = "tables") {
  quoted = as.character(DBI::dbQuoteIdentifier(conn, name))
  if (kind == "views") return(DBI::dbExecute(conn, paste("DROP VIEW IF EXISTS", quoted)))
  backend = detect_backend(conn)
  DBI::dbExecute(conn, paste(if (isTRUE(backend[["is_databricks"]])) "DROP TEMP TABLE IF EXISTS" else
    "DROP TABLE IF EXISTS", quoted))
}

fitted_resource_valid_connection = function(conn) {
  !is.null(conn) && isTRUE(tryCatch(DBI::dbIsValid(conn), error = function(e) FALSE))
}

# Abandoned leases on caller sessions are retired without calling a driver.
# Database work runs only at explicit, synchronous package entry points.
.dbreg_retired_resources = new.env(parent = emptyenv())

retire_fitted_resources = function(lease) {
  session = lease[["session"]]
  key = basename(tempfile("retired_"))
  .dbreg_retired_resources[[key]] = list(conn = session[["conn"]],
    tables = lease[["tables"]], views = lease[["views"]])
  lease[["tables"]] = lease[["views"]] = character()
  lease[["released"]] = TRUE
  session[["references"]] = session[["references"]] - 1L
  if (session[["references"]] == 0L) {
    session[["closed"]] = TRUE
    session[["conn"]] = NULL
  }
  invisible(NULL)
}

flush_retired_resources = function(conn, quiet = TRUE) {
  errors = character()
  for (key in ls(.dbreg_retired_resources, all.names = TRUE)) {
    entry = .dbreg_retired_resources[[key]]
    if (!fitted_resource_valid_connection(entry[["conn"]])) {
      rm(list = key, envir = .dbreg_retired_resources)
      next
    }
    if (!fitted_resource_same_connection(conn, entry[["conn"]])) next
    # Detach before driver calls: GC may append new entries during a query.
    rm(list = key, envir = .dbreg_retired_resources)
    for (kind in c("views", "tables")) for (name in rev(entry[[kind]])) {
      error = tryCatch({ drop_owned_relation(conn, name, kind); NULL }, error = identity)
      if (is.null(error)) entry[[kind]] = setdiff(entry[[kind]], name) else
        errors = c(errors, conditionMessage(error))
    }
    if (length(c(entry[["views"]], entry[["tables"]]))) .dbreg_retired_resources[[key]] = entry
  }
  if (length(errors) && !quiet) stop("Could not release retired database state: ",
    paste(errors, collapse = "; "), call. = FALSE)
  invisible(NULL)
}

fitted_resource_session = function(conn, own_conn = FALSE) {
  flush_retired_resources(conn)
  if (isTRUE(own_conn) && identical(detect_backend(conn)[["name"]], "duckdb")) {
    DBI::dbExecute(conn, "SET threads = 1")
  }
  inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) {
    fitted_resource_connection(inherited)
    return(inherited[["session"]])
  }
  session = new.env(parent = emptyenv())
  session[["conn"]] = conn
  session[["own_conn"]] = isTRUE(own_conn)
  session[["references"]] = 0L
  session[["closed"]] = FALSE
  session
}

# A public connection handle keeps its result lease alive while used as input
# to another operation. The shared session itself keeps an untagged handle, so
# unrelated leases do not prevent garbage collection of abandoned results.
fitted_resource_public_connection = function(lease) {
  conn = fitted_resource_connection(lease)
  attr(conn, "dbreg_resource_lease") = lease
  conn
}

fitted_resource_same_connection = function(a, b) {
  if (!fitted_resource_valid_connection(a) || !fitted_resource_valid_connection(b)) return(FALSE)
  for (slot in c("conn_ref", "ptr")) {
    pointer = attr(a, slot, exact = TRUE)
    if (typeof(pointer) == "externalptr" && identical(pointer, attr(b, slot, exact = TRUE))) return(TRUE)
  }
  attr(a, "dbreg_resource_lease") = attr(b, "dbreg_resource_lease") = NULL
  identical(a, b)
}

fitted_resource_lease = function(session, tables = character(), views = character()) {
  if (isTRUE(session[["closed"]]) || !fitted_resource_valid_connection(session[["conn"]])) {
    stop("Retained database state requires its live owning session.", call. = FALSE)
  }
  lease = new.env(parent = emptyenv())
  lease[["version"]] = 1L
  lease[["session"]] = session
  lease[["tables"]] = unique(tables)
  lease[["views"]] = unique(views)
  lease[["released"]] = FALSE
  session[["references"]] = session[["references"]] + 1L
  reg.finalizer(lease, fitted_resource_finalize, onexit = TRUE)
  lease
}

# Transfer only explicitly tracked objects. The fit still drops all of its
# other intermediates and unregisters any in-memory source on normal cleanup.
retain_fit_resources = function(resources, tables, views = character()) {
  if (any(!tables %in% resources[["tables"]]) || any(!views %in% resources[["views"]])) {
    stop("Cannot retain relations that are not owned by this fit.", call. = FALSE)
  }
  session = resources[["retained_session"]]
  if (is.null(session)) {
    session = fitted_resource_session(resources[["conn"]], resources[["own_conn"]])
  }
  lease = fitted_resource_lease(session, tables, views)
  resources[["retained_session"]] = session
  resources[["own_conn"]] = FALSE
  resources[["tables"]] = setdiff(resources[["tables"]], tables)
  resources[["views"]] = setdiff(resources[["views"]], views)
  lease
}

fitted_resource_connection = function(lease) {
  if (!is.environment(lease) || !identical(lease[["version"]], 1L)) {
    stop("Unsupported retained-resource version.", call. = FALSE)
  }
  if (isTRUE(lease[["released"]])) stop("Retained database state has been released.", call. = FALSE)
  session = lease[["session"]]
  if (isTRUE(session[["closed"]]) || !fitted_resource_valid_connection(session[["conn"]])) {
    stop("Retained database state requires its live owning session; reconnecting or deserializing a handle cannot restore it.", call. = FALSE)
  }
  flush_retired_resources(session[["conn"]])
  session[["conn"]]
}

fitted_resource_release = function(lease, quiet = FALSE) {
  if (is.null(lease)) return(invisible(NULL))
  if (!is.environment(lease) || !identical(lease[["version"]], 1L)) {
    stop("Unsupported retained-resource version.", call. = FALSE)
  }
  if (isTRUE(lease[["released"]])) return(invisible(NULL))
  session = lease[["session"]]
  conn = session[["conn"]]
  flush_retired_resources(conn, quiet = quiet)
  errors = character()
  if (fitted_resource_valid_connection(conn)) {
    backend = detect_backend(conn)
    for (kind in c("views", "tables")) {
      for (name in rev(lease[[kind]])) {
        error = tryCatch({
          drop_owned_relation(conn, name, kind)
          NULL
        }, error = identity)
        if (is.null(error)) lease[[kind]] = setdiff(lease[[kind]], name) else
          errors = c(errors, conditionMessage(error))
      }
    }
  } else {
    lease[["tables"]] = lease[["views"]] = character()
  }
  # Keep failed drops owned and retryable; successful drops are never retried.
  if (length(errors)) {
    if (!quiet) stop("Could not release retained database state: ", paste(errors, collapse = "; "), call. = FALSE)
    return(invisible(NULL))
  }
  if (session[["references"]] == 1L && isTRUE(session[["own_conn"]]) && fitted_resource_valid_connection(conn)) {
    error = tryCatch({ DBI::dbDisconnect(conn, shutdown = TRUE); NULL }, error = identity)
    if (!is.null(error)) {
      if (!quiet) stop("Could not close the retained database session: ", conditionMessage(error), call. = FALSE)
      return(invisible(NULL))
    }
  }
  lease[["released"]] = TRUE
  session[["references"]] = session[["references"]] - 1L
  if (session[["references"]] == 0L) {
    session[["closed"]] = TRUE
    session[["conn"]] = NULL
  }
  invisible(NULL)
}

fitted_resource_finalize = function(lease) {
  if (!is.environment(lease) || isTRUE(lease[["released"]])) return(invisible(NULL))
  session = lease[["session"]]
  if (isTRUE(session[["own_conn"]]) && session[["references"]] == 1L) {
    try(fitted_resource_release(lease, quiet = TRUE), silent = TRUE)
  } else {
    retire_fitted_resources(lease)
  }
  invisible(NULL)
}

#' Release retained database state
#'
#' Releases only relations owned by this model. Caller connections and explicit
#' sample/prediction exports survive. Portable fitted state requires no release.
#' Calling `close()` repeatedly is safe. Copies of a model share their resource
#' lease, so closing one also releases the database state used by those copies.
#' A model-owned connection closes when its last model or prediction lease ends.
#' Abandoned state on a caller-owned connection is queued without executing SQL
#' during garbage collection. It is released at the next synchronous package
#' operation on that connection; call `close()` for immediate cleanup. Private
#' retained DuckDB sessions use one worker thread and close when abandoned.
#' Live database state cannot be restored by serializing its connection handle.
#'
#' @param con A fitted `dbreg` or `dbglm` model.
#' @param ... Additional arguments are not supported.
#' @return Invisibly returns `NULL`.
#' @export
close.dbreg = function(con, ...) {
  check_unused_args(list(...))
  fitted_resource_release(con[["resource_state"]])
}
