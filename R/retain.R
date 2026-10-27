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
  if (!is.null(resources[["conn"]]) && !DBI::dbIsValid(resources[["conn"]])) {
    resources[["conn"]] = NULL
    resources[["own_conn"]] = FALSE
    resources[["tables"]] = resources[["views"]] = character()
    resources[["registered_table"]] = NULL
  }
  if (!is.null(resources[["conn"]])) {
    for (kind in c("views", "tables")) for (name in rev(resources[[kind]])) {
      error = tryCatch(
        {
          drop_owned_relation(resources[["conn"]], name, kind)
          NULL
        },
        error = identity
      )
      if (is.null(error)) resources[[kind]] = setdiff(resources[[kind]], name) else
        errors = c(errors, conditionMessage(error))
    }
    if (!is.null(resources[["registered_table"]])) {
      error = tryCatch(
        {
          duckdb_unregister(resources[["conn"]], resources[["registered_table"]])
          NULL
        },
        error = identity
      )
      if (is.null(error)) resources[["registered_table"]] = NULL else errors = c(errors, conditionMessage(error))
    }
    if (isTRUE(resources[["own_conn"]])) {
      error = tryCatch(
        {
          DBI::dbDisconnect(resources[["conn"]], shutdown = TRUE)
          NULL
        },
        error = identity
      )
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
    if (view) return(paste0("(", sql, ")"))
    if (identical(backend_info[["name"]], "spark")) validate_databricks_ap_backend(backend_info)
    if (!is.null(target) && isTRUE(backend_info[["is_databricks"]])) {
      DBI::dbExecute(conn, paste("CREATE OR REPLACE TEMP TABLE", target, "AS", sql))
      return(target)
    }
    name = basename(tempfile("dbglm_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", quote_id(name), "AS", sql))
    quote_id(name)
  }
  drop = function(table) {
    if (is.null(table) || startsWith(table, "(")) return(invisible(NULL))
    name = gsub("[\"`]", "", table)
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
  .dbreg_retired_resources[[key]] = list(
    conn = session[["conn"]],
    tables = lease[["tables"]], views = lease[["views"]]
  )
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
      error = tryCatch(
        {
          drop_owned_relation(conn, name, kind)
          NULL
        },
        error = identity
      )
      if (is.null(error)) entry[[kind]] = setdiff(entry[[kind]], name) else
        errors = c(errors, conditionMessage(error))
    }
    if (length(c(entry[["views"]], entry[["tables"]]))) .dbreg_retired_resources[[key]] = entry
  }
  if (length(errors) && !quiet) stop("Could not release retired database state: ",
    paste(errors, collapse = "; "),
    call. = FALSE
  )
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
        error = tryCatch(
          {
            drop_owned_relation(conn, name, kind)
            NULL
          },
          error = identity
        )
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
    error = tryCatch(
      {
        DBI::dbDisconnect(conn, shutdown = TRUE)
        NULL
      },
      error = identity
    )
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

retain_prediction = function(retain) length(retain) == 1L && !is.na(retain) && retain %in% c("prediction", "all")

retain_sample = function(retain) length(retain) == 1L && !is.na(retain) && retain %in% c("inference", "all")

prepare_retained_data = function(data, table, retain) {
  r_input = retain_sample(retain) && is.null(table) && is.data.frame(data)
  key = NULL
  if (r_input) {
    data = as.data.frame(data)
    key = safe_identifiers("dbreg_row_key__", taken = names(data), fold_case = TRUE)
    data[[key]] = seq_len(nrow(data))
  }
  list(data = data, key = key, r_input = r_input)
}

materialize_retained_source = function(conn, from_statement, resources, key = NULL, r_input = FALSE) {
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Full fitted-state retention supports DuckDB and native Databricks.", call. = FALSE)
  }
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  schema = DBI::dbGetQuery(conn, paste("SELECT *", from_statement, "LIMIT 0"))
  columns = names(schema)
  generated = is.null(key)
  if (generated) key = safe_identifiers("dbreg_row_key__", taken = columns, fold_case = folds_identifier_case(conn))
  if (!generated && !key %in% columns) stop("The retained source is missing its positional key.", call. = FALSE)
  fields = q(columns)
  if (generated) fields = c(fields, paste(if (isTRUE(backend[["is_databricks"]])) "monotonically_increasing_id()" else
    "ROW_NUMBER() OVER ()", "AS", q(key)))
  name = basename(tempfile("dbreg_original_sample_"))
  resources[["conn"]] = conn
  resources[["tables"]] = c(resources[["tables"]], name)
  DBI::dbExecute(conn, sql_sample_table(conn, name, paste("SELECT", paste(fields, collapse = ","), from_statement), c(columns, if (generated) key)))
  n = as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", q(name)))[["n"]])
  classes = vapply(schema, function(value) class(value)[[1L]], character(1))
  if (generated) classes[[key]] = "integer64"
  list(table = name, key = key, r_input = r_input, columns = setdiff(columns, if (!generated) key), column_classes = classes, nobs_input = n)
}

sql_distinct_tuple_count = function(conn, columns) {
  columns = as.character(DBI::dbQuoteIdentifier(conn, columns))
  if (length(columns) == 1L) return(paste0("COUNT(DISTINCT ", columns, ")"))
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  paste0(
    "COUNT(DISTINCT ", if (!native) "(" else "",
    paste(columns, collapse = ","), if (!native) ")" else "", ")"
  )
}

# Append scores to the immutable original rows while the fit projector lives.
# Only this final augmented relation is retained, including excluded positions.
retain_score_sample = function(
  conn, original, resources, score_query, score_key, scores,
  eta, residual, mu = NULL, row_id = NULL, fe_components = list(), weights = NULL, materialize = TRUE
) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  private = safe_identifiers(
    c(
      "dbreg_estimated__", "dbreg_eta__", "dbreg_residual__",
      if (!is.null(mu)) "dbreg_mu__", if (length(scores)) paste0("dbreg_score", seq_along(scores), "__")
    ),
    taken = c(original[["columns"]], original[["key"]]), fold_case = folds_identifier_case(conn)
  )
  estimated_column = private[[1L]]
  eta_column = private[[2L]]
  residual_column = private[[3L]]
  mu_column = if (!is.null(mu)) private[[4L]] else NULL
  score_columns = if (length(scores)) utils::tail(private, length(scores)) else character()
  classes = original[["column_classes"]]
  classes[private] = "numeric"
  classes[[estimated_column]] = "logical"
  select = c(
    "r.*", paste0("s.", q(score_key), " IS NOT NULL AS ", q(estimated_column)),
    paste(eta, "AS", q(eta_column)), paste(residual, "AS", q(residual_column)),
    if (!is.null(mu)) paste(mu, "AS", q(mu_column)), if (length(scores)) paste(scores, "AS", q(score_columns))
  )
  name = basename(tempfile("dbreg_fitted_sample_"))
  kind = if (materialize) "tables" else "views"
  resources[[kind]] = c(resources[[kind]], name)
  sql = paste0(
    "SELECT ", paste(select, collapse = ","), " FROM ", q(original[["table"]]), " r LEFT JOIN (",
    score_query, ") s ON r.", q(original[["key"]]), "=s.", q(score_key)
  )
  statement = if (materialize) sql_sample_table(
    conn, name, sql,
    c(original[["columns"]], original[["key"]], private)
  ) else
    paste("CREATE TEMP VIEW", q(name), "AS", sql)
  DBI::dbExecute(conn, statement)
  w = if (is.null(weights)) "1.0" else paste0("CAST(", q(weights), " AS DOUBLE)")
  summary = c(
    "COUNT(*) AS n_original", paste0("SUM(CASE WHEN ", q(estimated_column), " THEN 1 ELSE 0 END) AS n_fit"),
    paste0("SUM(CASE WHEN ", q(estimated_column), " THEN POWER(", q(residual_column), ",2) ELSE 0.0 END) AS rss_unweighted"),
    paste0("SUM(CASE WHEN ", q(estimated_column), " THEN ", w, "*POWER(", q(residual_column), ",2) ELSE 0.0 END) AS rss")
  )
  row = DBI::dbGetQuery(conn, paste("SELECT", paste(summary, collapse = ","), "FROM", q(name)))
  if (as.numeric(row[["n_original"]]) != original[["nobs_input"]]) stop("Retained score rows do not match the original sample.", call. = FALSE)
  fe_sizes = if (length(fe_components)) {
    terms = vapply(seq_along(fe_components), function(i) paste(sql_distinct_tuple_count(conn, fe_components[[i]]), "AS", paste0("g", i)), character(1))
    counts = DBI::dbGetQuery(conn, paste("SELECT", paste(terms, collapse = ","), "FROM", q(name), "WHERE", q(estimated_column)))
    vapply(counts, as.numeric, numeric(1))
  } else numeric()
  list(
    version = 1L, table = name, key = original[["key"]], r_input = original[["r_input"]],
    row_id = if (length(row_id)) row_id else original[["key"]], estimated = estimated_column,
    eta = eta_column, residual = residual_column, mu = mu_column, scores = score_columns,
    nobs_input = original[["nobs_input"]], nobs = as.numeric(row[["n_fit"]]), rss = as.numeric(row[["rss"]]),
    rss_unweighted = as.numeric(row[["rss_unweighted"]]), fe_components = fe_components, fe_sizes = fe_sizes,
    weights = weights,
    column_classes = classes,
    resources = list(
      storage = "database", tables = if (materialize) name else character(),
      views = if (!materialize) name else character(), version = 1L
    )
  )
}

retained_inference = function(object) {
  state = object[["fitted_state"]][["inference"]]
  if (is.null(state)) stop("This operation requires fitting with retain='all' or retain='inference'.", call. = FALSE)
  if (!identical(state[["version"]], 1L)) stop("Unsupported retained-inference version.", call. = FALSE)
  fitted_resource_connection(object[["resource_state"]])
  state
}

predict_retained_sample = function(object, type = "link", sample = c("estimation", "original")) {
  sample = match.arg(sample)
  state = retained_inference(object)
  conn = fitted_resource_connection(object[["resource_state"]])
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  from = paste("FROM", q(state[["table"]]), if (sample == "estimation") paste("WHERE", q(state[["estimated"]])) else "")
  if (isTRUE(state[["r_input"]])) {
    values = DBI::dbGetQuery(conn, paste("SELECT", q(state[["eta"]]), "AS eta", from, "ORDER BY", q(state[["key"]])))[["eta"]]
    return(if (inherits(object, "dbglm") && type == "response") object[["family"]][["linkinv"]](values) else values)
  }
  resources = fit_resources()
  resources[["conn"]] = conn
  resources[["retained_session"]] = object[["resource_state"]][["session"]]
  on.exit(release_fit_resources(resources))
  prediction = safe_identifiers("prediction", taken = state[["row_id"]], fold_case = folds_identifier_case(conn))
  create = function(sql) {
    name = basename(tempfile("dbreg_sample_prediction_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, sql_sample_table(conn, name, sql, c(state[["row_id"]], prediction)))
    name
  }
  table = create(paste("SELECT", paste(c(q(state[["row_id"]]), paste(q(state[["eta"]]), "AS", q(prediction))), collapse = ","), from))
  if (inherits(object, "dbglm") && type == "response") table = create(sql_prediction_response(conn, table, state[["row_id"]], prediction, object[["family"]][["link"]]))
  lease = retain_fit_resources(resources, table)
  structure(list(
    conn = fitted_resource_public_connection(lease), table = table, row_id = state[["row_id"]],
    prediction_column = prediction, resource_state = lease
  ), class = "dbreg_prediction")
}

inference_coefficients = function(object, V) {
  ids = attr(V, "dbreg_coefficient_ids")
  if (!is.null(ids) && !is.null(object[["linear_design"]])) return(linear_coefficients_for_ids(object, ids))
  beta = stats::coef(object)
  if (!identical(names(beta), rownames(V))) beta = beta[rownames(V)]
  beta
}

inference_coeftable = function(object, V, distribution = NULL, coefficients = NULL) {
  beta = if (is.null(coefficients)) inference_coefficients(object, V) else coefficients
  se = sqrt(Matrix::diag(V))
  statistic = beta / se
  normal = if (is.null(distribution)) inherits(object, "dbglm") && object[["family"]][["family"]] %in% c("binomial", "poisson") else distribution == "normal"
  p = if (normal) 2 * stats::pnorm(-abs(statistic)) else 2 * stats::pt(-abs(statistic), attr(V, "df.t"))
  table = cbind(estimate = beta, std.error = se, statistic = statistic, p.value = p)
  attr(table, "dbreg_coefficient_ids") = attr(V, "dbreg_coefficient_ids")
  table
}

#' Summarize a database regression
#'
#' Without a covariance request, returns the fitted covariance and coefficient
#' table. Post-fit requests use the immutable sample saved by `retain="all"`.
#' Supported requests are IID, HC1, one-/two-way OLS and one-way GLM clustering.
#' @param object A `dbreg` or `dbglm` model.
#' @param vcov,cluster Optional covariance type or cluster formula/column names.
#' @param ssc Optional list with `K.adj`, `K.fixef` (`"full"`, `"none"`,
#'   `"nonnested"`), `G.adj`, `G.df` (`"min"`, `"conventional"`) and `t.df`
#'   (`"min"`, `"conventional"`, or a positive numeric df). `K.exact=TRUE` is
#'   unsupported. Defaults inherit the fitted native policy.
#' @param vcov_fix Repair nonpositive covariance eigenvalues when needed.
#' @param ... Additional arguments are not supported.
#' @return A fitted-model summary with the selected covariance and coefficient
#'   table. Point estimates and the fitted sample are unchanged.
#' @export
summary.dbreg = function(object, vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE, ...) {
  check_unused_args(list(...))
  changed = !is.null(vcov) || !is.null(cluster) || !is.null(ssc) || !identical(vcov_fix, TRUE)
  if (changed) {
    inference = requested_inference_result(object, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
    compact = function(value) if (inherits(value, "formula")) stats::as.formula(paste(deparse(value, width.cutoff = 500L), collapse = ""), env = baseenv()) else value
    if (!is.null(vcov) || !is.null(cluster)) object[["covariance_request"]] = list(vcov = compact(vcov), cluster = compact(cluster))
    object[["all_coefficients"]] = coef.dbreg(object, fe = TRUE)
    object = install_inference_result(object, inference, inference_coefficients(object, inference[["covariance"]]))
    object[["postfit_covariance"]] = TRUE
  }
  class(object) = unique(c(if (inherits(object, "dbglm")) "summary.dbglm", "summary.dbreg", class(object)))
  object
}

#' @rdname summary.dbreg
#' @export
summary.dbglm = summary.dbreg

#' @rdname summary.dbreg
#' @param x A database regression summary.
#' @export
print.summary.dbreg = function(x, ...) {
  if (inherits(x, "dbglm")) print.dbglm(x, ...) else print.dbreg(x, ...)
}

#' @rdname summary.dbreg
#' @export
print.summary.dbglm = print.summary.dbreg

inference_confint = function(object, V, parm, level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("level must lie between zero and one.", call. = FALSE)
  beta = inference_coefficients(object, V)
  indices = if (is.null(parm)) seq_along(beta) else if (is.numeric(parm)) seq_along(beta)[parm] else match(parm, names(beta))
  if (anyNA(indices)) stop("Unknown coefficient in parm.", call. = FALSE)
  normal = inherits(object, "dbglm") && object[["family"]][["family"]] %in% c("binomial", "poisson")
  probabilities = c((1 - level) / 2, (1 + level) / 2)
  critical = if (normal) stats::qnorm(probabilities) else stats::qt(probabilities, attr(V, "df.t"))
  answer = outer(sqrt(Matrix::diag(V))[indices], critical) + beta[indices]
  rownames(answer) = names(beta)[indices]
  colnames(answer) = paste0(format(100 * probabilities), " %")
  if (!is.null(attr(V, "dbreg_coefficient_ids"))) attr(answer, "dbreg_coefficient_ids") = attr(V, "dbreg_coefficient_ids")[indices]
  answer
}

materialize_retained_rows = function(object, query, value_column = NULL) {
  state = retained_inference(object)
  conn = fitted_resource_connection(object[["resource_state"]])
  resources = fit_resources()
  resources[["conn"]] = conn
  resources[["retained_session"]] = object[["resource_state"]][["session"]]
  on.exit(release_fit_resources(resources))
  table = basename(tempfile("dbreg_retained_rows_"))
  resources[["tables"]] = table
  DBI::dbExecute(conn, sql_sample_table(conn, table, query, c(state[["row_id"]], value_column)))
  lease = retain_fit_resources(resources, table)
  structure(list(
    conn = fitted_resource_public_connection(lease), table = table, row_id = state[["row_id"]],
    value_column = value_column, resource_state = lease
  ), class = "dbreg_rows")
}

#' Retained database regression rows
#'
#' `case.names()` returns original R positions for R input, or a keyed database relation
#' for SQL input. `case.names()` uses the same observation identities. These
#' operations require fitting with `retain="all"` and its live owning session.
#' @param x,object A fitted database regression.
#' @param ... Additional arguments are not supported.
#' @return Original R positions, or a `dbreg_rows` relation. Database collection
#'   is explicit with `as.data.frame()`; `head()` previews and `close()` releases
#'   the result. Caller connections remain open.
#' @importFrom stats case.names fitted residuals
#' @name retained-observations
retained_observations = function(x, ...) {
  check_unused_args(list(...))
  state = retained_inference(x)
  conn = fitted_resource_connection(x[["resource_state"]])
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  from = paste("FROM", q(state[["table"]]), "WHERE", q(state[["estimated"]]))
  if (isTRUE(state[["r_input"]])) {
    positions = DBI::dbGetQuery(conn, paste("SELECT", q(state[["key"]]), "AS position", from, "ORDER BY", q(state[["key"]])))[["position"]]
    return(if (state[["nobs_input"]] <= .Machine[["integer.max"]]) as.integer(positions) else as.numeric(positions))
  }
  materialize_retained_rows(x, paste("SELECT", paste(q(state[["row_id"]]), collapse = ","), from))
}

#' @rdname retained-observations
#' @export
case.names.dbreg = function(object, ...) retained_observations(object, ...)

#' @rdname retained-observations
#' @export
as.data.frame.dbreg_rows = function(x, ...) as.data.frame.dbreg_prediction(x, ...)

#' @rdname retained-observations
#' @param n Number of rows to preview.
#' @export
head.dbreg_rows = function(x, n = 6L, ...) head.dbreg_prediction(x, n, ...)

#' @rdname retained-observations
#' @param con A `dbreg_rows` result.
#' @export
close.dbreg_rows = function(con, ...) close.dbreg_prediction(con, ...)

#' @rdname retained-observations
#' @export
print.dbreg_rows = function(x, ...) {
  check_unused_args(list(...))
  cat("Database rows in ", x[["table"]], "; keys: ", paste(x[["row_id"]], collapse = ", "), ".\n", sep = "")
  invisible(x)
}

#' Fitted values and residuals on the saved sample
#' @param object A fitted `dbreg` or `dbglm` model with `retain="all"`.
#' @param type Residual type: response, deviance, Pearson, or GLM working.
#' @param na.rm Return only estimation rows (`TRUE`), or include original
#'   excluded positions as missing (`FALSE`).
#' @param ... Additional arguments are not supported.
#' @return A numeric vector for R inputs, or a keyed database relation for SQL
#'   inputs. These methods require the live retained sample.
#' @name retained-values
NULL

#' @rdname retained-values
#' @export
fitted.dbreg = function(object, ...) {
  check_unused_args(list(...))
  predict_retained_sample(object, type = "response")
}

#' @rdname retained-values
#' @export
fitted.dbglm = fitted.dbreg

# Stable Poisson deviance around y=mu. The odd-power series avoids subtracting
# nearly equal O(y) values when the desired deviance is O((y-mu)^2/y).
sql_poisson_bd0 = function(y, mu) {
  v = paste0("((", y, "-", mu, ")/(", y, "+", mu, "))")
  polynomial = "0.0"
  for (j in 12:1) polynomial = paste0("(", dbglm_number(1 / (2 * j + 1)), "+POWER(", v, ",2)*", polynomial, ")")
  series = paste0("(", y, "-", mu, ")*", v, "+2*", y, "*POWER(", v, ",3)*", polynomial)
  paste0(
    "CASE WHEN ", y, "=0 THEN ", mu, " WHEN ABS(", y, "-", mu, ")<0.1*(", y, "+", mu,
    ") THEN (", series, ") ELSE ", y, "*(LN(", y, ")-LN(", mu, "))-(", y, "-", mu, ") END"
  )
}

#' @rdname retained-values
#' @export
residuals.dbreg = function(object, type = c("response", "deviance", "pearson", "working"), na.rm = TRUE, ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) stop("na.rm must be TRUE or FALSE.", call. = FALSE)
  state = retained_inference(object)
  conn = fitted_resource_connection(object[["resource_state"]])
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  e = q(state[["residual"]])
  eta = q(state[["eta"]])
  w = if (is.null(state[["weights"]])) "1.0" else paste0("CAST(", q(state[["weights"]]), " AS DOUBLE)")
  glm = inherits(object, "dbglm")
  if (!glm && type == "working") stop("Working residuals require a GLM.", call. = FALSE)
  value = e
  if (!glm && type %in% c("deviance", "pearson")) value = paste0(e, "*SQRT(", w, ")")
  if (glm && type != "response") {
    mu = q(state[["mu"]])
    gaussian = object[["family"]][["family"]] == "gaussian"
    poisson = object[["family"]][["family"]] == "poisson"
    if (type == "working") {
      derivative = if (gaussian) "1.0" else if (poisson) paste0("GREATEST(EXP(", eta, "),", dbglm_number(.Machine[["double.eps"]]), ")") else
        dbglm_family_sql(object[["family"]][["link"]], eta)[["derivative"]]
      value = paste0(e, "/(", derivative, ")")
    } else if (type == "pearson") {
      variance = if (gaussian) "1.0" else if (poisson) mu else paste0(mu, "*(1.0-", mu, ")")
      value = paste0(e, "*SQRT(", w, ")/SQRT(", variance, ")")
    } else {
      y = paste0("CAST(", q(object[["yvar"]]), " AS DOUBLE)")
      deviance = if (gaussian) paste0("POWER(", e, ",2)/2.0") else if (poisson) sql_poisson_bd0(y, mu) else paste0("CASE WHEN ", y, "=0 THEN -LN(1.0-", mu, ") ELSE -LN(", mu, ") END")
      value = paste0("(CASE WHEN ", e, "<0 THEN -1.0 ELSE 1.0 END)*SQRT(GREATEST(2.0*", w, "*(", deviance, "),0.0))")
    }
  }
  value = paste0("CASE WHEN ", q(state[["estimated"]]), " THEN (", value, ") ELSE CAST(NULL AS DOUBLE) END")
  from = paste("FROM", q(state[["table"]]), if (na.rm) paste("WHERE", q(state[["estimated"]])) else "")
  if (isTRUE(state[["r_input"]])) return(DBI::dbGetQuery(conn, paste("SELECT", value, "AS residual", from, "ORDER BY", q(state[["key"]])))[["residual"]])
  column = safe_identifiers("residual", taken = state[["row_id"]], fold_case = folds_identifier_case(conn))
  query = paste("SELECT", paste(c(q(state[["row_id"]]), paste(value, "AS", q(column))), collapse = ","), from)
  materialize_retained_rows(object, query, column)
}

#' @rdname retained-values
#' @export
residuals.dbglm = residuals.dbreg

#' Residual degrees of freedom for a database regression
#' @param object A fitted `dbreg` model.
#' @param ... Additional arguments are not supported.
#' @return Model residual degrees of freedom.
#' @export
df.residual.dbreg = function(object, ...) {
  check_unused_args(list(...))
  object[["df_residual"]]
}

#' Gaussian likelihood for database OLS
#' Uses the unweighted residual sum of squares, including for weighted fits,
#' matching `fixest::feols()`. Weighted models require `retain="all"` to save
#' that statistic. Cached mathematical likelihood metadata survives release.
#' @param object A fitted OLS `dbreg` model.
#' @param ... Additional arguments are not supported.
#' @return A `logLik` object with observation and parameter counts.
#' @export
logLik.dbreg = function(object, ...) {
  check_unused_args(list(...))
  if (inherits(object, "dbivreg")) stop("OLS likelihood is not supported for IV models.", call. = FALSE)
  state = object[["fitted_state"]][["inference"]]
  if (!is.null(object[["weights"]]) && is.null(state)) stop("Weighted OLS likelihood requires retain='all' or retain='inference'.", call. = FALSE)
  n = object[["nobs_orig"]]
  rss = if (!is.null(state)) state[["rss_unweighted"]] else attr(object[["vcov"]], "rss")
  if (is.null(rss) || !is.finite(rss) || rss < 0) stop("Accurate OLS likelihood requires retain='all' or retain='inference'.", call. = FALSE)
  K = if (!is.null(state)) state[["K_model"]] else n - object[["df_residual"]]
  value = -n / 2 - n / 2 * log(rss / n) - n / 2 * log(2 * pi)
  structure(value, class = "logLik", df = K + 1, nobs = n)
}

retain_linear_score_frame = function(fit, frame, y, inputs) {
  if (!retain_sample(inputs[["retain"]])) return(invisible(NULL))
  conn = inputs[["conn"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  residual = fit_residual_sql(fit, frame, y, conn)
  x = vapply(fit[["x"]], function(name) if (name == "(Intercept)") "1.0" else
    fit_column_sql(fit, frame, name, conn), character(1))
  w = frame[["weights_expr"]]
  if (is.null(w)) w = "1.0"
  scale = fit[["weight_scale"]]
  if (is.null(scale)) scale = 1
  w = paste0("(", w, ")/(", dbglm_number(scale), ")")
  scores = sprintf("CAST((%s) AS DOUBLE)*(%s)*(%s)", w, residual, x)
  score_names = if (length(x)) paste0("dbreg_fit_score", seq_along(x)) else character()
  select = c(
    paste(inputs[["retained_key_safe"]], "AS dbreg_fit_key__"),
    paste(residual, "AS dbreg_fit_residual__"), if (length(scores)) paste(scores, "AS", score_names)
  )
  query = paste(frame[["cte_sql"]], "SELECT", paste(select, collapse = ","), "FROM", frame[["cte_name"]])
  original = inputs[["retained_source"]]
  public_y = restore_names(y, inputs[["source_map"]])
  weights = restore_names(inputs[["weights"]], inputs[["source_map"]])
  components = lapply(inputs[["prediction_fe_components"]], function(columns) restore_names(unique(columns), inputs[["source_map"]]))
  state = retain_score_sample(conn, original, inputs[["sample_resources"]], query, "dbreg_fit_key__",
    if (length(scores)) paste0("s.", score_names) else character(),
    paste0("CAST(r.", q(public_y), " AS DOUBLE)-s.dbreg_fit_residual__"), "s.dbreg_fit_residual__",
    row_id = inputs[["row_id"]], fe_components = components, weights = weights
  )
  if (!is.finite(state[["rss"]]) || abs(state[["rss"]] - fit[["rss"]]) > 1e-6 * max(1, abs(fit[["rss"]]))) {
    stop("Retained OLS scores do not reproduce the fitted RSS.", call. = FALSE)
  }
  state[["bread"]] = unname(fit[["bread"]]) * scale
  if (!is.null(fit[["sandwich_bread"]])) state[["score_projection"]] = t(unname(fit[["sandwich_bread"]]) * scale)
  state[["coefficient_names"]] = fit[["x"]]
  state[["kind"]] = "ols"
  state[["weight_scale"]] = scale
  inputs[["inference_state"]] = state
  invisible(NULL)
}

finalize_linear_inference = function(result, inputs) {
  if (!retain_sample(inputs[["retain"]])) return(result)
  state = inputs[["inference_state"]]
  if (is.null(state)) stop("This strategy has not retained its inference frame.", call. = FALSE)
  if (state[["nobs"]] != result[["nobs_orig"]]) stop("Retained score rows differ from the fitted sample.", call. = FALSE)
  layout = result[["linear_design"]]
  if (!is.null(layout) && is.null(state[["coefficient_ids"]])) {
    positions = match(state[["coefficient_names"]], layout[["columns"]][["alias"]])
    if (anyNA(positions)) stop("Retained scores do not match fitted coordinate IDs.", call. = FALSE)
    state[["coefficient_ids"]] = layout[["columns"]][["id"]][positions]
  }
  map = coefficient_name_map(inputs[["name_map"]], inputs[["source_map"]])
  names = standardize_coef_names(restore_names(state[["coefficient_names"]], map))
  public_names = if (!is.null(state[["coefficient_ids"]])) {
    names(linear_coefficients_for_ids(result, state[["coefficient_ids"]]))
  } else names(coef.dbreg(result))
  if (is.null(public_names)) public_names = character()
  if (!identical(names, public_names)) stop("Retained inference coefficients do not match the fitted design.", call. = FALSE)
  state[["coefficient_names"]] = names
  dimnames(state[["bread"]]) = list(names, names)
  state[["fe_df"]] = fe_intercept_df(state[["fe_sizes"]])
  state[["K_model"]] = length(names) + linear_fe_df(inputs, state[["fe_sizes"]])
  state[["default_ssc"]] = inputs[["ssc"]]
  state[["default_vcov"]] = if (inputs[["vcov_type_req"]] == "cluster")
    stats::reformulate(restore_names(inputs[["cluster_var"]], inputs[["source_map"]]), env = baseenv()) else inputs[["vcov_type_req"]]
  result = attach_inference_state(result, state)
  result
}

attach_inference_state = function(result, inference) {
  state = result[["fitted_state"]]
  if (is.null(state)) state = list(
    version = 1L, specification = NULL, coefficients = NULL,
    sample = list(nobs = inference[["nobs"]]), inference = NULL, resources = NULL, prediction = NULL
  )
  state[["inference"]] = inference
  state[["resources"]] = list(
    storage = "database", version = 1L,
    tables = unique(c(state[["resources"]][["tables"]], inference[["resources"]][["tables"]])),
    views = unique(c(state[["resources"]][["views"]], inference[["resources"]][["views"]]))
  )
  result[["fitted_state"]] = state
  result
}

retain_glm_score_frame = function(
  result, fit, conn, base, x, fe, original, resources,
  components, row_id, weights, weight_scale, counts, ssc
) {
  if (isFALSE(result[["finite_mle"]])) stop("Full inference retention requires a finite GLM fit.", call. = FALSE)
  prefix = ""
  if (!length(fe)) {
    prefix = dbglm_working_sql(base, x, fit[["beta"]], result[["family"]][["link"]], conn = conn)
    from = "FROM working b"
    projected = paste0("b.", x)
  } else if (!is.null(fit[["raw_projected"]])) {
    from = paste("FROM", fit[["raw_projected"]], "b")
    projected = paste0("b.", x, "_tilde")
  } else {
    from = paste("FROM", fit[["working"]], "b", if (length(x)) paste("JOIN", fit[["projected"]], "p USING (__dbglm_key)") else "")
    projected = paste0("p.", x)
  }
  score_names = if (length(x)) paste0("dbreg_fit_score", seq_along(x)) else character()
  terms = c(
    "b.__dbglm_key AS dbreg_fit_key__", "b.eta AS dbreg_fit_eta__", "b.mu AS dbreg_fit_mu__",
    "(b.__dbglm_y-b.mu) AS dbreg_fit_residual__", if (length(x)) paste0(projected, "*b.score AS ", score_names)
  )
  query = paste(prefix, "SELECT", paste(terms, collapse = ","), from)
  state = retain_score_sample(conn, original, resources, query, "dbreg_fit_key__",
    if (length(x)) paste0("s.", score_names) else character(), "s.dbreg_fit_eta__", "s.dbreg_fit_residual__",
    mu = "s.dbreg_fit_mu__", row_id = row_id, fe_components = components, weights = weights
  )
  if (state[["nobs"]] != result[["nobs_orig"]]) stop("Retained GLM scores differ from the fitted sample.", call. = FALSE)
  state[["bread"]] = if (length(x)) dbglm_solve(fit[["H"]])[["inverse"]] else matrix(numeric(), 0L, 0L)
  state[["coefficient_names"]] = names(stats::coef(structure(result, class = c("dbglm", "dbreg"))))
  dimnames(state[["bread"]]) = rep(list(state[["coefficient_names"]]), 2L)
  state[["kind"]] = if (result[["family"]][["family"]] == "gaussian") "gaussian" else "glm"
  state[["dispersion_scaled"]] = result[["dispersion_scaled"]]
  state[["weight_scale"]] = weight_scale
  state[["K_model"]] = counts[["model"]]
  state[["fe_df"]] = fe_intercept_df(state[["fe_sizes"]])
  state[["default_ssc"]] = ssc
  state[["default_vcov"]] = if (result[["vcov_type"]] == "cluster") stats::reformulate(result[["cluster_var"]], env = baseenv()) else result[["vcov_type"]]
  attach_inference_state(result, state)
}

# Restricted, typed sample predicates shared by the native and adapter APIs.
subset_variables = function(subset) {
  if (is.null(subset)) return(character())
  if (is.character(subset) && length(subset) == 1L && !is.na(subset) && nzchar(subset)) return(subset)
  if (!inherits(subset, "formula") || length(subset) != 2L) stop("`subset` must name a logical column or be a one-sided predicate formula.", call. = FALSE)
  all.vars(subset)
}

sql_subset_predicate = function(conn, subset, schema, to_safe = identity) {
  if (is.null(subset)) return(NULL)
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Subset predicates support DuckDB and native Databricks.", call. = FALSE)
  }
  subset_variables(subset)
  expr = if (is.character(subset)) as.name(subset) else subset[[2L]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  missing = function(value) paste0(
    "(", value[["sql"]], " IS NULL",
    if (value[["type"]] == "numeric") paste0(" OR isnan(CAST(", value[["sql"]], " AS DOUBLE))"), ")"
  )
  node = function(sql, type) list(sql = sql, type = type)
  literal = function(value) {
    if (length(value) != 1L || !(is.numeric(value) || is.character(value) || is.logical(value))) stop("Unsupported subset literal.", call. = FALSE)
    node(sql_value_literals(conn, value), if (is.character(value)) "character" else if (is.logical(value)) "logical" else "numeric")
  }
  constant = function(expr) {
    if (is.call(expr)) {
      op = as.character(expr[[1L]])
      args = as.list(expr)[-1L]
      if (op == "c") return(unlist(lapply(args, constant)))
      if (op %in% c("+", "-") && length(args) == 1L) {
        value = constant(args[[1L]])
        if (length(value) == 1L && is.numeric(value)) return(if (op == "-") -value else value)
      }
      stop("Subset %in% requires literal values.", call. = FALSE)
    }
    if (is.symbol(expr)) stop("Subset %in% requires literal values.", call. = FALSE)
    expr
  }
  visit = function(expr) {
    if (is.symbol(expr)) {
      column = to_safe(as.character(expr))
      if (!column %in% names(schema)) stop("Unknown subset column: ", as.character(expr), call. = FALSE)
      value = schema[[column]]
      type = if (is.logical(value)) "logical" else if (is.numeric(value) || inherits(value, "integer64")) "numeric" else if (is.character(value) || is.factor(value)) "character" else "other"
      return(node(q(column), type))
    }
    if (!is.call(expr)) return(literal(expr))
    op = as.character(expr[[1L]])
    args = as.list(expr)[-1L]
    if (op == "(" && length(args) == 1L) return(visit(args[[1L]]))
    if (op == "is.na" && length(args) == 1L) return(node(missing(visit(args[[1L]])), "logical"))
    if (op == "!" && length(args) == 1L) {
      value = visit(args[[1L]])
      if (value[["type"]] != "logical") stop("Subset ! requires a logical operand.", call. = FALSE)
      return(node(paste0("(NOT ", value[["sql"]], ")"), "logical"))
    }
    if (op == "%in%" && length(args) == 2L) {
      value = visit(args[[1L]])
      constants = constant(args[[2L]])
      if (!length(constants)) return(node("FALSE", "logical"))
      present = constants[!is.na(constants)]
      type = if (is.character(constants)) "character" else if (is.logical(constants)) "logical" else "numeric"
      if (length(present) && value[["type"]] != type) stop("Subset %in% requires matching column and literal types.", call. = FALSE)
      membership = if (length(present)) paste0("COALESCE(", value[["sql"]], " IN (", paste(sql_value_literals(conn, present), collapse = ","), "),FALSE)") else "FALSE"
      has_nan = is.numeric(constants) && any(is.nan(constants))
      has_na = any(is.na(constants) & !(if (is.numeric(constants)) is.nan(constants) else FALSE))
      return(node(paste0(
        "(CASE WHEN ", value[["sql"]], " IS NULL THEN ", if (has_na) "TRUE" else "FALSE",
        if (value[["type"]] == "numeric") paste0(" WHEN isnan(CAST(", value[["sql"]], " AS DOUBLE)) THEN ", if (has_nan) "TRUE" else "FALSE"),
        " ELSE ", membership, " END)"
      ), "logical"))
    }
    if (op %in% c("+", "-") && length(args) == 1L) {
      value = visit(args[[1L]])
      if (value[["type"]] != "numeric") stop("Subset arithmetic requires numeric operands.", call. = FALSE)
      return(node(paste0("(", op, value[["sql"]], ")"), "numeric"))
    }
    if (length(args) == 2L && op %in% c("+", "-", "*", "/", "^", "==", "!=", ">", "<", ">=", "<=", "&", "|")) {
      left = visit(args[[1L]])
      right = visit(args[[2L]])
      if (op %in% c("&", "|")) {
        if (left[["type"]] != "logical" || right[["type"]] != "logical") stop("Subset boolean operators require logical operands.", call. = FALSE)
        return(node(paste0("(", left[["sql"]], if (op == "&") " AND " else " OR ", right[["sql"]], ")"), "logical"))
      }
      if (op %in% c("+", "-", "*", "/", "^")) {
        if (left[["type"]] != "numeric" || right[["type"]] != "numeric") stop("Subset arithmetic requires numeric operands.", call. = FALSE)
        return(node(if (op == "^") paste0("POWER(", left[["sql"]], ",", right[["sql"]], ")") else paste0("(", left[["sql"]], op, right[["sql"]], ")"), "numeric"))
      }
      if (left[["type"]] != right[["type"]] && !left[["sql"]] %in% c("NULL", "NA") && !right[["sql"]] %in% c("NULL", "NA")) stop("Subset comparisons require matching operand types.", call. = FALSE)
      comparison = paste0("(", left[["sql"]], if (op == "==") "=" else op, right[["sql"]], ")")
      return(node(paste0(
        "(CASE WHEN ", missing(left), " OR ", missing(right),
        " THEN CAST(NULL AS BOOLEAN) ELSE ", comparison, " END)"
      ), "logical"))
    }
    stop("Unsupported subset expression: ", paste(deparse(expr), collapse = ""), call. = FALSE)
  }
  result = visit(expr)
  if (result[["type"]] != "logical") stop("Subset formulas must return a logical predicate.", call. = FALSE)
  result[["sql"]]
}

filter_source_predicate = function(from_statement, predicate) {
  if (is.null(predicate)) return(from_statement)
  paste0("FROM (SELECT * ", from_statement, ") dbreg_subset_source WHERE (", predicate, ")")
}

linear_prune_source = function(inputs, from_statement, fe, y) {
  conn = inputs[["conn"]]
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("OLS singleton pruning supports DuckDB and native Databricks.", call. = FALSE)
  }
  resources = inputs[["sample_resources"]]
  resources[["conn"]] = conn
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, backend[["name"]]))
  flags = safe_identifiers(c("__dbglm_single", "__dbglm_infinite"), taken = names(schema), fold_case = folds_identifier_case(conn))
  manager = fit_table_manager(conn, resources, backend)
  source = paste0("(SELECT * ", from_statement, ")")
  pruning = dbglm_prune(conn, source, fe, y, "singletons", manager[["create"]], manager[["drop"]], flags = flags)
  if (pruning[["rounds"]] > 0 && as.numeric(DBI::dbGetQuery(
    conn,
    paste("SELECT COUNT(*) AS n FROM", pruning[["table"]])
  )[["n"]]) == 0) {
    stop("No observations remain after singleton removal.", call. = FALSE)
  }
  pruning
}

# Source keys are typed tuples. Never infer their alignment from fetch order.
validate_row_id_arg = function(row_id) {
  if (!is.null(row_id) && (!is.character(row_id) || !length(row_id) ||
    anyNA(row_id) || any(!nzchar(row_id)) || anyDuplicated(row_id))) {
    stop("`row_id` must name unique source key columns.", call. = FALSE)
  }
  invisible(NULL)
}

validate_source_row_id = function(conn, from_statement, row_id) {
  if (is.null(row_id)) return(invisible(NULL))
  quoted = as.character(DBI::dbQuoteIdentifier(conn, row_id))
  counter = safe_identifiers("dbreg_key_count__",
    taken = row_id,
    fold_case = folds_identifier_case(conn)
  )
  count = sql_count_expr(detect_backend(conn)[["name"]])
  missing = paste0(quoted, " IS NULL")
  backend = detect_backend(conn)
  if (backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) {
    text = sql_text_type(conn, backend)
    missing = c(missing, paste0(
      "(LOWER(typeof(", quoted, ")) IN ('double','float','real') AND LOWER(CAST(",
      quoted, " AS ", text, "))='nan')"
    ))
  }
  # Compare counts in SQL so even counts outside double's exact integer range
  # cannot turn a duplicate into a false successful validation in R.
  info = DBI::dbGetQuery(conn, paste0(
    "SELECT CASE WHEN COALESCE(SUM(", counter,
    "),0)=", count, " THEN 0 ELSE 1 END AS duplicate_keys,COALESCE(SUM(CASE WHEN ",
    paste(missing, collapse = " OR "), " THEN ", counter,
    " ELSE 0 END),0) AS n_null FROM (SELECT ", paste(quoted, collapse = ","),
    ",", sql_count(conn, counter), " ", from_statement,
    " GROUP BY ", paste(quoted, collapse = ","), ") dbreg_keys"
  ))
  if (info[["n_null"]] > 0 || info[["duplicate_keys"]] > 0) {
    stop("Source row_id must be unique and nonmissing.", call. = FALSE)
  }
  invisible(NULL)
}

# Both supported export backends accept this CREATE form. Deliberately avoid
# OR REPLACE: an export must never overwrite a caller-owned relation.
sample_table_exists = function(conn, table) {
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    # ODBC's SQLTables metadata can omit session tables. Do not use SHOW LIKE:
    # table names may contain literal pattern metacharacters such as * or |.
    tables = DBI::dbGetQuery(conn, "SHOW TABLES")[["tableName"]]
    return(any(tolower(tables) == tolower(table)))
  }
  DBI::dbExistsTable(conn, table)
}

sql_sample_table = function(conn, table, select_sql, columns, replace = FALSE) {
  # Delta's default Parquet naming constraints also apply to session tables.
  # Enable mapping only on our new tables when preserving such source names.
  # https://docs.databricks.com/aws/en/tables/features/column-mapping
  mapping = if (isTRUE(detect_backend(conn)[["is_databricks"]]) &&
    any(grepl("[ ,;{}()\n\r\t=]", columns))) {
    "TBLPROPERTIES ('delta.columnMapping.mode'='name')"
  } else ""
  paste(
    "CREATE", if (replace) "OR REPLACE" else "", "TEMP TABLE",
    DBI::dbQuoteIdentifier(conn, table), mapping, "AS", select_sql
  )
}

export_retained_sample = function(conn, table, keys, from_statement, public_names = keys) {
  columns = paste(
    DBI::dbQuoteIdentifier(conn, keys), "AS",
    DBI::dbQuoteIdentifier(conn, public_names)
  )
  sql = paste("SELECT", paste(columns, collapse = ","), from_statement)
  DBI::dbExecute(conn, sql_sample_table(conn, table, sql, public_names))
  invisible(NULL)
}

# Explicit sample exports use an owned snapshot for the whole fit. Default
# compact fits keep their existing query and resource behavior.
prepare_linear_sample = function(inputs, from_statement, variables) {
  keys = inputs[["row_id"]]
  table = inputs[["sample_table"]]
  full = retain_sample(inputs[["retain"]])
  prediction = retain_prediction(inputs[["retain"]]) && length(inputs[["fe_labels"]]) > 0L
  if (is.null(keys) && is.null(table) && !prediction && !full) return(from_statement)
  conn = inputs[["conn"]]
  if (full && !is.null(table) && sample_table_exists(conn, table)) stop("Sample export table already exists.", call. = FALSE)
  if (full) {
    original = materialize_retained_source(
      conn, from_statement, inputs[["sample_resources"]],
      inputs[["retained_input"]][["key"]], inputs[["retained_input"]][["r_input"]]
    )
    inputs[["retained_source"]] = original
    from_statement = paste("FROM", DBI::dbQuoteIdentifier(conn, original[["table"]]))
  }
  if ((!is.null(table) || prediction) && !full) {
    backend = detect_backend(conn)
    if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
      stop("Sample snapshots currently support DuckDB and native Databricks.", call. = FALSE)
    }
    if (!is.null(table) && sample_table_exists(conn, table)) stop("Sample export table already exists.", call. = FALSE)
    resources = inputs[["sample_resources"]]
    resources[["conn"]] = conn
    columns = unique(c(variables, keys))
    snapshot = basename(tempfile("dbreg_sample_"))
    # Track before CREATE so partially failed setup is cleaned by on.exit too.
    resources[["tables"]] = c(resources[["tables"]], snapshot)
    sql = paste("SELECT", paste(DBI::dbQuoteIdentifier(conn, columns), collapse = ","), from_statement)
    DBI::dbExecute(conn, sql_sample_table(conn, snapshot, sql, columns))
    from_statement = paste("FROM", DBI::dbQuoteIdentifier(conn, snapshot))
  }
  validate_source_row_id(conn, from_statement, keys)
  from_statement
}

finalize_linear_sample = function(result, inputs) {
  result[["row_id"]] = inputs[["row_id"]]
  result[["sample_table"]] = inputs[["sample_table"]]
  if (!is.null(inputs[["sample_table"]])) {
    export_retained_sample(
      inputs[["conn"]], inputs[["sample_table"]],
      inputs[["row_id_safe"]], inputs[["from_statement"]], inputs[["row_id"]]
    )
  }
  result
}

# Reuse the fitted cell design and inverse information. A compact lookup gives
# each original row its fitted mean and coefficient-space score multiplier;
# no second WLS or observation-level FE projection is needed.
retain_compressed_score_sample = function(
  fit, score_X, yhat, cells,
  groups, from_statement, layout, inputs
) {
  conn = inputs[["conn"]]
  resources = inputs[["sample_resources"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  selected_ids = if (length(inputs[["fe"]])) layout[["slope_ids"]] else layout[["fitted_ids"]]
  selected = match(selected_ids, layout[["fitted_ids"]])
  names = fit[["x"]][selected]
  scale = fit[["weight_scale"]]
  projection = fit[["sandwich_bread"]][selected, , drop = FALSE] * scale
  multipliers = as.matrix(score_X %*% t(projection))
  private = safe_identifiers(c("dbreg_cell_fit__", if (length(selected)) paste0("dbreg_cell_score", seq_along(selected), "__")),
    taken = groups, fold_case = folds_identifier_case(conn)
  )
  lookup = as.data.frame(lapply(cells[groups], function(value) {
    if (!is.factor(value)) return(value)
    raw = attr(value, "dbreg_key_values")
    if (is.null(raw)) as.character(value) else raw[as.integer(value)]
  }), check.names = FALSE)
  lookup[[private[[1L]]]] = yhat
  for (i in seq_along(selected)) lookup[[private[[i + 1L]]]] = multipliers[, i]
  types = if (length(groups)) {
    row = DBI::dbGetQuery(conn, paste(
      "SELECT", paste0("typeof(", q(groups), ") AS t", seq_along(groups), collapse = ","),
      from_statement, "LIMIT 1"
    ))
    vapply(row, function(value) as.character(value[[1L]]), character(1))
  } else character()
  table = prediction_upload_frame(
    conn, lookup, resources,
    c(types, rep("DOUBLE", length(private)))
  )
  join = if (length(groups)) paste0(
    " JOIN ", table, " c ON ",
    paste0("f.", q(groups), "=c.", q(groups), collapse = " AND ")
  ) else paste(" CROSS JOIN", table, "c")
  y = inputs[["yvar"]]
  residual = paste0(
    "(", sql_cast_double(paste0("f.", q(y)), conn), "-(",
    dbglm_number(fit[["centers"]][[y]]), ")-c.", q(private[[1L]]), ")"
  )
  weight = if (is.null(inputs[["weights"]])) "1.0" else sql_cast_double(paste0("f.", q(inputs[["weights"]])), conn)
  scores = if (length(selected)) paste0(
    "(", weight, ")/(", dbglm_number(scale), ")*", residual,
    "*c.", q(private[-1L])
  ) else character()
  score_names = if (length(selected)) paste0("dbreg_fit_score", seq_along(selected)) else character()
  query = paste0(
    "SELECT f.", q(inputs[["retained_key_safe"]]), " AS dbreg_fit_key__,",
    residual, " AS dbreg_fit_residual__",
    if (length(scores)) paste0(",", paste(scores, "AS", score_names, collapse = ",")),
    " FROM (SELECT * ", from_statement, ") f", join
  )
  original = inputs[["retained_source"]]
  public_y = restore_names(y, inputs[["source_map"]])
  components = lapply(inputs[["prediction_fe_components"]], function(columns) restore_names(unique(columns), inputs[["source_map"]]))
  materialize = inputs[["retain"]] == "all"
  state = retain_score_sample(conn, original, resources, query, "dbreg_fit_key__",
    if (length(scores)) paste0("s.", score_names) else character(),
    paste0("CAST(r.", q(public_y), " AS DOUBLE)-s.dbreg_fit_residual__"), "s.dbreg_fit_residual__",
    row_id = inputs[["row_id"]], fe_components = components,
    weights = restore_names(inputs[["weights"]], inputs[["source_map"]]), materialize = materialize
  )
  if (!is.finite(state[["rss"]]) || abs(state[["rss"]] - fit[["rss"]]) > 1e-6 * max(1, abs(fit[["rss"]]))) {
    stop("Retained compressed scores do not reproduce the fitted RSS.", call. = FALSE)
  }
  if (!materialize) {
    # The filtered source can refer to an owned singleton-pruning relation.
    # Keep every current dependency until the lazy view is explicitly closed.
    state[["resources"]][["tables"]] = resources[["tables"]]
    state[["resources"]][["views"]] = resources[["views"]]
  }
  state[["bread"]] = unname(fit[["bread"]][selected, selected, drop = FALSE]) * scale
  state[["score_projection"]] = diag(length(selected))
  state[["coefficient_names"]] = names
  state[["coefficient_ids"]] = selected_ids
  state[["kind"]] = "ols"
  state[["weight_scale"]] = scale
  state[["storage"]] = if (materialize) "rows" else "cells"
  state[["n_cells"]] = nrow(cells)
  inputs[["inference_state"]] = state
  invisible(NULL)
}
