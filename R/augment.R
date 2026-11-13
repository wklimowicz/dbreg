#' Collect the retained estimation sample
#'
#' Collects the original fitted rows and appends response-scale fitted values
#' and residuals. Fit with `retain = "inference"` or `"all"` and keep the
#' owning database session open.
#' @param x A retained database OLS or GLM fit.
#' @param data Optional original R input, supplied in its original row order.
#'   This option requires automatic R row positions and preserves the original R
#'   column classes and row names.
#' @param ... Unused.
#' @return A data frame ordered by the original row position or native row key.
#'   Databricks BIGINT source columns are collected as exact decimal text. Internal generated
#'   keys are excluded from source columns. The appended names \code{.fitted}
#'   and \code{.resid} are reserved: conflicting source columns raise an error
#'   rather than being overwritten. With automatic R row positions, renamed
#'   original columns can be supplied through \code{data}.
#' @importFrom generics augment
#' @export
augment.dbreg = function(x, data = NULL, ...) {
  if (length(list(...))) stop("Unused arguments to augment().", call. = FALSE)
  require_inference_rows(x)
  conn = retained_fit_connection(x[["resource_state"]])
  state = x[["retained_inference"]]
  if (is.null(state)) stop("The fit has no retained source rows.", call. = FALSE)
  rows = x[["retained_rows"]]
  glm = inherits(x, "dbglm")
  fields = if (glm) state[["cluster_fields"]] else {
    columns = names(state[["source_classes"]])
    stats::setNames(columns, columns)
  }
  public = if (is.null(data)) names(fields) else names(data)
  conflicts = intersect(c(".fitted", ".resid"), public)
  if (length(conflicts)) stop("Source columns conflict with augmentation fields: ",
    paste(conflicts, collapse = ", "), ". Rename them or supply renamed original data.", call. = FALSE)
  source = if (glm) state[["source"]] else state[["original"]]
  source_keys = state[["source_keys"]]
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  keys = paste0("s.", rows[["keys"]])
  fitted = if (glm) "s.__dbglm_saved_mu" else "s.__dbreg_saved_fit"
  residual = if (glm) "s.__dbglm_saved_y-s.__dbglm_saved_mu" else "s.__dbreg_saved_residual"
  aliases = safe_identifiers(c("dbreg_augment_position", "dbreg_augment_fitted", "dbreg_augment_resid"),
    taken = names(fields), fold_case = folds_identifier_case(conn))
  key_out = names(fields)[fields %in% source_keys]
  if (is.null(data)) {
    expressions = paste0("b.", q(fields))
    backend = detect_backend(conn)
    if (isTRUE(backend[["is_databricks"]]) && length(fields)) {
      types = DBI::dbGetQuery(conn, paste("SELECT", paste0("typeof(", q(fields),
        ") AS ", q(names(fields)), collapse = ","), "FROM", source, "LIMIT 1"))
      exact = names(fields)[vapply(types, function(type) length(type) && tolower(type[1L]) == "bigint", logical(1))]
      i = match(exact, names(fields))
      expressions[i] = sql_cast_text(conn, expressions[i])
    }
    select = c(paste(expressions, "AS", q(names(fields))), paste(c(fitted, residual), "AS", q(aliases[2:3])))
    join = paste0(" JOIN ", source, " b ON ", paste0(keys, "=b.", q(source_keys), collapse = " AND "))
  } else {
    if (!isTRUE(rows[["r_input"]]) || !is.data.frame(data) || !isTRUE(nrow(data) == rows[["input_n"]])) {
      stop("Supply the original R input in its original order, or omit data to collect the retained sample.",
        call. = FALSE)
    }
    select = paste(c(keys, fitted, residual), "AS", q(aliases))
    join = ""
  }
  sql = paste0("SELECT ", paste(select, collapse = ","), " FROM ", rows[["table"]],
    " s", join, " WHERE s.__dbreg_estimated ORDER BY ", paste(keys, collapse = ","))
  values = read_keyed_query(conn, sql, key_out, exact_keys = TRUE)
  if (is.null(data)) data = values[names(fields)]
  else data = data[as.integer(values[[aliases[1L]]]), , drop = FALSE]
  data[[".fitted"]] = values[[aliases[2L]]]
  data[[".resid"]] = values[[aliases[3L]]]
  data
}
