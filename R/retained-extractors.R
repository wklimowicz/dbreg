materialize_retained_rows = function(object, query, value_column = NULL) {
  state = retained_inference(object)
  conn = fitted_resource_connection(object[["resource_state"]])
  resources = dbglm_resources()
  resources[["conn"]] = conn
  resources[["retained_session"]] = object[["resource_state"]][["session"]]
  on.exit(dbglm_release_resources(resources))
  table = basename(tempfile("dbreg_retained_rows_"))
  resources[["tables"]] = table
  DBI::dbExecute(conn, sql_sample_table(conn, table, dbglm_sql(conn, query), c(state[["row_id"]], value_column)))
  lease = retain_fit_resources(resources, table)
  structure(list(conn = fitted_resource_public_connection(lease), table = table, row_id = state[["row_id"]],
    value_column = value_column, resource_state = lease), class = "dbreg_rows")
}

#' Retained database regression rows
#'
#' `obs()` returns original R positions for R input, or a keyed database relation
#' for SQL input. `case.names()` uses the same observation identities. These
#' operations require fitting with `retain="all"` and its live owning session.
#' @param x,object A fitted database regression.
#' @param ... Additional arguments are not supported.
#' @return Original R positions, or a `dbreg_rows` relation. Database collection
#'   is explicit with `as.data.frame()`; `head()` previews and `close()` releases
#'   the result. Caller connections remain open.
#' @importFrom stats case.names fitted residuals
#' @export
obs = function(x, ...) {
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

#' @rdname obs
#' @export
case.names.dbreg = function(object, ...) obs(object, ...)

#' @rdname obs
#' @export
as.data.frame.dbreg_rows = function(x, ...) as.data.frame.dbreg_prediction(x, ...)

#' @rdname obs
#' @param n Number of rows to preview.
#' @export
head.dbreg_rows = function(x, n = 6L, ...) head.dbreg_prediction(x, n, ...)

#' @rdname obs
#' @param con A `dbreg_rows` result.
#' @export
close.dbreg_rows = function(con, ...) close.dbreg_prediction(con, ...)

#' @rdname obs
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
  for (j in 12:1) polynomial = paste0("(", dbglm_number(1/(2*j+1)), "+POWER(", v, ",2)*", polynomial, ")")
  series = paste0("(", y, "-", mu, ")*", v, "+2*", y, "*POWER(", v, ",3)*", polynomial)
  paste0("CASE WHEN ", y, "=0 THEN ", mu, " WHEN ABS(", y, "-", mu, ")<0.1*(", y, "+", mu,
    ") THEN (", series, ") ELSE ", y, "*(LN(", y, ")-LN(", mu, "))-(", y, "-", mu, ") END")
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
  if (isTRUE(state[["r_input"]])) return(DBI::dbGetQuery(conn, dbglm_sql(conn,
    paste("SELECT", value, "AS residual", from, "ORDER BY", q(state[["key"]]))))[["residual"]])
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
  if (!is.null(object[["weights"]]) && is.null(state)) stop("Weighted OLS likelihood requires retain='all'.", call. = FALSE)
  n = object[["nobs_orig"]]
  rss = if (!is.null(state)) state[["rss_unweighted"]] else attr(object[["vcov"]], "rss")
  if (is.null(rss) || !is.finite(rss) || rss < 0) stop("Accurate OLS likelihood requires retain='all'.", call. = FALSE)
  K = if (!is.null(state)) state[["K_model"]] else n-object[["df_residual"]]
  value = -n/2 - n/2*log(rss/n) - n/2*log(2*pi)
  structure(value, class = "logLik", df = K+1, nobs = n)
}
