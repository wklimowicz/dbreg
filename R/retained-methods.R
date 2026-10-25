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
    DBI::dbExecute(conn, sql_sample_table(conn, name, backend_sql(conn, sql), c(state[["row_id"]], prediction)))
    name
  }
  table = create(paste("SELECT", paste(c(q(state[["row_id"]]), paste(q(state[["eta"]]), "AS", q(prediction))), collapse = ","), from))
  if (inherits(object, "dbglm") && type == "response") table = create(sql_prediction_response(conn, table, state[["row_id"]], prediction, object[["family"]][["link"]]))
  lease = retain_fit_resources(resources, table)
  structure(list(conn = fitted_resource_public_connection(lease), table = table, row_id = state[["row_id"]],
    prediction_column = prediction, resource_state = lease), class = "dbreg_prediction")
}

inference_coeftable = function(object, V) {
  beta = stats::coef(object)[rownames(V)]
  se = sqrt(diag(V))
  statistic = beta/se
  normal = inherits(object, "dbglm") && object[["family"]][["family"]] %in% c("binomial", "poisson")
  p = if (normal) 2*stats::pnorm(-abs(statistic)) else 2*stats::pt(-abs(statistic), attr(V, "df.t"))
  cbind(estimate = beta, std.error = se, statistic = statistic, p.value = p)
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
    V = vcov.dbreg(object, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
    compact = function(value) if (inherits(value, "formula")) stats::as.formula(paste(deparse(value, width.cutoff = 500L), collapse = ""), env = baseenv()) else value
    if (!is.null(vcov) || !is.null(cluster)) object[["covariance_request"]] = list(vcov = compact(vcov), cluster = compact(cluster))
    object[["all_coefficients"]] = coef.dbreg(object, fe = TRUE)
    object[["coeftable"]] = inference_coeftable(object, V)
    object[["vcov"]] = V
    object[["fitted_state"]][["inference"]][["default_ssc"]] = attr(V,"ssc")
    object[["postfit_covariance"]] = TRUE
    object[["vcov_type"]] = attr(V, "vcov_type")
    for (field in c("cluster_var", "n_clusters", "ssc", "n_params", "df_t")) {
      object[[field]] = attr(V, field)
    }
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
  beta = stats::coef(object)[rownames(V)]
  if (is.null(parm)) parm = names(beta)
  if (is.numeric(parm)) parm = names(beta)[parm]
  if (anyNA(parm) || any(!parm %in% names(beta))) stop("Unknown coefficient in parm.", call. = FALSE)
  normal = inherits(object, "dbglm") && object[["family"]][["family"]] %in% c("binomial", "poisson")
  probabilities = c((1-level)/2, (1+level)/2)
  critical = if (normal) stats::qnorm(probabilities) else stats::qt(probabilities, attr(V, "df.t"))
  answer = outer(sqrt(diag(V))[parm], critical) + beta[parm]
  colnames(answer) = paste0(format(100*probabilities), " %")
  answer
}
