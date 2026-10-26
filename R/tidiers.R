# Keep glance() to one row when covariance metadata has several components.
glance_metadata = function(value, missing = NA_character_) {
  if (!length(value)) return(missing)
  if (is.list(value)) return(paste(paste0(names(value), "=", unlist(value)), collapse = ";"))
  if (length(value) == 1L) value else paste(value, collapse = ";")
}

#' tidy results
#' @export
generics::tidy

#' Tidiers for `dbreg` objects
#' @description Provides \code{broom::tidy} and \code{broom::glance} methods for
#' "dbreg" objects.
#' @importFrom generics tidy
#' @importFrom stats confint
#' @param x a model of class `dbreg` produced by the \code{\link[dbreg]{dbreg}}
#' function.
#' @param conf.int Logical indicating whether to include confidence intervals.
#' Default is `FALSE`.
#' @param conf.level Confidence level for intervals. Default is 0.95.
#' @param fe Should the fixed effects be tidied too? Default is `FALSE`.
#' @param collin Include regressors removed for collinearity with `NA` values?
#' @param vcov,cluster,ssc Optional covariance request, as in [summary.dbreg()].
#'   `vcov` may also be a numeric covariance matrix for the fitted coefficients.
#' @param ... Additional arguments to tidying method. Currently unused except to
#' handle superseded arguments.
#' @examples
#' mod = dbreg(Temp ~ Wind | Month, data = airquality)
#' tidy(mod, conf.int = TRUE)
#' tidy(mod, conf.int = TRUE, fe = TRUE)
#' glance(mod)
#' @name tidiers
#' @rdname tidiers
#' @export
tidy.dbreg = function(x, conf.int = FALSE, conf.level = 0.95, fe = FALSE,
                      collin = FALSE, vcov = NULL, cluster = NULL, ssc = NULL, ...) {
  custom = is.matrix(vcov)
  if (custom) {
    if (!is.null(cluster) || !is.null(ssc) || isTRUE(fe)) stop("A supplied covariance matrix cannot be combined with cluster, ssc, or fe=TRUE.", call. = FALSE)
    names = names(stats::coef(x))
    if (!is.numeric(vcov) || !identical(dim(vcov), c(length(names), length(names))) ||
      any(!is.finite(vcov)) || !isSymmetric(vcov) || any(diag(vcov) < 0)) {
      stop("`vcov` must be a finite symmetric covariance matrix for the fitted coefficients.", call. = FALSE)
    }
    if (is.null(dimnames(vcov))) dimnames(vcov) = list(names, names)
    if (!setequal(rownames(vcov), names) || !setequal(colnames(vcov), names)) stop("Covariance names must match fitted coefficients.", call. = FALSE)
    df = attr(vcov, "df.t")
    vcov = vcov[names, names, drop = FALSE]
    attr(vcov, "df.t") = if (is.null(df)) x[["df_t"]] else df
    x[["coeftable"]] = inference_coeftable(x, vcov)
  } else if (!is.null(vcov) || !is.null(cluster) || !is.null(ssc)) {
    x = summary.dbreg(x, vcov = vcov, cluster = cluster, ssc = ssc)
  }
  # superseded args handled through ...
  dots = list(...)
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning(
        "The `fes` argument has been superseded by `fe` (without the \"s\") and will be deprecated in a future `dbreg` release.\n"
      )
    }
  }
  ct = x[["coeftable"]]
  if (!isTRUE(fe) && !is.null(x[["fe"]])) {
    coef_names = x[["coef_names"]]
    if (!is.null(coef_names)) {
      ct = ct[coef_names, , drop = FALSE]
    } else {
      xvars = x[["xvars"]]
      ct = ct[xvars, , drop = FALSE]
    }
  }
  if (isTRUE(collin)) {
    terms = names(coef(x, fe = fe, collin = TRUE))
    out = data.frame(
      term = terms,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      row.names = NULL
    )
    matched = match(terms, rownames(ct))
    present = !is.na(matched)
    out[present, colnames(ct)] = ct[matched[present], , drop = FALSE]
  } else {
    out = data.frame(term = rownames(ct), ct, row.names = NULL)
  }

  if (isTRUE(conf.int)) {
    ci = if (custom) inference_confint(x, vcov, NULL, conf.level) else confint(x, level = conf.level, fe = fe)
    ci_match = match(out[["term"]], rownames(ci))
    out[["conf.low"]] = ifelse(is.na(ci_match), NA_real_, ci[ci_match, 1])
    out[["conf.high"]] = ifelse(is.na(ci_match), NA_real_, ci[ci_match, 2])
  }

  return(out)
}

#' glance results
#' @export
generics::glance

#' @importFrom generics glance
#' @rdname tidiers
#' @export
glance.dbreg = function(x, ...) {
  gof_vals = gof(x)
  fe_levels = x[["n_fe_levels"]]
  if (is.null(fe_levels)) {
    fe_levels_text = NA_character_
  } else {
    fe_names = x[["fe_labels"]]
    if (is.null(fe_names)) fe_names = x[["fe"]]
    if (is.null(names(fe_levels)) && length(fe_names) == length(fe_levels)) {
      names(fe_levels) = fe_names
    }
    fe_levels_text = paste(paste0(names(fe_levels), "=", fe_levels), collapse = ";")
  }
  data.frame(
    r.squared = gof_vals["r2"],
    adj.r.squared = gof_vals["adj_r2"],
    rmse = gof_vals["rmse"],
    nobs = x[["nobs_orig"]],
    df.residual = x[["df_residual"]],
    nobs.input = x[["nobs_input"]],
    nobs.dropped = x[["nobs_dropped"]],
    vcov.type = x[["vcov_type"]],
    cluster = glance_metadata(x[["cluster_var"]]),
    n.clusters = glance_metadata(x[["n_clusters"]], NA_integer_),
    ssc = glance_metadata(x[["ssc"]]),
    n.params = if (is.null(x[["n_params"]])) NA_integer_ else x[["n_params"]],
    n.fe = length(x[["fe"]]),
    fe.levels = fe_levels_text,
    strategy = x[["strategy"]],
    ap.converged = x[["ap_converged"]],
    ap.iterations = if (is.null(x[["ap_iterations"]])) NA_integer_ else x[["ap_iterations"]],
    ap.tolerance = x[["ap_tolerance"]],
    ap.max.abs = if (is.null(x[["ap_max_abs"]])) NA_real_ else x[["ap_max_abs"]],
    row.names = NULL
  )
}

select_coefficient_names = function(names, keep = NULL, drop = NULL, order = NULL) {
  for (value in list(keep, drop, order)) if (!is.null(value) && (!is.character(value) || anyNA(value))) stop("Coefficient selectors must be nonmissing regular expressions.", call. = FALSE)
  match = function(pattern, names) {
    pattern = sub("^!%", "%!", pattern)
    pattern = sub("^%", "", pattern)
    negative = startsWith(pattern, "!")
    if (negative) pattern = substring(pattern, 2L)
    found = grepl(pattern, names, perl = TRUE)
    if (negative) !found else found
  }
  if (length(keep)) names = names[Reduce(`|`, lapply(keep, match, names = names))]
  if (length(drop)) names = names[!Reduce(`|`, lapply(drop, match, names = names))]
  for (pattern in rev(order)) {
    first = match(pattern, names)
    names = c(names[first], names[!first])
  }
  names
}

#' Coefficient summaries through fixest generics
#' @param object A native database regression.
#' @param keep,drop,order Optional regular-expression coefficient selectors.
#' @param ... Covariance controls passed to [summary.dbreg()].
#' @return A coefficient matrix or named vector.
#' @name coefficient-summaries
NULL

#' @rdname coefficient-summaries
#' @exportS3Method fixest::coeftable
coeftable.dbreg = function(object, keep = NULL, drop = NULL, order = NULL, ...) {
  value = summary.dbreg(object, ...)[["coeftable"]]
  names = select_coefficient_names(names(stats::coef(object)), keep, drop, order)
  value = as.matrix(value[names, , drop = FALSE])
  normal = inherits(object, "dbglm") && object[["family"]][["family"]] %in% c("binomial", "poisson")
  colnames(value) = c("Estimate", "Std. Error", if (normal) "z value" else "t value", if (normal) "Pr(>|z|)" else "Pr(>|t|)")
  value
}

#' @rdname coefficient-summaries
#' @exportS3Method fixest::se
se.dbreg = function(object, ...) {
  value = coeftable.dbreg(object, ...)
  stats::setNames(value[, 2L], rownames(value))
}

#' @rdname coefficient-summaries
#' @exportS3Method fixest::tstat
tstat.dbreg = function(object, ...) {
  value = coeftable.dbreg(object, ...)
  stats::setNames(value[, 3L], rownames(value))
}

#' @rdname coefficient-summaries
#' @exportS3Method fixest::pvalue
pvalue.dbreg = function(object, ...) {
  value = coeftable.dbreg(object, ...)
  stats::setNames(value[, 4L], rownames(value))
}

#' Augment the immutable fitted sample
#' Collects the saved estimation rows into R, including SQL inputs, and appends
#' response-scale `.fitted` and `.resid` columns. Requires `retain="all"` or
#' `retain="inference"` and a live owning session.
#' @param x A fitted database regression.
#' @param data Optional original R input, in its original row order.
#' @param ... Additional arguments are not supported.
#' @return A data frame of estimation rows.
#' @importFrom generics augment
#' @exportS3Method generics::augment
augment.dbreg = function(x, data = NULL, ...) {
  check_unused_args(list(...))
  state = retained_inference(x)
  conn = fitted_resource_connection(x[["resource_state"]])
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  private = unlist(state[c("key", "estimated", "eta", "residual", "mu", "scores")], use.names = FALSE)
  columns = setdiff(names(state[["column_classes"]]), private)
  aliases = safe_identifiers(c("dbreg_augment_position__", "dbreg_augment_fitted__", "dbreg_augment_resid__"),
    taken = names(state[["column_classes"]]), fold_case = folds_identifier_case(conn)
  )
  fitted = if (inherits(x, "dbglm")) state[["mu"]] else state[["eta"]]
  values = backend_get_query(conn, paste(
    "SELECT", paste(c(
      q(columns),
      paste(q(state[["key"]]), "AS", q(aliases[[1L]])),
      paste(q(fitted), "AS", q(aliases[[2L]])), paste(q(state[["residual"]]), "AS", q(aliases[[3L]]))
    ), collapse = ","),
    "FROM", q(state[["table"]]), "WHERE", q(state[["estimated"]]), "ORDER BY", q(state[["key"]])
  ))
  if (is.null(data)) data = values[columns] else {
    if (!isTRUE(state[["r_input"]]) || !is.data.frame(data) || nrow(data) != state[["nobs_input"]]) {
      stop("`data` must be the original R input in its original row order; omit it to collect the saved sample.", call. = FALSE)
    }
    data = data[as.numeric(values[[aliases[[1L]]]]), , drop = FALSE]
  }
  data[[".fitted"]] = values[[aliases[[2L]]]]
  data[[".resid"]] = values[[aliases[[3L]]]]
  data
}

#' @rdname augment.dbreg
#' @export
#' @name augment
generics::augment

fixed_effect_labels = function(frame) {
  if (ncol(frame) == 1L) return(as.character(frame[[1L]]))
  vapply(seq_len(nrow(frame)), function(i) paste0("(", paste(vapply(frame, function(column)
    encodeString(as.character(column[[i]]), quote = "\""), character(1)), collapse = ", "), ")"), character(1))
}

#' Extract saved fixed effects
#' Returns FE intercepts and varying-slope coefficients in the saved identifying
#' normalization. The model constant is assigned to the first FE intercept.
#' Values are collected into R, including database-backed dictionaries.
#' @param object A native model fitted with `retain="prediction"` or `retain="all"`.
#' @param ... Additional arguments are not supported.
#' @return A named list of named numeric vectors. Composite levels use quoted
#'   tuples to keep their component values distinguishable.
#' @exportS3Method fixest::fixef
fixef.dbreg = function(object, ...) {
  check_unused_args(list(...))
  if (!length(object[["fe"]])) return(list())
  state = object[["fitted_state"]]
  if (is.null(state[["prediction"]])) stop("Fixed-effect extraction requires retain='prediction' or retain='all'.", call. = FALSE)
  specification = state[["specification"]]
  dictionaries = specification[["dictionaries"]]
  group_names = vapply(dictionaries, function(dictionary) paste(dictionary[["columns"]], collapse = "^"), character(1))
  database = identical(state[["prediction"]][["storage"]], "database")
  if (database) {
    conn = fitted_resource_connection(object[["resource_state"]])
    q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
    values = lapply(dictionaries, function(dictionary) backend_get_query(
      conn,
      paste(
        "SELECT", paste(c(paste0(
          "CAST(", q(dictionary[["keys"]]), " AS VARCHAR) AS ",
          paste0("k", seq_along(dictionary[["keys"]]))
        ), "gamma"), collapse = ","),
        "FROM", dictionary[["table"]], "ORDER BY id"
      )
    ))
    out = lapply(values, function(value) stats::setNames(
      as.numeric(value[["gamma"]]),
      fixed_effect_labels(value[setdiff(names(value), "gamma")])
    ))
    names(out) = group_names
    beta = state[["coefficients"]][["slopes"]]
  } else {
    operations = specification[["operations"]]
    labels = lapply(dictionaries, function(dictionary) fixed_effect_labels(prediction_decode_keys(dictionary)))
    coefficients = state[["coefficients"]]
    out = lapply(seq_along(operations), function(i) {
      group = operations[[i]][["group"]]
      indices = coefficients[["starts"]][[i]] + match(dictionaries[[group]][["id"]], coefficients[["ids"]][[i]])
      stats::setNames(coefficients[["joint"]][indices], labels[[group]])
    })
    names(out) = vapply(operations, function(operation) {
      group = group_names[[operation[["group"]]]]
      if (is.null(operation[["variable"]])) group else paste0(group, "[[", operation[["variable"]], "]]")
    }, character(1))
    for (i in seq_along(operations)) {
      operation = operations[[i]]
      if (is.null(operation[["variable"]]) || is.null(operation[["center"]])) next
      intercept = which(vapply(operations, function(other) is.null(other[["variable"]]) &&
        other[["group"]] == operation[["group"]], logical(1)))[[1L]]
      out[[intercept]] = out[[intercept]] - operation[["center"]] * out[[i]]
    }
    beta = coefficients[["joint"]][coefficients[["fe_n"]] + seq_along(specification[["design_center"]])]
  }
  constant = state[["prediction"]][["constant"]]
  if (is.null(constant)) constant = 0
  out[[1L]] = out[[1L]] + constant - sum(beta * specification[["design_center"]])
  out
}

#' Covariance matrices for insight and modelsummary
#' @param x A native database regression.
#' @param vcov Covariance name or cluster formula, as in [summary.dbreg()].
#'   `"HC1"` and `"vcovCL"` are accepted for downstream compatibility.
#' @param vcov_args Named list containing `cluster`, `ssc`, or `vcov_fix`.
#' @param ... Downstream display arguments, which do not change covariance.
#' @return A numeric variance-covariance matrix.
#' @exportS3Method insight::get_varcov
get_varcov.dbreg = function(x, vcov = NULL, vcov_args = NULL, ...) {
  if (is.null(vcov_args)) vcov_args = list()
  if (!is.list(vcov_args) || (length(vcov_args) && (is.null(names(vcov_args)) || anyDuplicated(names(vcov_args))))) {
    stop("`vcov_args` must be a named list.", call. = FALSE)
  }
  unknown = setdiff(names(vcov_args), c("cluster", "ssc", "vcov_fix"))
  if (length(unknown)) stop("Unsupported covariance arguments: ", paste(unknown, collapse = ", "), call. = FALSE)
  if (is.character(vcov) && length(vcov) == 1L) {
    if (tolower(vcov) == "vcovcl") vcov = if (is.null(vcov_args[["cluster"]])) "cluster" else NULL
    if (identical(tolower(vcov), "vcov")) vcov = NULL
  }
  do.call(vcov.dbreg, c(list(object = x, vcov = vcov), vcov_args))
}
