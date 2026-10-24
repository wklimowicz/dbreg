#' Small-sample correction policy
#' @param K.adj,G.adj Apply parameter-count/group-count corrections?
#' @param K.fixef Count full, none, or nonnested FE intercept parameters.
#' @param K.exact Exact FE rank counting is unsupported; must be FALSE.
#' @param G.df Group correction: min or conventional for two-way OLS.
#' @param t.df Testing df: min, conventional, or a number at least one.
#' @param ... Unsupported arguments error.
#' @return A validated policy list compatible with supported fixest SSC requests.
#' @export
ssc = function(K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE, G.adj = TRUE,
  G.df = "min", t.df = "min", ...) {
  check_unused_args(list(...))
  value = inference_ssc(list(K.adj=K.adj, K.fixef=K.fixef, K.exact=K.exact, G.adj=G.adj, G.df=G.df, t.df=t.df))
  class(value) = "ssc_type"
  value
}

compat_coefficient_names = function(names, keep = NULL, drop = NULL, order = NULL) {
  for (value in list(keep,drop,order)) if (!is.null(value) && (!is.character(value) || anyNA(value))) stop("Coefficient selectors must be nonmissing regular expressions.", call. = FALSE)
  match = function(pattern, names) {
    pattern = sub("^!%", "%!", pattern)
    pattern = sub("^%", "", pattern)
    negative = startsWith(pattern,"!")
    if (negative) pattern = substring(pattern,2L)
    found = grepl(pattern,names,perl=TRUE)
    if (negative) !found else found
  }
  if (length(keep)) names = names[Reduce(`|`,lapply(keep,match,names=names))]
  if (length(drop)) names = names[!Reduce(`|`,lapply(drop,match,names=names))]
  for (pattern in rev(order)) {
    first = match(pattern,names)
    names = c(names[first],names[!first])
  }
  names
}

#' Select compatibility coefficients
#' @param object A compatibility model.
#' @param keep,drop,order Regular-expression selectors, including leading !
#'   negation. Matching uses public coefficient names.
#' @param collin Include omitted coefficients as missing values?
#' @param fe Native explicit FE coefficients when available.
#' @param ... Unsupported arguments error.
#' @return Named coefficient vector.
#' @export
coef.dbreg_fixest = function(object, keep = NULL, drop = NULL, order = NULL, collin = FALSE, ..., fe = FALSE) {
  check_unused_args(list(...))
  if (!is.logical(collin) || length(collin) != 1L || is.na(collin)) stop("collin must be TRUE or FALSE.", call. = FALSE)
  if (isTRUE(fe) && isTRUE(object[["compatibility"]][["lean"]]) && length(object[["fe"]])) stop("Explicit FE coefficients are unavailable after lean release.", call. = FALSE)
  value = coef.dbreg(object,fe=fe,collin=collin)
  value[compat_coefficient_names(names(value),keep,drop,order)]
}

#' Coefficient tables and components
#' @param object A fitted database regression.
#' @param keep,drop,order Public-name regular-expression selectors.
#' @param ... Supported summary covariance controls.
#' @return Matrix or named numeric vector. Native tables keep native column
#'   labels; compatibility tables use fixest's public labels.
#' @export
coeftable = function(object, ...) UseMethod("coeftable")

#' @rdname coeftable
#' @export
coeftable.default = function(object, ...) {
  if (inherits(object,"fixest") && requireNamespace("fixest",quietly=TRUE)) return(fixest::coeftable(object,...))
  stop("No coefficient-table method for this model class.", call. = FALSE)
}

#' @rdname coeftable
#' @export
#' @exportS3Method fixest::coeftable
coeftable.dbreg = function(object, keep = NULL, drop = NULL, order = NULL, ...) {
  value = summary(object,...)[["coeftable"]]
  value[compat_coefficient_names(rownames(value),keep,drop,order),,drop=FALSE]
}

#' @rdname coeftable
#' @export
#' @exportS3Method fixest::coeftable
coeftable.dbreg_fixest = function(object, keep = NULL, drop = NULL, order = NULL, ...) {
  value = coeftable.dbreg(object,keep,drop,order,...)
  normal = inherits(object,"dbglm") && object[["family"]][["family"]] %in% c("binomial","poisson")
  colnames(value) = c("Estimate","Std. Error",if(normal) "z value" else "t value",if(normal) "Pr(>|z|)" else "Pr(>|t|)")
  value
}

#' @rdname coeftable
#' @export
se = function(object, ...) UseMethod("se")
#' @rdname coeftable
#' @export
#' @exportS3Method fixest::se
se.dbreg = function(object, ...) {
  value = coeftable(object,...)
  stats::setNames(value[,2L],rownames(value))
}
#' @rdname coeftable
#' @export
se.default = se.dbreg
#' @rdname coeftable
#' @export
pvalue = function(object, ...) UseMethod("pvalue")
#' @rdname coeftable
#' @export
#' @exportS3Method fixest::pvalue
pvalue.dbreg = function(object, ...) {
  value = coeftable(object,...)
  stats::setNames(value[,4L],rownames(value))
}
#' @rdname coeftable
#' @export
pvalue.default = pvalue.dbreg
#' @rdname coeftable
#' @export
tstat = function(object, ...) UseMethod("tstat")
#' @rdname coeftable
#' @export
#' @exportS3Method fixest::tstat
tstat.dbreg = function(object, ...) {
  value = coeftable(object,...)
  stats::setNames(value[,3L],rownames(value))
}
#' @rdname coeftable
#' @export
tstat.default = tstat.dbreg

#' Covariance and testing degrees of freedom
#' @param x A fitted database regression.
#' @param type k (covariance parameter count), resid, or t (testing df).
#' @param vars Optional coefficient names to count instead.
#' @param vcov,cluster,ssc Supported post-fit covariance controls.
#' @param ... Unsupported arguments error.
#' @return A parameter count or degrees of freedom. fixest::degrees_freedom is
#'   not generic; use dbreg::degrees_freedom for these models.
#' @export
degrees_freedom = function(x, type, vars = NULL, vcov = NULL, cluster = NULL, ssc = NULL, ...) {
  check_unused_args(list(...))
  if (!is.null(vars)) {
    if (!is.character(vars) || anyNA(vars)) stop("vars must contain coefficient names.", call. = FALSE)
    return(length(intersect(vars,names(stats::coef(x)))))
  }
  type = match.arg(type,c("k","resid","t"))
  V = vcov.dbreg(x,vcov=vcov,cluster=cluster,ssc=ssc)
  K = attr(V,"df.K")
  if (is.null(K)) K = nobs(x)-df.residual(x)
  switch(type,k=K,resid=max(nobs(x)-K,1),t=if(!is.null(attr(V,"df.t"))) attr(V,"df.t") else max(nobs(x)-K,1))
}
