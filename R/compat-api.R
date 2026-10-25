compat_formula = function(fml, fixef = NULL) {
  if (is.character(fml) && length(fml) == 1L && !is.na(fml)) fml = stats::as.formula(fml, env = baseenv())
  fml = Formula::Formula(fml)
  if (length(fml)[[1L]] != 1L || length(fml)[[2L]] > 2L) stop("Use one outcome and at most one fixed-effect block.", call. = FALSE)
  if (!is.null(fixef)) {
    if (!is.character(fixef) || !length(fixef) || anyNA(fixef) || any(!nzchar(fixef))) stop("fixef must be a vector of column names.", call. = FALSE)
    if (length(fml)[[2L]] > 1L) stop("Specify fixed effects in the formula or fixef, not both.", call. = FALSE)
    rhs = Reduce(function(a,b) call("+",a,b), lapply(unique(fixef), as.name))
    fml = Formula::Formula(stats::as.formula(call("~",fml[[2L]],call("|",fml[[3L]],rhs)), env = baseenv()))
  }
  environment(fml) = baseenv()
  fml
}

compat_family = function(family) {
  if (is.character(family) && length(family) == 1L && !is.na(family)) {
    family = switch(family, gaussian = stats::gaussian(), binomial = stats::binomial(), poisson = stats::poisson(),
      stop("Supported family names are gaussian, binomial and poisson; pass a family object to choose its link.", call. = FALSE))
  }
  if (is.function(family)) family = family()
  dbglm_family(family)
}

compat_covariance = function(vcov, cluster, se) {
  if (!is.null(se)) {
    if (!is.null(vcov)) stop("Specify vcov or se, not both.", call. = FALSE)
    vcov = se
  }
  if (!is.null(cluster)) {
    if (!is.null(vcov) && !(is.character(vcov) && length(vcov) == 1L && vcov %in% c("cluster", "twoway"))) stop("Conflicting vcov and cluster requests.", call. = FALSE)
    vcov = NULL
  } else if (is.null(vcov)) vcov = "iid"
  list(vcov = vcov, cluster = cluster)
}

compat_run = function(expr, warn, notes) {
  for (value in list(warn,notes)) if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("warn and notes must be TRUE or FALSE.", call. = FALSE)
  run = function() if (warn) force(expr) else suppressWarnings(force(expr))
  if (notes) run() else suppressMessages(run())
}

compat_fit = function(estimator, fml, data, family, vcov, weights, offset, subset, cluster, se, ssc,
  fixef, fixef.rm, lean, conn, table, path, strategy, factor_refs, row_id, fixef.tol, fixef.iter,
  glm.tol = NULL, glm.iter = NULL, start = NULL, verbose = 0) {
  if (!is.logical(lean) || length(lean) != 1L || is.na(lean)) stop("lean must be TRUE or FALSE.", call. = FALSE)
  if (estimator == "feols" && identical(strategy,"mundlak")) stop("Mundlak is a native dbreg() strategy, not an exact feols compatibility strategy.",call.=FALSE)
  if (!(is.numeric(verbose) || is.logical(verbose)) || length(verbose) != 1L || is.na(verbose) || verbose < 0 || !is.finite(verbose)) stop("verbose must be a nonnegative number or logical value.", call. = FALSE)
  fixef.rm = match.arg(fixef.rm, c("perfect_fit", "singletons", "infinite_coef", "none"))
  covariance = compat_covariance(vcov, cluster, se)
  policy = inference_ssc(ssc)
  fml = compat_formula(fml, fixef)
  if (inherits(data, "tbl_lazy")) {
    if (!is.null(table)) stop("Specify one input relation.", call. = FALSE)
    table = data
    data = NULL
  }
  source = compat_input_vectors(data, table, weights, offset, subset)
  normalized = dbglm_formula(fml)
  parsed = parse_regression_formula(normalized[["fml"]], allow_empty = TRUE, allow_i = TRUE)
  request = inference_vcov_request(list(fitted_state = list(inference = list(fe_components = parsed[["fe_components"]]))),
    covariance[["vcov"]], covariance[["cluster"]])
  if (length(request[["clusters"]]) > 2L || estimator == "feglm" && length(request[["clusters"]]) > 1L) {
    stop("Supported clustering is one-/two-way OLS or one-way GLM.", call. = FALSE)
  }
  if (length(request[["clusters"]])) {
    predicate = if (is.null(source[["subset"]])) TRUE else if (is.character(source[["subset"]])) as.name(source[["subset"]]) else source[["subset"]][[2L]]
    for (column in unique(unlist(request[["clusters"]]))) predicate = call("&",predicate,call("!",call("is.na",as.name(column))))
    source[["subset"]] = stats::as.formula(call("~",predicate), env = baseenv())
  }
  gaussian_ols = estimator == "feols" && (length(normalized[["factors"]]) || length(normalized[["slopes"]]) || !length(parsed[["xvars"]]))
  if (gaussian_ols && strategy != "auto") stop("These OLS designs use the shared Gaussian projector; strategy must be auto.", call. = FALSE)
  args = list(fml = fml, conn = conn, table = table, data = source[["data"]], path = path,
    weights = source[["weights"]], offset = source[["offset"]], subset = source[["subset"]],
    factor_refs = factor_refs, row_id = row_id, retain = if (lean) "inference" else "all", vcov = "iid", verbose = verbose > 0)
  if (estimator == "feols" && !gaussian_ols) {
    args[["strategy"]] = strategy
    args[["fixef_rm"]] = if (fixef.rm %in% c("perfect_fit", "singletons")) "singletons" else "none"
    if (!is.null(fixef.tol)) args[["fe_tol"]] = fixef.tol
    if (!is.null(fixef.iter)) args[["fe_max_iter"]] = fixef.iter
    fit = do.call(dbreg, args)
  } else {
    if (estimator == "feglm" && strategy != "auto") stop("GLM strategy selection uses the shared WLS projector; strategy must be auto.", call. = FALSE)
    args[["family"]] = if (gaussian_ols) stats::gaussian() else compat_family(family)
    args[["fixef_rm"]] = fixef.rm
    controls = list()
    for (pair in list(c("fe_tol", "fixef.tol"), c("fe_max_iter", "fixef.iter"), c("glm_tol", "glm.tol"), c("glm_max_iter", "glm.iter"))) {
      value = get(pair[[2L]])
      if (!is.null(value)) controls[[pair[[1L]]]] = value
    }
    args[["control"]] = do.call(dbglm_control, controls)
    args[["start"]] = start
    fit = do.call(dbglm, args)
    if (gaussian_ols) {
      fit[["fitted_state"]][["inference"]][["kind"]] = "ols"
      fit[["linear_model"]] = TRUE
    }
  }
  success = FALSE
  on.exit(if (!success) close(fit), add = TRUE)
  state = fit[["fitted_state"]][["inference"]]
  fit[["fitted_state"]][["inference"]][["default_ssc"]] = policy
  fit[["n_params_model"]] = state[["K_model"]]
  fit[["nobs_input"]] = state[["nobs_input"]]
  fit[["nobs_dropped"]] = state[["nobs_input"]]-state[["nobs"]]
  fit[["df_residual"]] = nobs(fit)-state[["K_model"]]
  fit[["compatibility"]] = list(version = 1L, estimator = estimator, lean = lean,
    vector_offset = source[["vector_offset"]], offset_column = source[["offset"]], fixef_rm = fixef.rm)
  class(fit) = unique(c("dbreg_fixest", class(fit)))
  fit[["compatibility"]][["loglik"]] = if (estimator == "feols") logLik.dbreg(fit) else logLik.dbglm(fit)
  fit[["loglik"]] = as.numeric(fit[["compatibility"]][["loglik"]])
  classes = class(fit)
  fit = summary.dbreg(fit, vcov = covariance[["vcov"]], cluster = covariance[["cluster"]], ssc = policy)
  class(fit) = classes
  if (lean) {
    close(fit)
    fit[["resource_state"]] = NULL
    fit[["fitted_state"]][c("prediction", "specification", "coefficients", "resources")] = list(NULL,NULL,NULL,NULL)
    fit[["data"]] = NULL
    fit[["all_coefficients"]] = NULL
    fit[["composite_levels"]] = NULL
  }
  success = TRUE
  fit
}

#' Fixest-shaped database regression entry points
#'
#' Adapters preserve native SQL fitting conventions and normalize the supported
#' fixest user API. Native [dbreg()] and [dbglm()] defaults are unchanged.
#' @param fml Regression formula; absorbed effects follow a pipe.
#' @param data R data frame or lazy database relation. The second positional
#'   argument is data, matching fixest. SQL inputs can use conn/table/path instead.
#' @param vcov,cluster,se Covariance request: IID/standard, HC1/hetero, or cluster
#'   formula/columns. OLS supports two-way clustering; nonlinear GLMs one-way.
#' @param weights,offset Column name, one-column formula, or numeric R-data vector
#'   of original length (a scalar is recycled). Positional vectors are rejected
#'   for SQL inputs. Vector offsets require an explicit prediction offset for
#'   newdata; default fitted-sample prediction uses the saved offset.
#' @param subset Logical column/predicate formula, or original-length logical R
#'   vector. R integer indices must be negative exclusions or increasing unique
#'   positive positions. Duplicate/reordered positions are explicitly unsupported.
#' @param ssc Small-sample policy list, compatible with fixest::ssc().
#' @param fixef Optional vector of absorbed column names, as an alternative to |.
#' @param fixef.rm Recursive exclusion policy; default perfect_fit. OLS removes
#'   singletons, while nonlinear families also remove perfect-prediction groups.
#' @param lean Release row and learned-FE state after fitting/inference. Default
#'   FALSE retains an immutable database sample. Cached likelihood/covariance
#'   survive release; sample extraction/post-fit recomputation need live state.
#' @param fixef.tol,fixef.iter Optional native projection tolerance/iteration limit.
#'   NULL preserves native numerical defaults.
#' @param conn,table,path Native database input arguments; see [dbreg()].
#' @param strategy Native OLS strategy. Gaussian varying-slope/factor-coercion
#'   designs use the shared Gaussian projector; other OLS keeps its fast strategies.
#' @param factor_refs,row_id Native factor references and typed source keys.
#' @param warn,notes Emit warnings/messages? FALSE suppresses the corresponding
#'   conditions throughout fitting and covariance construction.
#' @param verbose Nonnegative verbosity; positive values enable native diagnostics.
#' @param ... Unsupported arguments error before fitting.
#' @return A dbreg_fixest model inheriting its native dbreg/dbglm class. It is
#'   not a fixest object and does not emulate fixest internal slots.
#' @export
feols = function(fml, data = NULL, vcov = NULL, weights = NULL, offset = NULL, subset = NULL,
  cluster = NULL, se = NULL, ssc = NULL, fixef = NULL, fixef.rm = "perfect_fit", lean = FALSE,
  fixef.tol = NULL, fixef.iter = NULL, warn = TRUE, notes = TRUE, verbose = 0,
  ..., conn = NULL, table = NULL, path = NULL,
  strategy = "auto", factor_refs = NULL, row_id = NULL) {
  check_unused_args(list(...))
  compat_run(compat_fit("feols", fml, data, NULL, vcov, weights, offset, subset, cluster, se, ssc,
    fixef, fixef.rm, lean, conn, table, path, strategy, factor_refs, row_id, fixef.tol, fixef.iter, verbose = verbose), warn, notes)
}

#' @rdname feols
#' @param family Built-in family object/function, or gaussian/binomial/poisson
#'   name. Default Gaussian identity matches fixest's entry point.
#' @param glm.tol,glm.iter Optional native deviance tolerance/IRLS iteration limit.
#' @param start Optional native starting coefficient vector.
#' @export
feglm = function(fml, data = NULL, family = "gaussian", vcov = NULL, offset = NULL, weights = NULL,
  subset = NULL, cluster = NULL, se = NULL, ssc = NULL, fixef = NULL, fixef.rm = "perfect_fit",
  lean = FALSE, fixef.tol = NULL, fixef.iter = NULL, glm.tol = NULL, glm.iter = NULL, start = NULL,
  warn = TRUE, notes = TRUE, verbose = 0, ..., conn = NULL, table = NULL, path = NULL, factor_refs = NULL, row_id = NULL) {
  check_unused_args(list(...))
  compat_run(compat_fit("feglm", fml, data, family, vcov, weights, offset, subset, cluster, se, ssc,
    fixef, fixef.rm, lean, conn, table, path, "auto", factor_refs, row_id, fixef.tol, fixef.iter, glm.tol, glm.iter, start, verbose), warn, notes)
}
