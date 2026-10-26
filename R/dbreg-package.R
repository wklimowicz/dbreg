#' @keywords internal
#' @importFrom utils globalVariables
"_PACKAGE"

# Silence R CMD check NOTEs for variables it cannot see
globalVariables(c(
  # GLM setup and execution fields copied through the inputs environment.
  "all_names",
  "all_x",
  "base",
  "beta",
  "cell",
  "compression",
  "control",
  "count_input",
  "covariance",
  "create",
  "design_center",
  "display",
  "drop",
  "factor_refs",
  "family",
  "fixef_rm",
  "from",
  "has_intercept",
  "keep",
  "n",
  "n_eligible",
  "normalized",
  "offset",
  "original_base",
  "parsed",
  "precomputed_fe_sizes",
  "prediction_keys",
  "pruning",
  "public",
  "resources",
  "retain",
  "retained_source",
  "row",
  "row_id",
  "sample_table",
  "separation",
  "setup",
  "snap_info",
  "spec",
  "start",
  "subset",
  "table",
  "trend_aliases",
  "vcov",
  "weight_scale",
  "x",
  # NSE in dbbinsreg / plot
  "bin",
  "x_mid",
  "u",
  "u2",
  "V_beta",
  "cb_lwr",
  "cb_upr",
  # Fields pulled from the `inputs` environment into local scope via env2env()
  "any_continuous",
  "caller",
  "cluster_var",
  "compress_nmax",
  "compress_ratio",
  "compression_ratio_est",
  "conn",
  "data",
  "data_only",
  "drop_missings",
  "endog_names",
  "endog_vars",
  "exog_names",
  "fe",
  "fml",
  "from_statement",
  "has_interactions",
  "instr_names",
  "instr_vars",
  "intercept",
  "path",
  "sql_only",
  "ssc",
  "strategy",
  "structural_fml",
  "union_names",
  "union_sql",
  "vcov_type_req",
  "verbose",
  "weights",
  "x_all_names",
  "xvars",
  "yvar",
  "z_all_names"
))

#' Operations on package load
#' @param libname library name
#' @param pkgname package name
#' @keywords internal
#' @noRd
.onLoad = function(libname, pkgname) {
  op = options()

  dbreg_op = list(
    dbreg.verbose = FALSE
  )

  to_set = dbreg_op[!(names(dbreg_op) %in% names(op))]
  if (length(to_set)) options(to_set)

  invisible()
}
