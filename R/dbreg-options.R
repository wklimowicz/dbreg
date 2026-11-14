#' Global options
#'
#' @description
#' **dbreg** reads a few global options, which can be set with [options()],
#' e.g. `options(dbreg.verbose = TRUE)`.
#'
#' @section General:
#' - `dbreg.verbose`: Default for the `verbose` argument of [dbreg()],
#'   [dbivreg()] and [dbbinsreg()]. Default `FALSE`.
#' - `dbreg.i_max_columns`: Positive whole-number cap on expanded `i()`
#'   columns; default `1000`.
#' - `dbreg.prediction_max_entries`: Positive entry budget for portable
#'   learned-effect dictionaries kept in R (`retain = "prediction"` or
#'   `"all"`); default `2e6`. Oversized one- or two-family additive
#'   dictionaries may remain in the active database session.
#'
#' @section Alternating projections:
#' Alternating projections (AP) recover exact fixed effects when the
#' `"demean"` strategy cannot demean analytically, e.g. for unbalanced or
#' weighted two-way FE, or three or more FE. Used by [dbreg()] and [dbivreg()].
#'
#' - `dbreg.ap_solver`: `"auto"` (default), `"moments"`, or `"iterative"`.
#'   `"auto"` tries bounded sufficient-statistic FE projection first on DuckDB
#'   and native Databricks, falling back to iterative AP (or SQL CG on DuckDB)
#'   when a fit exceeds internal capacity limits.
#' - `dbreg.ap_accelerate`: Use Irons-Tuck extrapolation after two AP cycles.
#'   Default `TRUE`. Set to `FALSE` for plain alternating projections.
#' - `dbreg.ap_max_iter`: Maximum number of AP cycles. Default `100`.
#' - `dbreg.ap_tol`: Convergence tolerance for the largest absolute and
#'   relative (weighted) group mean. Default `1e-10`. Large datasets may
#'   benefit from a looser tolerance (e.g. `1e-7`).
#' - `dbreg.ap_method`: How each iterative AP cycle is computed in SQL:
#'   `"join"` writes group means to a table and joins them back, while
#'   `"window"` chains cycles using window functions. Default `"window"` on
#'   Spark/Databricks and `"join"` elsewhere.
#' - `dbreg.ap_batch_size`: Number of cycles per table write for
#'   `"window"`. Default `3`.
#'
#' @section Fixed-effect projection & SQL CG:
#' Bounded fixed-effect projection and coefficient-space SQL conjugate gradient
#' (CG) solvers respect additional tuning and resource-budget options:
#'
#' - `dbreg.cg_max_iter`: Maximum inner SQL CG iterations (also the default for
#'   `dbglm_control(fe_max_iter = ...)`). Default `10000L`.
#' - `dbreg.cg_restart`: Restart interval for SQL CG iterations. Default `500L`.
#' - `dbreg.databricks_packet_bytes`: Target byte size for packetized Databricks
#'   ODBC fixed-effect summary transport (1 to 32000). Default `20000L`.
#' - `dbreg.projection_max_bytes`: Memory budget in bytes for bounded FE
#'   projection structures. Default `1024^3` (1 GiB).
#' - `dbreg.projection_max_rows`: Maximum summary rows transferred from SQL.
#' - `dbreg.projection_max_entries`: Maximum dense FE design entries (`2e6`).
#' - `dbreg.projection_max_cells`: Maximum scalar summary cells transferred
#'   (`16e6`).
#' - `dbreg.projection_max_upload_bytes`: Maximum byte size for uploaded FE
#'   coefficient lookup statements (`8 * 1024^2`).
#' - `dbreg.fe_dense_max`: Maximum number of levels for dense Schur elimination
#'   (`4096`).
#' - `dbreg.fe_cross_max`: Maximum cross-block entries in bounded FE systems
#'   (`5e7`).
#'
#' @examples
#' # Set an option for a block of code, then restore the previous value
#' op = options(dbreg.ap_tol = 1e-8)
#' options(op)
#' @name dbreg-options
#' @md
NULL
