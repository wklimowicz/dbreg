#' Global options
#'
#' @description
#' **dbreg** reads a few global options, which can be set with [options()],
#' e.g. `options(dbreg.verbose = TRUE)`.
#'
#' @section General:
#' - `dbreg.verbose`: Default for the `verbose` argument of [dbreg()],
#'   [dbivreg()] and [dbbinsreg()]. Default `FALSE`.
#'
#' @section Alternating projections:
#' Alternating projections (AP) recover exact fixed effects when the
#' `"demean"` strategy can't demean analytically, e.g. for unbalanced or
#' weighted two-way FE, or three or more FE. Used by [dbreg()] and [dbivreg()].
#'
#' - `dbreg.ap_accelerate`: Use Irons-Tuck extrapolation after two AP cycles.
#'   Default `TRUE`. The same absolute/relative certificate applies; a step
#'   that does not reduce the residual is followed by plain AP. `FALSE` restores
#'   ordinary AP. Fitted objects report `ap_iterations`, `ap_accelerated_steps`,
#'   `ap_max_abs`, `ap_max_relative` and `ap_method` when iterative AP is used.
#'   Iterations count projection cycles, including two per extrapolation.
#' - `dbreg.ap_max_iter`: Maximum number of AP cycles. Default `100`.
#' - `dbreg.ap_tol`: Convergence tolerance for both the largest absolute
#'   (weighted) group mean and that mean relative to each projected column's
#'   weighted RMS. Columns below original-scale floating-point roundoff are
#'   treated as absorbed. Default `1e-10`. The relative check prevents a small
#'   but unfinished FE-only column from appearing identified after rescaling.
#'   Large datasets may need a looser tolerance (e.g. `1e-7`).
#' - `dbreg.ap_method`: How each AP cycle is computed in SQL. `"join"` writes
#'   group means to a table and joins them back, once per fixed effect.
#'   `"window"` uses window functions instead, chaining several cycles into a
#'   single table write. Default `"window"` on Spark/Databricks, where each
#'   table write is costly, and `"join"` elsewhere.
#' - `dbreg.ap_batch_size`: Number of cycles per table write for
#'   `"window"`. Convergence is checked after each batch, so up to
#'   `dbreg.ap_batch_size - 1` extra cycles may run. Default `3`.
#'
#' @examples
#' # Set an option for a block of code, then restore the previous value
#' op = options(dbreg.ap_tol = 1e-8)
#' options(op)
#' @name dbreg-options
#' @md
NULL
