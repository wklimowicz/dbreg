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
#' - `dbreg.ap_max_iter`: Maximum number of AP cycles. Default `100`.
#' - `dbreg.ap_tol`: Convergence tolerance, i.e. the largest absolute
#'   (weighted) group mean left in any demeaned variable. Default `1e-10`. On
#'   large datasets, AP can stall around `1e-8` because of floating-point
#'   noise, so a looser tolerance (e.g. `1e-7`) or `strategy = "mundlak"` may
#'   be needed.
#' - `dbreg.ap_method`: How each AP cycle is computed in SQL. `"join"` writes
#'   group means to a table and joins them back, once per fixed effect.
#'   `"window"` uses window functions instead, chaining several cycles into a
#'   single table write. Default `"window"` on Spark/Databricks, where each
#'   table write is costly, and `"join"` elsewhere.
#' - `dbreg.ap_batch_size`: Number of cycles per table write for
#'   `"window"`. Convergence is checked after each batch, so up to
#'   `dbreg.ap_batch_size - 1` extra cycles may run. Default `3`.
#'
#' @section Projection solvers and budgets:
#' Budgets bound compact FE summaries and matrix allocations. Exceeding a
#' solver budget selects another supported projection method; it does not
#' collect the observation-level dataset into R.
#'
#' - `dbreg.ap_solver`: `"auto"` (default), `"moments"`, or `"iterative"` selects
#'   the compact-moment solver or iterative SQL projection.
#' - `dbreg.duckdb_ap_method`: DuckDB SQL projection method, `"join"` (default)
#'   or `"cg"` for the bounded two-FE PCG solver. The general SQL AP cycle
#'   method is controlled separately by `dbreg.ap_method`.
#' - `dbreg.duckdb_ap_max_edges`: Maximum compact summary rows. Default `1e6`.
#' - `dbreg.duckdb_ap_max_entries`: Maximum FE design entries. Default `2e6`.
#' - `dbreg.ap_moment_max_cells`: Maximum transferred summary cells.
#'   Default `8e6`.
#' - `dbreg.glm_fe_dense_max`: Maximum levels in a dense FE system.
#'   Default `4096`.
#' - `dbreg.glm_fe_cross_max`: Maximum dense FE cross-product entries.
#'   Default `5e7`.
#' - `dbreg.projection_max_bytes`: Projection matrix memory budget in bytes.
#'   Default `512 * 1024^2` (512 MiB).
#' - `dbreg.spark_ap_chunk_cycles`: SQL AP cycles per Spark batch. Default `3`.
#'
#' @section Design and retained prediction:
#' - `dbreg.i_max_columns`: Maximum total columns expanded by `i()` terms.
#'   Default `1000`.
#' - `dbreg.prediction_max_entries`: Matrix and dictionary entry budget for
#'   portable learned-effect prediction. Default `2e6`. Supported oversized
#'   one- or two-FE states can instead remain in their owning database session.
#'
#' @section GLM diagnostics and Databricks transport:
#' - `dbreg.glm_databricks_lookup`: `"json"` (default) or `"values"` for bounded
#'   SQL coefficient/key lookups. Both preserve double precision and UTF-8 keys.
#' - `dbreg.glm_databricks_array_max_entries`: Maximum entries in a dense
#'   coefficient-array lookup. Default `65536`; larger lookups use joins.
#' - `dbreg.glm_databricks_materialize_intermediates`: Materialize intermediate
#'   GLM relations for diagnostic comparisons. Default `FALSE`.
#' - `dbreg.glm_trace`: Emit periodic inner projection diagnostics.
#'   Default `FALSE`.
#'
#' @examples
#' # Set an option for a block of code, then restore the previous value
#' op = options(dbreg.ap_tol = 1e-8)
#' options(op)
#' @name dbreg-options
#' @md
NULL
