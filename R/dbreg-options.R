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
#' DuckDB and native Databricks first try bounded FE marginal/pair summaries
#' with a shared Schur solve. Original typed keys stay in SQL. Reconstructed
#' rows must satisfy the same absolute/relative group equations, finite-value
#' and row-count checks as the iterative path. Capacity/numerical declines are
#' reported in `projection_diagnostics`; unexpected database failures propagate.
#'
#' - `dbreg.ap_solver`: `"auto"` (default), `"moments"` or `"iterative"`.
#'   Eligible additive FE OLS/IV fits use sufficient-statistic FWL over one
#'   owned source snapshot, keeping transformed observations inline. Global
#'   cross-products share the transfer/memory caps; near-absorption cancellation
#'   triggers an actual projected-product aggregate. Numerical/capacity declines
#'   retain ordinary AP fallback and the unchanged row certificates.
#'   Auto tries bounded moments unless `dbreg.ap_method` is explicitly set.
#'   Moments also tries them with an explicit method; an expected decline uses
#'   iterative AP. Iterative bypasses the bounded solver. SQL-only output shows
#'   the equivalent ordinary AP plan.
#' - `dbreg.projection_max_bytes`: Planned matrix/buffer budget; default 1 GiB.
#'   This bounds representations and allocations, not process peak RSS.
#' - `dbreg.projection_max_rows`: Cumulative marginal/pair summary rows;
#'   default `floor(dbreg.projection_max_bytes / 512)`.
#' - `dbreg.projection_max_entries`: FE-level by design entries; default 2e6.
#' - `dbreg.projection_max_cells`: Cumulative decoded scalar cells, including
#'   overflow probes; default 16e6. SQL LIMIT probes enforce runtime caps.
#' - `dbreg.projection_max_upload_bytes`: Coefficient SQL byte budget;
#'   default 8 MiB, also bounded by one eighth of the matrix/buffer budget.
#' - `dbreg.databricks_packet_bytes`: Target uncompressed summary packet bytes;
#'   default 20000, whole numbers from 1 to 32000. Small-entry packets stay
#'   below twice the target. Larger entries get individual packets; entries
#'   at or above 65535 bytes signal a capacity decline before ODBC transfer.
#'   Server byte/count metadata detects truncation, including a lost final
#'   numeric digit. Malformed/truncated packets fail instead of fitting.
#' - `dbreg.fe_dense_max`: Non-absorbed FE levels permitted in a dense Schur
#'   system; default 4096. Larger two-family systems try bounded Schur PCG.
#' - `dbreg.fe_cross_max`: Cross-block/fill entry limit; default 5e7.
#'   Dense/cross limits may be zero; other limits must be positive integers.
#'   Legacy aliases `dbreg.duckdb_ap_max_edges`, `dbreg.duckdb_ap_max_entries`,
#'   `dbreg.ap_moment_max_cells`, `dbreg.glm_fe_dense_max` and
#'   `dbreg.glm_fe_cross_max` are accepted, respectively, when their canonical
#'   options are unset.
#'
#' FE marginals and pairs share one summary statement. DuckDB transports
#' compact list columns; native Databricks transports bounded uncompressed
#' packets and uploads coefficients as numeric JSON rows. Original typed
#' source keys stay in SQL. `ap_transport`, `ap_summary_statements`,
#' `ap_wire_rows` and (for packets) `ap_wire_bytes` describe the representation;
#' logical transfer limits still apply to decoded values, not just wire rows.
#'
#' Successful bounded fits report `ap_method="bounded_moments"`, row/cell/
#' upload counts, planned bytes and the solver used. `ap_iterations` counts
#' solve/refinement steps, capped by `fe_max_iter` and eight steps. Inner PCG
#' has a separate 2,000-iteration cap per solve and reports total iterations
#' and its preconditioner. Explicit methods preserve ordinary AP cycle counts.
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
