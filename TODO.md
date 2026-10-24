# wojtek-dev TODO

## Avoid the wide-LM fixed-effect transfer limit

**Status:** proposed; not implemented. Identified in the
[5M-row branch comparison](_benchmarks/branch_matrix_results.txt), with
[raw matrix results](_benchmarks/branch_matrix_summary.csv).

The ten-slope, six-FE LM produces 671,488 FE-pair rows and 2,389 FE-level
rows. The current transfer pads every pair row with eleven unused
outcome/regressor fields. That creates 10,108,155 rectangular cells,
exceeding `dbreg.ap_moment_max_cells = 8e6` and triggering slow iterative
demeaning. The benchmark confirmed the fallback on DuckDB; the same limit
applies to native Databricks.

The tested workaround for this model is:

```r
options(dbreg.ap_moment_max_cells = 16e6)
```

With that setting, the completed LM medians were 7.1 seconds on DuckDB and
54.9 seconds on Databricks, with strict fixest parity. This is a temporary
configuration workaround; the durable fix should reduce the transfer size.

### Proposed implementation

- [ ] Transfer the FE-pair weights separately from the FE-level
  outcome/regressor summaries. Pair rows need only their set identifier,
  two level indices and weight sum. For this fixture, the compact layout
  would contain about 2.72M cells rather than 10.11M, before any extra
  bookkeeping. Evaluate separate bounded result sets or a packed query
  using the existing SQL helpers; avoid one query per FE or FE pair.
- [ ] Unpack directly into the existing FE system and right-hand-side
  matrices. Do not recreate the wide, mostly-NULL rectangle in R. Reuse
  the shared projection machinery for OLS/IV and preserve the existing
  GLM fitting path, numerical tolerances and row certification.
- [ ] Enforce the combined transfer budget before fetching, including
  bookkeeping and overflow-probe rows. Preserve the independent row,
  FE-by-design, dense-system and coefficient-upload limits. Keep the
  default 8M-cell budget and the iterative fallback for genuinely
  oversized or numerically rejected systems.
- [ ] Make a capacity-triggered fallback explain its reason and applicable
  limit through the existing verbose output/diagnostics.

Start with [R/ap-moments.R](R/ap-moments.R),
[R/glm-projection-moments.R](R/glm-projection-moments.R) and
[R/glm-databricks.R](R/glm-databricks.R).

### Acceptance checks

- [ ] Add a small synthetic six-FE, ten-or-more-slope regression case with
  a budget between the compact and padded sizes. Assert that the bounded
  solver is selected and that coefficients, full clustered covariance,
  retained observations and projection certification match the reference.
  Exercise the shared helper on DuckDB and native Databricks.
- [ ] Check transfer budgets just below, exactly at and above the required
  size, including safe fallback and cleanup after failure. Extend existing
  weighted OLS/IV and identifier-type tests where the changed helper is
  exercised; retain GLM reference coverage. Do not add wall-clock assertions
  to the unit suite.
- [ ] Rerun the same synthetic 5M-row, ten-slope LM on both backends with
  the **default 8M cap**, checking strict fixest parity and bounded-solver
  selection. Rerun cloglog to check for regressions. Run engines sequentially
  and keep capped diagnostics distinct from completed fits.
- [ ] Record fresh timings and resource usage in a new benchmark report.
  Preserve the original branch-matrix results as historical evidence, and
  update the option documentation if its accounting changes.
