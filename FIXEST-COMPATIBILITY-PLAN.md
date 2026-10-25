# Fixest regression replacement: implementation plan and handoff

Prepared 2026-10-24. This document is the implementation handoff for the Todoist
task **“Make dbreg a drop-in replacement for fixest OLS and GLM regression calls”**
(`6hqR4rv693M8FH44`). It describes the full target, existing foundations, work
already started, remaining implementation, and acceptance tests.

## Final supported scope — 2026-10-24 23:58 UTC

The supported regression subset is implemented. Current source passes **5,418
full local assertions**. Actual native Databricks passes **120 formula/prediction
checks and 118 adapter checks**, with zero failures. Four additional current-source
preflight checks prove that unsupported `i()` interactions reject an unusable
connection before database setup. Native source precedes only that final guard
and two local-only regressions; `audit-source.R` verifies the exact delta and all
archived sources. The guard changes no supported native model calculation.

The implementation includes fixest-shaped OLS/GLM entry points; fractional weights,
offsets and restricted subsets; sample/pruning identities; supported `i()` and
varying-slope formulas; learned-FE R/SQL prediction; retained-sample post-fit
inference; two-way OLS clustering; extraction/default conventions; and resource
lifetime/serialization contracts. Native compact defaults and ordinary OLS fast
strategies remain available. Mundlak stays a native alternative and is rejected
by the exact compatibility adapter.

The support/difference table is in [the introduction](vignettes/intro.qmd), with
an [installed contract](inst/FIXEST-COMPATIBILITY.txt) and
[executable example](inst/examples/fixest-compatibility.R). Explicit boundaries
include categorical/binning/interactive `i()` forms, slope-only FEs, external new
cluster mappings, arbitrary R expressions, unrelated fixest consumers and
oversized three-or-more-FE/varying-slope prediction recovery. These are declared
limits of the replacement subset. The separate compact wide-AP-transfer TODO
is unchanged.

The final million-row three-level-factor fixture takes median
**0.103s native compact,
0.740s adapter lean and
0.762s adapter full**, with at most three
rows returned per query. Shared-Gaussian `i()` OLS and full row retention have a
substantial cost for simple models. Peak RSS includes embedded DuckDB, and small
returned R objects do not imply small total memory. Full measurements, repeated
cleanup checks and limitations are in
[`_benchmarks/formula_i_20261024/results.txt`](_benchmarks/formula_i_20261024/results.txt).

Next is draft-PR review against `wklimowicz/dbreg:wojtek-dev` and CI. Final
package check has zero errors/warnings and the existing optional `binsreg`
cross-reference NOTE; its complete log is recorded with that evidence. The continuation checkpoints
below are historical; their remaining-work lists no longer describe current scope.

## 0. Continuation checkpoints — historical status

The fractional-weight increment has completed local and actual native Databricks
validation: **710 assertions on each engine, zero failures**. Its full local
package run passed **3,549 assertions**. Package check passed with zero errors,
zero warnings and the existing optional `binsreg` cross-reference NOTE. Native
execution took about 44 minutes. Durable evidence, including source identity,
is in `_benchmarks/fixest_fractional_20261024/`; see `inst/GLM-VALIDATION.txt`.

This increment is committed as **`b9a8614`**. A subsequent argument/compact-state
follow-up passes **3,594 local assertions across 46 nonempty files**, with the
same package-check result. Unsupported OLS offsets and post-fit covariance
arguments now error; legacy aliases and shared GLM tidying remain supported.
Saved linear formulas no longer capture caller/strategy environments. Evidence
for that follow-up is in `_benchmarks/fixest_argument_contract_20261024/`.

The observation-identity foundation now passes **59 shared assertions on each
engine** and **3,653 full local assertions across 47 nonempty files**, with the
same package-check result. OLS and GLM share typed key validation/export helpers;
requested OLS exports use an owned immutable fit snapshot. Evidence is in
`_benchmarks/sample_identity_20261024/`. Default fits do not retain these tables.
Performance benchmarking of default fitting versus optional identity/export
paths is next, before broader retained-state work.

That performance checkpoint is now complete on a synthetic million-row DuckDB
panel: 210 measured warmed fits passed value checks. Normal fits show no large
consistent slowdown; explicit key validation/export can roughly double simple
OLS cost and increased peak client RSS about 16% in this fixture. Evidence and
limitations are in `_benchmarks/compat_performance_20261024/results.txt`. Keep
native defaults compact; prediction dictionaries and optional row-sized inference
state must be separate capabilities. Learned-effect recovery is the next milestone.

Bounded portable GLM learned-effect recovery/newdata prediction is now implemented
behind `retain="prediction"`, with native retention off by default. The latest
full local run passes **3,773 assertions**; package check has zero errors/warnings
and the existing optional `binsreg` NOTE. New tests cover four GLM families,
varying slopes, composite/three-way/disconnected effects, joint-design
estimability, exact keys, serialization and extreme group-specific weights.
The first native run passed 108 checks and failed eight analytic tiny-weight
checks, exposing a fixed-decimal score floor. That summation path is corrected;
the completed native rerun passes **116 assertions with zero failures** in about
18 minutes. The implementation is committed as **`a22b086`**; all 116 R/test source
files in the local/native manifests match it. Evidence:
`_benchmarks/prediction_state_20261024/`.

A 60-fit warmed million-row follow-up verified coefficients, covariance and
likelihood before timing comparisons. Portable prediction retention added about
0.26s/14% to the two-FE cloglog fit (1.811s to 2.071s pooled medians); saved state
was 107,896 bytes. Peak client RSS was about 15% higher during fitting. These
are local measurements with material timing variation, not native guarantees.
Database-backed large prediction state, original-sample/default prediction
and post-fit inference remain unfinished; the OLS recovery checkpoint follows.

Bounded portable OLS learned-effect prediction is now implemented behind the
same explicit `retain="prediction"` option. It preserves the original estimator's
slopes, recovers nuisance effects from the immutable fit frame, and certifies
training RSS. Demean/compress/auto paths, exact BIGINT and adjacent-double keys,
disconnected/collinear designs, source removal and serialization are covered.
Mundlak retention is explicitly unsupported. The shared recovery and typed-key
helpers also preserve GLM behavior.

This checkpoint passes **4,087 full local assertions** and **149 native OLS
assertions** (434.71 seconds). The native GLM regression run passes **118**
prediction checks plus **3** snapshot checks. Its frozen source predates a final
OLS-only compressed-aggregate conversion fix and local-runner expansion; those
two file differences are documented in the evidence. Latest full/local/native
OLS source manifests match. Package check has zero errors/warnings and the
existing optional `binsreg` NOTE. Evidence: `_benchmarks/linear_prediction_20261024/`.

New recovery checks exposed existing compressed-fit precision problems: numeric
factor labels could merge adjacent doubles; default numeric DBI conversion could
round BIGINT FE keys; integer64 sufficient statistics could truncate fractional
RSS arithmetic. These are corrected, including exact compressed cluster-cell
matching. Both DuckDB BIGINT conversion modes and native ODBC are validated.

A 90-fit million-row follow-up matches names, coefficients, full covariance, N
and strategy against the pre-compatibility baseline. Current default medians
are 0.0905s for OLS, 0.171s for two-FE OLS and 0.0255s for compressed OLS.
Explicit retained prediction raises the latter two to 0.422s and 0.1865s;
portable dictionaries occupy 107,296 and 6,712 bytes. Peak process RSS rises
about 25-30% with retention. These warmed local measurements do not establish
native or disk-bound performance. Compact defaults remain unchanged.

The next retained-state increment adds a database-backed dictionary for one
absorbed FE (including typed composite keys), without varying slopes. It keeps
original fitted slopes, recovers weighted group residual means, and uses a
bounded within-X null system to check new-row estimability. Mathematical state
and live resource leases are separate. `close(model)` releases owned relations;
caller connections and sample exports survive. Small portable state is unchanged.
Large multi-FE recovery and keyed database prediction output remain pending.

The completed increment passes **4,174 local assertions** and **60 native
Databricks assertions** (382.63 seconds), with zero package-check errors/warnings
and the existing optional `binsreg` NOTE. Final local/native/performance source
manifests match. Native validation first caught a text-literal issue (59/60
passed); bounded JSON key transport fixes it, and the source/failure/rerun are
preserved. Checked resource release also retains failed drops for retry rather
than accepting the old best-effort helper's suppressed error.

The million-row/400,000-level fixture now keeps the dictionary in SQL under the
default budget: mathematical state is 5,792 bytes; prediction recovery queries
return at most one row during fitting. Six measured fits match coefficients,
full covariance and N. Retention adds about 0.91s (0.243s to 1.154s medians) and
whole-process peak RSS rises about 29% in this warmed local fixture. This
includes embedded DuckDB execution memory and is not a native performance gate.
Evidence: `_benchmarks/database_prediction_20261024/`.

Two-FE database-backed prediction is now implemented as well, without varying
slopes. It fixes the original slopes, recovers nuisance effects/regressor
projections in SQL, and uses component labels plus a bounded within-X null
system to distinguish estimable new combinations. It supports FE-only GLMs,
checks eta/RSS recovery, and cleans up nonconverged requests. Small portable
state and compact native defaults retain their existing behavior.

This increment passes **4,258 local assertions**, **77 native two-FE checks**
(1,717.93 seconds) and **60 native single-FE regression checks** (529.97 seconds),
with zero package-check errors/warnings and the existing optional `binsreg`
NOTE. All five final local/native/performance source manifests match. Evidence:
`_benchmarks/two_fe_prediction_20261024/`.

A SQL-only million-row, 50,000/49,997-level two-FE fixture has a known analytical
predictor using paired opposite noise. Six fits match beta/full covariance/N;
retained predictions match the analytical target. Mathematical state is 8,120
bytes; recovery queries return at most one row during fitting. Retention adds
about 5.1s (4.588s to 9.713s medians) and 18% whole-process peak RSS locally.
Native large-state performance/disk-history storage remain unmeasured.

Keyed database prediction output now has a public API: explicit prediction
row_id plus conn/table/path or newdata returns a dbreg_prediction relation.
It supports portable and database-backed learned effects, fitted factor coding,
all four GLM link/response kernels, typed keys and output leases. Collection is
explicit; outputs survive model close, while caller-owned exports survive their
result's close. The final full local gate passes **4,427 assertions**. Broad
native checks pass **137 + 4** assertions; a final response-helper/Poisson tail
correction passes **16** focused native checks. The broad run predates only
that response-helper change and fixture registration; archived source and
audit-source.R make the scope explicit. Package check has zero errors/warnings
and the existing binsreg NOTE. Evidence: `_benchmarks/sql_prediction_20261024/`.
The million-row keyed output matches its analytical target to 1.4e-10, with at
most one row collected per setup query. Median output creation is 0.155s in
this warmed local fixture; backend row storage remains required and native
large-output performance/history storage are unmeasured.

Retained sample and post-fit inference now pass **4,793 full local assertions**,
**268 broad native checks** and **35 focused final native numerical checks**.
Package check has zero errors/warnings and the existing binsreg NOTE. Explicit
`retain="all"` supports default/original-sample prediction, fitted values,
residuals, observation identities and post-fit IID/HC1/one-/two-way OLS or
one-way GLM covariance. Compact native defaults are unchanged. The broad native
source preceded final numerical/key fixes and namespace-qualified test calls;
archived source plus final focused coverage are documented in
`_benchmarks/retained_inference_20261024/`. Its million-row fixture transfers
only bounded score/covariance summaries; full retention increases whole-process
peak RSS about 21% versus prediction-only in this noisy local comparison.

Native OLS offsets now pass **4,849 full local assertions** and **44 native
checks**, with the same package-check result. All existing strategies share a
derived outcome; R/SQL prediction adds the offset once and full retention keeps
original R positions. Evidence: `_benchmarks/linear_offset_20261024/`.

Gaussian identity support now passes **5,071 full local assertions** and **204
native checks**, with zero errors/warnings in package check and the existing
binsreg NOTE. It reuses the shared FE/varying-slope solver with one certified
WLS update and matches Gaussian feglm reporting. The OLS IID `K.adj=FALSE`
denominator edge is also corrected. Evidence: `_benchmarks/gaussian_20261024/`;
its no-FE GLM route is materially slower than OLS moments, so keep existing OLS
fast paths. Compatibility entry points and shared subset/vector normalization
are the next milestone; prototypes are in SCRATCH pending integration.

Compatibility entry points now pass **5,215 full local assertions**, **99 native
checks** and package check with zero errors/warnings plus the existing binsreg
NOTE. They normalize vectors/subsets before filtering, preserve original IDs,
support varying-slope OLS and fit-time two-way linear inference, use response
prediction defaults, and expose selected extraction methods. Native defaults
remain unchanged. Lean uses inference-only state transiently, avoiding learned
effect recovery. The million-row TwoFE fixture takes 3.608s native compact,
4.111s lean and 8.474s full retention; full compatibility has a real resource
cost. Evidence: `_benchmarks/compat_api_20261024/`. The next major implementation
is `i()` normalization, followed by integrated acceptance/support documentation.
External new cluster mapping data and oversized 3+ FE/varying-slope prediction
remain explicitly unsupported boundaries; original-source cluster columns work.

Final adapter source passes **5,236 local**, **118 broad native** and **11
focused native** checks, with zero package-check errors/warnings and the existing
binsreg NOTE. Extra probes caught/fixed an SSC class-tag typo and ODBC Boolean
literal typing in default cluster predicates. Exact earlier sources are archived;
final manifests match and audit-source.R scopes historical performance evidence.
The entry-point milestone is complete for its declared subset. `g[x]` supports
group intercepts and slopes; slope-only `g[[x]]` remains explicitly unsupported.

The original handoff and interrupted-run evidence in sections 1 and 8 below are
historical, preserved for provenance. They no longer describe the latest gate.
Remaining dependencies are `i()` formula support and final integrated
acceptance/support documentation. Adapters provide fit-time two-way OLS; native
compact fit-time clustering remains one-way. New external cluster mappings are
explicitly unsupported; later original-source cluster columns work. Three or more large FE blocks
and large varying-slope designs still need a supported recovery boundary. Use dbreg's
existing conventions; consult fixest for unresolved user-facing compatibility
choices. Concrete API/lifetime decisions are recorded in
`inst/FIXEST-COMPATIBILITY.txt`. The complete replacement remains unfinished.

## 1. Original handoff state (historical)

- Working branch: **`feat/fixest-compat`**, created from `wojtek-dev` at
  **`00f2868`** (`docs: track compact fixed-effect transfer follow-up`).
- The working tree was clean before this work. The implementation and documents
  described below are currently **uncommitted**. No push or PR was made.
- The first increment adds fractional GLM weights and likelihood conventions.
  It is locally validated, but **native Databricks validation is incomplete**.
- No compatibility wrappers, learned-FE prediction, two-way clustering, or
  post-estimation clustering have been implemented in this increment.
- Implementation was stopped at the user's request for this handoff. Continue
  from these changes; do not restart from upstream main or replace the existing
  GLM/projection machinery.

### Changes already present

| File | Change |
|---|---|
| `R/dbglm.R` | Accept fractional nonnegative original weights; preserve retained-row counts; normalize uniformly small weights for numerical fitting; restore reported likelihood/deviance/IID covariance to original units. |
| `R/glm-family.R` | SQL aggregation for fixest/R binomial likelihood reporting, including ties-to-even rounding. |
| `R/glm-methods.R` | `logLik(model, type = "weighted")` exposes the unrounded objective; default `type = "fixest"` returns the compatibility likelihood. |
| `inst/tinytest/helper_glm_fractional.R` | Shared DuckDB/native fixtures: reference comparisons, exact retained IDs, repeated weights, SQL fallback, cleanup, rounding, analytic small-weight checks. |
| `inst/tinytest/test_glm_fractional.R` | Local runner for the new fixtures. |
| `inst/tinytest/test_glm_fractional_live.R` | Native Databricks runner for the same fixtures. |
| `.github/scripts/run-databricks-tests.R` | Includes the new native suite. |
| `inst/tinytest/test_glm_reference.R`, `test_glm_sample.R` | Replace obsolete fractional-weight rejection assertions with value comparisons. |
| `man/dbglm.Rd`, `man/dbglm-methods.Rd`, `NEWS.md` | Documentation for the new behavior. |
| `inst/FIXEST-COMPATIBILITY.txt`, `inst/GLM-DEVELOPMENT.txt`, `inst/GLM-VALIDATION.txt` | Compatibility contract and maintainer documentation. |

The weight normalization fixes a demonstrated numerical problem: uniform weights
of `1e-16` or `1e-30` could satisfy the old absolute convergence criteria after one
iteration, with coefficient errors around 0.07–0.14 in simple binomial examples.
When all retained weights are below one, the new implementation divides them by
a power of two before fitting. It records this divisor as `weight_scale`, restores
the original scale of likelihood/deviance and IID covariance, and leaves the
scale-invariant sandwich covariance unchanged. Power-of-two scaling also preserves
the binary weight values used for likelihood rounding. Score diagnostics refer
to the internal scale. Analytic tests cover both absorbed and unabsorbed fits.

### Validation actually completed

- R 4.6.0; fixest 0.14.0; DuckDB R package 1.5.2; tinytest 1.4.3.
- New local fractional-weight suite: **606 assertions passed**.
- Full local suite: **3,445 assertions across 44 nonempty files**, zero failures.
  The 606 new assertions are included in this total.
- Package build and `R CMD check --no-manual --no-vignettes --no-tests` passed
  with **zero errors, zero warnings, one NOTE** for the existing unavailable
  optional `binsreg` documentation cross-reference. Tests were run separately.
- Supplemental analytic probes recovered binomial solutions at uniform weights
  down to `1e-30`; the committed-to-be suite uses `1e-20` as its tiny-weight case.
- Native credentials and the ODBC driver were present. The sandboxed attempt
  could not resolve the warehouse host. An unsandboxed attempt then exposed a
  fixture-upload issue: the coefficient VALUES helper deliberately rejects NA.
  The fixture now inserts its missing values in SQL without weakening that helper.
- The subsequent native run reached the cloglog comparison cases, but was
  interrupted for this handoff. **There is no completed native pass or native
  assertion total to claim.** Progress messages alone do not establish parity.

Session-local evidence is in `/tmp/dbreg-fractional-{local,full,native,build,check}.log`,
`/tmp/dbreg-fractional-{local,full}.rds`, and `/tmp/dbreg.Rcheck/00check.log`.
These are temporary artifacts, not durable release evidence. Rerun checks as
needed and record completed native results before declaring this increment ready.

## 2. Required end state and scope boundary

The supported replacement calls must provide:

1. OLS, including linear probability models, and GLMs, especially binomial
   complementary-log-log, with absorbed fixed effects on DuckDB and Databricks.
2. Fixest-style formulas, categorical reference levels, interactions,
   group-specific slopes, offsets and weights, preserving specification and
   public coefficient names.
3. Arbitrary nonnegative fractional observation weights, including externally
   supplied posterior probabilities, with different weights on successive calls.
4. Link and response predictions that include learned absorbed effects, align
   with supplied observations, and include estimable rows omitted solely for
   zero fitting weights. Newdata prediction must not require outcomes.
5. Two-way clustered OLS covariance alongside existing single-cluster inference.
6. Coefficients, full covariance, SEs, coefficient tables, likelihood, observation
   counts, retained-observation information, degrees of freedom and convergence
   status, with supported prediction/summary/extraction and post-fit clustering.
7. Matching statistical conventions for supported calls, with explicit errors
   for unsupported calls and documented unavoidable compatibility differences.
8. Substitution without rewriting surrounding supported analysis code, while
   preserving existing dbreg functionality.
9. Database inputs and large prediction outputs without full panel collection
   into R; explicit resource requirements and object lifetime/persistence limits.

This is a **regression replacement only**. Data construction, weight calculation,
EM or other outer iterations, posterior updates, initialization selection,
acceleration, checkpointing and outer convergence remain caller-owned. Use
synthetic fixtures; there are no application-specific dependencies.

## 3. Existing foundations to reuse

The current branch already has substantial functionality beyond the older package
overview. Consult source and `inst/WOJTEK-DEV.txt`, not historical WIP claims alone.

- OLS frames, sufficient statistics, score passes and strategy dispatch in
  `R/dbreg.R`, `R/frames.R`, `R/strategies.R`, `R/fit.R`, and `R/vcov.R`.
- SQL Fisher-scoring/IRLS for binomial cloglog/logit/probit and Poisson log in
  `R/dbglm.R`, `R/glm-fit.R`, `R/glm-fe-fit.R`, and `R/glm-family.R`.
- Absorbed categorical/composite FEs, GLM group-specific slopes, sample pruning,
  rank decisions, numerical boundaries and convergence certification.
- Bounded FE-moment and sparse projection solvers, shared with OLS/IV through
  `R/ap-moments.R`; SQL fallback paths and native Databricks transport.
- Restricted formula expansion and fitted categorical metadata in
  `R/formula-terms.R`, `R/linear-design.R`, `R/glm-design.R`,
  `R/glm-sample.R`, and `R/sql_model_matrix.R`.
- GLM retained-key export via `row_id`/`sample_table`, resource ownership and
  cleanup, and compact objects that currently discard fit-owned relations.
- Independent fixest reference refinement and shared local/native test helpers.

Reuse these layers. Do not add a second independent formula engine, collect the
panel to call fixest internally, or replace the stable GLM core wholesale.

## 4. Compatibility/API contract — settle before broad implementation

### Entry points

Add explicitly documented fixest-shaped entry points, for example
`dbreg::feols()` and `dbreg::feglm()` or clearly named equivalent adapters.
The exact names are an implementation decision. Their argument order must support
ordinary `function(fml, data, ...)` calls: native `dbreg()` currently takes `conn`
as its second argument, so a simple function rename is insufficient.

Adapters should normalize supported arguments into the existing fit pipeline.
Keep native dbreg defaults and supported legacy arguments compatible. Decide
explicitly how adapters accept database relations and how in-memory numeric
weight/offset vectors map to stable row identity. Never attach a positional
vector to an unordered database relation without an explicit alignment contract.

Inventory supported downstream calls: `coef`, `vcov`, `se`, `coeftable`,
`summary`, `confint`, `logLik`, `AIC`/`BIC`, `nobs`, `df.residual`, retained-row
extraction, `predict`, and any selected fixest-specific extraction generics.
Check their actual dispatch behavior. Do not assign a fake `fixest` class to an
incomplete object, or promise every fixest consumer/internal slot will work.

Unknown arguments must fail before fitting. Current probes demonstrated silent
ignoring of OLS `offset = ~o` and `vcov(model, cluster = ~g)`; these are priority
compatibility hazards. Preserve intentional legacy aliases explicitly.

### Defaults and versioning

- Pin a primary fixest reference version and record it in validation artifacts;
  test a supported version range separately. Current development uses 0.14.0.
- Characterize ordinary default calls as well as explicitly aligned settings.
- Adapter prediction defaults should match fixest (`type = "response"`); native
  `predict.dbglm()` currently defaults to `"link"`.
- Match fixest's supported `fixef.rm` and small-sample policies. OLS currently
  uses a different `ssc` interface/default from GLM/fixest.
- Distinguish model parameter count, covariance-adjustment K and testing df.
- Enumerate supported formulas, families, inference requests and extraction
  methods. Broad fixest features such as arbitrary R functions, multi-model
  formulas, panel operators, arbitrary VCOV functions, and unrelated estimators
  need explicit support decisions; they are not implied by a working basic fit.

## 5. Implementation work packages

### A. Complete the fractional-weight increment

1. Review current changes and finish native validation before committing/releasing.
2. Verify all supported families, with/without FEs, IID/HC1/cluster covariance,
   offsets, zero weights, missing weights, tiny weights and repeated calls.
3. Extend sample tests as needed for weight-dependent recursive exclusions and
   loss/recovery of regressor or FE identification between calls.
4. Confirm changes to weights in a reused relation cannot reuse stale summaries,
   rank decisions or retained-row masks. Earlier fits must remain immutable.
5. Check scale invariance against analytic examples, not only fixest: tiny
   weights can also stress the oracle's stopping rules.
6. Keep fitting and reporting objectives separate where necessary.

**Likelihood decision already implemented:** for binary responses, fixest 0.14.0
uses R's binomial `family$aic`, which rounds weights for reported likelihood.
`logLik(model)` follows that convention; `logLik(model, type = "weighted")`
returns the unrounded weighted Bernoulli objective. Fitting, deviance and scores
never round weights. Standard AIC/BIC use the default reported likelihood.
For Poisson and integer-weight binomial fits the two likelihoods coincide.

**Acceptance:** complete named coefficient/full covariance/sample/objective
comparisons on both engines, independent scaling checks, unchanged existing
integer/unweighted results, and cleanup on invalid/nonconvergent fits.

### B. Shared fitted-state and observation-identity design

This is the main dependency for prediction and post-estimation inference.

Define one internal fitted-state contract containing:

- Original model specification and immutable fitted coding: public terms,
  internal aliases, factor levels/references, intercept policy, transformations,
  offsets, composite FE components and varying slopes.
- Learned slope coefficients and FE dictionaries/coefficients, with a documented
  normalization and internal-centering corrections.
- Retained observation identity, exclusions, sample counts and required fitting
  metadata, without accidentally capturing the original R data in environments.
- Sufficient information or owned relations for supported residual/score passes.
- Backend/resource ownership, cleanup, serialization and persistence metadata.

Recover/save effects on **all** solver paths: analytic demeaning, bounded moments,
SQL/AP fallbacks and compression. The GLM moment path already computes FE
coefficients during fitting, but the public result currently discards them.
Fallback paths may need a final reconstruction solve. OLS needs corresponding
recovery of the absorbed component from the training fit.

Use in-memory dictionaries when bounded; support database-backed dictionaries
when cardinality makes collection inappropriate. Specify which operations are
portable after disconnect and which require a live owned snapshot. Provide
explicit cleanup and, where supported, export/reopen behavior. A serialized DBI
connection pointer is not persistence. A source relation changing after fitting
must not silently change the fitted model.

**Acceptance:** fitted-state invariants across strategies, numerical FE-sum
reconstruction, bounded transfer, object independence across repeated calls,
resource cleanup on success/error, and tested lifetime behavior.

### C. Formula, argument and sample compatibility

Normalize into the existing shared design representation:

- Factors and references, including the supported fixest `i()` forms.
- Main effects, numeric/factor interactions, public term names and ordering.
- Supported no-intercept coding and restricted transformations.
- Additive and composite FEs; group-specific slopes, including OLS support.
  Decide and document support for `g[x]` versus slope-only `g[[x]]`.
- Offsets for both OLS and GLM, including supported formula/vector/column forms.
- Weight/subset forms, missingness, zero weights and recursive FE exclusions.
- Collinear terms: retained basis and optional NA-filled extraction behavior.

Capture the specification independently of global contrast options and later
data changes. Preserve typed keys and integer64 identifiers. Keep SQL identifiers
and public coefficient names separate; cover quoting and generated-name collisions.

For OLS offsets, residualize the outcome appropriately for fitting and add the
offset back for predictions; use the correct residuals for likelihood/covariance.
Do not silently choose a Mundlak/CRE specification as an exact absorbed-FE
replacement where equivalence has not been established. Keep its existing native
functionality, but gate compatibility strategy choices on statistical equivalence.

OLS retained IDs and default pruning need implementation, not merely aggregate
diagnostic counts. Maintain separate estimation and prediction eligibility.

**Acceptance:** supported calls use unchanged formulas/arguments apart from the
entry point/data binding, match names and samples exactly, and reject unsupported
forms before execution. Test default fixest calls and `fixef.rm="none"` separately.

### D. Learned-FE prediction and large outputs

Implement prediction from the immutable fitted specification and learned state.
Current absorbed OLS predictions recompute effects from supplied outcomes;
absorbed GLM prediction is unsupported. Both need substantive changes.

- Return link and response predictions including offsets and learned effects.
- Do not require outcomes or fitting weights in ordinary newdata prediction.
- Return identical values when an observation is predicted alone, in a batch,
  or after permuting other prediction rows; altering outcomes must have no effect.
- Predict zero-weight fitting omissions when their required components are
  estimable. Weight-zero exclusion does not imply prediction nonestimability.
- Define missing/unseen-level and disconnected-component behavior explicitly.
  Knowing each marginal FE level is insufficient for all new combinations;
  unidentified combinations must not acquire normalization-dependent predictions.
- Handle varying-slope identification and rank-deficient specifications.
- Distinguish `predict(model)` on the fitting sample, original-sample alignment,
  and `predict(model, newdata=...)`; characterize fixest's behavior for each.
- Preserve positions/NA values for in-memory inputs. For SQL inputs, return
  stable keys with the output; SQL does not supply an implicit ordering.
- Offer a lazy/keyed relation or database destination for large outputs, with
  explicit collection to R and documented output ownership.

Do not make full FE prediction intervals a hidden prerequisite: implement and
document supported requests, and reject unsupported uncertainty calculations.

**Acceptance:** observation-by-observation fixest link/response parity, without
outcomes, including zero-weight rows and complex FEs; row/batch invariance;
well-defined nonestimability; and no observation-sized fetch for DB output.

### E. Two-way OLS clustering and post-estimation inference

Extend the shared score pass and covariance machinery. Two-way cluster meat is
the contribution for cluster A plus B minus their intersection. Build intersection
groups from typed tuples, not collision-prone pasted strings.

Match supported fixest conventions for:

- Model/covariance parameter counts, FE nesting and rank policies.
- `K.adj`, `K.fixef`, `G.adj`, `G.df`, `t.df`, and any selected `K.exact` support.
- Cluster counts, small numbers of clusters, covariance symmetry and indefinite
  matrix repair/warnings. Do not silently clip matrices by a different rule.
- Testing distributions, p-values and confidence intervals.
- Cluster missingness and the distinction between fit-time and post-fit requests.

Reuse one covariance engine from `vcov`, `summary`, SE, table and CI methods.
Post-estimation requests need enough saved state to calculate scores on the
original fitted sample. If new cluster columns are supplied later, require a
validated key mapping to that sample or another explicit supported mechanism.
Investigate fixest's missing-cluster behavior and match the selected policy;
do not silently refit on a different sample.

The task requires two-way clustering for linear models. GLM two-way clustering
can remain explicitly unsupported unless deliberately added to the contract.

**Acceptance:** full covariance matrices and downstream tests/CIs match fixest
for nested/crossed clusters, repeated/intersection clusters, weighted fits,
unbalanced panels and supported post-fit requests. Existing single-cluster
OLS/GLM/IV behavior remains covered.

### F. Results, extraction and failure semantics

Provide consistent metadata and methods for the agreed public contract:

- Coefficients, public names, full covariance, SEs and coefficient tables.
- `summary`, `confint`, `logLik`, AIC/BIC, `nobs`, `df.residual` and selected
  fixest-style extraction methods. Native OLS currently lacks a useful summary
  method and a `logLik` method.
- Input, eligible, pruned and retained counts; exact retained identities;
  original-sample alignment where requested.
- Model df, covariance-adjustment df, testing df and collinearity information.
- Outer/inner convergence, iterations, score/projector diagnostics and finite-MLE
  status. Failed/separated fits must not masquerade as ordinary successful fits.

Preserve existing `tidy`, `glance`, printing, IV and binscatter behavior where
shared helpers change. Define compact/lean modes truthfully: methods requiring
discarded state should error clearly, not silently ignore their arguments.

### G. Resources, backend integration and documentation

Keep observation-sized computation and output in SQL. Budget R memory for full
slope covariance (at least O(p²)), FE dictionaries, FE graph edges and bounded
summary/solver matrices. Bounded transfer does not automatically mean constant
memory. DuckDB is embedded: its memory contributes to total client-process RSS.

Document backend memory, spill/temp storage, round trips, retained tables and
connection lifetime. Exercise resource caps just below/at/above boundaries and
force genuine fallback paths. Do not bypass budgets by raising defaults without
evidence. The existing compact FE-transfer follow-up in `TODO.md` remains relevant
to wide OLS; coordinate changes to shared moment transport with that task.

Provide examples for a supported replacement call, repeated externally supplied
weights, prediction without outcomes, keyed database output, post-fit clustering,
and save/reopen/release behavior. Include a support/difference table.

## 6. Validation architecture and acceptance gates

### Primary oracle and independent checks

Use fixest as the primary oracle. Reuse `helper_glm_reference.R`, which refines
fixest from its own fitted predictor until coefficients and covariance stabilize.
Never initialize the oracle from dbreg results. Record reference version,
settings, convergence and any failure to stabilize.

Use explicit absolute-plus-relative tolerances, named before reviewing failures.
Existing helpers use coefficient tolerances around `1e-7 + 1e-6*abs(reference)`
and covariance tolerances around `1e-8 + 1e-5*abs(reference)`; assess them for each
fixture rather than assuming a universal threshold. Match IDs/names/policies
exactly. Investigate discrepancies before changing tolerances.

Independently check small explicit-dummy `lm`/`glm` fits, analytic saturated-cell
solutions, weighted likelihood/deviance, information/score equations and direct
sandwich calculations. Compare estimable predictions/FE sums when nuisance FE
normalizations differ. Do not compare iteration counts as numerical equality.

### Fixture dimensions

| Dimension | Required coverage |
|---|---|
| Estimator | OLS, LPM, cloglog, logit, probit, existing Poisson |
| Fixed effects | None, one, two, several; balanced/unbalanced; composite; varying slopes; disconnected/nested structures |
| Design | Continuous/categorical, references, interactions, no intercept where supported, collinearity, unsafe names, integer64 keys |
| Weights | Unweighted, integer, fractional, posterior-like, zeros, tiny scales, invalid inputs, changing columns/values across calls |
| Outcomes/sample | Missingness, rare outcomes, singleton chains, pure groups, separation/quasi-separation, finite-MLE controls |
| Inference | IID, HC1, one cluster; two-way OLS; relevant SSC options; fit-time and post-fit requests |
| Prediction | Link/response, no outcomes, zero-weight rows, missing/unseen/nonestimable rows, reordered/subset/batched input, DB destination |
| Execution | Each supported strategy, bounded solvers and forced fallbacks, DuckDB and actual Databricks |
| Lifetime | Repeated calls, source mutation, cleanup on failure, serialization, disconnect, release/reopen |

Use several fixed seeds and a deliberate selection of interacting cases, rather
than an unmanageably large Cartesian product. Add invariance checks for row
permutation, bijective FE relabeling, regressor units and uniform weight scaling.

### Test tiers

1. Small analytic and formula/contract tests: fast local execution.
2. Medium differential fixtures: full values, sample IDs and prediction parity.
3. Native Databricks: shared fixtures, actual SQL/transport/lifetime behavior;
   credentials missing or zero assertions is not a pass.
4. Large release checks: synthetic million-row panels and large prediction
   outputs; measure total client RSS, database resources, transfer volume,
   query count and repeated-call growth. Compare completed fits only.
5. Existing whole-package tests and package checks, including shared OLS/IV code.

A feature is complete only when its supported API, default statistical behavior,
native execution, failure semantics and documentation pass their relevant gates.
Existing historical benchmark results do not validate newly changed code.

## 7. Suggested delivery order and effort

| Workstream | Focused engineering days, including targeted tests |
|---|---:|
| API/support contract and adapters | 2–3 |
| Formula and sample compatibility | 4–7 |
| Fractional weights and likelihood audit | 3–5 |
| Fitted effects, lifecycle and prediction | 7–12 |
| Two-way and post-estimation inference | 4–7 |
| Result/extraction methods | 3–5 |
| Integrated native/resource validation and documentation | 6–10 |

Original total estimate: **about 30–50 engineering days (6–10 focused weeks)**
for one experienced developer on this integration branch. It covers the supported
regression subset, not unrestricted fixest emulation. Upstream discussion/review,
warehouse availability and release coordination are additional uncertainties.
The first weight increment has now been started; do not interpret the estimate
as a fresh 30–50 days after this handoff.

Deliver as reviewable increments:

1. Finish/review/checkpoint fractional weights and the compatibility contract.
2. Define and implement shared fitted state and stable observation identities.
3. Add formula/argument/sample normalization and compatibility entry points.
4. Add learned-FE prediction on both engines, including database output.
5. Add two-way OLS and supported post-estimation covariance.
6. Finish summary/extraction/default compatibility and lifecycle documentation.
7. Run integrated acceptance/resource checks and prepare focused PRs.

Items 2–3 can inform each other; prediction and post-fit covariance must share
the resulting state design. Keep existing native APIs compatible throughout.
Follow the repository's PR-only contribution policy; do not push directly to main.

## 8. Immediate next actions and reproduction

1. Read this plan and inspect `git diff` plus untracked files on
   `feat/fixest-compat`. Preserve the existing working changes.
2. Finish the native fractional suite, diagnosing any slow/hung query or fixture
   issue before broadening implementation. Its interrupted run is not evidence
   of a completed native pass. Consider adding per-phase diagnostics and a
   bounded timeout to make expensive fallback validation easier to monitor.
3. Review the likelihood convention and weight normalization carefully; inspect
   the full matrix and analytic tests, not just the assertion total.
4. Save durable validation evidence and commit the first increment once reviewed.
5. Begin the shared fitted-state design and explicit API compatibility tests.

Development must use `pkgload::load_all()`, not an installed `library(dbreg)`.

```r
pkgload::load_all()
tinytest::run_test_file("inst/tinytest/test_glm_fractional.R")
tinytest::test_all(".")

# Configured native-query connection/driver required. This can be slow because
# the suite exercises repeated calls and the forced SQL projection fallback.
Sys.setenv(DBREG_TEST_DATABRICKS = "true")
tinytest::run_test_file("inst/tinytest/test_glm_fractional_live.R")
```

```sh
# All native suites, with nonempty-result checks and saved artifacts:
make test-databricks

# Standard documentation/check workflow:
make document
make check
```

The local roxygen2 version used during this work was older than the repository's
recorded documentation tool version and regenerated unrelated files. Those
unrelated changes were removed. Keep future generated documentation diffs scoped.

Useful references: [FE OLS](https://lrberge.github.io/fixest/reference/feols.html),
[FE GLM](https://lrberge.github.io/fixest/reference/feglm.html),
[prediction](https://lrberge.github.io/fixest/reference/predict.fixest.html),
[small-sample corrections](https://lrberge.github.io/fixest/reference/ssc.html),
and [covariance methods](https://lrberge.github.io/fixest/reference/vcov.fixest.html).
Use the installed pinned version as the executable reference; web documentation
may describe a newer release.
