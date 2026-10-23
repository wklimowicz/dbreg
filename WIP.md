# Headline Cloglog Benchmark & Databricks Optimization Analysis (`WIP.md`)

## 1. The Headline Model (`~/research/teachers-leaving-wltl/`)

From [`scripts/03-annual/core-03a-main-model-estimate.R`](file:///home/wojtek/research/teachers-leaving-wltl/scripts/03-annual/core-03a-main-model-estimate.R#L45-L68) and [`R/annual_sample.R`](file:///home/wojtek/research/teachers-leaving-wltl/R/annual_sample.R#L70-L76), the headline discrete-time exit hazard is the Britton–Propper (BP) specification:

```r
feglm(
  wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK,
  data = hdl,
  family = binomial("cloglog"),
  vcov = ~district
)
```

- **Sample scale**: **4,594,532** eligible teacher-year observations (**4,591,507** retained after pruning 3,025 observations across 119 singleton / all-zero outcome `spine_band_year` cells).
- **Fixed-effect dimensions**: 6 additive categorical fixed effects totaling **2,389 active levels** (**2,385** model FE parameters; $K_{\text{vcov}} = 2,090$):
  - `spine_band_year`: 2,001 retained levels (2,120 raw)
  - `district`: 296 levels (also the cluster variable)
  - `Experience`: 44 levels
  - `AGE`: 44 levels
  - `SEX`: 2 levels
  - `FTPTWK`: 2 levels
- **Bipartite edge count across all $\binom{6}{2} = 15$ FE pairs**: **305,235** non-empty pairs total (dominated by `spine_band_year × district` with 132,843 edges, `spine_band_year × Experience` with 67,660, and `spine_band_year × AGE` with 67,679).

---

## 2. Current Execution Times Across Engines

> [!NOTE]
> **Update (2026-10-22, after implementation):** the headline now runs in **44.8 s on RAM-only DuckDB**: 67 SQL statements, 9 IRLS iterations, method `moments`, 0 N-row tables written during IRLS. Results: β = 0.891456572 (Δ 7.2e-8 vs feglm), SE = 0.096624867 (Δ 2.8e-7), N = 4,591,507, K_model = 2,385, K_vcov = 2,090. Before this change it was ~60–90 s and ~330 statements. See Section 5 below for what changed and what remains. Databricks has not been re-timed because benchmarking is stopped.

| Engine / Checkout | Solver Path | Elapsed Time | Convergence & Numerical Parity |
| :--- | :--- | :--- | :--- |
| **`fixest::feglm`** | In-memory C++ alternating projections + IRLS ([`core-03a`](file:///home/wojtek/research/teachers-leaving-wltl/scripts/03-annual/core-03a-main-model-estimate.R)) | **~3.5 – 5.0 min** (~210–300 s) | Converges in 6 outer iterations in standard production run (all 6 models in `core-03a` took [1,390 s total](file:///home/wojtek/research/teachers-leaving-wltl/logs/scripts/03-annual/core-03a-main-model-estimate.log)). *(With 25 strict oracle linear-predictor refinements, `feglm` took ~1 h 18 min).* |
| **DuckDB (`dbreg::dbglm`)** | RAM-only SQL IRLS + `dbglm_project_moments` ([`full-bp-safe`](file:///home/wojtek/projects/dbreg/dbreg-glm/investigation/implementation-fe/full-bp-safe/research-summary.csv)) | **~1.0 – 1.5 min** (~60–90 s solver; ~2.5 min w/ Parquet load) | Converges in **10 IRLS iterations** (max score $3.73 \times 10^{-9}$). Matches `fixest` to $\Delta\hat{\beta} = 3.21 \times 10^{-9}$, $\Delta\text{SE} = 6.17 \times 10^{-10}$, 0 key mismatches. |
| **Databricks (`dbreg::dbglm`)** | Remote SQL Warehouse + fallback `dbglm_project_window` | **~45 – 90+ min (est.)** *(Halted / not run on 4.59M real panel)* | On synthetic 2-FE data ([`RESULTS.md`](file:///home/wojtek/projects/dbreg/_benchmarks/glm-databricks/RESULTS.md)): 100k cloglog took **303.1 s** (5.1 min); 1M logit took **325.1 s**, 1M probit took **647.8 s**, and 1M cloglog was interrupted after iteration 3. |

---

## 3. Why Is Databricks (and Multi-FE `dbglm`) So Slow Right Now?

Inspecting [`R/wls-state.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/wls-state.R), [`R/glm-projection-moments.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-projection-moments.R), and [`R/glm-backend.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-backend.R) reveals **four compounding bottlenecks**:

### Bottleneck 1: Databricks disables the bounded moments solver and falls back to full-row Window AP for $>2$ FEs
In [`R/wls-state.R:8-55`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/wls-state.R#L8-L55):
```r
if (is.null(cache$moment_enabled) && !isTRUE(backend_info$is_databricks)) {
  ...
  cache$moment_enabled = sum(vapply(count, as.numeric, numeric(1))*(1+lengths(slopes))) <= 4096
} else if (is.null(cache$moment_enabled)) {
  cache$moment_enabled = FALSE
}
```
- Even though `sum(fe_sizes) = 2,389 <= 4,096` and [`dbglm_project_moments`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-projection-moments.R#L56-L60) already contains Databricks JSON lookup code, `wls-state.R` explicitly sets `cache$moment_enabled = FALSE` on Databricks!
- And because `length(fe) == 6 != 2`, Databricks skips `bounded_pcg` and falls into [`dbglm_project_window()`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-backend.R#L89-L127).
- Inside `dbglm_project_window()`, **every IRLS step** runs an inner loop of 3-cycle window chunks. Each chunk executes **18 sequential `OVER (PARTITION BY fe_k)` full-table window shuffles** (6 FEs $\times$ 3 cycles) over 4.59M rows and materializes a new 4.59M-row temporary table!

### Bottleneck 2: Materializing $N$-row temporary tables at every IRLS iteration
In [`dbglm_fe_fit`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/wls-state.R#L152-L266):
- The linear predictor `__dbglm_eta` is stored as a physical column on the $N$-row table `working`.
- At **every outer IRLS iteration** (10–14 iterations), `dbglm_fe_working` executes `CREATE OR REPLACE TEMP TABLE ... AS SELECT * FROM working`, writing millions of rows to remote storage.
- On Databricks SQL Warehouse, just 10–15 `CREATE TEMP TABLE` calls alone add 100–200+ seconds of storage/planning overhead (even on 1,000 rows, physical table creation took 154s out of 211s in [`PERFORMANCE_LOG.md`](file:///home/wojtek/projects/dbreg/PERFORMANCE_LOG.md#L1760-L1775)).

### Bottleneck 3: 200+ Sequential ODBC Round-Trips per Fit
Even when `dbglm_project_moments` is used (on DuckDB, or if enabled on Databricks):
- For $M = 6$ FEs, [`dbglm_project_moments`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-projection-moments.R#L4-L27) issues:
  - **6 separate** `SELECT ... GROUP BY fe_i` queries,
  - **$\binom{6}{2} = 15$ separate** `SELECT ... GROUP BY fe_i, fe_j` queries,
  - **1 full-table `CREATE TEMP TABLE projected`** joining 6 CTEs back to 4.59M rows,
  - **6 separate** `SELECT ... FROM projected GROUP BY fe_i` residual check queries,
  - **1 query** in `dbglm_fe_fit` to compute `SUM(w * x * x)` and `SUM(w * x * z)` on `projected` (even though by Frisch–Waugh–Lovell this is simply `raw_moments - crossprod(B, coefficients)` already in R!),
  - **1 full-table `CREATE TEMP TABLE next_working`** joining `working` and `projected` on `__dbglm_key` (a 4.59M-row hash join!),
  - **1 evaluation query** in `dbglm_fe_evaluate` containing 6 `GROUP BY fe_i` subqueries in a `UNION ALL`.
- That is **~30 SQL queries (37 full-table scans) + 2 full 4.59M-row table materializations per IRLS iteration** (~350+ queries across the fit). On a remote Databricks warehouse where each round-trip costs 0.5–3.0s of planning/network latency, 300 sequential queries automatically cost **5 to 15+ minutes** before counting data movement.

---

## 4. A Much Faster Architecture for Databricks (and DuckDB)

We can make the headline 6-FE cloglog run in **~25–45 seconds on Databricks** and **~8–12 seconds on DuckDB** (faster than in-memory `fixest::feglm`) without weakening any numerical tolerances:

### Optimization 1: Parameter-Space IRLS (Zero $N$-Row Table Writes During IRLS)
When the fixed effects are solved via moments / normal equations (or Schur complement / PCG on FE levels), the linear predictor $\eta_i^{(t)}$ on row $i$ is **completely determined** by the parameter vector $\theta^{(t)} = (\beta^{(t)}, \gamma_1^{(t)}, \dots, \gamma_M^{(t)}) \in \mathbb{R}^{p + K_{\text{FE}}}$:
$$\eta_i^{(t)} = o_i + x_i' \beta^{(t)} + \sum_{m=1}^M \gamma_{m, g_m(i)}^{(t)}$$
- For the headline model, $p + K_{\text{FE}} = 1 + 2,389 = 2,390$ doubles (~19 KB).
- In each WLS step, solving the FE normal equations $H_{\text{FE}} \alpha = B$ gives the FE coefficients $\alpha_X$ (for $X$) and $\alpha_z$ (for $z$).
- By Frisch–Waugh–Lovell, the demeaned cross-products are already known in R without scanning or materializing a `projected` table:
  $$H_{\text{step}} = X' W X - B_X' \alpha_X, \qquad b_{\text{step}} = X' W z - B_X' \alpha_z$$
  $$\beta^{\text{cand}} = H_{\text{step}}^{-1} b_{\text{step}}, \qquad \gamma^{\text{cand}} = \alpha_z - \alpha_X \beta^{\text{cand}}$$
- Even step-halving ($\lambda \in (0, 1]$) is exact in parameter space:
  $$\beta^{(t+1)} = (1 - \lambda)\beta^{(t)} + \lambda\beta^{\text{cand}}, \qquad \gamma^{(t+1)} = (1 - \lambda)\gamma^{(t)} + \lambda\gamma^{\text{cand}}$$
- **Result**: `base` is materialized **once** before IRLS. During IRLS, **zero $N$-row temporary tables are ever created or written**. At each iteration, R broadcasts the 2,389-row $\gamma^{(t)}$ lookup table (`/*+ BROADCAST(g1, ..., g6) */`) and computes $\eta_i, \mu_i, w_i, z_i, s_i, \text{dev}_i$ on the fly in a streaming read of `base`.

### Optimization 2: Single-Pass `GROUPING SETS` Sufficient Statistics (28 Queries $\to$ 1 Query per Step)
Instead of issuing $6 + \binom{6}{2} + 6 + 1 = 28$ separate `GROUP BY` queries per IRLS step, combine evaluation and moment collection into **one single-pass `GROUPING SETS` query**:
```sql
SELECT
  GROUPING_ID(f1, f2, f3, f4, f5, f6) AS __gid,
  f1, f2, f3, f4, f5, f6,
  SUM(dev) AS deviance,
  MIN(mu) AS mu_min, MAX(mu) AS mu_max,
  SUM(bad) AS bad,
  SUM(w) AS h,
  SUM(w * x) AS bx,
  SUM(w * z) AS bz,
  SUM(score) AS s,
  SUM(w * x * x) AS xwx,
  SUM(w * x * z) AS xwz,
  SUM(x * score) AS xs
FROM working_view
GROUP BY GROUPING SETS (
  (),
  (f1), (f2), (f3), (f4), (f5), (f6),
  (f1, f2), (f1, f3), (f1, f4), (f1, f5), (f1, f6),
  (f2, f3), (f2, f4), (f2, f5), (f2, f6),
  (f3, f4), (f3, f5), (f3, f6),
  (f4, f5), (f4, f6),
  (f5, f6)
)
```
- **Verified benchmark on the actual 4,594,532-row teacher dataset**:
  Running this exact 22-grouping-set query in RAM-only DuckDB returns all **307,624** moment rows in **0.787 seconds**!
- Assembling the $2,389 \times 2,389$ matrix $H_{\text{FE}}$ and solving via Cholesky in R takes **~0.08 seconds**.
- Total queries for a 10-iteration cloglog fit drops from **~330 queries to ~12 queries** (1 per IRLS iteration + 1 final clustered `vcov` query).

### Optimization 3: Enable `moment_enabled` (and Analytical First-FE Absorption) on Databricks
1. **Enable `dbglm_project_moments` on Databricks** when $\sum L_m$ is within the cap (e.g., $\le 10,000$) and the total bivariate edges fit within `dbreg.duckdb_ap_max_edges`.
2. **Schur-Complement Absorption of the Largest FE**:
   When one FE has very high cardinality $L_1$ (e.g., $L_1 = 50,000$ schools or `spine_band_year` cells) while the remaining FEs ($f_2, \dots, f_M$) have small total cardinality $K_{-1} = \sum_{m=2}^M L_m$ (for BP, $K_{-1} = 296 + 44 + 44 + 2 + 2 = 388$!):
   - Aggregate by $f_1$ first or absorb $f_1$ analytically via $D_1' W D_1$ (diagonal!), forming the Schur complement only on the remaining $388 \times 388$ block:
     $$S = H_{-1, -1} - H_{-1, 1} H_{1, 1}^{-1} H_{1, -1} \in \mathbb{R}^{388 \times 388}$$
   - Then R only solves a **$388 \times 388$ dense linear system** (instantaneous) and back-substitutes $\gamma_1 = H_{1,1}^{-1}(B_1 - H_{1,-1}\gamma_{-1})$!
3. **Pre-IRLS Sufficient-Statistic Cell Collapse (When Applicable)**:
   If the regressors $X$, fixed effects, and cluster variables take $C \ll N$ distinct combinations (e.g., when outside wages are defined at the `district × AGE × SEX × FTPTWK × YEAR` level rather than individual postcodes), collapsing the $N$ binary rows once before IRLS to sufficient statistics `(cell_id, n_trials = COUNT(*), y_sum = SUM(y), x, fe1..feM, cluster)` reduces the working table size by $5\times$–$50\times$ with mathematically identical likelihood, gradients, Hessian, and cluster scores.


---

## 5. Implemented (uncommitted, `feat/glm`)

| Change | File |
| :--- | :--- |
| Single-pass `GROUPING SETS` FE moments; residual and refinement computed in R; returns FE coefficients, `cross_reduction` and `extra_reduction`; NULL under caps | [`glm-projection-moments.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-projection-moments.R) |
| Parameter-space IRLS: `working` is a view on `base` plus a 1-column FE lookup; evaluation and next FE moments in one query | [`wls-state.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/wls-state.R) |
| Moments projector enabled on Databricks for >2 FEs (2-FE still `bounded_pcg`; `sql`/cap fallback still `window`) | [`wls-state.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/wls-state.R) |
| FE ranking on distinct-key dictionaries; `retained` dropped early | [`dbglm.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/dbglm.R) |
| Pruning via grouped per-FE summaries instead of N-row windows | [`glm-sample.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/R/glm-sample.R) |
| New 4-FE Databricks-fixture parity and cap-fallback assertions | [`test_glm_databricks.R`](file:///home/wojtek/projects/dbreg/dbreg-glm-pr/inst/tinytest/test_glm_databricks.R) |

**Tests:** 1,259 GLM assertions pass, 0 failures.

**Remaining headline cost (44.8 s):** 18 fused evaluations × ~1.9 s ≈ 34 s (9 iterations plus line-search trials), plus ~10 s in R.

Next candidates:
1. Line-search trials could compute only deviance and per-FE scores (cheap), and collect pair moments only for the accepted step.
2. Absorb `spine_band_year` analytically and solve only the 388×388 Schur block in R.
3. Re-time on Databricks once benchmarking is authorised.


---

## 6. Remote Benchmark Evidence: fixest vs DuckDB vs Databricks SQL Warehouse

Tested on 2026-10-22 on `input_100000` (100,000 observations):
$$\text{y\_cloglog} \sim x_1 + x_2 \mid g + h + k_1 + k_2, \quad \text{family} = \text{binomial}(\text{"cloglog"}), \quad \text{vcov} = \sim h$$
- **Complexity**: 4 fixed effects totaling **5,230 levels** ($g=5,000, h=200, k_1=20, k_2=10$).
- **Previous Databricks status**: Halted / unrunnable because $5,230 > 4,096$ triggered window AP fallback, requiring hundreds of distributed window shuffles over 13 IRLS iterations.
- **New Databricks status with raised cap (15,000) and broadcast JSON lookup**: **118.6 s** (13 iterations, method `moments`, zero physical tables created during IRLS).

### Benchmark Comparison (100k rows, 4 FEs, 5,230 levels)

| Engine | Solver Path | Elapsed Time | Iterations | $\hat{\beta}_1, \hat{\beta}_2$ | SE($\hat{\beta}_1$), SE($\hat{\beta}_2$) | Agreement vs `fixest` |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`fixest::feglm`** | In-memory C++ AP + IRLS | **0.52 s** | 6 | `0.269809, -0.187194` | `0.005053, 0.004338` | *Reference* |
| **`DuckDB (dbglm)`** | Single-pass Moments + Parameter-Space IRLS | **78.2 s** | 13 | `0.269824, -0.187204` | `0.005053, 0.004338` | $\Delta\beta \le 1.4 \times 10^{-5}$ |
| **`Databricks SQL`** | Raised Cap Moments + Broadcast JSON Lookup | **118.6 s** | 13 | `0.269824, -0.187204` | `0.005053, 0.004338` | $\Delta\beta \le 1.4 \times 10^{-5}$ |

*(Databricks and DuckDB produce bit-for-bit identical estimates to 8 decimal places: $\hat{\beta}_1 = 0.2698236$, $\hat{\beta}_2 = -0.1872039$).*

---

## 7. Conclusions: Sustainability & Architectural Ceiling

### Why This Paradigm Has Hit an Asymptotic Ceiling
Orchestrating iterative mathematical solvers (IRLS) step-by-step from R over remote SQL queries has a hard physical ceiling:
1. **The Catalyst Planning / Network Tax**: Every remote SQL execution against a Databricks SQL Warehouse carries an irreducible **~1.0 to 1.5 second floor** (TLS transit, query parsing, Catalyst compilation, Photon plan generation, executor barrier synchronization). Across 13 IRLS iterations with 2 queries each, **~30–40 seconds is spent purely in latency** before computing any floating-point arithmetic.
2. **Dense Bipartite Pair Scaling**: Computing all $\binom{K}{2}$ pairs in SQL scales quadratically in the number of fixed effects. When total levels reach 5,000+, assembling and factoring a dense $5,230 \times 5,230$ Cholesky matrix in R adds 6–10 seconds of CPU time.
3. **Array Operations vs Hash Aggregation**: In-memory C++ (`fixest`) updates contiguous memory arrays in L1/L2 CPU cache in **0.52 seconds**. A SQL database must scan rows, evaluate hash functions, update aggregation hash tables, and join back—inherently orders of magnitude heavier per iteration.

### The Two Sustainable Long-Term Paths
For further order-of-magnitude improvements:
1. **Push the Entire IRLS Loop into the Cluster (Databricks Serverless / Python / Spark)**:
   - Instead of 25 round trips over ODBC, R makes **one HTTP call** submitting an isolated script or job to Databricks.
   - The cluster executes all 10 IRLS iterations in parallel in cluster RAM using native C++/Python/NumPy, returning only the final 3 KB coefficients and covariance matrix. Expected time: **< 5 seconds**.
2. **Cell-Collapsed Sufficient Statistics (The Pure SQL Win)**:
   - For binary exit hazard / logistic / Poisson models, 4.59 million rows with discrete covariates collapse into **200,000–300,000 unique covariate cells** with weights $N_{\text{cell}}$ and counts $Y_{\text{cell}}$.
   - Running IRLS on the collapsed table preserves exact mathematical equivalence while cutting per-iteration scan and aggregation time by $5\times$ to $10\times$.
