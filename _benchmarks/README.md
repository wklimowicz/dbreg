# Benchmark Suite (`benchmark-suite` branch)

Reproducible scripts for the DuckDB, Databricks, `fixest`, and upstream `dbreg` comparisons. All scripts run out of the box on DuckDB and optionally compare against upstream `dbreg` and/or Databricks when environment variables are set:

| Script | Workload | Key Envvars |
|:---|:---|:---|
| `bench_fixest_dbreg_10m.R` | 10M-row `fixest` + `dbreg` suite (Discrete 2-FE, Easy & Difficult 1–3 way FE OLS, Poisson/Logit 2-FE) | `DBREG_BENCH_N`, `DBREG_BENCH_THREADS`, `DBREG_UPSTREAM_LIB`, `DATABRICKS_SQL_PATH` |
| `bench_comprehensive.R` | 14-case DuckDB vs. Databricks parity suite (OLS, 1–3 FE, composite FE, varying slopes, compress, IV, GLM, binscatter) | `DBREG_BENCH_THREADS`, `DATABRICKS_SQL_PATH` |
| `bench_12gb.R` | 100M-row (12 GB raw) disk-backed out-of-core DuckDB vs. Databricks suite | `DBREG_BENCH_N`, `DBREG_BENCH_MEM_LIMIT`, `DBREG_BENCH_DIR`, `DBREG_UPSTREAM_LIB`, `DATABRICKS_SQL_PATH` |

### Quick Example (with Upstream Comparison)
```bash
# Install upstream dbreg into a separate library directory
mkdir -p /tmp/dbreg_upstream_src /tmp/dbreg_upstream_lib
git archive ef061a3 | tar -x -C /tmp/dbreg_upstream_src
R CMD INSTALL --no-multiarch --no-help --no-docs -l /tmp/dbreg_upstream_lib /tmp/dbreg_upstream_src

# Run 10M-row fixest + dbreg benchmark
DBREG_UPSTREAM_LIB=/tmp/dbreg_upstream_lib Rscript _benchmarks/bench_fixest_dbreg_10m.R
```
