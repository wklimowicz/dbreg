# Running `dbreg` & `dbglm` on Databricks Free Edition

This guide documents how to connect and run regressions with `dbreg` and `dbglm` against a free-tier Databricks SQL Warehouse (such as the Serverless Starter Warehouse on AWS/Azure/GCP).

---

## 1. Prerequisites & Environment Setup

### ODBC Driver
The official Databricks ODBC driver (64-bit) must be installed. On this Linux system:
* Driver path: `/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so`

### Environment Variables
Store credentials in `~/.Renviron` or export them in your shell:

```bash
export DATABRICKS_HOST="your-workspace.cloud.databricks.com"
export DATABRICKS_TOKEN="dapi..."
export DATABRICKS_SQL_PATH="/sql/1.0/warehouses/YOUR_WAREHOUSE_ID"
```

---

## 2. Connecting from R

Use `odbc::databricks()` with `useNativeQuery = TRUE`:

```r
library(DBI)
library(odbc)
pkgload::load_all() # Load dbreg

con <- dbConnect(
  odbc::databricks(),
  driver = "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so",
  httpPath = Sys.getenv("DATABRICKS_SQL_PATH"),
  useNativeQuery = TRUE,
  DefaultStringColumnLength = 65536
)
```

Alternatively, connect with standard `odbc::odbc()`:

```r
con <- dbConnect(
  odbc::odbc(),
  Driver = "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so",
  Host = Sys.getenv("DATABRICKS_HOST"),
  Port = 443,
  HTTPPath = Sys.getenv("DATABRICKS_SQL_PATH"),
  SSL = 1,
  ThriftTransport = 2,
  AuthMech = 3,
  UID = "token",
  PWD = Sys.getenv("DATABRICKS_TOKEN"),
  useNativeQuery = TRUE,
  DefaultStringColumnLength = 65536
)
```

Verify backend detection:
```r
dbreg:::detect_backend(con)
# Should return: $name = "spark", $is_databricks = TRUE, $supports_physical_temp_tables = TRUE
```

---

## 3. Serverless Starter Warehouse Management & Free-Tier Gotchas

### Auto-Stop (10-minute timeout)
Free edition warehouses stop automatically after 10 minutes of inactivity. When a query is issued while stopped, startup takes ~5–15 seconds.

### Programmatic Warehouse Startup
To wake the warehouse before executing benchmarks:

```bash
curl -s -X POST -H "Authorization: Bearer $DATABRICKS_TOKEN" \
  "https://$DATABRICKS_HOST/api/2.0/sql/warehouses/YOUR_WAREHOUSE_ID/start"
```

Check state:
```bash
curl -s -H "Authorization: Bearer $DATABRICKS_TOKEN" \
  "https://$DATABRICKS_HOST/api/2.0/sql/warehouses/YOUR_WAREHOUSE_ID" | jq .state
# Returns "STARTING" or "RUNNING"
```

### Resource Gatekeeper: `DENY_NEW_AND_EXISTING_RESOURCES` (`denyReason: INACTIVE`)
If the REST API or ODBC connection fails with:
`HTTP Response code: 400, BAD_REQUEST: Cannot create the resource, please try again later.`
and the REST API details state:
`"reason": "DENY_NEW_AND_EXISTING_RESOURCES", "denyReason": "INACTIVE"`

* **Root Cause**: Platform-level suspension by Databricks resource-gatekeeper (e.g. trial expiration or pending account verification).
* **Fix**: Log in to the Databricks Account Console at [https://accounts.cloud.databricks.com](https://accounts.cloud.databricks.com) (not just the workspace URL) and confirm/reactivate the workspace or update verification details.

---

## 4. Recommended R Options for Databricks

When running regressions on Databricks, set:

```r
# AP convergence floor: avoid stalling at float serialization floor (~1e-8)
options(dbreg.ap_tol = 1e-7)

# Window AP chunking (cycles per query)
options(dbreg.spark_ap_chunk_cycles = 3L)
```

---

## 5. Usage Examples

### Single-Pass Linear Regressions (`moments` & `mundlak`)
Fastest execution path (~2.0–2.7s even on 5M rows):

```r
# No-FE OLS (single aggregation pass)
fit_ols <- dbreg(y ~ x1 + x2, conn = con, table = "dbreg_bench_100k")

# 2-FE Mundlak / Correlated Random Effects (single pass CTE)
fit_mun <- dbreg(y ~ x1 + x2 | fe1 + fe2, conn = con, table = "dbreg_bench_100k", strategy = "mundlak")
```

### Multi-Way Demeaning (`demean`)
Iterative alternating projections with cluster-robust standard errors:

```r
fit_fe <- dbreg(
  y ~ x1 + x2 | fe1 + fe2,
  conn = con,
  table = "dbreg_bench_100k",
  strategy = "demean",
  vcov = ~cl
)
```

### Generalized Linear Models (`dbglm`)
Parameter-space streaming IRLS without intermediate $N$-row table materialization:

```r
fit_glm <- dbglm(
  y_bin ~ x1 + x2 | fe1 + fe2,
  conn = con,
  table = "dbreg_bench_100k_glm",
  family = binomial("logit"),
  vcov = ~cl
)
```
