# DuckDB 2.0 alpha: CTE-heavy dbreg regression

Benchmark run on 2026-09-04 with
[`duckdb_cte_regression.py`](duckdb_cte_regression.py).

## Workload

- 10 million synthetic balanced-panel observations (500,000 units by 20
  periods), generated lazily from DuckDB's `range` function.
- Two continuous regressors and two fixed effects.
- The CTE-heavy sufficient-statistics SQL emitted for
  `dbreg(y ~ x1 + x2 | unit + time, strategy = "demean", vcov = "iid")`.
- Four DuckDB threads, one untimed warm-up, and two blocks of five and three
  timed repetitions. The second block reversed the engine order.
- In-memory execution with no file-system I/O. The timer includes SQL parsing,
  planning, execution, and fetching the single sufficient-statistics row. It
  excludes Python/DuckDB startup and dbreg's small R-side matrix solve.

Host: AMD Ryzen 7 PRO 6850U (8 cores/16 threads), 30 GiB RAM, Linux
7.1.9-200.fc44.x86_64.

## Results

| Engine | Python package | Median | Mean | SD | Min-max | Relative median |
|---|---:|---:|---:|---:|---:|---:|
| DuckDB v1.5.2 | 1.5.2 | 4.234 s | 4.237 s | 0.116 s | 4.069-4.374 s | 1.000 |
| DuckDB v2.0.0-alpha38615 | 1.6.0.dev365 | 4.357 s | 4.441 s | 0.229 s | 4.241-4.904 s | 1.029 |

Raw times, in seconds:

```text
v1.5.2:             4.0911 4.2531 4.0694 4.2143 4.2063 4.3117 4.3743 4.3728
v2.0.0-alpha38615:  4.9042 4.3503 4.3323 4.3634 4.6545 4.2411 4.4430 4.2424
```

The alpha was 2.9% slower by pooled median. This is small enough to describe as
no material speedup on this workload, especially for a laptop benchmark with
eight observations per engine.

The returned sufficient statistics agree to floating-point tolerance and imply
coefficients of 0.499999 and -0.300000 in v1.5.2. The alpha results imply the
same coefficients.

## Plan observation

The engines do not use identical physical plans. DuckDB v1.5.2 keeps the final
sums and two `COUNT(DISTINCT ...)` operations in one ungrouped aggregate. The
2.0 alpha introduces an internal aggregate CTE and derives the two distinct
counts through separate grouped branches. That rewrite did not reduce runtime
for this high-cardinality fixed-effect workload.

The headline v2.0 CTE improvement is a rewritten **recursive** CTE engine.
dbreg's strategy SQL uses ordinary, non-recursive CTEs, so the large recursive
CTE gains should not be expected to transfer automatically. Other v2.0
optimizer changes can still affect dbreg plans, as the internal aggregate CTE
here demonstrates.

## Reproduction

Use isolated Python environments containing the two indicated DuckDB package
versions, then run:

```sh
/path/to/duckdb-1.5.2/python \
  _benchmarks/duckdb_cte_regression.py \
  --rows 10000000 --threads 4 --warmup 1 --iterations 5

/path/to/duckdb-2.0-alpha/python \
  _benchmarks/duckdb_cte_regression.py \
  --rows 10000000 --threads 4 --warmup 1 --iterations 5
```

At the time of this benchmark, the official DuckDB R development package still
embedded v1.5.5. Therefore the v2.0 comparison uses the Python client that
embeds v2.0.0-alpha38615 and runs dbreg's exact regression SQL, rather than an
end-to-end `dbreg()` call.
