#!/usr/bin/env python3
"""Benchmark dbreg's CTE-heavy two-way demean regression SQL.

Run this script with separate Python environments to compare DuckDB engines:

    /path/to/python duckdb_cte_regression.py --rows 20000000

The query mirrors the SQL emitted by::

    dbreg(y ~ x1 + x2 | unit + time, strategy = "demean", vcov = "iid")

for a balanced panel. The synthetic input is a lazy ``range`` view, so the
benchmark needs no external data and does not include file-system I/O.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import statistics
import time

import duckdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=20_000_000)
    parser.add_argument("--periods", type=int, default=20)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--iterations", type=int, default=5)
    return parser.parse_args()


REGRESSION_SQL = """
WITH base AS (
    SELECT unit, time, y, x1 AS x1, x2 AS x2 FROM bench
),
unit_means AS (
    SELECT unit,
           AVG(y) AS y_u, AVG(x1) AS x1_u, AVG(x2) AS x2_u
    FROM base
    GROUP BY unit
),
time_means AS (
    SELECT time,
           AVG(y) AS y_t, AVG(x1) AS x1_t, AVG(x2) AS x2_t
    FROM base
    GROUP BY time
),
overall AS (
    SELECT AVG(y) AS y_o, AVG(x1) AS x1_o, AVG(x2) AS x2_o
    FROM base
),
demeaned AS (
    SELECT
        b.unit,
        b.time,
        (b.y - um.y_u - tm.y_t + o.y_o) AS y_tilde,
        (b.x1 - um.x1_u - tm.x1_t + o.x1_o) AS x1_tilde,
        (b.x2 - um.x2_u - tm.x2_t + o.x2_o) AS x2_tilde
    FROM base b
    JOIN unit_means um ON b.unit = um.unit
    JOIN time_means tm ON b.time = tm.time
    CROSS JOIN overall o
),
moments AS (
    SELECT
        CAST(COUNT(*) AS BIGINT) AS n_total,
        CAST(COUNT(DISTINCT unit) AS BIGINT) AS n_fe1,
        CAST(COUNT(DISTINCT time) AS BIGINT) AS n_fe2,
        SUM(CAST(y_tilde AS FLOAT) * CAST(y_tilde AS FLOAT)) AS sum_y_sq,
        SUM(CAST(x1_tilde AS FLOAT) * CAST(y_tilde AS FLOAT)) AS sum_1_y,
        SUM(CAST(x1_tilde AS FLOAT) * CAST(x1_tilde AS FLOAT)) AS sum_1_1,
        SUM(CAST(x2_tilde AS FLOAT) * CAST(y_tilde AS FLOAT)) AS sum_2_y,
        SUM(CAST(x2_tilde AS FLOAT) * CAST(x2_tilde AS FLOAT)) AS sum_2_2,
        SUM(CAST(x1_tilde AS FLOAT) * CAST(x2_tilde AS FLOAT)) AS sum_1_2
    FROM demeaned
)
SELECT * FROM moments
"""


def main() -> None:
    args = parse_args()
    if args.rows <= 0 or args.periods <= 0 or args.rows % args.periods:
        raise SystemExit("--rows must be positive and divisible by --periods")
    if args.threads <= 0 or args.warmup < 0 or args.iterations <= 0:
        raise SystemExit("threads/iterations must be positive and warmup non-negative")

    con = duckdb.connect(":memory:")
    con.execute(f"SET threads = {args.threads}")
    con.execute(
        f"""
        CREATE VIEW bench AS
        SELECT
            i // {args.periods} AS unit,
            i % {args.periods} AS time,
            sin(i * 0.001) + ((i % 17) - 8) * 0.01 AS x1,
            cos(i * 0.0007) + ((i % 23) - 11) * 0.01 AS x2,
            0.5 * x1 - 0.3 * x2
                + ((unit % 1000) - 500) * 0.001
                + (time - {args.periods / 2}) * 0.02
                + ((i % 29) - 14) * 0.005 AS y
        FROM range({args.rows}) AS t(i)
        """
    )

    result = None
    for _ in range(args.warmup):
        result = con.execute(REGRESSION_SQL).fetchone()

    timings = []
    for _ in range(args.iterations):
        started = time.perf_counter()
        result = con.execute(REGRESSION_SQL).fetchone()
        timings.append(time.perf_counter() - started)

    assert result is not None
    numeric_result = [float(value) for value in result]
    if not all(math.isfinite(value) for value in numeric_result):
        raise RuntimeError("Regression query returned a non-finite statistic")

    report = {
        "python_package": duckdb.__version__,
        "duckdb_engine": con.execute("SELECT version()").fetchone()[0],
        "platform": platform.platform(),
        "rows": args.rows,
        "periods": args.periods,
        "threads": args.threads,
        "warmup": args.warmup,
        "iterations": args.iterations,
        "timings_seconds": timings,
        "median_seconds": statistics.median(timings),
        "min_seconds": min(timings),
        "max_seconds": max(timings),
        "result": numeric_result,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
