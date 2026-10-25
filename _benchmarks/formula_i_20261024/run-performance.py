from pathlib import Path
import csv
import re
import statistics
import subprocess
import sys

root = Path('_benchmarks/formula_i_20261024')
database = sys.argv[1]
previous = root / (sys.argv[2] if len(sys.argv) > 2 else 'performance-before-eligibility')
previous.mkdir(exist_ok=True)
for path in root.glob('large-*'):
    if path.name != 'large-fit.R':
        path.rename(previous / path.name)
records = []
for mode in ['baseline', 'lean', 'full']:
    with (root / f'large-{mode}-run.txt').open('w') as log:
        result = subprocess.run([
            '/usr/bin/time', '-v', 'Rscript', str(root / 'large-fit.R'), database, mode
        ], stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(f'{mode}: {result.returncode}')
    with (root / f'large-{mode}.csv').open() as handle:
        data = list(csv.DictReader(handle))
    text = (root / f'large-{mode}-run.txt').read_text()
    peak = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', text)[1])
    records.append(dict(
        mode=mode,
        median_elapsed_seconds=statistics.median(float(row['elapsed']) for row in data),
        min_elapsed_seconds=min(float(row['elapsed']) for row in data),
        max_elapsed_seconds=max(float(row['elapsed']) for row in data),
        peak_rss_kib=peak,
        model_object_bytes=int(data[0]['model_object_bytes']),
        fit_read_calls=int(data[0]['fit_read_calls']),
        fit_execute_calls=int(data[0]['fit_execute_calls']),
        max_rows=int(data[0]['max_rows']),
        query_object_bytes=int(data[0]['query_object_bytes']),
        first_rss_after_close_kib=int(data[0]['rss_after_close_kib']),
        last_rss_after_close_kib=int(data[-1]['rss_after_close_kib'])
    ))
    print(f'{mode} completed', flush=True)
with (root / 'performance-summary.csv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=records[0].keys())
    writer.writeheader()
    writer.writerows(records)
