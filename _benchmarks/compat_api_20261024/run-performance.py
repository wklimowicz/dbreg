from pathlib import Path
import subprocess
root=Path('_benchmarks/compat_api_20261024')
for dimension in ['none','two']:
 for mode in ['baseline','lean','full']:
  print(f'Starting {dimension}/{mode}',flush=True)
  with (root/f'{dimension}-{mode}-run.txt').open('w') as log:
   result=subprocess.run(['/usr/bin/time','-v','Rscript',str(root/'large-fit.R'),'/tmp/dbreg-two-fe-prediction-perf.duckdb',dimension,mode],stdout=log,stderr=subprocess.STDOUT)
  if result.returncode:raise SystemExit(f'Failed {dimension}/{mode}: {result.returncode}')
  print(f'Completed {dimension}/{mode}',flush=True)
