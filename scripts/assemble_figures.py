"""Map fresh analysis outputs to the historical figure input names, without corrections."""
from pathlib import Path
import sys,csv,shutil
root=Path(__file__).resolve().parents[1]
if len(sys.argv)!=3:raise SystemExit('Usage: assemble_figures.py FULL_RUN_DIR NEW_INPUT_DIR')
run=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
for expected in sorted((root/'reference/figure_inputs').rglob('*')):
 if not expected.is_file():continue
 rel=expected.relative_to(root/'reference/figure_inputs');target=out/rel;target.parent.mkdir(parents=True,exist_ok=True)
 if rel.parts[0]=='dgp':source=run/'local-evidence'/expected.name
 elif 'sim_dgp' in rel.parts:source=root/'02_revision-simulation/sim_dgp'/expected.name
 elif rel.parts[0]=='03_hyperparameter-sensitivity':
  stage='sensitivity' if 'mixed_hyperparameter_path' in rel.parts else 'path-terminal'
  source=run/stage/'tables/sensitivity_results_raw.csv'
  with source.open() as f:reader=csv.DictReader(f);rows=list(reader);fields=reader.fieldnames
  # Fixed DGP 4.3 has 60 target cells among 400, as asserted by generator checks.
  for row in rows:row['global_mlsas_observed']=(float(row['true_local_mean'])*60+float(row['false_local_mean'])*340)/400
  with target.open('w') as f:w=csv.DictWriter(f,fieldnames=fields+['global_mlsas_observed']);w.writeheader();w.writerows(rows)
  continue
 elif rel.parts[0]=='04_current-statistics-comparison':source=run/'comparison/tables'/expected.name
 else:
  name=expected.name.replace('_results_corrected.csv','_results.csv')
  found=list((run/'representative').rglob(name))
  if len(found)!=1:raise RuntimeError(f'Expected exactly one fresh {name}; found {found}')
  source=found[0]
 shutil.copy2(source,target)
print('Prepared fresh figure inputs; no reference result CSV was substituted')
