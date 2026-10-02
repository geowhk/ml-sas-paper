"""End-to-end reproduction. Full mode includes the slow original experiments."""
from pathlib import Path
import argparse,subprocess,sys,time,json,hashlib,os,shutil
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--mode',choices=['saved','full'],required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--workers',type=int,default=1)
p.add_argument('--parallelism',type=int,default=1)
p.add_argument('--resume',action='store_true')
a=p.parse_args();out=a.out.resolve()
if a.workers<1 or a.parallelism<1:p.error('Worker counts must be positive')
if out.exists() and not a.resume:p.error('Output exists; use a new path or --resume')
if not out.exists() and a.resume:p.error('Cannot resume a missing output')
# Bind resumption to all code/config/reference bytes, not just a completion filename.
source_dirs=['analysis','02_revision-simulation','03_hyperparameter-sensitivity','04_current-statistics-comparison','05_empirical-study','scripts','tests','reference']
files=sorted(f for directory in source_dirs for f in (root/directory).rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix in ('.py','.R','.sh','.yml','.json','.csv','.rds','.gz') and 'runs' not in f.relative_to(root).parts and out not in f.parents)+[root/'renv.lock']
fingerprint=hashlib.sha256(''.join(str(f.relative_to(root))+hashlib.sha256(f.read_bytes()).hexdigest() for f in files).encode()).hexdigest()
statefile=out/'reproduction.json'
if a.resume:
 state=json.loads(statefile.read_text())
 if state['code_sha256']!=fingerprint or state['mode']!=a.mode:p.error('Code/input/mode changed; use a new run directory')
else:
 out.mkdir(parents=True);state={'mode':a.mode,'code_sha256':fingerprint,'steps':{},'status':'RUNNING'}
state['status']='RUNNING'
def save():statefile.write_text(json.dumps(state,indent=2)+'\n')
save()
def step(name,cmd):
 if state['steps'].get(name,{}).get('exit_code')==0:
  for relative,expected in state['steps'][name].get('output_sha256',{}).items():
   f=out/relative
   if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=expected:
    raise SystemExit(f'Completed output changed or missing: {f}; use a new output directory')
  return
 print(name,flush=True)
 target_name={'export-local':'local-evidence','assemble':'figure-inputs','dgp-figures':'dgp-figures'}.get(name,name)
 target=out/target_name
 if a.resume and name not in ['repeated','sensitivity','path-terminal'] and target.is_dir():
  failed=out/'incomplete';failed.mkdir(exist_ok=True)
  shutil.move(str(target),str(failed/(target_name+'-'+str(time.time_ns()))))
 state['steps'][name]={'command':[str(v) for v in cmd],'exit_code':None};save()
 start=time.monotonic()
 with (out/(name+'.log')).open('a') as log:r=subprocess.run([str(v) for v in cmd],cwd=root,stdout=log,stderr=subprocess.STDOUT)
 state['steps'][name]={'command':[str(v) for v in cmd],'elapsed_seconds':time.monotonic()-start,'exit_code':r.returncode,'output_sha256':{str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(target.rglob('*')) if f.is_file()} if target.is_dir() and r.returncode==0 else {}};save()
 if r.returncode:
  state['status']='FAILED';save();raise SystemExit(f'{name} failed; see {out/(name+".log")}')
py=sys.executable
step('preflight',['Rscript',root/'analysis/repeated/R/preflight_production.R'])
step('empirical-inputs',['Rscript',root/'scripts/prepare_empirical.R',root,out/'empirical-inputs'])
if a.mode=='full':
 for stage in ['representative','sensitivity','path-terminal','comparison','empirical','repeated']:
  cmd=[py,root/'scripts/run_analysis.py',stage,'--out',out/stage,'--workers',a.workers,'--parallelism',a.parallelism]
  if stage in ['repeated','sensitivity','path-terminal'] and (out/stage).exists():cmd+=['--resume']
  step(stage,cmd)
 step('f1-compare',[py,root/'scripts/verify_saved.py',out/'repeated/summary'])
 step('null-summary-check',[py,root/'scripts/compare_csv.py',root/'reference/legacy_summary/dgp2_null_operating_characteristics.csv',out/'repeated/summary/dgp2_null_operating_characteristics.csv'])
 for name in ['global_mlsas','classification_summary','continuous_vs_mixed_local']:
  step('empirical-check-'+name,[py,root/'scripts/compare_csv.py',root/f'reference/empirical/tables/{name}.csv',out/f'empirical/conditional/tables/{name}.csv'])
 step('export-local',['Rscript',root/'scripts/export_local_evidence.R',root,out/'repeated',out/'local-evidence'])
 evidence=out/'local-evidence/repeated_local_evidence.csv'
 step('assemble',[py,root/'scripts/assemble_figures.py',out,out/'figure-inputs'])
 step('figure-input-check',[py,root/'scripts/compare_figure_inputs.py',out/'figure-inputs'])
 inputs=out/'figure-inputs'
 empirical=out/'empirical/conditional'
else:
 step('saved-f1',['Rscript',root/'scripts/replay_saved_f1.R',root,out/'saved-f1'])
 step('f1-compare',[py,root/'scripts/verify_saved.py',out/'saved-f1'])
 evidence=root/'reference/repeated_local_evidence.csv.gz'
 inputs=root/'reference/figure_inputs'
 empirical=root/'reference/empirical'
step('appendix',[py,root/'scripts/analyze_appendix.py',evidence,out/'appendix'])
for name in ['dgp3_pvalue_tail_rates','category_diagnostics_summary','category_pair_auc_summary']:
 step('check-'+name,[py,root/'scripts/compare_csv.py',root/f'reference/appendix/{name}.csv',out/f'appendix/{name}.csv'])
step('section4',['Rscript',root/'scripts/render_section4.R',root,out/'section4',inputs])
step('dgp-figures',[py,root/'scripts/render_dgp_figures.py',out/'dgp-figures',inputs])
step('figure8',[py,root/'scripts/render_empirical_inputs.py',out/'empirical-inputs/empirical_inputs.geojson',out/'figure8'])
step('figure9',['Rscript',root/'scripts/render_empirical_results.R',empirical,out/'empirical-inputs/analysis_inputs.rds',out/'figure9'])
state['status']='PASS';save();print(f'PASS: {a.mode} reproduction. Timing and logs: {out}')
