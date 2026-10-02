"""Explicit output locations for the paper's full analyses; no result corrections."""
import argparse, subprocess, os, time, json, datetime
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('analysis',choices=['repeated','representative','sensitivity','path-terminal','comparison','lovo','empirical'])
p.add_argument('--out',type=Path,required=True)
p.add_argument('--resume',action='store_true',help='Resume repeated or sensitivity stages after checking code/configuration')
p.add_argument('--workers',type=int,default=1,help='Workers inside non-repeated analyses; does not change analysis dimensions')
p.add_argument('--parallelism',type=int,default=1,help='Repeated-DGP concurrent processes only')
a=p.parse_args();out=a.out.resolve()
if out.exists() and not (a.resume and a.analysis in ['repeated','sensitivity','path-terminal']):p.error('Output must not exist')
if a.resume and a.analysis not in ['repeated','sensitivity','path-terminal']:p.error('--resume is supported for repeated and sensitivity analyses')
if a.workers<1:p.error('Workers must be positive')
os.environ['REDUCED_WORKERS']=str(a.workers)
os.environ['MLSAS_WORKERS']=str(a.workers)
os.environ['MLSAS_RESUME']='1' if a.resume else '0'
if a.parallelism<1:p.error('Parallelism must be positive')
if out==root or root/'reference' in out.parents:p.error('Choose an output outside reference inputs')
def run(*cmd):
 start=time.monotonic();result=subprocess.run([str(x) for x in cmd],cwd=root)
 if out.exists():
  with (out/'timing.jsonl').open('a') as f:f.write(json.dumps({'command':[str(x) for x in cmd],'elapsed_seconds':time.monotonic()-start,'exit_code':result.returncode,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})+'\n')
 result.check_returncode()
if a.analysis=='repeated':
 run('bash',root/'analysis/repeated/run_production_batch.sh',out,a.parallelism,*(['--resume'] if a.resume else []))
elif a.analysis=='representative':
 run('Rscript','02_revision-simulation/R/run_all_dgps.R','--config','02_revision-simulation/config/revision_v2.yml','--output_root',out,'--workers',a.workers)
elif a.analysis=='lovo':
 run('Rscript','02_revision-simulation/R/run_dgp43_lovo.R',root,out,root/'02_revision-simulation/config/dgp43_lovo_revision_v2.R')
elif a.analysis=='empirical':
 out.mkdir(parents=True)
 run('Rscript','05_empirical-study/R/run_reduced_urbanity_analysis.R',root,out/'base')
 run('Rscript','05_empirical-study/R/upgrade_reduced_inference_99999.R',root,out/'base',out/'global')
 run('Rscript','05_empirical-study/R/rerun_conditional_local_inference.R',root,out/'global',out/'conditional',out/'base')
else:
 module,config = {
  'sensitivity':('03_hyperparameter-sensitivity','mixed_hyperparameter_path'),
  'path-terminal':('03_hyperparameter-sensitivity','mixed_path_vs_terminal'),
  'comparison':('04_current-statistics-comparison','current_statistics_20x20')
 }[a.analysis]
 script='run_current_statistics_comparison.R' if a.analysis=='comparison' else 'run_mixed_sensitivity.R'
 run('Rscript',f'{module}/R/{script}',f'{module}/config/{config}.yml',out)
