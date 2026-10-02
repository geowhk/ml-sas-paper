"""Check repeated runner's fail-closed resume and overwrite guards without fitting models."""
from pathlib import Path
import subprocess,tempfile
root=Path(__file__).resolve().parents[1];script=root/'analysis/repeated/run_production_batch.sh'
with tempfile.TemporaryDirectory() as tmp:
 out=Path(tmp)/'run';out.mkdir()
 for args in [[out,'1'],[out,'1','--resume'],[out,'1','invalid']]:
  result=subprocess.run(['bash',str(script),*[str(v) for v in args]],cwd=root,capture_output=True)
  assert result.returncode!=0
 assert not list(out.iterdir()), 'Rejected calls must not write to the existing directory'
print('PASS: overwrite, incomplete resume metadata and invalid flag rejected before analysis')
