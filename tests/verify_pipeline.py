"""Run finite saved-result checks; never run the full permutation experiment."""
from pathlib import Path
import csv, subprocess, tempfile, sys
root = Path(__file__).resolve().parents[1]
def run(*args): subprocess.run(args, cwd=root, check=True)
run('Rscript', str(root/'tests/check_v6.R'), str(root))
run('Rscript', str(root/'analysis/repeated/R/preflight_production.R'))
with tempfile.TemporaryDirectory(prefix='mlsas-v6-') as tmp:
    tmp = Path(tmp)
    replay = tmp/'replay'
    run('Rscript', str(root/'scripts/replay_saved_f1.R'), str(root), str(replay))
    run(sys.executable, str(root/'scripts/verify_saved.py'), str(replay))
    with (root/'reference/v6/all_replication_metrics.csv').open() as f:
        rows = list(csv.DictReader(f))
    for dgp in range(2,5):
        stage = tmp/f'dgp{dgp}'; stage.mkdir()
        (stage/'STAGE_VALIDATED.txt').write_text('SAVED-METRIC TEST FIXTURE ONLY\n')
        with (stage/f'dgp{dgp}_stage_metrics.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader()
            writer.writerows(row for row in rows if row['dgp']==str(dgp))
    script = str(root/'analysis/repeated/R/summarize_production.R')
    run('Rscript', script, f'--run_dir={tmp}')
    run(sys.executable, str(root/'scripts/verify_saved.py'), str(tmp/'summary'))
    result=subprocess.run(['Rscript',script,f'--run_dir={tmp}','--correct_legacy_local_geary_auc=true'],capture_output=True,text=True)
    assert result.returncode != 0 and 'forbidden' in result.stderr
print('PASS: saved replay, integrated summary, metric edge case, DGP generator, AUC guard')
