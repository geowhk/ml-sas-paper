"""Check every numeric/label figure input against v6, excluding legacy correction annotations."""
from pathlib import Path
import csv,math,sys
root=Path(__file__).resolve().parents[1];actual=Path(sys.argv[1])
for expected in sorted((root/'reference/figure_inputs').rglob('*.csv')):
 with expected.open() as f:a=list(csv.DictReader(f))
 with (actual/expected.relative_to(root/'reference/figure_inputs')).open() as f:b=list(csv.DictReader(f))
 assert len(a)==len(b),expected
 ignore={'geary_p_source_less','geary_p_adj_source_less','geary_flag_source_less','local_geary_alternative','local_geary_interpretation'} if expected.name.endswith('lisa_results_corrected.csv') else set()
 for x,y in zip(a,b):
  assert set(x)-ignore==set(y)-ignore,(expected,x.keys(),y.keys())
  for key in set(x)-ignore:
   if x[key]!=y[key]:
    try:ok=math.isclose(float(x[key]),float(y[key]),rel_tol=0,abs_tol=1e-12)
    except ValueError:ok=False
    assert ok,(expected.name,key,x[key],y[key])
print('PASS: all figure CSV inputs match v6 numeric values and labels')
