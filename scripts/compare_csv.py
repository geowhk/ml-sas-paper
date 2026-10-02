"""Strict schema/row-order CSV comparison with absolute numeric tolerance 1e-12."""
import csv,math,sys
from pathlib import Path
if len(sys.argv)!=3:raise SystemExit('Usage: compare_csv.py EXPECTED_CSV ACTUAL_CSV')
with Path(sys.argv[1]).open() as f:a=list(csv.DictReader(f))
with Path(sys.argv[2]).open() as f:b=list(csv.DictReader(f))
assert len(a)==len(b),(len(a),len(b))
for i,(x,y) in enumerate(zip(a,b)):
 assert x.keys()==y.keys(),(x.keys(),y.keys())
 for k in x:
  if x[k]!=y[k]:
   try:ok=math.isclose(float(x[k]),float(y[k]),rel_tol=0,abs_tol=1e-12)
   except ValueError:ok=False
   assert ok,(i,k,x[k],y[k])
print('PASS',Path(sys.argv[1]).name,len(a),'rows')
