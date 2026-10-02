"""Compare replayed CSVs with frozen v6 values (absolute tolerance 1e-12)."""
from pathlib import Path
import csv, math, sys
root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1])
cells = 0
for expected in sorted((root / 'reference/v6').glob('*.csv')):
    with expected.open() as f: a = list(csv.DictReader(f))
    with (out / expected.name).open() as f: b = list(csv.DictReader(f))
    assert len(a) == len(b), expected.name
    for x, y in zip(a, b):
        assert x.keys() == y.keys(), expected.name
        for key in x:
            v, w = x[key], y[key]
            cells += 1
            if v != w:
                assert math.isclose(float(v), float(w), rel_tol=0, abs_tol=1e-12), (expected.name,v,w)
print(f'PASS: five v6 CSVs, {cells} cells, absolute tolerance 1e-12')
