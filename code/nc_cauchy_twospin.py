# Exact two spin model of Proposition cm:prop:twospin: X = Y1 + Y2, P(Y = +-1) = q < 1/4.
# cosh a = (1-2q)/(2q); the cosine tilt is ferromagnetic exactly for w >= w*, cosh 2w* = cosh^2 a.
import json, os
from mpmath import mp, mpf, acosh, asinh, sinh, sqrt
mp.dps = 30
rows = []
for q in ['0.2', '0.24', '0.249', '0.2499', '0.24999']:
    q = mpf(q); a = acosh((1 - 2*q)/(2*q)); w = acosh(mp.cosh(a)**2)/2
    assert abs(w - asinh(sinh(a)/sqrt(2))) < mpf(10)**-25 and a/sqrt(2) < w < a
    rows.append({'q': float(q), 'a': float(a), 'wstar': float(w), 'ratio': float(w/a)})
    print(f"q={float(q)}: a={float(a):.6f} w*={float(w):.6f} w*/a={float(w/a):.6f}")
print("1/sqrt2 =", float(1/sqrt(2)))
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'nc_cauchy_twospin.json')
json.dump(rows, open(out, 'w'), indent=1)
