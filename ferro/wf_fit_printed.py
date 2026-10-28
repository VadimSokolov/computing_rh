# Note W, computation N1 (side check): the distance TV(q, p) of row 2 of Table fm:tab:bessel at the parameters PRINTED by the
# note of Nicholas Polson (a = 5.484, c = 4.126, nu = 4.292, i.e. book units a = 5.484, l = c/2 = 2.063, lam = nu/2 = 2.146)
# and at the locally reoptimised member quoted by the review (5.483990, 2.062861, 2.146192), by the objective of wf_fit.py
# (Simpson, step 1e-4), against the optimum found by wf_fit.py (wf_fit.json). Double precision, not certified.
# Usage: python3 -I wf_fit_printed.py > wf_fit_printed.log
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = [sys.argv[0]]
import wf_fit as F
OB = F.Obj(1e-4)
opt = json.load(open('wf_fit.json'))['row2']
for name, (a, l, lam) in (('printed', (5.484, 2.063, 2.146)), ('review reoptimised', (5.483990, 2.062861, 2.146192)),
                          ('wf_fit.py decimals', (opt['a'], opt['l'], opt['lams'][0]))):
    print('%-20s a=%s l=%s lam=%s  TV = %.6g' % (name, a, l, lam, OB.tv(float(a), float(l), [float(lam)], 0.0)))
print('wf_fit.py optimum (double precision, unrounded) = %.6g' % opt['tv_opt'])
