# ht_hcm2exp.py: the HCM ladder of the two exponential density in Nick Polson's residue note (R:327, Theorem nohcm (i) and the paragraph after it):
# for f(x) = e^(-x) - e^(-2x) (the density of the sum of two exponentials of rates 1 and 2 up to the constant 2) and h_u(w) = f(u v) f(u/v),
# w = v + 1/v, the first order k with (-1)^k h_u^(k)(w)/h_u(w) < 0 at w = 2.1 for u = 3, 4, 5, 6, 7, 8, 10, and its value at u = 5.
# Convention. v = (w + sqrt(w^2 - 4))/2 > 1; h_u is analytic near w = 2.1 (singular only at w = 2); derivatives are taken from the Taylor series of h_u at w = 2.1
# computed by mpmath at dps 300 and again at dps 450 (the number printed is the one of the lower precision; the relative difference of
# the printed value is logged). The chapter does not say at which w the cases u = 6, 7, 8 were tested; w = 2.1 is taken for all.
# Options: --quick (dps 120/160, orders to 20), --out FILE. Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, argparse
from mpmath import mp, mpf, exp, sqrt, taylor, nstr
ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/ht_hcm2exp.json')
a = ap.parse_args()
LOG = open(a.out.replace('.json', '.log'), 'w')
def lg(*x):
    s = ' '.join(str(t) for t in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
K = 20 if a.quick else 40
def run(dps, u, w0):
    mp.dps = dps
    u = mpf(u); w0 = mpf(w0)
    def h(w):
        v = (w + sqrt(w * w - 4)) / 2
        f = lambda x: exp(-x) - exp(-2 * x)
        return f(u * v) * f(u / v)
    c = taylor(h, w0, K)
    res = []
    for k in range(1, K + 1):
        res.append((-1) ** k * c[k] * mp.factorial(k) / c[0])
    return res
R = {}
for u in (3, 4, 5, 6, 7, 8, 10):
    r1 = run(120 if a.quick else 300, u, '2.1'); r2 = run(160 if a.quick else 450, u, '2.1')
    neg = [k + 1 for k, x in enumerate(r1) if x < 0]
    first = neg[0] if neg else None
    mp.dps = 30
    d = max(abs(r1[k] / r2[k] - 1) for k in range(K)) if True else 0
    lg('u=%d w=2.1: first negative order %s (value %s); negative orders %s; max relative difference of the 40 values between precisions %s' % (u, first, nstr(r1[first - 1], 4) if first else None, neg[:6], nstr(d, 3)))
    R['u%d' % u] = [first, nstr(r1[first - 1], 4) if first else None]
import json; json.dump(R, open(a.out, 'w'))
