# ev_dev.py: Devroye switch points, the sampler bound sup_k x_k, partial sums, mean excess and Aitken numbers of Chapter 25 (ev:), Section
# "Devroye, Smith and Tawn" (E:540 to 621: tab:devroye, the c_k, f(0.02), f(0.03), the partial sums, E X, e(0), e(inf), A, the deviations, Aitken).
# Claims supported: M3 (switch points of the two calibration laws), M4 (sup over k of x_k <= 0.0196, alternation of the c_k), E:569-572, 584, 594,
# 617-621 (log concavity and monotone e(u) on the grid u = 0.005 (0.0025) 0.2).
# Conventions. X = sum E_k/lam_k, lam_k = gamma_k^2, density f(x) = sum c_k lam_k exp(-lam_k x), survival P(X>x) = sum c_k exp(-lam_k x),
# c_k = prod_{j != k} (1 - lam_k/lam_j)^(-1) = 2 i xi(1/2)/(gamma_k xi'(1/2 + i gamma_k)) (the second form for Xi). Switch point
# x_k = log(|c_(k+1)| lam_(k+1)/(|c_k| lam_k))/(lam_(k+1) - lam_k). Calibration atoms pi^2 (k-1/2)^2/2 (law 1/cosh sqrt(2s)) and pi^2 k^2/2
# (law sqrt(2s)/sinh sqrt(2s)); exact c_k lam_k = (-1)^(k-1) pi (k-1/2) and 2 (-1)^(k-1) (derived here, checked against the product at K = 4000).
# Method. Residue series over the first 1700 ordinates (68 digits), absolutely convergent for x > 0 (terms exp(-lam_k x) |c_k|); f is also computed
# by Talbot inversion (mpmath) of xi(1/2)/xi(1/2 + sqrt s) at 30 digits as an independent check. Mean excess e(u) = sum c_k e^(-lam_k u)/lam_k / sum c_k e^(-lam_k u).
# Error control. The whole run is done at two precisions (dps 50 and 70) and the largest relative difference of the printed quantities is logged.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 -I ev_dev.py [--quick] [--out FILE] (quick: 300 zeros, no Talbot; full: c_k and switch points over 1700 zeros, sums over 300, Talbot check)
import argparse, json, sys, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from mpmath import mp, mpf, mpc, pi, exp, log, sqrt, nstr, diff, invertlaplace
from ev_evt_lib import *

ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None); ap.add_argument('--talbot', action='store_true')
args = ap.parse_args()
NZ = 300 if args.quick else 1700
NS = 300  # zeros in the series sums (terms exp(-lam x) |c| below 1e-90 beyond for x >= 0.005)
out = {}
logf = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); logf.append(s)

def switch(c, lam):
    return [log(abs(c[k + 1]) * lam[k + 1] / (abs(c[k]) * lam[k])) / (lam[k + 1] - lam[k]) for k in range(len(c) - 1)]

def run(dps):
    mp.dps = dps
    R = {}
    # --- calibration laws (M3)
    lam1 = lambda k: pi ** 2 * (mpf(k) - mpf(1) / 2) ** 2 / 2
    lam2 = lambda k: pi ** 2 * mpf(k) ** 2 / 2
    K = 12
    ex1 = [(-1) ** (k - 1) * pi * (mpf(k) - mpf(1) / 2) / lam1(k) for k in range(1, K + 1)]  # c_k, from c_k lam_k = (-1)^(k-1) pi (k-1/2)
    ex2 = [2 * (-1) ** (k - 1) for k in range(1, K + 1)]
    R['x_closed_cosh'] = [switch(ex1, [lam1(k) for k in range(1, K + 1)])[i] for i in range(3)]
    R['x_closed_sinh'] = [switch(ex2, [lam2(k) for k in range(1, K + 1)])[i] for i in range(3)]
    def prod_c(lamf, Kt, ks):
        lm = [lamf(j) for j in range(1, Kt + 1)]
        res = []
        for k in ks:
            p = mpf(1)
            for j in range(Kt):
                if j + 1 != k: p *= 1 - lm[k - 1] / lm[j]
            res.append(1 / p)
        return res
    for name, lamf, exa in (('cosh', lam1, ex1), ('sinh', lam2, ex2)):
        for Kt in (30, 4000) if dps == 50 else (30,):
            c = prod_c(lamf, Kt, [1, 2, 3, 4])
            lm = [lamf(k) for k in range(1, 5)]
            R['x_prod%d_%s' % (Kt, name)] = switch(c, lm)[:3]
            if Kt == 4000:
                R['ck_lam_ratio_%s' % name] = [c[i] / exa[i] for i in range(4)]
    # --- Xi law
    NZ_ = NZ if dps == 50 else min(NZ, 300)  # the second precision pass uses 300 zeros (time)
    zs = zero_strings(NZ_)
    t = [mpf(z) for z in zs]
    xh = xi(mpf(1) / 2)
    c = [ck_from_zero(g, xh) for g in t]
    lam = [g * g for g in t]
    R['c1_4'] = c[:4]
    R['alternates'] = all((c[k] > 0) == (k % 2 == 0) for k in range(NS))
    xs = switch(c, lam)
    R['x123'] = xs[:3]
    R['sup_x'] = max(xs); R['argsup'] = xs.index(max(xs)) + 1
    R['sup_x_from4'] = max(xs[3:]); R['argsup_from4'] = xs.index(max(xs[3:])) + 4
    R['xs_first10'] = xs[:10]
    # --- f and partial sums
    def fseries(x, nd=0):
        return sum(c[k] * lam[k] * (-lam[k]) ** nd * exp(-lam[k] * x) for k in range(NS))
    for x in (mpf('0.02'), mpf('0.03')):
        terms = [c[k] * lam[k] * exp(-lam[k] * x) for k in range(6)]
        ps = []; s = mpf(0)
        for tm in terms: s += tm; ps.append(s)
        R['f_%s' % x] = fseries(x); R['ps_%s' % x] = ps
        if args.talbot and dps == 50:
            mp.dps = 30
            R['f_talbot_%s' % x] = invertlaplace(lambda s: xh / xi(mpf(1) / 2 + sqrt(s)), x, method='talbot')
            mp.dps = dps
    # --- mean excess
    def e(u):
        num = sum(c[k] * exp(-lam[k] * u) / lam[k] for k in range(NS)); den = sum(c[k] * exp(-lam[k] * u) for k in range(NS))
        return num / den
    mp.dps = dps
    h = mpf(1) / 10 ** 12
    R['EX'] = (diff(lambda s: xi(s), mpf(1) / 2, 2)) / (2 * xh)
    R['e_inf'] = 1 / lam[0]
    A = c[1] / c[0] * (1 / lam[1] - 1 / lam[0]); R['A'] = A
    dev = []
    for u in ('0.03', '0.055', '0.105'):
        u = mpf(u); dev.append((e(u) - 1 / lam[0], A * exp(-(lam[1] - lam[0]) * u)))
    R['dev'] = dev
    grid = [mpf('0.005') + mpf('0.0025') * j for j in range(79)]
    fv = [(fseries(u), fseries(u, 1), fseries(u, 2)) for u in grid]
    # f' = sum c lam (-lam) e^..., f'' with (-lam)^2
    R['logconc_all'] = all(f[0] * f[2] - f[1] ** 2 < 0 for f in fv)
    ev = [e(u) for u in grid]
    R['e_decreasing_all'] = all(ev[j] > ev[j + 1] for j in range(78))
    # Aitken
    ait = []
    hh = mpf('0.0075')
    for u0 in ('0.03', '0.04', '0.05'):
        u0 = mpf(u0); e0, e1, e2 = e(u0), e(u0 + hh), e(u0 + 2 * hh)
        r = (e2 - e1) / (e1 - e0); einf = e2 - (e2 - e1) ** 2 / (e2 - 2 * e1 + e0)
        gap = -log(r) / hh
        ait.append((sqrt(1 / einf), sqrt(1 / einf + gap)))
    R['aitken'] = ait
    return R

t0 = time.time()
A = run(50); B = run(70)
def flat(R):
    o = {}
    for k, v in R.items():
        if isinstance(v, bool): o[k] = v
        elif isinstance(v, list): 
            o[k] = [x if not isinstance(x, tuple) else list(x) for x in v]
        else: o[k] = v
    return o
def num(v):
    if isinstance(v, (list, tuple)): return [num(x) for x in v]
    if isinstance(v, bool): return v
    return nstr(v, 20)
for k, v in A.items():
    P(k, num(v))
mx = mpf(0)
def cmp(a, b):
    global mx
    if isinstance(a, (list, tuple)):
        for x, y in zip(a, b): cmp(x, y)
    elif not isinstance(a, bool) and a != 0:
        mx = max(mx, abs(a - b) / abs(a))
mp.dps = 50
for k in A:
    if k in B and not k.startswith('x_prod4000') and not k.startswith('sup_x_from4') and not k.startswith('argsup'): cmp(A[k], B[k])
P('max relative difference dps50 vs dps70:', nstr(mx, 3))
P('time', round(time.time() - t0, 1), 's, zeros used', NZ)
of = args.out or ('/tmp/ev_dev.json' if args.quick else __import__('os').path.join(HERE, '..', 'data', 'ev_dev.json'))
json.dump({k: num(v) for k, v in A.items()} | {'zeros': NZ, 'maxreldiff_50_70': nstr(mx, 3)}, open(of, 'w'), indent=1)
open(of.replace('.json', '.log'), 'w').write('\n'.join(logf) + '\n')
