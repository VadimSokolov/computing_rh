# code/hx_edge.py (October 2026, from misc/reviews/nick-b4/reconcile-P.md, task N4, Proposition cp:prop:count (iii) and Problem cp:prob:rvm (b)): the right edge Re zeta_t(2 + i y), written from the mathematics (the collaborator's scripts, notes and pickles were not read, run or copied).
# Definition. zeta_t = xi_t/(g e^{t l^2/4}) is evaluated by the Gauss Hermite rule of Lemma cp:lem:shift (hx_shift.py; for s = 2 + i y with y small the formula Psi = xi(s+h)/(g(s) e^{h l}) is still valid, it is an identity of entire functions divided by g(s) != 0); at t = 0 it is zeta itself.
# Computation. The minimum of Re zeta_t(2 + i y) over 0 <= y <= 200 on the grid y = 0.01 k, for t = 0, 0.05, 0.1, 0.15, 0.2 (K = 40, and K = 60 at the same points; 30 digits), the grid minimiser refined by golden section search on [y - 0.01, y + 0.01] (K = 60), the value at y = 0, and Re zeta_t(2 + i y) at y = 10^3, 10^4, 10^6 (K = 40 and 60, 30 and 40 digits). The edge hypothesis of the proposition is Re zeta_t(2 + i y) > 0 on 0 <= y <= T; the minimum is [N4a] at y = [N4b] for t = 0.2.
# Error control. K = 40 against K = 60 over the whole grid (the largest difference of Re zeta_t is logged), t = 0 against mpmath zeta, the two working precisions at the three large heights.
# Output: hx_edge.json and hx_edge.log; python3 hx_edge.py [--quick] [--out FILE] [--procs N]. Hopper: a few minutes on 48 cores. Quick: y <= 20, step 0.1.
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hx_shift as hs

EV = {}


def ev(dps=30):
    if dps not in EV:
        EV[dps] = hs.Ev(dps)
    return EV[dps]


def rz(y, t, K, dps=30):
    e = ev(dps)
    return e.to_c(e.zt(e.cx('2', repr(float(y))), t, K))


def chunk(args):
    ys, t, Ks = args
    return [[rz(y, t, K).real for K in Ks] for y in ys]


def golden(t, lo, hi, K):
    g = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = rz(c, t, K).real, rz(d, t, K).real
    for _ in range(40):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = rz(c, t, K).real
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = rz(d, t, K).real
    y = (a + b) / 2
    return y, rz(y, t, K).real


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='hx_edge.json')
    ap.add_argument('--procs', type=int, default=0)
    a = ap.parse_args()
    hs.LOGF = os.path.splitext(a.out)[0] + '.log'
    open(hs.LOGF, 'w').close()
    say = hs.say
    procs = a.procs or (os.cpu_count() or 4)
    step, ymax = (0.1, 20.0) if a.quick else (0.01, 200.0)
    n = int(round(ymax / step))
    ys = [step * k for k in range(n + 1)]
    Ks = (40,) if a.quick else (40, 60)
    res = {'script': 'hx_edge.py', 'quick': a.quick, 'step': step, 'ymax': ymax, 'grid': {}}
    say('hx_edge.py, backend', 'flint' if hs.HAVE_FLINT else 'mpmath', 'step', step, 'y <=', ymax)
    with Pool(procs) as pool:
        for t in ['0', '0.05', '0.1', '0.15', '0.2']:
            out = []
            KK = Ks if t != '0' else (40,)
            for r in pool.imap(chunk, [(ys[i:i + 100], t, KK) for i in range(0, len(ys), 100)]):
                out.extend(r)
            kp = len(out[0]) - 1
            vals = [o[kp] for o in out]
            j = min(range(len(vals)), key=lambda i: vals[i])
            dK = max(abs(o[0] - o[kp]) for o in out)
            yr, vr = golden(t, max(0.0, ys[j] - step), ys[j] + step, KK[-1])
            say('t = %s: Re zeta_t(2 + i y) at y = 0: %.10f; grid minimum %.10f at y = %.2f; refined %.10f at y = %.6f; largest K=40 vs 60 difference over the grid %.2e' % (t, vals[0], vals[j], ys[j], vr, yr, dK))
            res['grid'][t] = {'value_at_0': vals[0], 'grid_min': vals[j], 'grid_argmin': ys[j], 'refined_min': vr, 'refined_argmin': yr, 'max_K40_vs_K60': dK}
            if t == '0':
                mp.mp.dps = 20
                d = max(abs(vals[i] - float(mp.zeta(mp.mpc(2, ys[i])).real)) for i in range(0, len(ys), 37))
                say('   t = 0 against mpmath zeta on a sample: largest difference %.2e' % d)
                res['grid'][t]['vs_mpmath_sample'] = d
        res['N4a'] = res['grid']['0.2']['refined_min']
        res['N4b'] = res['grid']['0.2']['refined_argmin']
        res['large'] = {}
        for y in ([100.0] if a.quick else [1e3, 1e4, 1e6]):
            row = {}
            for dps in (30, 40):
                for K in Ks:
                    row['dps%d_K%d' % (dps, K)] = rz(y, '0.2', K, dps).real
            say('Re zeta_0.2(2 + i %g): %s' % (y, {k: '%.12f' % v for k, v in row.items()}))
            res['large']['%g' % y] = row
    res['seconds'] = round(time.time() - hs.T0, 1)
    json.dump(res, open(a.out, 'w'), indent=1)
    say('wrote', a.out)


if __name__ == '__main__':
    main()
