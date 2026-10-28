# code/hx_scan.py (October 2026, from misc/reviews/nick-b4/reconcile-P.md, task N2, remark cp:rem:minflow and the proposition cp:prop:count): the zeros of xi_{0.2} on the critical line, 0 <= theta <= 4200, written from the mathematics (the collaborator's scripts, notes and pickles were not read, run or copied).
# Definitions. Z_t(theta) = e^{i(vartheta(theta) + (t/4) Im l(1/2+i theta)^2)} zeta_t(1/2+i theta), with zeta_t from the Gauss Hermite evaluator of Lemma cp:lem:shift (hx_shift.py: K nodes, Arb balls with python-flint, midpoints used); Z_t is real and its sign changes are the zeros of xi_t(1/2+i theta) of odd multiplicity. Main term of the count: M_t(T) = 1 + vartheta(T)/pi + (t/(4 pi)) Im l(1/2+iT)^2 (Proposition cp:prop:count (i)); N_t(T) = M_t(T) + S_t(T).
# Method. Grid theta_j = 0.05 j, 0 <= j <= 84000, K = 40 at 25 digits; at every sign change (bracket) of Re Z both endpoints are recomputed with K = 60 and the sign must agree (otherwise the bracket is flagged); brackets are bisected (K = 40) to width 1e-6 (at most 16 halvings). The located zeros are the midpoints. A pair of zeros closer than the grid step with no sign change on the grid is invisible: the counts below are counts of located sign changes, a lower bound for N_t(T) when all zeros are on the line. The flags: sign mismatches between K = 40 and K = 60, the maximum |Im Z|/|Re Z| over the grid (rounding of the phase), the number of brackets.
# Report: the first zero, the largest gap between consecutive located zeros on [theta_1, 4200] with its location, the counts below T = 500, 1000, ..., 4000 and their differences from M_t(T) (N2d is the largest absolute difference).
# Error control: the sign recheck with K = 60; a second run at 40 digits and K = 60 on a random sample of 400 grid points (relative difference of Re Z logged).
# Output: hx_scan.json and hx_scan.log (the zeros are in the json as 15 digit strings); python3 hx_scan.py [--quick] [--out FILE] [--tmax 4200] [--procs N]. Hopper: one node, 48 cores, about 10 minutes. Quick: tmax = 100.
import argparse
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hx_shift as hs

T = '0.2'
EV = {}


def ev(dps):
    if dps not in EV:
        EV[dps] = hs.Ev(dps)
    return EV[dps]


def zval(theta, K, dps=25):
    e = ev(dps)
    z = e.to_c(e.Z(repr(float(theta)), T, K))
    return z


def chunk(args):
    lo, hi, step, K = args
    out = []
    for j in range(lo, hi):
        z = zval(j * step, K)
        out.append((j, z.real, abs(z.imag)))
    return out


def recheck(args):
    ths, K = args[:2]
    dps = args[2] if len(args) > 2 else 25
    return [zval(th, K, dps) for th in ths]


def bisect(args):
    a, b, fa = args
    for _ in range(16):
        m = 0.5 * (a + b)
        fm = zval(m, 40).real
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
        if b - a < 1e-6:
            break
    return 0.5 * (a + b), a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='hx_scan.json')
    ap.add_argument('--tmax', type=float, default=4200.0)
    ap.add_argument('--procs', type=int, default=0)
    a = ap.parse_args()
    hs.LOGF = os.path.splitext(a.out)[0] + '.log'
    open(hs.LOGF, 'w').close()
    say = hs.say
    tmax = 100.0 if a.quick else a.tmax
    step = 0.05
    n = int(round(tmax / step))
    procs = a.procs or (os.cpu_count() or 4)
    say('hx_scan.py: t = 0.2, grid step', step, 'up to', tmax, '(%d points), K = 40, 25 digits, %d processes, backend %s' % (n + 1, procs, 'flint' if hs.HAVE_FLINT else 'mpmath'))
    cs = 200
    chunks = [(lo, min(lo + cs, n + 1), step, 40) for lo in range(0, n + 1, cs)]
    with Pool(procs) as pool:
        rows = []
        for r in pool.imap(chunk, chunks):
            rows.extend(r)
        say('grid done')
        rows.sort()
        vals = [r[1] for r in rows]
        maxim = max(r[2] / abs(r[1]) for r in rows if r[1] != 0)
        say('max |Im Z|/|Re Z| over the grid: %.3e' % maxim)
        br = [j for j in range(n) if (vals[j] > 0) != (vals[j + 1] > 0)]
        say('brackets (sign changes of Re Z on the grid):', len(br))
        # K = 60 recheck of both endpoints
        pts = sorted(set([j * step for j in br] + [(j + 1) * step for j in br]))
        parts = [(pts[i:i + 100], 60) for i in range(0, len(pts), 100)]
        z60 = {}
        for ps, r in zip([p[0] for p in parts], pool.imap(recheck, parts)):
            for th, z in zip(ps, r):
                z60[th] = z
        mism = []
        for j in br:
            for jj in (j, j + 1):
                if (z60[jj * step].real > 0) != (vals[jj] > 0):
                    mism.append(jj * step)
        say('K = 60 recheck: endpoints %d, sign mismatches with K = 40: %d' % (len(pts), len(mism)), mism[:10])
        mx = max(abs(z60[jj * step].real - vals[jj]) / abs(z60[jj * step].real) for j in br for jj in (j, j + 1))
        say('max relative difference of Re Z, K = 40 vs 60, at the bracket endpoints: %.3e' % mx)
        zs = list(pool.imap(bisect, [(j * step, (j + 1) * step, vals[j]) for j in br], chunksize=20))
        zeros = [z[0] for z in zs]
        say('zeros located:', len(zeros), 'first', zeros[0] if zeros else None)
        # second precision on a random sample
        random.seed(1)
        samp = sorted(random.sample(range(n + 1), 40 if a.quick else 400))
        sp = [([j * step for j in samp[i:i + 20]], 60, 40) for i in range(0, len(samp), 20)]
        d = 0.0
        for ps, r in zip(sp, pool.imap(recheck, sp)):
            for th, z in zip(ps[0], r):
                jj = int(round(th / step))
                d = max(d, abs(z.real - vals[jj]) / max(abs(z.real), 1e-300))
        say('sample of %d grid points, K = 60 at 40 digits against K = 40 at 25 digits: max relative difference of Re Z %.3e' % (len(samp), d))
    res = {'script': 'hx_scan.py', 'tmax': tmax, 'step': step, 'K': 40, 'digits': 25, 'quick': a.quick}
    res['N2a'] = tmax
    res['first_zero'] = zeros[0]
    res['n_zeros'] = len(zeros)
    res['brackets'] = len(br)
    res['k60_mismatches'] = mism
    res['max_imag_over_real'] = maxim
    gaps = [(zeros[i + 1] - zeros[i], zeros[i]) for i in range(len(zeros) - 1)]
    g, loc = max(gaps)
    res['largest_gap'] = g
    res['largest_gap_after_zero'] = loc
    res['largest_gap_between'] = [loc, loc + g]
    res['smallest_gap'] = min(gaps)
    res['gaps_over_9_2306'] = [(x, y) for x, y in gaps if x > 9.2306]
    say('N2b largest gap %.6f between zeros %.6f and %.6f' % (g, loc, loc + g))
    say('smallest gap %.6f at %.6f' % min(gaps))
    # counts and main term
    mp.mp.dps = 30
    e = hs.Ev(30)
    counts = {}
    mx = 0.0
    Ts = [500 * k for k in range(1, 9) if 500 * k <= tmax] if not a.quick else [50, 100]
    for TT in Ts:
        cnt = sum(1 for z in zeros if z <= TT)
        s = e.cx('0.5', TT)
        mt = 1 + e.to_mp(e.theta_rs(TT)) / mp.pi + mp.mpf(T) / (4 * mp.pi) * (e.to_mp(e.ell(s)) ** 2).imag
        diff = cnt - mt
        mx = max(mx, abs(float(diff)))
        counts[str(TT)] = {'count': cnt, 'main_term': mp.nstr(mt, 10), 'count_minus_main': mp.nstr(diff, 6)}
        say('T = %d: sign changes %d, main term %s, difference %s' % (TT, cnt, mp.nstr(mt, 10), mp.nstr(diff, 5)))
    res['counts'] = counts
    res['N2d_max_abs_difference'] = mx
    res['N2d_rounded_up'] = math.ceil(mx * 100) / 100
    say('N2d largest |difference| %.4f, rounded up %.2f' % (mx, res['N2d_rounded_up']))
    res['zeros'] = ['%.8f' % z for z in zeros]
    res['seconds'] = round(time.time() - hs.T0, 1)
    json.dump(res, open(a.out, 'w'))
    say('wrote', a.out)


if __name__ == '__main__':
    main()
