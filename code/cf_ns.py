# cf_ns.py: Table tab:ns of the source note L (Chapter 34, prefix cf, Section cf:sec:ns): the scores s(t) at the four heights of the table (and
# at 1004606.8, the plain Monte Carlo point of the text). Group CF-XI, task T8; supports R15 (the score) and N13. Written from the mathematics
# (the collaborator's scripts were not read, run or copied).
# Score (L:1571-1573): s(t) = u_{sigma0}(t)/rho(t), sigma0 = 0.6, u_theta(v) = pi^(-1) Re (xi'/xi)(theta + i v) (the density of the chapter's measure on the
# line; the source's m_sigma), rho(t) = (2 pi)^(-1) log(t/2 pi). The score is a single value of xi'/xi, so no sum over zeros and no truncation is needed:
# xi'/xi(s) = 1/s + 1/(s-1) - (log pi)/2 + psi(s/2)/2 + zeta'/zeta(s) in Arb (acb_series jet of zeta, which uses the Riemann Siegel expansion at
# large height), as enclosures, at two precisions (the difference is logged).
# Also computed around each height: the local minimum of s (ternary search in [t-0.05, t+0.05] at 64 bits and a scan of [t-0.3, t+0.3] with step 0.01),
# because the table prints heights to 0.01 and the minimum of the nested sampling run sits at such a height; and the normalised gap
# rho(t) (gamma_+ - gamma_-) between the consecutive zeros bracketing t, found by the sign changes of Hardy's Z function Z(t) = Re(e^(i theta(t)) zeta(1/2 + i t))
# on [t-3, t+3] (step 0.04, 64 bits) and bisected to 1e-12 in Arb (a zero is certified by Z of opposite sign at the ends of an interval of 1e-12;
# a missed pair of zeros closer than 0.04 would be missed, and the number of zeros found in the window is logged against the expected 6 rho).
# Usage: python3 -I cf_ns.py [--quick] [--out FILE]; --quick = the first height only, 64 and 96 bits (about 15 s); default: five heights, 128 and 192 bits,
# writes data/cf_ns.json and data/cf_ns.log. Cluster: bash misc/tools/hopper_run.sh -c 6 -t 120 -g cf_ns.json -g cf_ns.log code/cf_ns.py code/cf_eval.py -- --out cf_ns.json
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cf_eval import acb, arb, ctx, xi_logder, fl, density  # noqa: E402

HEIGHTS = [('1e6 band', '1009611.61', '0.175', '2.44'), ('1e7 band', '10005411.18', '0.182', '2.62'),
           ('1e9 band', '1000004287.61', '0.252', '2.54'), ('1e10 band', '10000009878.55', '0.292', '2.38'),
           ('1e6 band, plain Monte Carlo point', '1004606.8', '0.142', '2.87')]
SIG0 = arb(3) / 5


def rho(t):
    return (t / (2 * arb.pi())).log() / (2 * arb.pi())


def score(t):
    return density(xi_logder, SIG0, t) / rho(t)


def sfloat(t):
    return float(score(arb(t)).mid().str(17, radius=False))


def Z(t):
    th = (acb(1) / 4 + acb(0, t / 2)).lgamma().imag - t / 2 * arb.pi().log()
    return (acb(0, th).exp() * acb(arb(1) / 2, t).zeta()).real


def task(args):
    label, ts, smin, gap, p1, p2, quick = args
    out = {'label': label, 't': ts, 'printed_min_s': smin, 'printed_gap': gap}
    T0 = time.time()
    t = arb(ts)
    vals = []
    for p in (p1, p2):
        ctx.prec = p
        vals.append(score(arb(ts)))
    assert vals[0].overlaps(vals[1])
    out['s'] = vals[1].str(12)
    out["precision_diff"] = abs(float(vals[0].mid().str(25, radius=False)) - float(vals[1].mid().str(25, radius=False)))
    out['radii'] = [float(vals[0].rad()), float(vals[1].rad())]
    ctx.prec = 64
    tf = float(ts)
    # local minimum
    a, b = tf - 0.05, tf + 0.05
    for _ in range(40):
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        if sfloat(m1) < sfloat(m2):
            b = m2
        else:
            a = m1
    out['local_min_t'] = (a + b) / 2
    out['local_min_s'] = sfloat((a + b) / 2)
    n = 6 if quick else 60
    scan = [(sfloat(tf - 0.3 + 0.01 * i * (60 // n)), tf - 0.3 + 0.01 * i * (60 // n)) for i in range(n + 1)]
    out['scan_min'] = min(scan)[0]
    out['scan_min_t'] = min(scan)[1]
    # zeros bracketing t from Hardy Z
    ctx.prec = 64
    step = 0.04
    grid = [tf - 3 + step * i for i in range(int(6 / step) + 1)]
    zv = [Z(arb(g)) for g in grid]
    zeros = []
    ctx.prec = 128
    for i in range(len(grid) - 1):
        if (zv[i] > 0 and zv[i + 1] < 0) or (zv[i] < 0 and zv[i + 1] > 0):
            lo, hi = arb(grid[i]), arb(grid[i + 1])
            zl = Z(lo)
            while float(hi - lo) > 1e-12:
                mid = (lo + hi) / 2
                zm = Z(mid)
                if (zm > 0) == (zl > 0):
                    lo, zl = mid, zm
                else:
                    hi = mid
            zeros.append(float(lo.mid().str(20, radius=False)))
    out['zeros_in_window'] = zeros
    out['expected_zero_count'] = float(6 * rho(arb(tf)))
    below = [z for z in zeros if z < tf]
    above = [z for z in zeros if z > tf]
    if below and above:
        out['gamma_minus'], out['gamma_plus'] = below[-1], above[0]
        out['normalised_gap'] = float(rho(arb(tf)).mid()) * (above[0] - below[-1])
    out['seconds'] = time.time() - T0
    return out


def main():
    import multiprocessing as mpr
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or ('/tmp/cf_ns_quick.json' if a.quick else os.path.join(here, '..', 'data', 'cf_ns.json'))
    logf = os.path.splitext(out)[0] + '.log'
    hs = HEIGHTS[:1] if a.quick else HEIGHTS
    p1, p2 = (64, 96) if a.quick else (128, 192)
    with mpr.Pool(min(len(hs), os.cpu_count() or 1)) as pool:
        res = pool.map(task, [(h[0], h[1], h[2], h[3], p1, p2, a.quick) for h in hs])
    lines = []
    for r in res:
        lines.append('%s t = %s: s(t) = %s (difference of the two precisions %.1e as doubles; enclosure radii %.1e at %d bits, %.1e at %d bits); printed min s %s; local minimum of s near t: %.4f at t = %.4f; scan min %.4f at %.2f'
                     % (r['label'], r['t'], r['s'], r['precision_diff'], r['radii'][0], p1, r['radii'][1], p2, r['printed_min_s'], r['local_min_s'], r['local_min_t'], r['scan_min'], r['scan_min_t']))
        lines.append('   zeros found in [t-3, t+3]: %d (expected %.1f); bracketing zeros %s and %s; normalised gap %s (printed %s); %.0f s'
                     % (len(r['zeros_in_window']), r['expected_zero_count'], r.get('gamma_minus'), r.get('gamma_plus'),
                        '%.3f' % r['normalised_gap'] if 'normalised_gap' in r else 'n/a', r['printed_gap'], r['seconds']))
    print('\n'.join(lines))
    json.dump(res, open(out, 'w'), indent=1)
    open(logf, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
