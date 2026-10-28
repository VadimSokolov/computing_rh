# cf_xi.py: the exact values for the window of Section cf:sec:u34 and for Table tab:negmass of the source note L (Chapter 34, prefix cf).
# Group CF-XI, task T6; supports the window text of Proposition cf:prop:determine and Theorem cf:thm:negmass. Written from the mathematics
# (the collaborator's scripts were not read, run or copied).
# What is computed. Density u_theta(v) = pi^(-1) Re (L'/L)(theta + i v) of the measure of the chapter on the line (no factor 1/2; the source's m_sigma
# is this function), L = xi (zeta) or Lambda (Davenport and Heilbronn, DH), both from cf_eval.py in Arb.
#  (A) Level theta = 3/4 on the window [78,94]: the total increment I = int_78^94 u_{3/4} dv (acb.integral, an enclosure, at two precisions) for
#      xi and for DH; the source's 8.76 is stated for DH (L:758-760 and the text of the window), its 6.16 for zeta; the minimum of u_{3/4} on the
#      grid of spacing 0.02 and (DH) the refined spike depth near 85.70; the negative mass of u_{3/4} for DH on the window.
#  (B) Table tab:negmass: the negative mass of u_theta for DH, theta in {0.9, 0.79, 0.7, 0.6}, in intervals of length 3 centred at the ordinates
#      85.6993484854 and 114.1633427308 of the two lowest zeros off the line, (i) as the source does, a Riemann sum of max(-u,0) on a grid of
#      spacing 0.001, and (ii) exactly: the sign changes are located by bisection in Arb and u is integrated by acb.integral over the intervals
#      where it is negative (an enclosure of the mass up to the root locations, which are bisected to 1e-25).
#  (C) The minimum of u_theta over [60,120] for DH at theta = 0.8 and 0.6 (the source's -36.6 and -4.6), by a ternary search near each ordinate.
#  (D) The curve u_{3/4} on [78,94], grid 0.02, for xi and DH (data for a redrawn figure), 12 digits.
# Error control: (A) at 128 and 192 bits with tolerance 1e-25 and 1e-38 (difference logged); (B)(ii) at 128 bits with the integral enclosure radius logged.
# Usage: python3 -I cf_xi.py [--quick] [--out FILE]; --quick = grid 0.01 for (B), no (C), curve step 0.1, 64 bits (about 25 s); default writes
# data/cf_xi.json and data/cf_xi.log. Cluster: bash misc/tools/hopper_run.sh -c 8 -t 120 -g cf_xi.json -g cf_xi.log code/cf_xi.py code/cf_eval.py -- --out cf_xi.json
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cf_eval import acb, arb, ctx, xi_logder, dh_logder, fl, density  # noqa: E402

G1 = arb('85.6993484854')     # the source's ordinates of the two lowest DH zeros off the line (10 digits)
G2 = arb('114.1633427308')
LEVELS = {'0.9': arb(9) / 10, '0.79': arb(79) / 100, '0.7': arb(7) / 10, '0.6': arb(3) / 5, '0.75': arb(3) / 4, '0.8': arb(4) / 5}


def uf(fn, theta, t):
    g = xi_logder if fn == 'xi' else dh_logder
    return density(g, theta, t)


def ufloat(fn, theta, t):
    return float(uf(fn, theta, arb(t)).mid().str(17, radius=False))


def integral(fn, theta, a, b, tolexp):
    g = xi_logder if fn == 'xi' else dh_logder
    tol = arb(10) ** (-tolexp)
    I = acb.integral(lambda z, an: g(acb(theta) + acb(0, 1) * z), acb(a), acb(b), abs_tol=tol, rel_tol=tol, eval_limit=10 ** 8, depth_limit=2000)
    return I.real / arb.pi()


def bisect_root(fn, theta, lo, hi, iters):
    flo = uf(fn, theta, lo)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = uf(fn, theta, mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def negmass(args):
    """(label, level string, centre string, step, prec): grid Riemann sum and exact negative mass of DH u_theta in [c-1.5, c+1.5]"""
    label, lv, cs, h, prec = args
    ctx.prec = prec
    theta = LEVELS[lv]
    c = G1 if cs == '1' else G2
    lo, hi = c - arb(3) / 2, c + arb(3) / 2
    n = int(round(3 / h))
    lo_f = float(lo.mid().str(17, radius=False))
    ts = [lo_f + i * h for i in range(n + 1)]
    vals = [ufloat('dh', theta, t) for t in ts]
    grid = sum(max(-v, 0.0) for v in vals) * h
    umin = min(vals)
    brk = [i for i in range(n) if (vals[i] > 0) != (vals[i + 1] > 0)]
    roots = [bisect_root('dh', theta, arb(ts[i]), arb(ts[i + 1]), 85) for i in brk]
    pts = [lo] + roots + [hi]
    exact = arb(0)
    nneg = 0
    for a, b in zip(pts[:-1], pts[1:]):
        m = (a + b) / 2
        if uf('dh', theta, m) < 0:
            exact -= integral('dh', theta, a, b, 20)
            nneg += 1
    return label, grid, exact.str(10, radius=True), float(exact.mid().str(15, radius=False)), len(roots), nneg, umin


def tern(args):
    fn, lv, a, b, prec = args
    ctx.prec = prec
    theta = LEVELS[lv]
    a, b = float(a), float(b)
    for _ in range(80):
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        if ufloat(fn, theta, m1) < ufloat(fn, theta, m2):
            b = m2
        else:
            a = m1
    t = (a + b) / 2
    return fn, lv, t, ufloat(fn, theta, t)


def curve(args):
    t0, t1, h, prec = args
    ctx.prec = prec
    out = []
    n0 = int(round(t0 / h))
    for i in range(n0, int(round(t1 / h)) + 1):
        t = arb(i) * arb(str(h)) if False else arb(i) / int(round(1 / h))
        out.append((fl(t, 10), fl(uf('xi', arb(3) / 4, t), 12), fl(uf('dh', arb(3) / 4, t), 12)))
    return out


def main():
    import multiprocessing as mpr
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or ('/tmp/cf_xi_quick.json' if a.quick else os.path.join(here, '..', 'data', 'cf_xi.json'))
    logf = os.path.splitext(out)[0] + '.log'
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    res = {}
    ncpu = os.cpu_count() or 1
    T0 = time.time()
    prec1, tol1, prec2, tol2 = (64, 10, 96, 16) if a.quick else (128, 25, 192, 38)
    # (A) total increments
    th = arb(3) / 4
    A = {}
    for fn in ('xi', 'dh'):
        vals = []
        for p, tl in ((prec1, tol1), (prec2, tol2)):
            ctx.prec = p
            vals.append(integral(fn, arb(3) / 4, arb(78), arb(94), tl))
        d = abs(vals[0].mid() - vals[1].mid())
        assert vals[0].overlaps(vals[1])
        log('(A) %s: int_78^94 u_3/4 dv = %s  (precision difference %.1e)' % (fn, vals[1].str(15), float(d.mid())))
        A[fn] = {'total_increment': vals[1].str(15), 'precision_diff': float(d.mid())}
    ctx.prec = 64
    step = 0.1 if a.quick else 0.02
    for fn in ('xi', 'dh'):
        n = int(round(16 / step))
        g = [(ufloat(fn, th, 78 + i * step), 78 + i * step) for i in range(n + 1)]
        m = min(g)
        log('(A) %s: minimum of u_3/4 on the grid of spacing %g: %.4f at t = %.2f' % (fn, step, m[0], m[1]))
        A[fn]['grid_min'] = m[0]
        A[fn]['grid_min_t'] = m[1]
    res['A'] = A
    # (B) tab:negmass
    h = 0.01 if a.quick else 0.001
    tasks = [('%s@%s' % (lv, cs), lv, cs, h, 128) for cs in ('1', '2') for lv in ('0.9', '0.79', '0.7', '0.6')]
    tasks.append(('0.75@win', '0.75', '1', h, 128))
    if a.quick:
        tasks = [tasks[1], tasks[8]]
    with mpr.Pool(min(len(tasks) + 2, ncpu)) as pool:
        nm = pool.map(negmass, tasks)
        B = {}
        for label, grid, es, ev, nroots, nneg, umin in nm:
            log('(B) DH level@zero %s: negative mass, grid %g Riemann sum %.4f; exact %s; sign changes %d, negative pieces %d; min u %.3f' % (label, h, grid, es, nroots, nneg, umin))
            B[label] = {'grid_sum': grid, 'exact': es, 'exact_mid': ev, 'sign_changes': nroots, 'negative_pieces': nneg, 'min_u': umin, 'step': h}
        res['B'] = B
        # (C)
        if not a.quick:
            tt = [('dh', '0.8', 85.6, 85.8, 128), ('dh', '0.6', 85.6, 85.8, 128), ('dh', '0.6', 114.05, 114.28, 128), ('dh', '0.8', 114.05, 114.28, 128)]
            C = []
            for fn, lv, t, v in pool.map(tern, tt):
                log('(C) DH u_%s local minimum near t = %.6f: %.4f' % (lv, t, v))
                C.append({'level': lv, 't': t, 'u': v})
            res['C'] = C
        # (D)
        h2 = 0.1 if a.quick else 0.02
        cs = [(78 + 2 * j, 78 + 2 * j + (2 - h2 if j < 7 else 2), h2, 96) for j in range(8)]
        parts = pool.map(curve, [(t0 if False else 78 + 2 * j, 78 + 2 * j + 2, h2, 96) for j in range(8)])
    rows = []
    seen = set()
    for p in parts:
        for r in p:
            if r[0] not in seen:
                seen.add(r[0])
                rows.append(r)
    res['D'] = {'description': 'rows [t, u_3/4 for xi, u_3/4 for DH], 12 digits', 'rows': rows}
    log('(D) curve rows written: %d' % len(rows))
    log('total seconds %.1f' % (time.time() - T0))
    json.dump(res, open(out, 'w'), indent=1)
    open(logf, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
