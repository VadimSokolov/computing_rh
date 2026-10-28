# cf_moments.py: Table tab:dl of the source note L (Chapter 34, prefix cf, Section cf:sec:gbc): the contour moments of rectangles and the zeros
# recovered from them, as Arb interval enclosures. Group CF-XI, task T5; supports Proposition cf:prop:moments and the table of contour moments.
# Computed. For a rectangle R = [s0,s1] x [T1,T2] with L(s) = xi(s) (the zeta rows) or Lambda(s) (Davenport and Heilbronn, completion with
# Gamma((s+1)/2), see cf_eval.py): M_k(R) = (2 pi i)^(-1) oint_{dR} (s-c)^k L'/L(s) ds, k = 0..K, c the centre of R, the boundary positively
# oriented, four straight edges, each integrated by acb.integral (rigorous Gauss-Legendre with ball bounds, so every M_k is an enclosure,
# poles of L'/L being detected by non finite balls). K is the integer in the enclosure of M_0. Newton's identities (Proposition cf:prop:moments),
# j e_j = sum_{r=1..j} (-1)^(r-1) e_(j-r) M_r, give the polynomial z^K - e_1 z^(K-1) + ... whose roots are rho - c. Root enclosures: for an
# approximate root r (a root of the polynomial of the midpoints), some root of every polynomial of the ball family lies within K |P(r)|/|P'(r)|
# of r (P'/P = sum 1/(r - z_j)), with P, P' evaluated in ball arithmetic; if the K discs are pairwise disjoint each contains exactly one root.
# The approximate roots of the midpoint polynomial are polished by ten Newton steps in the working precision before the discs are formed.
# Rectangles. (1) xi [0.6,3]x[80,92]; (2) xi [0.55,3]x[60,120] (zeros at distance 0.05 from the left edge: the source's M_0 = 0.006 is quadrature
# error); (3) DH [0.6,3]x[80,92]; extra, not in the source: (4) DH [0.55,3]x[60,120] (K=2), (5) DH [0.55,3]x[80,180] (K=4: the four zeros
# off the line below 200 lying in it), (6) xi [0.55,3]x[80,92].
# Error control: every rectangle is run at two precisions (the second with tighter tolerance); the logged number is the largest distance between
# the midpoints, and the enclosures must overlap (asserted).
# written from the mathematics (the collaborator's scripts were not read, run or copied)
# Usage: python3 -I cf_moments.py [--quick] [--rects 0,1,2] [--out FILE]; the table rectangles are --rects 0,1,2 (data/cf_moments.json), the extra ones --rects 3,4,5 (data/cf_moments_extra.json); --quick = rectangles (1) and a short DH rectangle [0.6,3]x[84.5,87] at 53 and 64 bits (parallel, about 20 s);
# default writes data/cf_moments.json and data/cf_moments.log. Cluster: bash misc/tools/hopper_run.sh -c 1 -t 120 -g cf_moments.json
# -g cf_moments.log code/cf_moments.py code/cf_eval.py -- --out cf_moments.json
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cf_eval import acb, arb, ctx, xi_logder, dh_logder, fl, rad  # noqa: E402
from flint import acb_poly  # noqa: E402

KNOWN_DH = {1: ('0.8085171825', '85.6993484854'), 2: ('0.6508300806', '114.1633427308'),
            3: ('0.5743560505', '166.4793059132'), 4: ('0.7242576946', '176.7024612429')}

RECTS = [  # name, function, s0, s1, T1, T2 (rationals as strings)
    ('xi [0.6,3]x[80,92]', 'xi', '3/5', '3', '80', '92'),
    ('xi [0.55,3]x[60,120]', 'xi', '11/20', '3', '60', '120'),
    ('DH [0.6,3]x[80,92]', 'dh', '3/5', '3', '80', '92'),
    ('DH [0.55,3]x[60,120]', 'dh', '11/20', '3', '60', '120'),
    ('DH [0.55,3]x[80,180]', 'dh', '11/20', '3', '80', '180'),
    ('xi [0.55,3]x[80,92]', 'xi', '11/20', '3', '80', '92'),
]


QUICK_RECTS = [RECTS[0], ('DH [0.6,3]x[84.5,87]', 'dh', '3/5', '3', '84.5', '87')]


def q(x):
    if '/' in x:
        a, b = x.split('/')
        return arb(int(a)) / int(b)
    return arb(int(x)) if '.' not in x else arb(int(float(x) * 2)) / 2


def moment(g, c, k, corners, tol):
    tot = acb(0)
    for a, b in zip(corners[:-1], corners[1:]):
        tot += acb.integral(lambda z, an: (z - c) ** k * g(z), a, b, abs_tol=tol, rel_tol=tol,
                            eval_limit=10 ** 8, depth_limit=2000)
    return tot / (2 * arb.pi() * acb(0, 1))


def root_discs(es, K):
    """polynomial z^K - e1 z^(K-1) + ... with ball coefficients es[1..K]; returns list of (centre acb, radius float)"""
    coef = [acb(1)] + [(-1) ** j * es[j] for j in range(1, K + 1)]       # highest first
    P = acb_poly(coef[::-1])                                              # lowest first
    Pm = acb_poly([c.mid() for c in coef[::-1]])
    cents = []
    dPm = Pm.derivative()
    for r in Pm.roots():
        r = r.mid()
        for _ in range(10):      # Newton polish: roots() returns only a few digits
            r = (r - Pm(r) / dPm(r)).mid()
        cents.append(r)
    dP = P.derivative()
    out = []
    for r in cents:
        num = abs(P(r))
        den = abs(dP(r))
        R = K * num / den
        out.append((r, R))
    return out


def run_rect(name, fn, s0, s1, T1, T2, prec, tolexp, log):
    ctx.prec = prec
    g = xi_logder if fn == 'xi' else dh_logder
    s0, s1, T1, T2 = q(s0), q(s1), q(T1), q(T2)
    c = acb((s0 + s1) / 2, (T1 + T2) / 2)
    corners = [acb(s0, T1), acb(s1, T1), acb(s1, T2), acb(s0, T2), acb(s0, T1)]
    tol = arb(10) ** (-tolexp)
    t0 = time.time()
    M = [moment(g, c, 0, corners, tol)]
    K = int(round(float(M[0].real.mid())))
    assert M[0].real.contains(K) or abs(float(M[0].real.mid()) - K) < 1e-6, 'M_0 not near an integer'
    for k in range(1, K + 1):
        M.append(moment(g, c, k, corners, tol))
    es = [acb(1)]
    for j in range(1, K + 1):
        es.append(sum(((-1) ** (r - 1) * es[j - r] * M[r] for r in range(1, j + 1)), acb(0)) / j)
    zeros = []
    if K > 0:
        discs = root_discs(es, K)
        for r, R in discs:
            zeros.append((c + r, R))
        zeros.sort(key=lambda z: float(z[0].imag.mid()))
        for i in range(K):
            for j in range(i + 1, K):
                d = abs(zeros[i][0] - zeros[j][0])
                assert float(d.mid()) > float(zeros[i][1].mid() + zeros[j][1].mid()), 'discs not disjoint'
    return {'name': name, 'K': K, 'prec': prec, 'M': M, 'zeros': zeros, 'seconds': time.time() - t0}


def process_rect(args):
    r, precs = args
    lines = []
    log = lines.append
    runs = [run_rect(*r, p, t, log) for p, t in precs]
    r1, r2 = runs
    log('== %s: K = %d (prec %d: %.1f s; prec %d: %.1f s)' % (r[0], r1['K'], precs[0][0], r1['seconds'], precs[1][0], r2['seconds']))
    assert r1['K'] == r2['K']
    dmax = 0.0
    ent = {'K': r1['K'], 'M': [], 'zeros': []}
    for k, (m1, m2) in enumerate(zip(r1['M'], r2['M'])):
        d = abs(m1.mid() - m2.mid())
        dmax = max(dmax, float(d.mid()))
        assert (m1.real.overlaps(m2.real) and m1.imag.overlaps(m2.imag)), 'enclosures disjoint'
        log('   M_%d = %s   (enclosure radius %.1e at the higher precision; %s)' % (k, fl(m2, 12), rad(m2), m2.real.str(8) + ' ' + m2.imag.str(8) + 'i'))
        ent['M'].append({'k': k, 'value_str': m2.str(15), 'mid12': fl(m2, 12), 'radius': rad(m2)})
    for (z1, R1), (z2, R2) in zip(r1['zeros'], r2['zeros']):
        d = abs(z1 - z2)
        dmax = max(dmax, float(d.mid()))
        log('   zero = %s   enclosure radius %.2e' % (fl(z2, 14), float(R2.mid())))
        ent['zeros'].append({'mid': fl(z2, 14), 'radius': float(R2.mid())})
    log('   max midpoint difference between the two precisions: %.2e' % dmax)
    ent['precision_diff'] = dmax
    return r[0], lines, ent


def main():
    import multiprocessing as mpr
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--rects', default=None, help='comma separated indices of RECTS to run (default all)')
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or ('/tmp/cf_moments_quick.json' if a.quick else os.path.join(here, '..', 'data', 'cf_moments.json'))
    logf = os.path.splitext(out)[0] + '.log'
    rects = QUICK_RECTS if a.quick else RECTS
    if a.rects:
        rects = [RECTS[int(i)] for i in a.rects.split(',')]
    precs = [(53, 6), (64, 8)] if a.quick else [(128, 25), (192, 38)]
    with mpr.Pool(min(len(rects), os.cpu_count() or 1)) as pool:
        results = pool.map(process_rect, [(r, precs) for r in rects])
    lines = []
    res = {}
    for name, ls, ent in results:
        lines += ls
        res[name] = ent
    for k, (b, t) in KNOWN_DH.items():
        lines.append('known DH zero %d (source): %s + %s i' % (k, b, t))
    print('\n'.join(lines))
    json.dump(res, open(out, 'w'), indent=1)
    open(logf, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
