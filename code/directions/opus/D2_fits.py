#!/usr/bin/env python3
"""D2: refined fit of the margin law on all grid cells, the accuracy law for the zeros of the ground state,
the proved bound against every cell, and the comparison with Zhu's fitted law.

Run: timeout 120 python3 -I D2_fits.py   (after D2_law_vs_grid.py)
Status: floating point; the offset constants are fitted, the leading law has no free parameter.
"""
import json, math
import numpy as np
from scipy import integrate
ROOT = '.'
OPUS = ROOT + '/misc/research/directions/opus'


def params(Nx, L):
    T = math.pi * Nx / L
    lo, hi = 1e-12, math.pi / 2
    g = lambda p: (4 * L / math.pi) * p - math.log(T / (4 * math.pi * math.sin(p)))
    if g(hi) <= 0:
        return math.pi / 2, 4 * math.pi * math.exp(2 * L)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if g(m) > 0 else (m, hi)
    p = 0.5 * (lo + hi)
    return p, T / math.sin(p)


def E(Nx, L):
    p, b = params(Nx, L)
    return b + 4 * Nx * math.log(1 / math.tan(p / 2)), b, p


def U_closed(t, Nx, L):
    """U(t) of the equilibrium, t < b: int_t^b ( (1/2) arccos(s/b) + (L/pi) log|sin(phi+th)/sin(phi-th)| ) ds."""
    p, b = params(Nx, L)
    if t >= b:
        return 0.0

    def f(s):
        th = math.asin(min(s / b, 1.0))
        val = 0.5 * math.acos(min(s / b, 1.0))
        if p < math.pi / 2 - 1e-12:
            val += (L / math.pi) * math.log(abs(math.sin(p + th) / math.sin(p - th)))
        return val
    pts = [b * math.sin(p)] if t < b * math.sin(p) < b else None
    return integrate.quad(f, t, b, points=pts, limit=400)[0]


def main():
    rows = json.load(open(OPUS + '/D2_law_vs_grid.json'))
    # 1. bound of Theorem 2 against every cell
    bad = [(r['X'], r['N']) for r in rows if r['bound_ln'] is not None and -r['bound_ln'] > r['minus_ln_eps'] + 1e-6]
    usable = [r for r in rows if r['bound_ln'] is not None]
    print('Theorem 2 cells', len(usable), 'violations', bad, 'median ratio bound/-ln eps',
          round(float(np.median([-r['bound_ln'] / r['minus_ln_eps'] for r in usable if r['minus_ln_eps'] > 5])), 3))
    # 2. refined fit
    A, y, lab = [], [], []
    for r in rows:
        if r['N'] < 8 or r['minus_ln_eps'] < 1:
            continue
        L = r['L']
        Eh, b, p = E(r['N'] + 0.5, L)
        A.append([math.log(b), 1.0, math.log(L)])
        y.append(r['minus_ln_eps'] - Eh)
        lab.append((r['X'], r['N'], r['minus_ln_eps'], Eh))
    A, y = np.array(A), np.array(y)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    print('fit -ln eps - E(N+1/2,L) = a ln b + c + d ln L: a,c,d =', np.round(coef, 3), 'rms', round(float(res.std()), 3), 'max', round(float(np.abs(res).max()), 3), 'cells', len(y))
    A2 = A[:, :2]
    c2, *_ = np.linalg.lstsq(A2, y, rcond=None)
    res2 = y - A2 @ c2
    print('fit without ln L: a,c =', np.round(c2, 3), 'rms', round(float(res2.std()), 3), 'max', round(float(np.abs(res2).max()), 3))
    y3 = y + 4.5 * A[:, 0]
    print('with a = -9/2 fixed: c mean', round(float(y3.mean()), 3), 'rms', round(float(y3.std()), 3), 'range', round(float(y3.min()), 2), round(float(y3.max()), 2))
    # leading order only: relative error
    rel = [abs(r['minus_ln_eps'] - r['E']) / r['minus_ln_eps'] for r in rows if r['minus_ln_eps'] > 200]
    print('relative error of the bare law on cells with -ln eps > 200: median', round(float(np.median(rel)), 4), 'max', round(float(np.max(rel)), 4), 'n', len(rel))
    # 3. saturation: the book's cells against N_* = 2 X ln X
    print('saturation: X, N_*, cells (N, -ln eps)')
    for X in (2, 5, 13, 50, 200):
        cells = sorted((r['N'], round(r['minus_ln_eps'], 2)) for r in rows if r['X'] == X and r['N'] >= 20)
        print(' ', X, round(2 * X * math.log(X), 1), cells[:9])
    # 4. accuracy of the zeros of the ground state
    acc = []
    gam = [float(x) for x in open(ROOT + '/data/zeros_1700.txt') if x.strip()]
    p = json.load(open(ROOT + '/misc/reviews/nick-runs/vadim_runs/pilot/p13_40.json'))
    errs = [float(e) for e in p['err']]
    L = p['L']
    print('accuracy at X=13, N=40: k, gamma_k, -ln err, 2U(gamma_k)')
    for k in (1, 5, 10, 15, 20, 25, 28):
        print('  ', k, round(gam[k - 1], 2), round(-math.log(errs[k - 1]), 1), round(2 * U_closed(gam[k - 1], 40.5, L), 1))
    wl = json.load(open(ROOT + '/data/wl_acc.json'))
    xs, ys, labs, zz = [], [], [], []
    for e in wl:
        if e.get('kind') != 'zeta' or e.get('parity') != 'even':
            continue
        errlist = e.get('err_log10')
        if not errlist:
            continue
        L = 0.5 * math.log(e['X'])
        pb, bb = params(e['N'] + 0.5, L)
        for k, lg in enumerate(errlist, start=1):
            if lg is None or lg > -3:
                continue
            xs.append(2 * U_closed(gam[k - 1], e['N'] + 0.5, L))
            ys.append(-lg * math.log(10))
            zz.append(math.log(bb / gam[k - 1]))
            labs.append((e['X'], e['N'], k))
    if xs:
        xs, ys = np.array(xs), np.array(ys)
        slope, icpt = np.polyfit(xs, ys, 1)
        r = ys - (slope * xs + icpt)
        print('accuracy law over wl_acc: entries', len(xs), 'fit -ln err = %.4f * 2U + %.2f' % (slope, icpt), 'rms %.2f (nats) = %.2f digits' % (r.std(), r.std() / math.log(10)))
        zz = np.array(zz)
        A = np.column_stack([xs, zz, np.ones_like(xs)])
        co, *_ = np.linalg.lstsq(A, ys, rcond=None)
        r2 = ys - A @ co
        print('   fit -ln err = a 2U + c ln(b/gamma) + d: a,c,d =', np.round(co, 4), 'rms %.2f nats = %.2f digits, max %.2f digits' % (r2.std(), r2.std() / math.log(10), np.abs(r2).max() / math.log(10)))
        d = ys - xs
        print('   no fit: -ln err - 2U: mean %.2f, rms %.2f nats' % (d.mean(), d.std()))
    else:
        print('wl_acc keys', list(wl[0].keys()))
    # 5. Zhu's law against 4 pi X
    print('a, X, 4piX - 4.5 ln(4piX) - c, Zhu 20.13 N/log N, Zhu 2pi^2 N/log N')
    cfit = float(y3.mean())
    for a in (0.8, 1.0, 1.2, 1.6, 2.0, 2.3026, 2.649, 3.0, 3.454, 5.0):
        X = math.exp(2 * a)
        Tst = 2 * math.pi * X
        Nst = Tst / (2 * math.pi) * math.log(Tst / (2 * math.pi * math.e)) + 7 / 8
        law = 4 * math.pi * X - 4.5 * math.log(4 * math.pi * X) + cfit
        print(f"  {a:6.3f} {X:10.1f} {law:12.1f} {20.13 * Nst / math.log(Nst):12.1f} {2 * math.pi ** 2 * Nst / math.log(Nst):12.1f}  N*={2 * X * math.log(X):.0f}")


if __name__ == '__main__':
    main()
