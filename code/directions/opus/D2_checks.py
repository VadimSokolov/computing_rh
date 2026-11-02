#!/usr/bin/env python3
"""D2: attempts to refute Corollary 2 and Theorem 5 on the grid, and the lower bound of Theorem 1(c) (under RH) at sample cells.

Run: timeout 300 python3 -I D2_checks.py   (after D2_law_vs_grid.py)
Status: floating point in log space. The lower bound assumes RH and uses the first N+1 ordinates, or the N+1 best of the first 2N.
"""
import json, math
import numpy as np
ROOT = '.'
OPUS = ROOT + '/misc/research/directions/opus'
gam = np.array([float(x) for x in open(ROOT + '/data/zeros_8000.txt') if x.strip()])
rows = [r for r in json.load(open(OPUS + '/D2_law_vs_grid.json')) if r['minus_ln_eps'] >= 1]


def law(Nx, L):
    T = math.pi * Nx / L
    if T >= 4 * math.pi * math.exp(2 * L):
        return 4 * math.pi * math.exp(2 * L)
    lo, hi = 1e-12, math.pi / 2
    g = lambda p: (4 * L / math.pi) * p - math.log(T / (4 * math.pi * math.sin(p)))
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if g(m) > 0 else (m, hi)
    p = 0.5 * (lo + hi)
    return T / math.sin(p) + 4 * Nx * math.log(1 / math.tan(p / 2))


def lower_bound(N, L, S):
    """ln( sum ||l_gamma||^2 / 2 ), so that eps_N >= exp(-this) under RH; S = N+1 distinct ordinates, not in the lattice."""
    y = (math.pi * np.arange(0, N + 1) / L) ** 2
    kap = np.ones(N + 1)
    kap[0] = 2
    # log prod_{j != n} |y_n - y_j|
    d = np.abs(y[:, None] - y[None, :])
    np.fill_diagonal(d, 1.0)
    logd = np.log(d).sum(1)
    s2 = S ** 2
    terms = []
    for i, g in enumerate(S):
        others = np.delete(s2, i)
        logP = np.log(np.abs(y[:, None] - others[None, :])).sum(1)
        lognorm = np.logaddexp.reduce(np.log(L * kap) + 2 * logP - 2 * logd)
        logc = np.log(np.abs(g * g - y)).sum() - math.log(2 * g) - math.log(abs(math.sin(L * g))) - np.log(np.abs(g * g - others)).sum()
        terms.append(2 * logc + lognorm)
    return np.logaddexp.reduce(np.array(terms)) - math.log(2)


def main():
    print('Corollary 2(b): cells with N <= X, theta = log N/log X: check -ln eps >= c_theta N')
    bad = []
    ratio = []
    for r in rows:
        X, N = r['X'], r['N']
        if N > X:
            continue
        th = math.log(N) / math.log(X)
        c = 4 * math.log(1 / math.tan(math.pi * th / 4))
        ratio.append(r['minus_ln_eps'] / (c * N))
        if r['minus_ln_eps'] < c * N:
            bad.append((X, N, round(th, 3), round(c * N, 1), round(r['minus_ln_eps'], 1)))
    print('  cells', len(ratio), 'violations', len(bad), bad[:12], 'min ratio', round(min(ratio), 3), 'median', round(float(np.median(ratio)), 3))
    print('Corollary 2(a): cells with N <= X: -ln eps - c_*(X) N, min')
    v = []
    for r in rows:
        X, N = r['X'], r['N']
        if N <= X and X >= 13:
            v.append((r['minus_ln_eps'] - law(X, r['L']) / X * N, X, N))
    print('  ', min(v), 'number negative', sum(1 for t in v if t[0] < 0), 'of', len(v))
    print('Theorem 5: largest (E - (-ln eps)) / log^2(2+N), the least C the data allow')
    w = [((r['E'] - r['minus_ln_eps']) / math.log(2 + r['N']) ** 2, r['X'], r['N']) for r in rows]
    print('  ', max(w), 'violations of -ln eps >= E - C log^2 with C = 4:', sum(1 for t in w if t[0] > 4))
    print('Theorem 1(c), under RH: X, N, -ln eps(grid), LB first N+1, LB best N+1 of first 2N by |sin L gamma|, E')
    for X, N in [(13, 8), (13, 20), (50, 20), (1000, 8), (1000, 16), (1000, 30), (1000, 100), (1000000, 16), (1000000, 100)]:
        rr = [r for r in rows if r['X'] == X and r['N'] == N]
        if not rr:
            continue
        L = 0.5 * math.log(X)
        lb1 = lower_bound(N, L, gam[:N + 1])
        cand = gam[:2 * N + 2]
        keep = np.sort(cand[np.argsort(-np.abs(np.sin(L * cand)))[:N + 1]])
        lb2 = lower_bound(N, L, keep)
        print(f"  {X:8d} {N:4d} {rr[0]['minus_ln_eps']:10.2f} {lb1:10.2f} {lb2:10.2f} {rr[0]['E']:10.2f}")


if __name__ == '__main__':
    main()
