#!/usr/bin/env python3
"""D2: the lower bound of Theorem 1(c) with equilibrium nodes: S = all ordinates up to b(N+1,L), then in the band the ordinates nearest the
quantiles of the equilibrium density q, then a greedy exchange that lowers the bound. Under RH, eps_N >= exp(-LB) for every such S.
Run: timeout 600 python3 -I D2_lower.py
Status: floating point in log space, assumes RH for the inequality; S uses the first 8000 ordinates only (nodes capped below height 8148).
"""
import json, math, sys
import numpy as np
ROOT = '.'
OPUS = ROOT + '/misc/research/directions/opus'
gam = np.array([float(x) for x in open(ROOT + '/data/zeros_8000.txt') if x.strip()])
grid = {(r['X'], r['N']): r for r in json.load(open(OPUS + '/D2_law_vs_grid.json'))}


def params(Nx, L):
    T = math.pi * Nx / L
    if T >= 4 * math.pi * math.exp(2 * L):
        return math.pi / 2, 4 * math.pi * math.exp(2 * L)
    lo, hi = 1e-12, math.pi / 2
    g = lambda p: (4 * L / math.pi) * p - math.log(T / (4 * math.pi * math.sin(p)))
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if g(m) > 0 else (m, hi)
    p = 0.5 * (lo + hi)
    return p, T / math.sin(p)


def q_band(x, phi, b, L):
    r = x / b
    e = r + np.sqrt(np.maximum(r * r - 1, 0.0))
    w = e ** -2
    tanh = (1 - w) / (1 + w)
    tp = math.tan(min(phi, math.pi / 2 - 1e-14))
    d = np.arctan(tp * (2 * w / (1 + w)) / (1 + tp * tp * tanh))
    return (0.5 * np.log1p(w) + (2 * L / math.pi) * d) / math.pi


def nodes(N, L):
    """N+1 distinct ordinates: all up to b, then nearest the quantiles of q."""
    phi, b = params(N + 1, L)
    K = int(np.searchsorted(gam, b, side='right'))
    if K >= N + 1:
        return gam[:N + 1].copy()
    u = np.linspace(0, 40, 400001)
    s = b * np.exp(u)
    qs = q_band(s, phi, b, L) * s
    du = u[1] - u[0]
    tail = np.concatenate([np.cumsum(0.5 * (qs[1:] + qs[:-1])[::-1] * du)[::-1], [0.0]])
    A = (N + 1) - tail
    J = N + 1 - K
    targets = K + np.arange(1, J + 1) - 0.5
    betas = np.interp(targets, A, s)
    chosen = list(gam[:K])
    idx = K
    for bj in betas:
        i = int(np.searchsorted(gam, bj))
        i = max(idx, min(i, len(gam) - 1 - (J - 1)))
        i = min(i, len(gam) - 1)
        if i < idx:
            i = idx
        chosen.append(gam[i])
        idx = i + 1
    return np.array(chosen)


def lb_terms(N, L, S):
    y = (math.pi * np.arange(0, N + 1) / L) ** 2
    kap = np.ones(N + 1)
    kap[0] = 2
    d = np.abs(y[:, None] - y[None, :])
    np.fill_diagonal(d, 1.0)
    logd = np.log(d).sum(1)
    s2 = S ** 2
    terms = np.zeros(len(S))
    logPall = np.log(np.abs(y[:, None] - s2[None, :]))
    sumlogP = logPall.sum(1)
    D = np.abs(s2[:, None] - s2[None, :])
    np.fill_diagonal(D, 1.0)
    logD = np.log(D).sum(1)
    logY = np.log(np.abs(s2[:, None] - y[None, :])).sum(1)
    for i, g in enumerate(S):
        logP = sumlogP - logPall[:, i]
        lognorm = np.logaddexp.reduce(np.log(L * kap) + 2 * logP - 2 * logd)
        logc = logY[i] - math.log(2 * g) - math.log(abs(math.sin(L * g))) - logD[i]
        terms[i] = 2 * logc + lognorm
    return terms


def LB(N, L, S):
    t = lb_terms(N, L, S)
    return np.logaddexp.reduce(t) - math.log(2), t


def improve(N, L, S, rounds=40):
    """greedy: move the worst node to a neighbouring free ordinate if the bound drops."""
    S = S.copy()
    cur, t = LB(N, L, S)
    for _ in range(rounds):
        i = int(np.argmax(t))
        best = None
        k = int(np.searchsorted(gam, S[i]))
        used = set(np.round(S, 9))
        for j in list(range(max(0, k - 4), min(len(gam), k + 5))):
            if round(gam[j], 9) in used:
                continue
            T = S.copy()
            T[i] = gam[j]
            T.sort()
            v, tt = LB(N, L, T)
            if v < cur - 1e-6 and (best is None or v < best[0]):
                best = (v, T, tt)
        if best is None:
            break
        cur, S, t = best
    return cur, S, t


def main():
    print('    X    N  K  b  -ln eps(grid)  LB(first N+1)  LB(equilibrium)  LB(after exchange)  E  worst node, its height')
    out = []
    cells = [(13, 8), (13, 20), (13, 40), (50, 20), (50, 40), (50, 100), (200, 60), (200, 100), (1000, 8), (1000, 16), (1000, 30), (1000, 100), (1000, 200), (100000, 100), (1000000, 16), (1000000, 100)]
    for X, N in cells:
        if (X, N) not in grid:
            continue
        L = 0.5 * math.log(X)
        S0 = gam[:N + 1]
        v0, _ = LB(N, L, S0)
        S = nodes(N, L)
        v1, t1 = LB(N, L, S)
        rounds = 40 if N <= 40 else 12
        v2, S2, t2 = improve(N, L, S, rounds)
        phi, b = params(N + 1, L)
        K = int(np.searchsorted(gam, b, side='right'))
        e = grid[(X, N)]
        out.append(dict(X=X, N=N, K=K, b=b, grid=e['minus_ln_eps'], LB_first=v0, LB_eq=v1, LB_ex=v2, E=e['E']))
        print(f"{X:8d} {N:4d} {K:4d} {b:7.1f} {e['minus_ln_eps']:9.2f} {v0:9.2f} {v1:9.2f} {v2:9.2f} {e['E']:9.2f}  {S2[int(np.argmax(t2))]:8.2f}", flush=True)
    json.dump(out, open(OPUS + '/D2_lower.json', 'w'), indent=0)


if __name__ == '__main__':
    main()
