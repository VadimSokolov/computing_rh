#!/usr/bin/env python3
"""D2: two sided bracket of -ln eps_N from the zeros alone, beyond the grid.
Below: -ln RQ of the test function of Theorem 5 (zero side, first 8000 ordinates, density with mean sin^2 = 1/2 above height 8148: a model).
Above: LB of Theorem 1(c) with equilibrium nodes taken from the first 8000 ordinates (under RH; rigorous apart from floating point).
Nodes whose equilibrium position lies above about 0.6 of the top are spread geometrically below the top, then a short greedy exchange.
Run: timeout 590 python3 -I D2_bracket.py
"""
import json, math, sys
import numpy as np
OPUS = './misc/research/directions/opus'
sys.path.insert(0, OPUS)
import D2_lower as DL
import D2_testfunction as DT
gam = DL.gam
TOP = gam[-1]


def nodes(N, L):
    phi, b = DL.params(N + 1, L)
    K = int(np.searchsorted(gam, b, side='right'))
    if K >= N + 1:
        return gam[:N + 1].copy(), K, b
    u = np.linspace(0, 40, 400001)
    s = b * np.exp(u)
    qs = DL.q_band(s, phi, b, L) * s
    du = u[1] - u[0]
    tail = np.concatenate([np.cumsum(0.5 * (qs[1:] + qs[:-1])[::-1] * du)[::-1], [0.0]])
    A = (N + 1) - tail
    J = N + 1 - K
    betas = np.interp(K + np.arange(1, J + 1) - 0.5, A, s)
    cut = 0.6 * TOP
    over = betas > cut
    if over.any():
        j0 = int(np.argmax(over))
        m = J - j0
        lo = max(betas[j0 - 1] if j0 > 0 else b, cut * 0.9)
        betas[j0:] = lo * (TOP * 0.999 / lo) ** ((np.arange(m) + 1) / m)
    chosen = list(gam[:K])
    idx = K
    for bj in betas:
        i = max(idx, int(np.searchsorted(gam, bj)))
        chosen.append(gam[min(i, len(gam) - 1)])
        idx = min(i, len(gam) - 1) + 1
    S = np.array(chosen)
    if len(np.unique(S)) < N + 1 or idx > len(gam):
        return None, K, b
    return S, K, b


def main():
    print('    X     N     K       b   -ln RQ      LB   gap    E(N)   law fit   grid   worst node')
    res = []
    for X, N in [(5, 16), (13, 66), (50, 391), (1000, 600), (100, 921), (200, 1000), (200, 1500), (200, 2119), (1000, 1000), (1000, 2000), (1000, 4000), (500, 6000)]:
        L = 0.5 * math.log(X)
        S, K, b = nodes(N, L)
        if S is None:
            print(X, N, 'no node set inside the list')
            continue
        v, t = DL.LB(N, L, S)
        if N <= 1100:
            v, S, t = DL.improve(N, L, S, 6)
        rq = DT.rayleigh(N, L, gam)[0]
        Eh = DL.params(N + 0.5, L)
        import D2_fits as DF
        E0 = DF.E(N, L)[0]
        Ef, bf, pf = DF.E(N + 0.5, L)
        fit = Ef - 4.552 * math.log(bf) - 5.157 - 0.675 * math.log(L)
        g = DL.grid.get((X, N))
        gv = g['minus_ln_eps'] if g else float('nan')
        res.append(dict(X=X, N=N, K=K, b=b, rq=rq, LB=v, E=E0, fit=fit, grid=gv))
        print(f"{X:6d} {N:5d} {K:5d} {b:7.1f} {rq:8.2f} {v:8.2f} {v - rq:5.2f} {E0:8.2f} {fit:8.2f} {gv:8.2f} {S[int(np.argmax(t))]:8.1f}", flush=True)
    json.dump(res, open(OPUS + '/D2_bracket.json', 'w'), indent=0)


if __name__ == '__main__':
    main()
