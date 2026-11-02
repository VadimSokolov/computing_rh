#!/usr/bin/env python3
"""D2: LB of Theorem 1(c) for N > N_* (and the law values at three saturated cells). Equilibrium for N >= N_*: A = N_0 up to b = 4 pi X, then
A' = h' + q_inf, q_inf = log(1 + omega)/(2 pi), so the nodes follow the lattice with a small excess. U_N = U_inf for every N >= N_*.
Run: timeout 500 python3 -I D2_bracket3.py /path/to/zeros1
Status: floating point, under RH, published ordinates; heuristic choice of nodes (the bound holds for any distinct nodes).
"""
import json, math, sys
import numpy as np
OPUS = './misc/research/directions/opus'
sys.path.insert(0, OPUS)
import D2_bracket2 as B2
import D2_lower as DL
import D2_fits as DF
gam = B2.gam
TOP = gam[-1]


def nodes_big(N, L):
    X = math.exp(2 * L)
    b = 4 * math.pi * X
    K = int(np.searchsorted(gam, b, side='right'))
    N0b = b / (2 * math.pi) * math.log(b / (2 * math.pi * math.e))
    u = np.linspace(0, 40, 800001)
    s = b * np.exp(u)
    r = s / b
    e = r + np.sqrt(np.maximum(r * r - 1, 0.0))
    w = e ** -2
    hp = np.where(s <= math.pi * (N + 1) / L, L / math.pi, 0.0)
    dens = (hp + np.log1p(w) / (2 * math.pi)) * s
    du = u[1] - u[0]
    A = N0b + np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * du)])
    J = N + 1 - K
    betas = np.interp(K + np.arange(1, J + 1) - 0.5, A, s)
    cut = 0.6 * TOP
    over = betas > cut
    nover = int(over.sum())
    if nover:
        j0 = int(np.argmax(over))
        m = J - j0
        lo = max(betas[j0 - 1], cut * 0.9)
        betas[j0:] = lo * (TOP * 0.999 / lo) ** ((np.arange(m) + 1) / m)
    chosen = list(gam[:K])
    idx = K
    for bj in betas:
        i = max(idx, int(np.searchsorted(gam, bj)))
        chosen.append(gam[min(i, len(gam) - 1)])
        idx = min(i, len(gam) - 1) + 1
    return np.array(chosen), K, b, nover


def main():
    print('law at (1000,6000), (500,6214), (1000,13800):')
    for X, N in [(1000, 6000), (500, 6214), (1000, 13800)]:
        L = 0.5 * math.log(X)
        Ef, bf, pf = DF.E(N + 0.5, L)
        print('  ', X, N, 'E', round(DF.E(N, L)[0], 2), 'fit', round(Ef - 4.552 * math.log(bf) - 5.157 - 0.675 * math.log(L), 2))
    print('    X     N   N_*   K  squeezed   LB   grid   worst node')
    for X, N in [(5, 100), (5, 600), (13, 100), (13, 200), (13, 600), (50, 600), (50, 2000), (100, 1500), (100, 4000), (200, 3000), (200, 6000), (1000, 20000)]:
        L = 0.5 * math.log(X)
        S, K, b, nover = nodes_big(N, L)
        v, t = B2.lb_light(N, L, S)
        g = DL.grid.get((X, N))
        print(f"{X:6d} {N:6d} {2 * X * math.log(X):8.1f} {K:6d} {nover:4d} {v:10.2f} {(g['minus_ln_eps'] if g else float('nan')):9.2f} {S[int(np.argmax(t))]:9.1f}", flush=True)


if __name__ == '__main__':
    main()
