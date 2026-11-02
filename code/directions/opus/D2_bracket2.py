#!/usr/bin/env python3
"""D2: the upper end LB of the bracket (Theorem 1(c), equilibrium nodes) with the first 100000 ordinates (Odlyzko's table, about 9 decimals),
at cells near saturation. The table is not in the repository: pass its path.
Run: timeout 590 python3 -I D2_bracket2.py /path/to/zeros1
Status: floating point; under RH; ordinates from a published table, not certified here.
"""
import json, math, sys
import numpy as np
OPUS = './misc/research/directions/opus'
sys.path.insert(0, OPUS)
import D2_lower as DL
gam8 = DL.gam
tab = np.array([float(x) for x in open(sys.argv[1]) if x.strip()])
gam = np.concatenate([gam8, tab[8000:]])
DL.gam = gam
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
    nover = int(over.sum())
    if nover:
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
    return np.array(chosen), K, b, nover


def lb_light(N, L, S):
    y = (math.pi * np.arange(0, N + 1) / L) ** 2
    kap = np.ones(N + 1)
    kap[0] = 2
    s2 = S ** 2
    logd = np.zeros(N + 1)
    sumlogP = np.zeros(N + 1)
    logY = np.zeros(N + 1)
    logD = np.zeros(N + 1)
    for i0 in range(0, N + 1, 400):
        i1 = min(i0 + 400, N + 1)
        d = np.abs(y[i0:i1, None] - y[None, :])
        d[np.arange(i1 - i0), np.arange(i0, i1)] = 1.0
        logd[i0:i1] = np.log(d).sum(1)
        P = np.log(np.abs(y[i0:i1, None] - s2[None, :]))
        sumlogP[i0:i1] = P.sum(1)
        logY[i0:i1] = np.log(np.abs(s2[i0:i1, None] - y[None, :])).sum(1)
        D = np.abs(s2[i0:i1, None] - s2[None, :])
        D[np.arange(i1 - i0), np.arange(i0, i1)] = 1.0
        logD[i0:i1] = np.log(D).sum(1)
    base = np.log(L * kap) - 2 * logd
    terms = np.zeros(N + 1)
    for i, g in enumerate(S):
        x = base + 2 * (sumlogP - np.log(np.abs(y - s2[i])))
        m = x.max()
        lognorm = m + math.log(np.exp(x - m).sum())
        terms[i] = 2 * (logY[i] - math.log(2 * g) - math.log(abs(math.sin(L * g))) - logD[i]) + lognorm
    m = terms.max()
    return m + math.log(np.exp(terms - m).sum()) - math.log(2), terms


def main():
    print('    X     N     K      b   nodes_capped      LB   worst node')
    for X, N in [(50, 391), (100, 921), (200, 1000), (200, 1500), (200, 2119), (1000, 600), (1000, 2000), (500, 6214), (1000, 6000), (1000, 13800)]:
        L = 0.5 * math.log(X)
        S, K, b, nover = nodes(N, L)
        v, t = lb_light(N, L, S)
        print(f"{X:6d} {N:5d} {K:5d} {b:8.1f} {nover:4d} {v:10.2f} {S[int(np.argmax(t))]:9.1f}", flush=True)


if __name__ == '__main__':
    main()
