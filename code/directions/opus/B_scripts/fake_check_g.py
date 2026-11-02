#!/usr/bin/env python3
"""Checks of Theorem 1(a) and 1(g) of B.md. mpmath, floating point, NOT certified. Usage: python3 -I fake_check_g.py"""
import time
import mpmath as mp
mp.mp.dps = 60
T0 = time.time()

def closed(N):
    r = 4 * N; q = r * r; K = N * r
    y = [mp.cos(2 * mp.pi * j / q) for j in range(K)]
    W = lambda x: mp.fprod(x - yj for yj in y)
    W1 = W(mp.mpf(-1))
    c = [(-1) ** n * W(mp.cos(2 * mp.pi * n / q)) / W1 if min(n, q - n) >= K else mp.mpf(0) for n in range(q)]
    return r, q, K, c
def lagrange(N):
    r = 4 * N; q = r * r; K = N * r; half = q // 2
    y = [mp.cos(2 * mp.pi * n / q) for n in range(K, half + 1)]
    w = []
    for i, n in enumerate(range(K, half + 1)):
        p = mp.mpf(1)
        for j in range(len(y)):
            if j != i: p *= y[i] - y[j]
        w.append(1 / ((1 if n == half else 2) * p))
    return [w[i] / w[-1] for i in range(len(w))]
def ratio_formula(j, q, K):
    v = mp.cos(mp.pi * j / q)
    for i in range(j):
        v *= mp.cos(mp.pi * (K + i) / q) / mp.cos(mp.pi * (K - 1 - i) / q)
    return v

print("(a) closed form against the Lagrange weights, the ratio formula, the Gaussian bound")
for N in (2, 3, 4, 5, 6, 8):
    r, q, K, c = closed(N)
    half = q // 2
    lg = lagrange(N)
    d1 = max(abs(c[K + i] - lg[i]) for i in range(len(lg)))
    d2 = max(abs(abs(c[half - j]) - ratio_formula(j, q, K)) for j in range(K + 1))
    gm = max(abs(c[half - j]) * mp.exp(mp.pi * j * j / mp.mpf(q)) for j in range(K + 1))
    sign_ok = all(mp.sign(c[n]) == (-1) ** n for n in range(K, half + 1))
    print(f" N={N} q={q}: |closed-Lagrange|={mp.nstr(d1,3)}, |ratio formula error|={mp.nstr(d2,3)}, max_j |c_(q/2-j)| exp(pi j^2/q)={mp.nstr(gm,6)} (bound 2), signs {sign_ok}, {time.time()-T0:.0f}s")

def G_u(mu):
    u = mp.findroot(lambda u: u * (1 + u * u) / 2 - mu, mp.mpf(mu), tol=1e-40, verify=False)
    return (1 - u * u) / 2 - mu * (mp.pi / 2 - 2 * mp.atan(u)), u
print("G(mu):", {m: mp.nstr(G_u(mp.mpf(m))[0], 6) for m in (0.01, 0.05, 0.125, 0.25, 0.5, 0.75, 0.9, 0.99)})
print(f" limit G(0)=1/2, G(1)=0. Check of dG/dmu=-omega at mu=1/2: {mp.nstr(mp.diff(lambda x: G_u(x)[0], mp.mpf(0.5)),6)} against {mp.nstr(-(mp.pi/2-2*mp.atan(G_u(mp.mpf(0.5))[1])),6)}")
T = mp.mpf("3e12")
for mu in (mp.mpf(1) / 8, mp.mpf(1) / 2):
    NN = T / (8 * mp.pi * mu)
    print(f" T=3e12, mu={mp.nstr(mu,3)}: N={mp.nstr(mp.sqrt(NN),5)}, exponent 4 pi N^2 G(mu)={mp.nstr(4*mp.pi*NN*G_u(mu)[0],4)}")

def rL(s, r, q, c):
    tot = mp.mpf(0)
    for b in range(1, q + 1):
        v = c[b % q]
        if v != 0: tot += v * mp.zeta(s, mp.mpf(b) / q)
    return (mp.mpf(r) / q) ** s * tot
def M(s, r, q, c):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * rL(s, r, q, c)
def bound_g(s, N, q):
    t = s.imag; mu = t / (8 * mp.pi * N * N)
    G, u = G_u(mu); S = (1 - u * u) / (1 + u * u)
    return abs(s * (s - 1)) * (q + 3 + 2 / (5 * mp.pi * S)) / (mp.pi * N * N * S) * mp.exp(-mp.pi * t / 4 - 4 * mp.pi * N * N * G)
print("(g) |M| over the bound of Theorem 1(g), sigma=-1,1/2,2")
for N in (2, 3, 4):
    r, q, K, c = closed(N)
    for mu in (0.1, 0.5, 0.9):
        t = mp.mpf(mu) * 8 * mp.pi * N * N
        row = [mp.nstr(abs(M(mp.mpc(sg, t), r, q, c)) / bound_g(mp.mpc(sg, t), N, q), 3) for sg in (-1, 0.5, 2)]
        print(f" N={N} mu={mu} t={float(t):.2f}: {row} {time.time()-T0:.0f}s")
print(f"total {time.time()-T0:.0f}s")
