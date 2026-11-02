#!/usr/bin/env python3
"""Cluster fake (Theorem 2 of B.md): replace the zeros in [H, H+ell] by the roots of P + b.

Multiprecision floating point (mpmath), nothing here is certified.
Usage: python3 -I cluster_fake.py H T xi   (T = top of the data window, xi = relative excess of b)
Reads the table data/zeros_8000.txt of the book (ordinates up to 8148).
"""
import json, os, sys, time
from math import e as E, log, pi, exp
import numpy as np
import mpmath as mp

H = float(sys.argv[1]) if len(sys.argv) > 1 else 500.0
T = float(sys.argv[2]) if len(sys.argv) > 2 else H / 2
xi = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", "..", "..", "..", ".."))
zeros = [line.strip() for line in open(os.path.join(root, "data", "zeros_8000.txt"))]
ell = (H - T) / E
cl = [z for z in zeros if H <= float(z) <= H + ell]
n = len(cl)
mp.mp.dps = int(1.0 * n) + 120
gam = [mp.mpf(z) for z in cl]
c = (gam[0] + gam[-1]) / 2
tk = [g - c for g in gam]
tf = np.array([float(t) for t in tk])
print(f"H={H} T={T} ell={ell:.3f} n={n} first={float(gam[0]):.4f} last={float(gam[-1]):.4f} dps={mp.mp.dps}")

def logabsP(x):   # log |P(x)| at real x (double)
    return float(np.log(np.abs(x - tf)).sum())
ext = []
for j in range(n - 1):
    a, b = tf[j], tf[j + 1]
    lo, hi = a + 1e-12, b - 1e-12
    for _ in range(80):        # bisection on P'/P, decreasing on each gap
        m = 0.5 * (lo + hi)
        if (1.0 / (m - tf)).sum() > 0: lo = m
        else: hi = m
    x = 0.5 * (lo + hi)
    ext.append((logabsP(x), x, j))
ext.sort()
lmin, xstar, jstar = ext[0]
print(f"log e* = {lmin:.3f} at Y={xstar:.5f} (height {float(c) + xstar:.5f}), gap {float(gam[jstar+1]-gam[jstar]):.5f}; next log extremum {ext[1][0]:.3f}")

def P(t):
    r = mp.mpf(1)
    for z in tk: r = r * (t - z)
    return r
xs = mp.mpf(xstar)
for _ in range(6):             # refine the critical point of P in multiprecision
    d1 = sum(1 / (xs - z) for z in tk); d2 = -sum(1 / (xs - z) ** 2 for z in tk)
    xs = xs - d1 / d2
Pxs = P(xs)
b = -Pxs * (1 + mp.mpf(xi))
print(f"|b|=exp({float(mp.log(abs(b))):.3f}), n log(ell)={n*log(ell):.1f}, n log(ell/2)={n*log(ell/2):.1f}")
# the off line root near xs: start from the second order estimate
d2 = -sum(1 / (xs - z) ** 2 for z in tk)
delta0 = mp.sqrt(2 * mp.mpf(xi) / abs(d2))
t = mp.mpc(xs, delta0)
for _ in range(60):
    Pt = P(t); dlog = sum(1 / (t - z) for z in tk)
    t = t - (1 + b / Pt) / (dlog * 1) * 1 if False else t - (Pt + b) / (Pt * dlog)
print(f"off line root: ordinate {float(c + t.real):.8f}, delta = {float(t.imag):.8f}, residual |P+b|/|b| = {float(abs(P(t)+b)/abs(b)):.2e}")
# data difference on the window, at level theta, for heights 0..T
res = {"H": H, "T": T, "ell": ell, "n": n, "log_estar": lmin, "ordinate": float(c + t.real), "delta": float(t.imag)}
out = []
for theta in (1.0, 0.75, 0.6):
    vs = np.linspace(0, T, 41)
    mx = mp.mpf(0); arg = 0
    for v in vs:
        tf_ = mp.mpc(v, -(theta - 0.5))      # t_full = -i (s - 1/2) = v - i (theta - 1/2)
        tot = mp.mpc(0)
        for sign in (1, -1):
            Y = sign * tf_ - c
            Pv = P(Y); dlog = sum(1 / (Y - z) for z in tk)
            # d/ds log(1 + b/P(Y)) with Y = +-t_full - c, d/ds = -i d/dt_full
            tot += -1j * sign * (-(b / (Pv + b)) * dlog)
        du = abs(tot.real) / mp.pi
        if du > mx: mx, arg = du, v
    bound = 4 * n / (H - T) * mp.e ** (-n) / (1 - 2 * mp.e ** (-n)) / mp.pi
    print(f"level {theta}: sup over the grid of |Delta u| = {mp.nstr(mx, 4)} at v={arg:.1f}; bound of Theorem 2 = {mp.nstr(bound, 4)}; log10 sup = {float(mp.log10(mx)):.1f}, bits to see it = {float(-mp.log(mx, 2)):.0f}")
    out.append((theta, float(mp.log10(mx))))
res["log10_sup"] = out

# the prime side: D_tot(x) = sum over the moved zeros of the change of 2 Re e^{i t x}, by a contour integral
def D_tot(x, A=0.75):
    R = A * ell
    x = mp.mpf(x)
    M = int(2 * n + 2 * R * abs(x) + 200)
    s = mp.mpc(0)
    for k in range(M):
        z = R * mp.expjpi(mp.mpf(2 * k + 1) / M)
        s += mp.exp(1j * z * x) * mp.log(1 + b / P(z)) * z
    # D(x) = (1/2 pi i) oint e^{itx} d log(1+b/P) = -(x/2 pi) oint e^{itx} log(1+b/P) dt, dt = i z dphi
    Dx = -(x / (2 * mp.pi)) * s * 1j * (2 * mp.pi / M)
    return 2 * (mp.expj(c * x) * Dx).real
rho = log(H / (2 * pi)) / (2 * pi)
print(f"density rho(H) = {rho:.4f}, Nyquist point pi rho = {pi*rho:.3f}, e^(pi rho) = {exp(pi*rho):.2f} = (H/2 pi)^(1/2)")
tab = []
t0 = time.time()
mp.mp.dps = int(1.0 * n) + 120 + int(0.75 * ell * 6 / 2.3)
for X in (2, 3, 5, 7, 11, 13, 17, 23, 31, 47, 71, 101, 151, 251):
    x = log(X)
    val = D_tot(x)
    tab.append((X, float(mp.log10(abs(val)))))
    print(f"X={X:4d} x={x:.3f}  log10|D_tot(x)| = {tab[-1][1]:.1f}   ({time.time()-t0:.0f}s)")
res["D_tot"] = tab
json.dump(res, open(os.path.join(here, f"cluster_fake_H{int(H)}.json"), "w"), indent=0)
