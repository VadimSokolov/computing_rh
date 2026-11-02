#!/usr/bin/env python3
"""Exact fakes of xi (Theorem 1 of B.md): F_eps = xi + eps*M with the same primes below N.

Floating point only, nothing here is certified. Usage: python3 -I fake_exact.py N Tmax dt [Tscan]
Writes fake_exact_N<N>.json next to this script.
"""
import json, os, sys, time
import numpy as np
from math import pi

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
Tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
dt = float(sys.argv[3]) if len(sys.argv) > 3 else 0.02
Tscan = float(sys.argv[4]) if len(sys.argv) > 4 else max(1.0, Tmax - 40.0)
r = 4 * N; q = r * r; K = N * r            # K = q/4 constraints, eigenspace dimension q/4 + 1
idx = np.arange(q)
C = np.cos(2 * pi * np.outer(idx, idx) / q) / np.sqrt(q)   # the DFT on even vectors, 0 on odd vectors
# Exact vector, Theorem 1(a) of B.md: the Lagrange weights at the nodes cos(2 pi n/q), K <= n <= q/2.
# A floating point null vector is NOT usable: the K conditions have many singular values below rounding
# (null_dim below), so an SVD returns an arbitrary mix of near solutions.
import mpmath as mp
mp.mp.dps = 120
half = q // 2
nodes = [mp.cos(2 * mp.pi * n / q) for n in range(K, half + 1)]
wts = []
for i, n in enumerate(range(K, half + 1)):
    p = mp.mpf(1)
    for j in range(len(nodes)):
        if j != i: p *= nodes[i] - nodes[j]
    wts.append(1 / ((1 if n == half else 2) * p))
c = np.zeros(q)
for i, n in enumerate(range(K, half + 1)):
    c[n] = float(wts[i] / wts[-1]); c[(q - n) % q] = c[n]
w, V = np.linalg.eigh(C)
B = V[:, w > 0.5]
mplus = B.shape[1]
S = np.linalg.svd(B[:K, :], compute_uv=False)
null_dim = mplus - int(np.sum(S > 1e-10))
resid = np.abs(C @ c - c).max()
print(f"N={N} r={r} q={q} K={K} dim(+1)={mplus} null_dim={null_dim} sv_min={S[-1]:.2e} |Fc-c|={resid:.2e}")
neff = {th: int(np.argmax(np.abs(c) > th)) for th in (1e-6, 1e-3, 1e-1)}
print("first nonzero c_n:", [(int(n), float(f"{c[n]:.2e}")) for n in range(K, K + 4)], "first n with |c_n|>th:", neff, "lambda_eff(1e-3)=", neff[1e-3] / r, "2 pi lambda_eff^2=", round(2 * pi * (neff[1e-3] / r) ** 2, 1), "2 pi N^2=", round(2 * pi * N * N, 1))

# Hurwitz zeta at a=b/q, b=1..q, by Euler and Maclaurin, for one complex s
from math import factorial
BERN = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6]
A = (idx + 1) / q                            # a = 1/q, ..., 1
cb = np.roll(c, -1)                          # c_b for b = 1..q (index b-1); c_q = c_0
def hurwitz_vec(s, Mterms):
    k = np.arange(Mterms)[None, :]
    lg = np.log(A[:, None] + k)
    head = np.exp(-s * lg).sum(axis=1)
    X = A + Mterms
    lX = np.log(X)
    tail = np.exp((1 - s) * lX) / (s - 1) + 0.5 * np.exp(-s * lX)
    poch = s
    for j, b in enumerate(BERN, start=1):
        tail = tail + b / factorial(2 * j) * poch * np.exp((-s - 2 * j + 1) * lX)
        poch = poch * (s + 2 * j - 1) * (s + 2 * j)
    return head + tail
def zeta_and_M(s):
    Hz = hurwitz_vec(s, int(abs(s.imag)) + 40)
    qs = np.exp(-s * np.log(q))
    return qs * Hz.sum(), np.exp(s * np.log(r)) * qs * (cb * Hz).sum()

t0 = time.time()
# 1. reality of M/zeta on the critical line (functional equation test) and collision thresholds
ts = np.arange(1.0, Tmax, dt)
ratio = np.empty(ts.size, dtype=complex)
for i, t in enumerate(ts):
    z, m = zeta_and_M(complex(0.5, t))
    ratio[i] = z / m
imag_rel = float(np.median(np.abs(ratio.imag) / (np.abs(ratio.real) + 1e-300)))
Rr = -ratio.real                             # zeros of xi+eps*M on the line are the solutions of Rr(t)=eps
d = np.diff(Rr)
ext = []
for i in range(1, d.size):
    if d[i - 1] * d[i] < 0 and abs(Rr[i]) < 1e3 and abs(d[i]) < 0.5 and abs(d[i-1]) < 0.5:
        # parabolic refinement of the local extremum
        a, b_, c_ = Rr[i - 1], Rr[i], Rr[i + 1]
        den = a - 2 * b_ + c_
        tv = ts[i] + 0.5 * dt * (a - c_) / den
        val = b_ - 0.125 * (a - c_) ** 2 / den
        ext.append((float(tv), float(val)))
ext.sort(key=lambda e: e[0])
# running minimum of |extremum|: the smallest eps > 0 or < 0 at which two line zeros can collide below t
pos = [(t, v) for t, v in ext]
print(f"median |Im ratio|/|Re ratio| on the line = {imag_rel:.2e}; local extrema found: {len(ext)}; {time.time()-t0:.1f}s")
absM = np.abs(1.0 / ratio)
blocks = [(float(ts[i]), float(absM[i:i + int(5 / dt)].max())) for i in range(0, ts.size, int(5 / dt))]
print("max |M/zeta| on the line in blocks of length 5:", [(round(a, 0), float(f"{b:.1e}")) for a, b in blocks][:: max(1, len(blocks) // 24)])
rec = []; best = np.inf
for t, v in pos:
    if abs(v) < best:
        best = abs(v); rec.append((round(t, 4), float(f"{v:.5e}")))
print("record small |local extremum of -xi/M| (height, value):", rec)

# 2. data sensitivity: C_W = L2 norm over v in [0,W] of d/d eps of u_1^{F_eps} at eps=0
out = {"neff": {str(k): v for k, v in neff.items()}, "blocks": blocks, "N": N, "r": r, "q": q, "K": K, "dim_plus": mplus, "null_dim": null_dim, "imag_rel": imag_rel, "records": rec}
h = 1e-4
for W in (50.0, 120.0):
    vs = np.arange(0.0, W + 0.05, 0.1)
    g = np.empty(vs.size)
    for i, v in enumerate(vs):
        zp, mp_ = zeta_and_M(complex(1 + h, v)); zm, mm_ = zeta_and_M(complex(1 - h, v)) if v > 0.5 else zeta_and_M(complex(1 + 2 * h, v))
        if v > 0.5:
            g[i] = ((mp_ / zp - mm_ / zm) / (2 * h)).real / pi
        else:
            g[i] = ((mm_ / zm - mp_ / zp) / h).real / pi
    CW = float(np.sqrt(np.trapezoid(g ** 2, vs))); SW = float(np.abs(g).max())
    print(f"window [0,{W:.0f}]: L2 sensitivity C_W={CW:.3e}, sup={SW:.3e}")
    out[f"C_W_{int(W)}"] = CW; out[f"sup_W_{int(W)}"] = SW

# 3. zeros with Re s > 1/2 of zeta + eps*M: scan |zeta/M + eps| for a list of eps on a coarse grid
sig = np.arange(0.55, 1.6, 0.05); tg = np.arange(Tscan, Tmax, 0.1)
Z = np.empty((sig.size, tg.size), dtype=complex); Mv = Z.copy()
for i, sg in enumerate(sig):
    for j, t in enumerate(tg):
        Z[i, j], Mv[i, j] = zeta_and_M(complex(sg, t))
found = {}
for eps in (10.0, 1.0, 0.1, 0.01, -0.01, -0.1, -1.0, -10.0):
    f = Z + eps * Mv
    zs = []
    for i in range(sig.size - 1):
        for j in range(tg.size - 1):
            # argument principle on the cell, by the change of arg of f along its boundary
            corners = [f[i, j], f[i + 1, j], f[i + 1, j + 1], f[i, j + 1], f[i, j]]
            turn = sum(np.angle(corners[k + 1] / corners[k]) for k in range(4))
            if abs(turn) > 3.0:
                s = complex(sig[i] + 0.025, tg[j] + 0.05)
                for _ in range(30):
                    z1, m1 = zeta_and_M(s); val = z1 + eps * m1
                    z2, m2 = zeta_and_M(s + 1e-5); der = ((z2 + eps * m2) - val) / 1e-5
                    if der == 0: break
                    s = s - val / der
                z1, m1 = zeta_and_M(s)
                if abs(z1 + eps * m1) < 1e-8 and 0.51 < s.real < 2.0 and 0 < s.imag < Tmax + 20:
                    if not any(abs(s - complex(*w_)) < 1e-6 for w_ in zs):
                        zs.append((round(s.real, 8), round(s.imag, 8)))
    found[str(eps)] = sorted(zs, key=lambda w_: w_[1])
    print(f"eps={eps:+.2f}: zeros with Re s>0.51, t<{tg[-1]:.0f}: {found[str(eps)][:8]}")
out["offline_scan"] = found
# 4. functional equation test by mp-free direct formula is done through the reality test above
out["c"] = [float(x_) for x_ in c[K:half + 1]]
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"fake_exact_N{N}.json"), "w"), indent=0)
print(f"total {time.time()-t0:.1f}s")
