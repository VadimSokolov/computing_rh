#!/usr/bin/env python3
"""The local sign test of Proposition 5 of B.md: eta < 0 line doublet, eta > 0 off line pair.

Floating point, not certified; the noise is a modelling device. Usage: python3 -I sign_test.py
Part A: real u_1 (mpmath) on a window around the Lehmer pair near 7005, Gauss and Newton fit.
Part B: design only, se(eta)/tau for all adjacent pairs in windows of real zeros, up to 3e12.
"""
import json, os, sys, time
import numpy as np
import mpmath as mp
from math import pi, erf, sqrt

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", "..", "..", "..", ".."))
Phi = lambda z: 0.5 * (1 + erf(z / sqrt(2)))
A0 = 0.5     # level theta = 1

def psi(a, eta, x):        # doublet at distance x, level offset a, parameter eta (real part of 2(a+ix)/((a+ix)^2-eta) / pi)
    w = a + 1j * x
    return (2 * w / (w * w - eta)).real / pi
def dpsi_eta(a, eta, x):
    w = a + 1j * x
    return (2 * w / (w * w - eta) ** 2).real / pi
def dpsi_x(a, eta, x):     # derivative in x; d/dg is minus this
    w = a + 1j * x
    return (2j * (-(w * w) - eta) / (w * w - eta) ** 2).real / pi
def P1(a, x):
    return a / (a * a + x * x) / pi
def dP1(a, x):             # derivative of P_a(x) in x
    return -2 * a * x / (a * a + x * x) ** 2 / pi

def design(vs, zeros, k, a=A0, deg=5):
    """Jacobian columns for the doublet (k,k+1) at the current positions: returns t, T' (columns)."""
    g = 0.5 * (zeros[k] + zeros[k + 1]); s = zeros[k + 1] - zeros[k]; eta = -s * s / 4
    x = vs - g
    t = dpsi_eta(a, eta, x)
    cols = [-dpsi_x(a, eta, x)]
    for m, z in enumerate(zeros):
        if m in (k, k + 1): continue
        cols.append(-dP1(a, vs - z))
    w = (vs - vs.mean()) / (0.5 * (vs.max() - vs.min()))
    for d in range(deg + 1):
        cols.append(np.polynomial.chebyshev.chebval(w, [0] * d + [1]))
    return t, np.column_stack(cols), g, eta

def tperp(t, Tp):
    scale = np.linalg.norm(Tp, axis=0); scale[scale == 0] = 1
    Q, _ = np.linalg.qr(Tp / scale)
    r = t - Q @ (Q.T @ t)
    r = r - Q @ (Q.T @ r)      # reorthogonalise
    return np.linalg.norm(r), np.linalg.norm(t)

def part_b(name, zeros_all, c, half, dv=0.02, nwin=None):
    zs = np.array([z for z in zeros_all if abs(z - c) <= half + 4.0])
    vs = np.arange(c - half, c + half + dv / 2, dv)
    inner = [k for k in range(len(zs) - 1) if abs(0.5 * (zs[k] + zs[k + 1]) - c) <= half * 0.5]
    rows = []
    for k in inner:
        t, Tp, g, eta = design(vs, zs, k)
        tp, tn = tperp(t, Tp)
        rows.append((g, zs[k + 1] - zs[k], tp, tn))
    rows = np.array(rows)
    s = rows[:, 1]; tp = rows[:, 2]
    # z at tau per sample: |eta_0| / se = (s^2/4) * ||t_perp|| / tau
    tau_for_z4 = (s * s / 4) * tp / 4
    i = int(np.argmin(tau_for_z4))
    print(f"{name}: height {c:.0f}, {len(zs)} zeros, {vs.size} samples, {len(inner)} adjacent pairs: mean gap {s.mean():.4f}, min gap {s.min():.5f}")
    print(f"   ||t|| median {np.median(rows[:,3]):.2f}, ||t_perp|| median {np.median(tp):.3e}, min {tp.min():.3e}; tau needed for z=4: median {np.median(tau_for_z4):.2e}, worst {tau_for_z4.min():.2e} (gap {s[i]:.5f} at {rows[i,0]:.3f})")
    return dict(name=name, height=c, n=len(zs), mean_gap=float(s.mean()), min_gap=float(s.min()), tperp_median=float(np.median(tp)), tperp_min=float(tp.min()), tau_z4_median=float(np.median(tau_for_z4)), tau_z4_worst=float(tau_for_z4.min()))

t0 = time.time()
Z = np.array([float(l) for l in open(os.path.join(root, "data", "zeros_8000.txt"))])
res = {}
# Part A: Lehmer pair
k0 = int(np.argmin(np.diff(Z[Z > 7000])[:10])) + int(np.sum(Z <= 7000))
g0, g1 = Z[k0], Z[k0 + 1]; c = 0.5 * (g0 + g1); s0 = g1 - g0
print(f"pair {g0:.6f} {g1:.6f} gap {s0:.6f} eta0=-s^2/4={-s0*s0/4:.4e}")
half = 6.0; dv = 0.05
vs = np.arange(c - half, c + half + dv / 2, dv)
mp.mp.dps = 30
def u1(v):
    s = mp.mpc(1, v)
    val = 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)
    return float(val.real / mp.pi)
cache = os.path.join(here, "lehmer_u1.npy")
if os.path.exists(cache):
    data = np.load(cache)
else:
    data = np.array([u1(v) for v in vs]); np.save(cache, data)
print(f"real u_1 on {vs.size} points in {time.time()-t0:.1f}s")
zs = Z[np.abs(Z - c) <= half + 4.0].copy()
k = int(np.argmin(np.abs(zs - g0)))
def fit(y, eta_start, others, g, deg=5, iters=40):
    """Gauss and Newton in (eta, g, free ordinates, background); ordinates beyond the data window stay fixed."""
    free = np.abs(others - c) <= half - 1.0
    pars_z = others.copy(); eta = eta_start
    w = (vs - vs.mean()) / half
    bg = [np.polynomial.chebyshev.chebval(w, [0] * d + [1]) for d in range(deg + 1)]
    for it in range(iters):
        x = vs - g
        model = psi(A0, eta, x) + sum(P1(A0, vs - z) for z in pars_z)
        cols = [dpsi_eta(A0, eta, x), -dpsi_x(A0, eta, x)] + [-dP1(A0, vs - z) for z in pars_z[free]] + bg
        J = np.column_stack(cols)
        scale = np.linalg.norm(J, axis=0)
        sol, *_ = np.linalg.lstsq(J / scale, y - model, rcond=1e-14)
        sol = sol / scale
        nf = int(free.sum())
        lim = max(np.abs(sol[1:2 + nf]).max() / 0.05, abs(sol[0]) / 2e-3, 1.0)   # step limit
        sol = sol / lim
        eta += sol[0]; g += sol[1]; pars_z[free] = pars_z[free] + sol[2:2 + nf]
        if lim == 1.0 and abs(sol[0]) < 1e-16: break
    t, Tp = J[:, 0], J[:, 1:]
    tp, _ = tperp(t, Tp)
    return eta, g, tp
others = np.delete(zs, [k, k + 1])
rng = np.random.default_rng(1)
xg = vs - c
rows = []
for label, y0 in (("zeta", data), ("alt: same doublet moved off line, eta=+s^2/4", data - psi(A0, -s0 * s0 / 4, xg) + psi(A0, +s0 * s0 / 4, xg)), ("alt: eta=+s^2/400", data - psi(A0, -s0 * s0 / 4, xg) + psi(A0, s0 * s0 / 400, xg))):
    for tau in (1e-6, 1e-8, 1e-10, 1e-12):
        y = y0 + tau * rng.standard_normal(vs.size)
        eta_hat, g_hat, tp = fit(y, -1e-4, others, c)
        se = tau / tp
        rows.append((label, tau, eta_hat, se, Phi(eta_hat / se)))
        print(f"  {label[:44]:44s} tau={tau:.0e}: eta_hat={eta_hat:+.4e} se={se:.2e} P(eta>0|Y, flat)={Phi(eta_hat/se):.4f}")
res["lehmer"] = dict(g0=g0, g1=g1, gap=s0, rows=rows)
# Part B: design only
res["design"] = []
res["design"].append(part_b("zeros_8000", Z, 1000.0, 15.0))
res["design"].append(part_b("zeros_8000", Z, 7005.0, 15.0))
for fn, label in (("ex_zeros3e12.json", "3e12"), ("es_zeros1e13.json", "1e13")):
    d = json.load(open(os.path.join(root, "data", fn)))
    from decimal import Decimal
    cc = float(d["center"])
    zz_rel = np.array(sorted(float(Decimal(i["ordinate"]) - Decimal(str(d["center"]))) for i in d["zeros"]))
    print(f"{fn}: center {cc}, half window {d['half_window']}, {zz_rel.size} zeros, first offsets {[round(float(x), 3) for x in zz_rel[:3]]}")
    half_b = float(d["half_window"]) - 2.0
    r = part_b(label, zz_rel, 0.0, half_b)
    r["height"] = cc
    res["design"].append(r)
json.dump(res, open(os.path.join(here, "sign_test.json"), "w"), indent=0, default=float)
print(f"total {time.time()-t0:.1f}s")
