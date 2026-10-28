# bsy_levy.py: data of Figure bsy:fig:dhlevy of Chapter 35 (Balazard, Saias and Yor meet Thorin), the 2 sin^2 Levy density W_DH of the Davenport and Heilbronn function, group BSY-A,
# task C8 of misc/reviews/nick-b4/reconcile-M.md. Claims supported: M2:262 to M2:277 (W_DH computed by the zeros up to height 300 and by the explicit formula, agreement 6e-10
# relative for 0.1 <= u <= 0.5 and 2e-6 for u <= 2, off line terms never above 4 percent of W, W e^(gamma_1 u) tends to 1 from above) and the zero count M2:257 to M2:259 (204 sign
# changes of Lambda(1/2 + it) and 5 quadruples up to height 300, heights 85.70, 114.16, 166.48, 176.70, 240.40, count 204 + 2*5 = 214 by the argument principle). Data only: no figure.
# Normalisation: W is the trace of the chapter, W(u) = sum_(gamma' > 0) e^(-gamma' u) + 2 sum_(quadruples) e^(-gamma u) cos(delta u), u > 0, for Lambda(1/2 + x) with zeros
# +-i gamma' and quadruples +-delta +- i gamma (x = s - 1/2); equivalently the sum over the zeros of Lambda with positive imaginary part, e^(-gamma u) cos(delta u) (no factor 1/2).
# Conventions. f the Davenport and Heilbronn function (coefficients 1, kappa, -kappa, -1, 0 mod 5), Lambda(s) = (5/pi)^(s/2) Gamma((s+1)/2) f(s) = Lambda(1-s), real on Re s = 1/2,
# Lambda'/Lambda(s) = (1/2) log(5/pi) + (1/2) psi((s+1)/2) + f'/f(s). f and f' by Euler and Maclaurin (bsy_zdh.py, relative error 1e-13 against mpmath).
# Zeros. Line zeros: sign changes of Re(e^(iB) f(1/2 + it)) (B = Im of the logarithm of the gamma and power factors, so Lambda = e^A times this real number) on a grid of step 0.01
# in (0, 300], refined by 45 steps of bisection (error 1e-14). Zeros off the line: interior local minima of |f| on the grid sigma = 0.52(0.02)1.40, t step 0.05, refined by Newton's
# method on f; each is a zero of Lambda at 1/2 + delta + i gamma (and 1/2 - delta + i gamma by the functional equation). Completeness: the argument principle for Lambda on the rectangle
# [-1, 2] x [0, 300] (Lambda has no zero outside the strip -1 < Re s < 2, none on the real axis since f > 0 there; the left half uses Lambda(s) = Lambda(1 - s)), phase increments of
# Lambda(z_(k+1))/Lambda(z_k) on a polygon of step 0.005 and 0.0033, which must agree and be an integer; it must equal S + 2K, S line zeros and K quadruples. Explicit formula (derived
# here from the contour: the upper zeros of F(x) = Lambda(1/2 + x) are the poles of e^(izu) F'/F(z), the rectangle -alpha <= Re z <= alpha, Im z > 0, with F'/F odd and F(conj z) =
# conj F(z)): W(u) = (1/pi) Re[e^(i alpha u) A(u)] + (1/pi) int_0^alpha (Lambda'/Lambda)(1/2 + v) sin(v u) dv, A(u) = int_0^infinity (Lambda'/Lambda)(c + it) e^(-tu) dt, c = 3,
# alpha = c - 1/2. A is computed by composite Gauss and Legendre on [0, 450] (panels 0.5 with 20 nodes; second run 0.25 with 24), the second integral on [0, 0.5] and [0.5, 2.5].
# Output: data/bsy_levy.json (arrays on u = 0.1(0.005)2 and the tail curve W e^(gamma_1 u) on [0.1, 10], the zero lists, the maxima) and the log data/bsy_levy.log.
# Usage: python3 -I bsy_levy.py [--quick] [--out FILE]; --quick (heights up to 100, about 10 s) writes to /tmp unless --out is given; the full case on the cluster:
#   bash misc/tools/hopper_run.sh -c 4 -t 60 -g bsy_levy.json -g bsy_levy.log code/bsy_levy.py code/bsy_zdh.py -- --out bsy_levy.json
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse
import json
import math
import os
import sys
import time
import numpy as np
from scipy.special import loggamma, digamma

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_zdh as Z

LOG5PI = math.log(5 / math.pi)


def logfac(s):
    return 0.5 * s * LOG5PI + loggamma((s + 1) / 2)


def Lam(s, **kw):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    r = s.real < 0.5
    w = np.where(r, 1 - s, s)
    return np.exp(logfac(w)) * Z.dh(w, **kw)


def line_scan(T, step, kw=None):
    kw = kw or {}
    t = np.arange(step, T + step / 2, step)
    s = 0.5 + 1j * t
    v = np.exp(1j * logfac(s).imag) * Z.dh(s, **kw)
    return t, v


def line_val(t, kw=None):
    s = 0.5 + 1j * np.asarray(t)
    return (np.exp(1j * logfac(s).imag) * Z.dh(s, **(kw or {}))).real


def line_zeros(T, step=0.01):
    t, v = line_scan(T, step)
    imag = np.abs(v.imag).max() / np.abs(v).max()
    g = v.real
    idx = np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]
    a, b = t[idx].copy(), t[idx + 1].copy()
    ga = line_val(a)
    for _ in range(45):
        m = (a + b) / 2
        gm = line_val(m)
        left = np.sign(gm) == np.sign(ga)
        a = np.where(left, m, a)
        ga = np.where(left, gm, ga)
        b = np.where(left, b, m)
    return (a + b) / 2, imag


def offline_zeros(T):
    ts = np.arange(0.05, T, 0.05)
    sg = np.arange(0.52, 1.40 + 1e-9, 0.02)
    F = np.empty((len(ts), len(sg)))
    for j, sj in enumerate(sg):
        F[:, j] = np.abs(Z.dh(sj + 1j * ts))
    c = (F[1:-1, 1:-1] < F[:-2, 1:-1]) & (F[1:-1, 1:-1] < F[2:, 1:-1]) & (F[1:-1, 1:-1] < F[1:-1, :-2]) & (F[1:-1, 1:-1] < F[1:-1, 2:]) & (F[1:-1, 1:-1] < 0.5)
    ii, jj = np.nonzero(c)
    seeds = sg[jj + 1] + 1j * ts[ii + 1]
    found = []
    for z in seeds:
        s = np.array([z])
        for _ in range(30):
            f, d = Z.dh(s, deriv=True)
            s = s - f / d
        f, d = Z.dh(s, deriv=True)
        if abs(f[0]) < 1e-11 and 0.5 + 1e-6 < s[0].real < 1.45 and 0 < s[0].imag < T and not any(abs(s[0] - q) < 1e-6 for q in found):
            found.append(s[0])
    return sorted(found, key=lambda z: z.imag), len(seeds)


def count_rect(T, step, kw=None):
    kw = kw or {}
    n = int(round(T / step))
    t = np.linspace(0, T, n + 1)
    m = int(round(3 / step))
    sg = np.linspace(-1, 2, m + 1)
    pts = np.concatenate([2 + 1j * t, sg[::-1] + 1j * T, -1 + 1j * t[::-1]])  # the bottom edge (real axis, Lambda > 0) adds no phase; s = 1 is avoided
    L = Lam(pts, **kw)
    d = np.angle(L[1:] / L[:-1])
    d = np.append(d, np.angle(L[0] / L[-1]))
    return d.sum() / (2 * math.pi), float(np.abs(d).max())


def dlogLam(s, **kw):
    s = np.asarray(s, dtype=complex)
    return 0.5 * LOG5PI + 0.5 * digamma((s + 1) / 2) + Z.dh(s, deriv=True, **kw)[1] / Z.dh(s, **kw)


def explicit(us, panel, nn, kw=None):
    kw = kw or {}
    c, alpha, TA = 3.0, 2.5, 450.0
    x, w = np.polynomial.legendre.leggauss(nn)
    edges = np.arange(0, TA + 1e-9, panel)
    t = ((edges[1:, None] - edges[:-1, None]) / 2 * x[None, :] + (edges[1:, None] + edges[:-1, None]) / 2).ravel()
    wt = ((edges[1:, None] - edges[:-1, None]) / 2 * w[None, :]).ravel() * np.ones_like(t)
    Q = dlogLam(c + 1j * t, **kw)
    xs, ws = np.polynomial.legendre.leggauss(32)
    vs, vw = [], []
    for a_, b_ in ((0.0, 0.5), (0.5, alpha)):
        vs.append((b_ - a_) / 2 * xs + (b_ + a_) / 2)
        vw.append((b_ - a_) / 2 * ws)
    v = np.concatenate(vs)
    vw = np.concatenate(vw)
    R = dlogLam(0.5 + v, **kw)
    out = np.empty(len(us))
    for k, u in enumerate(us):
        A = (wt * Q * np.exp(-t * u)).sum()
        out[k] = (np.exp(1j * alpha * u) * A).real / math.pi + (vw * R.real * np.sin(v * u)).sum() / math.pi
    return out


def Wzeros(us, line, quads):
    u = np.asarray(us)[:, None]
    wl = np.exp(-line[None, :] * u).sum(1)
    if len(quads):
        g = np.array([z.imag for z in quads])[None, :]
        dl = np.array([z.real - 0.5 for z in quads])[None, :]
        off = 2 * (np.exp(-g * u) * np.cos(dl * u)).sum(1)
        offabs = 2 * (np.exp(-g * u) * np.abs(np.cos(dl * u))).sum(1)
    else:
        off = offabs = np.zeros(len(us))
    return wl + off, off, offabs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    A = ap.parse_args()
    T = 100.0 if A.quick else 300.0
    outp = A.out or ("/tmp/bsy_levy_quick.json" if A.quick else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "bsy_levy.json"))
    lf = open(os.path.splitext(outp)[0] + ".log", "w")

    def log(m):
        print(m, flush=True)
        lf.write(m + "\n")
        lf.flush()

    t00 = time.time()
    log(f"bsy_levy.py {'quick' if A.quick else 'full'}: height {T:g}")
    import mpmath as mp
    for s in (0.5 + 85.7j, 2.5 + 3j, 3 + 200j):
        m = Z.mp_h("dh", 0.5, s)
        log(f"  f vs mpmath (30 digits) at {s}: relative difference {abs(Z.dh(np.array([s]))[0] - m) / abs(m):.1e}")
    s1 = np.array([0.7 + 40j])
    l1 = np.exp(logfac(s1)) * Z.dh(s1)
    l2 = np.exp(logfac(1 - s1)) * Z.dh(1 - s1)
    log(f"  functional equation, direct evaluation at 0.7+40i and 0.3-40i: relative difference |Lambda(s) - Lambda(1-s)|/|Lambda(s)| = {abs(l1[0] - l2[0]) / abs(l1[0]):.1e}")
    line, imag = line_zeros(T)
    line2, _ = line_zeros(T, 0.0123)
    S = len(line)
    log(f"line zeros: S = {S} sign changes of Lambda(1/2 + it) on (0, {T:g}] (grid 0.01; grid 0.0123 gives {len(line2)}); max |Im|/max|Re| of e^(iB) f on the line {imag:.1e}; "
        f"first ordinate {line[0]:.10f} (chapter gamma_1 = 5.0942), last {line[-1]:.6f}")
    quads, nseeds = offline_zeros(T)
    K = len(quads)
    log(f"zeros off the line: {K} found from {nseeds} interior minima of |f| (sigma 0.52 to 1.40, step 0.02), by Newton:")
    for z in quads:
        log(f"   {z.real:.12f} + {z.imag:.10f} i  (delta = {z.real - 0.5:.6f})")
    cnt = []
    for st in (0.005, 0.0033):
        c, dm = count_rect(T, st)
        cnt.append(c)
        log(f"argument principle, Lambda on [-1, 2] x [0, {T:g}], polygon step {st}: count {c:.6f} (max |increment| {dm:.2f})")
    log(f"  S + 2K = {S} + 2*{K} = {S + 2 * K}; counts agree with it: {all(abs(c - (S + 2 * K)) < 0.01 for c in cnt)}")
    us = np.round(np.arange(0.1, 2.0 + 1e-9, 0.005), 6)
    Wz, off, offabs = Wzeros(us, line, [q for q in quads])
    t1 = time.time()
    We = explicit(us, 0.5, 20)
    We2 = explicit(us, 0.25, 24, {"nf": 0.45, "M": 36})
    log(f"explicit formula on {len(us)} values of u in [0.1, 2] ({time.time() - t1:.0f} s); second run (panels 0.25, 24 nodes, nf 0.45, M 36): max |difference| {np.abs(We - We2).max():.1e}, "
        f"max relative {np.abs((We - We2) / We).max():.1e}")
    rel = np.abs(We - Wz) / Wz
    for lo, hi in ((0.1, 0.5), (0.1, 2.0)):
        m = (us >= lo - 1e-9) & (us <= hi + 1e-9)
        log(f"  zeros against explicit formula, 0.1 <= u <= {hi}: max relative difference {rel[m].max():.2e} at u = {us[m][rel[m].argmax()]:.3f}; max absolute {np.abs(We - Wz)[m].max():.1e}")
    log(f"  W by zeros: min over the grid {Wz.min():.3e} (u = {us[Wz.argmin()]:.3f}), W(0.1) = {Wz[0]:.6f}, W(1) = {Wz[np.argmin(np.abs(us - 1))]:.6f}, W(2) = {Wz[-1]:.6e}; all positive: {bool((Wz > 0).all())}")
    ufine = np.round(np.arange(0.1, 2.0 + 1e-9, 0.001), 6)
    Wf, of, oa = Wzeros(ufine, line, quads)
    sh, sa = np.abs(of) / Wf, oa / Wf
    log(f"  off line share on 0.1 <= u <= 2 (step 0.001): max |off line sum|/W = {sh.max():.4f} at u = {ufine[sh.argmax()]:.3f}; max (sum of |terms|)/W = {sa.max():.4f} at u = {ufine[sa.argmax()]:.3f}")
    ut = np.round(np.arange(0.1, 10.0 + 1e-9, 0.05), 6)
    Wt, ot, oat = Wzeros(ut, line, quads)
    g1 = line[0]
    ratio = Wt * np.exp(g1 * ut)
    m5 = ut <= 5.0
    log(f"  W e^(gamma_1 u) - 1 on [0.1, 5]: min {(ratio[m5] - 1).min():.3e}, max {(ratio[m5] - 1).max():.3e}; above 1 on [0.1, 5]: {bool((ratio[m5] > 1).all())}; at u = 10: {ratio[-1] - 1:.3e} (rounding level 1e-15)")
    uw = np.round(np.arange(0.01, 2.0 + 1e-9, 0.001), 6)
    Ww, ow, oaw = Wzeros(uw, line, quads)
    log(f"  for information, u from 0.01 (the zero sum cut at height {T:g} is not converged below u = 0.1): max |off line sum|/W = {(np.abs(ow) / Ww).max():.4f} at u = {uw[(np.abs(ow) / Ww).argmax():][0]:.3f}; "
        f"max (sum of |terms|)/W = {(oaw / Ww).max():.4f} at u = {uw[(oaw / Ww).argmax()]:.3f}; share at u = 0.05: {abs(ow[np.argmin(np.abs(uw - 0.05))]) / Ww[np.argmin(np.abs(uw - 0.05))]:.4f}")
    # off-line terms relative to W on the explicit-formula grid, first quadruple alone
    json.dump({"height": T, "S": S, "K": K, "line_zeros": [float(x) for x in line], "offline_zeros": [[z.real, z.imag] for z in quads], "rect_counts": cnt,
               "u": us.tolist(), "W_zeros": Wz.tolist(), "W_explicit": We.tolist(), "W_explicit_second": We2.tolist(), "offline_sum": off.tolist(), "offline_abs": offabs.tolist(),
               "max_share_signed": float(sh.max()), "u_at_max_signed": float(ufine[sh.argmax()]), "max_share_abs": float(sa.max()), "u_at_max_abs": float(ufine[sa.argmax()]),
               "tail_u": ut.tolist(), "tail_ratio": ratio.tolist()}, open(outp, "w"))
    log(f"done in {time.time() - t00:.0f} s; wrote {outp}")


if __name__ == "__main__":
    main()
