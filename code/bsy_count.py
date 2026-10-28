# bsy_count.py: Theorem bsy:thm:thorincount of Chapter 35 (Balazard, Saias and Yor meet Thorin), "Counting with the Thorin density alone", group BSY-A, task C4 of
# misc/reviews/nick-b4/reconcile-M.md. Claims supported: the zeros of f are counted by the horizontal logarithmic derivative d_f = pi^(-1) Re(f'/f) alone,
#   sum_rho phi(rho) - sum_p phi(p) = -(1/2) int int d_f(u,t) K_phi(u,t) du dt,   K_phi = d_u phi + d_t^2 Phi,   Phi(u,t) = int_(-inf)^u phi(sigma,t) d sigma,
# on the box of the chapter: for zeta (no zero in the box) the functional returns 0; for the Davenport and Heilbronn function f it converges to phi(rho) = 0.1005312239 (M2:763 to
# M2:766: error 8.6e-6 at 512 nodes in t, deflated 1.3e-14); the deflation of Remark (iv): J_f = J_(fR) + phi(rho_hat) - phi(rho_hat - 1), R(s) = (s - rho_hat + 1)/(s - rho_hat), and the
# tail beyond u_max (Remark (iii): the integral of d_f over (u_max, infinity) is -pi^(-1) log|F(u_max + i t)|, F = f or f R).
# Normalisation: the Thorin density here is d_f = pi^(-1) Re f'/f at the level u, the measure without the factor 1/2 (Chapters 30 to 35, decisions.md).
# Conventions. Box 0.6 <= sigma <= 1.4, 83.7 <= t <= 87.7, phi(sigma, t) = b(x) b(y), x = (sigma - 1)/0.4, y = (t - 85.7)/2, b(x) = exp(-1/(1 - x^2)) on |x| < 1 (M2:618, 656).
# With B(x) = int_(-1)^x b: Phi = 0.4 B(x) b(y), K = b'(x) b(y)/0.4 + 0.1 B(x) b''(y), and for u >= 1.4 K = K_inf(t) = 0.1 B(1) b''(y), independent of u. f'/f is exact:
# Euler and Maclaurin for the Hurwitz zeta function and its derivative (bsy_zdh.py, direct partial sum plus the Bernoulli series, relative error 1e-13 against mpmath at 30
# digits, logged), f = zeta (pole at 1 outside the box) or the Davenport and Heilbronn f. The zero rho of f in the box is found by Newton's method on f (double precision).
# Method. t: periodic (midpoint) trapezoid rule with N_t = 128, 256, 512, 1024 nodes on [83.7, 87.7] (K vanishes with all derivatives at both ends); u: composite Gauss and
# Legendre, panels [0.6, 0.8, 1.0, 1.2, 1.4] with nu nodes each, then [1.4, 3] and [3, 6] (u_max = 6, 30, 45, 60) with nu/2 nodes each (nu = 24 and 48 in the second run); tail beyond u_max by the identity
# above. Variants: no deflation, deflation (Davenport and Heilbronn only) at rho_hat = rho + eps (1 + i)/sqrt 2 with eps = 0, 1e-6, 1e-4, 1e-2; for each, the tail left out, the tail of R only (the chapter's rule),
# and the full tail of f R. Error control: the u rule is run at two node counts and the t rule at four; the logged error is the difference to the target (0 for zeta, phi(rho) for f),
# and the u doubling difference separately; the Euler and Maclaurin values at N_t = 256 are recomputed with nf = 0.45, M = 36 (largest difference logged).
# Output: data/bsy_count.json, the log data/bsy_count.log. Usage: python3 -I bsy_count.py [--quick] [--out FILE]; --quick (N_t = 64, 128 and nu = 8, 16; about 3 s) writes to
# /tmp unless --out is given; the full case takes about a minute and fits the laptop limit only with --quick; on the cluster:
#   bash misc/tools/hopper_run.sh -c 4 -t 30 -g bsy_count.json -g bsy_count.log code/bsy_count.py code/bsy_zdh.py -- --out bsy_count.json
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse
import json
import math
import os
import sys
import time
import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_zdh as Z
import mpmath as mp

T0, T1, TC, TW = 83.7, 87.7, 85.7, 2.0
S0, S1, SC, SW = 0.6, 1.4, 1.0, 0.4
UMAXS = [6.0, 30.0, 45.0, 60.0]


def b(x):
    x = np.asarray(x, dtype=float)
    g = 1 - x * x
    out = np.zeros_like(x)
    m = g > 0
    out[m] = np.exp(-1 / g[m])
    return out


def bp(x):
    x = np.asarray(x, dtype=float)
    g = 1 - x * x
    out = np.zeros_like(x)
    m = g > 0
    out[m] = np.exp(-1 / g[m]) * (-2 * x[m] / g[m] ** 2)
    return out


def bpp(x):
    x = np.asarray(x, dtype=float)
    g = 1 - x * x
    out = np.zeros_like(x)
    m = g > 0
    xm, gm = x[m], g[m]
    out[m] = np.exp(-1 / gm) * (4 * xm ** 2 / gm ** 4 - 2 / gm ** 2 - 8 * xm ** 2 / gm ** 3)
    return out


_BC = {}


def _b1(z):
    return math.exp(-1 / (1 - z * z)) if abs(z) < 1 else 0.0


def Bint(x):
    out = []
    for xx in x:
        xx = float(xx)
        if xx not in _BC:
            _BC[xx] = quad(_b1, -1, min(xx, 1.0), epsabs=1e-17, epsrel=1e-14, limit=400)[0] if xx > -1 else 0.0
        out.append(_BC[xx])
    return np.array(out)


BETA = Bint(np.array([1.0]))[0]


def phi(s):
    return float(b(np.array([(s.real - SC) / SW]))[0] * b(np.array([(s.imag - TC) / TW]))[0])


def unodes(nu, umax):
    pts = [0.6, 0.8, 1.0, 1.2, 1.4] + [p for p in (2.0, 3.0, 4.5, 6.0, 9.0, 13.0, 20.0, 30.0, 45.0, 60.0) if p < umax] + [umax]
    ns = [nu] * 4 + [max(nu // 2, 8)] * (len(pts) - 5)
    xs, ws = [], []
    for a, c, n in zip(pts[:-1], pts[1:], ns):
        x, w = np.polynomial.legendre.leggauss(n)
        xs.append((c - a) / 2 * x + (c + a) / 2)
        ws.append((c - a) / 2 * w)
    return np.concatenate(xs), np.concatenate(ws)


def Kfun(u, t):
    x = (u - SC) / SW
    y = (t - TC) / TW
    Bx = Bint(x)
    return np.outer(bp(x), b(y)) / SW + 0.1 * np.outer(Bx, bpp(y)), 0.1 * BETA * bpp(y)


def dlog(kind, s, **kw):
    if kind == "zeta":
        z, d = Z.zeta(s, deriv=True, **kw)
        return d / z
    f, d = Z.dh(s, deriv=True, **kw)
    return d / f


def logabs(kind, s, **kw):
    return np.log(np.abs(Z.zeta(s, **kw) if kind == "zeta" else Z.dh(s, **kw)))


def zeta_plain(s, **kw):  # zeta without the factor s-1: f itself, pole at 1 is far from the box
    return Z.zeta(s, **kw)


def find_rho():
    s = np.array([0.8085171825 + 85.6993484854j])
    for _ in range(8):
        f, d = Z.dh(s, deriv=True)
        s = s - f / d
    return s[0]


def mp_dlog(kind, s, dps=30):
    mp.mp.dps = dps
    s = mp.mpmathify(s)
    if kind == "zeta":
        return complex(mp.zeta(s, 1, 1) / mp.zeta(s))
    k = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    c = {1: 1, 4: -1, 2: k, 3: -k}
    f = sum(c[r] * mp.zeta(s, mp.mpf(r) / 5) for r in c)
    df = sum(c[r] * mp.zeta(s, mp.mpf(r) / 5, 1) for r in c)
    return complex(df / f - mp.log(5))


def run(kind, rho, nts, nus, umax, log, kw=None):
    kw = kw or {}
    res = {}
    for nu in nus:
        u, wu = unodes(nu, umax)
        for nt in nts:
            h = (T1 - T0) / nt
            t = T0 + (np.arange(nt) + 0.5) * h
            K, Kinf = Kfun(u, t)
            S = u[:, None] + 1j * t[None, :]
            dl = dlog(kind, S.ravel(), **kw).reshape(S.shape)
            la_u = logabs(kind, umax + 1j * t, **kw)
            variants = {}
            for eps in ([None] + ([0.0, 1e-6, 1e-4, 1e-2] if kind == 'dh' else [])):
                if eps is None:
                    d = dl
                    J = lambda tailv: -0.5 * h * (((d.real / math.pi) * K * wu[:, None]).sum() + (Kinf * tailv).sum())
                    variants["none"] = {"no tail": J(0 * t), "full tail": J(-la_u / math.pi)}
                    continue
                rh = rho + eps * (1 + 1j) / math.sqrt(2)
                dR = 1 / (S - rh + 1) - 1 / (S - rh)
                d = dl + dR
                R_u = np.log(np.abs((umax + 1j * t - rh + 1) / (umax + 1j * t - rh)))
                J = lambda tailv: -0.5 * h * (((d.real / math.pi) * K * wu[:, None]).sum() + (Kinf * tailv).sum())
                add = phi(rh) - phi(rh - 1)
                variants[f"deflated eps={eps:g}"] = {"no tail": J(0 * t) + add, "R tail only": J(-R_u / math.pi) + add,
                                                     "full tail": J(-(la_u + R_u) / math.pi) + add}
            res[(nu, nt)] = variants
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    A = ap.parse_args()
    nts = [64, 128] if A.quick else [128, 256, 512, 1024]
    nus = [8, 16] if A.quick else [24, 48]
    outp = A.out or ("/tmp/bsy_count_quick.json" if A.quick else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "bsy_count.json"))
    lf = open(os.path.splitext(outp)[0] + ".log", "w")

    def log(m):
        print(m, flush=True)
        lf.write(m + "\n")
        lf.flush()

    t00 = time.time()
    log(f"bsy_count.py {'quick' if A.quick else 'full'}: N_t {nts}, nu {nus}, u_max {UMAXS}; beta = int b = {BETA:.15f}")
    for kind in ("zeta", "dh"):
        for s in (0.7 + 85.69j, 1.2 + 84.1j, 0.6 + 87.0j, 3 + 85j):
            v = dlog(kind, np.array([s]))[0]
            m = mp_dlog(kind, s)
            log(f"  f'/f vs mpmath (30 digits), {kind} at {s}: relative difference {abs(v - m) / abs(m):.1e}")
    rho = find_rho()
    f, d = Z.dh(np.array([rho]), deriv=True)
    tgt = phi(rho)
    log(f"rho = {rho.real:.12f} + {rho.imag:.12f} i (|f(rho)| = {abs(f[0]):.1e}); phi(rho) = {tgt:.10f}; 2 pi phi(rho) = {2 * math.pi * tgt:.8f}")
    out = {"rho": [rho.real, rho.imag], "phi_rho": tgt, "beta": BETA, "nts": nts, "nus": nus, "results": {}}
    for kind in ("zeta", "dh"):
        target = 0.0 if kind == "zeta" else tgt
        out["results"][kind] = {}
        for umax in UMAXS:
            res = run(kind, rho, nts, nus, umax, log)
            for nu in nus:
                for nt in nts:
                    out["results"][kind][f"umax{umax:g}_nu{nu}_nt{nt}"] = res[(nu, nt)]
            log(f"== {kind}, u_max = {umax:g}: target {target:.10f}; value (error against the target), nu = {nus[-1]}")
            names = list(res[(nus[0], nts[0])].keys())
            for name in names:
                for tv in res[(nus[-1], nts[0])][name].keys():
                    row = "; ".join(f"N_t={nt}: {res[(nus[-1], nt)][name][tv]:+.12f} ({res[(nus[-1], nt)][name][tv] - target:+.1e})" for nt in nts)
                    log(f"  {name:22s} {tv:12s} {row}")
            nt = nts[-1] if A.quick else 512
            log(f"  u rule doubling (nu {nus[0]} -> {nus[1]}) at N_t = {nt}, full tail, per variant: " + "; ".join(f"{name} {abs(res[(nus[1], nt)][name]['full tail'] - res[(nus[0], nt)][name]['full tail']):.1e}" for name in names))
    # second parameter set of the evaluator
    nt = nts[1]
    u, wu = unodes(nus[0], 6.0)
    t = T0 + (np.arange(nt) + 0.5) * (T1 - T0) / nt
    S = (u[:, None] + 1j * t[None, :]).ravel()
    for kind in ("zeta", "dh"):
        a1, a2 = dlog(kind, S), dlog(kind, S, nf=0.45, M=36)
        log(f"  evaluator, second parameters (nf 0.45, M 36), {kind}, {len(S)} points: max |d(f'/f)| = {np.abs(a1 - a2).max():.1e}")
    json.dump(out, open(outp, "w"), indent=1)
    log(f"done in {time.time() - t00:.0f} s; wrote {outp}")


if __name__ == "__main__":
    main()
