#!/usr/bin/env python3
"""D10: floating point checks of the lower bound (Theorem G of D10.md): sampling of every lemma, margins, the chain
of the proof on random functions, and the grid behind the conjecture kappa_isolated = kappa_1, now with lower pairs.
Nothing here is certified.  Usage: timeout 590 python3 -I D10_lower.py [outer_scale]
"""
import sys
import numpy as np

pi = np.pi
rng = np.random.default_rng(11)
T_C = 2.798386045783887
RC = np.cos(pi - T_C)
K1 = 0.3927607377918873


def kap(d, y, h):
    return np.arctan((y + h) / d) + np.arctan((y - h) / d)


def kc(t, Y, h):  # args in (0,pi), Y>h
    return np.arctan2(Y - h, t) + np.arctan2(Y + h, t)


def Kf(y, d):
    return y / (d * d + y * y)


def slope(d, w, h):
    return Kf(w + h, d) + Kf(w - h, d)


def sigma(t, h):
    return Kf(1 + h, t) + Kf(1 - h, t)


def E0(d, w, h):
    return kap(d, w, h) - w * kap(d, 1.0, h)


def Hb(w):
    return (1 - w ** 2) * (1 + 5 * w ** 2) / (4 * (1 - 3 * w ** 2))


def B1(w):
    return w * (1 + 5 * w ** 2) / (2 * (1 - 3 * w ** 2))


def Afun(u, w):
    return np.angle((u + 1j * w) ** 2 + 1)


def Ax(u, w):
    z = u + 1j * w
    return (2 * z / (z * z + 1)).imag


def lemmas(N=4000000):
    print("lemma checks by sampling, N =", N, "(max violation, positive would be a failure)")
    w = rng.uniform(1e-3, 0.3, N)
    d = np.exp(rng.uniform(np.log(1e-3), np.log(1e3), N))
    h = w * rng.uniform(0, 1, N)
    print(" A  E0(h)-E0(0), lower pairs   ", np.max(E0(d, w, h) - E0(d, w, 0 * h)))
    print(" B  E0(0)/(2K_w)-(1-w^2)(1+w^2/d^2)/4 ", np.max(E0(d, w, 0 * h) / (2 * Kf(w, d)) - (1 - w ** 2) * (1 + w ** 2 / d ** 2) / 4))
    print(" C  2w/(d^2+4w^2)-s            ", np.max(2 * w / (d ** 2 + 4 * w ** 2) - slope(d, w, h)))
    print(" E  -E0, lower pairs           ", np.max(-E0(d, w, h)))
    t = rng.normal(size=N) * np.exp(rng.uniform(-5, 5, N))
    hh = rng.uniform(0, 0.3, N)
    print(" F  sigma-2/((1-h^2)(1+t^2))   ", np.max(sigma(t, hh) - 2 / ((1 - hh ** 2) * (1 + t ** 2))))
    # far pairs: d>=4 max(1,h)
    hf = np.exp(rng.uniform(np.log(0.01), np.log(100), N))
    M = np.maximum(1, hf)
    tj = 5 * M * np.exp(rng.exponential(0.7, N))
    u = rng.uniform(-1, 1, N)
    dL = tj + u
    print(" far: 2w/(3d^3) bound on E0(0) ", np.max(E0(dL, w, 0 * w) - 2 * w / (3 * dL ** 3)))
    print(" far: E0(h)-E0(0)              ", np.max(E0(dL, w, hf) - E0(dL, w, 0 * w)))
    print(" far: -E0                      ", np.max(-E0(dL, w, hf)))
    print(" far: (24/17)w/d^2-s           ", np.max(24 / 17 * w / dL ** 2 - slope(dL, w, hf)))
    print(" far: sigma-(17/8)/t^2 (t>=4M) ", np.max(sigma(dL, hf) - (17 / 8) / dL ** 2))
    # the ratio bounds used in the two cases, lower pairs at the allowed distance
    a2 = 1 - w ** 2
    u = np.sqrt(a2) * rng.uniform(-1, 1, N)
    d = np.sqrt((1 - 3 * w ** 2) / 2) * np.exp(rng.exponential(1.0, N))
    v = np.maximum(-u, 0)
    s = slope(d, w, h)
    EL = kap(d, w, h) - w * kc(d - u, 1.0, h)
    ER = -kap(d, w, h) + w * kc(d + u, 1.0, h)
    CL = np.where(u >= 0, Hb(w), Hb(w) + v / (1 - w ** 2))
    CR = np.where(u >= 0, 0.0, v * (1 + v + 3 * w ** 2) / (1 - w ** 2))
    print(" ratio left  E_L-C s           ", np.max(EL - CL * s))
    print(" ratio right E_R-C s           ", np.max(ER - CR * s))
    print(" largest true E_L/s at u>=0    ", np.max((EL / s)[u >= 0]), " bound Hb at w=.3:", Hb(0.3))
    print(" v A_x - w/(2a)                ", np.max(v * Ax(u, w) - w / (2 * np.sqrt(a2))))


def margins():
    print("margins of the proof")
    print("   w    arccos(w)-pi w  B1(w)   case1      (pi/2)(1-2w)  case2 bound (left,far)  case2 bound (right)")
    for w in (0.05, 0.1, 0.15, 0.2, 0.25, 0.28, 0.29, 0.30, 0.302, 0.305, 0.31):
        a = np.sqrt(1 - w * w)
        vA = w / (2 * a)
        L = B1(w) + max(1 / (1 - w * w), 1.51) * vA
        R = max((2 + 3 * w * w) / (1 - w * w), 2.17) * vA
        print(f"  {w:.3f}  {np.arccos(w)-pi*w:.5f}  {B1(w):.5f}  {np.arccos(w)-pi*w-B1(w):+.5f}  {pi/2*(1-2*w):.5f}  {L:.5f}  {R:.5f}")


def nonreal(z, phi, w, tol=1e-7):
    co = (np.poly(z - 1j * w) * np.exp(1j * phi * w)).real
    co = co / np.max(np.abs(co))
    k = 0
    while k < len(co) - 1 and abs(co[k]) < 1e-13:
        k += 1
    co = co[k:]
    if len(co) <= 1:
        return False
    r = np.roots(co)
    return np.any(np.abs(r.imag) > tol * np.maximum(1.0, np.abs(r)))


def wstar(z, phi):
    last = 1.0
    for w in np.linspace(0.998, 0.002, 500):
        if nonreal(z, phi, w):
            lo, hi = w, last
            for _ in range(30):
                m = 0.5 * (lo + hi)
                if nonreal(z, phi, m):
                    lo = m
                else:
                    hi = m
            return hi
        last = w
    return 0.0


def chain():
    print("the chain of the proof on random F = (z^2+1) prod lower pairs prod reals exp(phi z) just above w*")
    worst_ratio, viol_level, viol_top, viol_ineq, n = 9.0, 0, 0, 0, 0
    worstD = np.inf
    for trial in range(1200):
        nr = rng.integers(0, 6)
        reals = rng.uniform(-5, 5, nr)
        npair = rng.integers(0, 2)
        pm = rng.uniform(-4, 4, npair)
        ph = rng.uniform(0.02, 0.3, npair)
        phi = rng.uniform(-3.2, 3.2)
        z = np.concatenate(([1j, -1j], reals.astype(complex), pm + 1j * ph, pm - 1j * ph))
        ws = wstar(z, phi)
        if ws <= 0 or (npair and ws <= ph.max()):
            continue
        w = ws + 0.01 * (1 - ws)
        worst_ratio = min(worst_ratio, ws)
        a = np.sqrt(1 - w * w)
        us = np.linspace(-a, a, 4001)[1:-1]
        Tall = phi * w + np.sum([np.angle(us - r + 1j * w) for r in reals], axis=0) if nr else phi * w * np.ones_like(us)
        if npair:  # continuous in x because y > h (the principal Arg of the product would jump on the ray above c_j)
            Tall = Tall + np.arctan2(w - ph[0], us - pm[0]) + np.arctan2(w + ph[0], us - pm[0])
        Phi = Tall + Afun(us, w)
        lev = np.floor((Phi - pi / 2) / pi)
        idx = np.where(lev[1:] > lev[:-1])[0]
        if len(idx) == 0:
            viol_level += 1
            continue
        u = us[idx[0]]
        sgn = 1.0
        rs, pms, ph_ = reals.copy(), pm.copy(), ph.copy()
        ph0 = phi

        def build(u, rs, pms, ph_, ph0):
            def Tn(x, y):
                val = ph0 * y
                for r in rs:
                    val = val + np.arctan2(y, x - r) - (pi if r > u else 0)
                for mj, hj in zip(pms, ph_):
                    val = val + np.arctan2(y - hj, x - mj) + np.arctan2(y + hj, x - mj) - (2 * pi if mj > u else 0)
                return val
            return Tn

        Tn = build(u, rs, pms, ph_, ph0)
        L = Tn(u, w) + Afun(u, w)
        if L < 0:
            u, rs, pms, ph0 = -u, -rs, -pms, -ph0
            Tn = build(u, rs, pms, ph_, ph0)
            L = Tn(u, w) + Afun(u, w)
        n += 1
        if abs(((L - pi / 2) / pi) - round((L - pi / 2) / pi)) > 2e-2:
            viol_level += 1
        Ttop = Tn(0.0, 1.0)
        if Ttop > L + pi / 2 + 1e-2 or Ttop < L - pi / 2 - 1e-2:
            viol_top += 1
        D0 = pi / 2 * (1 - 2 * w) - Afun(u, w)
        Diff = Tn(u, w) - w * Ttop
        S = 0.0
        for r in rs:
            S += w / ((u - r) ** 2 + w ** 2)
        for mj, hj in zip(pms, ph_):
            S += slope(abs(u - mj), w, hj)
        if Diff < D0 - 1e-2 or S > Ax(u, w) + 1e-2:
            viol_ineq += 1
        worstD = min(worstD, Diff - D0)
    print(f"  functions {n}; least w* {worst_ratio:.4f}; level off pi/2+pi Z {viol_level}; T(i) outside [L-pi/2,L+pi/2] {viol_top}; Diff<D0 or slope bound broken {viol_ineq}; least Diff-D0 {worstD:.4f}")


def grid(nw, nxp, nxq, hq, xq_top=False, wmax=K1 * 0.9999):
    worst = (np.inf, None)
    dist = np.geomspace(1e-3, 1e4, 240)
    hs = np.linspace(0, 1, 9)
    per_w = []
    for w in np.linspace(0.02, wmax, nw):
        lam = w / hq
        ap = np.sqrt(1 - w * w)
        aq = np.sqrt(max(1 - hq * hq, 0.0))
        best_w = np.inf
        H = (w * hs)[:, None]
        for xp in np.linspace(-ap, ap, nxp + 2)[1:-1]:
            Au = Ax(xp, w)
            xqs = [0.0] if xq_top else np.linspace(-aq, aq, nxq)
            for xq in xqs:
                Aq = -pi / 2 if xq_top else Afun(xq, hq)
                D = (pi / 2 - Afun(xp, w)) - lam * (pi / 2 - Aq)
                d = dist[None, :]
                s = slope(d, w, H)
                ok = s <= Au
                ok_real = Kf(w, d) <= Au
                EL = kap(d, w, H) - lam * kc(xq - xp + d, hq, H)
                ER = -kap(d, w, H) - lam * (kc(xq - xp - d, hq, H) - 2 * pi)
                muP = max(np.max((EL / s)[ok], initial=-9), np.max((ER / s)[ok], initial=-9))
                muR = max(np.max((EL[0] / s[0])[ok_real[0]], initial=-9), np.max((ER[0] / s[0])[ok_real[0]], initial=-9))
                mu = max(muP, muR, xq - xp, 0.0)
                margin = D - mu * Au
                best_w = min(best_w, margin)
                if margin < worst[0]:
                    worst = (margin, w, xp, xq, D, mu, Au)
        per_w.append((w, best_w))
    return worst, per_w


def n14():
    """High precision roots for the law n=14, beta=2.6 (double precision chebroots gave a false 0.7606)."""
    import mpmath as mp
    mp.mp.dps = 90
    n, A, beta = 14, mp.mpf(1) / 2, mp.mpf(26) / 10
    delta = beta / (n * mp.sqrt(3))
    eps = mp.sinh(delta) * mp.sqrt(1 - A ** 2 / mp.cosh(delta) ** 2)
    # power coefficients of q(y) = y^n ((y+A)^2+eps^2), then to Chebyshev by T_k recurrence
    N = n + 2
    q = [mp.mpf(0)] * (N + 1)
    q[n] = A ** 2 + eps ** 2
    q[n + 1] = 2 * A
    q[n + 2] = mp.mpf(1)
    T = [[mp.mpf(1)] + [mp.mpf(0)] * N, [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (N - 1)]
    for k in range(2, N + 1):
        new = [mp.mpf(0)] * (N + 1)
        for j in range(N + 1):
            if j > 0:
                new[j] += 2 * T[k - 1][j - 1]
            new[j] -= T[k - 2][j]
        T.append(new)
    c = [mp.mpf(0)] * (N + 1)
    rem = q[:]
    for k in range(N, -1, -1):  # triangular solve, top down
        c[k] = rem[k] / T[k][k]
        for j in range(N + 1):
            rem[j] -= c[k] * T[k][j]
    print("n=14, beta=2.6, delta =", mp.nstr(delta, 8), ": mpmath roots at 90 digits, largest |Im root| and largest |Re root| for T_wF in y = cos z")
    for r in ("0.30", "0.45", "0.48", "0.49", "0.50", "0.60", "0.70", "0.76", "0.80", "0.95"):
        w = mp.mpf(r) * delta
        co = [mp.mpf(0)] * (N + 1)
        for k in range(N + 1):
            ck = c[k] * mp.cosh(k * w)
            for j in range(N + 1):
                co[j] += ck * T[k][j]
        roots = mp.polyroots(co[::-1], maxsteps=2000, extraprec=800)
        print(f"  w/delta={r}: max|Im|={mp.nstr(max(abs(mp.im(x)) for x in roots), 4)}  max|Re|={mp.nstr(max(abs(mp.re(x)) for x in roots), 8)}", flush=True)


if __name__ == "__main__":
    scale = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    if scale == 0:
        chain()
        n14()
        sys.exit()
    lemmas()
    margins()
    chain()
    print("grid: q at the height r_c, real zeros and lower pairs (h = w*j/8), m = 0 worst case")
    worst, per_w = grid(int(25 * scale), int(41 * scale), int(9 * scale), RC)
    print(" worst margin", np.round(worst, 5))
    print(" per w:", " ".join(f"{w:.3f}:{m:.3f}" for w, m in per_w))
    print("grid: q at the top (x=0,y=1), A(q) replaced by -pi/2: margin per w, to w=0.38")
    worst, per_w = grid(int(19 * scale), int(41 * scale), 1, 1.0 - 1e-12, xq_top=True, wmax=0.38)
    print(" per w:", " ".join(f"{w:.3f}:{m:.3f}" for w, m in per_w))
