#!/usr/bin/env python3
"""Track C, script 5. The Davenport and Heilbronn law as a control: a genuine positive even density,
known not to be ferromagnetic, and still close in total variation to explicit ferromagnetic laws.

Run: timeout 500 python3 -I C_dh.py
Status: computed in mpmath and double precision, not certified.
"""
import numpy as np
import mpmath as mp
from scipy import optimize, special

mp.mp.dps = 25
KAPPA = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [mp.mpf(1), KAPPA, -KAPPA, mp.mpf(-1)]
Af = [float(a) for a in A] + [0.0]


def f_dh(s):
    return mp.mpf(5) ** (-s) * sum(a * mp.zeta(s, mp.mpf(r + 1) / 5) for r, a in enumerate(A))


def Lam(s):
    return (5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * f_dh(s)


LAM0 = Lam(mp.mpf(1) / 2)


def LDH(z):
    return Lam(mp.mpf(1) / 2 + z) / LAM0


def theta_dh(y):
    y = np.asarray(y, dtype=float)
    s = np.zeros_like(y)
    for n in range(1, 40):
        a = Af[(n - 1) % 5]
        if a != 0:
            s = s + a * n * np.exp(-np.pi * n**2 * y / 5)
    return s


def kernel(t):
    """e^{3t/2} theta(e^{2t}) for t >= 0, extended evenly; the series is used at e^{2|t|} >= 1."""
    t = np.abs(np.asarray(t, dtype=float))
    return np.exp(1.5 * t) * theta_dh(np.exp(2 * t))


def main():
    print("kappa =", mp.nstr(KAPPA, 15), " Lam(1/2) =", mp.nstr(LAM0, 15))
    # evenness of the raw kernel, positivity
    for t in (0.2, 0.5, 0.9):
        a = np.exp(1.5 * t) * theta_dh(np.exp(2 * t))
        b = np.exp(-1.5 * t) * theta_dh(np.exp(-2 * t))
        print(f"  evenness at t={t}: {a:.12e} vs {b:.12e}")
    xs = np.arange(0, 4.2001, 0.0005)
    k = kernel(xs)
    print(f"  min kernel on [0,3.0]: {k[xs<=3.0].min():.3e} (positive means a density), at 0: {k[0]:.6f}")
    norm = 2 * (np.sum(k) - 0.5 * k[0] - 0.5 * k[-1]) * 0.0005
    p = k / norm
    var = 2 * (np.sum(xs**2 * p)) * 0.0005
    k4 = 2 * (np.sum(xs**4 * p)) * 0.0005 - 3 * var**2
    print(f"  normaliser {norm:.10f}, Var={var:.8f}, kappa_4={k4:.5e}")
    # the transform of the density is Lam(1/2+z)/Lam(1/2)
    for z in (0.3 + 0j, 0.2 + 3j, 0.0 + 5.0j):
        direct = (2 * np.sum(p * np.cosh(z * xs)) - p[0]) * 0.0005
        print(f"  z={z}: quadrature {direct:.10f}, Lam ratio {complex(LDH(mp.mpc(z))):.10f}")
    # ferromagnetic approximants
    xg = np.arange(0, 3.6, 0.001)
    pg = kernel(xg) / norm

    def tv_q(q):
        d = np.abs(q - pg)
        return float(np.sum(d) - 0.5 * d[0] - 0.5 * d[-1]) * 0.001

    def bessel(par):
        lam, a, l = par
        xx = np.arange(0, 6, 0.0005)
        q = np.cosh(lam * xx) * np.exp(-abs(a) * np.cosh(abs(l) * xx))
        nn = 2 * (np.sum(q) - 0.5 * q[0]) * 0.0005
        return np.cosh(lam * xg) * np.exp(-abs(a) * np.cosh(abs(l) * xg)) / nn

    obj = lambda par: tv_q(bessel(par))
    r = optimize.minimize(obj, [1.8, 1.0, 2.0], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-11, "maxiter": 4000})
    r = optimize.minimize(obj, r.x, method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-12, "maxiter": 4000})
    print(f"  Bessel law cosh(lam t) exp(-a cosh(l t)), ferromagnetic by the book's theorem: TV={r.fun:.4e}, lam,a,l={np.round(r.x,6)}")

    # harmonic uniform sums on the DH law
    TH = np.arange(0, 45.0001, 0.02)
    Wt = np.ones_like(TH) * 0.02
    Wt[0] = Wt[-1] = 0.01
    XS = np.arange(0, 2.6001, 0.001)
    COS = np.cos(np.outer(XS, TH))
    P = kernel(XS) / norm

    def tv_phi(phi):
        fq = COS @ (Wt * phi) / np.pi
        d = np.abs(fq - P)
        return float(np.sum(d) - 0.5 * d[0] - 0.5 * d[-1]) * 0.001

    def tailphi(b, Kstart, K=500):
        ks = np.arange(Kstart, K + 1)
        out = np.ones_like(TH)
        for c in 0.5 / (ks + b):
            out = out * np.sinc(c * TH / np.pi)
        s2 = 0.25 * special.polygamma(1, K + 1 + b)
        s4 = 0.0625 * special.polygamma(3, K + 1 + b) / 6
        return out * np.exp(-(TH**2) * s2 / 6 - (TH**4) * s4 / 180)

    bstar = optimize.brentq(lambda b: special.digamma(1 + b) - (0.5772156649015329 - np.log(5)), -0.95, 3)
    print(f"  tail constant pi/5 needs psi(1+b)=gamma-log 5: b={bstar:.6f}; first zero 2 pi (1+b)={2*np.pi*(1+bstar):.4f}")
    print(f"  TV(harmonic uniform law with that b, DH law) = {tv_phi(tailphi(bstar,1)):.4e}")
    for J in (1, 2, 3):
        def o2(v):
            out = tailphi(v[J], J + 1)
            for c in np.abs(v[:J]):
                out = out * np.sinc(c * TH / np.pi)
            return tv_phi(out)
        v0 = np.concatenate([0.5 / (np.arange(1, J + 1) + bstar), [bstar]])
        res = optimize.minimize(o2, v0, method="Powell", options={"xtol": 1e-6, "ftol": 1e-9, "maxiter": 25})
        print(f"  J={J}: TV={res.fun:.4e}, pi/c_k={np.round(np.pi/np.abs(res.x[:J]),4)}, b={res.x[J]:.5f}")

    # the known zero off the line
    z0 = mp.findroot(f_dh, mp.mpc("0.8085", "85.6993"))
    print("  zero off the line: s =", mp.nstr(z0, 12), " |f(s)| =", mp.nstr(abs(f_dh(z0)), 3))
    zc = z0 - mp.mpf(1) / 2
    print("  size of the law at that height on the critical line: |L_DH(i 85.7)| =", mp.nstr(abs(LDH(mp.mpc(0, mp.im(z0)))), 4))
    rmin = {}
    for r in (0.05, 0.1, 0.2, 0.3):
        vals = [abs(LDH(zc + r * mp.expjpi(2 * mp.mpf(j) / 72))) for j in range(72)]
        rmin[r] = min(vals)
        print(f"  circle radius {r}: min |L_DH| = {mp.nstr(rmin[r],4)}")
    best_r = max(rmin, key=lambda u: rmin[u])
    R = float(mp.re(zc)) + best_r
    sig = np.sqrt(var) * 1.001
    # solve 5 eps exp(R sig sqrt(2 log 1/eps)) = min |L_DH|
    tgt = float(rmin[best_r])
    if tgt > 0:
        eps = optimize.brentq(lambda e: np.log(5) + e + R * sig * np.sqrt(2 * (-e)) - np.log(tgt), np.log(1e-120), np.log(1e-5))
        print(f"  lower bound for delta(DH law) from the zero: eps = {np.exp(eps):.3e} (R={R:.3f}, sigma={sig:.4f})")

    # lobes of the DH transform on the critical line: which zeros are real below 85, amplitude
    zs = []
    ts = np.arange(0.5, 90, 0.1)
    vals = [float(mp.re(LDH(mp.mpc(0, t)))) for t in ts]
    for i in range(len(ts) - 1):
        if vals[i] * vals[i + 1] < 0:
            zs.append(ts[i] - vals[i] * (ts[i + 1] - ts[i]) / (vals[i + 1] - vals[i]))
    zs = np.array(zs)
    print(f"  sign changes on the line below 90: {len(zs)}; first: {np.round(zs[:8],3)}")


if __name__ == "__main__":
    import builtins
    _p = builtins.print
    builtins.print = lambda *a, **k: _p(*a, flush=True, **k)
    main()
