#!/usr/bin/env python3
"""Track C, script 3. Sums of scaled uniforms against the Polya law.

S_b = sum_{k>=1} U_k / (2(k+b)), the harmonic uniform law, and fits with the first J
coefficients free. Every law here is ferromagnetic with no computation of zeros and no
positivity proof. Controls: the Gaussian law and the second row of the Bessel table of the
book, whose total variation distances are stated in the book.

Run: timeout 550 python3 -I C_uniform_fit.py
Status: computed in double precision on a grid, not certified.
"""
import numpy as np
from scipy import optimize, special

VARX = 0.0462099862308379
KAPPA4 = -4.46071191423236e-4
KAPPA6 = 3.46017435362336e-5
XI0 = 0.4971207781883141
EG = 0.5772156649015329

TH = np.arange(0, 45.0001, 0.02)
XS = np.arange(0, 1.8001, 0.001)
COS = np.cos(np.outer(XS, TH))
W = np.ones_like(TH) * 0.02
W[0] = W[-1] = 0.01


def p_density(x):
    x = np.abs(x)
    s = np.zeros_like(x)
    for n in range(1, 8):
        y = np.pi * n**2 * np.exp(2 * x)
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * x) - 3 * np.pi * n**2 * np.exp(2.5 * x)) * np.exp(-y)
    return 2 * s / XI0


P = p_density(XS)


def tv_of_phi(phi):
    """TV between the density with even real transform phi and p, on the half line times two / 2."""
    f = COS @ (W * phi) / np.pi
    d = np.abs(f - P)
    return float(np.sum(d) - 0.5 * (d[0] + d[-1])) * 0.001


def sinc_prod(cs):
    """Product of sin(c th)/(c th) over the array cs, on TH."""
    out = np.ones_like(TH)
    for c in cs:
        out = out * np.sinc(c * TH / np.pi)
    return out


def tail_phi(a, b, Kstart, K=500):
    """Product over k = Kstart..K of sinc(th a/(k+b)), then the cumulant series for k > K."""
    ks = np.arange(Kstart, K + 1)
    out = sinc_prod(a / (ks + b)) if len(ks) else np.ones_like(TH)
    # log sinc x = -x^2/6 - x^4/180 - x^6/2835 - ..., with sums over k > K from polygamma
    s2 = a**2 * special.polygamma(1, K + 1 + b)
    s4 = a**4 * special.polygamma(3, K + 1 + b) / 6
    s6 = a**6 * special.polygamma(5, K + 1 + b) / 120
    return out * np.exp(-(TH**2) * s2 / 6 - (TH**4) * s4 / 180 - (TH**6) * s6 / 2835)


def cum_of(cs, a, b, Kstart):
    """kappa_2, kappa_4, kappa_6 for free coefficients cs plus the harmonic tail."""
    def pw(m):
        # sum over the tail of (a/(k+b))^(2m), k >= Kstart
        return a ** (2 * m) * float(special.zeta(2 * m, Kstart + b))
    k2 = (np.sum(cs**2) + pw(1)) / 3
    k4 = -2 / 15 * (np.sum(cs**4) + pw(2))
    k6 = 16 / 63 * (np.sum(cs**6) + pw(3))  # kappa_6(U) = 2^6 B_6/6 = 16/63
    return k2, k4, k6


def main():
    print("controls")
    gauss = np.exp(-VARX * TH**2 / 2)
    print(f"  TV(Gaussian, p) = {tv_of_phi(gauss):.4e}  (book 1.48e-2)")
    # Bessel row: cosh(lam t) exp(-a cosh(l t)), lam=2.146184, a=5.483952, l=2.062866, on the x grid directly
    lam, a_, l_ = 2.146184, 5.483952, 2.062866
    xx = np.arange(0, 6, 0.0005)
    q = np.cosh(lam * xx) * np.exp(-a_ * np.cosh(l_ * xx))
    norm = 2 * (np.sum(q) - 0.5 * q[0]) * 0.0005
    d = np.abs(q / norm - p_density(xx))
    print(f"  TV(Bessel row 2, p) = {(np.sum(d)-0.5*d[0])*0.0005:.4e}  (book < 6.1e-6)")

    bstar = optimize.brentq(lambda b: special.digamma(1 + b) - EG, 0.5, 3)
    print(f"harmonic law S_b with psi(1+b)=Euler gamma: b*={bstar:.8f}, first zero 2 pi (1+b*)={2*np.pi*(1+bstar):.5f}")
    print(f"  tail constant pi e^(psi(1+b)-gamma) = pi; zero count (T/2pi)log(T/2pi e) + o(T)")
    k2, k4, k6 = cum_of(np.array([]), 0.5, bstar, 1)
    print(f"  Var={k2:.8f} (X {VARX:.8f}, rel {k2/VARX-1:+.2e}); k4={k4:.6e} (X {KAPPA4:.6e}, rel {k4/KAPPA4-1:+.2e}); k6={k6:.6e} (X {KAPPA6:.6e}, rel {k6/KAPPA6-1:+.2e})")
    print(f"  TV(S_b*, p) = {tv_of_phi(tail_phi(0.5, bstar, 1)):.4e}")

    bvar = optimize.brentq(lambda b: special.polygamma(1, 1 + b) / 12 - VARX, 0.5, 3)
    print(f"b matched to the variance: b={bvar:.8f}, TV = {tv_of_phi(tail_phi(0.5, bvar, 1)):.4e}")
    r = optimize.minimize_scalar(lambda b: tv_of_phi(tail_phi(0.5, b, 1)), bounds=(1.0, 1.5), method="bounded", options={"xatol": 1e-7})
    print(f"b minimising TV: b={r.x:.8f}, TV={r.fun:.4e}")

    best = {}
    x0 = np.array([])
    for J in (1, 2, 3, 4, 6, 8, 12):
        def obj(v):
            cs = np.abs(v[:J])
            b = v[J]
            a = 0.5
            if b <= -J + 0.05:
                return 1.0
            return tv_of_phi(sinc_prod(cs) * tail_phi(a, b, J + 1))
        start = np.concatenate([0.5 / (np.arange(1, J + 1) + bstar), [bstar]])
        if J > 1 and "v" in best:
            prev = best["v"]
            start[: len(prev) - 1] = prev[:-1]
            start[-1] = prev[-1]
        res = None
        v = start
        for it in range(4):
            res = optimize.minimize(obj, v, method="Powell", options={"xtol": 1e-6, "ftol": 1e-9, "maxiter": 40})
            v = res.x
        cs = np.abs(v[:J])
        k2, k4, k6 = cum_of(cs, 0.5, v[J], J + 1)
        print(f"J={J:2d}: TV={res.fun:.4e} b={v[J]:.5f} Var rel {k2/VARX-1:+.1e}  pi/c_k={np.array2string(np.pi/cs, precision=3, max_line_width=200)}")
        best["v"] = v
    # the last fit, with its coefficients in full
    v = best["v"]
    print("last fit coefficients c_k:", np.array2string(np.abs(v[:-1]), precision=8, max_line_width=200), "b=", v[-1])


if __name__ == "__main__":
    main()
