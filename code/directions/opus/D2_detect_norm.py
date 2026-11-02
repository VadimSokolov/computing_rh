#!/usr/bin/env python3
"""D2: norm level detection of a planted quadruple 1/2 +- delta +- i gamma, added to the first 8000 ordinates.

Two windows thresholds in a, on a grid of step 0.05:
 a_band: least eigenvalue of Weil's form on a band of even modes, frequencies within 45 of gamma, is below -1e-9 (float64).
 a_f: the explicit function of Theorem 6, f(u) = sinh(delta |u|) sin(gamma |u|) on [-a,a], has W(f) < 0.
Compared with a_s, the root of (1/pi)(sinh(2 a delta)/(2 delta) - a) = rho_0(gamma), with the horizon (1/2)log(gamma/2pi),
and with (1/2)log(gamma/4pi).
Run: timeout 600 python3 -I D2_detect_norm.py
Status: float64, tail above height 8148 from the density (a model for the band; an upper bound for the explicit function).
"""
import math
import numpy as np
ROOT = '.'
gam = np.array([float(x) for x in open(ROOT + '/data/zeros_8000.txt') if x.strip()])
TOP = gam[-1]


def rho0(t):
    return np.log(t / (2 * math.pi)) / (2 * math.pi)


def band_gram(g, a, margin=45.0):
    L = a
    nlo = max(1, int(math.floor((g - margin) * L / math.pi)))
    nhi = int(math.ceil((g + margin) * L / math.pi))
    n = np.arange(nlo, nhi + 1)
    y = (math.pi * n / L) ** 2
    sgn = (-1.0) ** n
    M = (np.sin(L * gam) * 2 * gam)[:, None] / (gam[:, None] ** 2 - y[None, :]) * sgn
    A = 2 * M.T @ M
    u = np.linspace(0, 12, 2401)
    t = TOP * np.exp(u)
    w = rho0(t) * t * (u[1] - u[0])
    w[0] *= 0.5
    w[-1] *= 0.5
    Fm = (2 * t)[:, None] / (t[:, None] ** 2 - y[None, :]) * sgn
    A += (Fm.T * w) @ Fm      # 2 * mean sin^2 (1/2) * density
    return A, n, y, sgn


def band_lmin(A, n, y, sgn, g, d, L):
    z = complex(g, -d)
    v = np.sin(L * z) * 2 * z / (z * z - y) * sgn
    Q = 4 * (np.outer(v.real, v.real) - np.outer(v.imag, v.imag))
    return np.linalg.eigvalsh((A + Q) / L)[0]


def H(x, d, a):
    """H(x) = -(1/pi) int_0^a sinh(d u) sin(x u) du, x complex, x^2 + d^2 != 0."""
    return -(d * math.cosh(d * a) * np.sin(x * a) - x * math.sinh(d * a) * np.cos(x * a)) / (d * d + x * x) / math.pi


def explicit_ratio(g, d, a):
    A0 = (math.sinh(2 * a * d) / (4 * d) - a / 2) / math.pi
    Fz = 1j * A0 - H(complex(2 * g, -d), d, a)
    F = H(gam - g, d, a) - H(gam + g, d, a)
    # tail above TOP: upper bound with |H(x)| <= 2 sinh(a d)/(pi |x|)
    u = np.linspace(0, 12, 2401)
    t = TOP * np.exp(u)
    tail = np.sum(2 * rho0(t) * (2 * 2 * math.sinh(a * d) / (math.pi * (t - g))) ** 2 * t) * (u[1] - u[0])
    W = 2 * np.sum(F * F) + tail + 4 * (Fz * Fz).real
    uu = np.linspace(0, a, 400001)
    val = (np.sinh(d * uu) * np.sin(g * uu) / math.pi) ** 2
    nrm = 2 * (np.sum(val) - 0.5 * (val[0] + val[-1])) * (uu[1] - uu[0])
    return W / nrm


def a_sinh(g, d):
    r = float(rho0(g))
    lo, hi = 0.01, 40.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if (math.sinh(2 * m * d) / (2 * d) - m) / math.pi > r:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def main():
    book = {10: (1.20, 1.15, 1.15, 1.10), 79: (1.85, 1.80, 1.75, 1.75), 649: (2.65, 2.65, 2.60, 2.50), 1517: (3.05, 3.00, 2.95, 2.45)}
    deltas = (0.01, 0.05, 0.2, 0.3, 0.45)
    agrid = np.arange(0.60, 4.60, 0.05)
    print(' k gamma delta a_band a_f a_s horizon (1/2)log(gamma/4pi) book(i) lmin_at_a_band')
    for k in (10, 79, 649, 1517, 3000, 5000, 7000):
        g = 0.5 * (gam[k - 1] + gam[k])
        grams = {}
        res = {d: None for d in deltas}
        lm = {}
        for a in agrid:
            a = round(float(a), 2)
            if all(res[d] is not None for d in deltas):
                break
            A, n, y, sgn = band_gram(g, a)
            for d in deltas:
                if res[d] is None:
                    l = band_lmin(A, n, y, sgn, g, d, a)
                    if l < -1e-9:
                        res[d] = a
                        lm[d] = l
        for i, d in enumerate(deltas):
            af = None
            for a in np.arange(0.60, 12.0, 0.05):
                if explicit_ratio(g, d, float(a)) < 0:
                    af = round(float(a), 2)
                    break
            bk = ''
            if k in book and d != 0.3:
                bk = book[k][(0, 1, 2, 3)[[0.01, 0.05, 0.2, 0.45].index(d)]]
            print(f"{k:5d} {g:9.2f} {d:5.2f} {str(res[d]):>6} {str(af):>6} {a_sinh(g, d):6.2f} {0.5 * math.log(g / (2 * math.pi)):6.3f} {0.5 * math.log(g / (4 * math.pi)):6.3f} {str(bk):>5} {lm.get(d, float('nan')):.2e}")


if __name__ == '__main__':
    main()
