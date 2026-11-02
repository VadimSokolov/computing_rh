#!/usr/bin/env python3
"""D1 exploration in double precision (NOT a certificate).

Weil form on L^2[-L,L] in the periodic cosine and sine modes, from the closed
forms derived in D1.md, Section 2. Prints least eigenvalues of
    W_kappa = P + A - kappa T
for the even and the odd sector, with T the prime 2 term.

Usage: timeout 300 python3 -I D1_explore.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sici
import scipy.linalg as sla
from functools import lru_cache
from scipy.special import roots_legendre

GAMMA = 0.57721566490153286061
K2 = np.sqrt(2.0) * np.log(2.0)   # kappa_2 = 2 Lambda(2) / sqrt(2)
LOG2 = np.log(2.0)


@lru_cache(maxsize=8)
def _gl(nq):
    return roots_legendre(nq)


def h_kernel(x):
    """h(x) = exp(x/2)/sinh(x) - 1/x, by series near 0."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x < 1e-3
    xs = x[small]
    # exp(x/2)/sinh(x) = 1/x + 1/2 + (1/8 - 1/6) x + ... ; three terms suffice below 1e-3
    out[small] = 0.5 + (1.0 / 8.0 - 1.0 / 6.0) * xs + (1.0 / 48.0 - 1.0 / 12.0) * xs ** 2
    xb = x[~small]
    out[~small] = np.exp(xb / 2) / np.sinh(xb) - 1.0 / xb
    return out


def integrals(L, N, nq=None):
    """I_n, and C_n = int_0^{2L} h(x)(1-x/2L) cos(w_n x) dx, for n=0..N."""
    if nq is None:
        nq = max(4000, 8 * N)
    t, w = _gl(nq)
    x = L * (t + 1.0)
    w = L * w
    hx = h_kernel(x)
    om = np.pi * np.arange(N + 1) / L
    In = np.zeros(N + 1)
    Cn = np.zeros(N + 1)
    tri = 1.0 - x / (2 * L)
    # blocks to save memory
    for start in range(0, N + 1, 256):
        idx = np.arange(start, min(start + 256, N + 1))
        ph = np.outer(om[idx], x)
        Cn[idx] = (np.cos(ph) * (hx * tri * w)).sum(axis=1)
        In[idx] = (np.sin(ph) * (hx * w)).sum(axis=1)
    si, _ = sici(2 * np.pi * np.arange(1, N + 1))
    In[1:] += si
    return In, Cn


def matrices(L, N, c_list=((LOG2, K2),)):
    """Return dict with even and odd arch (A), pole (P), prime (T) matrices.

    Even basis: e_0=(2L)^-1/2, e_n=L^-1/2 cos(pi n x/L), n=1..N, size N+1.
    Odd basis: s_n=L^-1/2 sin(pi n x/L), n=1..N, size N.
    c_list is a list of (log n, 2 Lambda(n)/sqrt(n)) with log n < 2L.
    """
    In, Cn = integrals(L, N)
    n = np.arange(N + 1, dtype=float)
    om = np.pi * n / L
    _, ci = sici(2 * np.pi * n[1:])
    # even sector
    Ae = np.zeros((N + 1, N + 1))
    Ao = np.zeros((N, N))
    Te = np.zeros((N + 1, N + 1))
    To = np.zeros((N, N))
    base = np.log(n[1:] / (2 * L)) - ci - Cn[1:]
    d_e = base + In[1:] / (2 * np.pi * n[1:])
    d_o = base - In[1:] / (2 * np.pi * n[1:])
    Ae[0, 0] = -GAMMA - np.log(np.pi) - np.log(4 * L) + 1.0 - Cn[0]
    Ae[np.arange(1, N + 1), np.arange(1, N + 1)] = d_e
    Ao[np.arange(N), np.arange(N)] = d_o
    m_, n_ = np.meshgrid(n[1:], n[1:], indexing="ij")
    sgn = (-1.0) ** (m_ + n_)
    denom = np.pi * (m_ ** 2 - n_ ** 2)
    np.fill_diagonal(denom, 1.0)
    I_m = np.outer(In[1:], np.ones(N))
    I_n = I_m.T
    off_e = -sgn * (n_ * I_n - m_ * I_m) / denom
    off_o = -sgn * (m_ * I_n - n_ * I_m) / denom
    np.fill_diagonal(off_e, 0.0)
    np.fill_diagonal(off_o, 0.0)
    Ae[1:, 1:] += off_e
    Ao += off_o
    # e_0 with e_n: -int k g_{0n} = +(-1)^n I_n/(sqrt2 pi n)
    v = ((-1.0) ** n[1:]) * In[1:] / (np.sqrt(2.0) * np.pi * n[1:])
    Ae[0, 1:] = v
    Ae[1:, 0] = v
    # pole
    sh = np.sinh(L / 2)
    ap = np.empty(N + 1)
    ap[0] = (2 * L) ** -0.5 * sh / 0.25
    ap[1:] = L ** -0.5 * ((-1.0) ** n[1:]) * sh / (om[1:] ** 2 + 0.25)
    am = -L ** -0.5 * ((-1.0) ** n[1:]) * 2 * om[1:] * sh / (om[1:] ** 2 + 0.25)
    Pe = 2 * np.outer(ap, ap)
    Po = -2 * np.outer(am, am)
    # prime power terms
    for c, wgt in c_list:
        if c >= 2 * L:
            continue
        S = np.sin(om * c)
        Cc = np.cos(om * c)
        ge = np.zeros((N + 1, N + 1))
        go = np.zeros((N, N))
        s_m = np.outer(S[1:], np.ones(N))
        s_n = s_m.T
        oe = sgn * (n_ * s_n - m_ * s_m) / denom
        oo = sgn * (m_ * s_n - n_ * s_m) / denom
        np.fill_diagonal(oe, 0.0)
        np.fill_diagonal(oo, 0.0)
        ge[1:, 1:] = oe
        go = oo
        tri = 1.0 - c / (2 * L)
        de = tri * Cc[1:] - S[1:] / (2 * np.pi * n[1:])
        do = tri * Cc[1:] + S[1:] / (2 * np.pi * n[1:])
        ge[np.arange(1, N + 1), np.arange(1, N + 1)] = de
        go[np.arange(N), np.arange(N)] = do
        ge[0, 0] = tri
        v0 = -((-1.0) ** n[1:]) * S[1:] / (np.sqrt(2.0) * np.pi * n[1:])
        ge[0, 1:] = v0
        ge[1:, 0] = v0
        Te += wgt * ge
        To += wgt * go
    return dict(Ae=Ae, Pe=Pe, Te=Te, Ao=Ao, Po=Po, To=To, ap=ap, am=am)


def lam(M, k=0):
    M = 0.5 * (M + M.T)
    return sla.eigvalsh(M, subset_by_index=[k, k])[0]


def lams(M, k):
    M = 0.5 * (M + M.T)
    return sla.eigvalsh(M, subset_by_index=[0, k - 1])


def kappa_interval(A, T, tol=1e-12):
    """Interval of kappa with A - kappa T >= 0 by bisection, if A - T >= 0."""
    if lam(A - T) <= 0:
        return (np.nan, np.nan)
    # upper end
    lo, hi = 1.0, 2.0
    while lam(A - hi * T) > 0 and hi < 1e6:
        lo, hi = hi, 2 * hi
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        if lam(A - mid * T) > 0:
            lo = mid
        else:
            hi = mid
    up = lo
    lo2, hi2 = 1.0, 0.0
    step = 1.0
    while lam(A - hi2 * T) > 0 and hi2 > -1e6:
        lo2, hi2 = hi2, hi2 - step
        step *= 2
    for _ in range(48):
        mid = 0.5 * (lo2 + hi2)
        if lam(A - mid * T) > 0:
            lo2 = mid
        else:
            hi2 = mid
    return (lo2, up)


if __name__ == "__main__":
    Ltop = 0.5 * np.log(3.0)
    Lbot = 0.5 * np.log(2.0)
    print("kappa_2 =", K2, " norm of T =", K2 / 2)
    print("check against the book at a=0.5: even 9.335e-7, odd 1.937e-4 (extrap.)")
    for N in (26, 160, 400, 800):
        d = matrices(0.5, N)
        print(N, lam(d["Ae"] + d["Pe"] - d["Te"]), lam(d["Ao"] + d["Po"] - d["To"]))
    print()
    print("L  e^{2L}  even: lam1 lam2 | arch only A+P: lam1 lam2 | odd: lam1 lam2 | arch only odd lam1")
    Ls = [Lbot, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, Ltop]
    N = 800
    out = []
    for L in Ls:
        d = matrices(L, N)
        Wxe = d["Ae"] + d["Pe"]
        Wxo = d["Ao"] + d["Po"]
        We = Wxe - d["Te"]
        Wo = Wxo - d["To"]
        le = lams(We, 2)
        la = lams(Wxe, 2)
        lo = lams(Wo, 2)
        lao = lams(Wxo, 2)
        print("%.5f %.4f  %.4e %.4e | %.4e %.4e | %.4e %.4e | %.4e %.4e" % (L, np.exp(2 * L), le[0], le[1], la[0], la[1], lo[0], lo[1], lao[0], lao[1]))
    print()
    print("feasible interval of the weight kappa of the prime 2 term at N=300"); N = 300
    for L in [0.36, 0.40, 0.44, 0.48, 0.50, 0.52, 0.54, Ltop]:
        d = matrices(L, N)
        ke = kappa_interval(d["Ae"] + d["Pe"], d["Te"])
        ko = kappa_interval(d["Ao"] + d["Po"], d["To"])
        print("%.5f even [%.8g, %.8g]  odd [%.8g, %.8g]" % (L, ke[0], ke[1], ko[0], ko[1]))
