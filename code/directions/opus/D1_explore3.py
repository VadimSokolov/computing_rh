#!/usr/bin/env python3
"""D1 exploration, part 3 (double precision, NOT a certificate): thresholds of the
archimedean form, weight interval of the prime 2, reduced form in the Legendre basis,
and the Davenport Heilbronn rung. Sections chosen by argv: thr wei red dh.
Usage: OMP_NUM_THREADS=1 timeout 500 python3 -I -u D1_explore3.py thr wei red dh
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.linalg as sla
from scipy.special import spherical_jn, digamma, sici
from numpy.polynomial.legendre import leggauss
import D1_explore as E
Ltop = 0.5 * np.log(3.0)
Lbot = 0.5 * np.log(2.0)
K2 = E.K2
SEC = set(sys.argv[1:])
def bis(fun, a, b, it=30):
    fa = fun(a)
    for _ in range(it):
        c = 0.5 * (a + b)
        fc = fun(c)
        if (fc > 0) == (fa > 0):
            a, fa = c, fc
        else:
            b = c
    return 0.5 * (a + b)

def _A(L, sec, k=0):
    d = E.matrices(L, 200)
    return E.lam(d["A" + sec] + d["P" + sec], k)
if "thr" in SEC:
    print("=== thresholds of mu1(A+P)=0 and of mu2(A+P)=|T|, N=200")
    print("even mu1=0 at L=", bis(lambda L: _A(L, "e"), 0.38, 0.45, 16))
    print("odd  mu1=0 at L=", bis(lambda L: _A(L, "o"), 0.3466, 0.40, 16))
    print("even mu2=|T| at L=", bis(lambda L: _A(L, "e", 1) - K2 / 2, 0.36, 0.5493, 14))
    print("odd  mu2=|T| at L=", bis(lambda L: _A(L, "o", 1) - K2 / 2, 0.36, 0.5493, 14))
if "wei" in SEC:
    print("=== weight interval N=200: L even[low,up] odd[low,up]")
    for L in [0.36, 0.40, 0.44, 0.48, 0.50, 0.52, 0.54, Ltop]:
        d = E.matrices(L, 200)
        ke = E.kappa_interval(d["Ae"] + d["Pe"], d["Te"])
        ko = E.kappa_interval(d["Ao"] + d["Po"], d["To"])
        print("%.5f even [%.9g, %.9g] odd [%.9g, %.9g]" % (L, ke[0], ke[1], ko[0], ko[1]))
def mult(t):
    return digamma(0.25 + 0.5j * t).real - np.log(np.pi)


def reduced(L, T, Nmax, panels_per_unit=1, ng=40):
    """Least eigenvalues of beta I + pole + band form in the Legendre basis, both sectors."""
    tg, wg = leggauss(ng)
    npan = int(np.ceil(T * panels_per_unit))
    edges = np.linspace(0, T, npan + 1)
    ts = np.concatenate([0.5 * (edges[i + 1] - edges[i]) * tg + 0.5 * (edges[i + 1] + edges[i]) for i in range(npan)])
    ws = np.concatenate([0.5 * (edges[i + 1] - edges[i]) * wg for i in range(npan)])
    beta = np.log(T / (2 * np.pi)) - 1.2 / T ** 2 - K2
    psi = mult(ts) - K2 * np.cos(ts * np.log(2.0))
    q = np.array([np.sqrt(2 * L * (2 * n + 1)) * spherical_jn(n, L * ts) for n in range(Nmax + 1)])
    # modified spherical Bessel at L/2 via series
    def inb(n, z):
        s = 0.0
        term = z ** n
        for k in range(1, n + 1):
            term /= (2 * k + 1)
        s = term
        for k in range(1, 40):
            term *= (z * z / 2) / (k * (2 * n + 2 * k + 1))
            s += term
        return s
    a = np.array([np.sqrt(2 * L * (2 * n + 1)) * inb(n, L / 2) for n in range(Nmax + 1)])
    res = {}
    for sec, idx, sg in (("e", np.arange(0, Nmax + 1, 2), 2.0), ("o", np.arange(1, Nmax + 1, 2), -2.0)):
        s = (-1.0) ** (idx // 2)
        Q = (q[idx].T * s).T
        K = (Q * (ws * (psi - beta) / np.pi)) @ Q.T
        an = a[idx]
        M = beta * np.eye(len(idx)) + K + sg * np.outer(an, an)
        res[sec] = sla.eigvalsh(0.5 * (M + M.T), subset_by_index=[0, 1])
    return beta, res

if "red" in SEC:
    print("=== reduced form (lower bound) in the Legendre basis at L*")
    for T in (20, 25, 30, 40, 50, 70, 100, 140, 200, 300, 400, 600):
        Nmax = int(1.5 * Ltop * T * np.e / 2) + 30
        beta, r = reduced(Ltop, T, Nmax)
        print("T %4d orders<=%4d beta %.4f even %.4e %.4e odd %.4e %.4e" % (T, Nmax, beta, r["e"][0], r["e"][1], r["o"][0], r["o"][1]))
    print("=== reduced form at a=0.5")
    for T in (50, 100, 200, 400):
        Nmax = int(1.5 * 0.5 * T * np.e / 2) + 30
        beta, r = reduced(0.5, T, Nmax)
        print("T %4d even %.4e odd %.4e" % (T, r["e"][0], r["o"][0]))
def dh_matrices(L, N):
    """Arch form with kernel exp(-x/2)/sinh x, constant log(5/pi), no pole, weight sqrt2 kap log2 at n=2."""
    nq = max(4000, 8 * N)
    t, w = E._gl(nq)
    x = L * (t + 1.0)
    w = L * w
    h = np.exp(-x / 2) / np.sinh(x) - 1.0 / x
    n = np.arange(N + 1, dtype=float)
    om = np.pi * n / L
    tri = 1 - x / (2 * L)
    In = np.zeros(N + 1)
    Cn = np.zeros(N + 1)
    for start in range(0, N + 1, 256):
        idx = np.arange(start, min(start + 256, N + 1))
        ph = np.outer(om[idx], x)
        Cn[idx] = (np.cos(ph) * (h * tri * w)).sum(axis=1)
        In[idx] = (np.sin(ph) * (h * w)).sum(axis=1)
    si, _ = sici(2 * np.pi * n[1:])
    In[1:] += si
    _, ci = sici(2 * np.pi * n[1:])
    base = np.log(5 * n[1:] / (2 * L)) - ci - Cn[1:]
    Ae = np.zeros((N + 1, N + 1))
    Ao = np.zeros((N, N))
    Ae[0, 0] = -E.GAMMA + np.log(5 / np.pi) - np.log(4 * L) + 1 - Cn[0]
    Ae[np.arange(1, N + 1), np.arange(1, N + 1)] = base + In[1:] / (2 * np.pi * n[1:])
    Ao[np.arange(N), np.arange(N)] = base - In[1:] / (2 * np.pi * n[1:])
    m_, n_ = np.meshgrid(n[1:], n[1:], indexing="ij")
    sgn = (-1.0) ** (m_ + n_)
    denom = np.pi * (m_ ** 2 - n_ ** 2)
    np.fill_diagonal(denom, 1.0)
    I_m = np.outer(In[1:], np.ones(N))
    I_n = I_m.T
    oe = -sgn * (n_ * I_n - m_ * I_m) / denom
    oo = -sgn * (m_ * I_n - n_ * I_m) / denom
    np.fill_diagonal(oe, 0)
    np.fill_diagonal(oo, 0)
    Ae[1:, 1:] += oe
    Ao += oo
    v = ((-1.0) ** n[1:]) * In[1:] / (np.sqrt(2.0) * np.pi * n[1:])
    Ae[0, 1:] = v
    Ae[1:, 0] = v
    return Ae, Ao

if "dh" in SEC:
    print("=== Davenport Heilbronn rung, cos modes N=300")
    kap = (np.sqrt(10 - 2 * np.sqrt(5)) - 2) / (np.sqrt(5) - 1)
    print("kappa_DH", kap, "weight of n=2:", np.sqrt(2) * kap * np.log(2))
    zeta = E.matrices(Ltop, 300)
    for L in [0.2, 0.3, Lbot, 0.4, 0.5, Ltop]:
        Ae, Ao = dh_matrices(L, 300)
        cl = [(np.log(2.0), np.sqrt(2) * kap * np.log(2.0))]
        d = E.matrices(L, 300, cl) if L > Lbot else None
        Te = d["Te"] if d else 0 * Ae
        To = d["To"] if d else 0 * Ao
        print("DH L=%.4f even lam1 %.4e arch-only %.4e | odd lam1 %.4e arch-only %.4e" % (L, E.lam(Ae - Te), E.lam(Ae), E.lam(Ao - To), E.lam(Ao)))
