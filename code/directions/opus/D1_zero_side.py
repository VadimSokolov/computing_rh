#!/usr/bin/env python3
"""D1 zero side pilot in double precision. NOT a certificate.

(1) Checks the identity of D1.md, Section 3.7, on the window L=0.5 with a cut-off a(t).
(2) Shows which zeros carry the least eigenvalue at the top of the rung.
(3) Least eigenvalue of the zero side lower form, as a function of the number K of zeros kept.
(4) Last negative frequency of the symbol on the rung; norm of the prime operator against A_L.

Reads only the ordinates in data/zeros_8000.txt. Usage: timeout 500 python3 -I -u D1_zero_side.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.linalg as sla
from scipy.special import digamma, spherical_jn, roots_legendre
import D1_explore as E

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))
GAM = np.loadtxt(os.path.join(ROOT, "data", "zeros_8000.txt"))
K2 = E.K2
LOG2 = np.log(2.0)
Ltop = 0.5 * np.log(3.0)
Lbot = 0.5 * np.log(2.0)


def Fcos(t, c, L):
    """Transform of sum c_n e_n (even sector basis), real, at real t > 0."""
    t = np.asarray(t, dtype=float)
    N = len(c) - 1
    om = np.pi * np.arange(N + 1) / L
    out = np.zeros_like(t)
    for start in range(0, t.size, 4096):
        tt = t[start:start + 4096]
        den = tt[:, None] ** 2 - om[None, :] ** 2
        den[np.abs(den) < 1e-9] = 1e-9
        sgn = (-1.0) ** np.arange(N + 1)
        fac = 2 * tt[:, None] * sgn / den
        fac[:, 0] *= 1 / np.sqrt(2.0)
        out[start:start + 4096] = np.sin(tt * L) * (fac @ c) / np.sqrt(L)
    return out


def Fsin(t, c, L):
    """|transform| up to sign of sum c_n s_n (odd sector basis), n=1.., at real t > 0."""
    t = np.asarray(t, dtype=float)
    N = len(c)
    om = np.pi * np.arange(1, N + 1) / L
    out = np.zeros_like(t)
    sgn = (-1.0) ** np.arange(1, N + 1)
    for start in range(0, t.size, 4096):
        tt = t[start:start + 4096]
        den = tt[:, None] ** 2 - om[None, :] ** 2
        den[np.abs(den) < 1e-9] = 1e-9
        fac = 2 * om[None, :] * sgn / den
        out[start:start + 4096] = np.sin(tt * L) * (fac @ c) / np.sqrt(L)
    return out


def symbol(t):
    return digamma(0.25 + 0.5j * t).real - np.log(np.pi) - K2 * np.cos(t * LOG2)


def nodes(a, b, width, ng):
    tg, wg = roots_legendre(ng)
    n = int(np.ceil((b - a) / width))
    ed = np.linspace(a, b, n + 1)
    t = np.concatenate([0.5 * (ed[i + 1] - ed[i]) * tg + 0.5 * (ed[i + 1] + ed[i]) for i in range(n)])
    w = np.concatenate([0.5 * (ed[i + 1] - ed[i]) * wg for i in range(n)])
    return t, w


def cut(t, eta, r):
    return np.sinc(eta * t / (r * np.pi)) ** r


def cut_i2(eta, r):
    v = eta / (2 * r)
    return (np.sinh(v) / v) ** r


print("=== (1) identity on L=0.5, eta=0.04, r=8; 8000 zeros")
L, eta, r = 0.5, 0.04, 8
d = E.matrices(L, 200)
Me = d["Ae"] + d["Pe"] - d["Te"]
Mo = d["Ao"] + d["Po"] - d["To"]
we, Ve = sla.eigh(0.5 * (Me + Me.T))
wo, Vo = sla.eigh(0.5 * (Mo + Mo.T))
rng = np.random.default_rng(1)
crand = rng.standard_normal(201) / (1 + np.arange(201.0)) ** 2
crand /= np.linalg.norm(crand)
ts, ws = nodes(0.0, 12000.0, 1.0, 16)
tests = [("even ground", Ve[:, 0], "e"), ("odd ground", Vo[:, 0], "o"), ("even random", crand, "e")]
for name, c, sec in tests:
    if sec == "e":
        Ft = Fcos(ts, c, L)
        Fz = Fcos(GAM, c, L)
        pole = 2 * (d["ap"] @ c) ** 2
        Wf = c @ Me @ c
    else:
        Ft = Fsin(ts, c, L)
        Fz = Fsin(GAM, c, L)
        pole = -2 * (d["am"] @ c) ** 2
        Wf = c @ Mo @ c
    a2t = cut(ts, eta, r) ** 2
    a2z = cut(GAM, eta, r) ** 2
    Sz = 2 * np.sum(a2z * Fz ** 2)
    Sd = cut_i2(eta, r) ** 2 * pole + np.sum(ws * a2t * Ft ** 2 * symbol(ts)) / np.pi
    high = Wf - pole - np.sum(ws * a2t * Ft ** 2 * symbol(ts)) / np.pi  # (1/pi) int (1-a^2) F^2 Psi
    rhs = Sz + high + (1 - cut_i2(eta, r) ** 2) * pole
    print("%-12s W(f)=%.10e zeros=%.10e direct=%.10e |zeros-direct|=%.1e  identity residual=%.1e  high part=%.3e" % (name, Wf, Sz, Sd, abs(Sz - Sd), abs(rhs - Wf), high))

print("=== (2) share of the least eigenvalue by zero at L* (N=300, no cut-off)")
d = E.matrices(Ltop, 300)
Me = d["Ae"] + d["Pe"] - d["Te"]
Mo = d["Ao"] + d["Po"] - d["To"]
we, Ve = sla.eigh(0.5 * (Me + Me.T))
wo, Vo = sla.eigh(0.5 * (Mo + Mo.T))
Fe = 2 * Fcos(GAM, Ve[:, 0], Ltop) ** 2
Fo = 2 * Fsin(GAM, Vo[:, 0], Ltop) ** 2
print("even lam1=%.5e first terms/lam1: %s" % (we[0], np.array2string(Fe[:6] / we[0], precision=4)))
print("odd  lam1=%.5e first terms/lam1: %s" % (wo[0], np.array2string(Fo[:6] / wo[0], precision=4)))
for K in (1, 2, 3, 5, 10, 20, 50, 100, 300, 1000, 8000):
    print("K=%5d gamma_K=%9.3f even partial/lam1=%.5f odd partial/lam1=%.5f" % (K, GAM[K - 1], Fe[:K].sum() / we[0], Fo[:K].sum() / wo[0]))
te, wte = nodes(GAM[-1], 4e6, 4.0, 4)
tail_e = np.sum(wte * Fcos(te, Ve[:, 0], Ltop) ** 2 * np.log(te / (2 * np.pi))) / np.pi
print("mean density tail beyond gamma_8000 even/lam1 about %.4f" % (tail_e / we[0]))

print("=== (3) zero side lower form in the Legendre basis, modes to order 110, T=104, beta=m(T)-K2")


def legendre_q(t, L, Nmax):
    return np.array([np.sqrt(2 * L * (2 * n + 1)) * spherical_jn(n, L * t) for n in range(Nmax + 1)])


def inb(n, z):
    term = z ** n
    for k in range(1, n + 1):
        term /= 2 * k + 1
    s = term
    for k in range(1, 40):
        term *= (z * z / 2) / (k * (2 * n + 2 * k + 1))
        s += term
    return s


def zero_form(L, eta, r, Klist, T=104, Nmax=110):
    t1, w1 = nodes(0.0, float(T), 1.0, 40)
    t2, w2 = nodes(float(T), 12000.0, 4.0, 24)
    q1 = legendre_q(t1, L, Nmax)
    q2 = legendre_q(t2, L, Nmax)
    qz = legendre_q(GAM, L, Nmax)
    beta = digamma(0.25 + 0.5j * T).real - np.log(np.pi) - K2
    a_vec = np.array([np.sqrt(2 * L * (2 * n + 1)) * inb(n, L / 2) for n in range(Nmax + 1)])
    a2i = cut_i2(eta, r) ** 2
    out = {}
    for sec, ix, sg in (("e", np.arange(0, Nmax + 1, 2), 2.0), ("o", np.arange(1, Nmax + 1, 2), -2.0)):
        s = (-1.0) ** (ix // 2)
        Q1 = (q1[ix].T * s).T
        Q2 = (q2[ix].T * s).T
        Qz = (qz[ix].T * s).T
        a1 = cut(t1, eta, r) ** 2
        a2 = cut(t2, eta, r) ** 2
        az = cut(GAM, eta, r) ** 2
        base = beta * np.eye(len(ix))
        base += (Q1 * (w1 * ((1 - a1) * symbol(t1) - beta) / np.pi)) @ Q1.T
        base -= (Q2 * (w2 * beta * a2 / np.pi)) @ Q2.T
        an = a_vec[ix]
        base += (1 - a2i) * sg * np.outer(an, an)
        # direct reduced form at the same T for comparison
        direct = beta * np.eye(len(ix)) + (Q1 * (w1 * (symbol(t1) - beta) / np.pi)) @ Q1.T + sg * np.outer(an, an)
        ldirect = sla.eigvalsh(0.5 * (direct + direct.T), subset_by_index=[0, 0])[0]
        res = []
        for K in Klist:
            Z = (Qz[:, :K] * (2 * az[:K])) @ Qz[:, :K].T
            M = base + Z
            res.append(sla.eigvalsh(0.5 * (M + M.T), subset_by_index=[0, 0])[0])
        out[sec] = (ldirect, res)
    return out


Klist = [5, 10, 20, 40, 80, 150, 300, 600, 1000, 2000, 4000, 8000]
for (L, eta, r) in ((0.5, 0.04, 8), (0.54, 0.004, 8)):
    z = zero_form(L, eta, r, Klist)
    print("L=%.3f eta=%g r=%d" % (L, eta, r))
    for sec in "eo":
        print("  sector %s direct reduced form lam=%.4e; zero side by K: %s" % (sec, z[sec][0], " ".join("%d:%.3e" % (K, v) for K, v in zip(Klist, z[sec][1]))))

print("=== (4) last negative frequency of Psi on the rung; prime operator norm against A_L")
t = np.linspace(0.0, 40.0, 400001)
psi = symbol(t)
neg = t[psi < 0]
print("Psi<0 only on t<= %.4f ; min Psi = %.4f at t=%.3f; connected pieces end at %s" % (neg.max(), psi.min(), t[np.argmin(psi)], np.array2string(t[1:][(psi[:-1] < 0) & (psi[1:] >= 0)], precision=4)))


def comb_list(L):
    out = []
    n = 2
    while np.log(n) < 2 * L:
        q = n
        p = None
        for f in range(2, n + 1):
            if q % f == 0:
                p = f
                while q % f == 0:
                    q //= f
                break
        if q == 1:
            out.append((np.log(n), 2 * np.log(p) / np.sqrt(n)))
        n += 1
    return out


for L in (0.5493, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0):
    cl = comb_list(L)
    A_L = sum(w for _, w in cl)
    # sup V
    xs = np.linspace(0, L, 20001)
    V = np.array([0.5 * sum(w * (c <= L + x) + w * (c <= L - x) for c, w in cl) for x in xs[::20]])
    d = E.matrices(L, 1000, cl)
    nT = max(abs(sla.eigvalsh(d["Te"], subset_by_index=[0, 0])[0]), abs(sla.eigvalsh(d["Te"], subset_by_index=[1000, 1000])[0]), abs(sla.eigvalsh(d["To"], subset_by_index=[0, 0])[0]), abs(sla.eigvalsh(d["To"], subset_by_index=[999, 999])[0]))
    print("L=%.4f primes<=%d A_L=%.4f A_L/2=%.4f sup V=%.4f |T| at N=1000=%.4f  log10 T1(A_L)=%.2f log10 T1(supV)=%.2f log10 T1(|T|)=%.2f" % (L, int(np.exp(2 * L)), A_L, A_L / 2, V.max(), nT, np.log10(2 * np.pi) + A_L / np.log(10), np.log10(2 * np.pi) + V.max() / np.log(10), np.log10(2 * np.pi) + nT / np.log(10)))
