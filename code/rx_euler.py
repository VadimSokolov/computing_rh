# rx_euler.py: the numerical claims of the Brownian truncated Euler product (Theorem ht:thm:euler and Table ht:tab:euler, Section 3.2 of
# misc/reviews/nick-b4/reconcile-N.md, task E1; source: the reciprocal Xi note, cited N:n). Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Claim. For U > 0 and s with Re(s^2) > -gamma_1^2,
#   xi(1/2)/xi(1/2 + s) = exp(-A_U(s)) prod_p E_{p,U}(s) exp(-R_U(s)),   E_{p,U}(s) = exp(sum_m p^(-m/2)/m B_U(m log p, s)),
# B_U(x,s) = erfc(x/sqrt(2U)) - (1/2)[e^(-xs) erfc((x - sU)/sqrt(2U)) + e^(xs) erfc((x + sU)/sqrt(2U))],
# A_U(s) = int_0^U (1 - e^(-s^2 u/2)) k_inf(u) du/u,  k_inf(u) = e^(u/8) - log(pi)/(2 sqrt(2 pi u)) + (1/2pi) int_0^inf e^(-u r^2/2) Re psi(1/4 + i r/2) dr,
# R_U(s) = sum_{gamma > 0} [E_1(gamma^2 U/2) - E_1((gamma^2 + s^2) U/2)], with the bound |R_U| <= (2/U) e^(-m_s U/2) sum (1/gamma^2 + 1/Re(gamma^2 + s^2)).
# Table ht:tab:euler rows (s, U): (0.1, 0.05), (0.1, 1), (0.49, 0.05), (0.49, 1), (2, 1), (0.3 + 3i, 0.25). Columns printed: log xi(1/2+s)/xi(1/2),
# A_U - sum_p log E_{p,U}, R_U, the bound, the residual LHS - (A_U - sum_p log E) - R_U.
# Conventions. xi(z) = z (z-1) pi^(-z/2) Gamma(z/2) zeta(z)/2; the left side is the principal logarithm of xi(1/2+s)/xi(1/2) (its imaginary part is checked small).
# The zeros are the first 1700 ordinates of data/zeros_1700.txt (RH holds there, gamma = ordinate); the rest of the zeros enter only through the bound.
# Method. (i) The three parts of A_U by closed forms: int_0^U (e^(-au) - e^(-bu)) du/u = Ein(bU) - Ein(aU) for e^(u/8) (a = -1/8, b = s^2/2 - 1/8) and
# for the digamma part (a = r^2/2, b = (r^2+s^2)/2, then mp.quad in r over [0, inf) with Ein(bU) - Ein(aU) written as log(1 + s^2/r^2) + E_1(bU) - E_1(aU)
# for r >= 6 so that no digits are lost), and int_0^U (1 - e^(-cu)) u^(-3/2) du = -2 (1 - e^(-cU))/sqrt(U) + 2 sqrt(pi c) erf(sqrt(cU)), c = s^2/2.
# (ii) The prime part sum_{n = p^m <= X} Lambda(n)/(sqrt(n) log n) B_U(log n, s): terms with (log n)^2/(2U) <= 50 in mpmath at the working precision,
# the remaining terms (each below e^-50) in double precision (numpy/scipy erfcx), whose total error is below 1e-30. The cutoff X is the smallest in a list
# for which the tail bound (1 + e^(|s|^2 U/2)) sum_{n > X} Lambda(n) n^(-1/2) (log n)^(-1) e^(-(log n)^2/(2U)) <= 1.04 X g(X) + 1.04 int_X^inf g, with
# g(y) = y^(-1/2)/log y e^(-(log y)^2/(2U)) and psi(y) <= 1.04 y (partial summation), is below 1e-34; the bound is logged. (iii) R_U from the zeros with
# E_1 of mpmath (terms with Re(gamma^2 U/2) > 700 are below 1e-300 and dropped); the rest of the zeros is bounded as in the theorem.
# Error control. Every row is computed at two working precisions (dps D1 < D2) and the largest difference is logged (it is the error of the column);
# the closed form of B_U is also checked against the defining expectation int_0^U (1 - e^(-s^2 u/2)) f_x(u) du by quadrature at x = 1.5.
# Output: data/rx_euler.json (numbers as strings) and the log rx_euler.log next to it.
# Usage: python3 rx_euler.py [--quick] [--out FILE] [--zeros FILE]; --quick uses 20 and 26 digits and primes up to 2e5 (U <= 0.25 rows complete, the
# U = 1 rows get the logged tail bound only), writes to /tmp unless --out is given; the full case (30 and 40 digits, primes up to 2e6) takes about
# 25 seconds; if it does not fit the local limit, on the cluster: bash misc/tools/hopper_run.sh -c 1 -t 30 -g rx_euler.json -g rx_euler.log code/rx_euler.py
# data/zeros_1700.txt -- --out rx_euler.json
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import mpmath as mp
from scipy import special as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = [("0.1", "0.05"), ("0.1", "1"), ("0.49", "0.05"), ("0.49", "1"), ("2", "1"), ("0.3+3j", "0.25")]
ZEROS = os.path.join(HERE, "..", "data", "zeros_1700.txt")


def sieve_lambda(X):
    """Return arrays n (prime powers <= X) and Lambda(n) = log p."""
    isp = np.ones(X + 1, dtype=bool)
    isp[:2] = False
    for q in range(2, int(X ** 0.5) + 1):
        if isp[q]:
            isp[q * q::q] = False
    pr = np.nonzero(isp)[0].astype(np.int64)
    ns, ls = [pr], [np.log(pr.astype(float))]
    pw = pr.copy()
    while True:
        pw = pw * pr
        m = pw <= X
        if not m.any():
            break
        pr, pw = pr[m], pw[m]
        ns.append(pw)
        ls.append(np.log(pr.astype(float)))
    n = np.concatenate(ns)
    lam = np.concatenate(ls)
    o = np.argsort(n)
    return n[o], lam[o]


def Ein(z):
    z = mp.mpmathify(z)
    if z == 0:
        return mp.mpf(0)
    if abs(z) <= 8:
        with mp.workdps(mp.mp.dps + 12):
            t, term, k = mp.mpf(0), z, 1
            while True:
                add = term / k
                t += add
                if abs(add) < mp.mpf(10) ** (-mp.mp.dps - 2) and k > abs(z):
                    break
                k += 1
                term = -term * z / k
            r = t
        return +r
    return mp.e1(z) + mp.log(z) + mp.euler


def Bu_mp(x, s, U):
    r2U = mp.sqrt(2 * U)
    return mp.erfc(x / r2U) - (mp.exp(-x * s) * mp.erfc((x - s * U) / r2U) + mp.exp(x * s) * mp.erfc((x + s * U) / r2U)) / 2


def Bu_np(x, s, U):
    """Double precision B_U(x, s) on an array x via erfcx (no overflow)."""
    r2U = math.sqrt(2 * U)
    z0 = x / r2U
    z1 = (x - s * U) / r2U
    z2 = (x + s * U) / r2U
    t0 = np.exp(-z0 * z0) * sp.erfcx(z0)
    t1 = np.exp(-x * s - z1 * z1) * sp.erfcx(z1)
    t2 = np.exp(x * s - z2 * z2) * sp.erfcx(z2)
    return t0 - (t1 + t2) / 2


def tail_bound(X, U, s):
    """Rigorous bound for the omitted primes n > X: (1+e^{|s|^2 U/2}) (1.04 X g(X) + 1.04 int_X^inf g)."""
    X = mp.mpf(X)
    L = mp.log(X)
    g = lambda y: 1 / (mp.sqrt(y) * mp.log(y)) * mp.exp(-mp.log(y) ** 2 / (2 * U))
    # int_X^inf g dy <= (1/log X) int_L^inf e^{x/2 - x^2/(2U)} dx = (1/L) sqrt(pi U/2) e^{U/8} erfc((L - U/2)/sqrt(2U))
    ig = mp.sqrt(mp.pi * U / 2) * mp.exp(U / 8) * mp.erfc((L - U / 2) / mp.sqrt(2 * U)) / L
    return (1 + mp.exp(abs(s) ** 2 * U / 2)) * mp.mpf("1.04") * (X * g(X) + ig)


def compute_row(args):
    srow, Urow, dps, Xmax, zeros_path, quick = args
    mp.mp.dps = dps
    s = mp.mpmathify(srow)
    U = mp.mpf(Urow)
    t0 = time.time()
    # left side
    xi = lambda z: z * (z - 1) * mp.pi ** (-z / 2) * mp.gamma(z / 2) * mp.zeta(z) / 2
    half = mp.mpf(1) / 2
    lhs = mp.log(xi(half + s) / xi(half))
    # A_U
    a1, b1 = -mp.mpf(1) / 8, s * s / 2 - mp.mpf(1) / 8
    T1 = Ein(b1 * U) - Ein(a1 * U)
    c = s * s / 2
    J = -2 * (1 - mp.exp(-c * U)) / mp.sqrt(U) + 2 * mp.sqrt(mp.pi * c) * mp.erf(mp.sqrt(c * U))
    T2 = -mp.log(mp.pi) / (2 * mp.sqrt(2 * mp.pi)) * J

    def f3(r):
        if r < 6:
            br = Ein((r * r + s * s) / 2 * U) - Ein(r * r / 2 * U)
        else:
            br = mp.log1p(s * s / (r * r)) + mp.e1((r * r + s * s) / 2 * U) - mp.e1(r * r / 2 * U)
        return mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2)) * br

    T3, e3 = mp.quad(f3, [0, 1, 3, 6, 15, 50, 300, 3000, mp.inf], error=True)
    T3 = T3 / (2 * mp.pi)
    AU = T1 + T2 + T3
    # primes
    n, lam = sieve_lambda(Xmax)
    Xc = Xmax
    for cand in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 2 * 10 ** 6, 5 * 10 ** 6):
        if cand <= Xmax and tail_bound(cand, U, s) < mp.mpf(10) ** (-34):
            Xc = cand
            break
    tb = tail_bound(Xc, U, s)
    m = n <= Xc
    n, lam = n[m], lam[m]
    x = np.log(n.astype(float))
    thr = x * x / (2 * float(U)) <= 50
    S = mp.mpf(0)
    nn = n[thr]
    for ni, li in zip(nn, lam[thr]):
        xi_ = mp.log(int(ni))
        S += mp.log(int(round(math.exp(li)))) / (mp.sqrt(int(ni)) * xi_) * Bu_mp(xi_, s, U)
    rest = ~thr
    xr = x[rest]
    wr = lam[rest] / (np.sqrt(n[rest].astype(float)) * xr)
    sc = complex(s)
    Sr = np.sum(wr * Bu_np(xr, sc, float(U)))
    S = S + mp.mpmathify(complex(Sr))
    arith = AU - S
    # zeros
    gam = [mp.mpf(line.strip()) for line in open(zeros_path) if line.strip()]
    R = mp.mpf(0)
    for g in gam:
        a = g * g * U / 2
        if a > 700:
            break
        R += mp.e1(a) - mp.e1((g * g + s * s) * U / 2)
    ms = gam[0] ** 2 + min(0, mp.re(s * s))
    K = len(gam)
    sumz = sum(1 / g ** 2 + 1 / mp.re(g * g + s * s) for g in gam)
    gK = gam[-1]
    tailz = 2 * ((mp.log(gK / (2 * mp.pi)) + 1) / (2 * mp.pi * gK)) * mp.mpf("1.1")  # Riemann-von Mangoldt average, 10% margin
    bound = 2 / U * mp.exp(-ms * U / 2) * (sumz + tailz)
    resid = lhs - arith - R
    return dict(s=srow, U=Urow, dps=dps, X=Xc, lhs=lhs, arith=arith, R=R, bound=bound, resid=resid, tailbound=tb, quad_err=e3,
                imlhs=mp.im(lhs), secs=time.time() - t0, Kzeros=K, quick=quick)


def check_B(dps):
    mp.mp.dps = dps
    s, U, x = mp.mpmathify("0.3+3j"), mp.mpf("0.25"), mp.mpf("1.5")
    f = lambda u: (1 - mp.exp(-s * s * u / 2)) * x * mp.exp(-x * x / (2 * u)) / mp.sqrt(2 * mp.pi * u ** 3)
    q = mp.quad(f, [0, U / 8, U / 2, U])
    return Bu_mp(x, s, U), q


def fmt(z, d=17):
    z = mp.mpmathify(z)
    if isinstance(z, mp.mpc) and mp.im(z) != 0:
        return mp.nstr(mp.re(z), d) + (" + " if mp.im(z) >= 0 else " - ") + mp.nstr(abs(mp.im(z)), d) + "i"
    return mp.nstr(mp.re(z), d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--zeros", default=ZEROS)
    a = ap.parse_args()
    d1, d2, Xmax = (20, 26, 200000) if a.quick else (30, 40, 2000000)
    out = a.out or ("/tmp/rx_euler.json" if a.quick else os.path.join(HERE, "..", "data", "rx_euler.json"))
    logp = os.path.splitext(out)[0] + ".log"
    log = open(logp, "w")

    def P(*t):
        line = " ".join(str(v) for v in t)
        print(line)
        log.write(line + "\n")
        log.flush()

    t0 = time.time()
    jobs = [(sr, Ur, d, Xmax, a.zeros, a.quick) for (sr, Ur) in ROWS for d in (d1, d2)]
    with Pool(min(len(jobs), os.cpu_count() or 1)) as pool:
        res = pool.map(compute_row, jobs)
    P("rx_euler.py: Theorem ht:thm:euler, Table ht:tab:euler; dps", d1, d2, "primes up to", Xmax, "quick" if a.quick else "full")
    data = []
    for i, (sr, Ur) in enumerate(ROWS):
        r1, r2 = res[2 * i], res[2 * i + 1]
        dif = {k: abs(r1[k] - r2[k]) for k in ("lhs", "arith", "R", "resid")}
        P("row s =", sr, "U =", Ur, "primes to", r2["X"], "tail bound", mp.nstr(r2["tailbound"], 3), "quad err", mp.nstr(r2["quad_err"], 3), "secs", round(r2["secs"], 1))
        P("  log xi(1/2+s)/xi(1/2) =", fmt(r2["lhs"]), "(Im part before: ", mp.nstr(r2["imlhs"], 5), ")")
        P("  A_U - sum_p log E     =", fmt(r2["arith"]))
        P("  R_U                   =", fmt(r2["R"], 6), " bound of the theorem =", mp.nstr(r2["bound"], 4))
        P("  residual lhs-(arith+R) =", mp.nstr(abs(r2["resid"]), 4), " | difference of the two precisions: lhs", mp.nstr(dif["lhs"], 3),
          "arith", mp.nstr(dif["arith"], 3), "R", mp.nstr(dif["R"], 3), "resid", mp.nstr(dif["resid"], 3))
        data.append(dict(s=sr, U=Ur, X=r2["X"], lhs=fmt(r2["lhs"], 25), arith=fmt(r2["arith"], 25), R=fmt(r2["R"], 10), bound=mp.nstr(r2["bound"], 6),
                         resid=mp.nstr(abs(r2["resid"]), 4), tailbound=mp.nstr(r2["tailbound"], 4), prec_diff_arith=mp.nstr(dif["arith"], 3),
                         prec_diff_lhs=mp.nstr(dif["lhs"], 3)))
    cf, q = check_B(d1)
    P("closed form of B_U at x=1.5, s=0.3+3i, U=0.25:", fmt(cf, 20), " expectation by quad:", fmt(q, 20), " |difference|", mp.nstr(abs(cf - q), 3))
    P("total seconds", round(time.time() - t0, 1))
    json.dump(dict(rows=data, B_check_diff=mp.nstr(abs(cf - q), 3), dps=[d1, d2], primes_to=Xmax), open(out, "w"), indent=1)
    log.close()


if __name__ == "__main__":
    main()
