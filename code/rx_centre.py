# rx_centre.py: the numerical claims of Remark bp:rem:centre (the centre alpha = 1/2: finite part F(1/2) and the identity at every t > 0; Section 3.10
# of misc/reviews/nick-b4/reconcile-N.md, task E2; source: the reciprocal Xi note, cited N:n). Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Claim. With W(t) = sum_{gamma > 0} e^(-t gamma^2) (zeros 1/2 + i gamma), F(1/2) = 2 - zeta'/zeta(1/2) = 2 - (1/2) log pi + (1/2) psi(1/4), and
# Abar(t) = 2 sqrt(pi t) pi^(-1) int_0^inf e^(-t r^2) sum_{k >= 1} q_{1/2 + 2k}(r) dr, q_c(r) = r^2/(c (c^2 + r^2)),
#   2 sqrt(pi t) W(t) = Abar(t) + F(1/2) - int_[1,inf) e^(-(log y)^2/4t) y^(-1/2) d(psi(y) - y),     psi = Chebyshev function,
# the last integral being sum_n Lambda(n) n^(-1/2) e^(-(log n)^2/4t) - int_1^inf e^(-(log y)^2/4t) y^(-1/2) dy (both divergent separately in the limit X -> inf
# only through the pole, the difference converges). Rows: t = 0.01, 0.025, 0.125, 0.5, 2.
# Conventions. Left side from the first 1700 ordinates of data/zeros_1700.txt (the rest is below the logged bound); Lambda(n) = log p for n = p^m.
# Method. F(1/2) from both closed forms (mp.zeta derivative, mp.digamma). Abar(t) three ways: (a) mp.quad of the closed form
# sum_{k >= 1} q_{1/2+2k}(r) = (1/2)[Re psi(1/4 + i r/2) - psi(1/4)] - 2 + (1/2)/(1/4 + r^2) (the closed form is also compared with the series at r = 1.3),
# (b) the termwise integrals int_0^inf e^(-t r^2) q_c(r) dr = (1/c)[(1/2) sqrt(pi/t) - (pi c/2) e^(t c^2) erfc(c sqrt t)] summed over c = 1/2 + 2k (directly up to k = 500, then the tail from the asymptotic series of the terms summed with Hurwitz zeta functions),
# (c) the same sum with a second splitting point (directly summed up to k = 1000 instead of 500); the tail beyond is the same asymptotic series. The compensator integral up to X: int_1^X e^(-(log y)^2/4t) y^(-1/2) dy = e^(t/4) sqrt(pi t)
# [erf((log X - t)/(2 sqrt t)) + erf(sqrt t/2)]. The prime sum over prime powers n <= X = 1e7 in double precision (math.fsum of numpy terms; rounding error ~1e-15),
# checked against mpmath at 25 digits for all n <= 2e5. Omitted tail beyond X: by partial summation with psi(y) <= 1.04 y,
# sum_{n > X} Lambda(n) h(n) <= 1.04 (X h(X) + int_X^inf h), h(y) = y^(-1/2) e^(-(log y)^2/4t) (decreasing for log y >= t), and the omitted part of the compensator is the exact
# erfc integral; the sum of the two is the logged tail error of the difference. Do not use t beyond about 3 (the integrand peaks near log y = t).
# Error control. Everything but the prime sum is computed at two precisions (30 and 40 digits) and the largest difference is logged; the prime sum is also
# computed with X = 1e6 and X = 1e7 and the change is logged.
# Output: data/rx_centre.json (numbers as strings) and the log rx_centre.log next to it.
# Usage: python3 rx_centre.py [--quick] [--out FILE] [--zeros FILE]; --quick uses X = 1e6 and 20/26 digits and writes to /tmp unless --out is given; the full case
# takes about 15 seconds; on the cluster: bash misc/tools/hopper_run.sh -c 1 -t 30 -g rx_centre.json -g rx_centre.log code/rx_centre.py data/zeros_1700.txt -- --out rx_centre.json
import argparse
import json
import math
import os
import time

import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ZEROS = os.path.join(HERE, "..", "data", "zeros_1700.txt")
TS = ["0.01", "0.025", "0.125", "0.5", "2"]


def prime_powers(X):
    isp = np.ones(X + 1, dtype=bool)
    isp[:2] = False
    for q in range(2, int(X ** 0.5) + 1):
        if isp[q]:
            isp[q * q::q] = False
    pr = np.nonzero(isp)[0].astype(np.int64)
    ns, ls = [pr], [np.log(pr.astype(float))]
    pw, p0 = pr.copy(), pr.copy()
    while True:
        pw = pw * p0
        m = pw <= X
        if not m.any():
            break
        p0, pw = p0[m], pw[m]
        ns.append(pw)
        ls.append(np.log(p0.astype(float)))
    n = np.concatenate(ns)
    lam = np.concatenate(ls)
    o = np.argsort(n)
    return n[o], lam[o]


def prime_sum_double(n, lam, t, X):
    m = n <= X
    nn = n[m].astype(float)
    x = np.log(nn)
    terms = lam[m] / np.sqrt(nn) * np.exp(-x * x / (4 * t))
    return math.fsum(terms.tolist())


def comp_int(t, X):
    L = mp.log(X)
    return mp.exp(t / 4) * mp.sqrt(mp.pi * t) * (mp.erf((L - t) / (2 * mp.sqrt(t))) + mp.erf(mp.sqrt(t) / 2))


def tail_pieces(t, X):
    L = mp.log(X)
    h = mp.exp(-L / 2 - L * L / (4 * t))
    ih = mp.exp(t / 4) * mp.sqrt(mp.pi * t) * mp.erfc((L - t) / (2 * mp.sqrt(t)))
    return mp.mpf("1.04") * (X * h + ih), ih


def Abar_all(t, dps):
    mp.mp.dps = dps
    t = mp.mpf(t)
    psi14 = mp.digamma(mp.mpf(1) / 4)
    sumq = lambda r: (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2)) - psi14) / 2 - 2 + mp.mpf(1) / 2 / (mp.mpf(1) / 4 + r * r)
    # (a) quadrature of the closed form
    pts = [mp.mpf(0), 1, 3, 8, 20] + [k / mp.sqrt(t) for k in (3, 6, 10, 16) if k / mp.sqrt(t) > 20]
    ia = mp.quad(lambda r: mp.exp(-t * r * r) * sumq(r), pts + [mp.inf])
    Aa = 2 * mp.sqrt(mp.pi * t) / mp.pi * ia

    # (b), (c) termwise integrals; for x = c sqrt(t) >= 12 the bracket is the asymptotic series of e^(x^2) erfc(x), which avoids the 1/c^2 cancellation:
    # (pi c/2) e^(t c^2) erfc(c sqrt t) = (1/2) sqrt(pi/t) sum_{n >= 0} (-1)^n (2n-1)!!/(2x^2)^n, with a smallest term below e^(-x^2)
    def term(k):
        c = mp.mpf(1) / 2 + 2 * k
        x = c * mp.sqrt(t)
        if x >= 12:
            y = 1 / (2 * x * x)
            tot, tm = mp.mpf(0), mp.mpf(1)
            for nn in range(1, 60):
                tm = -tm * (2 * nn - 1) * y
                tot -= tm
                if abs(tm) < mp.mpf(10) ** (-dps - 5):
                    break
            return mp.sqrt(mp.pi / t) / 2 * tot / c
        return (mp.sqrt(mp.pi / t) / 2 - mp.pi * c / 2 * mp.exp(t * c * c) * mp.erfc(x)) / c

    def tailsum(K):
        # sum_{k >= K} term(k) from the asymptotic series of term: with c = 1/2 + 2k, term = (1/2) sqrt(pi/t) sum_{n >= 1} (-1)^(n+1) (2n-1)!!/(2t)^n c^(-2n-1),
        # and sum_{k >= K} c^(-m) = 2^(-m) zeta(m, K + 1/4) (Hurwitz); valid since x = c sqrt(t) >= 12 at k = K, and the series is stopped at 10^(-dps-5) relative
        tot, df = mp.mpf(0), mp.mpf(1)
        for nn in range(1, 80):
            df = df * (2 * nn - 1)
            m = 2 * nn + 1
            piece = (-1) ** (nn + 1) * df / (2 * t) ** nn * mp.zeta(m, K + mp.mpf(1) / 4) / mp.mpf(2) ** m
            tot += piece
            if abs(piece) < abs(tot) * mp.mpf(10) ** (-dps - 5):
                break
        return mp.sqrt(mp.pi / t) / 2 * tot

    sb = mp.fsum(term(k) for k in range(1, 500)) + tailsum(500)
    Ab = 2 * mp.sqrt(mp.pi * t) / mp.pi * sb
    sc = mp.fsum(term(k) for k in range(1, 1000)) + tailsum(1000)
    Ac = 2 * mp.sqrt(mp.pi * t) / mp.pi * sc
    return Aa, Ab, Ac


def check_sumq(dps):
    mp.mp.dps = dps
    r = mp.mpf("1.3")
    psi14 = mp.digamma(mp.mpf(1) / 4)
    closed = (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2)) - psi14) / 2 - 2 + mp.mpf(1) / 2 / (mp.mpf(1) / 4 + r * r)
    ser = mp.nsum(lambda k: r * r / ((mp.mpf(1) / 2 + 2 * k) * ((mp.mpf(1) / 2 + 2 * k) ** 2 + r * r)), [1, mp.inf])
    return closed, ser


def run(dps, X1, X2, n, lam, zeros_path):
    mp.mp.dps = dps
    Fa = 2 - mp.zeta(mp.mpf(1) / 2, derivative=1) / mp.zeta(mp.mpf(1) / 2)
    Fb = 2 - mp.log(mp.pi) / 2 + mp.digamma(mp.mpf(1) / 4) / 2
    gam = [mp.mpf(line.strip()) for line in open(zeros_path) if line.strip()]
    rows = []
    for ts in TS:
        t = mp.mpf(ts)
        W = mp.fsum(mp.exp(-t * g * g) for g in gam if t * g * g < 800)
        wtail = mp.exp(-t * gam[-1] ** 2) * (gam[-1] ** 2 + 1 / t)  # crude: N(T) <= T^2 for T >= 14
        lhs = 2 * mp.sqrt(mp.pi * t) * W
        Aa, Ab, Ac = Abar_all(ts, dps)
        res = {}
        for X in (X1, X2):
            S = prime_sum_double(n, lam, float(t), X)
            C = comp_int(t, mp.mpf(X))
            I = mp.mpf(S) - C
            res[X] = I
        I2 = res[X2]
        rhs = Aa + Fa - I2
        tb, ih = tail_pieces(t, mp.mpf(X2))
        rows.append(dict(t=ts, lhs=lhs, Aa=Aa, Ab=Ab, Ac=Ac, I=I2, I_X1=res[X1], rhs=rhs, resid=lhs - rhs, tail=tb + ih, wtail=wtail))
    return Fa, Fb, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--zeros", default=ZEROS)
    a = ap.parse_args()
    d1, d2, X1, X2 = (20, 26, 10 ** 5, 10 ** 6) if a.quick else (30, 40, 10 ** 6, 10 ** 7)
    out = a.out or ("/tmp/rx_centre.json" if a.quick else os.path.join(HERE, "..", "data", "rx_centre.json"))
    log = open(os.path.splitext(out)[0] + ".log", "w")

    def P(*v):
        line = " ".join(str(x) for x in v)
        print(line)
        log.write(line + "\n")
        log.flush()

    t0 = time.time()
    n, lam = prime_powers(X2)
    P("rx_centre.py: Remark bp:rem:centre; dps", d1, d2, "primes up to", X2, "(and", X1, ") quick" if a.quick else "full")
    # float64 prime sum against mpmath for n <= 2e5
    mp.mp.dps = 25
    m = n <= 200000
    nn, ll = n[m], lam[m]
    smp = mp.mpf(0)
    smp_terms = {}
    for t_s in ("0.5", "2"):
        t = mp.mpf(t_s)
        smp = mp.fsum(mp.log(int(round(math.exp(l)))) / mp.sqrt(int(k)) * mp.exp(-mp.log(int(k)) ** 2 / (4 * t)) for k, l in zip(nn, ll))
        sd = prime_sum_double(n, lam, float(t), 200000)
        P("prime sum n <= 2e5 at t =", t_s, ": mpmath 25 digits", mp.nstr(smp, 18), " double", repr(sd), " |diff|", mp.nstr(abs(smp - sd), 3))
        smp_terms[t_s] = mp.nstr(abs(smp - sd), 3)
    cl, se = check_sumq(d1)
    P("sum_k q_{1/2+2k}(1.3): closed form", mp.nstr(cl, 20), "series", mp.nstr(se, 20), "|diff|", mp.nstr(abs(cl - se), 3))
    Fa1, Fb1, r1 = run(d1, X1, X2, n, lam, a.zeros)
    Fa2, Fb2, r2 = run(d2, X1, X2, n, lam, a.zeros)
    P("F(1/2) = 2 - zeta'/zeta(1/2) =", mp.nstr(Fa2, 25), "; 2 - log(pi)/2 + psi(1/4)/2 =", mp.nstr(Fb2, 25), "; |diff|", mp.nstr(abs(Fa2 - Fb2), 3),
      "; |diff of the two precisions|", mp.nstr(abs(Fa1 - Fa2), 3))
    data = []
    for x1, x2 in zip(r1, r2):
        P("t =", x2["t"])
        P("  2 sqrt(pi t) W(t)  =", mp.nstr(x2["lhs"], 15), " (zeros beyond the 1700th: crude bound", mp.nstr(x2["wtail"], 3), ")")
        P("  Abar(t)  quad / nsum / nsum(r+s+e):", mp.nstr(x2["Aa"], 15), mp.nstr(x2["Ab"], 15), mp.nstr(x2["Ac"], 15), " max diff", mp.nstr(max(abs(x2["Aa"] - x2["Ab"]), abs(x2["Aa"] - x2["Ac"])), 3))
        P("  I(t) = int e^(-(log y)^2/4t) y^(-1/2) d(psi(y)-y) (X =", X2, ") =", mp.nstr(x2["I"], 15), " change from X =", X1, ":", mp.nstr(abs(x2["I"] - x2["I_X1"]), 3))
        P("  Abar + F - I =", mp.nstr(x2["rhs"], 15), " residual lhs - rhs =", mp.nstr(x2["resid"], 4), " tail error bound of I beyond X:", mp.nstr(x2["tail"], 3),
          " | precision diff of Abar", mp.nstr(abs(x1["Aa"] - x2["Aa"]), 3))
        data.append(dict(t=x2["t"], lhs=mp.nstr(x2["lhs"], 20), Abar=mp.nstr(x2["Aa"], 20), I=mp.nstr(x2["I"], 20), rhs=mp.nstr(x2["rhs"], 20), resid=mp.nstr(x2["resid"], 4),
                         tail=mp.nstr(x2["tail"], 4)))
    P("total seconds", round(time.time() - t0, 1))
    json.dump(dict(F_half=mp.nstr(Fa2, 25), F_half_other=mp.nstr(Fb2, 25), rows=data, X=X2, prime_sum_double_vs_mp=smp_terms), open(out, "w"), indent=1)
    log.close()


if __name__ == "__main__":
    main()
