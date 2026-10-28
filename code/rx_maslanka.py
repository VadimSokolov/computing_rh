# rx_maslanka.py: the numerical claims of Proposition ch19:maslanka (Maslanka's series as Laguerre moments of Williams' variable; Section 3.6 of
# misc/reviews/nick-b4/reconcile-N.md, task E4; source: the reciprocal Xi note, cited N:n; chapter: ch/ch19.tex after line 86). Added for the book in October 2026,
# written from the mathematics (the collaborator's scripts were not read, run or copied).
# Claim. With A_k = sum_{j=0}^k (-1)^j C(k,j) (2j+1) zeta(2j+2), W = sum_n Gamma_n/n^2 (Gamma_n of shape 2) and E e^(-lambda W) = phi(lambda) = (pi sqrt(lambda)/sinh(pi sqrt(lambda)))^2:
# (i) A_k = E[W L_k^(1)(W)]/(2(k+1)); (ii) sum_k A_k y^k = (1 - phi(y/(1-y)))/(2y), |y| < 1; (iii) for 0 < Re(sigma) < 2,
# (sigma-1) zeta(sigma) = sin(pi sigma/2)/(2 pi) int_0^inf r^(-sigma/2-1) (1 - phi(r)) dr, equal for real sigma to E[(1 - E[e^(-W R)|B])/(2(1-B))] with B ~ Beta(sigma/2, 1 - sigma/2), R = (1-B)/B;
# and Maslanka's expansion (sigma-1) zeta(sigma) = sum_k A_k (1 - sigma/2)_k/k!.
# Method. A_k from the definition. The power series of (z/sinh z)^2 = sum_m a_m z^(2m) from z/sinh z = sum_n (2 - 2^(2n)) B_(2n) z^(2n)/(2n)! (Bernoulli numbers, no zeta values),
# so phi(lambda) = sum_m a_m pi^(2m) lambda^m, E W^m = (-1)^m m! a_m pi^(2m); then (ii) gives A_k = -(1/2) sum_{m=1}^{k+1} a_m pi^(2m) C(k, k+1-m) (the coefficient of y^(k+1) in 1 - phi(y/(1-y)) is
# -sum_m a_m pi^(2m) C(k, k+1-m)), and (i) gives A_k = sum_j (-1)^j C(k+1, k-j) E W^(j+1)/(j! 2 (k+1)); so (i) and (ii) are compared with the definition, which uses zeta values: three computations that share only the
# identity under test. (iii) by mp.quad on [0, 1, 10, 100, 1000, 1e4, inf) (1 - phi evaluated at +60 digits near 0) at sigma = 0.5, 1.3, 0.5 + 3i against (sigma-1) zeta(sigma); the beta form at sigma = 0.5, 1.3 by quadrature in b;
# the partial sums of Maslanka's series up to k = 15 against (sigma-1) zeta(sigma).
# Error control. Everything at two working precisions (70 and 110 digits), the largest difference logged (the printed digits are those of the second run).
# Output: data/rx_maslanka.json (numbers as strings) and the log rx_maslanka.log next to it.
# Usage: python3 rx_maslanka.py [--quick] [--out FILE]; --quick uses 40 and 60 digits and writes to /tmp unless --out is given; the full case takes about 10 seconds; on the cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g rx_maslanka.json -g rx_maslanka.log code/rx_maslanka.py -- --out rx_maslanka.json
import argparse
import json
import os
import time
from multiprocessing import Pool

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
K = 16  # k = 0..15
SIGMAS = ["0.5", "1.3", "0.5+3j"]


def a_coeffs(M):
    """a_m, m = 0..M, of (z/sinh z)^2 = sum a_m z^(2m)."""
    b = [(2 - mp.mpf(2) ** (2 * n)) * mp.bernoulli(2 * n) / mp.factorial(2 * n) for n in range(M + 1)]
    return [mp.fsum(b[i] * b[m - i] for i in range(m + 1)) for m in range(M + 1)]


def phi(r):
    x = mp.pi * mp.sqrt(r)
    return (x / mp.sinh(x)) ** 2


def one_minus_phi(r):
    with mp.workdps(mp.mp.dps + 60):
        v = 1 - phi(mp.mpf(r))
    return +v


def job(dps):
    mp.mp.dps = dps
    out = {}
    A_def = [mp.fsum((-1) ** j * mp.binomial(k, j) * (2 * j + 1) * mp.zeta(2 * j + 2) for j in range(k + 1)) for k in range(K)]
    a = a_coeffs(K + 2)
    A_ii = [-mp.fsum(a[m] * mp.pi ** (2 * m) * mp.binomial(k, k + 1 - m) for m in range(1, k + 2)) / 2 for k in range(K)]
    EW = [(-1) ** m * mp.factorial(m) * a[m] * mp.pi ** (2 * m) for m in range(K + 2)]
    A_i = [mp.fsum((-1) ** j * mp.binomial(k + 1, k - j) * EW[j + 1] / mp.factorial(j) for j in range(k + 1)) / (2 * (k + 1)) for k in range(K)]
    out["A_def"], out["A_ii"], out["A_i"] = A_def, A_ii, A_i
    iii = {}
    for sg in SIGMAS:
        s = mp.mpmathify(sg)
        z = s / 2
        f = lambda r: r ** (-z - 1) * one_minus_phi(r)
        q = mp.sin(mp.pi * z) / (2 * mp.pi) * mp.quad(f, [0, 1, 10, 100, 1000, 10000, mp.inf])
        target = (s - 1) * mp.zeta(s)
        ser = [mp.fsum(A_def[k] * mp.rf(1 - z, k) / mp.factorial(k) for k in range(n + 1)) for n in (5, 10, 15)]
        rec = dict(mellin=q, target=target, series=ser)
        if mp.im(s) == 0:
            zz = z
            g = lambda b: b ** (zz - 1) * (1 - b) ** (-zz) / (mp.pi / mp.sin(mp.pi * zz)) * one_minus_phi((1 - b) / b) / (2 * (1 - b))
            rec["beta"] = mp.quad(g, [0, 0.5, 0.9, 0.99, 1])
        iii[sg] = rec
    out["iii"] = iii
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    d1, d2 = (40, 60) if a.quick else (70, 110)
    out = a.out or ("/tmp/rx_maslanka.json" if a.quick else os.path.join(HERE, "..", "data", "rx_maslanka.json"))
    log = open(os.path.splitext(out)[0] + ".log", "w")

    def P(*v):
        line = " ".join(str(x) for x in v)
        print(line)
        log.write(line + "\n")
        log.flush()

    t0 = time.time()
    with Pool(2) as pool:
        r1, r2 = pool.map(job, [d1, d2])
    mp.mp.dps = d2
    P("rx_maslanka.py: Proposition ch19:maslanka; dps", d1, d2, "quick" if a.quick else "full")
    res = dict(A=[])
    for k in range(K):
        d = r2["A_def"][k]
        P("A_%d = %s | (ii) minus def %s | (i) minus def %s | precision diff %s" % (k, mp.nstr(d, 22), mp.nstr(abs(r2["A_ii"][k] - d), 3), mp.nstr(abs(r2["A_i"][k] - d), 3), mp.nstr(abs(d - r1["A_def"][k]), 3)))
        res["A"].append(mp.nstr(d, 50))
    P("A_0 - pi^2/6 =", mp.nstr(r2["A_def"][0] - mp.pi ** 2 / 6, 3))
    res["iii"] = {}
    for sg in SIGMAS:
        r = r2["iii"][sg]
        t = r["target"]
        P("sigma = %s: (sigma-1) zeta(sigma) = %s ; Mellin form %s (|diff| %s, precision diff %s)" % (sg, mp.nstr(t, 18), mp.nstr(r["mellin"], 18), mp.nstr(abs(r["mellin"] - t), 3), mp.nstr(abs(r["mellin"] - r1["iii"][sg]["mellin"]), 3)))
        if "beta" in r:
            P("    beta form %s (|diff| %s)" % (mp.nstr(r["beta"], 18), mp.nstr(abs(r["beta"] - t), 3)))
        P("    Maslanka partial sums k <= 5, 10, 15: |sum - target| = %s, %s, %s" % tuple(mp.nstr(abs(x - t), 3) for x in r["series"]))
        res["iii"][sg] = dict(target=mp.nstr(t, 20), mellin_diff=mp.nstr(abs(r["mellin"] - t), 3), series_err=[mp.nstr(abs(x - t), 3) for x in r["series"]])
    P("total seconds", round(time.time() - t0, 1))
    json.dump(res, open(out, "w"), indent=1)
    log.close()


if __name__ == "__main__":
    main()
