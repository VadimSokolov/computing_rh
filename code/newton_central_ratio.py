"""Reproduces the table of the Lambert result for the central ratio of the Laguerre expressions of the heat flow (Chapter 12 or 21 of the book).

With M_k(t) = int u^(2k) e^(t u^2) Phi(u) du (Phi the Polya kernel in the variable of the Newton chapter, Phi(u) = sum_m pi m^2 (2 pi m^2 e^(4u) - 3) e^(5u - pi m^2 e^(4u)), which is the kernel of eq:ch2:Phi at 2u up to a constant, so R_n here is one quarter of the value in the variable of eq:ch2:Phi), the central ratio is R_n(t) = A_(n+1) / (2 A_n), where A_n = sum_j binom(2n, 2j) M_j M_(n-j), which equals k_n U_(n+1)(0) / U_n(0) with k_n = (n+1)(2n+1) and is half the second moment of b = u - v under the law proportional to b^(2n) nu_t(u) nu_t(v). The moments are computed in logs on a midpoint grid of step 2e-5 on u > 0 (windows down to e^-32 of the maximum) and compared with (1/8) W(n/pi)^2. Output: code/data/newton_central_ratio.json.
"""
import json
import os
import time

import numpy as np
from scipy.special import gammaln, lambertw, logsumexp

H = 2e-5
UMAX = 6.0
CUT = 32.0
COARSE = 500


def log_phi(u):
    """log Phi(u) for u > 0, factoring out exp(-pi e^(2u)) to avoid underflow."""
    e2 = np.exp(4.0 * u)
    s = np.zeros_like(u)
    for m in range(1, 13):
        s += (np.pi * m * m * (2 * np.pi * m * m * e2 - 3) * np.exp(5 * u)
              * np.exp(-np.pi * (m * m - 1) * e2))
    return -np.pi * e2 + np.log(s)


U = (np.arange(int(UMAX / H)) + 0.5) * H
LOGU = np.log(U)
U2 = U * U
LOGPHI = log_phi(U)


def log_moment(k, t):
    """log M_k(t), M_k = 2 * int_0^inf u^(2k) e^(t u^2) Phi(u) du, by the midpoint rule."""
    f_c = 2 * k * LOGU[::COARSE] + t * U2[::COARSE] + LOGPHI[::COARSE]
    i0 = int(np.argmax(f_c)) * COARSE
    m = 2 * COARSE
    while True:
        lo, hi = max(0, i0 - m), min(len(U), i0 + m)
        f = 2 * k * LOGU[lo:hi] + t * U2[lo:hi] + LOGPHI[lo:hi]
        mx = f.max()
        ok_lo = lo == 0 or f[0] < mx - CUT - 2
        ok_hi = hi == len(U) or f[-1] < mx - CUT - 2
        if ok_lo and ok_hi:
            break
        m *= 2
    return np.log(2 * H) + mx + np.log(np.exp(f - mx).sum())


def log_A(n, t, lm, klo):
    j = np.arange(0, n + 1)
    ok = (j >= klo) & (j < klo + len(lm)) & (n - j >= klo) & (n - j < klo + len(lm))
    j = j[ok]
    terms = (gammaln(2 * n + 1) - gammaln(2 * j + 1) - gammaln(2 * (n - j) + 1)
             + lm[j - klo] + lm[n - j - klo])
    edge = max(terms[0], terms[-1]) - terms.max()
    return logsumexp(terms), edge


def central_ratio(n, t):
    w = int(12 * np.sqrt(n)) + 30
    klo = max(0, (n + 1) // 2 - w)
    khi = min(n + 1, (n + 1) // 2 + w + 1)
    lm = np.array([log_moment(k, t) for k in range(klo, khi + 1)])
    a0, e0 = log_A(n, t, lm, klo)
    a1, e1 = log_A(n + 1, t, lm, klo)
    full = klo == 0 and khi == n + 1
    return 0.5 * np.exp(a1 - a0), max(e0, e1), full


def main():
    rows = []
    t0 = time.time()
    print(f"{'n':>8} {'t':>4} {'R_n':>10} {'W form':>10} {'diff':>10} {'edge':>8}")
    for t in (0.5, 0.0):
        for n in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
            r, edge, full = central_ratio(n, t)
            lam = float((lambertw(n / np.pi).real ** 2) / 8)
            rows.append({"n": n, "t": t, "R_n": float(f"{r:.6g}"),
                         "lambert_form": float(f"{lam:.6g}"),
                         "difference": float(f"{r - lam:.6g}")})
            print(f"{n:>8} {t:>4} {r:10.5f} {lam:10.5f} {r - lam:10.5f} {edge:8.1f}"
                  f"{'' if full else ' (window)'}", flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(here, "data"), exist_ok=True)
    with open(os.path.join(here, "data", "newton_central_ratio.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    print(f"time {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
