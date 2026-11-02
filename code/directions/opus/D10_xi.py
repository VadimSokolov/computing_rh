#!/usr/bin/env python3
"""D10: Xi with planted pairs: the exact criterion of Theorem 5, the two sided bracket, and the book's fixed point.

For a planted pair c = m + i eta with G = F/((z-c)(z-cbar)), at height w < eta put a = sqrt(eta^2-w^2),
alpha = arcsin(w/eta), T(u) = arg G(m+u+iw) (unwrapped on |u| <= a), A(u) = arg((u+iw)^2+eta^2), Phi = T + A.
 exact: the zeros born from the pair are real iff Phi rises through a level pi/2 + k pi on (-a,a).
 gap upper bound (Theorem 6(a)): w <= w* whenever max_{|u|<=a} dist(T, pi Z) <= arccos(w/eta).
 gap lower bound (Theorem 6(b)): w >= w* whenever a level lies strictly between Phi(-a) and Phi(a) (for all larger w).
Floating point (mpmath, 20 digits), nothing certified.
Usage: timeout 590 python3 -I D10_xi.py
"""
import math
import mpmath as mp
import numpy as np

mp.mp.dps = 20
HALF = mp.mpf(1) / 2


def Xi(z):
    s = HALF + 1j * z
    return HALF * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


class Config:
    def __init__(self, k, pairs):
        self.g1 = mp.zetazero(k).imag
        self.g2 = mp.zetazero(k + 1).imag
        self.pairs = [(mp.mpf(m), mp.mpf(h)) for m, h in pairs]

    def G(self, z, j):
        """F with pair j removed."""
        v = Xi(z) / ((z * z - self.g1 ** 2) * (z * z - self.g2 ** 2))
        for i, (m, h) in enumerate(self.pairs):
            c = m + 1j * h
            cb = m - 1j * h
            if i == j:
                v *= (z + c) * (z + cb)
            else:
                v *= (z * z - c * c) * (z * z - cb * cb)
        return v

    def field(self, w, j):
        m, h = self.pairs[j]
        z = m + 1j * w
        d = mp.mpf(10) ** (-7)
        g = self.G(z, j)
        return (self.G(z + d, j) - self.G(z - d, j)) / (2 * d * g), g

    def profile(self, w, j, n=41):
        m, h = self.pairs[j]
        a = mp.sqrt(h * h - w * w)
        us = [(-a + 2 * a * i / (n - 1)) for i in range(n)]
        T = np.unwrap([float(mp.arg(self.G(m + u + 1j * w, j))) for u in us])
        A = np.array([float(mp.arg((u + 1j * w) ** 2 + h * h)) for u in us])
        return T, A


def levels_up(Phi):
    lev = np.floor((Phi - np.pi / 2) / np.pi)
    return np.any(lev > np.minimum.accumulate(lev))


def bisect(pred, lo, hi, steps=26):
    """pred(hi) True, pred(lo) False: smallest w with pred True on the connected piece below hi."""
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        if pred(m):
            hi = m
        else:
            lo = m
    return hi


def analyse(cfg, j, label, printed=None):
    m, h = cfg.pairs[j]
    h_f = float(h)

    def crit(w):
        T, A = cfg.profile(mp.mpf(w), j)
        return levels_up(T + A)

    def upper_gap_pred(w):  # True while NOT yet certainly nonreal: tau(w) > arccos(w/eta)
        T, A = cfg.profile(mp.mpf(w), j)
        d = np.abs((T + np.pi / 2) % np.pi - np.pi / 2)
        return d.max() > math.acos(w / h_f)

    def lower_gap_pred(w):
        T, A = cfg.profile(mp.mpf(w), j)
        return np.floor((T[-1] + A[-1] - np.pi / 2) / np.pi) > np.floor((T[0] + A[0] - np.pi / 2) / np.pi)

    lo = 0.55 * h_f
    hi = h_f * (1 - 1e-7)
    w_star = bisect(crit, lo, hi)
    w_f = bisect(upper_gap_pred, lo, hi)
    w_s = bisect(lower_gap_pred, lo, hi)
    # book's fixed point, iterated from eps = 0
    eps = 0.0
    for _ in range(40):
        w = h_f - eps
        Y, g = cfg.field(mp.mpf(w), j)
        arg = float(mp.arg(g))
        dist = abs((arg + math.pi / 2) % math.pi - math.pi / 2)
        b = -float(Y.imag)
        new = dist ** 2 / (4 * (b + 1 / (2 * h_f - eps)))
        if abs(new - eps) < 1e-12:
            eps = new
            break
        eps = new
    gap = h_f - w_star
    line = (f"{label:26s} gap {gap:.4e}  lower {h_f-w_s:.4e}  upper {h_f-w_f:.4e}  fp {eps:.4e}  "
            f"fp/gap {eps/gap:.4f}  book {printed if printed else ''}")
    print(line, flush=True)


def main():
    print("one planted pair at the midpoint of gamma_k, gamma_{k+1}")
    book = {(1, 0.1): 0.099810, (1, 0.3): 0.295090, (3, 0.1): 0.099746, (3, 0.3): 0.293544,
            (8, 0.1): 0.099711, (8, 0.3): 0.292783, (15, 0.1): 0.099619, (15, 0.3): 0.290699}
    for (k, h), val in book.items():
        g1 = mp.zetazero(k).imag
        g2 = mp.zetazero(k + 1).imag
        cfg = Config(k, [((g1 + g2) / 2, h)])
        analyse(cfg, 0, f"k={k} eta={h}", f"gap {h-val:.4e}")
    print("two planted pairs at M -+ d/2 replacing gamma_8, gamma_9 (Table cm:tab:pairs)")
    g1 = mp.zetazero(8).imag
    g2 = mp.zetazero(9).imag
    M = (g1 + g2) / 2
    table = {(3.0, 0.1): (5.88e-4, 8.92e-5), (2.0, 0.1): (1.13e-3, 6.47e-7), (0.8, 0.1): (4.18e-3, 1.21e-3),
             (0.3, 0.1): (1.31e-2, 8.62e-3), (1.2, 0.2): (1.53e-2, 1.93e-3)}
    for (d, h), (gl, gu) in table.items():
        cfg = Config(8, [(M - mp.mpf(d) / 2, h), (M + mp.mpf(d) / 2, h)])
        analyse(cfg, 0, f"d={d} eta={h} lower", f"gap {gl}")
        analyse(cfg, 1, f"d={d} eta={h} upper", f"gap {gu}")


if __name__ == "__main__":
    main()
