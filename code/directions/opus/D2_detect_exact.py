#!/usr/bin/env python3
"""D2: the heuristic exact sign criterion for a planted quadruple, 2 U_N(gamma) >= 2 log(1/delta) + (9/2) log b + c0,
with the potential U_N of Theorem 4 at the book's number of modes, N = ceil(a (gamma+45)/pi) + 2, and at N >= N_* (saturated).
For each entry a of Table zh:tab:det, criterion (ii), prints the margin m(a) = 2U - 2 log(1/delta) - (9/2) log b at a - 0.05 and at a:
the constant c0 must lie between them if the criterion is right. Then predicts a_det at the upper two heights.
Run: timeout 60 python3 -I D2_detect_exact.py
Status: heuristic, floating point.
"""
import math
from scipy import integrate


def params(Nx, L):
    T = math.pi * Nx / L
    if T >= 4 * math.pi * math.exp(2 * L):
        return math.pi / 2, 4 * math.pi * math.exp(2 * L)
    lo, hi = 1e-12, math.pi / 2
    g = lambda p: (4 * L / math.pi) * p - math.log(T / (4 * math.pi * math.sin(p)))
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if g(m) > 0 else (m, hi)
    p = 0.5 * (lo + hi)
    return p, T / math.sin(p)


def U(t, Nx, L):
    p, b = params(Nx, L)
    if t >= b:
        return 0.0, b

    def f(s):
        th = math.asin(min(s / b, 1.0))
        val = 0.5 * math.acos(min(s / b, 1.0))
        if p < math.pi / 2 - 1e-12:
            val += (L / math.pi) * math.log(abs(math.sin(p + th) / math.sin(p - th)))
        return val
    pts = [b * math.sin(p)] if t < b * math.sin(p) < b else None
    return integrate.quad(f, t, b, points=pts, limit=400)[0], b


def margin(g, d, a, mode):
    L = a
    if mode == 'book':
        N = math.ceil(a * (g + 45) / math.pi) + 2
    else:
        N = 4 * a * math.exp(2 * a) + 1
    u, b = U(g, N, L)
    return 2 * u - 2 * math.log(1 / d) - 4.5 * math.log(b), N, b


def main():
    table = {51.37: {0.01: 1.05, 0.05: 1.00, 0.2: 0.95, 0.45: 0.95}, 199.64: {0.01: 1.55, 0.05: 1.50, 0.2: 1.50, 0.45: 1.50}}
    print('gamma delta entry | book modes: N, b, m(a-.05), m(a) | saturated: b, m(a-.05), m(a)')
    for g, row in table.items():
        for d, a in row.items():
            m0, _, _ = margin(g, d, a - 0.05, 'book')
            m1, N, b = margin(g, d, a, 'book')
            s0, _, _ = margin(g, d, a - 0.05, 'sat')
            s1, _, bs = margin(g, d, a, 'sat')
            print(f"{g:8.2f} {d:5.2f} {a:5.2f} | {N:5d} {b:8.1f} {m0:8.1f} {m1:8.1f} | {bs:8.1f} {s0:8.1f} {s1:8.1f}")
    print('prediction of a_det on a fine grid, delta = 0.01 and 0.05, for c0 = 0, 6, 12')
    for g in (199.64, 1000.57, 1999.99):
        for d in (0.01, 0.05):
            out = []
            for mode in ('book', 'sat'):
                for c0 in (0, 6, 12):
                    a = 0.5 * math.log(g / (4 * math.pi))
                    while margin(g, d, a, mode)[0] < c0:
                        a += 0.002
                    out.append(round(a, 3))
            print(g, d, 'book modes', out[:3], 'saturated', out[3:])


if __name__ == '__main__':
    main()
