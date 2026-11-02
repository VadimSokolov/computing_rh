#!/usr/bin/env python3
"""D1: size of the direct certificate at the next rungs (plain double arithmetic, a planning table, not a result).

T from m(T) >= tau + 0.75, with tau = sup V, the Schur bound for the prime operator.
Nmax = 1.25 e L T / 2, modes = Nmax / 2, nodes = about 10 T.
Margin: Zhu's fitted law exp(-20.13 N / log N), N = number of zeros below 2 pi X: a conjecture, used only for the precision.
Bits: 1.27 Nmax (growth of the ball recurrence) + log2(1/margin) + 250.
Ops: nodes * modes^2 ball multiplications for the two sectors together, times 2.

Usage: timeout 30 python3 -I D1_scaling.py ../../../../data/zeros_8000.txt
"""
import math, sys
gam = sorted(float(x) for x in open(sys.argv[1]).read().split() if float(x) > 1)


def ppower(n):
    for p in range(2, n + 1):
        if n % p == 0:
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            return p if n == 1 else None
    return None


def comb(X):
    return [(math.log(n), 2 * math.log(ppower(n)) / math.sqrt(n)) for n in range(2, X) if ppower(n)]


def tau(L, cb):
    pts = sorted(set([0.0, L] + [abs(L - c) for c, _ in cb]))
    best = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-11:
            continue
        x = 0.5 * (a + b)
        v = 0.5 * sum(w for c, w in cb if c <= L + x) + 0.5 * sum(w for c, w in cb if c <= L - x)
        best = max(best, v)
    return best


print("X  L  A_L  tau  T  T_1=2pi e^A  horizon 2pi X  Nmax  modes  N(2piX)  log10 margin(law)  bits  ops")
for X in [3, 5, 7, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 49, 53, 64]:
    L = 0.5 * math.log(X)
    cb = comb(X)
    A = sum(w for _, w in cb)
    t = tau(L, cb)
    T = 2 * math.pi * math.exp(t + 0.75)
    Nm = 1.25 * math.e * L * T / 2
    modes = Nm / 2
    N = sum(1 for g in gam if g < 2 * math.pi * X)
    lm = -20.13 * N / math.log(N) / math.log(10) if N > 1 else 0
    bits = 1.27 * Nm - lm * math.log2(10) + 250
    ops = 2 * 10 * T * modes ** 2
    print("%d %.4f %.3f %.3f %.3g %.3g %.1f %.3g %.3g %d %.1f %.3g %.3g" % (X, L, A, t, T, 2 * math.pi * math.exp(A), 2 * math.pi * X, Nm, modes, N, lm, bits, ops))
