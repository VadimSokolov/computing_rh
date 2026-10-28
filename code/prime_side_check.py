# prime_side_check.py: the smoothed prime sum Pi(T;x) = sum_n Lambda(n) n^(-1/2) e^(-(log n)^2/8x) cos(T log n) at heights from 10^12 to 10^50, the weight bound
# M(x) = 2 sum_n Lambda(n) n^(-1/2) e^(-(log n)^2/8x) and the scale x_c(T) with M(x_c) = log(T/2 pi), for the unnumbered table of the subsection
# "The prime side of the ladder" (ch/ch12.tex). By (5.1.3) of Section 5.1, 2 Pi(T;x) < log(T/2 pi) is, to leading order, the positivity of the rung at the height T and
# the scale x. The phases T log n are computed in ball arithmetic at 300 bits (T log n must be known modulo 2 pi); the sums over prime powers n <= N use float weights.
# The omitted weight is the Prime Number Theorem estimate 2 int_{log N}^oo e^(L/2 - L^2/8x) dL of the sum of the absolute values beyond N (an estimate, not a bound).
# Run: python3 prime_side_check.py; writes prime_side_check.json.
import json
import math
import time

import mpmath as mp
import numpy as np
from flint import arb, ctx

ctx.prec = 300
T00 = time.time()


def prime_powers(N):  # (n, log p) for the prime powers n <= N
    sieve = np.ones(N + 1, bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    ps = np.nonzero(sieve)[0]
    for p in ps.tolist():
        q = p
        while q <= N:
            yield q, p
            q *= p


def tail(N, x):  # 2 int_{log N}^oo e^(L/2 - L^2/8x) dL
    f = lambda L: mp.e ** (L / 2 - L * L / (8 * x))
    return float(2 * mp.quad(f, [math.log(N), mp.inf]))


def Pi(T, x, N):
    Tb, tot = arb(T), 0.0
    logs = {}
    for n, p in prime_powers(N):
        if p not in logs:
            logs[p] = math.log(p)
        L = arb(n).log()
        Lf = float(L.mid())
        tot += logs[p] / math.sqrt(n) * math.exp(-Lf * Lf / (8 * x)) * float((Tb * L).cos().mid())
    return tot


# M(x) from the prime powers n <= 10^8 with the PNT estimate beyond (the weights are positive; the estimate is below 10^-3 for x <= 2.4)
NM = 10**8
PP = np.array([(n, p) for n, p in prime_powers(NM)], dtype=np.int64)
lg, lp = np.log(PP[:, 0].astype(float)), np.log(PP[:, 1].astype(float))
print('prime powers to 1e8', len(PP), round(time.time() - T00, 1), flush=True)


def M(x):
    return float(2 * np.sum(lp / np.sqrt(PP[:, 0].astype(float)) * np.exp(-lg * lg / (8 * x)))) + tail(NM, x)


def bisect(f, lo, hi):
    flo = f(lo)
    for _ in range(50):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


out = {'M': {str(x): M(x) for x in (0.5, 1.0, 1.5, 1.7, 2.0)}, 'x_c': {}, 'rows': []}
for T in (10**12, 3 * 10**12, 10**20):
    out['x_c'][str(T)] = bisect(lambda x: M(x) - math.log(T / (2 * math.pi)), 1.0, 2.4)
print(out, round(time.time() - T00, 1), flush=True)
rows = [(10**12, 0.5, 22), (10**12, 1.0, 22), (10**12, 1.7, 20), (10**12, 2.5, 16), (10**12, 3.0, 14), (3 * 10**12, 1.7, 20),
        (3 * 10**12, 3.0, 14), (10**20, 1.0, 22), (10**30, 1.0, 22), (10**50, 1.0, 22)]
for T, x, D in rows:
    N = int(math.exp(math.sqrt(8 * x * D)))
    P = Pi(T, x, N)
    r = dict(T=str(T), x=x, D=D, N=N, twoPi=2 * P, logT=math.log(T / (2 * math.pi)), detects=math.pi / (4 * x), omitted=tail(N, x))
    out['rows'].append(r)
    print(r, round(time.time() - T00, 1), flush=True)
json.dump(out, open('prime_side_check.json', 'w'), indent=1)
