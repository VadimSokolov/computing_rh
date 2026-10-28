# gl_pole_lock.py: the pole lock heuristic of Section "Landau Siegel zeros and the twisted heat trace" of Chapter 21 (label prefix gl:,
# the end of Section 21.6, Proposition gl:prop:polelock). If chi(p) = -1 for every prime p, the prime sum of the Gaussian explicit formula
#   S(x) = (2 sqrt(pi x))^-1 sum_p log p p^(-1/2) exp(-(log p)^2/(4x))
# is compared with exp(x/4) Phi(sqrt(x/2)), Phi the standard normal distribution function, at x = 0.5, 1, 2, 3 and primes p <= N = 6e6.
# Output: the ratio S/exp(x/4), Phi(sqrt(x/2)), the part of the integral (2 sqrt(pi x))^-1 int_0^(log 2) e^(u/2 - u^2/(4x)) du / exp(x/4)
# that lies below the first prime (the integral over u < log 2), and the remainder, the error of the prime number theorem at small primes.
# Usage: python3 gl_pole_lock.py [--N 6000000] [--out ../data/gl_pole_lock.json]. About 1 s on one core. Heuristic, not a theorem.
import argparse
import json
from math import erf, exp, log, pi, sqrt

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--N", type=int, default=6_000_000)
ap.add_argument("--out", default="../data/gl_pole_lock.json")
a = ap.parse_args()
N = a.N
s = np.ones(N + 1, bool)
s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
p = np.nonzero(s)[0].astype(float)
lp = np.log(p)
Phi = lambda z: 0.5 * (1 + erf(z / sqrt(2)))
rows = []
print("x   prime sum   exp(x/4)   ratio   Phi(sqrt(x/2))   below first prime   remainder")
for x in [0.5, 1, 2, 3]:
    S = float(np.sum(lp * p ** -0.5 * np.exp(-lp ** 2 / (4 * x))) / (2 * np.sqrt(np.pi * x)))
    E = exp(x / 4)
    u = np.linspace(0.0, log(2.0), 200001)
    f = np.exp(u / 2 - u ** 2 / (4 * x))
    low = float(np.sum((f[1:] + f[:-1]) / 2 * np.diff(u)) / (2 * sqrt(pi * x)) / E)
    ph = Phi(sqrt(x / 2))
    rows.append({"x": x, "prime_sum": S, "exp_x_over_4": E, "ratio": S / E, "Phi_sqrt_x_over_2": ph,
                 "integral_below_log2_over_exp": low, "remainder": ph - low - S / E})
    print(f"{x:<4} {S:9.4f} {E:9.4f} {S/E:7.3f} {ph:8.3f} {low:8.3f} {ph-low-S/E:8.3f}")
with open(a.out, "w") as fh:
    json.dump({"N": N, "rows": rows}, fh, indent=1)
