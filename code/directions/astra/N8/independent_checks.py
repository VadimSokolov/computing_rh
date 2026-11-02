"""Independent numerical audit. No coauthor files or code are used.

Run from this folder with: timeout 600 python3 -I independent_checks.py
All results are floating point computations, not interval certificates.
"""
import bisect
import json
import math
import platform
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar

mp.mp.dps = 30
OUT = Path(__file__).resolve().parent
START = time.monotonic()
result = {"precision_dps": mp.mp.dps, "python": platform.python_version(),
          "mpmath": mp.__version__, "numpy": np.__version__,
          "status": "computed, not certified"}


def save():
    result["elapsed_seconds"] = time.monotonic() - START
    (OUT / "independent_numbers.json").write_text(json.dumps(result, indent=2) + "\n")


def say(label, value):
    result[label] = value
    save()
    print(label, json.dumps(value), flush=True)


def xi(s):
    return s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


a = 2 * mp.pi
say("constants", {"twice_clock_mean": str(mp.diff(xi, mp.mpf('0.5'), 2) / xi(mp.mpf('0.5'))),
                  "rho_sum": str(2 + mp.euler - mp.log(4 * mp.pi))})
say("bessel_ratios", {str(v): str(mp.besselk(0, a) / mp.besselk(mp.sqrt(v), a)) for v in [1, 10, 50]})

# The Davenport and Heilbronn function is independently evaluated by Hurwitz zeta.
kappa = mp.sqrt(10 - 2 * mp.sqrt(5)) - 2
kappa /= mp.sqrt(5) - 1
coefficients = [1, kappa, -kappa, -1]


def dh(s):
    return mp.power(5, -s) * sum(c * mp.zeta(s, mp.mpf(j) / 5)
                               for j, c in enumerate(coefficients, 1))


dhzero = mp.findroot(dh, (mp.mpc('.80', '85.69'), mp.mpc('.81', '85.70')))
pole = (dhzero - mp.mpf('.5')) ** 2
say("davenport_heilbronn", {"kappa": str(kappa), "zero": str(dhzero),
                           "residual": str(abs(dh(dhzero))), "reciprocal_pole": str(pole)})

# Compute zeta zeros afresh. This cache contains only this script's own output.
zero_cache = OUT / "independent_zeta_zeros.json"
if zero_cache.exists():
    gammas = [mp.mpf(x) for x in json.loads(zero_cache.read_text())]
else:
    gammas = []
for n in range(len(gammas) + 1, 401):
    gammas.append(mp.im(mp.zetazero(n)))
    if n % 25 == 0:
        zero_cache.write_text(json.dumps([str(x) for x in gammas]) + "\n")
        print("computed_zeta_zeros", n, flush=True)
zero_cache.write_text(json.dumps([str(x) for x in gammas]) + "\n")
g = np.array([float(x) for x in gammas])
say("zeta_endpoints", [str(gammas[0]), str(gammas[-1])])
say("screw_partial_mean", float(np.sum(2 / g ** 2)))


def screw400(t):
    return float(np.sum(2 * (1 - np.cos(g * t)) / g ** 2))


grid = np.linspace(.1, 12, 119001)
values = np.empty(len(grid))
for j in range(0, len(grid), 1000):
    block = grid[j:j + 1000]
    values[j:j + 1000] = np.sum(2 * (1 - np.cos(block[:, None] * g)) / g ** 2, axis=1)
j = int(np.argmin(values))
minimum = minimize_scalar(screw400, bounds=(grid[max(0, j-1)], grid[min(len(grid)-1, j+1)]), method="bounded")
say("screw_minimum", {"grid_step": float(grid[1]-grid[0]), "grid_minimum": float(values[j]),
                       "grid_argmin": float(grid[j]), "refined_minimum": float(minimum.fun),
                       "refined_argmin": float(minimum.x),
                       "weyl_tail_mean_estimate": float((mp.log(gammas[-1] / (2*mp.pi)) + 1) / (mp.pi*gammas[-1]))})

# Counting phase scan without approximating reference roots. In theta coordinates,
# N_ref(t) differs from floor(theta(t)/pi - c) by a constant, irrelevant to span.
phases = [mp.siegeltheta(x) / mp.pi for x in gammas]
lower, upper = mp.mpf(30), mp.mpf(680)
phase_lower = mp.siegeltheta(lower)/mp.pi
phase_upper = mp.siegeltheta(upper)/mp.pi


def discrepancy_span(c):
    ds = [bisect.bisect_right(g, float(lower)) - int(mp.floor(phase_lower-c)),
          bisect.bisect_right(g, float(upper)) - int(mp.floor(phase_upper-c))]
    for n, (height, phase) in enumerate(zip(gammas, phases), 1):
        if lower < height < upper:
            reference = int(mp.floor(phase-c))
            ds.extend([n-1-reference, n-reference])
    return max(ds)-min(ds)


histogram = {}
for j in range(200):
    span = discrepancy_span(mp.mpf(j)/200)
    histogram[str(span)] = histogram.get(str(span), 0) + 1
breaks = sorted(set([mp.mpf(0), mp.mpf(1)] +
                    [p-mp.floor(p) for p,h in zip(phases,gammas) if lower < h < upper] +
                    [phase_lower-mp.floor(phase_lower), phase_upper-mp.floor(phase_upper)]))
cell_histogram = {}
for left, right in zip(breaks[:-1], breaks[1:]):
    span = discrepancy_span((left+right)/2)
    cell_histogram[str(span)] = cell_histogram.get(str(span), 0) + 1
say("phase_scan", {"grid_histogram": histogram, "open_cell_histogram": cell_histogram,
                   "breakpoint_count": len(breaks),
                   "note": "Open phase cells exhaust generic phases. Coincident events require simultaneous jumps."})
selected = [(n,p) for n,(h,p) in enumerate(zip(gammas,phases),1) if lower < h < upper]
say("S_extrema_on_30_680", {"left_limit_min": str(min(n-2-p for n,p in selected)),
                           "right_limit_max": str(max(n-1-p for n,p in selected))})

# Find Bessel zeros by sign brackets, using exp(pi*nu/2) to keep scale reasonable.
def base(nu):
    return mp.re(mp.besselk(1j*nu, a)) * mp.exp(mp.pi*nu/2)


def two_term(nu):
    return (2*mp.pi**2*mp.re(mp.besselk(mp.mpf(9)/4+1j*nu,a))
            -3*mp.pi*mp.re(mp.besselk(mp.mpf(5)/4+1j*nu,a))) * mp.exp(mp.pi*nu/2)


def roots_by_brackets(fun, count, start=mp.mpf(5)):
    roots = []
    left, fleft = start, fun(start)
    step = mp.mpf('.25')
    while len(roots) < count:
        right = left + step
        fright = fun(right)
        if fleft * fright < 0:
            root = mp.findroot(fun, (left,right), solver="anderson", tol=mp.mpf('1e-24'))
            roots.append(root)
        left, fleft = right, fright
    return roots


base_roots = roots_by_brackets(base, 490)
say("base_roots", {"first_three": [str(x) for x in base_roots[:3]],
                   "last": str(base_roots[-1]),
                   "phase_offsets": [str((r*mp.log(2*r/a)-r+mp.pi/4)/mp.pi-k)
                                     for k,r in [(1,base_roots[0]),(490,base_roots[-1])]]})

# Specify a reproducible leading Weyl tail correction, since the note does not.
cut = base_roots[-1]
tail = (mp.log(2*cut/a)+1)/(mp.pi*cut)
products = {}
for lam in [1,10,50]:
    log_product = -sum(mp.log1p(lam/(r*r)) for r in base_roots)
    exact = mp.besselk(0,a)/mp.besselk(mp.sqrt(lam),a)
    corrected = mp.exp(log_product-lam*tail)
    products[str(lam)] = {"raw": str(mp.exp(log_product)), "leading_weyl_corrected": str(corrected),
                          "exact": str(exact), "corrected_relative_error": str(corrected/exact-1)}
say("bessel_products", products)
two_roots = roots_by_brackets(two_term, 400)
gram_roots = []
for n in range(410):
    r = mp.grampoint(n)
    if r > gammas[-1]:
        break
    gram_roots.append(r)


def interval_summary(roots):
    roots = [float(x) for x in roots if x < gammas[-1]]
    counts = [bisect.bisect_left(g, right)-bisect.bisect_right(g,left)
              for left,right in zip(roots[:-1],roots[1:])]
    failures = [i for i,n in enumerate(counts) if n != 1]
    return {"intervals":len(counts), "failures":len(failures),
            "first_failure_left":roots[failures[0]] if failures else None,
            "first_failure_right":roots[failures[0]+1] if failures else None,
            "first_failure_count":counts[failures[0]] if failures else None}


say("reference_intervals", {"gram":interval_summary(gram_roots),
                            "two_term":interval_summary([2*x for x in two_roots]),
                            "base":interval_summary([2*x for x in base_roots])})
(OUT / "independent_bessel_zeros.json").write_text(json.dumps({"base": [str(x) for x in base_roots],
                                                            "two_term": [str(x) for x in two_roots]},indent=2)+"\n")
say("finished", True)
