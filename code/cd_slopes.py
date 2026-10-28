# cd_slopes.py: the numerical claims of the chapter on Cauchy deconvolution (ch/cd1.tex): the example of the signed source
# (Example cd:ex:signed and the critical width after it) and the stability experiment of Table cd:tab:slopes (Proposition cd:prop:finite).
# Written for the book in October 2026 from the mathematics in ch/cd1.tex.
# Part 1 (signed source). u = P_a - q P_a(. - d) with a = 1/2, d = 1. Positivity of P_(a-c)*mu holds iff q <= (r_c - d)/(r_c + d),
# r_c = sqrt(d^2 + 4 (a-c)^2), i.e. iff c <= c* = a - sqrt(m^2 - d^2)/2, m = d (1+q)/(1-q). Checked: q = 9 - 4 sqrt5 gives m = sqrt5/2 and
# c* = 1/4 (mpmath, 30 digits); the critical width is also found by bisection on the minimum of the density (grid plus golden section) and
# compared with the formula; q = 0.0557 gives c* = 0.25007; the datum bound (r - d)/(r + d) at c = 0 is 3 - 2 sqrt2.
# Part 2 (stability). Source: 40 Lorentzians of width w = 0.05 on a period L = 400, centres uniform on [0, L) and weights uniform on [0.5, 1.5]
# (seed 12345; the source is a test source, the slope does not depend on the choice), mu = sum_j m_j delta(x - x_j). The Fourier coefficients
# (frequencies omega_k = 2 pi k / L, k = -K..K) of P_s * mu are exp(-s |omega|) mu_k. Data at width a = 1/2: exp(-a|omega|) mu_k + eps e_k
# with e_k white noise (unit modulus, random phase, e_(-k) = conj e_k) at every frequency. Recovery of P_b * mu by the multiplier
# exp((a - b)|omega|) on |omega| <= Omega; the cutoff Omega is the oracle one (the grid value that minimises the L2 error against the truth,
# by Parseval over the periodic coefficients, per eps). The error is the L2 error of P_b * mu (averaged over 8 noise draws in the mean square).
# The fitted slope is the least squares slope of log error against log eps over eps = 1e-12 .. 1e-3 (19 values, two per decade;
# 7 with --quick; see EPS). The predicted slope is kappa = (b + w)/(a + w) of Proposition cd:prop:finite. The levels b of the table are not stated in
# the chapter; b = 3/8, 1/4, 1/8, 0 give the predicted 0.7727, 0.5455, 0.3182, 0.0909 and are used.
# Error control. K = 16384 coefficients: the largest cutoff needed is about ln(1e12)/(a + w) = 50, k = 50 L/(2 pi) = 3200 < K, and the truncated
# tail of the truth beyond K is below exp(-(b+w) 2 pi K / L) ~ 3e-6 in amplitude at b = 0 (K = 4096 gave 0.04 and biased the b = 0 slope). All exact quantities of Part 1
# are at 30 digits.
# Output: data/cd_slopes.json. Usage: python3 cd_slopes.py [--quick] [--out FILE]; --quick uses 2 noise draws and 7 eps values, and writes to
# /tmp unless --out is given. The full case takes a few seconds.
import argparse
import json
import math
import sys
import tempfile
import os
import numpy as np
import mpmath as mp
from scipy.optimize import minimize_scalar

ap = argparse.ArgumentParser()
ap.add_argument("--quick", action="store_true")
ap.add_argument("--out", default=None)
args = ap.parse_args()
out = args.out or (os.path.join(tempfile.gettempdir(), "cd_slopes.json") if args.quick else "data/cd_slopes.json")
res = {}

# ---- Part 1
mp.mp.dps = 30
a, d = mp.mpf(1) / 2, mp.mpf(1)
q = 9 - 4 * mp.sqrt(5)
cstar = lambda q, a=a, d=d: a - mp.sqrt(((1 + q) / (1 - q) * d) ** 2 - d ** 2) / 2
m = (1 + q) / (1 - q) * d
print("q = 9 - 4 sqrt5 =", mp.nstr(q, 12), " m =", mp.nstr(m, 12), " sqrt5/2 =", mp.nstr(mp.sqrt(5) / 2, 12),
      " c* =", mp.nstr(cstar(q), 12))
print("datum bound (r-d)/(r+d) at c=0:", mp.nstr((mp.sqrt(d * d + 4 * a * a) - d) / (mp.sqrt(d * d + 4 * a * a) + d), 10),
      " 3-2sqrt2 =", mp.nstr(3 - 2 * mp.sqrt(2), 10))
print("c* at q = 0.0557:", mp.nstr(cstar(mp.mpf("0.0557")), 8))

def minlevel(c, q):  # min over x of P_(a-c)(x) - q P_(a-c)(x-d) divided by P_(a-c)(x); the ratio f(x)
    s = a - c
    f = lambda x: (s * s + (x - d) ** 2) / (s * s + x * x)
    x0 = (d + mp.sqrt(d * d + 4 * s * s)) / 2
    return f(x0) - q
lo, hi = mp.mpf(0), a - mp.mpf("1e-9")
for _ in range(100):
    mid = (lo + hi) / 2
    if minlevel(mid, q) >= 0: lo = mid
    else: hi = mid
print("bisection critical width (ratio minimum) at q = 9-4sqrt5:", mp.nstr(lo, 15))
# independent check by brute minimisation on a fine grid at c slightly below and above 1/4
xs = np.linspace(-20, 20, 400001)
def gridmin(c, qq):
    s = 0.5 - c
    return float(np.min((s * s + (xs - 1) ** 2) / (s * s + xs ** 2)) - qq)
qf = float(q)
print("grid min of the ratio minus q at c = 0.2499, 0.2501:", gridmin(0.2499, qf), gridmin(0.2501, qf))
res["part1"] = {"q": mp.nstr(q, 25), "m": mp.nstr(m, 25), "cstar_formula": mp.nstr(cstar(q), 25), "cstar_bisection": mp.nstr(lo, 20),
                "cstar_q0.0557": mp.nstr(cstar(mp.mpf("0.0557")), 10), "bound_c0": mp.nstr((mp.sqrt(d * d + 4 * a * a) - d) / (mp.sqrt(d * d + 4 * a * a) + d), 10),
                "grid_min_0.2499": gridmin(0.2499, qf), "grid_min_0.2501": gridmin(0.2501, qf)}

# ---- Part 2
L, K, w, aa = 400.0, 16384, 0.05, 0.5
rng = np.random.default_rng(12345)
nat = 40
xj = rng.uniform(0, L, nat)
mj = rng.uniform(0.5, 1.5, nat)
k = np.arange(0, K + 1)
om = 2 * np.pi * k / L
muk = (mj[None, :] * np.exp(-1j * om[:, None] * xj[None, :])).sum(axis=1) / L   # Fourier coefficients k >= 0 (real source: k < 0 conjugates)
wt = np.where(k == 0, 1.0, 2.0)                                                  # both signs enter the norm
EPS = 10.0 ** np.arange(-12, -2.9, 0.5 if not args.quick else 1.5)
ndraw = 2 if args.quick else 8
levels = [0.375, 0.25, 0.125, 0.0]
res["part2"] = {"L": L, "K": K, "w": w, "a": aa, "atoms": {"x": xj.tolist(), "m": mj.tolist()}, "eps": EPS.tolist(), "levels": {}}
print("eps range", EPS[0], EPS[-1], "points", len(EPS), "draws", ndraw)
for b in levels:
    true_b = np.exp(-(b + w) * om) * muk       # coefficients of P_b * h, h = P_w * mu
    mult = np.exp((aa - b) * om)
    errs = []
    for eps in EPS:
        e2 = 0.0
        for _ in range(ndraw):
            ph = np.exp(2j * np.pi * rng.uniform(size=K + 1)); ph[0] = 1.0
            data = np.exp(-(aa + w) * om) * muk + eps * ph
            rec = mult * data
            diff = np.abs(rec - true_b) ** 2 * wt
            tail = np.abs(true_b) ** 2 * wt
            # error for cutoff index n: sum_{k<=n} diff + sum_{k>n} tail; take the minimum over n (oracle)
            c1 = np.cumsum(diff); c2 = np.append(np.cumsum(tail[::-1])[::-1][1:], 0.0)
            e2 += (c1 + c2).min()
        errs.append(math.sqrt(e2 / ndraw * L))   # L2 norm over one period: L * sum |coef|^2
    errs = np.array(errs)
    sl = np.polyfit(np.log(EPS), np.log(errs), 1)[0]
    pred = (b + w) / (aa + w)
    print("b = %.3f  predicted %.4f  fitted %.4f" % (b, pred, sl))
    res["part2"]["levels"][str(b)] = {"predicted": pred, "fitted": float(sl), "errors": errs.tolist()}
json.dump(res, open(out, "w"), indent=1)
print("written", out)
