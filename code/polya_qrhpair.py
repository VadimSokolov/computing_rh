# polya_qrhpair.py: the numbers behind Theorem pol:thm:qrhpair (the compensated clock under Q_{7/8}, paired form) of Section pol:sec:qrh of Chapter 16 (Polya and the Zeros of Entire Functions). Added for the book in October 2026; it does not check a claim of the chapter author, whose argument does not pair a zero with its mirror image.
# The theorem. Assume that zeta has no zero with Re s > alpha_0 (alpha_0 = 7/8 in the chapter). For alpha_0 < alpha <= 1 write a = alpha - 1/2, c_max = alpha + alpha_0 - 1 < 1, f(c, x) = c/(c^2 + x^2) and D_alpha(r) = Re (xi'/xi)(alpha + ir) - (xi'/xi)(alpha) = sum over all zeros rho of f(alpha - beta, r - gamma) - b_alpha, b_alpha = (xi'/xi)(alpha). Every zero has alpha - beta >= alpha - alpha_0 > 0, and the zeros with ordinate in [r - d, r + d] come in pairs rho, 1 - conj(rho) with the same ordinate and widths c, c' with c + c' = 2a, c, c' <= c_max, so that f(c, x) + f(c', x) >= 2a/(c_max^2 + d^2) for |x| <= d; a zero on the line is its own mirror and f(a, x) >= a/(c_max^2 + d^2). Hence D_alpha(r) >= n a/(c_max^2 + d^2) - b_alpha, n the number of zeros in the window, and b_alpha <= 0.04621 a (Lemma pol:lem:Bbound), so D_alpha(r) > 0 as soon as n > 0.04621 (c_max^2 + d^2); for d = 4 and c_max < 1 this is n >= 1 (0.04621 * 17 = 0.786). For |r| <= theta_0 = gamma_1 - alpha/sqrt(3) every term f(c, r - gamma) is strictly convex in r (|r - gamma| >= alpha/sqrt(3) > c/sqrt(3)), so D_alpha is strictly convex and even and positive for r != 0 (the argument of Proposition bp:prop:v1min).
# What this script checks. (1) The pair inequality f(c, x) + f(2a - c, x) >= 2a/(c_max^2 + d^2) for alpha_0 = 7/8, d = 4, on grids of alpha, c and x, and the margin n a/(c_max^2 + d^2) - b_alpha for n = 1 with b_alpha computed by mpmath (30 digits) at 60 values of alpha in (7/8, 1], and the sufficient test 0.04621 (c_max^2 + d^2) < 1. (2) The window of half width 4 around every r from gamma_1 - 4 to gamma_1700 + 4 contains an ordinate of data/zeros_1700.txt: the largest gap between consecutive ordinates is at most 8, and the convexity radius theta_0(1) = gamma_1 - 1/sqrt(3) is at least gamma_1 - 4. (3) Lemma pol:lem:window with d = 4, variant H of polya_qrh.py (s(t) = 2.5167 up to X_0 = 30 610 046 000, then 0.1095 log t + 0.2042 log log t + 3.0305): nu(r, 4) > 0 on a geometric grid of 4000 points per decade from r = 100 to 1e16, with the smallest r from which nu stays positive on the grid, the minimum of nu above it, and a recomputation of nu at ten heights with mpmath (siegeltheta, 40 digits); the theorem uses the lemma only from r = 2200, where nu = 2.43 (the grid has nu > 0 from about 327.5 on, with the smallest value 5.8e-4 just above that height); beyond 1e16 the first term (4/pi) log(r/2 pi) grows faster than the second (0.219 log r + 0.41 log log r + 6.07) and nu increases. (4) The window of Lemma pol:lem:window against the real zeros: nu(r, 4) never exceeds the count of ordinates in [r - 4, r + 4] for r from 40 to 2190.
# Output: data/polya_qrhpair.json (numbers as strings) and the log data/polya_qrhpair.log. Usage: python3 polya_qrhpair.py [--out FILE] [--zeros FILE]; it takes about 10 s in one process, so it was run locally and not on the Hopper cluster.
import argparse
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

T0 = time.time()
LOG = []
PI = math.pi
X0 = 30610046000.0
KS = 7
CK = [(1 - 2.0 ** (1 - 2 * k)) * abs(float(mp.bernoulli(2 * k))) / (4 * k * (2 * k - 1)) for k in range(1, KS + 1)]
D = 4.0
ALPHA0 = 7.0 / 8.0


def log(*a):
    s = '[%6.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(s, flush=True)
    LOG.append(s)


def dtheta(r, d):
    # theta(r + d) - theta(r - d) from the Stirling series, valid for r - d >= 30
    def c(t):
        return sum(CK[k - 1] * t ** (-(2 * k - 1)) for k in range(1, KS + 1))
    return 0.5 * (r * np.log1p(2 * d / (r - d)) + d * np.log((r * r - d * d) / (4 * PI * PI))) - d + c(r + d) - c(r - d)


def sH(t):
    t = np.asarray(t, dtype=float)
    big = 0.1095 * np.log(t) + 0.2042 * np.log(np.log(t)) + 3.0305
    return np.where(t <= X0, 2.5167, big)


def nu(r, d=D):
    return dtheta(r, d) / PI - sH(r + d) - sH(r - d)


def find_zeros(arg):
    here = Path(__file__).resolve().parent
    for c in ([arg] if arg else []) + [here / 'zeros_1700.txt', Path('data/zeros_1700.txt'), here.parent / 'data' / 'zeros_1700.txt']:
        if c and Path(c).is_file():
            return Path(c)
    raise SystemExit('zeros_1700.txt not found; pass --zeros')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None)
    ap.add_argument('--zeros', default=None)
    args = ap.parse_args()
    out = {}
    mp.mp.dps = 30

    # (1) the pair inequality and the margin for one zero in the window
    worst = 1e9
    for alpha in np.linspace(ALPHA0 + 1e-4, 1.0, 120):
        a = alpha - 0.5
        cmax = alpha + ALPHA0 - 1.0
        cmin = alpha - ALPHA0
        cs = np.linspace(cmin, cmax, 241)
        xs = np.linspace(0.0, D, 321)
        C, X = np.meshgrid(cs, xs)
        f1 = C / (C * C + X * X)
        C2 = 2 * a - C
        f2 = C2 / (C2 * C2 + X * X)
        g = (f1 + f2) / (2 * a / (cmax * cmax + D * D))
        worst = min(worst, float(g.min()))
    out['pair_inequality_min_ratio'] = repr(worst)
    log('pair inequality: smallest ratio (f(c,x) + f(2a - c,x)) / (2a/(c_max^2 + d^2)) over the grids =', worst)

    def xi(s):
        return 0.5 * s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

    def b_alpha(alpha):
        return mp.diff(lambda s: mp.log(xi(s)), mp.mpf(alpha))

    mins = []
    ratio_max = 0.0
    rows = []
    for alpha in [0.8751 + k * (1.0 - 0.8751) / 59 for k in range(60)]:
        a = alpha - 0.5
        cmax = alpha + ALPHA0 - 1.0
        b = b_alpha(alpha)
        margin = a / (cmax * cmax + D * D) - float(b)
        mins.append((margin, alpha))
        ratio_max = max(ratio_max, float(b) / a)
        if abs(alpha - 0.8751) < 1e-12 or abs(alpha - 1.0) < 1e-12 or len(rows) == 29:
            rows.append((alpha, float(b), margin))
    mm = min(mins)
    out['margin_n1_min'] = repr(mm[0])
    out['margin_n1_min_at_alpha'] = repr(mm[1])
    out['b_over_a_max'] = repr(ratio_max)
    out['sufficient_test_max'] = repr(0.04621 * ((1.0 + ALPHA0 - 1.0) ** 2 + D * D))
    out['sufficient_test_at_alpha_1'] = repr(0.04621 * ((1.0 + ALPHA0 - 1.0) ** 2 + D * D))
    out['rows_alpha_b_margin'] = [[repr(x) for x in r] for r in rows]
    log('one zero in the window: smallest of a/(c_max^2 + d^2) - b_alpha over 60 values of alpha in (7/8, 1] =', mm[0], 'at alpha =', mm[1])
    log('largest b_alpha/a over the same values =', ratio_max, '(Lemma pol:lem:Bbound: <= 0.04621)')
    log('sufficient test 0.04621 (c_max^2 + d^2) at alpha = 1 =', 0.04621 * ((1.0 + ALPHA0 - 1.0) ** 2 + D * D), '< 1')

    # (2) gaps between ordinates and the convexity radius
    zs = [mp.mpf(l.strip()) for l in find_zeros(args.zeros).read_text().split() if l.strip()]
    z = np.array([float(x) for x in zs])
    gaps = np.diff(z)
    out['n_zeros'] = len(z)
    out['gamma_1'] = repr(float(z[0]))
    out['gamma_last'] = repr(float(z[-1]))
    out['max_gap'] = repr(float(gaps.max()))
    out['max_gap_between'] = [repr(float(z[int(gaps.argmax())])), repr(float(z[int(gaps.argmax()) + 1]))]
    second = np.sort(gaps)[-2]
    out['second_gap'] = repr(float(second))
    theta0 = float(z[0]) - 1.0 / math.sqrt(3.0)
    out['theta0_at_alpha_1'] = repr(theta0)
    out['window_start'] = repr(float(z[0]) - D)
    log('zeros:', len(z), 'ordinates from', z[0], 'to', z[-1], '; largest gap', gaps.max(), 'between', z[int(gaps.argmax())], 'and', z[int(gaps.argmax()) + 1], '; second largest', second)
    log('convexity radius theta_0(1) = gamma_1 - 1/sqrt(3) =', theta0, '; windows of half width 4 start at', z[0] - D, '; need gap <= 8:', gaps.max() <= 2 * D, 'and theta_0 >= start:', theta0 >= z[0] - D)

    # (3) Lemma window with d = 4, variant H
    per = 4000
    lg = np.linspace(2.0, 16.0, int(14 * per) + 1)
    rs = 10 ** lg
    rs = rs[rs - D >= 30.0]
    vals = nu(rs)
    pos = vals > 0
    # smallest r from which nu stays positive on the grid
    bad = np.where(~pos)[0]
    start = float(rs[bad[-1] + 1]) if len(bad) else float(rs[0])
    out['nu_positive_from'] = repr(start)
    out['nu_min_above_start'] = repr(float(vals[bad[-1] + 1:].min()) if len(bad) else float(vals.min()))
    out['nu_min_from_2200'] = repr(float(vals[rs >= 2200.0].min()))
    out['nu_at_335'] = repr(float(nu(np.array([335.0]))[0]))
    out['nu_at_2200'] = repr(float(nu(np.array([2200.0]))[0]))
    out['nu_at_X0_minus_d'] = repr(float(nu(np.array([X0 - D - 1.0]))[0]))
    out['nu_at_X0_plus_d'] = repr(float(nu(np.array([X0 + D + 1.0]))[0]))
    out['nu_at_1e16'] = repr(float(nu(np.array([1e16]))[0]))
    log('nu(r, 4) > 0 on the grid from r =', start, '; smallest value above it', out['nu_min_above_start'], '; smallest value from r = 2200:', out['nu_min_from_2200'], '; nu(335) =', out['nu_at_335'], '; nu(2200) =', out['nu_at_2200'], '; just below and above X_0:', out['nu_at_X0_minus_d'], out['nu_at_X0_plus_d'], '; nu(1e16) =', out['nu_at_1e16'])
    # recomputation of theta(r + d) - theta(r - d) with mpmath
    mp.mp.dps = 40
    worstdiff = 0.0
    for r in [100.0, 335.0, 1000.0, 5000.0, 1e5, 1e7, 1e9, 1e11, 1e12, 1e14]:
        ex = (mp.siegeltheta(mp.mpf(r) + D) - mp.siegeltheta(mp.mpf(r) - D))
        ap_ = dtheta(r, D)
        worstdiff = max(worstdiff, abs(float(ex) - float(ap_)))
    out['theta_difference_max_error'] = repr(worstdiff)
    log('largest difference between the Stirling form of theta(r + 4) - theta(r - 4) and mpmath at ten heights:', worstdiff)

    # (4) the window lower bound against the actual zeros
    over = 0.0
    for r in np.arange(40.0, 2190.0, 0.5):
        cnt = int(((z >= r - D) & (z <= r + D)).sum())
        over = max(over, float(nu(np.array([r]))[0]) - cnt)
    out['nu_minus_count_max'] = repr(over)
    log('largest nu(r, 4) - (number of ordinates in [r - 4, r + 4]) for r from 40 to 2190 step 0.5:', over, '(must be <= 0)')

    out['seconds'] = repr(time.time() - T0)
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        Path(args.out).with_suffix('.log').write_text('\n'.join(LOG) + '\n')
        log('wrote', args.out)


if __name__ == '__main__':
    main()
