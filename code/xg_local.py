# code/xg_local.py (October 2026, from misc/reviews/nick-xiggc/: review-A1.md, the row prop:local, errors E8 and E9, item (d)11): the local model for the
# first failing order of Proposition xg:prop:local (Section sec:xg:local, the end of Section 5.2) against the certified first failing orders and times of
# Table tab:ch13:first. Written from the mathematics, not from the authors' scripts.
# Setting. The pair of consecutive zeros g -+ (B-A)/2 of a case of Table tab:ch13:first is replaced by the quadruple 1/2 +- delta +- i g, g the midpoint, and
# c_k(t) = sum a^k e^{-at} over one zero of each pair {rho, 1-rho}; the quadruple contributes 2 Re a^k e^{-at} with a = (g - i delta)^2, a zero on the line
# at height gamma contributes gamma^{2k} e^{-gamma^2 t}. Let g - Dm and g + Dp be the nearest zeros on the line below and above the moved pair (the zero
# numbers n1 - 1 and n2 + 1 of zeros_1700.txt) and D = max(Dm, Dp). With t = k(1+v)/g^2 and v = -eta D/g, eta in [0,1], the weights relative to g^{2k} e^{-g^2 t} are
#   quadruple: 2 cos(omega eta) (1 + O(delta/D)),  omega = 2 k delta D/g^2,   neighbour: exp(-kappa (1 - eta)) (1 + ...),  kappa = 2 k D^2/g^2,
# so the order k fails, to leading order, iff min over eta in [0,1] of h(eta) = 2 cos(omega eta) + exp(-kappa (1 - eta)) is negative. Put r = D/delta = kappa/omega.
# h decreases in omega on (0, pi] for every eta, so the least failing omega is the root of min h = 0, found here by bisection with a fine grid in eta (check 1).
# At the root h = 0 and h' = 0, hence tan(omega eta) = -r and exp(-kappa(1-eta)) = 2/sqrt(1+r^2), which gives the closed form
#   omega_c(r) = pi - arctan(r) + log(sqrt(1+r^2)/2)/r  (a second check, not used in the book),
# and its expansion for large r is omega* = pi/2 + (1 + log(r/2))/r + O(r^-3), the value stated in the Proposition (check 1: omega* against omega_c).
# The first failing order of the model is k_loc = omega* g^2/(2 delta D) and the time t_loc = omega*/(2 delta D) = k_loc/g^2. These are compared with the
# certified k* and t* of data/firstfail_arb.json, the file that certifies Table tab:ch13:first (the zero numbers, delta, k*, t* of its seven cases).
# (2) Truncation to three zeros, not used in the Proposition: the least order k_3 at which the exact sum S(k,t) = 2 Re a^k e^{-at} + gm^{2k} e^{-gm^2 t} +
# gp^{2k} e^{-gp^2 t} (quadruple and the two neighbours, exact weights, divided by g^{2k} e^{-g^2 t}) is negative for some t > 0. Failure is monotone in k,
# since c_{k-1}(t) = int_t^oo c_k(s) ds termwise, so k_3 is found by bisection on k; for each k the minimum over t = k/u^2 comes from a grid in u of step 0.001
# between g - 1.2 Dm and g + 1.2 Dp, refined by Brent's method at every local minimum (double precision, with the large exponents formed in mpmath first;
# the margin of S per unit of k is about 1e-7, far above the rounding). It separates the error of the expansion (k_loc against k_3) from the error of
# the truncation (k_3 against the certified k*). Not a certificate: the certified values are those of firstfail_arb.json.
# (3) The model drops the factor exp(-2 t delta^2) = exp(-omega delta/D) of the weight of the neighbour relative to the quadruple (the amplification
# e^{2t delta^2} of the quadruple in the Gaussian smoothing of Section sec:ch12:ray). Keeping it, the tangency analysis of (1) gives omega_c r^2/(1+r^2),
# checked here by bisection as in (1); it lowers the model times in every row, which explains the sign of the differences. Not used in the Proposition.
# Run from a folder that holds zeros_1700.txt and firstfail_arb.json: python3 xg_local.py (seconds, one core); writes xg_local.json there.
import json
import math
import os
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar

mp.mp.dps = 50
t0 = time.time()


def say(*a):
    print(*a, flush=True)


def find(name):
    for p in (name, 'data/' + name, '../data/' + name):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)


Z = [mp.mpf(s) for s in open(find('zeros_1700.txt')).read().split()]
cases = json.load(open(find('firstfail_arb.json')))['cases']
assert len(Z) == 1700 and len(cases) == 7


# ---------------------------------------------------------------- 1. the model function and its threshold
def omega_star(r):                      # the Proposition: pi/2 + (1 + log(r/2))/r
    return math.pi / 2 + (1 + math.log(r / 2)) / r


def omega_closed(r):                    # exact threshold of the model function (tangency), r > sqrt(3)
    return math.pi - math.atan(r) + math.log(math.sqrt(1 + r * r) / 2) / r


ETA = np.linspace(0.0, 1.0, 100001)


def hmin(omega, r, keep=False):         # min over eta in [0,1] of 2 cos(omega eta) + exp(-r omega (1-eta)) [times exp(-omega/r) if keep]
    return float(np.min(2 * np.cos(omega * ETA) + np.exp(-r * omega * (1 - ETA) - (omega / r if keep else 0.0))))


def omega_bisect(r, keep=False):        # least omega in (pi/2, pi] with min h < 0
    lo, hi = math.pi / 2, math.pi
    assert hmin(lo, r, keep) > 0 and hmin(hi, r, keep) < 0
    for _ in range(60):
        mid = (lo + hi) / 2
        if hmin(mid, r, keep) < 0:
            hi = mid
        else:
            lo = mid
    return hi


# ---------------------------------------------------------------- 2. the three zero model, exact weights
def make_smin(g, d, gm, gp):
    a = mp.mpc(g, -d) ** 2
    Dm, Dp = g - gm, gp - g
    ulo, uhi = float(g - mp.mpf('1.2') * Dm), float(g + mp.mpf('1.2') * Dp)
    ugrid = np.arange(ulo, uhi, 0.001)

    def smin(k):
        La = complex(k * mp.log(a / g ** 2))
        Lm, Lp = float(2 * k * mp.log(gm / g)), float(2 * k * mp.log(gp / g))
        da, dmm, dpp = complex(a - g ** 2), float(gm ** 2 - g ** 2), float(gp ** 2 - g ** 2)

        def S(u):
            t = k / (u * u)
            return 2 * np.real(np.exp(La - da * t)) + np.exp(Lm - dmm * t) + np.exp(Lp - dpp * t)

        with np.errstate(over='ignore', under='ignore'):
            s = S(ugrid)
            best = (float(np.min(s)), float(ugrid[int(np.argmin(s))]))
            for i in range(1, len(ugrid) - 1):
                if s[i] <= s[i - 1] and s[i] <= s[i + 1] and s[i] < 0.05:
                    res = minimize_scalar(S, bounds=(ugrid[i - 1], ugrid[i + 1]), method='bounded', options={'xatol': 1e-13})
                    if res.fun < best[0]:
                        best = (float(res.fun), float(res.x))
        return best[0], best[1], k / best[1] ** 2

    return smin


def first_fail(smin, k_guess):
    lo, hi = int(0.8 * k_guess), int(1.25 * k_guess)
    assert smin(lo)[0] > 0 and smin(hi)[0] < 0, (lo, hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if smin(mid)[0] < 0:
            hi = mid
        else:
            lo = mid
    return hi


# ---------------------------------------------------------------- 3. the seven rows of Table tab:ch13:first
rows = []
say('Local model of Proposition xg:prop:local against Table tab:ch13:first; D = larger gap from g to the zeros next to the moved pair')
say('row  zeros        g          delta   Dm      Dp      D      r=D/delta  omega*   omega_bis  omega_closed   t_loc    t*cert   diff t   k_loc      k*cert    diff k   k_3        diff k_3  diff k_3/k*')
for i, c in enumerate(cases, 1):
    n1, n2 = c['zero_numbers']
    assert n2 == n1 + 1
    fr = Fraction(c['delta'])
    d = mp.mpf(fr.numerator) / fr.denominator
    A, B = Z[n1 - 1], Z[n2 - 1]
    g = (A + B) / 2
    gm, gp = Z[n1 - 2], Z[n2]
    assert gm < A < B < gp
    Dm, Dp = g - gm, gp - g
    D = max(Dm, Dp)
    r = float(D / d)
    w_star, w_bis, w_cl = omega_star(r), omega_bisect(r), omega_closed(r)
    assert abs(w_bis - w_cl) < 1e-8, (w_bis, w_cl)
    w_ref, w_ref_bis = w_cl * r * r / (1 + r * r), omega_bisect(r, keep=True)      # the factor exp(-2 t delta^2) kept (check 3)
    assert abs(w_ref - w_ref_bis) < 1e-8, (w_ref, w_ref_bis)
    t_loc = w_star / float(2 * d * D)
    t_ref = w_ref / float(2 * d * D)
    k_loc = t_loc * float(g) ** 2
    t_cert, k_cert = float(mp.mpf(c['t'])), c['k']
    diff_t, diff_k, diff_ref = t_loc / t_cert - 1, k_loc / k_cert - 1, t_ref / t_cert - 1
    smin = make_smin(g, d, gm, gp)
    k3 = first_fail(smin, k_cert)
    s3, u3, tt3 = smin(k3)
    assert smin(k3 - 1)[0] > 0 and smin(k3 + 1)[0] < 0 and smin(int(k3 * 1.01))[0] < 0
    side = 'above g' if u3 > g else 'below g'
    rows.append(dict(
        row=i, window=c['window'], zero_numbers=[n1, n2], delta=float(d), delta_string=c['delta'], g=float(g),
        gap_below=float(Dm), gap_above=float(Dp), gap_larger=float(D), ratio_gap_over_delta=r,
        omega_star=w_star, omega_bisection=w_bis, omega_closed_form=w_cl, omega_star_minus_exact=w_star - w_cl,
        t_loc=t_loc, t_certified=t_cert, relative_difference_t=diff_t,
        omega_refined=w_ref, t_refined=t_ref, relative_difference_t_refined=diff_ref,
        k_loc=k_loc, k_certified=k_cert, relative_difference_k=diff_k,
        k_three_zero=k3, t_three_zero=tt3, u_three_zero_minus_g=float(u3 - g), three_zero_side=side,
        relative_difference_k_three_zero_to_certified=k3 / k_cert - 1, relative_difference_k_loc_to_three_zero=k_loc / k3 - 1,
        min_S_at_k3=s3))
    say('%d  %4d,%4d  %9.4f  %.2f  %.4f  %.4f  %.4f  %7.3f  %.6f  %.6f  %.6f  %.5f  %.5f  %+.4f%%  %9.0f  %9d  %+.4f%%  %9d  %+.4f%%  %+.4f%%' % (
        i, n1, n2, float(g), float(d), float(Dm), float(Dp), float(D), r, w_star, w_bis, w_cl, t_loc, t_cert, 100 * diff_t,
        k_loc, k_cert, 100 * diff_k, k3, 100 * (k_loc / k3 - 1), 100 * (k3 / k_cert - 1)))
    say('     three zero model: first failing order %d at t = %.5f, u - g = %+.4f (%s); S there %.3e; certified order %d at t* = %s' % (
        k3, tt3, float(u3 - g), side, s3, k_cert, c['t']))

diffs = [100 * x['relative_difference_t'] for x in rows]
say('t_loc against the certified t*: relative differences in per cent %s; smallest %.2f, largest %.2f, all positive: %s' % (
    ['%+.2f' % x for x in diffs], min(abs(x) for x in diffs), max(abs(x) for x in diffs), all(x > 0 for x in diffs)))
dk = [100 * x['relative_difference_k'] for x in rows]
say('k_loc = t_loc g^2 against the certified k*: per cent %s; largest modulus %.2f' % (['%+.2f' % x for x in dk], max(abs(x) for x in dk)))
d3 = [100 * x['relative_difference_k_three_zero_to_certified'] for x in rows]
say('three zero model against the certified k*: per cent %s; largest modulus %.4f' % (['%+.4f' % x for x in d3], max(abs(x) for x in d3)))
say('omega* against the exact threshold of the model function: largest difference %.2e' % max(abs(x['omega_star_minus_exact']) for x in rows))
dr = [100 * x['relative_difference_t_refined'] for x in rows]
say('with the factor exp(-2 t delta^2) kept, omega = omega_c r^2/(1+r^2): t against the certified t*, per cent %s; largest modulus %.2f' % (['%+.2f' % x for x in dr], max(abs(x) for x in dr)))

json.dump(dict(script='xg_local.py', source_certified='firstfail_arb.json (Table tab:ch13:first)', zeros='zeros_1700.txt',
               model='t_loc = omega*/(2 delta D), omega* = pi/2 + (delta/D)(1 + log(D/(2 delta))), D the larger gap from g to the zeros next to the moved pair',
               rows=rows, seconds=round(time.time() - t0, 1)), open('xg_local.json', 'w'), indent=1)
say('wrote xg_local.json in %.1f s' % (time.time() - t0))
