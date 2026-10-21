# yo_symbol.py: group Y4 of the numerical checks of Chapter 17 (Yoshida's Hermitian forms, Weil positivity and the Thorin program), Section yo:sec:damped (damped forms and sharp coercivity): the symbol v_alpha(tau) = Re xi'/xi(alpha + i tau), the Laguerre type quantity L(tau), the constants kappa_1 and b_alpha and the margins of the proof of Theorem yo:thm:coercive (ii). Added for the book in October 2026, written from the mathematics (the chapter's own scripts were neither read nor run).
# Claims checked (chapter line numbers of ch_yoshida.tex of 17 October 2026). (a) Figure yo:fig:symbol (line 329) and item 1 of Section yo:sec:path (line 394): the minimum of v_alpha on [20, 2000] for alpha = 1, 0.75, 0.6, 0.55 (printed 0.190, 0.097, 0.039 and no value for the fourth), the ratio to b_alpha = xi'/xi(alpha) ("about eight"), where each minimum is attained, and that on [0, 2000] the minimum is at tau = 0. (b) Lines 296 to 298 and 395: L(tau) = sum over all zeros of (tau - gamma)^(-2) = -(log Xi)''(tau), L(0) = 2 kappa_1 = 0.046210, whether the global minimum on [0, 1800] is at 0, the next local minimum 0.2301 at tau = 17.50 and the ratio 4.98 ("five times"). (c) Line 278: 1/(c^2 + d^2) at c = 1/2 for d = 4.15 (0.057233) and d = 3.4437 (0.082583), 2/gamma_1^2 = 0.010010, 2 kappa_1 = 0.046210 and the margins 0.001013 and 0.026362, and the margin at other c (line 294: "wide for c near 0"), with the algebraic lemma of line 274 (P_gamma(tau) >= P_gamma(0) iff tau^2 <= 3 gamma^2 - c^2) tested on 200000 random triples. (d) kappa_1 = sum over gamma > 0 of gamma^(-2) = 0.0231049931, xi'/xi(1) = 1 + gamma_E/2 - (1/2) log 4 pi = 0.0230957, and the floor row 2 b_alpha of Table yo:tab:lam (0.0046210, 0.0092418, 0.0231027, 0.0461914 at alpha = 0.55, 0.6, 0.75, 1), with the two sided bound of Theorem yo:thm:coercive (iv) on it.
# Method. xi'/xi(s) = 1/s + 1/(s - 1) - (1/2) log pi + (1/2) psi(s/2) + zeta'/zeta(s) is evaluated directly with mpmath (the tail free value); the Hadamard sums are evaluated in float64 over the first 1700 ordinates of data/zeros_1700.txt (68 digits, all on the line, to height 2197.26) with all zeros of both signs, v_alpha(tau) = sum over gamma > 0 of k_c(tau - gamma) + k_c(tau + gamma), k_c(x) = c/(c^2 + x^2), c = alpha - 1/2, and L(tau) the same with x^(-2). The zeros beyond the 1700th are replaced by their mean density, the Riemann von Mangoldt tail: the integral from T0 to infinity of nbar(g) (k_c(g - tau) + k_c(g + tau)) dg with nbar = theta'/pi = (Re psi(1/4 + i g/2) - log pi)/(2 pi) and T0 = the Gram point g_1699 (theta(T0) = 1699 pi, so that N(T0) = 1700 = theta(T0)/pi + 1 and S(T0) = 0; mpmath.nzeros confirms the count), integrated by 24 point Gauss Legendre panels (ratio 2.5) out to 1e14. v_alpha is scanned on the grid of step c/10 (c/4 with --divisor 4) over [0, tau_max]; the 30 lowest grid minima above 20 and the lowest interior one are refined by Brent's method on the bracket of the grid. L is convex between consecutive ordinates (each term is), so it has exactly one minimum in each gap, found by Brent's root finder on L'; the first interval (-gamma_1, gamma_1) has its minimum at 0 and L increases on [0, gamma_1); a dense grid scan of L (step 0.01) confirms. kappa_1, kappa_2 are the Taylor coefficients of xi'/xi at 1/2, 2 kappa_1 z - 2 kappa_2 z^3 + ..., extracted by Cauchy's formula on circles of radius 5 and 8 (xi'/xi is analytic in |z| < 14.13), and are also summed over the zeros with the tail.
# Error control. (1) The Hadamard sum with the tail is compared with the direct mpmath value (25 to 30 digits) at the located minima, at tau = 0 and at points up to tau_max, for every alpha; the largest difference is the error of the sums (the tail beyond the 1700th ordinate has size about 1e-3 c at height 2000 and its error, which is the integral of S(g) against the derivative of the kernel, is of order 1e-9 c). (2) The tail is rebuilt from the first 1200 ordinates and its Gram point g_1199 and the two totals are compared. (3) The quadrature of the tail is compared with mpmath.quad. (4) kappa_1 from the two circles, from the zero sum with the tail and from the book's V = (log xi)''(1/2) in data/dbeta_minimum.json; b_1 from the closed form, from the formula and from riesz_program/data/v1min.json; the sample of ordinates is checked against |zeta(1/2 + i gamma)|; v_1(17.6), the minimum of v_1 on the grid of [13.5, 2000], the largest gap, the largest later gap and the largest distance to an ordinate are compared with riesz_program/data/v1min.json. (5) Nothing here is a proof: the scan is a grid. What the book certifies is recorded in the key book_results (Propositions bp:prop:v1min and bp:prop:valpha, Theorem thm:dbeta:minall (ii), scripts riesz_program/code/v1min.py, repro/valpha_min.py, code/dbeta_minimum.py), and the claim of the figure that v_alpha > b_alpha for tau != 0 is a theorem of the book on the whole plotted range for the four alpha, whereas the values of the minima above tau = 20 are numerics.
# Output: data/yo_symbol.json (values as strings; parameters, zero counts, the Gram point, the minima, the direct checks, the margins, the data of the figure), data/yo_symbol.log and fig/yo_symbol.pdf (the figure, redrawn from the figure data of the JSON; the paperback is black and white, so the four curves differ by line style and grey level). With --quick the files go to /tmp unless --out and --fig are given.
# Usage: python3 yo_symbol.py [--quick] [--out FILE] [--fig FILE] [--zeros FILE] [--root DIR] [--stages constants,symbol,laguerre,figure] [--alphas 1,0.75,0.6,0.55] [--tau-max T] [--divisor D] [--update]. --quick runs all stages on [0, 300] with divisor 4 in about 4 seconds in one process. The full case (tau_max 2000, divisor 10, L on [0, 1800]) takes about 20 seconds in one process on a laptop (the scan of alpha = 0.55 with 400001 grid points is half of it), so no cluster job is needed; it writes data/yo_symbol.json, data/yo_symbol.log and fig/yo_symbol.pdf from the repository root. --stages with --update runs single stages and merges them into an existing JSON.
import argparse
import json
import math
import os
import platform
import re
import socket
import sys
import time
from fractions import Fraction

import numpy as np
import mpmath
from mpmath import mp, mpf, mpc
from scipy.optimize import brentq, minimize_scalar
from scipy.special import digamma

HERE = os.path.dirname(os.path.abspath(__file__))
T_START = time.time()
LOGFILE = None
ALPHAS_ALL = ['1', '0.75', '0.6', '0.55']
CHAPTER_MIN = {'1': '0.190', '0.75': '0.097', '0.6': '0.039', '0.55': None}
CHAPTER_FLOOR = {'0.55': '0.0046210', '0.6': '0.0092418', '0.75': '0.0231027', '1': '0.0461914'}


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T_START, msg)
    print(line, flush=True)
    if LOGFILE is not None:
        LOGFILE.write(line + '\n')
        LOGFILE.flush()


def fmt(x, n=10):
    if isinstance(x, (mpf, mpc)):
        return mpmath.nstr(x, n)
    return '%.*g' % (n, float(x))


def cval(alpha):
    return Fraction(alpha) - Fraction(1, 2)


# ---------------------------------------------------------------------------------------------------------------
# The zero sums with the Riemann von Mangoldt tail
# ---------------------------------------------------------------------------------------------------------------

def nbar(g):
    # mean density of ordinates, theta'(g)/pi with theta'(g) = (1/2) Re psi(1/4 + i g/2) - (1/2) log pi
    return (digamma(0.25 + 0.5j * np.asarray(g)).real - math.log(math.pi)) / (2.0 * math.pi)


def tail_nodes(T0, order=24):
    xs, ws = np.polynomial.legendre.leggauss(order)
    edges = [0.0, 100.0]
    while edges[-1] < 1e14:
        edges.append(edges[-1] * 2.5)
    G, W = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        G.append(T0 + 0.5 * (b - a) * xs + 0.5 * (b + a))
        W.append(0.5 * (b - a) * ws)
    g, w = np.concatenate(G), np.concatenate(W)
    return g, w * nbar(g)


class Model:
    # the first n ordinates summed exactly, the rest replaced by the mean density from the Gram point g_{n-1}
    def __init__(self, gam, n):
        self.n = n
        self.gam = np.asarray(gam[:n], dtype=float)
        self.T0 = float(mpmath.grampoint(n - 1))
        self.count_T0 = int(mpmath.nzeros(self.T0))
        self.nodes, self.wts = tail_nodes(self.T0)

    def zero_part(self, kind, tau, c=None):
        tau = np.atleast_1d(np.asarray(tau, dtype=float))
        out = np.empty(tau.shape)
        rows = max(1, 3000000 // self.n)
        g = self.gam[None, :]
        with np.errstate(divide='ignore'):
            for i in range(0, len(tau), rows):
                t = tau[i:i + rows, None]
                a, b = t - g, t + g
                if kind == 'v':
                    c2 = c * c
                    out[i:i + rows] = c * (1.0 / (c2 + a * a) + 1.0 / (c2 + b * b)).sum(axis=1)
                elif kind == 'L':
                    out[i:i + rows] = (1.0 / (a * a) + 1.0 / (b * b)).sum(axis=1)
                else:
                    out[i:i + rows] = -2.0 * (1.0 / (a * a * a) + 1.0 / (b * b * b)).sum(axis=1)
        return out

    def tail_part(self, kind, tau, c=None):
        tau = np.atleast_1d(np.asarray(tau, dtype=float))
        assert tau.max() < self.T0 - 100.0, 'tau too close to the cut of the tail'
        out = np.empty(tau.shape)
        rows = max(1, 3000000 // len(self.nodes))
        g, w = self.nodes[None, :], self.wts[None, :]
        for i in range(0, len(tau), rows):
            t = tau[i:i + rows, None]
            a, b = g - t, g + t
            if kind == 'v':
                c2 = c * c
                F = c / (c2 + a * a) + c / (c2 + b * b)
            elif kind == 'L':
                F = 1.0 / (a * a) + 1.0 / (b * b)
            else:
                F = 2.0 / (a * a * a) - 2.0 / (b * b * b)
            out[i:i + rows] = (F * w).sum(axis=1)
        return out

    def total(self, kind, tau, c=None):
        return self.zero_part(kind, tau, c) + self.tail_part(kind, tau, c)

    def kappa_tail(self):
        return float(np.sum(self.wts / self.nodes ** 2))


# ---------------------------------------------------------------------------------------------------------------
# Direct (tail free) evaluations with mpmath
# ---------------------------------------------------------------------------------------------------------------

def xi_logderiv(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s)


def v_direct(alpha, tau, dps=25):
    with mp.workdps(dps):
        return xi_logderiv(mpc(mpf(alpha), mpf(tau))).real


def b_alpha(alpha):
    a = mpf(alpha)
    if a == 1:
        return 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
    return 1 / a + 1 / (a - 1) - mp.log(mp.pi) / 2 + mp.digamma(a / 2) / 2 + mp.zeta(a, 1, 1) / mp.zeta(a)


def L_direct(tau, dps=30):
    # L(tau) = (xi'/xi)'(1/2 + i tau) = -1/s^2 - 1/(s-1)^2 + psi'(s/2)/4 + zeta''/zeta - (zeta'/zeta)^2
    with mp.workdps(dps):
        s = mpc(mpf(1) / 2, mpf(tau))
        z0, z1, z2 = mp.zeta(s), mp.zeta(s, 1, 1), mp.zeta(s, 1, 2)
        return -1 / s ** 2 - 1 / (s - 1) ** 2 + mp.psi(1, s / 2) / 4 + z2 / z0 - (z1 / z0) ** 2


def circle_coeffs(center, r, M, kmax):
    # Taylor coefficients a_0..a_kmax of xi'/xi at the centre from M equally spaced points of the circle of radius r
    r = mpf(r)
    pts = [(j, center + r * mp.expjpi(2 * (j + mpf(1) / 2) / M)) for j in range(M)]
    F = [xi_logderiv(z) for _, z in pts]
    return [sum(F[j] * mp.expjpi(-2 * k * (j + mpf(1) / 2) / M) for j in range(M)) / (M * r ** k) for k in range(kmax + 1)]


def load_zero_strings(path, n):
    with open(path) as fh:
        return [x.strip() for x in fh if x.strip()][:n]


def find_file(names, dirs):
    for d in dirs:
        for nme in names:
            p = os.path.join(d, nme)
            if os.path.exists(p):
                return p
    return None


def ball_mid(s):
    # the midpoint of a string '[x +/- r]' of the book's JSON files, or the string itself
    return mpf(s.split('[')[1].split('+/-')[0].strip()) if '[' in s else mpf(s)


# ---------------------------------------------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------------------------------------------

def stage_constants(args, model, model2, zs, gam, dirs):
    res = {}
    N = len(gam)
    g1, g2 = mpf(zs[0]), mpf(zs[1])
    res['zeros'] = {'file_count': N, 'first': zs[0], 'second': zs[1], 'last': zs[-1], 'gram_point_T0': fmt(model.T0, 12),
                    'N_of_T0_by_mpmath_nzeros': model.count_T0, 'S_of_T0_is_zero': bool(model.count_T0 == N),
                    'gram_point_T0_1200': fmt(model2.T0, 12), 'N_of_T0_1200': model2.count_T0, 'n_used_for_N_test': model2.n}
    gaps = np.diff(gam)
    order = np.argsort(gaps)[::-1]
    later = [i for i in order if i >= 1][0]
    res['gaps'] = {'gamma2_minus_gamma1': mpmath.nstr(g2 - g1, 15), 'largest_later_gap': fmt(gaps[later], 12),
                   'largest_later_gap_between_zero_numbers': [int(later + 1), int(later + 2)],
                   'largest_later_gap_ordinates': [zs[later][:16], zs[later + 1][:16]],
                   'max_distance_to_nearest_ordinate_on_13.5_tau_max': mpmath.nstr((g2 - g1) / 2, 15),
                   'five_largest_gaps': [fmt(gaps[i], 8) for i in order[:5]]}
    samp = [0, 1, 99] if args.quick else [0, 1, 2, 99, 999, N - 1]
    chk = {}
    with mp.workdps(50):
        for i in samp:
            chk['zero_number_%d' % (i + 1)] = fmt(mp.log10(abs(mp.zeta(mpc(mpf(1) / 2, mpf(zs[i]))))), 5)
    res['zero_file_sample_log10_abs_zeta_at_ordinate'] = chk
    log('zeros: N %d, Gram point T0 %.6f, nzeros(T0) %d, S(T0)=0 %s; N=%d test: T0 %.6f, nzeros %d' % (
        N, model.T0, model.count_T0, model.count_T0 == N, model2.n, model2.T0, model2.count_T0))
    # kappa_1, kappa_2, kappa_3 from the circles
    half = mpf(1) / 2
    kap = {}
    for r, M in ((5, 96), (8, 128)):
        a = circle_coeffs(mpc(half, 0), r, M, 5)
        kap[r] = a
    k1, k2, k3 = kap[5][1] / 2, -kap[5][3] / 2, kap[5][5] / 2
    k1b, k2b = kap[8][1] / 2, -kap[8][3] / 2
    ksum = math.fsum((1.0 / gam ** 2).tolist()) + model.kappa_tail()
    with mp.workdps(20):
        tq = mp.quad(lambda g: (mp.re(mp.digamma(mpf(1) / 4 + 0.5j * g)) - mp.log(mp.pi)) / (2 * mp.pi) / g ** 2, [model.T0, 4 * model.T0, 100 * model.T0, mp.inf])
    res['kappa'] = {'kappa1_circle_r5': mpmath.nstr(k1.real, 25), 'kappa1_circle_r8': mpmath.nstr(k1b.real, 25),
                    'kappa1_two_circles_differ_by': mpmath.nstr(abs(k1 - k1b), 3), 'kappa2_circle_r5': mpmath.nstr(k2.real, 20),
                    'kappa2_circle_r8': mpmath.nstr(k2b.real, 20), 'kappa3': mpmath.nstr(k3.real, 15),
                    'imag_part_of_kappa1_and_even_coefficients_a0_a2_a4': [mpmath.nstr(abs(k1.imag), 3)] + [mpmath.nstr(abs(kap[5][j]), 3) for j in (0, 2, 4)],
                    'kappa1_zero_sum_1700_plus_tail': fmt(ksum, 12), 'kappa1_zero_sum_minus_circle': fmt(ksum - float(k1.real), 3),
                    'tail_of_kappa1_gauss_legendre': fmt(model.kappa_tail(), 12), 'tail_of_kappa1_mpmath_quad': mpmath.nstr(tq, 12),
                    'two_kappa1': mpmath.nstr(2 * k1.real, 15), 'chapter_kappa1': '0.0231049931', 'chapter_two_kappa1': '0.046210'}
    log('kappa_1 circle %s, zero sum plus tail %s, tail GL %s quad %s' % (res['kappa']['kappa1_circle_r5'], res['kappa']['kappa1_zero_sum_1700_plus_tail'],
                                                                         res['kappa']['tail_of_kappa1_gauss_legendre'], res['kappa']['tail_of_kappa1_mpmath_quad']))
    # b_alpha, the floor row and Theorem yo:thm:coercive (iv)
    b1c = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
    bal = {}
    for a in ALPHAS_ALL:
        bb = b_alpha(a)
        cc = mpf(str(cval(a).numerator)) / mpf(str(cval(a).denominator))
        lo = 4 * k1.real * cc - 4 * k2.real * cc ** 3
        up = 4 * k1.real * cc
        bal[a] = {'c': str(cval(a)), 'b_alpha': mpmath.nstr(bb, 15), 'two_b_alpha': mpmath.nstr(2 * bb, 12), 'floor_row_chapter': CHAPTER_FLOOR[a],
                  'floor_row_agrees_to_printed_digits': bool(('%.7f' % float(2 * bb)) == CHAPTER_FLOOR[a]),
                  'iv_lower_4k1c_minus_4k2c3': mpmath.nstr(lo, 12), 'iv_upper_4k1c': mpmath.nstr(up, 12),
                  'iv_holds': bool(lo <= 2 * bb <= up), 'b_over_c_to_two_kappa1': mpmath.nstr(bb / cc / (2 * k1.real), 8)}
    res['b_alpha'] = bal
    res['xi_over_xi_at_1'] = {'closed_form_1_plus_gammaE_over_2_minus_half_log_4pi': mpmath.nstr(b1c, 25), 'chapter': '0.0230957',
                              'from_formula_at_alpha_1_plus_1e-12': mpmath.nstr(b_alpha(mpf(1) + mpf(10) ** -12), 15)}
    log('b_alpha: ' + ', '.join('%s %s' % (a, bal[a]['b_alpha']) for a in ALPHAS_ALL) + '; floor row ' + ', '.join('%s' % bal[a]['two_b_alpha'] for a in ALPHAS_ALL))
    # the margins of the proof of Theorem yo:thm:coercive (ii) (line 278)
    c = mpf(1) / 2
    t2 = 2 / g1 ** 2
    t3 = 2 * k1.real
    mar = {'gamma1': mpmath.nstr(g1, 15), 'two_over_gamma1_sq': mpmath.nstr(t2, 12), 'two_kappa1': mpmath.nstr(t3, 12),
           'chapter': {'two_over_gamma1_sq': '0.010010', 'two_kappa1': '0.046210', 'inverse_4.15': '0.057233', 'inverse_3.4437': '0.082583',
                       'margin_4.15': '0.001013', 'margin_3.4437': '0.026362'}}
    for dl in ('4.15', '3.4437'):
        d = mpf(dl)
        t1 = 1 / (c * c + d * d)
        mar['d_' + dl] = {'one_over_c2_plus_d2': mpmath.nstr(t1, 12), 'margin_chapter_formula': mpmath.nstr(t1 - t2 - t3, 10),
                          'margin_chapter_formula_from_printed_6_digits': str(round(float(mpmath.nstr(t1, 6)) - 0.010010 - 0.046210, 6)),
                          'margin_without_the_term_2_over_gamma1_sq': mpmath.nstr(t1 - t3, 10),
                          'largest_d_with_positive_margin_chapter_formula': mpmath.nstr(mp.sqrt(1 / (t2 + t3) - c * c), 8),
                          'largest_d_without_the_term': mpmath.nstr(mp.sqrt(1 / t3 - c * c), 8)}
        mar['d_' + dl]['margin_by_c_chapter_formula'] = {str(cc): mpmath.nstr(1 / (mpf(cc) ** 2 + d * d) - t2 - t3, 8) for cc in (0, 0.05, 0.1, 0.25, 0.5)}
    res['margins'] = mar
    rng = np.random.default_rng(20261017)
    cs, gs, ts = rng.uniform(1e-3, 0.5, 200000), rng.uniform(14.0, 120.0, 200000), rng.uniform(0.0, 400.0, 200000)
    Pk = lambda t: cs / (cs ** 2 + (t - gs) ** 2) + cs / (cs ** 2 + (t + gs) ** 2)
    dif = Pk(ts) - Pk(0.0)
    okp = (dif >= 0) == (ts ** 2 <= 3 * gs ** 2 - cs ** 2)
    near = np.abs(dif) < 1e-12 * Pk(0.0)
    res['lemma_P_gamma_tau_ge_P_gamma_0_iff_tau2_le_3gamma2_minus_c2'] = {'random_triples': 200000, 'mismatches': int(np.sum(~okp & ~near)), 'near_boundary_ignored': int(np.sum(~okp & near))}
    log('margins at c=1/2: d=4.15 %s, d=3.4437 %s (chapter formula); without the term 2/gamma_1^2: %s, %s' % (
        mar['d_4.15']['margin_chapter_formula'], mar['d_3.4437']['margin_chapter_formula'],
        mar['d_4.15']['margin_without_the_term_2_over_gamma1_sq'], mar['d_3.4437']['margin_without_the_term_2_over_gamma1_sq']))
    # the book's own certified data
    bk = {}
    p = find_file(['v1min.json'], [os.path.join(args.root, 'riesz_program', 'data')] + dirs)
    if p:
        d = json.load(open(p))
        bk['v1min_json'] = os.path.relpath(p, args.root)
        bk['b1_book'] = d['b1']
        bk['b1_difference'] = mpmath.nstr(b1c - ball_mid(d['b1']), 3)
        bk['v1min_zeros_certified'] = d['zeros']['zeros_used']
        bk['v1min_largest_certified_ordinate'] = d['zeros']['largest_ordinate']
        bk['gamma2_minus_gamma1_book'] = d['zeros']['max_gap_gamma2_minus_gamma1']
        bk['largest_later_gap_book'] = d['zeros']['largest_later_gap']
        bk['max_distance_to_ordinate_book'] = d['zeros']['max_distance_to_nearest_zero_on_13.5_2000']
        bk['profile_min_13.5_2000_book'] = d['profile_min_on_13.5_2000']
        bk['v1_at_17.6_book'] = d['v1_at']['17.6']
        bk['v1_at_17.6_direct_here'] = mpmath.nstr(v_direct('1', '17.6', 30), 20)
        bk['v1_at_17.5_direct_here'] = mpmath.nstr(v_direct('1', '17.5', 30), 20)
        bk['v1_at_17.6_difference'] = mpmath.nstr(v_direct('1', '17.6', 30) - ball_mid(d['v1_at']['17.6']), 3)
        bk['v1_at_17.5_minus_book_profile_min'] = mpmath.nstr(v_direct('1', '17.5', 30) - ball_mid(d['profile_min_on_13.5_2000']['value']), 3)
    p = find_file(['dbeta_minimum.json'], [os.path.join(args.root, 'data')] + dirs)
    if p:
        d = json.load(open(p))
        bk['dbeta_minimum_json'] = os.path.relpath(p, args.root)
        bk['V_book_log_xi_second_derivative_at_half'] = d['xi']['V']
        bk['V_difference_from_two_kappa1'] = mpmath.nstr(2 * k1.real - mpf(d['xi']['V']), 3)
        bk['b_alpha_book'] = {a: d['xi']['b_alpha'][a]['b'] for a in ('0.55', '0.6', '0.75', '1')}
        bk['H0_verified_height_for_xi'] = d['xi']['H0']
        bk['alpha_max_for_xi'] = d['delta_max']['xi']['sigma_max = 1 + delta_max']
    res['book_crosschecks'] = bk if bk else {'available': False}
    res['kappa1_value'] = float(k1.real)
    return res


def stage_symbol(args, model, model2, gam, alphas, out):
    sym = out.setdefault('symbol', {})
    fig = out.setdefault('figure', {'left': {'step': '0.01', 'tau_max': '60', 'v': {}}, 'right': {'window': '10', 'lo': '20', 'min': {}}, 'b': {}, 'global_min': {}})
    tmax = args.tau_max
    for a in alphas:
        c = float(cval(a))
        h = min(c / args.divisor, 0.05)
        K = int(round(tmax / h))
        tau = np.arange(K + 1) * h
        bmp = b_alpha(a)
        b = float(bmp)
        t0 = time.time()
        v = model.total('v', tau, c)
        log('alpha %s: c %s step %s, %d grid points, scan %.1f s' % (a, c, h, K + 1, time.time() - t0))
        k20 = int(round(20.0 / h))
        f = lambda t: float(model.total('v', [t], c)[0])

        def refine(i):
            try:
                r = minimize_scalar(f, bracket=(tau[i - 1], tau[i], tau[i + 1]), method='brent', tol=1e-13)
                x, y = float(r.x), float(r.fun)
            except Exception:
                r = minimize_scalar(f, bounds=(tau[i - 1], tau[i + 1]), method='bounded', options={'xatol': 1e-12})
                x, y = float(r.x), float(r.fun)
            if not (tau[i - 1] <= x <= tau[i + 1]) or y > v[i] + 1e-14:
                return float(tau[i]), float(v[i])
            return x, y

        interior = np.where((v[1:-1] <= v[:-2]) & (v[1:-1] < v[2:]))[0] + 1
        above = interior[interior >= k20]
        cand = above[np.argsort(v[above])[:30]]
        ref = sorted((refine(int(i)) for i in cand), key=lambda p: p[1])
        boundary = bool(v[k20] <= v[k20 + 1])
        if boundary:
            ref.append((float(tau[k20]), float(v[k20])))
            ref.sort(key=lambda p: p[1])
        tstar, vstar = ref[0]
        n_idx = int(np.searchsorted(gam, tstar))
        low_int = interior[np.argsort(v[interior])[:5]]
        li = sorted((refine(int(i)) for i in low_int), key=lambda p: p[1])[0]
        # direct checks (tail free) and the 1200 zero version at the minimum and at points of the range
        pts = [tstar, li[0]] + ([100.3, tmax] if args.quick else [100.3, 777.7, 1500.5, tmax])
        checks = []
        for t in pts:
            dv = v_direct(a, t)
            sv = float(model.total('v', [t], c)[0])
            s2 = float(model2.total('v', [t], c)[0]) if t < model2.T0 - 150.0 else None
            checks.append({'tau': fmt(t, 12), 'direct': mpmath.nstr(dv, 15), 'zero_sum_plus_tail': fmt(sv, 14), 'difference': fmt(sv - float(dv), 3),
                           'sum_with_1200_zeros_minus_1700': (fmt(s2 - sv, 3) if s2 is not None else 'above the cut of the 1200 zero tail')})
        v0 = float(model.total('v', [0.0], c)[0])
        beyond1 = tau >= 1.0
        w = int(round(10.0 / h))
        wmin = np.minimum.reduceat(v[k20:], np.arange(0, len(v) - k20, w))
        tl = np.arange(6001) * 0.01
        vl = model.total('v', tl, c)
        e = {'c': str(cval(a)), 'b_alpha': mpmath.nstr(bmp, 15), 'grid_step': fmt(h), 'grid_points': K + 1, 'tau_max': fmt(tmax),
             'v_at_0_zero_sum_plus_tail': fmt(v0, 14), 'v_at_0_minus_b_alpha': fmt(v0 - b, 3),
             'min_on_[20,tau_max]': {'value': fmt(vstar, 12), 'tau': fmt(tstar, 12), 'ratio_to_b_alpha': fmt(vstar / b, 8),
                                    'attained_at_boundary_20': bool(abs(tstar - 20.0) < 1e-9), 'between_zero_numbers': [n_idx, n_idx + 1],
                                    'between_ordinates': [fmt(gam[n_idx - 1], 12), fmt(gam[n_idx], 12)],
                                    'chapter': CHAPTER_MIN[a], 'direct_value_at_the_minimum': mpmath.nstr(v_direct(a, tstar), 15)},
             'lowest_5_refined_minima_above_20': [{'tau': fmt(t, 10), 'value': fmt(y, 10), 'ratio_to_b_alpha': fmt(y / b, 6)} for t, y in ref[:5]],
             'lowest_interior_local_minimum': {'value': fmt(li[1], 12), 'tau': fmt(li[0], 12), 'ratio_to_b_alpha': fmt(li[1] / b, 8)},
             'v_exceeds_b_alpha_at_every_grid_point_with_tau_ge_1': bool(np.all(v[beyond1] > b)),
             'smallest_v_minus_b_over_b_for_tau_ge_1': fmt(float(np.min((v[beyond1] - b) / b)), 6),
             'grid_local_minima_above_20': int(len(above)), 'direct_checks': checks,
             'largest_abs_difference_sum_minus_direct': fmt(max(abs(float(x['difference'])) for x in checks), 3)}
        sym[a] = e
        fig['left']['v'][a] = ['%.6g' % x for x in vl]
        fig['right']['min'][a] = ['%.6g' % x for x in wmin]
        fig['b'][a] = mpmath.nstr(bmp, 12)
        fig['global_min'][a] = {'tau': fmt(tstar, 10), 'value': fmt(vstar, 10)}
        log('alpha %s: b %s; min over [20, %g] = %s at tau %s (ratio %s, chapter %s), gap of zeros %d, %d; lowest interior %s at %s (ratio %s); max |sum - direct| %s' % (
            a, e['b_alpha'], tmax, e['min_on_[20,tau_max]']['value'], e['min_on_[20,tau_max]']['tau'], e['min_on_[20,tau_max]']['ratio_to_b_alpha'], CHAPTER_MIN[a],
            n_idx, n_idx + 1, e['lowest_interior_local_minimum']['value'], e['lowest_interior_local_minimum']['tau'], e['lowest_interior_local_minimum']['ratio_to_b_alpha'],
            e['largest_abs_difference_sum_minus_direct']))


def stage_laguerre(args, model, model2, gam, kappa1, out):
    lmax = args.l_max
    L0 = float(model.total('L', [0.0])[0])
    res = {'chapter': {'L0': '0.046210', 'next_local_minimum': '0.2301', 'at_tau': '17.50', 'ratio': '4.98'}, 'range': [0, lmax]}
    res['L0_zero_sum_plus_tail'] = fmt(L0, 12)
    res['two_kappa1'] = fmt(2 * kappa1, 12)
    res['L0_minus_two_kappa1'] = fmt(L0 - 2 * kappa1, 3)
    taus, Ls, ids = [], [], []
    fd = lambda t: float(model.total('dL', [t])[0])
    t0 = time.time()
    for n in range(len(gam) - 1):
        lo, hi = float(gam[n]), float(gam[n + 1])
        if lo > lmax:
            break
        eps = 1e-3 * (hi - lo)
        try:
            t = brentq(fd, lo + eps, hi - eps, xtol=1e-12, rtol=1e-14, maxiter=300)
        except ValueError:
            t = float(minimize_scalar(lambda x: float(model.total('L', [x])[0]), bounds=(lo + eps, hi - eps), method='bounded', options={'xatol': 1e-12}).x)
        taus.append(t)
        Ls.append(float(model.total('L', [t])[0]))
        ids.append(n + 1)
    taus, Ls, ids = np.array(taus), np.array(Ls), np.array(ids)
    order = np.argsort(Ls)
    res['gaps_scanned'] = int(len(Ls))
    res['minimum_of_L_in_each_gap_exceeds_L0'] = bool(Ls.min() > L0)
    res['global_minimum_on_[0,lmax]_is_at_0'] = bool(Ls.min() > L0)
    res['lowest_gap_minima'] = [{'between_zero_numbers': [int(ids[i]), int(ids[i]) + 1], 'tau': fmt(taus[i], 12), 'L': fmt(Ls[i], 12), 'ratio_to_L0': fmt(Ls[i] / L0, 8)} for i in order[:6]]
    dchk = []
    for i in order[:3]:
        z = L_direct(taus[i])
        dchk.append({'tau': fmt(taus[i], 12), 'direct_real': mpmath.nstr(z.real, 15), 'direct_imag': mpmath.nstr(z.imag, 3), 'zero_sum_plus_tail': fmt(Ls[i], 14),
                     'difference': fmt(Ls[i] - float(z.real), 3), 'sum_with_1200_zeros': fmt(float(model2.total('L', [taus[i]])[0]) - Ls[i], 3)})
    z0 = L_direct(0.0)
    dchk.append({'tau': '0', 'direct_real': mpmath.nstr(z0.real, 15), 'direct_imag': mpmath.nstr(z0.imag, 3), 'zero_sum_plus_tail': fmt(L0, 14), 'difference': fmt(L0 - float(z0.real), 3)})
    res['direct_checks'] = dchk
    # a dense grid as a check of the convexity argument
    step = 0.01
    lg = np.arange(int(round(args.l_grid_max / step)) + 1) * step
    Lg = model.total('L', lg)
    sel = lg >= float(gam[0]) + 0.5
    res['grid_check'] = {'step': fmt(step), 'range': [0, args.l_grid_max], 'grid_min_over_all': fmt(Lg.min(), 12), 'grid_argmin': fmt(lg[int(np.argmin(Lg))], 6),
                         'grid_min_beyond_gamma_1_plus_0.5': fmt(Lg[sel].min(), 12), 'at': fmt(lg[sel][int(np.argmin(Lg[sel]))], 6),
                         'grid_min_beyond_gamma_1_minus_lowest_gap_minimum': fmt(Lg[sel].min() - Ls.min(), 3)}
    t5 = brentq(lambda t: float(model.total('L', [t])[0]) - Ls.min(), 1.0, float(gam[0]) - 0.01, xtol=1e-10)
    res['L_reaches_the_next_local_minimum_at_tau'] = fmt(t5, 8)
    res['L_increases_on_[0,gamma_1)_checked_on_grid'] = bool(np.all(np.diff(Lg[lg < gam[0] - 0.01]) > 0))
    out['laguerre'] = res
    log('L(0) = %s (2 kappa_1 = %s); lowest gap minima: %s; global minimum at 0: %s; scan %.1f s' % (
        res['L0_zero_sum_plus_tail'], res['two_kappa1'], ', '.join('%s at %s (zeros %d, %d)' % (x['L'], x['tau'], x['between_zero_numbers'][0], x['between_zero_numbers'][1]) for x in res['lowest_gap_minima'][:3]),
        res['global_minimum_on_[0,lmax]_is_at_0'], time.time() - t0))
    log('ratio of the next local minimum to L(0): %s; direct checks of L: %s' % (res['lowest_gap_minima'][0]['ratio_to_L0'], ', '.join(x['difference'] for x in dchk)))


def make_figure(out, dst):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})
    fg = out.get('figure')
    if not fg or not fg['left']['v']:
        log('no figure data, figure skipped')
        return
    style = {'1': ('black', '-', 1.0), '0.75': ('black', (0, (5, 2)), 0.9), '0.6': ('0.35', (0, (6, 2, 1, 2)), 1.0), '0.55': ('0.15', (0, (1, 1.6)), 1.3)}
    alphas = [a for a in ALPHAS_ALL if a in fg['left']['v']]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.4, 2.7))
    tl = np.arange(len(fg['left']['v'][alphas[0]])) * float(fg['left']['step'])
    lo, wd = float(fg['right']['lo']), float(fg['right']['window'])
    for a in alphas:
        col, ls, lw = style[a]
        ax1.plot(tl, [float(x) for x in fg['left']['v'][a]], color=col, ls=ls, lw=lw, label=r'$\alpha=%s$' % a)
        wm = [float(x) for x in fg['right']['min'][a]]
        xc = lo + wd * (np.arange(len(wm)) + 0.5)
        ax2.plot(xc, wm, color=col, ls=ls, lw=lw, label=r'$\alpha=%s$' % a)
        b = float(fg['b'][a])
        for ax in (ax1, ax2):
            ax.axhline(b, color='0.78', lw=0.5, zorder=0)
        gm = fg['global_min'][a]
        ax2.plot([float(gm['tau'])], [float(gm['value'])], 'o', mfc='white', mec='black', ms=3.2, mew=0.7, zorder=5)
        ax1.text(60.8, b, r'$b_{%s}$' % a, fontsize=5.5, va='center', color='0.3', clip_on=False)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_xlim(0, 60)
    ax1.set_ylim(1.5e-3, 40)
    ax2.set_xlim(0, lo + wd * len(fg['right']['min'][alphas[0]]))
    ax2.set_ylim(1.5e-3, 12)
    ax1.set_xlabel(r'$\tau$')
    ax2.set_xlabel(r'$\tau$')
    ax1.set_ylabel(r'$\mathrm{Re}\,\xi^{\prime}/\xi(\alpha+i\tau)$')
    ax2.set_ylabel('smallest value in a window of length 10')
    ax2.legend(loc='upper left', frameon=False, fontsize=6.5, ncol=2, handlelength=2.6)
    ax2.set_title(r'circled: the minimum above $\tau=20$', fontsize=6.5)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.32)
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    fig.savefig(dst)
    plt.close(fig)
    log('wrote %s' % dst)


def agree(chapter, value):
    if chapter is None:
        return 'no value in the chapter'
    dec = len(chapter.split('.')[1]) if '.' in chapter else 0
    return 'agrees to the printed digits' if '%.*f' % (dec, float(value)) == chapter else 'differs'


def build_summary(out):
    S = []

    def add(claim, line, chapter, value):
        S.append({'claim': claim, 'chapter_line': line, 'chapter_value': chapter, 'computed': value, 'status': agree(chapter, value)})

    c, sym, lag = out.get('constants', {}), out.get('symbol', {}), out.get('laguerre', {})
    for a in ALPHAS_ALL:
        if a in sym:
            m = sym[a]['min_on_[20,tau_max]']
            add('minimum of v_alpha on [20, tau_max], alpha = ' + a, '329, 394', CHAPTER_MIN[a], m['value'])
            add('that minimum over b_alpha ("about eight"), alpha = ' + a, '329, 394', '8', m['ratio_to_b_alpha'])
            add('tau of that minimum, alpha = ' + a, '(not in the chapter)', None, m['tau'])
            add('lowest interior local minimum over b_alpha (minimum on (0, tau_max] is at 0 iff above 1), alpha = ' + a, '329', None, sym[a]['lowest_interior_local_minimum']['ratio_to_b_alpha'])
    if lag:
        g = lag['lowest_gap_minima'][0]
        add('L(0) = 2 kappa_1', '298', '0.046210', lag['two_kappa1'])
        add('next local minimum of L', '298', '0.2301', g['L'])
        add('tau of the next local minimum of L', '298', '17.50', g['tau'])
        add('next local minimum over L(0)', '298, 395', '4.98', g['ratio_to_L0'])
    if c:
        m = c['margins']
        add('2/gamma_1^2', '278', '0.010010', m['two_over_gamma1_sq'])
        add('2 kappa_1', '278', '0.046210', m['two_kappa1'])
        add('1/(c^2 + 4.15^2) at c = 1/2', '278', '0.057233', m['d_4.15']['one_over_c2_plus_d2'])
        add('1/(c^2 + 3.4437^2) at c = 1/2', '(review M10)', '0.082583', m['d_3.4437']['one_over_c2_plus_d2'])
        add('margin at d = 4.15', '278', '0.001013', m['d_4.15']['margin_chapter_formula'])
        add('margin at d = 3.4437', '(review M10)', '0.026362', m['d_3.4437']['margin_chapter_formula'])
        add('kappa_1', '(task)', '0.0231049931', c['kappa']['kappa1_circle_r5'])
        add('xi\'/xi(1)', '(task)', '0.0230957', c['xi_over_xi_at_1']['closed_form_1_plus_gammaE_over_2_minus_half_log_4pi'])
        for a in ('0.55', '0.6', '0.75', '1'):
            add('floor row 2 b_alpha, alpha = ' + a, '312', CHAPTER_FLOOR[a], c['b_alpha'][a]['two_b_alpha'])
    return S


def dumps_compact(obj):
    text = json.dumps(obj, indent=1)
    return re.sub(r'\[\s+([^\[\]{}]*?)\s+\]', lambda m: '[' + re.sub(r'\s+', ' ', m.group(1)) + ']', text)


BOOK_RESULTS = {
    'alpha_1': 'Proposition bp:prop:v1min (ch/basepoint.tex): v_1(theta) > v_1(0) = b_1 for every real theta != 0 with no hypothesis, from the first 1800 zeros certified to height 2304.3 by riesz_program/code/v1min.py (riesz_program/data/v1min.json) and the bounds of HSW22 above height 2000; Proposition bp:prop:valpha (certificate repro/valpha_min.py, repro/valpha_min.json) for every alpha >= 1; so the whole plotted range [0, 2000] (and every height) is covered.',
    'alpha_0.75_0.6_0.55': 'Theorem thm:dbeta:minall (ii) (ch/basepoint.tex, constants of Table tab:dbeta:minall from code/dbeta_minimum.py, data/dbeta_minimum.json): for 1/2 < alpha <= alpha_max = 7.35, v_alpha(theta) > b_alpha for 0 < |theta| <= H_0 - (1 + sqrt 2)/2 with H_0 = 3e12 (zeros verified on the line to that height), no hypothesis; part (vi) says that the same for every theta and every alpha in (1/2, 1) is equivalent to RH; so the plotted range [0, 2000] is covered for all three.',
    'only_numerics': 'On the plotted range nothing of the statement v_alpha > b_alpha (tau != 0) is only numerics for the four alpha; the numerics are the values and places of the minima above tau = 20 (the grid scan, refined) and the ratios to b_alpha, which no theorem of the book states. Beyond height 3e12 for alpha < 1 the statement is equivalent to RH.',
}


def main():
    global LOGFILE
    ap = argparse.ArgumentParser(description='Group Y4 checks of Chapter 17: the symbol, the Laguerre quantity, kappa_1, b_alpha and the margins.')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--fig')
    ap.add_argument('--root')
    ap.add_argument('--zeros')
    ap.add_argument('--stages', default='constants,symbol,laguerre,figure')
    ap.add_argument('--alphas', default=','.join(ALPHAS_ALL))
    ap.add_argument('--tau-max', type=float)
    ap.add_argument('--divisor', type=int)
    ap.add_argument('--n-zeros', type=int, default=1700)
    ap.add_argument('--update', action='store_true')
    args = ap.parse_args()
    parent = os.path.dirname(HERE)
    args.root = args.root or (parent if os.path.exists(os.path.join(parent, 'data', 'zeros_1700.txt')) else HERE)
    dirs = [HERE, os.path.join(HERE, 'data')]
    args.tau_max = args.tau_max or (300.0 if args.quick else 2000.0)
    args.divisor = args.divisor or (4 if args.quick else 10)
    args.l_max = args.l_grid_max = 300.0 if args.quick else 1800.0
    out_path = args.out or ('/tmp/yo_symbol_quick.json' if args.quick else os.path.join(args.root, 'data', 'yo_symbol.json'))
    fig_path = args.fig or ('/tmp/yo_symbol_quick.pdf' if args.quick else os.path.join(args.root, 'fig', 'yo_symbol.pdf'))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    LOGFILE = open(os.path.splitext(out_path)[0] + '.log', 'a' if args.update else 'w')
    zpath = args.zeros or find_file(['zeros_1700.txt'], [os.path.join(args.root, 'data')] + dirs)
    stages = args.stages.split(',')
    alphas = args.alphas.split(',')
    log('yo_symbol.py %s on %s, python %s, numpy %s, mpmath %s; stages %s, alphas %s, tau_max %g, divisor %d' % (
        'quick' if args.quick else 'full', socket.gethostname(), platform.python_version(), np.__version__, mpmath.__version__, args.stages, args.alphas, args.tau_max, args.divisor))
    mp.dps = 40
    out = json.load(open(out_path)) if (args.update and os.path.exists(out_path)) else {}
    out['script'] = 'code/yo_symbol.py'
    out.setdefault('runs', []).append({'mode': 'quick' if args.quick else 'full', 'host': socket.gethostname(), 'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'python': platform.python_version(),
                  'numpy': np.__version__, 'mpmath': mpmath.__version__, 'zeros_file': os.path.basename(zpath), 'zeros_used': args.n_zeros, 'tau_max': args.tau_max,
                  'grid_divisor': args.divisor, 'stages_of_this_run': stages, 'precision': 'float64 sums, mpmath 25 to 40 digits for direct values'})
    zs = load_zero_strings(zpath, args.n_zeros)
    gam = np.array([float(x) for x in zs])
    model = Model(gam, len(gam))
    model2 = Model(gam, min(1200, len(gam) - 1))
    out['book_results'] = BOOK_RESULTS
    kappa1 = None
    if 'constants' in stages:
        out['constants'] = stage_constants(args, model, model2, zs, gam, dirs)
        kappa1 = out['constants'].pop('kappa1_value')
    if 'symbol' in stages:
        stage_symbol(args, model, model2, gam, alphas, out)
    if 'laguerre' in stages:
        if kappa1 is None:
            mp.dps = 30
            kappa1 = float(circle_coeffs(mpc(mpf(1) / 2, 0), 5, 96, 1)[1].real / 2)
        stage_laguerre(args, model, model2, gam, kappa1, out)
    out['summary'] = build_summary(out)
    for e in out['summary']:
        log('%-95s chapter %-12s computed %-16s %s' % (e['claim'], e['chapter_value'], e['computed'][:14], e['status']))
    with open(out_path, 'w') as fh:
        fh.write(dumps_compact(out) + '\n')
    log('wrote %s' % out_path)
    if 'figure' in stages:
        make_figure(out, fig_path)
    log('done')


if __name__ == '__main__':
    main()
