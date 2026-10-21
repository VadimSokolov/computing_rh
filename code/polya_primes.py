# polya_primes.py: the primes in the compensated clock, Proposition pol:prop:primes, Lemma pol:lem:Bbound and
# the numerical sentences of Theorem pol:thm:almostall, and the direct computations of Section pol:sec:qrh,
# in Chapter 16 (Polya and the zeros of entire functions). Added for the book in October 2026, written from the mathematics (the chapter
# author's own scripts were not read, run or copied).
# What it checks. (1) The identity Re zeta'/zeta(alpha + i r) - zeta'/zeta(alpha) = sum_{n >= 2} Lambda(n) n^-alpha (1 - cos(r log n)): the left side
# by mpmath at 30 digits, the right side by a segmented sieve over n <= N (N = 1e10, eight processes, double precision with fsum over segments)
# at (alpha, r) = (3, 7), (3, 12), (1.5, 7), (1.5, 20). Rigorous tail: with psi(x) < 1.03883 x for x > 0 (Rosser and Schoenfeld 1962) and partial
# summation, sum_{n > N} Lambda(n) n^-alpha <= 1.03883 alpha N^(1 - alpha)/(alpha - 1) - psi(N) N^-alpha, psi(N) computed by the sieve, and the weight
# 1 - cos is at most 2. The tail is also completed by the prime number theorem main term int_N^oo x^-alpha (1 - cos(r log x)) dx in closed form (a
# check, not a bound). (2) alpha_1: the archimedean bracket Br(alpha, r) = Re A(alpha + i r) - A(alpha), A = 1/s + 1/(s-1) - (log pi)/2 + (1/2) psi(s/2),
# its minimum over r in (0, 200] by a grid (logarithmic from 1e-7 and linear step 0.1, about 4400 points) refined by bounded Brent minimisation at 40
# digits, and the threshold by bisection on the sign of that minimum (negative below minus 1e-25); also the root of A''(alpha) = 0 (the local
# form at r = 0, which is where the last dip closes), a scan of the minimum over alpha in [1.05, 12], positivity above the threshold up to r = 1e10, and the minima at alpha = 1.2, 2, 6.
# (3) Lemma Bbound: a_1 = sum 1/gamma^2 = (1/2)(xi'/xi)'(1/2) (mpmath diff at 40 digits, and the sum over data/zeros_1700.txt plus the smooth tail),
# B_alpha/(alpha - 1/2) at alpha = 0.55, 0.75, 1 (B_1 = 1 + gamma/2 - (log pi)/2 - log 2), 2 a_1, 4 a_1, and the arithmetic of the proof of Theorem
# almostall (kappa >= (4/5) c_min, 4 a_1 against the printed 0.0924, S = sum_{gamma > H} gamma^-2 against 1.6e-12). (4) D_alpha(r) = Re xi'/xi(alpha + i r)
# - xi'/xi(alpha) on 20000 heights r = 0.01 i in (0, 200] at fifteen equally spaced base points in [1, 7.8] and thirteen more; the minimum, where it is
# attained, the number of nonpositive values, D/r^2 at the first height, the smallest local minima, the distance of the argmin to the nearest ordinate.
# (5) The runs on (0, 600] (60000 heights, step 0.01) at alpha = 0.9, 0.95, 1, 1.001, 2, 4, 7.8: positivity, the minimum over (2, 600] with its location
# (local minima refined by Brent at 20 digits), D at r = 2 and 2.01.
# Error control. The sieve sum is compared with the zeta side at 12 or more digits where the rigorous tail allows it, the bracket minima are found
# on a grid and then refined (not by a coarse grid alone), a_1 is computed two ways, and each D value is a difference of 20 digit numbers. Nothing here
# is a proof: positivity on grids is not positivity for all r.
# Output: data/polya_primes.json (numbers as strings), the log is the standard output (data/polya_primes.log).
# Usage: python3 polya_primes.py [--quick] [--out FILE] [--workers W] [--zeros FILE]; --quick runs a small case in one process (about 25 s, writes to /tmp
# unless --out is given); the full case is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8 -g polya_primes.json code/polya_primes.py data/zeros_1700.txt -- --out polya_primes.json
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar

T0 = time.time()
RS = 1.03883   # Rosser and Schoenfeld (1962): psi(x) < 1.03883 x for all x > 0
H_VER = 3e12


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def default_out(name):
    return os.path.join('data', name) if os.path.isdir('data') else name


def find_zeros(path):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in ([path] if path else []) + ['data/zeros_1700.txt', 'zeros_1700.txt', os.path.join(here, 'zeros_1700.txt'), os.path.join(here, '..', 'data', 'zeros_1700.txt')]:
        if os.path.exists(c):
            return c
    raise FileNotFoundError('zeros_1700.txt')


# ---------------------------------------------------------------- mathematical ingredients
def A_arch(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2


def xi_logder(s):
    return A_arch(s) + mp.zeta(s, derivative=1) / mp.zeta(s)


def B_val(a):
    a = mp.mpf(a)
    if a == 1:
        return 1 + mp.euler / 2 - mp.log(mp.pi) / 2 - mp.log(2)
    return mp.re(xi_logder(a))


def zl(s):
    return mp.zeta(s, derivative=1) / mp.zeta(s)


def bracket(a, r):
    s = mp.mpc(a, r)
    return mp.re(1 / s + 1 / (s - 1) + mp.digamma(s / 2) / 2) - (1 / a + 1 / (a - 1) + mp.digamma(a / 2) / 2)


def make_grid(quick):
    if quick:
        u = np.linspace(-7, math.log10(200), 400)
        lin = np.arange(0.2, 200.001, 0.5)
    else:
        u = np.linspace(-7, math.log10(200), 2400)
        lin = np.arange(0.1, 200.0001, 0.1)
    return [float(x) for x in np.unique(np.concatenate([10 ** u, lin]))]


def bracket_min(alpha, grid):
    mp.mp.dps = 40
    a = mp.mpf(alpha)
    vals = [bracket(a, mp.mpf(r)) for r in grid]
    i = min(range(len(vals)), key=lambda j: vals[j])
    best, rbest = vals[i], grid[i]
    if 0 < i < len(grid) - 1:
        res = minimize_scalar(lambda x: float(bracket(a, mp.mpf(x))), bounds=(grid[i - 1], grid[i + 1]), method='bounded', options={'xatol': 1e-13})
        if res.fun < float(best):
            best, rbest = mp.mpf(res.fun), float(res.x)
    return best, rbest


# ---------------------------------------------------------------- sieve
_SP = {}


def small_primes(n):
    s = np.ones(n + 1, dtype=np.bool_)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def seg_task(args):
    lo, hi, nsmall, pairs = args
    if nsmall not in _SP:
        _SP[nsmall] = small_primes(nsmall)
    seg = np.ones(hi - lo, dtype=np.bool_)
    for p in _SP[nsmall]:
        p = int(p)
        if p * p >= hi:
            break
        seg[max(p * p, ((lo + p - 1) // p) * p) - lo::p] = False
    if lo < 2:
        seg[:2 - lo] = False
    lp = np.log((np.nonzero(seg)[0] + lo).astype(np.float64))
    return [float(np.sum(lp * np.exp(-a * lp) * (1.0 - np.cos(r * lp)))) for a, r in pairs], float(np.sum(lp))


def prime_sums(N, pairs, pmap, seg=1 << 24):
    nsmall = int(math.isqrt(N)) + 1
    res = pmap(seg_task, [(lo, min(lo + seg, N + 1), nsmall, pairs) for lo in range(0, N + 1, seg)])
    sums = [math.fsum(r[0][i] for r in res) for i in range(len(pairs))]
    theta = math.fsum(r[1] for r in res)
    psi_extra = 0.0
    for p in small_primes(int(math.isqrt(N))):
        p = int(p)
        lp = math.log(p)
        k, pk = 2, p * p
        while pk <= N:
            for i, (a, r) in enumerate(pairs):
                sums[i] += lp * math.exp(-a * k * lp) * (1.0 - math.cos(r * k * lp))
            psi_extra += lp
            k += 1
            pk *= p
    return sums, theta + psi_extra


# ---------------------------------------------------------------- tasks
def t_scan(args):
    alpha, grid = args
    best, rb = bracket_min(alpha, grid)
    return alpha, mp.nstr(best, 12), '%.9g' % rb


def t_D(args):
    alpha, step, idx, dps = args
    mp.mp.dps = dps
    a = mp.mpf(alpha)
    B = B_val(a)
    return [float(mp.re(xi_logder(mp.mpc(a, i * mp.mpf(step)))) - B) for i in idx]


def D_scalar(alpha, r, dps=20):
    mp.mp.dps = dps
    a = mp.mpf(alpha)
    return float(mp.re(xi_logder(mp.mpc(a, mp.mpf(r)))) - B_val(a))


def analyse(rs, D, zeros):
    n = len(D)
    imin = min(range(n), key=lambda i: D[i])
    lm = sorted((D[i], rs[i]) for i in range(1, n - 1) if D[i] < D[i - 1] and D[i] <= D[i + 1])
    far = [(D[i], rs[i]) for i in range(n) if rs[i] >= 1.0]
    mfar = min(far) if far else (None, None)
    sub = min((D[i], rs[i]) for i in range(n) if rs[i] >= 0.02 - 1e-12)
    return {
        'min_D_r_ge_0.02': '%.12g' % sub[0], 'argmin_r_ge_0.02': '%.6g' % sub[1],
        'min_D': '%.12g' % D[imin], 'argmin_r': '%.6g' % rs[imin],
        'argmin_is_first_grid_point': imin == 0,
        'nearest_ordinate_to_argmin': '%.4f' % min(abs(g - rs[imin]) for g in zeros),
        'n_nonpositive': sum(1 for x in D if x <= 0),
        'D_first': '%.12g' % D[0], 'D_first_over_r2': '%.12g' % (D[0] / rs[0] ** 2),
        'min_D_over_r_ge_1': '%.12g' % mfar[0] if mfar[0] is not None else None,
        'argmin_r_ge_1': '%.6g' % mfar[1] if mfar[1] is not None else None,
        'smallest_local_minima': [{'r': '%.6g' % r, 'D': '%.12g' % d} for d, r in lm[:4]],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--workers', type=int)
    ap.add_argument('--zeros')
    args = ap.parse_args()
    quick = args.quick
    workers = 1 if quick else (args.workers or min(8, os.cpu_count() or 1))
    out_path = args.out or ('/tmp/polya_primes_quick.json' if quick else default_out('polya_primes.json'))
    if workers > 1:
        pool = get_context('fork').Pool(workers)

        def pmap(f, tasks):
            return pool.map(f, tasks, chunksize=1)
    else:
        def pmap(f, tasks):
            return [f(t) for t in tasks]
    zeros_s = open(find_zeros(args.zeros)).read().split()
    zeros = [float(x) for x in zeros_s]
    result = {'script': 'polya_primes.py', 'quick': quick, 'workers': workers}

    # ---------------- (1) identity
    N = 10 ** 7 if quick else 10 ** 10
    pairs = [(3.0, 7.0), (3.0, 12.0), (1.5, 7.0), (1.5, 20.0)]
    log('identity: sieve to N = %d for %s' % (N, pairs))
    sums, psi_N = prime_sums(N, pairs, pmap)
    mp.mp.dps = 30
    ident = []
    for (a, r), sN in zip(pairs, sums):
        am, rm = mp.mpf(a), mp.mpf(r)
        lhs = mp.re(zl(mp.mpc(am, rm))) - zl(am)
        tail = 2 * (RS * a * N ** (1 - a) / (a - 1) - psi_N * N ** (-a))
        main_tail = mp.mpf(N) ** (1 - am) / (am - 1) - mp.re(mp.mpf(N) ** (1 - am + 1j * rm) / (am - 1 - 1j * rm))
        diff = lhs - sN
        ident.append({
            'alpha': a, 'r': r, 'N': N, 'lhs_zeta': mp.nstr(lhs, 20), 'rhs_partial_sum': '%.16g' % sN,
            'lhs_minus_rhs': mp.nstr(diff, 6), 'rigorous_tail_bound': '%.4g' % tail,
            'within_rigorous_bound_plus_float64_allowance': bool(abs(diff) <= tail + 1e-14 * abs(float(lhs))),
            'float64_allowance': '1e-14 relative (the sieve sum has about 4.5e8 terms accumulated in double precision)',
            'agree_digits_relative': round(float(-mp.log10(abs(diff) / abs(lhs))), 1) if diff != 0 else 99.0,
            'pnt_main_term_of_tail': mp.nstr(main_tail, 8),
            'lhs_minus_rhs_after_pnt_completion': mp.nstr(diff - main_tail, 6),
        })
        log('identity alpha=%g r=%g: lhs %s rhs %.15g diff %s (rigorous tail bound %.3g; after the PNT main term %s)' % (
            a, r, mp.nstr(lhs, 15), sN, mp.nstr(diff, 3), tail, mp.nstr(diff - main_tail, 3)))
    result['identity'] = {'N': N, 'psi_N': '%.10g' % psi_N, 'psi_N_over_N': '%.10g' % (psi_N / N), 'cases': ident}

    # ---------------- (2) alpha_1 and the bracket minima
    grid = make_grid(quick)
    log('bracket grid has %d points' % len(grid))
    scan_alphas = ['1.2', '2', '6'] + (['7.7', '7.8', '8', '12'] if quick else ['%.4f' % (1.05 + 0.05 * k) for k in range(0, 220)])
    scan = pmap(t_scan, [(a, grid) for a in scan_alphas])
    result['bracket_minima'] = {'claimed': {'1.2': {'min': '-4.83', 'r': '1.58'}, '2': {'min': '-0.74'}, '6': {'min': '-0.011'}},
                                'scan': [{'alpha': a, 'min': m, 'argmin_r': r} for a, m, r in scan]}
    for a, m, r in scan[:3]:
        log('bracket minimum at alpha=%s: %s at r=%s' % (a, m, r))
    neg = [float(a) for a, m, r in scan if mp.mpf(m) < mp.mpf('-1e-25')]
    pos = [float(a) for a, m, r in scan if mp.mpf(m) >= mp.mpf('-1e-25')]
    result['bracket_minima']['scan_summary'] = {'largest_alpha_with_negative_min': max(neg), 'smallest_alpha_with_nonnegative_min': min(pos),
                                                'sign_changes_once': min(pos) > max(neg)}
    log('scan: largest alpha with a negative minimum %g, smallest alpha scanned with nonnegative minimum %g' % (max(neg), min(pos)))
    lo, hi = mp.mpf(6), mp.mpf(9)
    assert bracket_min(lo, grid)[0] < mp.mpf('-1e-25') and bracket_min(hi, grid)[0] >= 0
    for it in range(25 if quick else 46):
        mid = (lo + hi) / 2
        if bracket_min(mid, grid)[0] < mp.mpf('-1e-25'):
            lo = mid
        else:
            hi = mid
    mp.mp.dps = 30
    a1 = (lo + hi) / 2
    aA = mp.findroot(lambda a: 2 / a ** 3 + 2 / (a - 1) ** 3 + mp.polygamma(2, a / 2) / 8, mp.mpf('7.789'))
    above = []
    big = [float(x) for x in np.logspace(math.log10(200), 10, 80)]
    for da in ['1e-6', '0.01', '0.2', '2', '12', '92']:
        a = a1 + mp.mpf(da)
        m, r = bracket_min(a, grid)
        mb = min(bracket(a, mp.mpf(x)) for x in big)
        above.append({'alpha': mp.nstr(a, 12), 'min_over_grid': mp.nstr(m, 6), 'argmin_r': '%.6g' % r, 'min_for_r_in_200_to_1e10': mp.nstr(mb, 6)})
    result['alpha_1'] = {'bisection': mp.nstr(a1, 14), 'bisection_width': mp.nstr(hi - lo, 3), 'root_of_A_second_derivative': mp.nstr(aA, 14),
                         'claimed': '7.789...', 'checks_above': above,
                         'note': 'the r grid reaches down to 1e-7, so the dips just below the threshold (at r about 1e-4) are seen; the sign test uses a depth of 1e-25, which limits the bisection to about 1e-11 in alpha; the root of A\'\'(alpha) = 0 is the exact threshold when the dip at r -> 0 is the last one to close'}
    log('alpha_1 by bisection %s (width %s); root of A\'\'(alpha) = 0: %s' % (mp.nstr(a1, 12), mp.nstr(hi - lo, 2), mp.nstr(aA, 12)))
    below = []
    for da in ['-1e-6', '-1e-3', '-0.1']:
        a = a1 + mp.mpf(da)
        m, r = bracket_min(a, grid)
        below.append({'alpha': mp.nstr(a, 12), 'min': mp.nstr(m, 6), 'argmin_r': '%.6g' % r})
    result['alpha_1']['checks_below'] = below

    # ---------------- (3) Lemma Bbound
    mp.mp.dps = 40
    a1sum_diff = mp.diff(xi_logder, mp.mpf(1) / 2, 1) / 2
    gam = [mp.mpf(x) for x in zeros_s]
    K = len(gam)
    partial = mp.fsum(1 / g ** 2 for g in gam)
    delta = 2 * mp.pi / mp.log(gam[-1] / (2 * mp.pi))
    Tc = gam[-1] + delta / 2
    tail = (mp.log(Tc / (2 * mp.pi)) + 1) / (2 * mp.pi * Tc)
    lem = {'a_1_from_derivative': mp.nstr(a1sum_diff, 20), 'chapter_a_1': '0.0231049931',
           'zero_sum_K': K, 'zero_sum_partial': mp.nstr(partial, 12), 'zero_sum_smooth_tail': mp.nstr(tail, 8),
           'zero_sum_plus_tail': mp.nstr(partial + tail, 12), 'relative_difference_of_the_two': mp.nstr(abs(partial + tail - a1sum_diff) / a1sum_diff, 3),
           '2_a_1': mp.nstr(2 * a1sum_diff, 12), '4_a_1': mp.nstr(4 * a1sum_diff, 12)}
    rows = []
    for a in ['0.501', '0.55', '0.6', '0.75', '0.9', '0.95', '1']:
        Bv = B_val(a)
        ratio = Bv / (mp.mpf(a) - mp.mpf(1) / 2)
        rows.append({'alpha': a, 'B_alpha': mp.nstr(Bv, 12), 'B_over_alpha_minus_half': mp.nstr(ratio, 9), 'below_2a1': bool(ratio < 2 * a1sum_diff)})
    lem['rows'] = rows
    lem['chapter_ratios'] = {'0.55': '0.046210', '0.75': '0.046205', '1': '0.046191'}
    mp.mp.dps = 80
    lem['B_1_limit_check'] = mp.nstr(mp.re(xi_logder(mp.mpf(1) + mp.mpf('1e-30'))), 25) + ' versus closed form ' + mp.nstr(B_val('1'), 25)
    mp.mp.dps = 40
    cmin_max = mp.mpf(1) / 4
    Hh = mp.mpf(H_VER)
    S_lead = (mp.log(Hh / (2 * mp.pi)) + 1) / (2 * mp.pi * Hh)
    S_chapter = (mp.log(Hh / (2 * mp.pi)) + 2) / (2 * mp.pi * Hh)
    lem['almost_all_arithmetic'] = {
        'kappa_over_c_min_at_alpha_1': mp.nstr(1 / (1 + cmin_max ** 2), 8), 'four_fifths': '0.8',
        'B_alpha_bound_coefficient_4a1': mp.nstr(4 * a1sum_diff, 10), 'printed': '0.0924', 'printed_is_below_true_bound': bool(mp.mpf('0.0924') < 4 * a1sum_diff),
        'margin_0.4_minus_4a1': mp.nstr(mp.mpf('0.4') - 4 * a1sum_diff, 8),
        'S_leading_term': mp.nstr(S_lead, 6), 'S_with_plus_2': mp.nstr(S_chapter, 6), 'chapter_S_bound': '1.6e-12'}
    result['lemma_bbound'] = lem
    log('a_1 = %s (derivative), zero sum + tail %s; 2a_1 = %s; ratios %s' % (mp.nstr(a1sum_diff, 12), mp.nstr(partial + tail, 12), mp.nstr(2 * a1sum_diff, 10), [r['B_over_alpha_minus_half'] for r in rows]))
    log('4 a_1 = %s against the printed 0.0924; S = %s (leading), %s (with +2)' % (mp.nstr(4 * a1sum_diff, 8), mp.nstr(S_lead, 5), mp.nstr(S_chapter, 5)))

    # ---------------- (4) D on a grid of base points, heights in (0, 200]
    if quick:
        step, n, chunk = '0.5', 400, 100
        alphas15 = ['1.000000000000', '4.400000000000', '7.800000000000']
        extra = ['1.001']
    else:
        step, n, chunk = '0.01', 20000, 500
        alphas15 = ['%.12f' % (1 + 6.8 * j / 14) for j in range(15)]
        extra = ['1.001', '1.01', '1.1', '1.2', '1.5', '2', '3', '4', '5', '6', '7', '7.789', '7.8']
    allalpha = alphas15 + extra
    idx = list(range(1, n + 1))
    rs = [i * float(step) for i in idx]
    tasks, owner = [], []
    for a in allalpha:
        for c in range(0, n, chunk):
            tasks.append((a, step, idx[c:c + chunk], 20))
            owner.append(a)
    log('D grid: %d base points, %d heights each (step %s), %d tasks' % (len(allalpha), n, step, len(tasks)))
    res = pmap(t_D, tasks)
    series = {a: [] for a in allalpha}
    for a, r in zip(owner, res):
        series[a].extend(r)
    dg = {}
    for a in allalpha:
        an = analyse(rs, series[a], zeros)
        dg[a] = an
    result['D_grid'] = {'step': step, 'heights': '(0, %g]' % rs[-1], 'fifteen_base_points': alphas15, 'extra_base_points': extra, 'cases': dg,
                        'all_argmin_at_first_grid_point': all(dg[a]['argmin_is_first_grid_point'] for a in allalpha),
                        'total_nonpositive_values': sum(dg[a]['n_nonpositive'] for a in allalpha)}
    log('D grid: argmin at the first grid point for every base point: %s; nonpositive values in total: %d' % (
        result['D_grid']['all_argmin_at_first_grid_point'], result['D_grid']['total_nonpositive_values']))
    for a in allalpha[:3]:
        log('  alpha=%s: min %s at r=%s, D/r^2 at first height %s, smallest later local minimum %s' % (a, dg[a]['min_D'], dg[a]['argmin_r'], dg[a]['D_first_over_r2'], dg[a]['smallest_local_minima'][:1]))

    # ---------------- (5) (0, 600]
    if quick:
        step6, n6, chunk6, al6 = '0.5', 400, 100, ['0.9']
    else:
        step6, n6, chunk6, al6 = '0.01', 60000, 500, ['0.9', '0.95', '1', '1.001', '2', '4', '7.8']
    idx6 = list(range(1, n6 + 1))
    rs6 = [i * float(step6) for i in idx6]
    tasks, owner = [], []
    for a in al6:
        for c in range(0, n6, chunk6):
            tasks.append((a, step6, idx6[c:c + chunk6], 20))
            owner.append(a)
    log('(0, 600] runs: %d base points, %d heights each, %d tasks' % (len(al6), n6, len(tasks)))
    res = pmap(t_D, tasks)
    s6 = {a: [] for a in al6}
    for a, r in zip(owner, res):
        s6[a].extend(r)
    d600 = {}
    for a in al6:
        D = s6[a]
        an = analyse(rs6, D, zeros)
        sel = [i for i in range(len(D)) if rs6[i] > 2.0]
        j = min(sel, key=lambda i: D[i])
        an['min_over_r_gt_2'] = '%.9g' % D[j]
        an['argmin_over_r_gt_2'] = '%.6g' % rs6[j]
        refined = []
        lms = sorted((D[i], i) for i in sel if 0 < i < len(D) - 1 and D[i] < D[i - 1] and D[i] <= D[i + 1])[:3]
        if rs6[j] == rs6[sel[0]]:
            lms.append((D[j], j))
        for d0, i in lms:
            if 0 < i < len(D) - 1:
                rr = minimize_scalar(lambda x: D_scalar(a, x), bounds=(rs6[i - 1], rs6[i + 1]), method='bounded', options={'xatol': 1e-9})
                refined.append({'r': '%.8g' % rr.x, 'D': '%.9g' % rr.fun})
        an['refined_local_minima_gt_2'] = refined
        i2 = idx6.index(int(round(2.0 / float(step6))))
        an['D_at_2'] = '%.9g' % D[i2]
        an['D_at_2_plus_step'] = '%.9g' % D[i2 + 1]
        d600[a] = an
        log('(0,600] alpha=%s: min %s at r=%s, nonpositive %d, min over r>2 %s at r=%s, D(2) %s' % (
            a, an['min_D'], an['argmin_r'], an['n_nonpositive'], an['min_over_r_gt_2'], an['argmin_over_r_gt_2'], an['D_at_2']))
    result['D_600'] = {'step': step6, 'heights': '(0, %g]' % rs6[-1], 'cases': d600,
                       'claimed_min_over_r_gt_2': {'0.9': '4e-4', '0.95': '4e-4', '1': '5e-4'}}
    result['seconds'] = round(time.time() - T0, 1)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=1)
    log('written', out_path)


if __name__ == '__main__':
    main()
