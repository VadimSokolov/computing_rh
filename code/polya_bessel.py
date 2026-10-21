# polya_bessel.py: the Bessel model f(nu) = K_nu(2 pi) and the compensated clock (cosine) condition, Table pol:tab:bessel in Section pol:sec:lemma
# (Is the 2 sin^2 lemma useful with regard to Polya's work?) and the numerical sentences after Theorem pol:thm:exactlemma (the cosine condition
# for the Bessel model and for xi, Figure pol:fig:two right) of Chapter 16 (Polya and the zeros of entire functions). Added for the book in
# October 2026, written from the mathematics (the chapter author's own scripts were not read, run or copied).
# What it checks. (1) Table bessel: E_j(nu) = (-1)^j d^(2+j)/dnu^(2+j) log K_nu(2 pi) for nu = 0.25, 1, 3, 5, 8 and j = 0..5, computed two ways:
# mpmath diff of log K_nu(2 pi) at 50 digits (central differences, step 2^-(prec+10), working precision (prec+20)(n+1)), and the Cauchy integral
# on the circle |z - nu| = 3 (trapezoid rule at N = 192 nodes, 60 digits, conjugate symmetry halves the work). K_z(2 pi) is entire in z with zeros
# only on the imaginary axis at distance at least 8.99 from the real line segment, so the rule converges like (3/8.99)^192; the logarithm is
# continued along the circle by unwrapping its argument (it must come back to 0 at the real axis, which is checked). The number of digits on
# which the two agree is written for every entry, with the chapter's printed value and whether our value rounds to it (4 decimals for j = 0, two
# significant digits for j >= 1). The signs claimed after the table (third entry negative for nu <= 5, a later entry negative at nu = 8, the
# j = 1 column positive) are tested on the computed values. (2) The compensated condition D_alpha(r) = Re g'(alpha + i r) - g'(alpha) >= 0,
# g = log f: for the Bessel model (g' = d/dnu log K_nu(2 pi), computed as (d/dnu K)/K with mpmath diff at complex order) at alpha = 0.25, 1, 3, 8
# on 0 < r <= 60 (step 0.02), and for xi (g' = xi'/xi = A + zeta'/zeta) at alpha = 1, 0.75, 0.6, 0.55 on 0 < r <= 100 (step 0.01, finer than the
# Lorentzian width 0.05 of the nearest zeros). For each case: the minimum over the grid and where it is attained, the number of grid points with
# D <= 0, D/r^2 at the first grid point against the Taylor prediction -g'''(alpha)/2 = E_1/2 (Bessel: from the table), the smallest local minima
# away from r = 0. A spot check of the Bessel derivative uses the integral d/dnu K_nu(x) = int_0^oo t sinh(nu t) exp(-x cosh t) dt (mp.quad at
# 80 digits) at a few complex orders. (3) The Figure pol:fig:two right series: D_alpha(r) for xi at alpha = 1, 0.75, 0.6 and r = 0.15 i, i = 1..400.
# For xi at alpha = 1 the pole cancels: B_1 = 1 + gamma/2 - (log pi)/2 - log 2 (also checked as the limit alpha -> 1).
# Error control. Every derivative is computed twice (Table), every Bessel derivative of the cosine runs uses working precision 30 digits against
# values of D as small as 1e-8, and the spot check recomputes three values by an integral. Nothing here is a proof: the cosine condition is
# checked on a grid, not for all r.
# Output: data/polya_bessel.json (numbers as strings), the log is the standard output (data/polya_bessel.log).
# Usage: python3 polya_bessel.py [--quick] [--out FILE] [--workers W]; --quick runs a small case in one process (about 25 s, writes to /tmp unless
# --out is given); the full case is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8 -g polya_bessel.json code/polya_bessel.py -- --out polya_bessel.json
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import mpmath as mp

T0 = time.time()
NUS = ['0.25', '1', '3', '5', '8']
ORDERS = list(range(2, 8))   # derivative orders 2..7 give j = 0..5
CHAPTER_TABLE = {
    '0.25': ['0.1482', '6.7e-4', '-2.7e-3', '-9.8e-5', '3.9e-4', '3.6e-5'],
    '1': ['0.1469', '2.6e-3', '-2.5e-3', '-3.7e-4', '3.2e-4', '1.3e-4'],
    '3': ['0.1373', '6.6e-3', '-1.3e-3', '-6.8e-4', '-1.4e-5', '1.5e-4'],
    '5': ['0.1223', '8.0e-3', '-1.6e-4', '-4.6e-4', '-1.6e-4', '6.9e-6'],
    '8': ['0.0991', '7.1e-3', '5.6e-4', '-6.9e-5', '-7.7e-5', '-3.5e-5'],
}


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def default_out(name):
    return os.path.join('data', name) if os.path.isdir('data') else name


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


def agree_digits(a, b):
    if a == b:
        return 999.0
    return float(-mp.log10(abs(a - b) / max(abs(a), abs(b))))


# ---------------------------------------------------------------- tasks (module level so that they can be pickled)
def t_diff(args):
    nu, n, dps = args
    mp.mp.dps = dps
    v = mp.diff(lambda t: mp.log(mp.besselk(t, 2 * mp.pi)), mp.mpf(nu), n)
    return mp.nstr(v, dps - 3)


def t_circle(args):
    nu, R, N, k, dps = args
    mp.mp.dps = dps
    z = mp.mpf(nu) + R * mp.expjpi(2 * (k + mp.mpf(1) / 2) / N)
    lg = mp.log(mp.besselk(z, 2 * mp.pi))
    return k, mp.nstr(mp.re(lg), dps - 3), mp.nstr(mp.im(lg), dps - 3)


def gprime_bessel(nu):
    return mp.diff(lambda t: mp.besselk(t, 2 * mp.pi), nu, 1) / mp.besselk(nu, 2 * mp.pi)


def t_bessel_D(args):
    alpha, step, idx, dps = args
    mp.mp.dps = dps
    a = mp.mpf(alpha)
    ga = mp.re(gprime_bessel(a))
    out = []
    for i in idx:
        v = gprime_bessel(mp.mpc(a, i * mp.mpf(step)))
        out.append(mp.nstr(mp.re(v) - ga, 18))
    return out


def t_xi_D(args):
    alpha, step, idx, dps = args
    mp.mp.dps = dps
    a = mp.mpf(alpha)
    B = B_val(a)
    out = []
    for i in idx:
        s = mp.mpc(a, i * mp.mpf(step))
        out.append(mp.nstr(mp.re(xi_logder(s)) - B, 18))
    return out


def t_integral_check(args):
    # d/dnu K_nu(x) / K_nu(x) at nu = alpha + i r: the mpmath diff against the integral int_0^oo t sinh(nu t) exp(-x cosh t) dt
    alpha, r, dps = args
    mp.mp.dps = dps
    nu = mp.mpc(mp.mpf(alpha), mp.mpf(r))
    x = 2 * mp.pi
    gp1 = gprime_bessel(nu)
    mp.mp.dps = 80
    nu = mp.mpc(mp.mpf(alpha), mp.mpf(r))
    x = 2 * mp.pi
    pts = [0] + [mp.mpf(j) / 2 for j in range(1, 17)]
    num = mp.quad(lambda t: t * mp.sinh(nu * t) * mp.exp(-x * mp.cosh(t)), pts)
    den = mp.quad(lambda t: mp.cosh(nu * t) * mp.exp(-x * mp.cosh(t)), pts)
    gp2 = num / den
    mp.mp.dps = dps
    return alpha, r, mp.nstr(gp1, 20), mp.nstr(gp2, 20), float(-mp.log10(abs(gp1 - gp2) / abs(gp2)))


# ---------------------------------------------------------------- analysis helpers
def unwrap(ims):
    out = []
    prev = mp.mpf(0)
    for x in ims:
        m = mp.nint((prev - x) / (2 * mp.pi))
        y = x + 2 * mp.pi * m
        out.append(y)
        prev = y
    return out


def printed_match(j, v, chapter):
    if j == 0:
        return '%.4f' % v == '%.4f' % float(chapter)
    return '%.1e' % v == '%.1e' % float(chapter)


def analyse(rs, D):
    n = len(D)
    imin = min(range(n), key=lambda i: D[i])
    lm = sorted((D[i], rs[i]) for i in range(1, n - 1) if D[i] < D[i - 1] and D[i] <= D[i + 1])
    far = [(D[i], rs[i]) for i in range(n) if rs[i] >= 1.0]
    mfar = min(far) if far else (None, None)
    return {
        'grid_points': n,
        'r_first': '%.6g' % rs[0], 'r_last': '%.6g' % rs[-1],
        'min_D': '%.12g' % D[imin], 'argmin_r': '%.6g' % rs[imin],
        'argmin_is_first_grid_point': imin == 0,
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
    args = ap.parse_args()
    quick = args.quick
    workers = 1 if quick else (args.workers or min(8, os.cpu_count() or 1))
    out_path = args.out or ('/tmp/polya_bessel_quick.json' if quick else default_out('polya_bessel.json'))
    if workers > 1:
        pool = get_context('fork').Pool(workers)

        def pmap(f, tasks):
            return pool.map(f, tasks, chunksize=1)
    else:
        def pmap(f, tasks):
            return [f(t) for t in tasks]

    if quick:
        nus, dps_d, dps_c, N, R = ['1', '5'], 30, 40, 64, 3
        b_alphas, b_step, b_n, b_chunk = ['1', '0.25'], '0.5', 20, 10
        x_alphas, x_step, x_n, x_chunk = ['0.75'], '0.25', 80, 20
        fig_n, fig_chunk = 40, 20
        dps_cos = 20
    else:
        nus, dps_d, dps_c, N, R = NUS, 50, 60, 192, 3
        b_alphas, b_step, b_n, b_chunk = ['0.25', '1', '3', '8'], '0.02', 3000, 25
        x_alphas, x_step, x_n, x_chunk = ['1', '0.75', '0.6', '0.55'], '0.01', 10000, 100
        fig_n, fig_chunk = 400, 50
        dps_cos = 30
    result = {'script': 'polya_bessel.py', 'quick': quick, 'workers': workers}

    # ---------------- (1) Table bessel
    log('Table bessel: mpmath diff at %d digits and Cauchy circle R=%s, N=%d at %d digits' % (dps_d, R, N, dps_c))
    res_d = pmap(t_diff, [(nu, n, dps_d) for nu in nus for n in ORDERS])
    res_c = pmap(t_circle, [(nu, R, N, k, dps_c) for nu in nus for k in range(N // 2)])
    mp.mp.dps = dps_c
    half = N // 2
    table = {}
    E1 = {}
    for ni, nu in enumerate(nus):
        vals = res_c[ni * half:(ni + 1) * half]
        re_ = [mp.mpf(v[1]) for v in vals]
        im_ = unwrap([mp.mpf(v[2]) for v in vals])
        closing_arg = float(im_[-1])
        assert abs(closing_arg) < 0.5, ('the continued argument does not return to 0', nu, closing_arg)
        row = []
        for j, n in enumerate(ORDERS):
            S = mp.mpf(0)
            for k in range(half):
                th = 2 * mp.pi * (k + mp.mpf(1) / 2) / N
                S += mp.re(mp.mpc(re_[k], im_[k]) * mp.expj(-n * th))
            cauchy = mp.factorial(n) / mp.mpf(R) ** n * 2 * S / N
            diff_v = mp.mpf(res_d[ni * len(ORDERS) + j])
            sgn = (-1) ** j
            e_c, e_d = sgn * cauchy, sgn * diff_v
            ag = agree_digits(e_c, e_d)
            ch = CHAPTER_TABLE.get(nu, [None] * 6)[j]
            ours = float(e_c)
            row.append({
                'j': j, 'derivative_order': n,
                'E_mpdiff': mp.nstr(e_d, dps_d - 3), 'E_cauchy': mp.nstr(e_c, dps_c - 3), 'agree_digits': round(ag, 1),
                'value': mp.nstr(e_c, 15), 'chapter': ch,
                'rounds_to_chapter': printed_match(j, ours, ch) if ch is not None else None,
                'sign': '+' if ours > 0 else '-',
            })
            if j == 1:
                E1[nu] = e_c
        table[nu] = {'closing_arg': '%.3g' % closing_arg, 'entries': row}
        log('nu=%s: ' % nu + ', '.join('%s%s' % (e['value'][:9], '' if e['rounds_to_chapter'] in (True, None) else ' (differs, chapter %s)' % e['chapter']) for e in row)
            + ' | agree digits min %.1f' % min(e['agree_digits'] for e in row))
    result['table_bessel'] = table
    signs = {nu: ''.join(e['sign'] for e in table[nu]['entries']) for nu in nus}
    claims = {'signs_j0_to_j5': signs}
    if all(nu in table for nu in ['0.25', '1', '3', '5']):
        claims['third_entry_negative_for_nu_le_5'] = all(table[nu]['entries'][2]['sign'] == '-' for nu in ['0.25', '1', '3', '5'])
    if '8' in table:
        claims['later_entry_negative_at_nu_8'] = [e['j'] for e in table['8']['entries'] if e['sign'] == '-' and e['j'] >= 3]
    claims['j1_column_positive_everywhere'] = all(table[nu]['entries'][1]['sign'] == '+' for nu in nus)
    claims['some_entry_negative_at_every_nu'] = all('-' in signs[nu] for nu in nus)
    claims['all_30_entries_round_to_chapter'] = all(e['rounds_to_chapter'] in (True, None) for nu in nus for e in table[nu]['entries'])
    claims['min_agree_digits_over_table'] = min(e['agree_digits'] for nu in nus for e in table[nu]['entries'])
    result['table_bessel_claims'] = claims
    log('claims:', json.dumps(claims))

    # ---------------- spot check of the complex order derivative by an integral
    spots = [('0.25', '60'), ('3', '25'), ('8', '5')] if not quick else [('1', '3')]
    chk = pmap(t_integral_check, [(a, r, dps_cos) for a, r in spots])
    result['bessel_derivative_integral_check'] = [{'alpha': a, 'r': r, 'mpdiff': g1, 'integral': g2, 'agree_digits': round(dg, 1)} for a, r, g1, g2, dg in chk]
    log('integral spot checks (alpha, r, agree digits):', [(a, r, round(dg, 1)) for a, r, g1, g2, dg in chk])
    t = time.time()
    t_bessel_D(('0.25', '1', [60], dps_cos))
    log('one Bessel D evaluation at nu = 0.25 + 60 i takes %.3f s' % (time.time() - t))

    # ---------------- (2) cosine condition, Bessel model
    def run_D(task, alphas, step, n, chunk):
        idx = list(range(1, n + 1))
        tasks, owner = [], []
        for a in alphas:
            for c in range(0, n, chunk):
                tasks.append((a, step, idx[c:c + chunk], dps_cos))
                owner.append(a)
        res = pmap(task, tasks)
        series = {a: [] for a in alphas}
        for a, r in zip(owner, res):
            series[a].extend(float(x) for x in r)
        return series

    log('cosine condition, Bessel model, alpha in %s, step %s, %d points' % (b_alphas, b_step, b_n))
    sb = run_D(t_bessel_D, b_alphas, b_step, b_n, b_chunk)
    rs_b = [i * float(b_step) for i in range(1, b_n + 1)]
    bes = {}
    for a in b_alphas:
        an = analyse(rs_b, sb[a])
        if a in E1:
            an['taylor_prediction_E1_over_2'] = mp.nstr(E1[a] / 2, 12)
        bes[a] = an
        log('Bessel alpha=%s: min D %s at r=%s (first point: %s), nonpositive %d, D/r^2 at first point %s, prediction %s' % (
            a, an['min_D'], an['argmin_r'], an['argmin_is_first_grid_point'], an['n_nonpositive'], an['D_first_over_r2'], an.get('taylor_prediction_E1_over_2')))
    result['bessel_cosine'] = {'step': b_step, 'range': '(0, %s]' % ('%g' % rs_b[-1]), 'cases': bes}

    # ---------------- (2b) cosine condition for xi
    log('cosine condition, xi, alpha in %s, step %s, %d points' % (x_alphas, x_step, x_n))
    mp.mp.dps = dps_cos
    sx = run_D(t_xi_D, x_alphas, x_step, x_n, x_chunk)
    rs_x = [i * float(x_step) for i in range(1, x_n + 1)]
    xis = {}
    for a in x_alphas:
        an = analyse(rs_x, sx[a])
        an['B_alpha'] = mp.nstr(B_val(a), 15)
        small = t_xi_D((a, '0.001', [1, 2, 5], 30))
        an['D_over_r2_at_r_0.001_0.002_0.005'] = ['%.10g' % (float(v) / (i * 0.001) ** 2) for v, i in zip(small, [1, 2, 5])]
        xis[a] = an
        log('xi alpha=%s: min D %s at r=%s (first point: %s), nonpositive %d, D/r^2 at first point %s, at r=0.001 %s' % (
            a, an['min_D'], an['argmin_r'], an['argmin_is_first_grid_point'], an['n_nonpositive'], an['D_first_over_r2'], an['D_over_r2_at_r_0.001_0.002_0.005'][0]))
    result['xi_cosine'] = {'step': x_step, 'range': '(0, %s]' % ('%g' % rs_x[-1]), 'cases': xis}
    mp.mp.dps = 50
    b1_closed = B_val('1')
    mp.mp.dps = 80
    s_near = mp.mpf(1) + mp.mpf('1e-30')
    b1_limit = mp.re(xi_logder(s_near))
    result['B_1_closed_form_and_limit'] = {'closed': mp.nstr(b1_closed, 30), 'limit_at_1_plus_1e-30': mp.nstr(b1_limit, 30)}

    # ---------------- (3) Figure 2 right series
    log('Figure 2 right series: alpha = 1, 0.75, 0.6, r = 0.15 i, i = 1..%d' % fig_n)
    mp.mp.dps = dps_cos
    sf = run_D(t_xi_D, ['1', '0.75', '0.6'], '0.15', fig_n, fig_chunk)
    result['figure2_right'] = {
        'r': ['%.2f' % (0.15 * i) for i in range(1, fig_n + 1)],
        'D': {a: ['%.12g' % v for v in sf[a]] for a in sf},
        'note': 'D_alpha(r) = Re xi\'/xi(alpha + i r) - xi\'/xi(alpha), r = 0.15 i; the chapter figure may use another 400 point grid',
    }
    result['seconds'] = round(time.time() - T0, 1)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=1)
    log('written', out_path)


if __name__ == '__main__':
    main()
