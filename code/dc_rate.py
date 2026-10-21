# dc_rate.py: the rate of the ladder in the chapter "Deconvolution as a Damped Weil Criterion" (label prefix dc:, Chapter 20 of the book): Table dc:tab:rate (Table 1.2), the ratios 1.19 and 1.04 of line 147, Problem dc:prob:rate (line 179, "a correction of order a^(-1)") and the trial function and the lower bound in the proof of Theorem dc:thm:ladder(b). Group D2 of the number checks.
# Added for the book in October 2026, written from the mathematics (no script of the chapter's author was read or run). It uses code/yo_matrix.py (matrices of the damped Weil form, group Y2) and code/dc_ladder.py (the ladder study: the sizes and the extrapolation to the converged value) as they are.
# Notation. sigma = 1/2 + y, f = xi'/xi, w_y(theta) = Re f(sigma + i theta), m(y) = 2 w_y(0) = 2 f(sigma), w_y''(0) = -f''(sigma), c_y = pi^2 w_y''(0)/4 (the upper constant of Theorem dc:thm:ladder(b)), S(a, y) = a^2 (lambda(a, y) - m(y)) the scaled excess, lambda(a, y) the least eigenvalue of Q_y on [-a, a] (see dc_ladder.py).
# Statements checked (chapter lines). R1 Table 1.2 (154 to 156) and the ratios of line 147: S(a, y) for a in {1, 2, 3, 4} and y in {0.5, 0.2, 0.1} from the converged lambda of the ladder study (Fourier sizes up to 768, Legendre up to 128, power law extrapolation, the residual change stated), against c_y, with the printed M = 24 entries reproduced for comparison and a replacement table to three digits; the constants c_y themselves (5.485, 2.200, 1.100 in units of 1e-4). R2 Problem 1.2: S(a, y) on a in {1.5, 2, ..., 6} (prime powers to e^12 = 1.6e5), least squares fits S = C + D/a + E/a^2 (and variants), C/c_y with an uncertainty from the data errors and from the choice of window and degree. R3 the cosine trial function h = cos(pi u/2a) of the proof of (b): Q_y(h)/||h||^2 = 2 <w_y>_h from the prime side (W of e^(-y|x|) (h * h), closed form of h * h, ||h||^2 = a) against m(y) + c_y/a^2, the remainder a^2 (T - m) - c_y, its product with a, and its own fit (informational only: for y below 1/2 the remainder oscillates in a because |h^(theta)|^2 ~ 4 k^2 cos^2(a theta)/theta^4 has a slowly decaying tail that meets the peaks of w_y at the zeros, of height 1/y). R4 a computable lower bound from the proof: with psi(theta) = w_y(theta) - w_y(0) and the largest eigenvalue Lambda_0(kappa) of the time and band limiting operator (Slepian and Pollak, computed here by Gauss Legendre discretisation of the sinc kernel on [-1, 1]), the chapter's version lambda - m >= 2 c_1 kappa^2 (1 - Lambda_0(kappa))/a^2 for kappa <= a with c_1 = inf psi/min(theta^2, 1), and a sharper layer cake version lambda - m >= 2 int psi~'(r) (1 - Lambda_0(a r)) dr with psi~ the nondecreasing minorant of psi (both use only P(|theta| <= r) <= Lambda_0(r a)).
# Method. const: f''(sigma) from the Taylor coefficients of f on the circle |s - sigma| = 1/5 (trapezoid rule with 64 points, 40 digits, geometric accuracy because the nearest singularity of f is a zero of xi at distance 14), checked by a four point real difference (h = 1e-6, 120 digits); b_alpha from the module. table and fit: lambda from dc_ladder.estimate (the three largest Fourier sizes give the power law value, the shift when the largest size is dropped is the uncertainty, the Legendre best value is an upper bound that must lie above it). trial: ym.weil_functional with the closed form (h * h)(x) = (2a - |x|) cos(kx)/2 + sin(k|x|)/(2k), k = pi/2a, checked against the Legendre matrix route (c^T B c / c^T c, c the Legendre coefficients of h), against 30 digits and against 1.5 times the nodes. lower: psi on a grid of 600 points in [0, 40] at 25 digits (the part above 20 has its minimum at height 27.67, Chapter 17, group Y4), Lambda_0 on a kappa grid of step 0.005 up to 16 (beyond it 1 - Lambda_0 is set to 0, below 1e-13, which keeps the bound a lower bound), right endpoint rule for the Stieltjes sum.
# Error control and honesty. All eigenvalues are variational upper bounds; the converged value is an extrapolation (the uncertainty is the shift described above, not a bound). The fit uncertainty is the standard deviation of C over 400 random perturbations of the data of the size of their uncertainties (floor 3e-10) and the spread over windows (a >= 1.5, 2, 3) and degrees (1, 2, 3). The lower bound of R4 is a floating point computation on a grid with the assumption that psi stays above its minimum on [20, 40] beyond 40 (the book proves the minimum at the origin up to a huge height but not this quantitative margin): it is a numerical bound, not a certificate. For y < 1/2 the identity Q_y = 2 <w_y> assumes y to the right of the zeros (checked here only through the numerics). Nothing here is a proof.
# Output: data/dc_rate.json (numbers as strings; keys const, table, fit, trial, lower, meta; the rows of lower carry the two sided table: lower bounds, converged value, trial upper bound; a part written again replaces only its own entries) and data/dc_rate.log (appended).
# Usage (repository root): timeout 60 python3 -I code/dc_rate.py --part const|table|fit|trial|lower|all [--quick] [--out FILE] [--procs P] [--y 0.5,0.2] [--a 1,2] [--afit 2,3,4] [--atrial 1,2] [--alower 1,2] [--fourier 96,192] [--legendre 24,32]. --quick runs a small case of every part in under 30 s in one process and writes to /tmp unless --out is given. Locally one value of y takes about 10 to 20 s per part; the complete run with the larger sizes (Fourier 1536, a up to 8) is meant for the Hopper cluster:
#   bash misc/tools/hopper_run.sh -c 16 -t 240 -m 48G -g dc_rate.json -g dc_rate.log code/dc_rate.py code/dc_ladder.py code/yo_matrix.py -- --part all --procs 16 --fourier 96,192,384,768,1536 --afit 1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7 --atrial 1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8 --out dc_rate.json
# (a = 8 with the Fourier size 1536 would need about 15 GB for the prime sum arrays, so the fits stop at a = 7; the trial function has no such limit.)
import argparse
import json
import math
import os
import platform
import socket
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import mpmath
import numpy as np

import yo_matrix as ym
import dc_ladder as dl

mpmath.mp.dps = 50

Y_ALL = ['0.5', '0.375', '0.2', '0.1', '0.05']
Y_RATE = ['0.5', '0.2', '0.1']
A_TABLE = ['1', '2', '3', '4']
A_FIT = ['1.5', '2', '2.5', '3', '3.5', '4', '4.5', '5', '5.5', '6']
A_TRIAL = ['1', '1.5', '2', '2.5', '3', '3.5', '4', '4.5', '5', '5.5', '6']
A_LOWER = ['1', '2', '3', '4', '5', '6']
FOU_RATE = [24, 48, 96, 192, 384, 768]
FOU_FIT = [96, 192, 384, 768]
LEG_RATE = [24, 32, 48, 64, 96, 128]
# Table 1.2 as printed in the chapter (lines 154 to 156), units of 1e-4, a = 1, 2, 3, 4, and the constants of column 2
PRINTED_T = {'0.5': ['6.518', '5.958', '5.798', '5.720'], '0.2': ['2.711', '2.446', '2.369', '2.325'], '0.1': ['1.384', '1.245', '1.203', '1.172']}
PRINTED_C = {'0.5': '5.485', '0.2': '2.200', '0.1': '1.100'}
PRINTED_RATIO = {('0.5', '1'): '1.19', ('0.5', '4'): '1.04'}


def stamp(res, part, seconds, mode, **extra):
    meta = res.setdefault('meta', {})
    meta['script'] = 'dc_rate.py'
    meta[part] = dict(seconds=round(seconds, 1), host=socket.gethostname(), mode=mode, date=time.strftime('%Y-%m-%d %H:%M'),
                      python=platform.python_version(), numpy=np.__version__, mpmath=mpmath.__version__, **extra)


def e4(x, n=5):
    """x in units of 1e-4, n significant digits."""
    return mpmath.nstr(mpmath.mpf(x) * 10 ** 4, n)


# ---------------------------------------------------------------------------------------------------------------------
# const: f = xi'/xi near sigma
# ---------------------------------------------------------------------------------------------------------------------

def xi_logder(s):
    """xi'/xi(s) = 1/s + 1/(s-1) - log(pi)/2 + psi(s/2)/2 + zeta'/zeta(s) (mpmath, current precision, complex or real s)."""
    return 1 / s + 1 / (s - 1) - mpmath.log(mpmath.pi) / 2 + mpmath.digamma(s / 2) / 2 + mpmath.zeta(s, derivative=1) / mpmath.zeta(s)


def taylor_coeffs(sigma, kmax, K, r, dps):
    """c_k = f^(k)(sigma)/k!, k = 0..kmax, by the trapezoid rule on the circle |s - sigma| = r with K points."""
    with mpmath.workdps(dps):
        sg, rr = mpmath.mpf(sigma), mpmath.mpf(r)
        vals = [xi_logder(sg + rr * mpmath.expjpi(2 * mpmath.mpf(m) / K)) for m in range(K)]
        out = []
        for k in range(kmax + 1):
            tot = mpmath.mpc(0)
            for m in range(K):
                tot += vals[m] * mpmath.expjpi(-2 * mpmath.mpf(k * m) / K)
            out.append(tot / K / rr ** k)
        return out


def f2_fourpoint(sigma, h, dps):
    """f''(sigma) from f(x +- 2h) and f(x +- h) only (so that the removable singularity at sigma = 1 is never evaluated): (f(2h) + f(-2h) - f(h) - f(-h)) / (3 h^2) + O(h^2)."""
    with mpmath.workdps(dps):
        sg, hh = mpmath.mpf(sigma), mpmath.mpf(h)

        def F(x):
            return xi_logder(sg + x)

        return (F(2 * hh) + F(-2 * hh) - F(hh) - F(-hh)) / (3 * hh ** 2)


def part_const(res, log, ys):
    t0 = time.time()
    store = res.setdefault('const', {})
    for y in ys:
        with mpmath.workdps(40):
            sigma = mpmath.mpf('0.5') + mpmath.mpf(y)
            cs = taylor_coeffs(sigma, 6, 64, '0.2', 40)
            f0 = cs[0].real
            f2 = 2 * cs[2].real
            f4 = 24 * cs[4].real
            imag = max(abs(c.imag) for c in cs)
            f2r = f2_fourpoint(sigma, '1e-6', 120)
            m = 2 * ym.b_alpha(sigma, dps=40)
            w2 = -f2
            cy = mpmath.pi ** 2 * w2 / 4
        printed = PRINTED_C.get(y)
        ent = dict(sigma=dl.fnum(sigma, 6), m=dl.fnum(m, 25), f0_circle_minus_b=dl.fnum(f0 - m / 2, 3), imaginary_part_of_coefficients=mpmath.nstr(imag, 3),
                   f2_circle=dl.fnum(f2, 25), f2_four_point=dl.fnum(f2r, 25), difference=dl.fnum(f2 - f2r, 3), w2=dl.fnum(w2, 25), c_y=dl.fnum(cy, 25),
                   f4=dl.fnum(f4, 12), c_y_units_1e4=e4(cy, 6))
        if printed:
            ent['c_y_printed_units_1e4'] = printed
            ent['agrees_with_printed'] = bool('%.3f' % float(cy * 10 ** 4) == printed)
        store['y=' + y] = ent
        log('const y = %-5s sigma = %s: m = %s, f\'\'(sigma) = %s (four point %s, difference %s), w\'\'(0) = %s, c_y = pi^2 w\'\'/4 = %se-4 (printed %s)' % (
            y, dl.fnum(sigma, 4), dl.fnum(m, 12), dl.fnum(f2, 14), dl.fnum(f2r, 14), dl.fnum(f2 - f2r, 3), dl.fnum(w2, 12), e4(cy, 6), printed))
    stamp(res, 'const', time.time() - t0, 'mpmath 40 digits')


def cy_of(res, log, y):
    if 'y=' + y not in res.get('const', {}):
        part_const(res, log, [y])
    return mpmath.mpf(res['const']['y=' + y]['c_y'])


# ---------------------------------------------------------------------------------------------------------------------
# table (R1)
# ---------------------------------------------------------------------------------------------------------------------

def part_table(res, log, ys, avals, legs, fous, procs):
    t0 = time.time()
    jobs = [(a, y, 'legendre', M, 15) for y in ys for a in avals for M in legs] + [(a, y, 'fourier', N, 15) for y in ys for a in avals for N in fous]
    vals = dl.collect(dl.run_jobs(dl.series_job, jobs, procs))
    store = res.setdefault('table', {})
    for y in ys:
        cy = cy_of(res, log, y)
        m = dl.m_of(y)
        row = store.setdefault('y=' + y, {})
        log('R1 y = %s: m = %s, c_y = %se-4 (units 1e-4 below)' % (y, dl.fnum(m, 9), e4(cy, 5)))
        log('   a    printed(M=24)  S(M=24)   S(M=32)   S converged (+- shift)   residual change   S best upper   ratio S/c_y (converged, M=24)')
        for a in avals:
            v = vals[(y, a)]
            e, conv, unc = dl.estimate(v['legendre'], v['fourier'])
            a2 = mpmath.mpf(a) ** 2
            ent = dict(lambda_converged=e.get('converged'), lambda_uncertainty=e.get('uncertainty'), lambda_legendre_best=e.get('legendre_best'), m=dl.fnum(m, 20))
            last = v['fourier'][max(v['fourier'])]
            best_up = min(list(v['fourier'].values()) + list(v['legendre'].values()))
            ent['S_converged'] = dl.fnum(a2 * (conv - m), 8)
            ent['S_uncertainty_shift'] = dl.fnum(a2 * unc, 3)
            ent['S_residual_change_last_size_to_converged'] = dl.fnum(a2 * (last - conv), 3)
            ent['S_best_upper_bound'] = dl.fnum(a2 * (best_up - m), 8)
            ent['ratio_converged'] = mpmath.nstr(a2 * (conv - m) / cy, 5)
            ent['ratio_upper_bound'] = mpmath.nstr(a2 * (best_up - m) / cy, 5)
            if 24 in v['legendre']:
                ent['S_M24'] = dl.fnum(a2 * (v['legendre'][24] - m), 8)
                ent['ratio_M24'] = mpmath.nstr(a2 * (v['legendre'][24] - m) / cy, 5)
            if 32 in v['legendre']:
                ent['S_M32'] = dl.fnum(a2 * (v['legendre'][32] - m), 8)
                ent['ratio_M32'] = mpmath.nstr(a2 * (v['legendre'][32] - m) / cy, 5)
            if y in PRINTED_T and a in A_TABLE:
                pr = PRINTED_T[y][A_TABLE.index(a)]
                ent['printed'] = pr
                ent['M24_agrees_with_printed'] = bool(24 in v['legendre'] and '%.3f' % float(a2 * (v['legendre'][24] - m) * 10 ** 4) == pr)
                ent['converged_three_digits'] = '%.3f' % float(a2 * (conv - m) * 10 ** 4)
            if (y, a) in PRINTED_RATIO:
                ent['ratio_printed'] = PRINTED_RATIO[(y, a)]
                ent['ratio_M24_two_digits'] = '%.2f' % float(ent.get('ratio_M24', '0'))
                ent['ratio_M32_two_digits'] = '%.2f' % float(ent.get('ratio_M32', '0'))
                ent['ratio_converged_two_digits'] = '%.2f' % float(ent['ratio_converged'])
            row['a=' + a] = ent
            log('%5s   %-8s     %-8s  %-8s  %-8s (%s)       %-9s         %-8s       %s, %s' % (
                a, ent.get('printed', '-'), e4(mpmath.mpf(ent.get('S_M24', '0')), 4), e4(mpmath.mpf(ent.get('S_M32', '0')), 4), e4(mpmath.mpf(ent['S_converged']), 4),
                e4(mpmath.mpf(ent['S_uncertainty_shift']), 2), e4(mpmath.mpf(ent['S_residual_change_last_size_to_converged']), 2), e4(mpmath.mpf(ent['S_best_upper_bound']), 4),
                ent['ratio_converged'][:6], ent.get('ratio_M24', '-')[:6]))
    stamp(res, 'table', time.time() - t0, 'double precision', legendre=legs, fourier=fous)


# ---------------------------------------------------------------------------------------------------------------------
# fit (R2)
# ---------------------------------------------------------------------------------------------------------------------

def lsq(points, deg, weights=None):
    """Least squares S = c0 + c1/a + ... + c_deg/a^deg on points (a, S); returns coefficients and residuals (S in any unit)."""
    x = 1.0 / np.array([p[0] for p in points])
    s = np.array([p[1] for p in points])
    V = np.vander(x, deg + 1, increasing=True)
    w = np.ones(len(x)) if weights is None else weights
    coef = np.linalg.lstsq(V * w[:, None], s * w, rcond=None)[0]
    return coef, s - V @ coef


def fit_block(points, cy_units, rng, draws=400):
    """points: (a, S, sigma) with S in units 1e-4. Variants over the window and the degree, and a Monte Carlo for the preferred one."""
    out = {'variants': {}}
    cs = []
    for amin in (1.5, 2.0, 3.0):
        sel = [p for p in points if p[0] >= amin - 1e-9]
        for deg in (1, 2, 3):
            if len(sel) < deg + 2:
                continue
            coef, resid = lsq([(p[0], p[1]) for p in sel], deg)
            key = 'amin=%g,degree=%d' % (amin, deg)
            ent = dict(points=len(sel), C=float(coef[0]), C_over_cy=float(coef[0] / cy_units), D=float(coef[1]), D_over_cy=float(coef[1] / cy_units),
                       max_abs_residual=float(np.abs(resid).max()))
            if deg >= 2:
                ent['E_over_cy'] = float(coef[2] / cy_units)
            if deg >= 3:
                ent['F_over_cy'] = float(coef[3] / cy_units)
            out['variants'][key] = ent
            if deg >= 2:
                cs.append(ent['C_over_cy'])
    sel = [p for p in points if p[0] >= 2.0 - 1e-9]
    if len(sel) >= 4:
        base, _ = lsq([(p[0], p[1]) for p in sel], 2)
        draws_c, draws_d = [], []
        for _ in range(draws):
            pert = [(p[0], p[1] + rng.standard_normal() * p[2]) for p in sel]
            c, _ = lsq(pert, 2)
            draws_c.append(c[0] / cy_units)
            draws_d.append(c[1] / cy_units)
        out['preferred'] = dict(model='S = C + D/a + E/a^2 on a >= 2', C_over_cy=float(base[0] / cy_units), C_over_cy_monte_carlo_std=float(np.std(draws_c)),
                                D_over_cy=float(base[1] / cy_units), D_over_cy_monte_carlo_std=float(np.std(draws_d)), E_over_cy=float(base[2] / cy_units),
                                spread_C_over_cy_over_variants_deg_2_3=[min(cs), max(cs)])
    return out


def part_fit(res, log, ys, afit, fous, procs):
    t0 = time.time()
    jobs = [(a, y, 'fourier', N, 15) for y in ys for a in afit for N in fous]
    vals = dl.collect(dl.run_jobs(dl.series_job, jobs, procs))
    store = res.setdefault('fit', {})
    rng = np.random.default_rng(20261017)
    for y in ys:
        cy = cy_of(res, log, y)
        cy_u = float(cy) * 1e4
        m = dl.m_of(y)
        pts, rows = [], {}
        for a in afit:
            e, conv, unc = dl.estimate({}, vals[(y, a)]['fourier'])
            if conv is None:
                continue
            a2 = mpmath.mpf(a) ** 2
            S = float(a2 * (conv - m)) * 1e4
            sig = float(a2 * max(unc if unc is not None else 0, mpmath.mpf('3e-10'))) * 1e4
            pts.append((float(a), S, sig))
            rows['a=' + a] = dict(lambda_converged=e['converged'], lambda_uncertainty=e.get('uncertainty'), S_units_1e4='%.6f' % S, sigma_units_1e4='%.2e' % sig, ratio='%.5f' % (S / cy_u))
        blk = fit_block(pts, cy_u, rng)
        blk['points'] = rows
        blk['c_y_units_1e4'] = '%.6f' % cy_u
        blk['fourier_sizes'] = fous
        store['y=' + y] = blk
        log('R2 y = %s, c_y = %.5fe-4: S(a)/c_y at a = %s' % (y, cy_u, ' '.join('%g:%.4f' % (p[0], p[1] / cy_u) for p in pts)))
        for key, v in blk['variants'].items():
            log('   fit %-20s points %2d  C/c_y = %.5f  D/c_y = %+.4f  max residual %.2e (units 1e-4)' % (key, v['points'], v['C_over_cy'], v['D_over_cy'], v['max_abs_residual']))
        if 'preferred' in blk:
            p = blk['preferred']
            log('   preferred (%s): C/c_y = %.5f +- %.5f (Monte Carlo), D/c_y = %+.4f +- %.4f, E/c_y = %+.3f; spread over degree 2, 3 variants %.5f to %.5f' % (
                p['model'], p['C_over_cy'], p['C_over_cy_monte_carlo_std'], p['D_over_cy'], p['D_over_cy_monte_carlo_std'], p['E_over_cy'], *p['spread_C_over_cy_over_variants_deg_2_3']))
    stamp(res, 'fit', time.time() - t0, 'double precision', fourier=fous)


# ---------------------------------------------------------------------------------------------------------------------
# trial (R3)
# ---------------------------------------------------------------------------------------------------------------------

def trial_T(a, y, dps=15, nodes=None):
    """T = Q_y(h)/||h||^2 for h = cos(pi u/2a) on [-a, a], from the prime side: W(e^(-y|x|) (h * h)) / a (||h||^2 = a)."""
    if dps <= 15:
        af = float(a)
        k = math.pi / (2 * af)
        X = 2 * af

        def g(x):
            x = float(x)
            return 0.5 * (X - x) * math.cos(k * x) + math.sin(k * x) / (2 * k)

        return mpmath.mpf(ym.weil_functional(g, X, y=y, dps=15, nodes=nodes)) / mpmath.mpf(a)
    with mpmath.workdps(dps + 15):
        am = mpmath.mpf(a)
        k = mpmath.pi / (2 * am)
        X = 2 * am

        def g(x):
            return (X - x) * mpmath.cos(k * x) / 2 + mpmath.sin(k * x) / (2 * k)

        return mpmath.mpf(ym.weil_functional(g, X, y=y, dps=dps, nodes=nodes)) / am


def trial_T_matrix(a, y, M=24):
    """The same Rayleigh quotient through the Legendre matrix: c^T B c / c^T c with c the Legendre coefficients of h."""
    af = float(a)
    xs, ws = np.polynomial.legendre.leggauss(80)
    pj = np.array([math.sqrt((2 * j + 1) / (2 * af)) * np.polynomial.legendre.Legendre.basis(j)(xs) for j in range(M)])
    c = (pj * ws * af).dot(np.cos(math.pi * xs / 2))
    B = ym.build_matrix(a, y=y, basis='legendre', size=M)
    return float(c @ B @ c / (c @ c))


def part_trial(res, log, ys, atrial, quick):
    t0 = time.time()
    store = res.setdefault('trial', {})
    chk = {}
    for a, y in (('1', '0.5'), ('4', '0.2')):
        t15, tm = trial_T(a, y), trial_T_matrix(a, y)
        t30 = trial_T(a, y, dps=30)
        t15n = trial_T(a, y, nodes=600)
        chk['a=%s,y=%s' % (a, y)] = dict(prime_side_double=dl.fnum(t15, 14), prime_side_30_digits=dl.fnum(t30, 20), legendre_matrix_route=dl.fnum(tm, 14),
                                         difference_double_30=dl.fnum(t15 - t30, 3), difference_matrix_route=dl.fnum(t15 - tm, 3), difference_more_nodes=dl.fnum(t15 - t15n, 3))
        log('R3 check a = %s y = %s: prime side %s (30 digits %s, difference %s), matrix route difference %s, more nodes difference %s' % (
            a, y, dl.fnum(t15, 13), dl.fnum(t30, 13), dl.fnum(t15 - t30, 3), dl.fnum(t15 - tm, 3), dl.fnum(t15 - t15n, 3)))
        if quick:
            break
    store['check'] = chk
    rng = np.random.default_rng(7)
    for y in ys:
        cy = cy_of(res, log, y)
        cy_u = float(cy) * 1e4
        m = dl.m_of(y)
        pts, rows = [], {}
        for a in atrial:
            T = trial_T(a, y)
            a2 = mpmath.mpf(a) ** 2
            S = a2 * (T - m)
            rho = S - cy
            ent = dict(T=dl.fnum(T, 14), S_trial_units_1e4=e4(S, 8), remainder_units_1e4=e4(rho, 6), remainder_over_cy='%.6f' % float(rho / cy), a_times_remainder_over_cy='%.5f' % float(mpmath.mpf(a) * rho / cy),
                       ratio_to_cy='%.6f' % float(S / cy))
            for src in ('table', 'fit'):
                sub = res.get(src, {}).get('y=' + y, {})
                got = sub.get('a=' + a) or sub.get('points', {}).get('a=' + a)
                if got:
                    sc = got.get('S_converged')
                    sc = mpmath.mpf(sc) if sc is not None else mpmath.mpf(got['S_units_1e4']) * mpmath.mpf('1e-4')
                    ent['S_converged_units_1e4'] = e4(sc, 8)
                    ent['trial_minus_converged_units_1e4'] = e4(S - sc, 5)
                    ent['trial_over_converged_excess'] = '%.5f' % float(S / sc)
                    break
            rows['a=' + a] = ent
            pts.append((float(a), float(S) * 1e4, 1e-6))
        blk = fit_block(pts, cy_u, rng, draws=50)
        store['y=' + y] = dict(points=rows, fit=blk)
        log('R3 y = %s: trial S/c_y at a = %s' % (y, ' '.join('%s:%.4f' % (a, float(rows['a=' + a]['ratio_to_cy'])) for a in atrial)))
        log('   a * remainder / c_y at a = %s' % ' '.join('%s:%.3f' % (a, float(rows['a=' + a]['a_times_remainder_over_cy'])) for a in atrial))
        if 'preferred' in blk.keys():
            p = blk['preferred']
            log('   trial fit on a >= 2 (informational, the trial remainder is not a clean series in 1/a for y < 1/2): C/c_y = %.5f, D/c_y = %+.4f, E/c_y = %+.3f; variants %.5f to %.5f' % (
                p['C_over_cy'], p['D_over_cy'], p['E_over_cy'], *p['spread_C_over_cy_over_variants_deg_2_3']))
    stamp(res, 'trial', time.time() - t0, 'double precision (30 digits for the check)')


# ---------------------------------------------------------------------------------------------------------------------
# lower (R4)
# ---------------------------------------------------------------------------------------------------------------------

_GL = {}
KAPPA_STEP = 0.005
KAPPA_MAX = 16.0
_OML = []


def prolate_top(kappa, n=64):
    """Largest eigenvalue Lambda_0(kappa) of the sinc kernel sin(kappa (t - s))/(pi (t - s)) on [-1, 1], by n point Gauss Legendre."""
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    t, w = _GL[n]
    d = t[:, None] - t[None, :]
    eye = np.eye(n, dtype=bool)
    safe = np.where(eye, 1.0, d)
    K = np.where(eye, kappa / math.pi, np.sin(kappa * d) / (math.pi * safe))
    sw = np.sqrt(w)
    return float(np.linalg.eigvalsh(sw[:, None] * K * sw[None, :])[-1])


def one_minus_L0_table():
    if not _OML:
        ks = np.arange(0, int(round(KAPPA_MAX / KAPPA_STEP)) + 1) * KAPPA_STEP
        vals = [1.0] + [max(1.0 - prolate_top(float(k)), 0.0) for k in ks[1:]]
        _OML.extend(vals)
    return _OML


def omL(kappa):
    """A lower bound for 1 - Lambda_0(kappa): the table value at the next grid point to the right (1 - Lambda_0 decreases), 0 beyond KAPPA_MAX."""
    tab = one_minus_L0_table()
    if kappa >= KAPPA_MAX:
        return 0.0
    return tab[min(int(math.ceil(kappa / KAPPA_STEP - 1e-12)), len(tab) - 1)]


def prolate_checks():
    out = {}
    out['Lambda_0(0.2)'] = prolate_top(0.2)
    out['small_kappa_series_2k_over_pi_times_(1-k^2/9)'] = 2 * 0.2 / math.pi * (1 - 0.04 / 9)
    for k in (0.5, 1.0, 2.0, 4.0):
        out['Lambda_0(%g) n=64' % k] = prolate_top(k, 64)
        out['Lambda_0(%g) n=96' % k] = prolate_top(k, 96)
    out['1-Lambda_0(8)'] = 1 - prolate_top(8.0)
    out['asymptotic 4 sqrt(pi k) exp(-2k) at 8'] = 4 * math.sqrt(math.pi * 8) * math.exp(-16)
    return {k: '%.10g' % v for k, v in out.items()}


def part_lower(res, log, ys, avals, quick):
    t0 = time.time()
    store = res.setdefault('lower', {})
    ks = np.arange(1, int(round(KAPPA_MAX / KAPPA_STEP)) + 1) * KAPPA_STEP
    tab = one_minus_L0_table()
    I1 = float(np.sum(0.5 * (ks[1:] * np.array(tab[2:]) + ks[:-1] * np.array(tab[1:-1]))) * KAPPA_STEP)
    store['prolate_checks'] = prolate_checks()
    store['I1_integral_kappa_times_1_minus_Lambda0'] = '%.8f' % I1
    store['single_kappa_max_kappa2_times_1_minus_Lambda0'] = None
    kg = np.arange(1, 600) * 0.02
    kv = np.array([k * k * omL(k) for k in kg])
    ib = int(np.argmax(kv))
    store['single_kappa_max_kappa2_times_1_minus_Lambda0'] = '%.6f at kappa = %.2f' % (kv[ib], kg[ib])
    log('R4 prolate: I1 = int kappa (1 - Lambda_0) = %.6f; max kappa^2 (1 - Lambda_0) = %.6f at kappa = %.2f; asymptotic layer cake constant over pi^2/4 (times w\'\'/2 each): 8 I1/pi^2 = %.4f; single kappa: 4 max/pi^2 = %.4f' % (
        I1, kv[ib], kg[ib], 8 * I1 / math.pi ** 2, 4 * kv[ib] / math.pi ** 2))
    rgrid = np.concatenate([np.arange(1, 101) * 0.005, 0.5 + np.arange(1, 141) * 0.025, 4.0 + np.arange(1, 361) * 0.1])
    if quick:
        rgrid = np.concatenate([np.arange(1, 41) * 0.0125, 0.5 + np.arange(1, 36) * 0.1, 4.0 + np.arange(1, 37) * 1.0])
    for y in ys:
        cy = cy_of(res, log, y)
        w2 = float(mpmath.mpf(res['const']['y=' + y]['w2']))
        with mpmath.workdps(25):
            al = mpmath.mpf('0.5') + mpmath.mpf(y)
            w0 = ym.b_alpha(al, dps=25)
            psi = np.array([float(ym.symbol_exact(al, float(r), dps=25) - w0) for r in rgrid])
        ratio = psi / np.minimum(rgrid ** 2, 1.0)
        c1 = min(float(ratio.min()), w2 / 2)
        psit = np.minimum.accumulate(psi[::-1])[::-1]
        mask20 = rgrid >= 20.0
        floor20 = float(psi[mask20].min()) if mask20.any() else float('nan')
        ent = dict(w2=res['const']['y=' + y]['w2'], c1_chapter='%.8g' % c1, c1_over_half_w2='%.6f' % (c1 / (w2 / 2)), psi_nonnegative_on_grid=bool((psi > 0).all()),
                   psi_floor_on_20_40='%.6g' % floor20, grid_points=len(rgrid), argmin_ratio='%.4f' % rgrid[int(np.argmin(ratio))])
        rows = {}
        for a in avals:
            af = float(a)
            kmax = min(af, 12.0)
            kk = np.arange(1, int(kmax / 0.02) + 1) * 0.02
            vals = np.array([k * k * omL(k) for k in kk])
            ib = int(np.argmax(vals))
            lb_ch = 2 * c1 * vals[ib] / af ** 2
            lb_lay = 0.0
            for j in range(len(rgrid) - 1):
                lb_lay += (psit[j + 1] - psit[j]) * omL(af * rgrid[j + 1])
            lb_lay = 2 * (lb_lay + psit[0] * omL(af * rgrid[0]))
            row = dict(S_lower_chapter_units_1e4='%.5f' % (af ** 2 * lb_ch * 1e4), kappa_opt='%.2f' % kk[ib], S_lower_layer_cake_units_1e4='%.5f' % (af ** 2 * lb_lay * 1e4),
                       lower_chapter_over_cy='%.5f' % (af ** 2 * lb_ch / float(cy)), lower_layer_over_cy='%.5f' % (af ** 2 * lb_lay / float(cy)))
            for src in ('table', 'fit'):
                sub = res.get(src, {}).get('y=' + y, {})
                got = sub.get('a=' + a) or sub.get('points', {}).get('a=' + a)
                if got:
                    sc = mpmath.mpf(got['S_converged']) if 'S_converged' in got else mpmath.mpf(got['S_units_1e4']) * mpmath.mpf('1e-4')
                    row['S_converged_units_1e4'] = e4(sc, 6)
                    row['lower_layer_over_S_converged'] = '%.4f' % (af ** 2 * lb_lay / float(sc))
                    break
            tr = res.get('trial', {}).get('y=' + y, {}).get('points', {}).get('a=' + a)
            if tr:
                row['S_trial_upper_units_1e4'] = tr['S_trial_units_1e4']
            rows['a=' + a] = row
        ent['rows'] = rows
        store['y=' + y] = ent
        log('R4 y = %s: c_1 = %s (%s of w\'\'/2), psi floor on [20, 40] = %s; S bounds in units of 1e-4 (c_y = %.4f):' % (y, ent['c1_chapter'], ent['c1_over_half_w2'], ent['psi_floor_on_20_40'], float(cy) * 1e4))
        for a in avals:
            r = rows['a=' + a]
            log('   a = %-3s lower (chapter, kappa = %s) %s   lower (layer cake) %s   converged %s   trial %s' % (
                a, r['kappa_opt'], r['S_lower_chapter_units_1e4'], r['S_lower_layer_cake_units_1e4'], r.get('S_converged_units_1e4', '-'), r.get('S_trial_upper_units_1e4', '-')))
    stamp(res, 'lower', time.time() - t0, 'double precision grid, psi at 25 digits')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', default='all', help='const, table, fit, trial, lower or all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--y', default=None)
    ap.add_argument('--a', default=None)
    ap.add_argument('--afit', default=None)
    ap.add_argument('--atrial', default=None)
    ap.add_argument('--alower', default=None)
    ap.add_argument('--fourier', default=None)
    ap.add_argument('--legendre', default=None)
    args = ap.parse_args()
    t_all = time.time()
    if args.quick:
        out = args.out or '/tmp/dc_rate_quick.json'
        ys, at, af, atr, alo = ['0.5'], ['1', '2'], ['2', '3', '4', '5'], ['1', '2', '3', '4'], ['2', '4']
        fous, fou_fit, legs = [24, 48, 96, 192], [48, 96, 192], [24, 32]
    else:
        out = args.out or dl.default_path('dc_rate.json')
        ys, at, af, atr, alo = Y_RATE, A_TABLE, A_FIT, A_TRIAL, A_LOWER
        fous, fou_fit, legs = FOU_RATE, FOU_FIT, LEG_RATE
    ys = dl.parse_list(args.y) or ys
    at = dl.parse_list(args.a) or at
    af = dl.parse_list(args.afit) or af
    atr = dl.parse_list(args.atrial) or atr
    alo = dl.parse_list(args.alower) or alo
    fous_t = dl.parse_list(args.fourier, int) or fous
    fou_fit = dl.parse_list(args.fourier, int) or fou_fit
    legs = dl.parse_list(args.legendre, int) or legs
    log = dl.Log(os.path.splitext(out)[0] + '.log')
    log('=== dc_rate.py %s part %s, %s ===' % (time.strftime('%Y-%m-%d %H:%M'), args.part, 'quick' if args.quick else 'full'))
    res = dl.load(out)
    parts = ['const', 'table', 'fit', 'trial', 'lower'] if args.part == 'all' else [args.part]
    for part in parts:
        if part == 'const':
            part_const(res, log, ys)
        elif part == 'table':
            part_table(res, log, ys, at, legs, fous_t, args.procs)
        elif part == 'fit':
            part_fit(res, log, ys, af, fou_fit, args.procs)
        elif part == 'trial':
            part_trial(res, log, ys, atr, args.quick)
        elif part == 'lower':
            part_lower(res, log, ys, alo, args.quick)
        else:
            raise SystemExit('unknown part ' + part)
        dl.save(out, res)
    log('done in %.1f s, written to %s' % (time.time() - t_all, out))


if __name__ == '__main__':
    main()
