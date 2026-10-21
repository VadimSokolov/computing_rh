# tf_visibility.py: the numbers of Section tf:sec:window, subsection tf:sec:curve ("The visibility curve") of the chapter "Three Families of Generalised Gamma Convolutions and the Deconvolution Obstruction" (lines 12, 184, 237 to 246, 263 and 265 of the chapter source) and Figure tf:fig:visibility, regenerated as fig/tf_visibility.pdf from the data this script writes. Added for the book in October 2026, written from the mathematics (the chapter author's scripts were neither read nor run).
# Definitions. At a level y > 0 a quadruple 1/2 +- delta +- i gamma* is seen at the ordinate (v_y(gamma*) < 0) iff delta < delta_det(y), delta_det(y)^2 = y^2 + 2 y/(pi v0_y(gamma*)) (Proposition tf:prop:curve), where pi v0_y(gamma*) is the sum over the other zeros of y/(y^2 + d^2), d the distance of the zero's ordinate from gamma*; for zeros on the line this is Re xi'/xi(1/2 + y + i gamma*). Three lattice models for the other zeros, with gap Delta: midpoint (ordinates gamma* +- (k + 1/2) Delta, k >= 0, the quadruple inserted in a gap), replacement (gamma* +- k Delta, k >= 1, the model of Section sec:ch13:visibility: an on line zero replaced by the quadruple) and count preserving (gamma* +- (k + 1/2) Delta, k >= 1: two on line zeros replaced by the quadruple, the perturbation of Theorem thm:ch13:perturb, finding M3 of review A). Their sums are (pi/Delta) tanh(pi y/Delta), (pi/Delta) coth(pi y/Delta) - 1/y and the first minus 2 y/(y^2 + Delta^2/4); the first order constants delta_c = sqrt(2/S_0) are sqrt(2) Delta/pi, sqrt(6) Delta/pi and sqrt(2/(pi^2 - 8)) Delta. The greedy count of Corollary tf:cor:levels is z_0 = delta_c, z_{k+1} = delta_det(z_k), K the least k with z_k >= 1/2 (K = 0 when delta_c >= 1/2), with Delta = 2 pi/log(h/2 pi).
# Checks. (1) The constants, and the heights at which delta_c = 1/2 (h = 2 pi exp(4 pi delta_c/Delta)): 1.131e5 and 1799 in the chapter. (2) The closed forms of the lattice sums against direct summation (mpmath nsum, 30 digits), the limits delta_det^2 -> 2 Delta^2/pi^2 and 6 Delta^2/pi^2 as y -> 0 and delta_det(y) - y -> Delta/pi as y/Delta -> infinity, the band width delta_det(y) - y from the line to above (monotone fall from delta_c to Delta/pi, claimed in lines 12 and 265). (3) The greedy counts at h = 1e5, 1e8, 3e12, 1e22 (and 1e50, 1e100, 1e300 for the size of the O(log log h) term) in the three models, against (1/4) log(h/2 pi). (4) The actual zeros at the midpoint of the 10000th gap (zeros 10000 and 10001 from mpmath.zetazero): pi v0_y(gamma*) = Re xi'/xi(1/2 + y + i gamma*) from xi'/xi = 1/s + 1/(s - 1) - (log pi)/2 + psi(s/2)/2 + zeta'/zeta; delta_det(0.19) (chapter 0.4503), the root of delta_det(y) = 0.45 (0.1895), the largest difference from the midpoint lattice with the same gap on 60 points of [0.002, 0.5] and, refined, over the interval (chapter 0.0134, 'within 0.015'). The formula for xi'/xi is checked against the Hadamard sum over the 1700 zeros of data/zeros_1700.txt plus the integral of the smooth density beyond them, at the midpoint of gap 100 and y = 0.19 (agreement expected to about 1e-7, set by the smooth tail). The images of the quadruple at minus gamma*, which Proposition tf:prop:curve counts in v0, are about 1e-9 and are omitted, with their size printed. (5) The figure.
# Error control. The lattice formulas are exact; the numerical checks are at 30 digits. The actual zeros use mpmath at 20 digits (the zeta function at height 9878 and its derivative are evaluated by mpmath, the derivative by Euler Maclaurin summation; the Hadamard check tests the whole chain), and the root and the maximum are found to 1e-10 and 1e-6. Nothing here is a proof about xi: the hypothetical quadruple is data added to the computed v, as in the chapter, and both lattice models (and the count preserving one) are models.
# Output: data/tf_visibility.json (numbers as strings, with the curve of the figure), fig/tf_visibility.pdf; the log is data/tf_visibility.log. The figure is black and white (solid, dashed and dotted lines with markers, grey shading). Usage: python3 tf_visibility.py [--quick] [--out FILE] [--fig FILE] [--zeros FILE] [--fig-only]; --quick uses the 100th gap (zeros from data/zeros_1700.txt) and 12 points (a few seconds in one process) and writes to /tmp unless --out and --fig are given; --fig-only redraws the figure from the JSON file given by --out (default data/tf_visibility.json) without computing. The full case (the 10000th gap, about 100 evaluations of the zeta function at height 9878) takes about 8 s on one core (run locally within the 30 s limit on 17 October 2026).
import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

T0 = time.time()
KINDS = ('midpoint', 'replacement', 'count preserving')
HEIGHTS = ['1e5', '1e8', '3e12', '1e22']
EXTRA = ['1e50', '1e100', '1e300']


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def sd(x, n=20):
    return mp.nstr(x, n)


def spacing(h):
    return 2 * mp.pi / mp.log(h / (2 * mp.pi))


def s0_const(kind):
    """S_0 times Delta^2."""
    if kind == 'midpoint':
        return mp.pi ** 2
    if kind == 'replacement':
        return mp.pi ** 2 / 3
    return mp.pi ** 2 - 8


def dc_over_gap(kind):
    return mp.sqrt(2 / s0_const(kind))


def lattice_sum(kind, y, D):
    x = mp.pi * y / D
    mid = mp.pi / D * mp.tanh(x)
    if kind == 'midpoint':
        return mid
    if kind == 'replacement':
        return mp.pi / D * mp.coth(x) - 1 / y
    return mid - 2 * y / (y * y + D * D / 4)


def direct_sum(kind, y, D):
    half = mp.mpf(1) / 2
    if kind == 'midpoint':
        return 2 * mp.nsum(lambda j: y / (y * y + (j + half) ** 2 * D * D), [0, mp.inf])
    if kind == 'replacement':
        return 2 * mp.nsum(lambda k: y / (y * y + k * k * D * D), [1, mp.inf])
    return 2 * mp.nsum(lambda j: y / (y * y + (j + half) ** 2 * D * D), [1, mp.inf])


def delta_det(kind, y, D):
    return mp.sqrt(y * y + 2 * y / lattice_sum(kind, y, D))


def greedy(kind, D):
    half = mp.mpf(1) / 2
    z = dc_over_gap(kind) * D
    edges = [z]
    K = 0
    while z < half:
        z = delta_det(kind, z, D)
        edges.append(z)
        K += 1
        if K > 100000:
            raise RuntimeError('greedy count did not stop')
    return K, edges


def read_zeros(path):
    here = Path(__file__).resolve().parent
    cands = ([path] if path else []) + ['data/zeros_1700.txt', 'zeros_1700.txt', str(here / 'zeros_1700.txt'), str(here.parent / 'data' / 'zeros_1700.txt')]
    for c in cands:
        if c and os.path.exists(c):
            with open(c) as fh:
                return [mp.mpf(line.split()[0]) for line in fh if line.strip()], c
    return None, None


def xi_logder(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)


def pi_v0(y, th):
    return xi_logder(mp.mpc(mp.mpf(1) / 2 + y, th)).real


def make_figure(d, dst):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})
    fg = d['figure']
    y = np.array([float(v) for v in fg['y']])
    act = np.array([float(v) for v in fg['delta_det_actual']])
    lat = np.array([float(v) for v in fg['delta_det_lattice']])
    yc = d['actual_zeros'].get('root_delta_det_equals_0.45')
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    ax.fill_between(y, y, act, color='0.88', lw=0)
    ax.axhline(0.5, color='0.75', lw=0.5)
    ax.text(0.497, 0.508, r'$\delta=\frac{1}{2}$', ha='right', va='bottom', fontsize=6.5, color='0.35')
    ax.plot(y, y, color='0.45', ls=':', lw=1.0)
    ax.plot(y, act, color='black', ls='-', marker='o', ms=2.4, markevery=5)
    ax.plot(y, lat, color='black', ls='--', marker='s', ms=2.4, mfc='white', mew=0.6, markevery=(2, 5))
    if yc is not None:
        yc = float(yc)
        ax.plot([0, yc], [0.45, 0.45], color='0.45', ls=(0, (1, 2)), lw=0.6)
        ax.plot([yc], [0.45], marker='D', color='black', ms=3.0, ls='none')
        ax.annotate(r'$\varepsilon=%.4f$, $\delta=0.45$' % yc, xy=(yc, 0.45), xytext=(0.215, 0.405), fontsize=6.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='0.3'))
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.72)
    ax.set_xlabel(r'level $\varepsilon$')
    ax.set_ylabel(r'distance $\delta$ of the quadruple from the line')
    handles = [Patch(facecolor='0.88', edgecolor='none', label=r'seen: $\varepsilon<\delta<\delta_{\rm det}(\varepsilon)$'),
               Line2D([0], [0], color='black', ls='-', marker='o', ms=2.4, label='actual zeros'),
               Line2D([0], [0], color='black', ls='--', marker='s', ms=2.4, mfc='white', mew=0.6, label='midpoint lattice, same gap'),
               Line2D([0], [0], color='0.45', ls=':', label=r'$\delta=\varepsilon$')]
    ax.legend(handles=handles, loc='upper left', frameon=False, fontsize=6.5)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    fig.savefig(dst)
    plt.close(fig)
    log('wrote', dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--fig', default=None)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--fig-only', action='store_true')
    args = ap.parse_args()
    quick = args.quick
    out = args.out or ('/tmp/tf_visibility_quick.json' if quick else 'data/tf_visibility.json')
    fig = args.fig or ('/tmp/tf_visibility_quick.pdf' if quick else 'fig/tf_visibility.pdf')
    if args.fig_only:
        with open(out) as fh:
            make_figure(json.load(fh), fig)
        return
    mp.mp.dps = 30
    log('mode', 'quick' if quick else 'full')
    res = {'script': 'tf_visibility.py', 'mode': 'quick' if quick else 'full', 'parameters': {'dps_lattice': 30, 'dps_actual_zeros': 20}}

    # (1) constants and the heights where delta_c reaches 1/2
    consts = {}
    for kind in KINDS:
        c = dc_over_gap(kind)
        h_half = 2 * mp.pi * mp.exp(4 * mp.pi * c)
        consts[kind] = {'S0_times_Delta2': sd(s0_const(kind), 18), 'delta_c_over_Delta': sd(c, 18), 'height_where_delta_c_is_one_half': sd(h_half, 10)}
        log('%-16s S0 Delta^2 = %s  delta_c/Delta = %s  delta_c = 1/2 at h = %s' % (kind, sd(s0_const(kind), 12), sd(c, 12), sd(h_half, 8)))
    res['lattice_constants'] = consts

    # (2) closed forms against direct sums, limits, band width
    pts = [('0.05', '0.3'), ('0.19', '0.872118'), ('0.4', '0.65'), ('1.3', '0.4')]
    worst = mp.mpf(0)
    chk = []
    for kind in KINDS:
        for ys, Ds in pts:
            y, D = mp.mpf(ys), mp.mpf(Ds)
            cl, dr = lattice_sum(kind, y, D), direct_sum(kind, y, D)
            rel = abs(cl - dr) / abs(cl)
            worst = max(worst, rel)
            chk.append({'kind': kind, 'y': ys, 'Delta': Ds, 'closed_form': sd(cl, 25), 'direct_sum': sd(dr, 25), 'rel_diff': sd(rel, 3)})
    log('closed form lattice sums against nsum: largest relative difference %s over %d cases' % (sd(worst, 3), len(chk)))
    res['lattice_sums_check'] = {'largest_rel_diff': sd(worst, 3), 'cases': chk}
    D = mp.mpf(1)
    lim = {}
    with mp.workdps(60):
        # the coth form loses about 24 digits at y = 1e-12 Delta (1/y against pi^2 y/3 Delta^2), hence 60 digits here
        for kind, target in (('midpoint', 2 / mp.pi ** 2), ('replacement', 6 / mp.pi ** 2), ('count preserving', 2 / (mp.pi ** 2 - 8))):
            y = mp.mpf(10) ** -12
            v = delta_det(kind, y, D) ** 2
            lim[kind] = {'delta_det2_at_y_1e-12_Delta': sd(v, 20), 'limit_c_Delta2': sd(target, 20), 'abs_diff': sd(abs(v - target), 3), 'expected_abs_diff_about': '1e-24 (the O(y^2) term)'}
    for kind in ('midpoint', 'replacement'):
        lim[kind]['chapter_limit'] = '2 Delta^2/pi^2' if kind == 'midpoint' else '6 Delta^2/pi^2'
    big = []
    for j in (1, 2, 3, 4):
        y = mp.mpf(10) ** j * D
        row = {'y_over_Delta': '1e%d' % j}
        for kind in KINDS:
            dd = delta_det(kind, y, D) - y
            row[kind] = {'delta_det_minus_y': sd(dd, 15), 'minus_Delta_over_pi_times_y_over_Delta2': sd((dd - D / mp.pi) * y / D ** 2, 8)}
        big.append(row)
    nxt = {'midpoint': -1 / (2 * mp.pi ** 2), 'replacement': 1 / (2 * mp.pi ** 2), 'count preserving': 3 / (2 * mp.pi ** 2)}
    lim['large_y'] = {'Delta_over_pi': sd(D / mp.pi, 15), 'next_term_coefficient_of_Delta2_over_y': {k: sd(v, 8) for k, v in nxt.items()},
                      'note': 'delta_det(y) = y + Delta/pi + c Delta^2/y + O(Delta^3/y^2) with c = -1/(2 pi^2), +1/(2 pi^2), +3/(2 pi^2) in the three models', 'rows': big}
    log('limit y -> 0 (y = 1e-12 Delta), delta_det^2 against the limit:', {k: lim[k]['abs_diff'] for k in KINDS})
    log('large y: (delta_det - y - Delta/pi) y/Delta^2 at y = 1e4 Delta:', {k: big[-1][k]['minus_Delta_over_pi_times_y_over_Delta2'] for k in KINDS}, 'expected', {k: sd(v, 8) for k, v in nxt.items()})
    res['limits'] = lim
    band = {}
    for kind in KINDS:
        rows = []
        mono = True
        prev = None
        amin, ymin = None, None
        for k in range(0, 601):
            y = D * mp.mpf(10) ** (mp.mpf(-3) + mp.mpf(6) * k / 600)
            a = delta_det(kind, y, D) - y
            if prev is not None and a > prev:
                mono = False
            prev = a
            if amin is None or a < amin:
                amin, ymin = a, y
        for r in ('0.001', '0.5', '1', '2', '5', '10', '100'):
            y = mp.mpf(r) * D
            rows.append({'y_over_Delta': r, 'band_width_over_Delta': sd((delta_det(kind, y, D) - y) / D, 8)})
        band[kind] = {'band_width_over_Delta_at_y_to_0': sd(dc_over_gap(kind), 8), 'band_width_over_Delta_at_infinity_1_over_pi': sd(1 / mp.pi, 8),
                      'decreasing_on_y_over_Delta_1e-3_to_1e3': mono, 'least_band_width_over_Delta_on_grid': sd(amin / D, 8),
                      'y_over_Delta_of_least_on_grid': sd(ymin / D, 5), 'grid': '601 points, y/Delta from 1e-3 to 1e3 geometrically', 'rows': rows}
        log('band width %-16s from %s to %s (1/pi); decreasing on a grid: %s; least %s at y/Delta = %s' % (kind, sd(dc_over_gap(kind), 5), sd(1 / mp.pi, 5), mono, sd(amin / D, 5), sd(ymin / D, 4)))
    res['band_width'] = band

    # (3) greedy counts
    levels = []
    for hs in HEIGHTS + EXTRA:
        h = mp.mpf(hs)
        Dl = spacing(h)
        quarter = mp.log(h / (2 * mp.pi)) / 4
        rec = {'h': hs, 'log_h_over_2pi': sd(mp.log(h / (2 * mp.pi)), 12), 'Delta': sd(Dl, 12), 'quarter_log_h_over_2pi': sd(quarter, 8)}
        for kind in KINDS:
            K, edges = greedy(kind, Dl)
            rec[kind] = {'K': K, 'K_minus_quarter_log': sd(K - quarter, 5), 'delta_c': sd(edges[0], 10)}
            if hs in HEIGHTS:
                rec[kind]['edges'] = [sd(e, 8) for e in edges]
        levels.append(rec)
        log('h = %-6s Delta = %s (1/4) log(h/2pi) = %s | K: midpoint %d, replacement %d, count preserving %d | delta_c (replacement) = %s' % (hs, sd(Dl, 6), sd(quarter, 8), rec['midpoint']['K'], rec['replacement']['K'], rec['count preserving']['K'], sd(dc_over_gap('replacement') * Dl, 6)))
    res['levels'] = levels
    res['level_definition'] = 'z_0 = delta_c (first order test covers (0, z_0)); z_{k+1} = delta_det(z_k) (the level y = z_k covers (z_k, z_{k+1})); K = least k with z_k >= 1/2; Delta = 2 pi/log(h/2 pi); midpoint: ordinates gamma* +- (k + 1/2) Delta, k >= 0; replacement: gamma* +- k Delta, k >= 1; count preserving: gamma* +- (k + 1/2) Delta, k >= 1'

    # (4) the actual zeros
    mp.mp.dps = 20
    zl, zpath = read_zeros(args.zeros)
    n = 100 if quick else 10000
    if quick:
        if zl is None:
            raise SystemExit('--quick needs data/zeros_1700.txt (pass --zeros FILE)')
        g1, g2 = zl[n - 1], zl[n]
    else:
        g1, g2 = mp.zetazero(n).imag, mp.zetazero(n + 1).imag
    gstar, gap = (g1 + g2) / 2, g2 - g1
    log('zeros %d and %d: %s, %s; gamma* = %s, gap %s' % (n, n + 1, sd(g1, 14), sd(g2, 14), sd(gstar, 14), sd(gap, 10)))
    act = {'n': n, 'gamma_n': sd(g1, 18), 'gamma_n_plus_1': sd(g2, 18), 'gamma_star': sd(gstar, 18), 'gap': sd(gap, 12)}
    if zl is not None:
        th = (zl[99] + zl[100]) / 2
        y = mp.mpf('0.19')
        direct = pi_v0(y, th)
        S = mp.fsum(y / (y * y + (th - g) ** 2) + y / (y * y + (th + g) ** 2) for g in zl)
        gl = zl[-1]
        Tc = gl + mp.pi / mp.log(gl / (2 * mp.pi))
        tail = mp.quad(lambda g: mp.log(g / (2 * mp.pi)) / (2 * mp.pi) * (y / (y * y + (th - g) ** 2) + y / (y * y + (th + g) ** 2)), [Tc, 10 * Tc, mp.inf])
        act['hadamard_check_gap_100_y_0.19'] = {'xi_logder_real_part': sd(direct, 15), 'sum_over_1700_zeros': sd(S, 15), 'smooth_tail_integral': sd(tail, 8),
                                                 'difference': sd(direct - S - tail, 4), 'zeros_file': zpath}
        log('xi\'/xi check at the midpoint of gap 100, y = 0.19: Re xi\'/xi = %s, Hadamard sum %s + tail %s, difference %s' % (sd(direct, 12), sd(S, 12), sd(tail, 6), sd(direct - S - tail, 3)))
    else:
        log('zeros file not found: Hadamard check skipped')
    ngrid = 12 if quick else 60
    lo, hi = mp.mpf('0.002'), mp.mpf('0.5')
    ys = [lo + k * (hi - lo) / (ngrid - 1) for k in range(ngrid)]
    dd_act, dd_lat, piv = [], [], []
    for y in ys:
        pv = pi_v0(y, gstar)
        piv.append(pv)
        dd_act.append(mp.sqrt(y * y + 2 * y / pv))
        dd_lat.append(delta_det('midpoint', y, gap))
    diffs = [abs(a - b) for a, b in zip(dd_act, dd_lat)]
    imax = max(range(ngrid), key=lambda i: diffs[i])
    act['grid_points'] = ngrid
    act['max_abs_diff_from_midpoint_lattice_on_grid'] = sd(diffs[imax], 8)
    act['y_of_max_on_grid'] = sd(ys[imax], 8)
    from scipy.optimize import minimize_scalar
    fun = lambda yf: -float(abs(mp.sqrt(mp.mpf(yf) ** 2 + 2 * mp.mpf(yf) / pi_v0(mp.mpf(yf), gstar)) - delta_det('midpoint', mp.mpf(yf), gap)))
    r = minimize_scalar(fun, bounds=(float(ys[max(imax - 1, 0)]), float(ys[min(imax + 1, ngrid - 1)])), method='bounded', options={'xatol': 1e-5})
    act['max_abs_diff_refined'] = '%.8f' % (-r.fun)
    act['y_of_max_refined'] = '%.6f' % r.x
    log('largest |actual - midpoint lattice| on the %d point grid: %s at y = %s; refined %s at y = %s' % (ngrid, sd(diffs[imax], 6), sd(ys[imax], 5), act['max_abs_diff_refined'], act['y_of_max_refined']))
    if not quick:
        y = mp.mpf('0.19')
        pv = pi_v0(y, gstar)
        d019 = mp.sqrt(y * y + 2 * y / pv)
        yc = mp.findroot(lambda t: t * t + 2 * t / pi_v0(t, gstar) - mp.mpf('0.2025'), (mp.mpf('0.18'), mp.mpf('0.20')), solver='secant')
        mirror = sum((y - b) / ((y - b) ** 2 + (2 * gstar) ** 2) for b in (mp.mpf('0.45'), -mp.mpf('0.45')))
        act['pi_v0_at_0.19'] = sd(pv, 12)
        act['delta_det_at_0.19'] = sd(d019, 10)
        act['delta_det_lattice_at_0.19'] = sd(delta_det('midpoint', y, gap), 10)
        act['root_delta_det_equals_0.45'] = sd(yc, 10)
        act['mirror_images_at_minus_gamma_star_at_y_0.19_delta_0.45'] = sd(mirror, 4)
        log('pi v0(0.19) = %s, delta_det(0.19) = %s (lattice %s), root of delta_det(y) = 0.45 at y = %s; mirror terms %s' % (sd(pv, 10), sd(d019, 8), sd(delta_det('midpoint', y, gap), 8), sd(yc, 8), sd(mirror, 3)))
    res['actual_zeros'] = act
    res['figure'] = {'description': 'delta_det(y) at the midpoint of gap %d (solid) and of the midpoint lattice with the same gap (dashed)' % n,
                     'y': [sd(v, 12) for v in ys], 'pi_v0': [sd(v, 12) for v in piv], 'delta_det_actual': [sd(v, 12) for v in dd_act], 'delta_det_lattice': [sd(v, 12) for v in dd_lat]}
    res['seconds'] = round(time.time() - T0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as fh:
        json.dump(res, fh, indent=1)
    log('wrote', out)
    make_figure(res, fig)


if __name__ == '__main__':
    main()
