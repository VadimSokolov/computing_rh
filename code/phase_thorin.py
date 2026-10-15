# phase_thorin.py: the numbers and the curves of Section cp:sec:see (Table cp:tab:valleys, Proposition cp:prop:see, Figure cp:fig:thorint)
# and of Figure cp:fig:phase of Chapter 14, "The Thorin phase diagram" (Theorem cp:thm:phase, Corollary cp:cor:known, Definition cp:def:G).
# Added for the book in October 2026 from the checks of the authors' chapter, misc/reviews/nick-ch14/ (review-B on cp3 and cp4, review-A M1);
# written from the mathematics, not from the chapter's thorin_t.py or figs.py.
# The shifted time t clock at the basepoint alpha is phi_{t,alpha}(u) = xi_t(alpha)/xi_t(alpha + sqrt(u)), xi_t(s) = 8 H_t(-2i(s - 1/2)). In the chart
# theta = sqrt(z) its Thorin density is pi^{-1} v_t(alpha, theta), v_t(alpha, theta) = Re (xi_t'/xi_t)(alpha + i theta), whenever (t, alpha) lies in the
# GGC region (Theorem cp:thm:transfer, Remark cp:rem:transfer); the script computes v_t whether or not it does (at t = 0 and alpha = 0.6 that is
# unknown, so there v_t/pi is a function, not yet a density). In the lambda = theta^2 chart of Section sec:maxwald the density is v_t/(2 pi theta).
# 1. v_t(alpha, theta) at t in {0, 0.2}, alpha in {0.6, 0.75}, theta = 0, 0.25, ..., 100 (the P\'olya route, no primes and no zeros), and log10 |xi_t|.
# 2. Table cp:tab:valleys and Proposition cp:prop:see: for each of the 28 gaps (gamma_k, gamma_{k+1}), k <= 28, of the zeros of zeta (the file
#    zeros_1700.txt) the grid minimum of v_t(alpha, .) on the step 0.25 and the true minimum (golden section from the grid minimum; the table
#    and the proposition use the true minima), at the four combinations; the relative change from t = 0 to t = 0.2; the ratio alpha = 0.75 over
#    alpha = 0.6; the rows of the table as printed, with the printed values stored in this file and compared cell by cell; (i) the minimum over the grid and positivity;
#    (ii) the range of the valleys and the ratios; (iii) the median and the maximum of the 28 changes; the zeros of xi_0.2 below 100 (and of xi_0), with
#    their shifts (mean, median, minimum, maximum, for the 14 zeros below 60 and for all 29); the mean of v_t over [50, 100] against the average
#    (1/2) log(theta/2 pi) and the least ratio of v_t to it (the paragraph after the proposition); the constants (alpha - 1/2)/(pi (2 alpha + 1)).
# 3. Figure cp:fig:phase: the boundaries alpha = 1/2 + sqrt((alpha_0 - 1/2)^2 - (t - t_0)/2) of the regions carried along the flow from
#    (0, 1), (0, 7/8) and (0.2, 0.6) by Theorem cp:thm:phase(i), their ends t_0 + 2 (alpha_0 - 1/2)^2 on alpha = 1/2 (exact fractions), the level
#    sets t + 2 (alpha - 1/2)^2 = 0.1, ..., 0.4 and the point 1/2 + sqrt(0.1) where the level set 0.2 meets t = 0.
# Method. xi_t(s) = 8 int_0^oo e^{t u^2} Phi(u) cosh(2 w u) du and xi_t'(s) = 8 int_0^oo e^{t u^2} Phi(u) 2 u sinh(2 w u) du, w = s - 1/2, Phi = Phi_dBN
# of Section sec:ch21:flow (the integrand is even and analytic in |Im u| < pi/8). Trapezoid rule with step h = 1/160 on [0, 1.6], where
# Phi(1.6) < 10^{-800}; with w = a + i theta the sums are sum_j A_j cos(2 theta u_j), sum_j B_j sin(2 theta u_j) (and the same for xi') with the weights
# of a and t computed once, and cos, sin(2 theta u_j) by the rotation recurrence. The cancellation |xi_t(alpha + i theta)| ~ e^{-pi theta/4} costs
# 0.341 theta digits, so the working precision is int(0.36 theta) + 35 digits. The aliasing error of the trapezoid rule is of the size
# M e^{-d (2 pi/h - 2 theta)} for every d < pi/8 (Poisson summation, analyticity strip), so that h = 1/160 leaves more than 90 digits at theta = 100.
# Error control, all reported in the JSON: (a) a second discretisation (h = 1/125 on [0, 1.9], working precision 15 digits higher, tables rebuilt at
# 20 digits more) at every point of the grid, the largest relative difference of v (check every 4th point in the quick case); (b) at t = 0
# the closed form xi(s) = s (s - 1) pi^{-s/2} Gamma(s/2) zeta(s)/2 and xi'/xi = 1/s + 1/(s - 1) - (log pi)/2 + psi(s/2)/2 + zeta'/zeta (s)
# at alpha in {0.6, 0.75}, theta in {0, 17.5, 50, 100}; (c) xi_t' against the central difference of xi_t at t = 0.2; (d) the zeros: those of
# xi_0 against zeros_1700.txt, those of xi_t for t = 0.2, 0, -0.4, -1 and k = 1, 2, 3, 5, 10 against Table tab:ch21:flowzeros, each by a scan of
# step 0.1 on [gamma_k - 0.4, gamma_k + 0.5] (exactly one sign change required) and Newton's method; (e) xi_t(sigma) on the real axis.
# Writes data/phase_thorin.json (values as strings, 22 significant digits for v, 28 for the zeros), fig/cp_thorin_t.pdf and fig/cp_phase.pdf (line
# styles and direct labels as well as colours, so that they survive grayscale). Needs mpmath; the figures need matplotlib and numpy.
# Run from the repository root (a folder holding zeros_1700.txt or data/zeros_1700.txt):
#   timeout 600 python3 code/phase_thorin.py [--quick] [--out PATH] [--workers N] [--figdir DIR] [--no-figure]
#   timeout 600 python3 code/phase_thorin.py --figure [--out PATH] [--figdir DIR]      (redraws the figures from an existing JSON, no computation)
# The full case (about 8000 evaluations of v_t at up to 71 digits, 2.4 ms each in the quick timing, so under a minute on one core) is meant for Hopper (compute-C1-run.sh); --quick (theta <= 50 and the gaps 24 and 27, the first 8 gaps
# and 10 zeros, one process, under 30 s) writes next to --out (by default in the temporary folder) and never to data/ or fig/.
import argparse
import json
import math
import multiprocessing
import os
import statistics
import sys
import tempfile
import time
from fractions import Fraction

import mpmath as mp

HALF = Fraction(1, 2)
NDIG = 22      # significant digits of v in the JSON
NDIGZ = 28     # significant digits of the zeros
H1_INV, U1 = 160, Fraction(8, 5)       # primary discretisation: step 1/160 on [0, 1.6]
H2_INV, U2 = 125, Fraction(19, 10)     # check discretisation: step 1/125 on [0, 1.9]
TIMES = [Fraction(0), Fraction(1, 5)]
ALPHAS = [Fraction(3, 5), Fraction(3, 4)]
COMBOS = [(t, al) for al in ALPHAS for t in TIMES]    # index ci = 2 * (alpha index) + (time index)
TSTR = {Fraction(0): '0', Fraction(1, 5): '0.2', Fraction(-2, 5): '-0.4', Fraction(-1): '-1'}
ASTR = {Fraction(3, 5): '0.6', Fraction(3, 4): '0.75'}
# the printed Table cp:tab:valleys: (row label, gap index k or None for the origin, [alpha 0.6: t = 0, t = 0.2, change; alpha 0.75: same]);
# a row is the true minimum of v_t in the gap k, its label the ordinate of that minimum at t = 0 and alpha = 0.6 to two decimals
SOURCE_TABLE = [('0', None, ['0.004621', '0.004640', '+0.4', '0.011551', '0.011599', '+0.4']),
                ('17.50', 1, ['0.0230', '0.0231', '+0.4', '0.0573', '0.0575', '+0.4']),
                ('23.01', 2, ['0.0582', '0.0578', '-0.7', '0.1439', '0.1429', '-0.7']),
                ('27.67', 3, ['0.0389', '0.0393', '+0.9', '0.0968', '0.0976', '+0.9']),
                ('35.27', 5, ['0.0520', '0.0526', '+1.1', '0.1291', '0.1304', '+1.1']),
                ('42.12', 7, ['0.1539', '0.1496', '-2.8', '0.3729', '0.3627', '-2.7']),
                ('45.63', 8, ['0.0569', '0.0579', '+1.8', '0.1412', '0.1438', '+1.8']),
                ('95.26', 27, ['0.5625', '0.4829', '-14.2', '1.2474', '1.0906', '-12.6'])]
# Table tab:ch21:flowzeros: theta_k(t) for k = 1, 2, 3, 5, 10
FLOWZEROS = {Fraction(1, 5): ['14.07681', '20.94669', '24.95221', '32.88801', '49.72988'],
             Fraction(0): ['14.13473', '21.02204', '25.01086', '32.93506', '49.77383'],
             Fraction(-2, 5): ['14.25083', '21.17365', '25.12774', '33.02595', '49.85213'],
             Fraction(-1): ['14.42563', '21.40341', '25.30189', '33.15277', '49.93350']}
FLOWK = [1, 2, 3, 5, 10]
CFG_FULL = dict(theta_max=100, step_inv=4, gaps=28, extra_gaps=[], nz=29, check_every=1, golden=22, cf_theta=['0', '35/2', '50', '100'],
                fd=[('3/5', '35/2'), ('3/4', '50')], window=(50, 100), flow_extra=True)
CFG_QUICK = dict(theta_max=50, step_inv=4, gaps=8, extra_gaps=[24, 27], nz=10, check_every=4, golden=14, cf_theta=['0', '35/2', '100'],
                 fd=[('3/5', '35/2')], window=(25, 50), flow_extra=True)
ST = {}        # shared state, filled before the worker processes are forked
PH = {}        # wall seconds of each parallel phase


def mpq(fr):
    return mp.mpf(fr.numerator) / fr.denominator


def dps_for(theta):
    return int(0.36 * float(theta)) + 35


def cname(ci):
    t, al = COMBOS[ci]
    return 't=%s,alpha=%s' % (TSTR[t], ASTR[al])


def phi_dbn(u):
    """Phi_dBN(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}) at the current precision."""
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    eps = mp.mpf(10) ** (-(mp.mp.dps + 5))
    tot = mp.mpf(0)
    n = 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n * n * e5) * mp.exp(-mp.pi * n * n * e4)
        tot += term
        if n > 1 and abs(term) < eps * abs(tot):
            return tot
        n += 1


class Kernel:
    """The trapezoid nodes u_j = j h, j <= J, and Phi(u_j); the weights of one (t, a), a = alpha - 1/2, are built once and cached."""

    def __init__(self, h_inv, umax, dps):
        self.dps = dps
        with mp.workdps(dps):
            self.h = mp.mpf(1) / h_inv
            self.J = int(umax * h_inv)
            self.u = [j * self.h for j in range(self.J + 1)]
            self.phi = [phi_dbn(x) for x in self.u]
            self.phi_end = self.phi[-1]
        self.cache = {}

    def arrays(self, t, a):
        key = (t, a)
        if key not in self.cache:
            with mp.workdps(self.dps):
                tt, aa = mpq(t), mpq(a)
                rest, a0 = [], None
                for j, (x, p) in enumerate(zip(self.u, self.phi)):
                    w = mp.exp(tt * x * x) * p
                    if j == 0:
                        a0 = w / 2
                        continue
                    ch, sh = mp.cosh(2 * aa * x), mp.sinh(2 * aa * x)
                    rest.append((w * ch, w * sh, 2 * x * w * sh, 2 * x * w * ch))
            self.cache[key] = (a0, rest)
        return self.cache[key]

    def sums(self, t, a, theta, dps):
        """(S0r, S0i, S1r, S1i) with xi_t = 8h (S0r + i S0i) and xi_t' = 8h (S1r + i S1i) at s = 1/2 + a + i theta."""
        a0, rest = self.arrays(t, a)
        with mp.workdps(dps):
            ang = 2 * theta * self.h
            c1, s1 = mp.cos(ang), mp.sin(ang)
            c, s = mp.mpf(1), mp.mpf(0)
            s0r, s0i, s1r, s1i = +a0, mp.mpf(0), mp.mpf(0), mp.mpf(0)
            for aj, bj, cj, dj in rest:
                c, s = c * c1 - s * s1, s * c1 + c * s1
                s0r += aj * c
                s0i += bj * s
                s1r += cj * c
                s1i += dj * s
        return s0r, s0i, s1r, s1i

    def v(self, t, a, theta, dps):
        """(v_t(alpha, theta), log10 |xi_t(alpha + i theta)|)."""
        s0r, s0i, s1r, s1i = self.sums(t, a, theta, dps)
        with mp.workdps(dps):
            den = s0r * s0r + s0i * s0i
            val = (s1r * s0r + s1i * s0i) / den
        with mp.workdps(20):
            lx = float(mp.log10(8 * self.h) + mp.log10(den) / 2)
        return val, lx

    def xi(self, t, a, theta, dps):
        """(xi_t(s), xi_t'(s)) as mpc at s = 1/2 + a + i theta."""
        s0r, s0i, s1r, s1i = self.sums(t, a, theta, dps)
        with mp.workdps(dps):
            f = 8 * self.h
            return f * mp.mpc(s0r, s0i), f * mp.mpc(s1r, s1i)


def par_map(fn, tasks, workers):
    if workers <= 1:
        return [fn(x) for x in tasks]
    with multiprocessing.get_context('fork').Pool(workers) as pool:
        return list(pool.imap_unordered(fn, tasks, chunksize=1))


def golden_min(f, lo, hi, iters):
    gr = (mp.sqrt(5) - 1) / 2
    x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
    f1, f2 = f(x1), f(x2)
    for _ in range(iters):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - gr * (hi - lo)
            f1 = f(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + gr * (hi - lo)
            f2 = f(x2)
    return (lo + hi) / 2


# ---------------------------------------------------------------- tasks (run in the worker processes)
def task_grid(arg):
    ci, idxs = arg
    t, al = COMBOS[ci]
    a = al - HALF
    kp, kc = ST['kp'], ST['kc']
    res = []
    for i in idxs:
        th = ST['theta'][i]
        dps = dps_for(th)
        t0 = time.time()
        v, lx = kp.v(t, a, th, dps)
        dt = time.time() - t0
        rel = None
        if i % ST['check_every'] == 0:
            v2, _ = kc.v(t, a, th, dps + 15)
            with mp.workdps(30):
                rel = float(abs(v - v2) / abs(v))
        res.append((i, mp.nstr(v, NDIG), lx, rel, dt))
    return ci, res


def task_zero(arg):
    t, k = arg
    kp = ST['kp']
    g0 = ST['gam'][k - 1]
    zero = Fraction(0)
    dps_s, dps_n = int(0.36 * float(g0)) + 20, dps_for(g0) + 5
    t0 = time.time()
    for per_unit in (10, 50):                            # scan step 0.1, and 0.02 if that does not give exactly one sign change
        with mp.workdps(50):
            xs = [g0 - mp.mpf(4) / 10 + mp.mpf(i) / per_unit for i in range(int(0.9 * per_unit) + 1)]
        fs = [kp.sums(t, zero, x, dps_s)[0] for x in xs]
        idx = [i for i in range(len(xs) - 1) if fs[i] * fs[i + 1] < 0]
        if len(idx) == 1:
            break
    if len(idx) != 1:
        return t, k, None, len(idx), None, time.time() - t0
    with mp.workdps(50):
        x = (xs[idx[0]] + xs[idx[0] + 1]) / 2
    last = None
    for _ in range(14):
        s0r, _, _, s1i = kp.sums(t, zero, x, dps_n)      # F' = -S1i, Newton: x - F/F' = x + S0r/S1i
        with mp.workdps(dps_n):
            dx = s0r / s1i
            x = x + dx
            last = abs(dx)
        if last < mp.mpf(10) ** -30:
            break
    return t, k, mp.nstr(x, 45), 1, float(last), time.time() - t0


def task_min(arg):
    ci, k, i_best = arg
    t, al = COMBOS[ci]
    a = al - HALF
    kp = ST['kp']
    step = mp.mpf(1) / ST['step_inv']
    thg = ST['theta'][i_best]

    def f(x):
        return kp.v(t, a, x, int(0.36 * float(x)) + 18)[0]

    x = golden_min(f, thg - step, thg + step, ST['golden'])
    v, _ = kp.v(t, a, x, dps_for(x))
    return ci, k, mp.nstr(x, 14), mp.nstr(v, NDIG)


def task_closed(arg):
    al, th = arg
    dps = dps_for(th) + 10
    t0 = time.time()
    with mp.workdps(dps):
        s = mp.mpc(mpq(al), mpq(th))
        xi_c = s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
        ld_c = 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)
    xi_p, xip_p = ST['kp'].xi(Fraction(0), al - HALF, mpq(th), dps)
    with mp.workdps(dps):
        e1 = abs(xi_p - xi_c) / abs(xi_c)
        e2 = abs(xip_p / xi_p - ld_c) / abs(ld_c)
        out = dict(alpha=ASTR[al], theta=mp.nstr(mpq(th), 8), xi_closed=mp.nstr(xi_c, 20), xi_polya=mp.nstr(xi_p, 20),
                   rel_error_xi=float(e1), digits_xi=float(-mp.log10(e1)) if e1 > 0 else None,
                   ld_closed=mp.nstr(ld_c, 20), ld_polya=mp.nstr(xip_p / xi_p, 20),
                   rel_error_ld=float(e2), digits_ld=float(-mp.log10(e2)) if e2 > 0 else None, seconds=round(time.time() - t0, 2))
    return out


def task_fd(arg):
    al, th = arg
    t, a, eps = Fraction(1, 5), al - HALF, Fraction(1, 10 ** 15)
    dps = dps_for(th) + 30
    kp = ST['kp']
    _, xip = kp.xi(t, a, mpq(th), dps)
    xp, _ = kp.xi(t, a + eps, mpq(th), dps)
    xm, _ = kp.xi(t, a - eps, mpq(th), dps)
    with mp.workdps(dps):
        fd = (xp - xm) / (2 * mpq(eps))
        e = abs(fd - xip) / abs(xip)
        return dict(t='0.2', alpha=ASTR[al], theta=mp.nstr(mpq(th), 8), xi_prime=mp.nstr(xip, 20), central_difference=mp.nstr(fd, 20),
                    rel_error=float(e), eps=1e-15)


def task_real(arg):
    t, sigma = arg
    x, _ = ST['kp'].xi(t, sigma - HALF, mp.mpf(0), 60)
    return TSTR[t], str(sigma), mp.nstr(x.real, 14)


def run_phase(label, fn, tasks, workers):
    t0 = time.time()
    res = par_map(fn, tasks, workers)
    PH[label] = round(time.time() - t0, 2)
    print('%s: %d tasks in %.1f s' % (label, len(tasks), time.time() - t0), flush=True)
    return res


# ---------------------------------------------------------------- the phase diagram numbers (exact)
def phase_block():
    F = Fraction
    bases = [('critical_strip', F(0), F(1)), ('half_plane_7_8_conditional', F(0), F(7, 8)), ('polymath15', F(1, 5), F(3, 5))]
    out = dict(curves={}, level_sets={})
    for name, t0, a0 in bases:
        t_end = t0 + 2 * (a0 - HALF) ** 2
        dev = 0.0
        for i in range(301):
            t = float(t0) + float(t_end - t0) * i / 300
            al = 0.5 + math.sqrt(max((float(a0) - 0.5) ** 2 - (t - float(t0)) / 2, 0.0))
            dev = max(dev, abs((al - 0.5) ** 2 + t / 2 - (float(a0 - HALF) ** 2 + float(t0) / 2)))
        out['curves'][name] = dict(t0=str(t0), alpha0=str(a0), t_end=str(t_end), t_end_decimal=float(t_end), bound_on_Lambda=str(t_end),
                                   invariant=str((a0 - HALF) ** 2 + t0 / 2), max_deviation_of_invariant_on_300_samples=dev)
    with mp.workdps(30):
        for lev in ('1/10', '1/5', '3/10', '2/5'):
            lv = F(lev)
            out['level_sets'][str(float(lv))] = dict(
                t_on_axis_alpha_half=str(lv), alpha_on_axis_t_zero=mp.nstr(mp.mpf(1) / 2 + mp.sqrt(mpq(lv) / 2), 15))
        out['alpha_where_level_set_0.2_meets_t_0'] = mp.nstr(mp.mpf(1) / 2 + mp.sqrt(mp.mpf(1) / 10), 15)
        out['margin_constant_(alpha-1/2)/(pi(2alpha+1))'] = {ASTR[al]: mp.nstr(mpq(al - HALF) / (mp.pi * (2 * mpq(al) + 1)), 8) for al in ALPHAS}
    return out


# ---------------------------------------------------------------- the figures
def draw_figures(path, figdir):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    data = json.load(open(path))
    os.makedirs(figdir, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                         'grid.color': '#e6e6e6', 'grid.linewidth': 0.6, 'pdf.fonttype': 42})
    dark, warm, grey = '#14284b', '#d95f02', '0.55'
    # Figure cp:fig:thorint
    th = np.array([float(x) for x in data['grid']['theta']])
    fig, axs = plt.subplots(2, 1, figsize=(6.2, 4.8), sharex=True)
    for ax, al in zip(axs, ('0.6', '0.75')):
        for tt, kw in (('0', dict(color=dark, lw=1.6, ls='-', label=r'$t=0$')), ('0.2', dict(color=warm, lw=1.2, ls='--', label=r'$t=0.2$'))):
            vv = np.array([float(x) for x in data['grid']['combos']['t=%s,alpha=%s' % (tt, al)]['v']])
            brk = np.nonzero(np.diff(th) > 0.3)[0] + 1            # a quick run has a gap in theta: break the line there
            ax.semilogy(np.insert(th, brk, np.nan), np.insert(vv, brk, np.nan), **kw)
        for g in data['grid']['zeta_zeros_in_range']:
            if np.min(np.abs(th - float(g))) < 0.5:
                ax.axvline(float(g), color=grey, lw=0.5, ls=':')
        ax.set_ylabel(r'$\mathrm{Re}\,\xi_t^{\prime}/\xi_t(\alpha+i\theta)$')
        ax.text(0.99, 0.06, r'$\alpha=%s$' % al, transform=ax.transAxes, ha='right')
        ax.set_xlim(0, th.max())
    axs[0].legend(frameon=False, fontsize=8, loc='upper left')
    axs[1].set_xlabel(r'$\theta$')
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, 'cp_thorin_t.pdf'))
    plt.close(fig)
    # Figure cp:fig:phase
    ph = data['phase']
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    top = 1.02
    T = np.linspace(0, 0.55, 551)
    strip = ph['curves']['critical_strip']
    pm = ph['curves']['polymath15']

    def branch(c, tt):
        t0, a0 = float(Fraction(c['t0'])), float(Fraction(c['alpha0']))
        return np.where(tt >= t0, 0.5 + np.sqrt(np.maximum((a0 - 0.5) ** 2 - (tt - t0) / 2, 0)), np.inf)

    low = np.minimum(branch(strip, T), branch(pm, T))
    ax.fill_between(T, low, top, color='0.90', zorder=0, lw=0)
    for lev in ph['level_sets']:
        lv = float(lev)
        tt = np.linspace(0, lv, 200)
        ax.plot(tt, 0.5 + np.sqrt(np.maximum((lv - tt) / 2, 0)), color=grey, lw=0.8, ls=':')
        ax.text(-0.012, 0.5 + math.sqrt(lv / 2), lev, color='0.35', fontsize=7, ha='right', va='center')
    spec = [('critical_strip', dict(color=dark, lw=2.0, ls='-'), 'o', 'critical strip, flowed (de Bruijn): $\\Lambda\\leq\\frac{1}{2}$'),
            ('polymath15', dict(color=warm, lw=1.6, ls='-'), 's', 'Polymath15 certificate $(0.2,0.6)$, flowed: $\\Lambda\\leq0.22$'),
            ('half_plane_7_8_conditional', dict(color='0.25', lw=1.3, ls='--'), '^',
             r'$\mathrm{Re}\,s>\frac{7}{8}$ (announced, conditional), flowed: $\Lambda\leq\frac{9}{32}$')]
    for name, kw, mk, lab in spec:
        c = ph['curves'][name]
        t0, a0, te = float(Fraction(c['t0'])), float(Fraction(c['alpha0'])), float(Fraction(c['t_end']))
        tt = np.linspace(t0, te, 300)
        ax.plot(tt, 0.5 + np.sqrt(np.maximum((a0 - 0.5) ** 2 - (tt - t0) / 2, 0)), label=lab, **kw)
        ax.plot([t0], [a0], mk, color=kw['color'], ms=5, mfc='white' if name.startswith('half') else kw['color'])
        ax.plot([te], [0.5], 'x', color=kw['color'], ms=5)
    ax.plot([0], [0.5], '*', color='black', ms=10)
    ax.text(0.008, 0.515, 'RH', fontsize=8)
    ax.set_xlim(-0.04, 0.55)
    ax.set_ylim(0.49, top)
    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'basepoint $\alpha$')
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor='0.90', edgecolor='0.6', lw=0.5))
    labels.append(r'known part of $\mathcal{G}$ (certified)')
    handles.append(Line2D([], [], color=grey, lw=0.8, ls=':'))
    labels.append(r'level sets of $t+2(\alpha-\frac{1}{2})^2$, labelled by their value')
    ax.legend(handles, labels, frameon=False, fontsize=7, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, 'cp_phase.pdf'))
    plt.close(fig)
    print('wrote %s and %s' % (os.path.join(figdir, 'cp_thorin_t.pdf'), os.path.join(figdir, 'cp_phase.pdf')), flush=True)


def find_zeros_file():
    here = os.path.dirname(os.path.abspath(__file__))
    for p in ('zeros_1700.txt', 'data/zeros_1700.txt', '../data/zeros_1700.txt', os.path.join(here, 'zeros_1700.txt'),
              os.path.join(here, '..', 'data', 'zeros_1700.txt')):
        if os.path.exists(p):
            return p
    sys.exit('zeros_1700.txt not found (looked in ., data/, ../data/ and next to the script)')


def sgn(x):
    return '%+.1f' % x


def main():
    ap = argparse.ArgumentParser(description='Thorin density along the flow (Section cp:sec:see) and the phase diagram (Figure cp:fig:phase)')
    ap.add_argument('--quick', action='store_true', help='reduced case (theta <= 50, first 8 gaps), one process, writes next to --out')
    ap.add_argument('--out', default=None, help='JSON path (default data/phase_thorin.json; with --quick the temporary folder)')
    ap.add_argument('--workers', type=int, default=1, help='worker processes (0: SLURM_CPUS_PER_TASK or all cores)')
    ap.add_argument('--figdir', default=None, help='folder of the figures (default fig; with --quick the folder of --out)')
    ap.add_argument('--no-figure', action='store_true')
    ap.add_argument('--figure', action='store_true', help='only redraw the figures from the JSON')
    args = ap.parse_args()
    out_path = args.out or (os.path.join(tempfile.gettempdir(), 'phase_thorin_quick.json') if args.quick else 'data/phase_thorin.json')
    figdir = args.figdir or (os.path.dirname(os.path.abspath(out_path)) if args.quick else 'fig')
    if args.quick:
        for p in (os.path.dirname(os.path.abspath(out_path)), os.path.abspath(figdir)):
            if os.path.basename(p) in ('data', 'fig'):
                sys.exit('a quick run never writes to data/ or fig/: %s' % p)
    if args.figure:
        draw_figures(out_path, figdir)
        return
    workers = args.workers if args.workers > 0 else int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1))
    cfg = CFG_QUICK if args.quick else CFG_FULL
    T0 = time.time()
    mp.mp.dps = 50
    here_out = dict(script='phase_thorin.py', quick=args.quick, mpmath=mp.__version__, backend=mp.libmp.BACKEND, workers=workers)
    with mp.workdps(120):
        gam = [mp.mpf(line.strip()) for line in open(find_zeros_file()) if line.strip()]
    assert len(gam) == 1700 and abs(gam[0] - mp.mpf('14.134725141734693790457251983562')) < mp.mpf('1e-28')
    si = cfg['step_inv']
    gap_list = list(range(1, cfg['gaps'] + 1)) + list(cfg['extra_gaps'])
    active = set(range(cfg['theta_max'] * si + 1))
    for k in cfg['extra_gaps']:                  # the quick case adds the interior grid points of the gap (gamma_k, gamma_{k+1})
        active |= {i for i in range(int(gam[k - 1] * si), int(gam[k] * si) + 2) if gam[k - 1] < mp.mpf(i) / si < gam[k]}
    active = sorted(active)
    nth = active[-1] + 1
    theta = [mp.mpf(i) / si for i in range(nth)]          # exact binary fractions
    kp = Kernel(H1_INV, U1, 115)
    kc = Kernel(H2_INV, U2, 135)
    ST.update(kp=kp, kc=kc, gam=gam, theta=theta, check_every=cfg['check_every'], step_inv=cfg['step_inv'], golden=cfg['golden'])
    for ci in range(4):                                                  # warm the caches before the fork
        t, al = COMBOS[ci]
        kp.arrays(t, al - HALF)
        kc.arrays(t, al - HALF)
    for t in TIMES + [Fraction(-2, 5), Fraction(-1)]:
        kp.arrays(t, Fraction(0))
    with mp.workdps(30):
        tail = mp.log10(8 * kp.phi_end * mp.exp(mp.mpf(1) / 5 * kp.u[-1] ** 2 + 2 * kp.u[-1] / 4))
    res = dict(**here_out)
    res['params'] = dict(h_inv=H1_INV, umax=str(U1), nodes=kp.J, check_h_inv=H2_INV, check_umax=str(U2), check_nodes=kc.J,
                         precision_digits='int(0.36 theta) + 35 (check: 15 more)', digits_at_theta_100=dps_for(100),
                         log10_tail_beyond_umax=float(tail), theta_max=cfg['theta_max'], grid_step=1.0 / cfg['step_inv'],
                         grid_points=len(active), gaps=cfg['gaps'], extra_gaps=cfg['extra_gaps'], zeros=cfg['nz'], check_every=cfg['check_every'],
                         asymptotic_aliasing_exponent_at_theta_100='d (2 pi/h - 2 theta) = %.0f d, d < pi/8 = 0.3927 (h = 1/160); %.0f d (h = 1/125)' % (
                             2 * math.pi * H1_INV - 200, 2 * math.pi * H2_INV - 200))
    print('kernels built in %.1f s (nodes %d and %d, backend %s)' % (time.time() - T0, kp.J, kc.J, mp.libmp.BACKEND), flush=True)

    # ---- checks
    closed = run_phase('closed form at t = 0', task_closed, [(al, Fraction(s)) for al in ALPHAS for s in cfg['cf_theta']], workers)
    closed.sort(key=lambda r: (r['alpha'], float(r['theta'])))
    fd = run_phase('central difference at t = 0.2', task_fd, [(Fraction(a), Fraction(s)) for a, s in cfg['fd']], workers)
    reals = run_phase('xi_t on the real axis', task_real,
                      [(t, Fraction(s)) for t in TIMES for s in ('1/2', '3/5', '9/10', '1', '3/2')], workers)
    res['checks'] = dict(closed_form_t0=closed, xi_prime_central_difference_t02=fd,
                         xi_real_axis={'%s,sigma=%s' % (a, b): c for a, b, c in sorted(reals)})

    # ---- zeros of xi_t on the critical line
    zks = list(range(1, cfg['nz'] + 1)) + sorted({j for k in cfg['extra_gaps'] for j in (k, k + 1)})
    ztasks = [(t, k) for t in TIMES for k in zks]
    if cfg['flow_extra']:
        ztasks += [(t, k) for t in (Fraction(-2, 5), Fraction(-1)) for k in FLOWK]
    zres = run_phase('zeros', task_zero, ztasks, workers)
    Z = {}
    bad = []
    for t, k, x, nsc, last, sec in zres:
        if x is None:
            bad.append((TSTR[t], k, nsc))
        else:
            Z[(t, k)] = (mp.mpf(x), last)
    if bad:
        print('WARNING: the zero scan did not find exactly one sign change for (t, k) in %s' % bad, flush=True)
    nz = cfg['nz']
    ks = [k for k in zks if (Fraction(0), k) in Z and (Fraction(1, 5), k) in Z]
    z0 = [Z[(Fraction(0), k)][0] for k in ks]
    z2 = [Z[(Fraction(1, 5), k)][0] for k in ks]
    d0 = [abs(Z[(Fraction(0), k)][0] - gam[k - 1]) for k in ks]
    shift = [a - b for a, b in zip(z0, z2)]

    def stats(xs):
        xs = [float(x) for x in xs]
        return dict(n=len(xs), min=min(xs), max=max(xs), mean=statistics.mean(xs), median=statistics.median(xs))
    below60 = [j for j, k in enumerate(ks) if gam[k - 1] < 60]
    flow_cmp = []
    for t, col in FLOWZEROS.items():
        for k, printed in zip(FLOWK, col):
            if (t, k) in Z:
                flow_cmp.append(dict(t=TSTR[t], k=k, printed=printed, computed=mp.nstr(Z[(t, k)][0], 12),
                                     abs_diff=float(abs(Z[(t, k)][0] - mp.mpf(printed)))))
    res['zeros'] = dict(
        k=ks, failed_scans=bad, zeta=[mp.nstr(gam[k - 1], NDIGZ) for k in ks], xi_t0=[mp.nstr(x, NDIGZ) for x in z0],
        xi_t02=[mp.nstr(x, NDIGZ) for x in z2], shift_t0_minus_t02=[mp.nstr(x, 12) for x in shift],
        xi_t0_vs_zeros_1700_max_abs_difference=float(max(d0)), newton_last_step_max=float(max(Z[(t, k)][1] for (t, k) in Z)),
        shift_stats_all=stats(shift), shift_stats_below_60=stats([shift[k] for k in below60]), first_zero_t0=mp.nstr(z0[0], 12),
        first_zero_t02=mp.nstr(z2[0], 12), table_ch21_flowzeros=flow_cmp,
        table_ch21_flowzeros_max_abs_diff=max(r['abs_diff'] for r in flow_cmp))
    print('zeros: max |zero(xi_0) - gamma| = %.2e, first zero of xi_0.2 = %s, shifts %s' % (
        max(d0), mp.nstr(z2[0], 10), stats(shift)), flush=True)

    # ---- the grid
    chunks = [(ci, active[j:j + 6]) for ci in range(4) for j in range(0, len(active), 6)]
    gres = run_phase('grid', task_grid, chunks, workers)
    V = [[None] * nth for _ in range(4)]
    LX = [[None] * nth for _ in range(4)]
    REL, DT = [], []
    for ci, rows in gres:
        for i, vs, lx, rel, dt in rows:
            V[ci][i], LX[ci][i] = mp.mpf(vs), lx
            DT.append((float(theta[i]), dt))
            if rel is not None:
                REL.append((rel, ci, i))
    worst = max(REL)
    res['grid'] = dict(theta=[mp.nstr(theta[i], 8) for i in active], zeta_zeros_in_range=[mp.nstr(g, 14) for g in gam if g < theta[active[-1]]],
                       combos={cname(ci): dict(v=[mp.nstr(V[ci][i], NDIG) for i in active], log10_abs_xi=[round(LX[ci][i], 4) for i in active])
                               for ci in range(4)})
    buckets = {}
    for th_, dt in DT:
        buckets.setdefault(int(th_ // 25) * 25, []).append(dt)
    res['two_discretisations'] = dict(points_checked=len(REL), max_relative_difference_of_v=worst[0], worst_combo=cname(worst[1]),
                                      worst_theta=mp.nstr(theta[worst[2]], 8),
                                      digits=float(-math.log10(worst[0])) if worst[0] > 0 else None)
    res['timing'] = dict(evals_of_v_in_grid=len(DT), mean_ms_per_eval_by_theta_block={'%d-%d' % (b, b + 25): round(1000 * statistics.mean(x), 2)
                                                                                        for b, x in sorted(buckets.items())},
                         grid_seconds_cpu=round(sum(dt for _, dt in DT), 1))
    print('grid done; two discretisations differ by at most %.2e (%s, theta = %s)' % (worst[0], cname(worst[1]), mp.nstr(theta[worst[2]], 6)), flush=True)
    json.dump(res, open(out_path, 'w'), indent=1)

    # ---- valleys of the zeta gaps: grid minima, then the true minima
    K = len(gap_list)
    grid_min = {}
    for k in gap_list:
        idx = [i for i in active if gam[k - 1] < theta[i] < gam[k]]
        for ci in range(4):
            grid_min[(ci, k)] = min(idx, key=lambda i: V[ci][i])
    mres = run_phase('true minima', task_min, [(ci, k, grid_min[(ci, k)]) for ci in range(4) for k in gap_list], workers)
    RM = {(ci, k): (mp.mpf(x), mp.mpf(v)) for ci, k, x, v in mres}
    valleys = []
    for k in gap_list:
        rec = dict(k=k, gamma_k=mp.nstr(gam[k - 1], 14), gamma_k1=mp.nstr(gam[k], 14), grid={}, refined={}, change_percent={}, ratio_075_over_06={})
        for ci in range(4):
            i = grid_min[(ci, k)]
            rec['grid'][cname(ci)] = dict(theta=mp.nstr(theta[i], 8), v=mp.nstr(V[ci][i], NDIG))
            rec['refined'][cname(ci)] = dict(theta=mp.nstr(RM[(ci, k)][0], 10), v=mp.nstr(RM[(ci, k)][1], NDIG))
        for ai, al in enumerate(ALPHAS):
            c0, c2 = 2 * ai, 2 * ai + 1
            rec['change_percent']['alpha=' + ASTR[al]] = dict(
                grid=float(100 * (V[c2][grid_min[(c2, k)]] / V[c0][grid_min[(c0, k)]] - 1)),
                refined=float(100 * (RM[(c2, k)][1] / RM[(c0, k)][1] - 1)))
        for ti, t in enumerate(TIMES):
            rec['ratio_075_over_06']['t=' + TSTR[t]] = dict(grid=float(V[2 + ti][grid_min[(2 + ti, k)]] / V[ti][grid_min[(ti, k)]]),
                                                            refined=float(RM[(2 + ti, k)][1] / RM[(ti, k)][1]))
        valleys.append(rec)
    res['valleys'] = valleys

    # ---- Table cp:tab:valleys as printed, and the comparison
    rows, ncell, nagree = [], 0, 0
    for label, k, printed in SOURCE_TABLE:
        if k is not None and k not in gap_list:
            continue
        cells = []
        for ai in range(2):
            c0, c2 = 2 * ai, 2 * ai + 1
            if k is None:
                a0_, a2_, th_lab = V[c0][0], V[c2][0], '0'
            else:
                a0_, a2_ = RM[(c0, k)][1], RM[(c2, k)][1]
                th_lab = '%.2f' % float(RM[(c0, k)][0])
            nd = 6 if k is None else 4
            cells += ['%.*f' % (nd, float(a0_)), '%.*f' % (nd, float(a2_)), sgn(float(100 * (a2_ / a0_ - 1)))]
            if ai == 0:
                lab06 = th_lab
        agree = [c == p for c, p in zip(cells, printed)] + ([lab06 == label] if k is not None else [])
        ncell += len(agree)
        nagree += sum(agree)
        rows.append(dict(row=label, k=k, theta_label_alpha06=lab06, theta_label_alpha075=th_lab, printed=printed, recomputed=cells, agree=agree))
    res['table_cp_tab_valleys'] = dict(rows=rows, cells=ncell, cells_agreeing_as_printed=nagree)

    # ---- Proposition cp:prop:see
    prop = dict(i={}, ii={}, iii={})
    for ci in range(4):
        order = sorted(active, key=lambda i: V[ci][i])
        prop['i'][cname(ci)] = dict(argmin_theta=mp.nstr(theta[order[0]], 8), min_v=mp.nstr(V[ci][order[0]], NDIG), v_at_origin=mp.nstr(V[ci][0], NDIG),
                                    second_smallest_theta=mp.nstr(theta[order[1]], 8), second_smallest_v=mp.nstr(V[ci][order[1]], NDIG),
                                    minimum_at_origin=order[0] == 0, positive_on_grid=bool(V[ci][order[0]] > 0))
    for ci in range(4):
        gv = [(V[ci][grid_min[(ci, k)]], k, theta[grid_min[(ci, k)]]) for k in gap_list]
        rv = [(RM[(ci, k)][1], k, RM[(ci, k)][0]) for k in gap_list]
        lo, hi = min(gv), max(gv)
        rlo, rhi = min(rv), max(rv)
        prop['ii'][cname(ci)] = dict(grid_smallest=dict(v=mp.nstr(lo[0], 8), k=lo[1], theta=mp.nstr(lo[2], 8)),
                                     grid_largest=dict(v=mp.nstr(hi[0], 8), k=hi[1], theta=mp.nstr(hi[2], 8)),
                                     refined_smallest=dict(v=mp.nstr(rlo[0], 8), k=rlo[1], theta=mp.nstr(rlo[2], 8)),
                                     refined_largest=dict(v=mp.nstr(rhi[0], 8), k=rhi[1], theta=mp.nstr(rhi[2], 8)))
    for ti, t in enumerate(TIMES):
        rg = [r['ratio_075_over_06']['t=' + TSTR[t]]['grid'] for r in valleys]
        rr = [r['ratio_075_over_06']['t=' + TSTR[t]]['refined'] for r in valleys]
        prop['ii']['ratio_075_over_06_t=' + TSTR[t]] = dict(grid=dict(min=min(rg), max=max(rg)), refined=dict(min=min(rr), max=max(rr)))
    for al in ALPHAS:
        key = 'alpha=' + ASTR[al]
        for kind in ('grid', 'refined'):
            ch = [r['change_percent'][key][kind] for r in valleys]
            kmax = max(range(K), key=lambda j: abs(ch[j]))
            prop['iii'][key + ',' + kind] = dict(valleys=K, median_abs_change_percent=statistics.median(abs(x) for x in ch),
                                                 max_abs_change_percent=abs(ch[kmax]), max_at_k=gap_list[kmax],
                                                 max_at_theta_t0=valleys[kmax]['grid']['t=0,' + key]['theta'],
                                                 signed_change_at_max=ch[kmax], changes_percent=[round(x, 3) for x in ch])
    prop['iii']['zeros_shift'] = dict(first_zero=[res['zeros']['first_zero_t0'], res['zeros']['first_zero_t02']],
                                      all=res['zeros']['shift_stats_all'], below_60=res['zeros']['shift_stats_below_60'])
    res['proposition_cp_prop_see'] = prop

    # ---- the average and the margin
    w0, w1 = cfg['window']
    avg = {}
    for ci in range(4):
        ii = [i for i in active if w0 <= theta[i] <= w1]
        mv = sum(V[ci][i] for i in ii) / len(ii)
        mh = sum(mp.log(theta[i] / (2 * mp.pi)) / 2 for i in ii) / len(ii)
        jj = [i for i in active if theta[i] >= 10]
        rr = min((V[ci][i] / (mp.log(theta[i] / (2 * mp.pi)) / 2), i) for i in jj)
        avg[cname(ci)] = dict(window=[w0, w1], mean_v=mp.nstr(mv, 8), mean_half_log_theta_over_2pi=mp.nstr(mh, 8), ratio=mp.nstr(mv / mh, 6),
                              least_ratio_v_to_half_log_above_10=mp.nstr(rr[0], 6), at_theta=mp.nstr(theta[rr[1]], 8))
    res['average_check'] = avg
    res['phase'] = phase_block()
    res['timing']['phase_seconds'] = dict(PH)
    res['seconds'] = round(time.time() - T0, 1)
    json.dump(res, open(out_path, 'w'), indent=1)
    print('wrote %s in %.1f s' % (out_path, time.time() - T0), flush=True)
    if not args.no_figure:
        try:
            draw_figures(out_path, figdir)
        except ImportError as e:
            print('figures not drawn (%s): run  python3 phase_thorin.py --figure --out %s  where matplotlib is installed' % (e, out_path), flush=True)


if __name__ == '__main__':
    main()
