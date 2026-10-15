# maxwald_thorin.py: the zeros as atoms of the smoothed Thorin density of the reciprocal Xi law, from xi'/xi with no zero as input, and the
# trigonometric analogue, for Section sec:maxwald (ch/maxwald2.tex): mw:sec:zeros with mw:tab:zeros, mw:fig:thorin and mw:fig:thorinM.
# Added for the book in October 2026 from the checks of the authors' chapter (misc/reviews/nick-maxwald/), recomputed from the mathematics
# (mpmath; nothing is taken from the chapter's scripts). With psi(lambda) = d/dlambda log xi(1/2 + sqrt(lambda)) = (xi'/xi)(1/2 + sqrt(lambda))/(2 sqrt(lambda))
# = sum_k 1/(lambda + gamma_k^2) and rho_eps(x) = pi^{-1} Im psi(-x - i eps) = sum_k pi^{-1} eps/((x - gamma_k^2)^2 + eps^2) (Poisson smoothing of
# the Thorin measure sum_k delta_{gamma_k^2} of T), and for M, psi_M(lambda) = d/dlambda log(sinh sqrt(pi lambda)/sqrt(pi lambda)) = sum_n 1/(lambda + pi n^2)
# = (pi/(2w)) coth(w) - 1/(2 lambda), w = sqrt(pi lambda):
#  1. rho_eps for eps = 3 on x in [50, 1700] (step 0.5, 30 digits per point, from xi'/xi = 1/s + 1/(s-1) - log(pi)/2 + psi_digamma(s/2)/2 + zeta'/zeta);
#     its local maxima, refined to 30 digits as roots of Im psi'(-x - i eps) (psi' in closed form with zeta''); the table tab:zeros: peak, sqrt(peak),
#     gamma_k (from data/zeros_1700.txt, read from the run folder or ../data), peak height against 1/(pi eps); the chapter's printed values compared.
#  2. The mass of rho_eps on [50, 1700] (trapezoid rule; and from the zeros, sum of arctan differences over the 1700 zeros), the number of
#     zeros in the window (7), the minimum of rho_eps on the window and where it is (the chapter says 5.6e-5 > 0), and the same minimum for rho^M at eps = 3.
#  3. Agreement of the grid values of rho_eps from xi'/xi with the Poisson smoothing of the 1700 zeros (largest difference).
#  4. rho^M_eps for eps = 0.3 on x in [0.5, 60] (step 0.02); its peaks near pi n^2, n = 1..4, with the heights against 1/(pi eps); psi_M against
#     sum_n 1/(lambda + pi n^2) at 30 digits.
#  5. fig/maxwald_thorin.pdf (rho_eps of T, eps = 3, dashed lines at gamma_k^2) and fig/maxwald_thorinM.pdf (rho^M_eps, eps = 0.3, dashed lines at pi n^2).
# Run from a folder holding this script and zeros_1700.txt: python3 maxwald_thorin.py [--quick]; writes maxwald_thorin.json and the figures to fig/
# (created if absent). The grid of rho_eps goes to a process pool with SLURM_CPUS_PER_TASK workers (one if unset). --quick is the laptop test
# (window [150, 500], step 1; the mass, minimum and figure of rho_eps then refer to that window). On Hopper (about 2 minutes, 8 cores):
# bash misc/tools/hopper_run.sh -c 8 -t 60 -g maxwald_thorin.json -g fig/maxwald_thorin.pdf -g fig/maxwald_thorinM.pdf code/maxwald_thorin.py data/zeros_1700.txt
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QUICK = '--quick' in sys.argv
DPS_GRID = 25
DPS = 30
EPS_T = 3
XLO, XHI = (150, 500) if QUICK else (50, 1700)
STEP_T = 1.0 if QUICK else 0.5
EPS_M = 0.3
MLO, MHI = 0.5, 60.0
STEP_M = 0.1 if QUICK else 0.02
T0 = time.time()
pi = mp.pi
HALF = mp.mpf(1) / 2
# the chapter's table (peak, sqrt(peak), gamma_k), k = 1..7, with the unit of the last printed digit
CH_PEAK = [199.790462, 441.926159, 625.542987, 925.673105, 1084.718263, 1412.720791, 1674.341571]
CH_SQRT = [14.1347254, 21.0220398, 25.0108574, 30.4248764, 32.9350613, 37.5861782, 40.9187191]
CH_GAMMA = [14.1347251, 21.0220396, 25.0108576, 30.4248761, 32.9350616, 37.5861782, 40.9187190]


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def load_zeros(n):
    here = os.path.dirname(os.path.abspath(__file__))
    for p in ('zeros_1700.txt', os.path.join(here, 'zeros_1700.txt'), os.path.join(here, '..', 'data', 'zeros_1700.txt')):
        if os.path.exists(p):
            z = [mp.mpf(l.strip()) for l in open(p) if l.strip()]
            if len(z) >= n:
                return z[:n]
    return [mp.zetazero(k).imag for k in range(1, n + 1)]


def Lfun(s):
    """(xi'/xi)(s)"""
    return 1 / s + 1 / (s - 1) - mp.log(pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)


def dLfun(s):
    z, z1, z2 = mp.zeta(s), mp.zeta(s, derivative=1), mp.zeta(s, derivative=2)
    return -1 / s ** 2 - 1 / (s - 1) ** 2 + mp.psi(1, s / 2) / 4 + (z2 * z - z1 ** 2) / z ** 2


def psi_T(lam):
    w = mp.sqrt(lam)
    return Lfun(HALF + w) / (2 * w)


def dpsi_T(lam):
    w = mp.sqrt(lam)
    s = HALF + w
    return dLfun(s) / (4 * lam) - Lfun(s) / (4 * lam * w)


def psi_M(lam):
    w = mp.sqrt(pi * lam)
    return (pi / (2 * w)) * mp.coth(w) - 1 / (2 * lam)


def _init(dps):
    mp.mp.dps = dps


def rho_T_task(x):
    return float(mp.im(psi_T(mp.mpc(-x, -EPS_T))) / mp.pi)


def rho_M(x, eps):
    return mp.im(psi_M(mp.mpc(-x, -eps))) / mp.pi


def main():
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    mp.mp.dps = DPS
    out = {'working_precision': 'mpmath %d digits (the grid of rho_eps for T at %d); numpy float64 for the zeros sum' % (DPS, DPS_GRID),
           'quick': QUICK, 'eps_T': EPS_T, 'window_T': [XLO, XHI], 'grid_step_T': STEP_T, 'eps_M': EPS_M}
    gam = load_zeros(1700)
    a = np.array([float(g) ** 2 for g in gam])
    # ------------------------------------------------------------ 1. grid of rho_eps for T from xi'/xi
    xs = np.arange(XLO, XHI + STEP_T / 2, STEP_T)
    if nproc > 1:
        with Pool(nproc, initializer=_init, initargs=(DPS_GRID,)) as pool:
            rho = np.array(pool.map(rho_T_task, [float(x) for x in xs], chunksize=8))
    else:
        mp.mp.dps = DPS_GRID
        rho = np.array([rho_T_task(float(x)) for x in xs])
        mp.mp.dps = DPS
    log('grid of rho_eps (T):', len(xs), 'points done')
    pk = [i for i in range(1, len(xs) - 1) if rho[i] > rho[i - 1] and rho[i] > rho[i + 1] and rho[i] > 0.01]
    rows = []
    for k, i in enumerate(pk):
        x0 = mp.mpf(float(xs[i]))
        xm = mp.findroot(lambda x: mp.im(dpsi_T(mp.mpc(-x, -EPS_T))), (x0, x0 + mp.mpf('0.05')), solver='secant', tol=mp.mpf(10) ** (-DPS + 8))
        height = rho_T_task(float(xm))
        kk = int(np.argmin(np.abs(a - float(xm))))
        rows.append({'k': kk + 1, 'peak': mp.nstr(xm, 15), 'sqrt_peak': mp.nstr(mp.sqrt(xm), 15), 'gamma_k': mp.nstr(gam[kk], 15),
                     'peak_minus_gamma_k^2': mp.nstr(xm - gam[kk] ** 2, 6), 'height': height, 'height_times_pi_eps': height * np.pi * EPS_T,
                     'sqrt_peak_minus_gamma_k': mp.nstr(mp.sqrt(xm) - gam[kk], 6)})
    out['table_zeros'] = rows
    if not QUICK:
        cmp_rows = []
        for r in rows:
            k = r['k'] - 1
            if k < 7:
                cmp_rows.append({'k': r['k'],
                                 'peak': {'chapter': CH_PEAK[k], 'ours': float(r['peak']), 'agree': bool(abs(float(r['peak']) - CH_PEAK[k]) <= 0.5e-6 * (1 + 1e-9))},
                                 'sqrt_peak': {'chapter': CH_SQRT[k], 'ours': float(r['sqrt_peak']), 'agree': bool(abs(float(r['sqrt_peak']) - CH_SQRT[k]) <= 0.5e-7 * (1 + 1e-9))},
                                 'gamma_k': {'chapter': CH_GAMMA[k], 'ours': float(r['gamma_k']), 'agree': bool(abs(float(r['gamma_k']) - CH_GAMMA[k]) <= 0.5e-7 * (1 + 1e-9))}})
        out['chapter_comparison_table'] = cmp_rows
        out['number_of_peaks_found'] = len(rows)
    log('peaks refined:', [r['peak'][:10] for r in rows])
    # ------------------------------------------------------------ 2. mass, minimum
    trap = getattr(np, 'trapezoid', None) or np.trapz
    mass_trap = float(trap(rho, xs))
    mass_zeros = float(np.sum((np.arctan((XHI - a) / EPS_T) - np.arctan((XLO - a) / EPS_T)) / np.pi))
    imin = int(np.argmin(rho))
    out['window'] = {'mass_of_rho_eps_trapezoid': mass_trap, 'mass_from_the_1700_zeros_(arctan differences)': mass_zeros,
                     'zeros_with_gamma^2_in_window': int(np.sum((a > XLO) & (a < XHI))), 'min_rho_eps': float(rho[imin]), 'x_at_min': float(xs[imin]),
                     'rho_eps_at_the_left_end_of_the_window': float(rho[0]), 'rho_eps_at_the_right_end_of_the_window': float(rho[-1])}
    rho_sum = np.sum(EPS_T / ((a[None, :] - xs[:, None]) ** 2 + EPS_T ** 2), axis=1) / np.pi
    out['zeros_sum_vs_xi_over_xi'] = {'max_abs_diff_of_grid_values': float(np.max(np.abs(rho_sum - rho))),
                                      'note': 'the sum over the first 1700 zeros omits atoms beyond gamma = 2197, whose Poisson tails on the window are below 1e-10'}
    xs_mp = [mp.mpf(float(x)) for x in xs]
    rhoM_T = np.array([float(rho_M(x, EPS_T)) for x in (xs_mp[::4] if not QUICK else xs_mp)])
    out['window']['min_of_rho^M_eps_(the trigonometric analogue)_at_eps_3_on_the_same_window_step_2'] = float(rhoM_T.min())
    if not QUICK:
        CHW = {'mass of rho_eps on [50,1700] (L557)': (6.95, 0.01, mass_trap), 'zeros in the window (L557)': (7, 1, out['window']['zeros_with_gamma^2_in_window']),
               'minimum of rho_eps on [50,1700] (L561)': (5.6e-5, 1e-6, float(rho[imin])), 'spike height 1/(pi eps) (L537)': (0.1061, 0.0001, 1 / (np.pi * EPS_T))}
        out['chapter_comparison_window'] = {k: {'chapter': v[0], 'unit': v[1], 'ours': v[2], 'agree': bool(abs(v[2] - v[0]) <= 0.5 * v[1] * (1 + 1e-9))} for k, v in CHW.items()}
    log('mass', mass_trap, mass_zeros, 'min', rho[imin], 'at', xs[imin])
    # ------------------------------------------------------------ 4. the trigonometric analogue
    lam = mp.mpc('7.3', '-0.4')
    out['psi_M_closed_form_vs_sum'] = {'lambda': '7.3 - 0.4i', 'closed form': mp.nstr(psi_M(lam), 25),
                                       'sum_n 1/(lambda + pi n^2)': mp.nstr(mp.nsum(lambda n: 1 / (lam + pi * n * n), [1, mp.inf]), 25)}
    xm_grid = np.arange(MLO, MHI + STEP_M / 2, STEP_M)
    mp.mp.dps = DPS_GRID
    rhoM = np.array([float(rho_M(mp.mpf(float(x)), mp.mpf(EPS_M))) for x in xm_grid])
    mp.mp.dps = DPS
    pkm = [i for i in range(1, len(xm_grid) - 1) if rhoM[i] > rhoM[i - 1] and rhoM[i] > rhoM[i + 1] and rhoM[i] > 0.05]
    rows = []
    for i in pkm:
        x0 = mp.mpf(float(xm_grid[i]))
        f = lambda x: rho_M(x, mp.mpf(EPS_M))
        xm = mp.findroot(lambda x: mp.diff(f, x), (x0, x0 + mp.mpf('0.01')), solver='secant')
        n = int(round(float(mp.sqrt(xm / pi))))
        h = float(f(xm))
        rows.append({'n': n, 'peak': mp.nstr(xm, 15), 'pi_n^2': mp.nstr(pi * n * n, 15), 'peak_minus_pi_n^2': mp.nstr(xm - pi * n * n, 6), 'height': h,
                     'height_times_pi_eps': h * np.pi * EPS_M})
    out['table_trigonometric_analogue'] = rows
    out['rho_M_min_on_grid'] = float(rhoM.min())
    # ------------------------------------------------------------ 5. figures
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
    os.makedirs('fig', exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.4, 2.7))
    ax.plot(xs, rho, color='black', lw=1.2)
    for g2 in a:
        if XLO <= g2 <= XHI:
            ax.axvline(g2, color='#777777', lw=0.7, ls='--')
    ax.set_xlabel(r'$x$ (squared ordinate scale)')
    ax.set_ylabel(r'$\rho_\varepsilon(x)$')
    ax.set_xlim(XLO, XHI)
    fig.tight_layout()
    fig.savefig('fig/maxwald_thorin.pdf')
    fig, ax = plt.subplots(figsize=(4.4, 2.7))
    ax.plot(xm_grid, rhoM, color='black', lw=1.2)
    for n in range(1, 5):
        ax.axvline(np.pi * n * n, color='#777777', lw=0.7, ls='--')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$\rho^M_\varepsilon(x)$')
    ax.set_xlim(MLO, MHI)
    fig.tight_layout()
    fig.savefig('fig/maxwald_thorinM.pdf')
    out['figures'] = {'thorin': {'file': 'fig/maxwald_thorin.pdf', 'x_range': [XLO, XHI], 'eps': EPS_T, 'dashed_lines': 'gamma_k^2 inside the window'},
                      'thorinM': {'file': 'fig/maxwald_thorinM.pdf', 'x_range': [MLO, MHI], 'eps': EPS_M, 'dashed_lines': 'pi n^2, n = 1..4'}}
    out['seconds'] = round(time.time() - T0, 1)
    with open('maxwald_thorin.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    log('done; wrote maxwald_thorin.json, fig/maxwald_thorin.pdf and fig/maxwald_thorinM.pdf')


if __name__ == '__main__':
    main()
