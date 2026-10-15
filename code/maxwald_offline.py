# maxwald_offline.py: detection of an off line zero in the z chart smoothed Thorin density rho_eps(x) = pi^{-1} Im psi(-x - i eps),
# psi(lambda) = d/dlambda log xi(1/2 + sqrt(lambda)) = sum_{Im rho > 0} 1/(lambda - lambda_rho), lambda_rho = (rho - 1/2)^2 = -alpha + i b,
# alpha = gamma^2 - (beta - 1/2)^2, b = 2 (beta - 1/2) gamma, for Section sec:maxwald (ch/maxwald.tex): mw:lem:detect (Lemma lem:detect
# of the authors' chapter, with the multiplicity m = 1 and gamma > 0 of its review) and mw:fig:offline. Added for the book in October 2026
# from the checks of the authors' chapter (misc/reviews/nick-maxwald/), recomputed from the mathematics (numpy float64 on the grids, mpmath
# at 30 digits on the line x = alpha), with the finer grid near the artificial off line zero that the review asks for.
# The model of the chapter: psi built from the first 60 zeros (the first 60 of data/zeros_1700.txt, read from the run folder or ../data);
# the model "off" replaces the term 1/(lambda + gamma_5^2) of the fifth zero by the two terms 1/(lambda - lambda_rho) + 1/(lambda - conj
# lambda_rho) of rho = 0.6 + i gamma_5 (the zero and its mirror 1 - conj rho both have positive imaginary part).
#  1. The chapter's grid (x in [1020, 1150] step 0.5, eps in [0.2, 14] step 0.1): the minimum of rho_eps for the zeros on the line, and for
#     the model off, with its location and the range of eps where rho_eps < 0.
#  2. A finer grid (x step 0.125, eps step 0.05): the same, and the figure fig/maxwald_offline.pdf (left: on the line; right: off the line;
#     grey scale of rho_eps clipped to [0, 0.05], hatched where rho_eps < 0, solid curve rho_eps = 0).
#  3. A zoom grid around (alpha, 0): x in [alpha - 8, alpha + 8] step 0.01, eps in (0, 8] step 0.005. The negative set of the pair alone
#     is the half disk (x - alpha)^2 + eps^2 < b^2 (the pair contributes (eps + b)/(d^2 + (eps + b)^2) - (b - eps)/(d^2 + (b - eps)^2),
#     d = x - alpha, which is < 0 exactly when d^2 < b^2 - eps^2); all negative points of the model lie in that half disk; the fraction of
#     the half disk where the model is negative.
#  4. Lemma lem:detect at 30 digits on the line x = alpha: pi rho_eps(alpha) from the full sum against eps (A_eps - 2/(b^2 - eps^2)) with A_eps
#     from its definition, for eps up to b (1 - 10^{-10}) and above b; A_eps against 2/b^2; A_0 = sum (gamma'^2 - alpha)^{-2}.
# Run from a folder holding this script and zeros_1700.txt: python3 maxwald_offline.py [--quick]; writes maxwald_offline.json and the
# figure to fig/ (created if absent). --quick is the laptop test (coarser zoom grid). On Hopper (about 1 minute, one core):
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g maxwald_offline.json -g fig/maxwald_offline.pdf code/maxwald_offline.py data/zeros_1700.txt
import json
import os
import sys
import time

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QUICK = '--quick' in sys.argv
DPS = 30
mp.mp.dps = DPS
NZ = 60
K0 = 4  # the fifth zero, index 4
DELTA = mp.mpf('0.1')  # beta - 1/2 for beta = 0.6
T0 = time.time()
out = {'working_precision': 'float64 on the grids; mpmath %d digits on the line x = alpha' % DPS, 'quick': QUICK, 'model_zeros': NZ,
       'moved_zero_index': K0 + 1, 'beta': '0.6'}


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


gam = load_zeros(NZ)
g5 = gam[K0]
alpha = g5 ** 2 - DELTA ** 2
b = 2 * DELTA * g5
out['gamma_5'] = mp.nstr(g5, 25)
out['gamma_5_squared'] = mp.nstr(g5 ** 2, 20)
out['alpha'] = mp.nstr(alpha, 20)
out['b'] = mp.nstr(b, 20)
fa, fb = float(alpha), float(b)
lam_on = [complex(-float(g) ** 2, 0.0) for g in gam]
lam_off = [lam for k, lam in enumerate(lam_on) if k != K0] + [complex(-fa, fb), complex(-fa, -fb)]


def rho(X, E, lams):
    """rho_eps(x) = pi^{-1} Im sum_j 1/(-x - i eps - lambda_j) on the grid (eps rows, x columns)."""
    z = -X[None, :] - 1j * E[:, None]
    tot = np.zeros(z.shape)
    for lam in lams:
        tot += np.imag(1.0 / (z - lam))
    return tot / np.pi


def pair_only(X, E):
    z = -X[None, :] - 1j * E[:, None]
    return (np.imag(1.0 / (z - complex(-fa, fb))) + np.imag(1.0 / (z - complex(-fa, -fb)))) / np.pi


def closed_off(X, E):
    """The same model from real formulas, as a check of the complex sum."""
    tot = np.zeros((len(E), len(X)))
    for k, g in enumerate(gam):
        if k == K0:
            continue
        g2 = float(g) ** 2
        tot += E[:, None] / ((g2 - X[None, :]) ** 2 + E[:, None] ** 2)
    d2 = (fa - X[None, :]) ** 2
    tot += (E[:, None] + fb) / (d2 + (E[:, None] + fb) ** 2) + (E[:, None] - fb) / (d2 + (E[:, None] - fb) ** 2)
    return tot / np.pi


def describe(X, E, H):
    i, j = np.unravel_index(H.argmin(), H.shape)
    neg_rows = np.any(H < 0, axis=1)
    d = {'min': float(H.min()), 'x_at_min': float(X[j]), 'eps_at_min': float(E[i]), 'grid_points': int(H.size),
         'x_step': float(X[1] - X[0]), 'eps_step': float(E[1] - E[0])}
    if neg_rows.any():
        d['eps_range_with_negative_values'] = [float(E[neg_rows].min()), float(E[neg_rows].max())]
        cols = np.any(H < 0, axis=0)
        d['x_range_with_negative_values'] = [float(X[cols].min()), float(X[cols].max())]
    return d


# ---------------------------------------------------------------- 1. the chapter's grid
X0 = np.linspace(1020, 1150, 261)
E0 = np.linspace(0.2, 14, 139)
H_on0, H_off0 = rho(X0, E0, lam_on), rho(X0, E0, lam_off)
out['chapter_grid'] = {'x': 'linspace(1020, 1150, 261)', 'eps': 'linspace(0.2, 14, 139)', 'on_line': describe(X0, E0, H_on0),
                       'off_line': describe(X0, E0, H_off0)}
jx = int(np.argmin(np.abs(X0 - 1084.5)))
je = int(np.argmin(np.abs(E0 - 6.4)))
out['chapter_grid']['off_line_at_x_1084.5_eps_6.4'] = float(H_off0[je, jx])
out['chapter_grid']['complex_sum_vs_real_formulas_max_abs_diff_off_line'] = float(np.max(np.abs(H_off0 - closed_off(X0, E0))))
log('chapter grid: on line min', out['chapter_grid']['on_line']['min'], ' off line min', out['chapter_grid']['off_line']['min'])

# ---------------------------------------------------------------- 2. finer grid and the figure
X1 = np.linspace(1020, 1150, 261 if QUICK else 1041)
E1 = np.linspace(0.05, 14, 70 if QUICK else 280)
H_on1, H_off1 = rho(X1, E1, lam_on), rho(X1, E1, lam_off)
out['fine_grid'] = {'on_line': describe(X1, E1, H_on1), 'off_line': describe(X1, E1, H_off1)}
log('fine grid: on line min', out['fine_grid']['on_line']['min'], ' off line min', out['fine_grid']['off_line']['min'])

plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
os.makedirs('fig', exist_ok=True)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True, constrained_layout=True)
vmax = 0.05
im = None
for ax, H, title in [(axs[0], H_on1, 'zeros on the line'), (axs[1], H_off1, r'$\gamma_5$ moved to $\beta=0.6$')]:
    im = ax.imshow(np.clip(H, 0, vmax), origin='lower', aspect='auto', cmap='Greys', vmin=0, vmax=vmax,
                   extent=[X1[0], X1[-1], E1[0], E1[-1]], interpolation='nearest')
    if H.min() < 0:
        ax.contourf(X1, E1, H, levels=[H.min() - 1, 0], colors='none', hatches=['////'])
        ax.contour(X1, E1, H, levels=[0], colors='black', linewidths=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(r'$x$')
    ax.grid(False)
axs[0].set_ylabel(r'$\varepsilon$')
axs[1].text(0.97, 0.97, r'hatched: $\rho_\varepsilon<0$', transform=axs[1].transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=0.5))
cb = fig.colorbar(im, ax=axs, shrink=0.85)
cb.set_label(r'$\rho_\varepsilon(x)$, clipped to $[0, 0.05]$')
fig.savefig('fig/maxwald_offline.pdf')
out['figure'] = {'file': 'fig/maxwald_offline.pdf', 'x_window': [1020, 1150], 'eps_window': [0.05, 14], 'clip': [0, vmax]}

# ---------------------------------------------------------------- 3. zoom: the negative set is inside the half disk of radius b
Xz = np.linspace(fa - 8, fa + 8, 321 if QUICK else 1601)
Ez = np.linspace(0.005, 8, 160 if QUICK else 1600)
Hz = rho(Xz, Ez, lam_off)
Pz = pair_only(Xz, Ez)
disk = (Xz[None, :] - fa) ** 2 + Ez[:, None] ** 2 < fb ** 2
neg = Hz < 0
out['zoom_grid'] = {'x_step': float(Xz[1] - Xz[0]), 'eps_step': float(Ez[1] - Ez[0]), 'min': float(Hz.min()),
                    'x_at_min': float(Xz[np.unravel_index(Hz.argmin(), Hz.shape)[1]]),
                    'eps_at_min': float(Ez[np.unravel_index(Hz.argmin(), Hz.shape)[0]]),
                    'negative_points': int(neg.sum()), 'negative_points_outside_half_disk': int((neg & ~disk).sum()),
                    'fraction_of_half_disk_negative': float((neg & disk).sum() / disk.sum()),
                    'pair_alone_sign_disagrees_with_half_disk_points': int(((Pz < 0) != disk).sum()),
                    'largest_eps_with_a_negative_value': float(Ez[np.any(neg, axis=1)].max()), 'b': float(fb)}
log('zoom grid: min', out['zoom_grid']['min'], 'outside half disk', out['zoom_grid']['negative_points_outside_half_disk'])

# ---------------------------------------------------------------- 4. Lemma lem:detect at 30 digits on x = alpha
lam_rho, lam_rho_bar = mp.mpc(-alpha, b), mp.mpc(-alpha, -b)
others = [-(g ** 2) for k, g in enumerate(gam) if k != K0]


def rho_alpha_full(eps):
    z = mp.mpc(-alpha, -eps)
    return mp.im(1 / (z - lam_rho) + 1 / (z - lam_rho_bar) + mp.fsum(1 / (z - l) for l in others)) / mp.pi


def A_eps(eps):
    z = mp.mpc(-alpha, -eps)
    return mp.fsum(mp.im(1 / (z - l)) for l in others) / eps


def rho_alpha_lemma(eps):
    return eps * (A_eps(eps) - 2 / (b ** 2 - eps ** 2)) / mp.pi


rows = []
elist = [('0.5', mp.mpf('0.5')), ('1', mp.mpf(1)), ('3', mp.mpf(3)), ('6', mp.mpf(6)), ('6.4', mp.mpf('6.4')), ('6.5', mp.mpf('6.5'))]
elist += [('b(1-1e-%d)' % j, b * (1 - mp.mpf(10) ** (-j))) for j in (2, 4, 6, 10)] + [('b(1+1e-2)', b * (1 + mp.mpf('0.01'))),
                                                                                       ('b(1+1e-4)', b * (1 + mp.mpf('0.0001')))]
for name, eps in elist:
    full = rho_alpha_full(eps)
    lem = rho_alpha_lemma(eps)
    rows.append({'eps': name, 'eps_value': mp.nstr(eps, 15), 'rho_eps(alpha) full sum': mp.nstr(full, 25),
                 'lemma formula': mp.nstr(lem, 25), 'rel_diff': mp.nstr(abs(full - lem) / abs(full), 3),
                 'A_eps': mp.nstr(A_eps(eps), 10), '2/b^2': mp.nstr(2 / b ** 2, 10)})
out['lemma_on_line_x_equals_alpha'] = rows
out['A_0_sum_(gamma^2-alpha)^-2'] = mp.nstr(mp.fsum(1 / (-l - alpha) ** 2 for l in others), 15)
out['A_eps_below_2_over_b2_for_all_eps_in_(0,b)_on_a_grid'] = bool(all(A_eps(b * mp.mpf(j) / 200) < 2 / b ** 2 for j in range(1, 200)))
out['seconds'] = round(time.time() - T0, 1)

# comparison with the numbers of the chapter (value, unit of its last printed digit)
CH = {'on line minimum over the window (L340)': (1.8e-5, 1e-6, out['chapter_grid']['on_line']['min']),
      'b for beta = 0.6 (L340)': (6.59, 0.01, float(b)),
      'location of the negative region, x (L340)': (1084.7, 0.1, float(alpha)),
      'gamma_5^2 (caption, L345)': (1084.7, 0.1, float(g5 ** 2)),
      'grid minimum with the model off (L340)': (-0.735, 0.001, out['chapter_grid']['off_line']['min']),
      'eps at the grid minimum (L340)': (6.4, 0.1, out['chapter_grid']['off_line']['eps_at_min']),
      'largest eps of the chapter grid with a negative value, "(0, 6.5)" (L340)': (6.5, 0.1, out['chapter_grid']['off_line']['eps_range_with_negative_values'][1])}
out['chapter_comparison'] = {k: {'chapter': v[0], 'unit': v[1], 'ours': v[2], 'agree': bool(abs(v[2] - v[0]) <= 0.5 * v[1] * (1 + 1e-9))}
                             for k, v in CH.items()}
with open('maxwald_offline.json', 'w') as f:
    json.dump(out, f, indent=1)
log('done; wrote maxwald_offline.json and fig/maxwald_offline.pdf;', json.dumps(out['chapter_comparison']))
