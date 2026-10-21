# polya_fig.py: Figures pol:fig:one and pol:fig:two of Chapter 16 (Polya and the Zeros of Entire Functions), redrawn from the data files of the book's own scripts.
# Added for the book in October 2026 from the figures of the authors' chapter (misc/reviews/nick-polya/), whose curves are not used: every curve here is read from a data file or computed below.
# Figure 1, left: the clock U_{3/4}(tau) on [0, 60] with the zero count N(tau) (dotted) and the zeros of Re xi(3/4 + i tau) on the half integer levels (data/polya_clock.json, key figure1_left, written by polya_clock.py). Right: theta(t_k)/pi - (k - 3/2) for the 138 zeros of Polya's Xi* below height 300 with the asymptote -(2 pi^2 - 3)/(pi t) of Proposition pol:prop:gram (data/polya_xistar.json, key figure1_right, written by polya_xistar.py).
# Figure 2, left: -(log f_X)''(x) on a log scale, from direct quadrature on [0, 0.95] (solid) and from the residue expansion with twelve zeros on [0.95, 4] (dashed) (data/polya_density.json, key log_concavity.series, written by polya_density.py). Right: the compensated clock condition D_alpha(r) = Re (xi'/xi)(alpha + i r) - (xi'/xi)(alpha) for alpha = 1, 0.75, 0.6 on 0 < r <= 60 at step 0.02. The right panel is computed here from the Hadamard sum over the first 1700 zeros (data/zeros_1700.txt, all on the line), D_alpha(r) = sum over gamma > 0 of [f(a, r - gamma) + f(a, r + gamma) - 2 f(a, gamma)], f(a, x) = a/(a^2 + x^2), a = alpha - 1/2, in which the tail beyond the 1700th zero contributes less than 1e-9 a r^2; the 400 point series of data/polya_bessel.json (step 0.15, key figure2_right) is too coarse for the peaks at alpha = 0.6 (width about 0.2) and is used only as a check, whose largest difference is printed.
# The paperback is black and white: black curves, gray guides, line styles instead of colours.
# Run from the repository root: python3 code/polya_fig.py [--root DIR]; writes fig/polya_fig1.pdf and fig/polya_fig2.pdf (about 3 s in one process).
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})


def load(root, name):
    with open(os.path.join(root, 'data', name)) as fh:
        return json.load(fh)


def floats(xs):
    return np.array([float(x) for x in xs])


def figure_one(root, dst):
    clock = load(root, 'polya_clock.json')['figure1_left']
    tau, U = floats(clock['tau']), floats(clock['U'])
    zeros, levels = floats(clock['zeros_of_Re_xi']), floats(clock['levels'])
    ords = np.sort(floats(clock['ordinates_below_height']))
    xs = load(root, 'polya_xistar.json')['figure1_right']
    t, dev = floats(xs['t']), floats(xs['theta_over_pi_minus_k_minus_3_over_2'])
    keep = t <= 300.0
    t, dev = t[keep], dev[keep]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.6, 2.4))
    for lev in levels:
        ax1.axhline(lev, color='0.8', lw=0.5)
    ax1.plot(tau, U, color='black', label=r'clock $U_{3/4}(\tau)$')
    steps = np.concatenate(([0.0], ords, [tau[-1]]))
    counts = np.arange(len(steps) - 1)
    ax1.step(steps[:-1], counts, where='post', color='black', ls=':', lw=1.0, label=r'$N(\tau)$')
    ax1.plot(zeros, levels, 'o', color='black', ms=2.8, label=r'zeros of $\mathrm{Re}\,\xi(\frac{3}{4}+i\tau)$')
    ax1.set_xlim(0, tau[-1])
    ax1.set_ylim(-0.4, 14)
    ax1.set_xlabel(r'$\tau$')
    ax1.legend(loc='upper left', frameon=False, fontsize=6.5)
    asym = -(2 * np.pi ** 2 - 3) / (np.pi * np.linspace(8.0, 300.0, 400))
    ax2.axhline(0.0, color='0.8', lw=0.5)
    ax2.plot(np.linspace(8.0, 300.0, 400), asym, color='0.45', label=r'$-(2\pi^2-3)/(\pi t)$')
    ax2.plot(t, dev, 'o', color='black', ms=1.8, label=r'$\vartheta(t_k^*)/\pi-(k-\frac{3}{2})$')
    ax2.set_xlim(0, 310)
    ax2.set_ylim(-0.6, 0.02)
    ax2.set_xlabel(r'$t_k^*$')
    ax2.legend(loc='lower right', frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(dst)
    plt.close(fig)
    print('wrote', dst, '(N at height 60 =', clock['N_at_height'], ', zeros of Xi* below 300 plotted =', len(t), ')')


def figure_two(root, dst):
    dens = load(root, 'polya_density.json')['log_concavity']['series']
    x = floats(dens['x'])
    quad = np.array([np.nan if v is None else float(v) for v in dens['solid_quadrature_x_le_0_95']])
    res = np.array([np.nan if v is None else float(v) for v in dens['residues12']])
    with open(os.path.join(root, 'data', 'zeros_1700.txt')) as fh:
        gam = np.array([float(line.split()[0]) for line in fh if line.strip()])
    def hadamard(alpha, rr):
        a = alpha - 0.5
        f = lambda xx: a / (a * a + xx * xx)
        D = np.zeros_like(rr)
        for g in gam:
            D += f(rr - g) + f(rr + g) - 2.0 * f(g)
        return D

    r = np.arange(1, 3001) * 0.02
    series = {alpha: hadamard(alpha, r) for alpha in (1.0, 0.75, 0.6)}
    chk = load(root, 'polya_bessel.json')['figure2_right']
    rc = floats(chk['r'])
    worst = 0.0
    for alpha, key in ((1.0, '1'), (0.75, '0.75'), (0.6, '0.6')):
        ref = floats(chk['D'][key])
        worst = max(worst, float(np.max(np.abs(hadamard(alpha, rc) - ref) / np.abs(ref))))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.6, 2.4))
    ax1.semilogy(x, quad, color='black', label='direct quadrature')
    ax1.semilogy(x, res, color='black', ls='--', label='residue expansion, 12 zeros')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r"$-(\log f_X)''(x)$")
    ax1.set_xlim(0, 4)
    ax1.legend(loc='lower left', frameon=False, fontsize=6.5)
    ax2.axhline(0.0, color='0.8', lw=0.5)
    ax2.plot(r, series[1.0], color='black', lw=0.8, label=r'$\alpha=1$')
    ax2.plot(r, series[0.75], color='0.4', lw=0.8, ls='--', label=r'$\alpha=0.75$')
    ax2.plot(r, series[0.6], color='black', lw=0.9, ls=':', label=r'$\alpha=0.6$')
    ax2.set_xlabel(r'$r$')
    ax2.set_ylabel(r"$\mathrm{Re}\,\frac{\xi'}{\xi}(\alpha+ir)-\frac{\xi'}{\xi}(\alpha)$")
    ax2.set_xlim(0, 60)
    ax2.set_ylim(-0.4, 10.6)
    ax2.legend(loc='upper left', frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(dst)
    plt.close(fig)
    print('wrote', dst, '(smallest D on the grid: %.3e, 1.000 / 0.75 / 0.6: %s; largest relative difference from the 400 point series of polya_bessel.json: %.2e)'
          % (min(float(s.min()) for s in series.values()), ' / '.join('%.3e' % float(series[a].min()) for a in (1.0, 0.75, 0.6)), worst))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    args = ap.parse_args()
    os.makedirs(os.path.join(args.root, 'fig'), exist_ok=True)
    figure_one(args.root, os.path.join(args.root, 'fig', 'polya_fig1.pdf'))
    figure_two(args.root, os.path.join(args.root, 'fig', 'polya_fig2.pdf'))


if __name__ == '__main__':
    main()
