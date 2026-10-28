# ht_figs.py: Figures ht:fig:ifpt and ht:fig:yor of Chapter 24 (ch/ht3.tex). They replace two interim figures of the first author's note.
#  ht_ifpt.pdf: (a) the solved inverse first passage boundary b(t) of the reciprocal Xi law on [0.003, 0.12] with the line of slope -sqrt(2) gamma_1 drawn through its last node
#    (from ht_curved_A_b_h4e-05.npy, ht_curved_A.json); (b) b(t) for t from 0.003 to 0.02 against log(1/t)/(2 sqrt 2) (the asymptote of Proposition smallb);
#    (c) normalised boundary changes for off line pairs at the zeros 6, 20, 60 (ht_curved_C2.json) against the shift response K(t) (ht_curved_D.json).
#  ht_yor.pdf: (a) the minimum of the Thorin densities N_alpha'(t) = Re xi'/xi(alpha + it)/pi on [0, 50], attained at t = 0, and (b) N_alpha(50) against N(50) = 10,
#    for alpha in {1, 0.75, 0.6, 0.55} (from ht_misc_F.json, key yor). The stored data hold these summary values only, not the curves on [0, 50], so the curves are not drawn.
# Reads only data files, nothing is recomputed. Run from the folder holding data/: python3 -I code/ht_figs.py; writes the two figures to fig/ (created if absent).
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
os.makedirs('fig', exist_ok=True)
DATA = 'data' if os.path.isdir('data') else '.'


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


# ---------------------------------------------------------------- 1. inverse first passage boundary
t, b = np.load(os.path.join(DATA, 'ht_curved_A_b_h4e-05.npy'))
c = load('ht_curved_A.json')['formula']['c']  # sqrt(2) gamma_1
fig, ax = plt.subplots(1, 3, figsize=(5.2, 2.5))
ax[0].plot(t, b, color='black', lw=1.2, label=r'$b(t)$')
ax[0].plot(t, b[-1] - c * (t - t[-1]), color='black', lw=1.0, ls='--', label=r'slope $-\sqrt{2}\,\gamma_1$')
ax[0].set_xlabel('$t$')
ax[0].set_title('(a)', fontsize=10)
ax[0].legend(frameon=False, fontsize=7, loc='lower left')
m = t <= 0.02
ax[1].plot(t[m], b[m], color='black', lw=1.2, label=r'$b(t)$')
ax[1].plot(t[m], np.log(1 / t[m]) / (2 * math.sqrt(2)), color='black', lw=1.0, ls=':', label=r'$\log(1/t)/2\sqrt{2}$')
ax[1].set_xlabel('$t$')
ax[1].set_title('(b)', fontsize=10)
ax[1].legend(frameon=False, fontsize=7, loc='center right')
C2 = load('ht_curved_C2.json')['resp']['normalised']
K = load('ht_curved_D.json')['resp']['K']
marks = {'6': ('o', 'white'), '20': ('s', 'black'), '60': ('^', 'gray')}
for g, (mk, fc) in marks.items():
    ax[2].plot([p[0] for p in C2[g]], [p[1] for p in C2[g]], ls='none', marker=mk, ms=4.5, mfc=fc, mec='black', label=f'$\\gamma={g}$')
kt = sorted(K, key=float)
ax[2].plot([float(x) for x in kt], [K[x][0] for x in kt], color='black', lw=1.0, ls='-', label='$K(t)$')
ax[2].set_xlabel('$t$')
ax[2].set_title('(c)', fontsize=10)
ax[2].set_xticks([0.02, 0.05, 0.08])
ax[2].legend(frameon=False, fontsize=7, loc='lower right')
fig.tight_layout()
fig.savefig('fig/ht_ifpt.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 2. Yor's densities and counting functions
Y = load('ht_misc_F.json')['yor']
al = sorted(Y, key=float)
x = [float(a) for a in al]
fig, ax = plt.subplots(1, 2, figsize=(5.2, 2.6))
ax[0].plot(x, [Y[a]['min_density'] for a in al], color='black', marker='o', ms=4.5, mfc='white', lw=1.0)
ax[0].set_xlabel(r'$\alpha$')
ax[0].set_ylabel(r"min of $N_\alpha'$")
ax[0].set_title('(a)', fontsize=10)
ax[1].plot(x, [Y[a]['N50_simpson'] for a in al], color='black', marker='s', ms=4.5, mfc='black', lw=1.0, label=r'$N_\alpha(50)$')
ax[1].axhline(10, color='black', lw=1.0, ls='--', label='$N(50)=10$')
ax[1].set_xlabel(r'$\alpha$')
ax[1].set_title('(b)', fontsize=10)
ax[1].legend(frameon=False, fontsize=8, loc='center left', bbox_to_anchor=(0.02, 0.75))
fig.tight_layout()
fig.savefig('fig/ht_yor.pdf')
plt.close(fig)
