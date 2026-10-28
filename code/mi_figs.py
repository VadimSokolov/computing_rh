# mi_figs.py: Figures mi:fig:escape and mi:fig:path of Chapter 22 (ch/mi4.tex, ch/mi5.tex). They replace two interim figures of the first author's note.
#  mi_zeros.pdf: the 31 zeros of g_{9/2} with real part below 160 in the plane of z (from mi_gb_g92_{a,b,c,d}.json), the escape law (eq:escape)
#    with the constant C_{9/2} stored in the same files, and the lines |Im z| = 9/2, the de Bruijn strip |Im w| <= 9/4 in the scale w = z/2 of g.
#  mi_thorin_path.pdf: the first twelve Thorin atoms gamma_k as functions of alpha = 1/2 + b on the seven stored values
#    alpha in {1, 0.9, 0.8, 0.7, 0.6, 0.55, 0.5} (from mi_ladder_path.json, table).
# Reads only data files, nothing is recomputed. Run from the folder holding the json files: python3 mi_figs.py; writes the two figures to fig/ (created if absent).
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


# ---------------------------------------------------------------- 1. zeros of g_{9/2}, escape law, strip
zs, Cb = [], None
for part in 'abcd':
    d = load(f'mi_gb_g92_{part}.json')
    v = next(iter(d.values()))
    zs += [tuple(p) for p in v['zeros (Re,Im)']]
    Cb = v['C_9/2']
zs.sort()
b = 4.5


def escape(x):
    L = math.log(x / (2 * math.pi))
    return b - 0.5 + 2 / L * (math.pi * (b - 2) * x / (4 * b) - 1.5 * math.log(x) + math.log(Cb) - 0.5 * math.log(4 * math.pi))


xs = np.linspace(15, 165, 300)
fig, ax = plt.subplots(figsize=(5.2, 3.4))
ax.plot(xs, [escape(x) for x in xs], color='black', lw=1.2, ls='-', label='escape law')
ax.plot([z[0] for z in zs], [z[1] for z in zs], ls='none', marker='o', ms=4.5, mfc='white', mec='black', label=f'zeros of $g_{{9/2}}$ ({len(zs)})')
ax.axhline(b, color='black', lw=1.0, ls='--', label=r'de Bruijn strip $|\mathrm{Im}\,z|=9/2$')
ax.axhline(-b, color='black', lw=1.0, ls='--')
ax.set_xlim(0, 170)
ax.set_xlabel(r'$\mathrm{Re}\,z$')
ax.set_ylabel(r'$\mathrm{Im}\,z$ (upper half plane)')
ax.legend(frameon=False, loc='upper left', fontsize=9)
fig.tight_layout()
fig.savefig('fig/mi_zeros.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 2. the Thorin path
tab = load('mi_ladder_path.json')['path']['table']
alphas = [r['c'] for r in tab]
gam = np.array([[float(g) for g in r['gamma']] for r in tab])
fig, ax = plt.subplots(figsize=(5.2, 3.6))
styles = ['-', '--', '-.', ':']
marks = ['o', 's', '^', 'v', 'D', 'x']
for k in range(gam.shape[1]):
    ax.plot(alphas, gam[:, k], color='black', lw=1.0, ls=styles[k % 4], marker=marks[k % 6], ms=3.5, mfc='white')
ax.invert_xaxis()
ax.set_xlabel(r'$\alpha=1/2+b$')
ax.set_ylabel(r'$\gamma_k(b)$')
ax.text(0.98, 0.04, r'$k=1,\dots,12$ from bottom to top', transform=ax.transAxes, ha="right", fontsize=9, bbox=dict(fc="white", ec="none", pad=1.5))
fig.tight_layout()
fig.savefig('fig/mi_thorin_path.pdf')
plt.close(fig)
