# maxwald_basepoint_fig.py: Figure mw:fig:basepoint of Section sec:maxwald, the Thorin density u_alpha of the basepoint law T_alpha.
# Added for the book in October 2026 from the checks of the authors' chapter, misc/reviews/nick-maxwald/ (its fig_basepoint.pdf, redrawn
# from curves computed by the book's own script). The four panels are u_alpha(lambda) = Re (xi'/xi)(alpha + i sqrt(lambda))/(2 pi sqrt(lambda)),
# computed on the line Re s = alpha with no zeros as input, for alpha = 1.5, 1, 0.7, 0.55 on 20 <= lambda <= 1750 (from
# maxwald_basepoint.json, where they are computed by maxwald_basepoint.py on Hopper); the dashed lines mark gamma_k^2. Under RH the term of
# the zero gamma_k is (lambda + gamma_k^2 + delta^2)/(2 gamma_k sqrt(lambda)) times a Cauchy density centred at gamma_k^2 - delta^2 with half
# width 2 delta gamma_k, delta = alpha - 1/2, so near its peak it is a Lorentzian of full width about 4 delta gamma_k and mass 1.
# The paperback is black and white: one black curve per panel, gray dashed marks.
# Run from a folder holding maxwald_basepoint.json: python3 maxwald_basepoint_fig.py [json [pdf]]; writes fig/maxwald_basepoint.pdf (fig/ is
# created if absent).
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})

src = sys.argv[1] if len(sys.argv) > 1 else 'maxwald_basepoint.json'
dst = sys.argv[2] if len(sys.argv) > 2 else 'fig/maxwald_basepoint.pdf'
os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)
D = json.load(open(src))['figure']
lam = np.array(D['lambda_grid'])
fig, axs = plt.subplots(4, 1, figsize=(6.0, 6.4), sharex=True)
for ax, key in zip(axs, D['u']):
    ax.plot(lam, np.array(D['u'][key]), color='black', lw=1.1)
    for g2 in D['gamma_squared']:
        ax.axvline(g2, color='0.55', lw=0.6, ls='--')
    ax.set_ylabel(r'$u_\alpha$')
    ax.text(0.98, 0.84, r'$\alpha=%s$' % key, transform=ax.transAxes, ha='right')
    ax.set_ylim(bottom=0)
axs[-1].set_xlabel(r'$\lambda$')
axs[-1].set_xlim(0, 1800)
fig.tight_layout()
fig.savefig(dst)
print('wrote', dst)
