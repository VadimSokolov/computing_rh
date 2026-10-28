# ev_figs.py: Figures ev:fig:st, ev:fig:tail and ev:fig:bspline of Chapter 25 (ch/ev1.tex). They replace three figures of the source chapter, whose drawing scripts were not run for the book.
#  ev_smith_tawn.pdf: (a) the zero of the Dirichlet polynomial of the Dyck chain near gamma_1 for n = 100 to 12800 (real and imaginary part, from ev_g1_tail.json, part dyck),
#    with gamma_1 = 14.1347 and the real part 1/2 drawn; (b) the threshold bias of the second atom estimated after deflating the first, against gamma_2 (part deflate).
#    The stored data have nothing for the logistic Markov chain with exceedance clusters of Smith, Tawn and Coles, so that panel of the source is not drawn.
#  ev_tail.pdf: the hazard ratio estimate sqrt(f(x)/P(T>x)) of gamma_1 at x = 0.18, 0.19, 0.20 from Talbot inversion at 100 and 130 nodes, as the error against gamma_1 (part talbot).
#    The stored data hold these estimates only, not the density curve or its three atom expansion, so the tail curve of the source is not drawn.
#  ev_bspline.pdf: the supremum distance between the mean corrected truncation of the reciprocal Xi law and the Laplace inversion density for K = 1, 2, 5, 20, 40 (part bspline, from ev_g1_tail.json
#    and the values printed in ch/ev1.tex for K = 1, 2), with the peak density 80 as a line. The stored data hold K = 5, 20, 40 only; the curves of the Gumbel and the reciprocal Xi limits are not drawn.
# Reads only data files, nothing is recomputed. Run from the folder holding data/: python3 -I code/ev_figs.py; writes the three figures to fig/ (created if absent).
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
os.makedirs('fig', exist_ok=True)
DATA = 'data' if os.path.isdir('data') else '.'
G1 = 14.134725141734693790457
G2 = 21.022039638771554993


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


tail = load('ev_g1_tail.json')

# ---------------------------------------------------------------- 1. Dyck chain and threshold bias
dy = np.array([[float(a), float(b), float(c)] for a, b, c in tail['tail_dyck.json']['dyck']])
de = np.array([[float(a), float(b)] for a, b in tail['tail_talbot.json']['deflate']])
fig, ax = plt.subplots(1, 2, figsize=(5.2, 2.5))
ax[0].plot(dy[:, 0], dy[:, 2], 'o-', color='black', lw=1, ms=3, label=r'$\operatorname{Im}$')
ax[0].axhline(G1, color='gray', lw=0.8, ls='--')
ax2 = ax[0].twinx()
ax2.plot(dy[:, 0], dy[:, 1], 's-', color='0.5', lw=1, ms=3, label=r'$\operatorname{Re}$')
ax2.axhline(0.5, color='0.5', lw=0.8, ls=':')
ax2.spines['right'].set_visible(True)
ax2.grid(False)
ax[0].set_xscale('log')
ax[0].set_xlabel(r'$n$')
ax[0].set_ylabel(r'$\operatorname{Im}$ (dashed: $\gamma_1$)')
ax2.set_ylabel(r'$\operatorname{Re}$ (dotted: $1/2$)')
ax[0].set_title('Dyck chain zero', fontsize=9)
ax[1].plot(de[:, 0], de[:, 1], 'o-', color='black', lw=1, ms=3)
ax[1].axhline(G2, color='gray', lw=0.8, ls='--')
ax[1].set_xlabel(r'threshold $x$')
ax[1].set_ylabel('second atom estimate')
ax[1].set_title(r'bias (dashed: $\gamma_2$)', fontsize=9)
fig.tight_layout()
fig.savefig('fig/ev_smith_tawn.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 2. hazard ratio estimate of gamma_1
tl = tail['tail_talbot.json']['talbot']
fig, ax = plt.subplots(figsize=(4.6, 2.5))
for key, mk in (('M100', 'o'), ('M130', 'x')):
    xs = [float(a) for a, _ in tl[key]]
    err = [abs(float(b) - G1) for _, b in tl[key]]
    ax.semilogy(xs, err, mk, color='black', ms=5, label=f'Talbot, {key[1:]} nodes')
ax.set_xlabel(r'threshold $x$')
ax.set_ylabel(r'$|\sqrt{f/S}-\gamma_1|$')
ax.set_xticks([0.18, 0.19, 0.20])
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig('fig/ev_tail.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 3. B spline sup distances
bs = tail['tail_bspline.json']['bspline']
K = [5, 20, 40]
sup = [float(bs['K%d' % k]['sup']) for k in K]
K = [1, 2] + K
sup = [128.0, 42.0] + sup  # K = 1, 2 as printed in ch/ev1.tex (not in the stored data)
fig, ax = plt.subplots(figsize=(4.6, 2.5))
ax.loglog(K[:2], sup[:2], 's', color='0.5', ms=5, label='printed in the text')
ax.loglog(K[2:], sup[2:], 'o', color='black', ms=5, label='stored data')
ax.axhline(80, color='gray', lw=0.8, ls='--')
ax.text(1.1, 88, 'peak density 80', fontsize=8)
ax.set_xlabel(r'$K$')
ax.set_ylabel('sup distance')
ax.legend(frameon=False, fontsize=8, loc='lower left')
fig.tight_layout()
fig.savefig('fig/ev_bspline.pdf')
plt.close(fig)
print('wrote fig/ev_smith_tawn.pdf, fig/ev_tail.pdf, fig/ev_bspline.pdf')
