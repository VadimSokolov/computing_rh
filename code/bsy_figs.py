# bsy_figs.py: the six figures of Chapter 35 (ch/bsy1.tex, ch/bsy2.tex), group BSY-C. They replace the figures that the author supplied. Reads only data files, nothing is recomputed.
#  fig/bsy_vertical.pdf   (bsy:fig:vertical)  left: P(f > lambda) (black) and P(f < -lambda) (grey), f = log|F(1/2 + i tau)| under the Cauchy law P_{1/2}(t0 - .), for zeta at t0 = 0 (solid)
#                         and the Davenport and Heilbronn function at t0 = 85.7 (dashed); right: the running balance A+(L) - A-(L) - log|h(1 + i t0)| for zeta at t0 = 0 and 85.7 and for DH at
#                         t0 = 85.7, with the prediction of data/bsy_defects.json for DH (data/bsy_vertical.json from code/bsy_vertical.py).
#  fig/bsy_yks_rates.pdf  (bsy:fig:yks)       root mean square error against N for plain, YKS and YKS with control, theta = 1/2 and 7/8, 400 replications, with the slope -2 of the YKS theorem
#                         and the slope -1/2 (data/bsy_yks_r400.json, code/bsy_yks.py).
#  fig/bsy_riesz2d_maps.pdf (bsy:fig:maps)    the integrand log|f| Delta phi on the box for zeta and for DH, the zero marked by a cross (data/bsy_riesz2d.json).
#  fig/bsy_riesz2d_rates.pdf (bsy:fig:rates2d) root mean square error against N of the estimators, zeta and DH (data/bsy_riesz2d.json), with the slopes -1/2 and -3/2.
#  fig/bsy_dhlevy.pdf     (bsy:fig:dhlevy)    left: K_DH by the explicit formula and by the zero sum on 0.1 <= u <= 2 with the off line terms in absolute value; right: K_DH(u) e^(gamma_1 u) - 1
#                         on a logarithmic scale (data/bsy_levy.json).
#  fig/bsy_halfdisk.pdf   (bsy:fig:halfdisk)  root mean square error against N of plain Monte Carlo, the sorted trapezoid rule and the rule after the change of variable Phi, for the half disc identity (zeta at theta = 3/4 and 1/2, DH at theta = 1/2), dotted slopes -1/2, -2, -5/2 (data/bsy_halfdisk.json from code/bsy_halfdisk.py, group BSY-E).
# Run from the book folder: python3 -I code/bsy_figs.py [--alt DIR] (DIR is searched for a data file that is not in data/, for a test with quick runs); writes to fig/ (created if absent).
import argparse, json, math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(); ap.add_argument('--alt', default=None); ap.add_argument('--figdir', default='fig')
args = ap.parse_args()
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
os.makedirs(args.figdir, exist_ok=True)
DATA = 'data' if os.path.isdir('data') else '.'


def load(name):
    for d in (DATA, args.alt):
        if d and os.path.exists(os.path.join(d, name)):
            with open(os.path.join(d, name)) as f:
                return json.load(f)
    raise FileNotFoundError(name)


# ------------------------------------------------------------------ 1. vertical curves
V = load('bsy_vertical.json')
lam = np.array(V['lambda']); hh = 'hh=' + str(V['hh'][0])
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
for key, ls, lab in (('zeta|0.0', '-', r'$\zeta$, $t_0=0$'), ('dh|85.7', '--', r'DH, $t_0=85.7$')):
    c = V['cases'][key][hh]
    up = np.array(c['P_up']); lo = np.array(c['P_lo'])
    ax[0].semilogy(lam[up > 0], up[up > 0], color='black', ls=ls, lw=1.2, label=lab + r': $P(f>\lambda)$')
    ax[0].semilogy(lam[lo > 0], lo[lo > 0], color='0.6', ls=ls, lw=1.6, label=lab + r': $P(f<-\lambda)$')
ax[0].set_xlabel(r'$\lambda$'); ax[0].set_xlim(0, 10); ax[0].set_ylim(1e-6, 1.5); ax[0].legend(frameon=False, fontsize=5.5, loc='upper right')
ax[0].set_title('(a)', fontsize=9)
styles = {'zeta|0.0': ('black', '-', r'$\zeta$, $t_0=0$'), 'zeta|85.7': ('0.5', '-', r'$\zeta$, $t_0=85.7$'), 'dh|85.7': ('black', '--', r'DH, $t_0=85.7$')}
for key, (col, ls, lab) in styles.items():
    ax[1].plot(lam, V['cases'][key][hh]['balance'], color=col, ls=ls, lw=1.3, label=lab)
try:
    pred = [r['prediction'] for r in load('bsy_defects.json')['rows'] if str(r['function']).lower().startswith('dh') and abs(float(r['theta']) - 0.5) < 1e-9 and abs(float(r['t0']) - 85.7) < 0.1][0]
    ax[1].axhline(float(pred), color='black', lw=0.8, ls=':', label='DH prediction')
except Exception:
    pass
ax[1].set_xlabel(r'level cutoff $\Lambda$'); ax[1].set_xlim(0, 10); ax[1].legend(frameon=False, fontsize=7, loc='center right')
ax[1].set_title('(b)', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(args.figdir, 'bsy_vertical.pdf')); plt.close(fig)

# ------------------------------------------------------------------ 2. YKS rates
Y = load('bsy_yks_r400.json')
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7), sharey=True)
mk = {'plain': ('o', 'white', '-'), 'yks': ('s', 'black', '-'), 'yks_ctrl': ('^', '0.5', '-')}
nm = {'plain': 'plain', 'yks': 'YKS', 'yks_ctrl': 'YKS with control'}
for a, (key, ttl) in zip(ax, (('zeta_0.5', r'$\theta=1/2$'), ('zeta_0.875', r'$\theta=7/8$'))):
    S = Y['summary'][key]; N = np.array(S['N'], dtype=float)
    for est, (m, fc, ls) in mk.items():
        a.loglog(N, S['rms'][est], color='black', marker=m, mfc=fc, ms=4.5, lw=1.0, ls=ls, label=f"{nm[est]} ({S['slope'][est]:.2f})")
    r0 = S['rms']['yks'][0]
    a.loglog(N, r0 * (N / N[0]) ** -2.0, color='black', ls=':', lw=1.0, label='slope $-2$')
    a.loglog(N, S['rms']['plain'][0] * (N / N[0]) ** -0.5, color='0.6', ls=':', lw=1.0, label='slope $-1/2$')
    a.set_xlabel('$N$'); a.set_title(ttl, fontsize=9); a.legend(frameon=False, fontsize=6, loc='lower left')
ax[0].set_ylabel('root mean square error')
fig.tight_layout(); fig.savefig(os.path.join(args.figdir, 'bsy_yks_rates.pdf')); plt.close(fig)

# ------------------------------------------------------------------ 3. maps
R = load('bsy_riesz2d.json')
x = np.array(R['map_x']); sg = 1.0 + 0.4 * x; tt = 85.7 + 2.0 * x
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.8))
mx = max(np.abs(np.array(R['maps'][k])).max() for k in ('zeta', 'dh'))
for a, (k, ttl) in zip(ax, (('zeta', r'$\zeta$: $I(\varphi)=0$'), ('dh', r'DH: $I(\varphi)=2\pi\varphi(\rho)=%.3f$' % R['exact']['dh']))):
    G = np.array(R['maps'][k])
    im = a.imshow(G.T, origin='lower', extent=(sg[0], sg[-1], tt[0], tt[-1]), aspect='auto', cmap='gray', vmin=-mx, vmax=mx)
    a.grid(False); a.set_xlabel(r'$\sigma$'); a.set_title(ttl, fontsize=8)
    if k == 'dh':
        a.plot([R['rho'][0]], [R['rho'][1]], marker='x', color='white', ms=7, mew=1.5)
        a.plot([R['rho'][0]], [R['rho'][1]], marker='+', color='black', ms=7, mew=0.8)
ax[0].set_ylabel(r'$t$')
fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, label=r'$\log|f|\,\Delta\varphi$')
fig.savefig(os.path.join(args.figdir, 'bsy_riesz2d_maps.pdf'), bbox_inches='tight'); plt.close(fig)

# ------------------------------------------------------------------ 4. rates in two dimensions
Ns = np.array(R['Ns'], dtype=float); S = R['summary']
ests = [('plain', 'plain', 'o', 'white'), ('voronoi', 'Voronoi', 's', 'white'), ('delaunay', 'Delaunay', '^', 'white'), ('jittered', 'jittered', 'v', '0.6'), ('grid', 'grid', 'D', 'black')]
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.9), sharey=True)
for a, (kind, ttl) in zip(ax, (('zeta', r'$\zeta$'), ('dh', 'Davenport and Heilbronn'))):
    lst = list(ests) + ([('grid_defl_1.4e-4', 'deflated grid', '*', 'black')] if kind == 'dh' else [])
    for key, lab, m, fc in lst:
        r = [S[f'{kind}|{int(N)}|{key}']['rms'] for N in Ns]
        a.loglog(Ns, r, color='black', marker=m, mfc=fc, ms=4.5, lw=0.9, ls='--' if key == 'grid_defl_1.4e-4' else '-', label=lab)
    a.loglog(Ns, S[f'{kind}|{int(Ns[0])}|plain']['rms'] * (Ns / Ns[0]) ** -0.5, color='0.6', ls=':', lw=1.0, label='slope $-1/2$')
    a.loglog(Ns, S[f'{kind}|{int(Ns[0])}|delaunay']['rms'] * (Ns / Ns[0]) ** -1.5, color='0.3', ls=':', lw=1.0, label='slope $-3/2$')
    a.set_xlabel('$N$'); a.set_title(ttl, fontsize=9); a.legend(frameon=False, fontsize=6, loc='lower left')
ax[0].set_ylabel('root mean square error of $I(\\varphi)$')
fig.tight_layout(); fig.savefig(os.path.join(args.figdir, 'bsy_riesz2d_rates.pdf')); plt.close(fig)

# ------------------------------------------------------------------ 5. Levy density of the Davenport and Heilbronn function
Lv = load('bsy_levy.json')
u = np.array(Lv['u']); sel = u <= 2.0001
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
ax[0].semilogy(u[sel], np.array(Lv['W_explicit'])[sel], color='black', lw=1.4, label='explicit formula')
ax[0].semilogy(u[sel], np.array(Lv['W_zeros'])[sel], color='0.6', lw=1.0, ls='--', label='sum over zeros')
ax[0].semilogy(u[sel], np.abs(np.array(Lv['offline_abs'])[sel]), color='black', lw=1.0, ls=':', label='off line terms (absolute value)')
ax[0].set_ylim(1e-8, None); ax[0].set_xlabel('$u$'); ax[0].set_ylabel(r'$K_{\mathrm{DH}}(u)$'); ax[0].legend(frameon=False, fontsize=6.5, loc='upper right'); ax[0].set_title('(a)', fontsize=9)
tu = np.array(Lv['tail_u']); tr = np.array(Lv['tail_ratio']) - 1.0
pos = tr > 0
ax[1].semilogy(tu[pos], tr[pos], color='black', lw=1.2)
if (~pos).any():
    ax[1].semilogy(tu[~pos], -tr[~pos], color='0.6', lw=1.0, ls='--', label='negative values (absolute)')
    ax[1].legend(frameon=False, fontsize=6.5)
ax[1].set_xlabel('$u$'); ax[1].set_ylabel(r'$K_{\mathrm{DH}}(u)e^{\gamma_1u}-1$'); ax[1].set_title('(b)', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(args.figdir, 'bsy_dhlevy.pdf')); plt.close(fig)
# ------------------------------------------------------------------ 6. half disc rates
Hd = load('bsy_halfdisk.json')
Nh = np.array(Hd['Ns'], dtype=float)
fig, ax = plt.subplots(1, 3, figsize=(6.4, 2.5), sharey=True)
for a, (key, ttl) in zip(ax, (('zeta|0.75', r'$\zeta$, $\theta=3/4$'), ('zeta|0.5', r'$\zeta$, $\theta=1/2$'), ('dh|0.5', r'DH, $\theta=1/2$'))):
    C = Hd['cases'][key]
    for est, lab, m, fc in (('plain', 'plain', 'o', 'white'), ('trap', 'sorted trapezoid', 's', 'white'), ('trap_phi', r'with $\Phi$', '^', 'black')):
        a.loglog(Nh, C['rmse'][est], color='black', marker=m, mfc=fc, ms=3.8, lw=0.9, label=f"{lab} ({C['slope'][est]:.2f})")
    a.loglog(Nh, C['rmse']['plain'][0] * (Nh / Nh[0]) ** -0.5, color='0.6', ls=':', lw=1.0)
    a.loglog(Nh, C['rmse']['trap'][0] * (Nh / Nh[0]) ** -2.0, color='0.3', ls=':', lw=1.0)
    a.loglog(Nh, C['rmse']['trap_phi'][0] * (Nh / Nh[0]) ** -2.5, color='0.3', ls=':', lw=1.0)
    a.set_xlabel('$N$'); a.set_title(ttl, fontsize=8); a.legend(frameon=False, fontsize=5.2, loc='lower left')
ax[0].set_ylabel('root mean square error')
fig.tight_layout(); fig.savefig(os.path.join(args.figdir, 'bsy_halfdisk.pdf')); plt.close(fig)
print('wrote six figures to', args.figdir)
