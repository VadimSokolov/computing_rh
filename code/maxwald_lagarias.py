# maxwald_lagarias.py: the arithmetic approximants A_X of Proposition mw:prop:AX of Section sec:maxwald (ch/maxwald3.tex), subsection
# "The arithmetic form below alpha = 1, and how far it sees": the identity with the pole terms corrected (review-3 M19), Table
# mw:tab:AX, Figure mw:fig:lagarias (fig/maxwald_lagarias.pdf) and the failure of the converse (review-3 M20). Added for the book
# in October 2026 from the checks of the authors' chapter, misc/reviews/nick-maxwald/. For alpha = 1/2 + delta, t real, X > 0,
# g_X(x) = e^(-delta |x|) cos(t x) (1 - |x|/X)_+, h_X(r) = int g_X(x) e^(irx) dx = Re[J(-delta + i(t + r)) + J(-delta + i(t - r))] with
# J(c) = int_0^X e^(cx) (1 - x/X) dx = -1/c + (e^(cX) - 1)/(X c^2), the script computes, from Weil's explicit formula,
#   A_X(delta, t) = 1/2 [ h_X(i/2) + h_X(-i/2) - log pi + Gamma_X - 2 sum_{n < e^X} Lambda(n) n^(-alpha) cos(t log n) (1 - log n/X) ],
# h_X(i/2) + h_X(-i/2) = 2 Re[J(1/2 - delta + it) + J(-1/2 - delta + it)] (each of the two terms is half of this, which is the factor 2
# of M19), and Gamma_X = (1/2 pi) int h_X(r) Re psi(1/4 + ir/2) dr, evaluated as the single integral over [0, 2X] that Gauss's
# formula psi(z) = -gamma + int_0^oo (e^(-u) - e^(-zu))/(1 - e^(-u)) du gives with g_X(u/2) in place of the Fourier transform:
#   Gamma_X = -gamma_E + int_0^(2X) [e^(-u) - e^(-u/4) g_X(u/2)]/(1 - e^(-u)) du - log(1 - e^(-2X)).
#  identity:  A_X against (1/2) sum_rho h_X(gamma_rho) = sum_{k >= 1} h_X(gamma_k) over the first 1700 zeros plus the integral of the
#             non-oscillating part of h_X against the density log(x/2 pi)/(2 pi) from the point where the smooth count is 1700.5,
#             at delta in {0.25, 0.1, 0.05}, t in {0, 5, gamma_1, 20, 33}, X in {2, 4, 8, 13};
#  tab_AX:    min over t in [0, 60] of A_X(alpha - 1/2, t) for X in {2, 4, 8, 13}, alpha in {0.75, 0.6, 0.55}: on the grid of step
#             0.05 (what the chapter tabulates) and refined; and Re xi'/xi(alpha), the limit of the minima as X grows;
#  samples:   A_X at alpha = 0.75 and 0.6, t in {5, gamma_1, 20}, X in {4, 7, 10, 13}, against Re xi'/xi(alpha + it);
#  poles:     the two pole terms against quadrature of int g_X(x) e^(+-x/2) dx;
#  converse:  for the artificial zeros 1/2 +- b +- i gamma_0 (delta = 0.1, b = 0.2, gamma_0 = 20) the contribution of the four
#             members to A_X at t = gamma_0 for X = 5, 10, 20, 40, against the limit 1/(delta - b) + 1/(delta + b) the proof uses
#             (review-3 M20: it grows like e^((b - delta)X)/(X (b - delta)^2)), and its minimum over t near gamma_0;
#  figure:    A_X(1/4, t) for X = 4, 8, 13 against Re xi'/xi(0.75 + it) on [0, 60], with the ordinates gamma_k marked.
# Run from a folder holding zeros_1700.txt: python3 maxwald_lagarias.py [--small]; writes maxwald_lagarias.json and
# fig/maxwald_lagarias.pdf there (the figure needs matplotlib; python3 maxwald_lagarias.py --figure redraws it from
# maxwald_lagarias_curves.npz). --small leaves out X = 13 for a quick test.
import json
import math
import os
import sys
import time

import numpy as np
import mpmath as mp
from scipy import integrate
from scipy.optimize import minimize_scalar

GAMMA_E = 0.5772156649015329


def J(c, X):
    return -1 / c + np.expm1(c * X) / (X * c * c)


def lam_arrays(X):
    nmax = int(math.exp(X))
    sv = np.ones(nmax + 1, dtype=bool)
    sv[:2] = False
    for i in range(2, int(nmax ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    ps = np.nonzero(sv)[0].astype(np.int64)
    ns, ls = [], []
    cur, idx = ps.copy(), np.arange(len(ps))
    while cur.size:
        ns.append(cur)
        ls.append(np.log(ps[idx].astype(float)))
        nxt = cur * ps[idx]
        keep = nxt <= nmax
        cur, idx = nxt[keep], idx[keep]
    n = np.concatenate(ns)
    return np.log(n.astype(float)), np.concatenate(ls)


_LAM = {}


def lam(X):
    if X not in _LAM:
        _LAM[X] = lam_arrays(X)
    return _LAM[X]


def prime_sum(delta, ts, X):
    logn, Lam = lam(X)
    wgt = Lam * np.exp(-(0.5 + delta) * logn) * (1 - logn / X)
    out = np.empty(len(ts))
    B = max(1, 4_000_000 // len(logn))
    for i in range(0, len(ts), B):
        out[i:i + B] = np.cos(np.outer(ts[i:i + B], logn)) @ wgt
    return out


def gamma_term(delta, ts, X, panel=0.2, order=12):
    xg, wg = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(0, 2 * X, int(math.ceil(2 * X / panel)) + 1)
    a, b = edges[:-1, None], edges[1:, None]
    u = ((a + b) / 2 + (b - a) / 2 * xg[None, :]).ravel()
    wu = ((b - a) / 2 * wg[None, :]).ravel()
    x = u / 2
    base = np.exp(-u / 4 - delta * x) * (1 - x / X)
    den = -np.expm1(-u)
    out = np.empty(len(ts))
    B = max(1, 2_000_000 // len(u))
    for i in range(0, len(ts), B):
        cm = np.cos(np.outer(ts[i:i + B], x))
        out[i:i + B] = ((np.exp(-u)[None, :] - base[None, :] * cm) / den[None, :]) @ wu
    return -GAMMA_E + out - math.log1p(-math.exp(-2 * X))


def pole_terms(delta, ts, X):
    return 2 * np.real(J(0.5 - delta + 1j * ts, X) + J(-0.5 - delta + 1j * ts, X))


def AX(delta, ts, X, parts=False):
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    pole, gam, ps = pole_terms(delta, ts, X), gamma_term(delta, ts, X), prime_sum(delta, ts, X)
    A = 0.5 * (pole - math.log(math.pi) + gam - 2 * ps)
    return (A, pole, gam, ps) if parts else A


def h_real(r, delta, t, X):
    return np.real(J(-delta + 1j * (t + r), X) + J(-delta + 1j * (t - r), X))


def zero_side(delta, t, X, gam):
    K = len(gam)
    # smooth count N(T) = (T/2pi) log(T/(2 pi e)) + 7/8 = K + 1/2 fixes the start of the tail
    Hc = float(mp.findroot(lambda T: T / (2 * mp.pi) * mp.log(T / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8 - K - mp.mpf(1) / 2, gam[-1] + 1))

    def smooth(r):
        cp, cm = -delta + 1j * (t + r), -delta + 1j * (t - r)
        return np.real(-1 / cp - 1 / (X * cp * cp) - 1 / cm - 1 / (X * cm * cm))
    tail = integrate.quad(lambda r: smooth(r) * math.log(r / (2 * math.pi)) / (2 * math.pi), Hc, np.inf, limit=400)[0]
    return float(np.sum(h_real(gam, delta, t, X)) + tail), tail


def xi_logder_re(alpha, t):
    s = mp.mpc(alpha, t)
    return float(mp.re(1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)))


def quad_contrib(delta, b, g0, t, X):
    def h(r):
        return 0.5 * sum(J(-delta + 1j * (s1 * t + s2 * r), X) for s1 in (1, -1) for s2 in (1, -1))
    return 0.5 * (h(g0 - 1j * b) + h(g0 + 1j * b) + h(-g0 - 1j * b) + h(-g0 + 1j * b))


def make_figure(npz):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print('no matplotlib, figure not drawn:', e)
        return
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False,
                         'axes.spines.right': False})
    os.makedirs('fig', exist_ok=True)
    t, true, gam = npz['t'], npz['true'], npz['gamma']
    fig, ax = plt.subplots(figsize=(6.4, 3.1))
    ax.plot(t, true, color='k', lw=2.0, label=r"$\mathrm{Re}\,\xi'/\xi(0.75+it)$")
    for X, ls, c in ((4, ':', '0.6'), (8, '--', '0.4'), (13, '-', '0.25')):
        if f'A{X}' in npz.files:
            ax.plot(t, npz[f'A{X}'], ls=ls, color=c, lw=1.3, label=rf'$A_X$, $X={X}$')
    g = gam[gam <= t[-1]]
    ax.plot(g, np.full(len(g), -0.17), '|', color='k', ms=5, mew=0.8, label=r'$\gamma_k$')
    ax.axhline(0, color='0.5', lw=0.6)
    ax.set_xlabel(r'$t$')
    ax.set_ylim(-0.3, 4.4)
    ax.legend(fontsize=8, frameon=False, loc='upper left', ncol=1)
    fig.tight_layout()
    fig.savefig('fig/maxwald_lagarias.pdf')
    print('wrote fig/maxwald_lagarias.pdf')


def main():
    a = sys.argv[1:]
    if a and a[0] == '--figure':
        make_figure(np.load('maxwald_lagarias_curves.npz'))
        return
    small = '--small' in a
    t0 = time.time()
    zf = next(p for p in ('zeros_1700.txt', 'data/zeros_1700.txt', '../data/zeros_1700.txt') if os.path.exists(p))
    mp.mp.dps = 40
    zs = [mp.mpf(x) for x in open(zf).read().split()]
    gam = np.array([float(x) for x in zs])
    Xs = [2, 4, 8] if small else [2, 4, 8, 13]
    res = {'working_precision': 'double precision (numpy) for the sums and integrals; mpmath at 40 digits for Re xi\'/xi and the zeros',
           'primes_below_e_X': {str(X): int(len(lam(X)[0])) for X in Xs}}
    ident = []
    for X in Xs:
        for delta in (0.25, 0.1, 0.05):
            tt = [0.0, 5.0, float(gam[0]), 20.0, 33.0]
            A, pole, gm, ps = AX(delta, tt, X, parts=True)
            for i, t in enumerate(tt):
                z, tail = zero_side(delta, t, X, gam)
                ident.append({'delta': delta, 't': t, 'X': X, 'A_X_primes': float(A[i]), 'A_X_zeros': z, 'zero_tail': tail,
                              'abs_difference': float(abs(A[i] - z)), 'pole_terms': float(pole[i]), 'gamma_term': float(gm[i]),
                              'prime_sum': float(ps[i])})
    res['identity'] = {'worst_abs_difference': max(d['abs_difference'] for d in ident), 'rows': ident}
    ts = np.arange(0, 60 + 1e-9, 0.05)
    tab = {}
    for alpha in (0.75, 0.6, 0.55):
        d = alpha - 0.5
        row = {'xi_logder_at_alpha': xi_logder_re(alpha, 0)}
        for X in Xs:
            A = AX(d, ts, X)
            imin = int(np.argmin(A))
            best = (float(A[imin]), float(ts[imin]))
            for i in range(len(ts)):
                if (i == 0 or A[i] <= A[i - 1]) and (i == len(ts) - 1 or A[i] <= A[i + 1]):
                    r = minimize_scalar(lambda t: float(AX(d, [t], X)[0]), bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                                        method='bounded', options={'xatol': 1e-9})
                    if r.fun < best[0]:
                        best = (float(r.fun), float(r.x))
            row[f'X={X}'] = {'min_on_grid_0.05': float(A[imin]), 't_at_grid_min': float(ts[imin]), 'min_refined': best[0],
                             't_at_refined_min': best[1], 'A_X_at_t_0': float(A[0])}
        tab[f'alpha={alpha}'] = row
    res['tab_AX'] = tab
    smp = {}
    for alpha in (0.75, 0.6):
        for t in (5.0, float(gam[0]), 20.0):
            smp[f'alpha={alpha},t={t:.6f}'] = {'xi_logder_re': xi_logder_re(alpha, t),
                                               **{f'X={X}': float(AX(alpha - 0.5, [t], X)[0]) for X in ([4, 7, 10] + ([13] if not small else []))}}
    res['samples'] = smp
    # the two pole terms against quadrature
    d, t, X = 0.25, 5.0, 4.0
    g = lambda x: math.exp(-d * abs(x)) * math.cos(t * x) * (1 - abs(x) / X)
    hp = integrate.quad(lambda x: g(x) * math.exp(x / 2), -X, X, limit=400)[0]
    hm = integrate.quad(lambda x: g(x) * math.exp(-x / 2), -X, X, limit=400)[0]
    R_ = float(np.real(J(0.5 - d + 1j * t, X) + J(-0.5 - d + 1j * t, X)))
    res['poles'] = {'delta': d, 't': t, 'X': X, 'h_X_of_minus_i_over_2_quadrature': hp, 'h_X_of_i_over_2_quadrature': hm,
                    'Re_J_plus_J': R_, 'twice_Re_J_plus_J': 2 * R_, 'sum_of_the_two_h_values': hp + hm,
                    'chapter_reading_each_equals_2Re_gives_sum': 4 * R_}
    dl, b, g0 = 0.1, 0.2, 20.0
    conv = {'delta': dl, 'b': b, 'gamma_0': g0, 'limit_1_over_delta_minus_b_plus_1_over_delta_plus_b': 1 / (dl - b) + 1 / (dl + b), 'X': {}}
    for X in (5, 10, 20, 40):
        c0 = quad_contrib(dl, b, g0, g0, X)
        tg = np.arange(g0 - 3, g0 + 3, 0.001)
        cc = np.real(np.array([quad_contrib(dl, b, g0, t, X) for t in tg]))
        conv['X'][str(X)] = {'contribution_at_t_gamma_0': float(np.real(c0)), 'imag_part': float(np.imag(c0)),
                             'min_over_t_in_gamma_0_pm_3': float(cc.min()), 't_of_min': float(tg[np.argmin(cc)]),
                             'growth_estimate_e_to_b_minus_delta_X_over_X_b_minus_delta_sq': math.exp((b - dl) * X) / (X * (b - dl) ** 2)}
    res['converse'] = conv
    tf = np.arange(0, 60 + 1e-9, 0.05)
    curves = {'t': tf, 'gamma': gam[:30]}
    curves['true'] = np.array([xi_logder_re(0.75, t) for t in tf])
    for X in (4, 8) + (() if small else (13,)):
        curves[f'A{X}'] = AX(0.25, tf, X)
    np.savez('maxwald_lagarias_curves.npz', **curves)
    res['seconds'] = round(time.time() - t0, 1)
    json.dump(res, open('maxwald_lagarias.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('identity',)}, indent=1)[:7000])
    print('identity: worst abs difference', res['identity']['worst_abs_difference'])
    make_figure(np.load('maxwald_lagarias_curves.npz'))


if __name__ == '__main__':
    main()
