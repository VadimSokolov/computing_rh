# yo_ladder.py: Theorem yo:thm:rungs and its stress test (Remark "A stress test", line 105) and Figure yo:fig:W (line 362) of Chapter 17
# (Yoshida's Hermitian Forms, Weil Positivity and the Thorin Program). Added for the book in October 2026, written from the mathematics
# (the chapter's own ladder.py was never read or run).
# Model. T_k(x) = sum over the zeros with positive ordinate of w^k e^(-w x), w = gamma_rho^2, gamma_rho = gamma minus i delta; on the line w = gamma^2.
# With x = k/y and the positive factor e^(-k (log y minus 1)) removed, F_k(y) = sum over the zeros of e^(k g(w, y)), g(w, y) = log(w/y) minus w/y + 1.
# For a real zero g = h(gamma^2/y) with h(u) = log u minus u + 1 <= 0; the hypothetical pair (gamma, delta) = (60.5, 0.45) and its mirror 1 minus conj(rho)
# give w and conj(w), so their contribution is 2 e^(k Re g) cos(k Im g) with Re g = log(sec theta) + h(y_*/y), Im g = theta + 2 gamma delta/y, theta = minus 2
# arctan(delta/gamma), y_* = gamma^2 minus delta^2, sec theta minus 1 = 2 delta^2/(gamma^2 minus delta^2). T_k(k/y) < 0 iff F_k(y) < 0. The true zeros are those of
# data/zeros_1700.txt with 40 < gamma < 85 (the window of the remark), or all 1700 when a check says so.
# What is checked. (1) sec theta_* minus 1 = 1.1065e-4, y_* = 3660.05. (2) The sign of F_k on the 3001 points y = y_* (1 + s), s = minus 0.06 to 0.06 in steps of 4e-5, for
# every order k from 1 to 9348: the first failing order K_0 (chapter: 6349), the number of failing orders among 6349 to 9348 (3000), among 3349 to 6348 (0), the 395
# orders 6349 + 997 j below 4e5 (all fail), and K_0 (sec theta_* minus 1) (0.7025). (3) The same on a grid ten times finer with the minimum over y polished by a bounded
# scalar minimisation, which gives the first order at which T_k takes a negative value for some x (the grid is only a lower bound for that order); and the review's check
# of T_k >= 0 for all x > 0 at k = K_0 minus 1 (by the initial segment property T_(k minus 1)(x) = integral from x to infinity of T_k(s) ds one order suffices): the scale
# free margin m(y) = F/(A + |B|), A the sum over the real zeros and B the pair, on a fine grid of 200001 points of |log(y/y_*)| <= 1/2, on a wide grid y = 1 to 1e8 with all
# 1700 zeros, with polishing of every local minimum of m, and the analytic dominance of the real zeros for y outside (gamma_1^2, gamma_1700^2). (4) The zeros outside the
# window: the maximum over the grid of the ratio of their contribution to the main terms, the sign comparison of F with the window and with all 1700 zeros, and the bound
# of order e^(minus 0.23 k) from h at the edges of the window. (5) The heuristic K_0 of order gamma^2/(2 delta^2) (lines 102, 466) on further pairs (gamma, delta): K_0 by
# bisection on the grid and a scan of 120 orders around it, K_0 (sec theta minus 1), K_0 divided by gamma^2/(2 delta^2), and K_0 with all 1700 zeros. (6) Figure yo:fig:W:
# W(T, x) = sum over all zeros of e^(minus (gamma_rho minus T)^2 x) = sum over gamma > 0 of [e^(minus (T minus gamma)^2 x) + e^(minus (T + gamma)^2 x)] on [0, 60] at
# x = 0.02, 0.1, 0.5 from the 1700 zeros (terms beyond the file are below 1e-300 for T <= 60 and x >= 0.02), drawn black and white as fig/yo_W.pdf in the style of
# polya_fig.py; also the zero side of Table yo:tab:W to 15 digits (mpmath, 30 digits) for comparison with a prime side computed elsewhere.
# Error control. The sums are float64 with h(u) = log1p(u minus 1) minus (u minus 1) and log(sec theta) = log1p(2 delta^2/(gamma^2 minus delta^2)) evaluated without
# cancellation; the phase k Im g is accurate to about 1e-12 radians at k = 4e5. The sign of the minimum found at the first failing order and one order below is confirmed
# in mpmath (60 digits) from the definition T_k(x) = sum w^k e^(minus w x) with the 68 digit ordinates, and the two computations are compared. The grid results are
# numerical evidence (a grid cannot exclude a negative value between its points); the polishing and the margin m(y) say how far from zero the minimum is.
# Output: data/yo_ladder.json (numbers as strings, parameters, counts, the claims table with the chapter's printed digits against the computed ones, the sampled W curves),
# data/yo_ladder.log (the log written next to the JSON), fig/yo_W.pdf.
# Usage: python3 yo_ladder.py [--quick] [--out FILE] [--fig FILE] [--zeros FILE]; --quick runs the reference pair around K_0, one more pair and a coarse figure (about 5 s,
# writes to /tmp unless --out and --fig are given); the full case is one process and was timed to stay inside the local limit of 30 seconds.
import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp
from scipy.optimize import minimize_scalar

T0 = time.time()
LOG = []
CLAIMS = []
CHECKS = {}


def log(*a):
    line = '[%6.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(line, flush=True)
    LOG.append(line)


def check(name, cond):
    CHECKS[name] = bool(cond)
    log(('ok   ' if cond else 'FAIL ') + name)
    return cond


def claim(key, line, printed, computed, ok, note=''):
    CLAIMS.append({'claim': key, 'chapter_line': line, 'chapter': printed, 'computed': computed, 'agrees_to_printed_digits': bool(ok), 'note': note})
    log('claim %-60s line %-4s chapter %-12s computed %-12s %s %s' % (key, line, printed, computed, 'agrees' if ok else 'DIFFERS', note))


def find_zeros(arg):
    here = Path(__file__).resolve().parent
    for c in ([Path(arg)] if arg else []) + [here.parent / 'data' / 'zeros_1700.txt', here / 'zeros_1700.txt', Path('zeros_1700.txt')]:
        if c.is_file():
            return c
    sys.exit('zeros_1700.txt not found')


def hfun(eps):
    # h(u) = log u minus u + 1 for u = 1 + eps, without cancellation
    return np.log1p(eps) - eps


class Model:
    # F_k(y) on a fixed array of y for the real zeros gam (array of ordinates) and the pair (gp, dp)
    def __init__(self, gam, gp, dp, y):
        self.y = np.asarray(y, dtype=float)
        self.gam2 = np.asarray(gam, dtype=float) ** 2
        self.gp, self.dp = gp, dp
        self.G = hfun(self.gam2[:, None] / self.y[None, :] - 1.0)
        wr = gp * gp - dp * dp
        self.S = math.log1p(2 * dp * dp / wr)
        self.RG = self.S + hfun((wr - self.y) / self.y)
        self.IG = -2.0 * math.atan(dp / gp) + 2.0 * gp * dp / self.y

    def F(self, k):
        return np.exp(k * self.G).sum(axis=0) + 2.0 * np.exp(k * self.RG) * np.cos(k * self.IG)

    def parts(self, k, mask=None):
        kg = k * self.G
        M = np.maximum(kg.max(axis=0), k * self.RG)
        e = np.exp(kg - M)
        A = e.sum(axis=0) if mask is None else e[mask].sum(axis=0)
        B = 2.0 * np.exp(k * self.RG - M) * np.cos(k * self.IG)
        return A, B


def parts_at(gam2, gp, dp, k, y):
    # scaled parts (A, B) at the points y for the real zeros with squares gam2; the common scale exp(minus M) is removed
    y = np.atleast_1d(np.asarray(y, dtype=float))
    wr = gp * gp - dp * dp
    S = math.log1p(2 * dp * dp / wr)
    RG = S + hfun((wr - y) / y)
    IG = -2.0 * math.atan(dp / gp) + 2.0 * gp * dp / y
    kg = k * hfun(gam2[:, None] / y[None, :] - 1.0)
    M = np.maximum(kg.max(axis=0), k * RG)
    A = np.exp(kg - M).sum(axis=0)
    B = 2.0 * np.exp(k * RG - M) * np.cos(k * IG)
    return A, B


def margin_m(gam2, gp, dp, k, y):
    A, B = parts_at(gam2, gp, dp, k, y)
    return (A + B) / (A + np.abs(B))


def scan(model, ks):
    neg = np.zeros(len(ks), dtype=bool)
    fmin = np.zeros(len(ks))
    amin = np.zeros(len(ks), dtype=int)
    for i, k in enumerate(ks):
        F = model.F(float(k))
        j = int(F.argmin())
        fmin[i], amin[i], neg[i] = F[j], j, F[j] < 0
    return neg, fmin, amin


def find_K0(model, S):
    # first failing order on the grid by bisection (the failing orders form an up set), then a scan of 121 orders around it
    def f(k):
        # the sign of F = A + B is decided on the scaled parts, so that terms far below the largest one at a point cannot be lost to underflow
        A, B = model.parts(float(k))
        return bool((A + B).min() < 0)
    lo, hi = max(1, int(1e-4 / S)), int(10.0 / S)
    assert not f(lo) and f(hi), (lo, hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if f(mid):
            hi = mid
        else:
            lo = mid
    pat = {k: f(k) for k in range(max(1, hi - 60), hi + 61)}
    first = min(k for k, v in pat.items() if v)
    anomalies = [k for k, v in pat.items() if (v and k < hi) or ((not v) and k > first)]
    return first, anomalies


def mp_F(k, y, zero_strs, gp, dp, dps=60):
    mp.mp.dps = dps
    x = mp.mpf(k) / mp.mpf(float(y))
    tot = mp.mpf(0)
    for t in zero_strs:
        w = mp.mpf(t) ** 2
        tot += mp.exp(k * mp.log(w) - w * x)
    w = (mp.mpf(gp) - 1j * mp.mpf(dp)) ** 2
    tot += 2 * mp.re(mp.exp(k * mp.log(w) - w * x))
    scale = mp.exp(-k * (mp.log(mp.mpf(float(y))) - 1))
    return tot, tot * scale


def W_curve(gam, T, x):
    T = np.asarray(T, dtype=float)
    out = np.zeros_like(T)
    for a in range(0, len(T), 400):
        t = T[a:a + 400, None]
        out[a:a + 400] = (np.exp(-x * (t - gam[None, :]) ** 2) + np.exp(-x * (t + gam[None, :]) ** 2)).sum(axis=1)
    return out


def draw_W(path, T, curves, zero_marks):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})
    fig, ax = plt.subplots(figsize=(5.6, 2.5))
    for g in zero_marks:
        ax.axvline(g, color='0.8', lw=0.5, ls=':')
    styles = {0.02: dict(color='black', ls='-', lw=1.1), 0.1: dict(color='black', ls=(0, (5, 2)), lw=0.9), 0.5: dict(color='0.35', ls=(0, (1, 1.2)), lw=1.3)}
    for x, w in curves:
        ax.plot(T, w, label=r'$x=%g$' % x, **styles[x])
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 4.8)
    ax.set_xlabel(r'window centre $u$')
    ax.set_ylabel(r'$\mathcal{G}(u,x/2)+\mathcal{G}(-u,x/2)$')
    ax.legend(loc='upper left', frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--fig')
    ap.add_argument('--zeros')
    args = ap.parse_args()
    quick = args.quick
    root = Path(__file__).resolve().parent.parent
    out_path = args.out or ('/tmp/yo_ladder_quick.json' if quick else str(root / 'data' / 'yo_ladder.json'))
    fig_path = args.fig or ('/tmp/yo_W_quick.pdf' if quick else str(root / 'fig' / 'yo_W.pdf'))
    zpath = find_zeros(args.zeros)
    zero_strs = [line.split()[0] for line in open(zpath) if line.strip()]
    gam = np.array([float(t) for t in zero_strs])
    out = {'quick': quick, 'zeros_file': 'data/zeros_1700.txt', 'zeros_used': len(gam), 'highest_ordinate': zero_strs[-1][:12]}

    # the reference pair and its derived numbers
    GP, DP = 60.5, 0.45
    mp.mp.dps = 40
    gpm, dpm = mp.mpf('60.5'), mp.mpf('0.45')
    y_star = gpm ** 2 - dpm ** 2
    sec1 = 2 * dpm ** 2 / (gpm ** 2 - dpm ** 2)
    theta_star = -2 * mp.atan(dpm / gpm)
    out['pair'] = {'gamma': '60.5', 'delta': '0.45', 'y_star': mp.nstr(y_star, 15), 'sec_theta_star_minus_1': mp.nstr(sec1, 12), 'theta_star': mp.nstr(theta_star, 15),
                   'log_sec_theta_star': mp.nstr(mp.log1p(sec1), 12), 'gamma_squared_over_2_delta_squared': mp.nstr(gpm ** 2 / (2 * dpm ** 2), 10)}
    log('pair gamma = 60.5, delta = 0.45: y_* = %s, sec theta_* minus 1 = %s, theta_* = %s, gamma^2/(2 delta^2) = %s' % (mp.nstr(y_star, 12), mp.nstr(sec1, 8), mp.nstr(theta_star, 10), mp.nstr(gpm ** 2 / (2 * dpm ** 2), 8)))
    claim('sec(theta_*) minus 1', '105', '1.1065e-4', '%.4e' % float(sec1), '%.4e' % float(sec1) == '1.1065e-04'.replace('e-04', 'e-04'))
    claim('y_*', '105', '3660.05', '%.2f' % float(y_star), '%.2f' % float(y_star) == '3660.05')

    win = gam[(gam > 40) & (gam < 85)]
    win_strs = [t for t, g in zip(zero_strs, gam) if 40 < g < 85]
    out['window'] = {'lo': 40, 'hi': 85, 'zeros': len(win), 'first': '%.6f' % win[0], 'last': '%.6f' % win[-1], 'nearest_to_pair_below_and_above': ['%.6f' % win[win < GP][-1], '%.6f' % win[win > GP][0]]}
    log('window (40, 85): %d true zeros from %.4f to %.4f; the pair sits between the ordinates %.4f and %.4f' % (len(win), win[0], win[-1], win[win < GP][-1], win[win > GP][0]))
    S_ref = math.log1p(2 * DP * DP / (GP * GP - DP * DP))
    ys = (GP * GP - DP * DP) * (1 + np.linspace(-0.06, 0.06, 3001))
    model = Model(win, GP, DP, ys)

    # (2) the scan on the chapter's grid
    if quick:
        K0, anomalies = find_K0(model, S_ref)
        ks_scan = np.arange(max(1, K0 - 150), K0 + 151)
    else:
        ks_scan = np.arange(1, 9349)
    neg, fmin, amin = scan(model, ks_scan)
    K0 = int(ks_scan[neg].min())
    log('scan of %d orders (%d to %d) on the 3001 point grid, window zeros: first failing order K_0 = %d' % (len(ks_scan), ks_scan[0], ks_scan[-1], K0))
    out['scan'] = {'orders_scanned': [int(ks_scan[0]), int(ks_scan[-1])], 'K0_first_failing_on_grid': K0}
    i0 = int(np.where(ks_scan == K0)[0][0])
    out['scan']['min_F_at_K0'] = '%.6e' % fmin[i0]
    out['scan']['y_of_min_at_K0'] = '%.6f' % ys[amin[i0]]
    i1 = int(np.where(ks_scan == K0 - 1)[0][0])
    out['scan']['min_F_at_K0_minus_1'] = '%.6e' % fmin[i1]
    out['scan']['y_of_min_at_K0_minus_1'] = '%.6f' % ys[amin[i1]]
    log('min of F over the grid: K_0 minus 1 = %d: %.6e at y = %.3f; K_0 = %d: %.6e at y = %.3f' % (K0 - 1, fmin[i1], ys[amin[i1]], K0, fmin[i0], ys[amin[i0]]))
    if not quick:
        sel_hi = (ks_scan >= 6349) & (ks_scan <= 9348)
        sel_lo = (ks_scan >= 3349) & (ks_scan <= 6348)
        n_hi, n_lo = int(neg[sel_hi].sum()), int(neg[sel_lo].sum())
        below = int(neg[ks_scan < K0].sum())
        out['scan'].update({'failing_orders_among_6349_to_9348': n_hi, 'orders_in_that_range': int(sel_hi.sum()), 'failing_orders_among_3349_to_6348': n_lo, 'orders_in_that_range_lo': int(sel_lo.sum()),
                            'failing_orders_below_K0': below, 'failing_orders_among_K0_to_9348': int(neg[ks_scan >= K0].sum()), 'orders_K0_to_9348': int((ks_scan >= K0).sum())})
        out['scan']['non_failing_orders_in_6349_to_9348'] = [int(k) for k in ks_scan[sel_hi & ~neg]][:60]
        out['scan']['failing_orders_in_3349_to_6348'] = [int(k) for k in ks_scan[sel_lo & neg]][:60]
        log('failing orders among 6349 to 9348: %d of %d; among 3349 to 6348: %d of %d; among all orders below K_0: %d; among K_0 to 9348: %d of %d' % (n_hi, sel_hi.sum(), n_lo, sel_lo.sum(), below, neg[ks_scan >= K0].sum(), (ks_scan >= K0).sum()))
        claim('first failing order K_0', '105', '6349', str(K0), K0 == 6349)
        claim('failing orders among the 3000 orders 6349 to 9348', '105', '3000', str(n_hi), n_hi == 3000)
        claim('failing orders among the 3000 orders 3349 to 6348', '105', '0', str(n_lo), n_lo == 0)
        samp = [6349 + 997 * j for j in range(400) if 6349 + 997 * j < 4e5]
        negs, fminS, _ = scan(model, samp)
        out['sampled_orders'] = {'count': len(samp), 'failing': int(negs.sum()), 'first': samp[0], 'last': samp[-1], 'non_failing': [int(k) for k, v in zip(samp, negs) if not v]}
        log('sampled orders 6349 + 997 j < 4e5: %d orders, %d fail (last order %d)' % (len(samp), negs.sum(), samp[-1]))
        claim('sampled orders 6349 + 997 j below 4e5 that fail (395 orders)', '105', '395', str(int(negs.sum())), len(samp) == 395 and int(negs.sum()) == 395)
        if K0 != 6349:
            samp2 = [K0 + 997 * j for j in range(400) if K0 + 997 * j < 4e5]
            negs2, _, _ = scan(model, samp2)
            out['sampled_orders_from_my_K0'] = {'count': len(samp2), 'failing': int(negs2.sum())}
            log('sampled orders K_0 + 997 j < 4e5 with my K_0: %d orders, %d fail' % (len(samp2), negs2.sum()))
    prod = float(K0 * sec1)
    out['K0_times_sec_theta_star_minus_1'] = '%.6f' % prod
    claim('K_0 (sec theta_* minus 1)', '105', '0.7025', '%.4f' % prod, '%.4f' % prod == '0.7025')
    log('K_0 (sec theta_* minus 1) = %.6f; K_0 divided by gamma^2/(2 delta^2) = %.5f; K_0 times log sec theta_* = %.6f' % (prod, K0 / float(gpm ** 2 / (2 * dpm ** 2)), K0 * S_ref))

    # (3) a finer grid with polishing: the first order at which T_k takes a negative value for some x, and k = K_0 minus 1 for all x
    nfine = 6001 if quick else 30001
    ysf = (GP * GP - DP * DP) * (1 + np.linspace(-0.06, 0.06, nfine))
    modelf = Model(win, GP, DP, ysf)
    gam2w = win ** 2

    def polished(k):
        F = modelf.F(float(k))
        j = int(F.argmin())
        lo_, hi_ = ysf[max(j - 1, 0)], ysf[min(j + 1, nfine - 1)]
        f1 = lambda yy: float(modelf_scalar(k, yy))
        r = minimize_scalar(f1, bounds=(lo_, hi_), method='bounded', options={'xatol': 1e-9})
        return min(float(F[j]), float(r.fun)), float(r.x)

    def modelf_scalar(k, yy):
        yy = np.array([yy])
        wr = GP * GP - DP * DP
        RG = S_ref + hfun((wr - yy) / yy)
        IG = -2.0 * math.atan(DP / GP) + 2.0 * GP * DP / yy
        return (np.exp(k * hfun(gam2w[:, None] / yy[None, :] - 1.0)).sum(axis=0) + 2.0 * np.exp(k * RG) * np.cos(k * IG))[0]

    rng_k = range(K0 - 80, K0 + 41)
    pol = {k: polished(k) for k in rng_k}
    kneg = [k for k in rng_k if pol[k][0] < 0]
    Kc = min(kneg)
    anomalies = [k for k in rng_k if (pol[k][0] < 0 and k < Kc) or (pol[k][0] >= 0 and k > Kc)]
    out['fine_grid'] = {'points': nfine, 'K0_first_failing_continuous_y': Kc, 'anomalies_in_scan_of_121_orders': anomalies,
                        'min_F_at_Kc_minus_1': '%.6e' % pol[Kc - 1][0], 'y_at_Kc_minus_1': '%.6f' % pol[Kc - 1][1], 'min_F_at_Kc': '%.6e' % pol[Kc][0], 'y_at_Kc': '%.6f' % pol[Kc][1],
                        'min_F_at_6348': '%.6e' % pol[6348][0] if 6348 in pol else None, 'min_F_at_6349': '%.6e' % pol[6349][0] if 6349 in pol else None}
    log('grid of %d points with polishing: first failing order for some y = %d (the 3001 point grid gave %d); min F at that order minus 1: %.6e (y = %.3f), at that order: %.6e (y = %.3f); anomalies: %s' % (
        nfine, Kc, K0, pol[Kc - 1][0], pol[Kc - 1][1], pol[Kc][0], pol[Kc][1], anomalies or 'none'))
    check('the failing orders near K_0 form an up set on the fine polished grid (no anomaly in 121 orders)', not anomalies)

    # k = K_0 minus 1 and k = K_0: all x > 0 (the pair model with the real zeros of the file)
    allx = {}
    gam2all = gam ** 2
    ymid = GP * GP - DP * DP
    nfx = 40001 if quick else 200001
    ly = np.linspace(-0.5, 0.5, nfx)
    yfine = ymid * np.exp(ly)
    ywide = np.exp(np.arange(0.0, math.log(1e8), 1e-3 if not quick else 4e-3))
    for tag, k in (('K0_minus_1', Kc - 1), ('K0', Kc), ('6348_chapter_K0_minus_1', 6348)):
        if quick and tag.startswith('6348'):
            continue
        # fine grid with the window zeros (a lower bound for the sum of all, since the omitted terms are positive)
        mfine = np.concatenate([margin_m(gam2w, GP, DP, k, yfine[a:a + 2000]) for a in range(0, nfx, 2000)])
        mwide = np.concatenate([margin_m(gam2all, GP, DP, k, ywide[a:a + 500]) for a in range(0, len(ywide), 500)])
        # polish each local minimum of m below 0.3
        locs = [i for i in range(1, nfx - 1) if mfine[i] <= mfine[i - 1] and mfine[i] <= mfine[i + 1] and mfine[i] < 0.3]
        best = float(mfine.min())
        for i in locs:
            r = minimize_scalar(lambda t: float(margin_m(gam2w, GP, DP, k, ymid * math.exp(t))[0]), bounds=(ly[i - 1], ly[i + 1]), method='bounded', options={'xatol': 1e-12})
            best = min(best, float(r.fun))
        iw = int(mwide.argmin())
        allx[tag] = {'k': int(k), 'min_m_fine_polished': '%.6e' % best, 'fine_grid_points': nfx, 'local_minima_polished': len(locs), 'min_m_wide_grid': '%.6e' % float(mwide.min()),
                     'wide_grid_points': len(ywide), 'y_of_wide_min': '%.4f' % ywide[iw], 'positive_everywhere_on_both_grids': bool(best > 0 and mwide.min() > 0)}
        log('all x check at k = %d (%s): minimum of m = F/(A + |B|) on the fine grid after polishing %.6e (%d local minima polished), on the wide grid y = 1 to 1e8 (%d points, all %d zeros) %.6e at y = %.2f' % (
            k, tag, best, len(locs), len(ywide), len(gam), float(mwide.min()), ywide[iw]))
    # analytic dominance outside (gamma_1^2, gamma_max^2): the largest real term against the bound 2 e^(k (S + h)) of the pair
    k_ = Kc - 1
    dom_hi = math.log(gam2all[-1] / (GP * GP - DP * DP)) - (0.0)   # h(u_max) minus h(u_pair) = log(u_max/u_pair) minus (u_max minus u_pair) with u_max <= 1 at y >= gamma_max^2
    dom_lo = math.log(gam2all[0] / (GP * GP - DP * DP))
    allx['dominance'] = {'log_gamma_max2_over_y_star': '%.4f' % dom_hi, 'log_gamma_1_squared_over_y_star': '%.4f' % dom_lo,
                         'statement': 'for y >= gamma_max^2 the largest real term e^(k h(gamma_max^2/y)) exceeds 2 e^(k (S + h(y_*/y))) because h(u_max) minus h(u_pair) >= log(gamma_max^2/y_*) minus 1 = %.3f > S + (log 2)/k; for y <= gamma_1^2 the term of gamma_1 exceeds the pair because h(u_1) minus h(u_pair) >= (y_*/gamma_1^2 minus 1) minus log(y_*/gamma_1^2) = %.3f > S + (log 2)/k' % (dom_hi - 1.0, (GP * GP - DP * DP) / gam2all[0] - 1 - math.log((GP * GP - DP * DP) / gam2all[0]))}
    out['all_x'] = allx
    out['all_x']['initial_segment_property'] = 'T_(k minus 1)(x) = integral from x to infinity of T_k(s) ds, so if T_k >= 0 on (0, infinity) then T_j >= 0 for every j <= k'
    check('T_k > 0 for all x > 0 on the grids at k = %d (one order below the first failing order for some y), real zeros of the file plus the pair' % (Kc - 1), allx['K0_minus_1']['positive_everywhere_on_both_grids'])
    check('T_k < 0 somewhere at k = %d (the first failing order for some y)' % Kc, not allx['K0']['positive_everywhere_on_both_grids'])

    # (4) the zeros outside the window
    gam2a_all = gam ** 2
    maskw = (gam > 40) & (gam < 85)
    modelall = Model(gam, GP, DP, ys)
    outside = {}
    ks_chk = [3349, Kc - 1, K0, 6349, 9348] if not quick else [Kc - 1, K0]
    if not quick:
        ks_chk += [samp[-1]]
    for k in ks_chk:
        A_all, B = modelall.parts(float(k))
        A_win, _ = modelall.parts(float(k), mask=maskw)
        Fa, Fw = A_all + B, A_win + B
        sign_diff = int(((Fa < 0) != (Fw < 0)).sum())
        kG = float(k) * modelall.G
        Mx = np.maximum(kG.max(axis=0), float(k) * modelall.RG)
        nexcl = int((~maskw).sum())
        # the excluded terms underflow in float64, so bound them in logs: count times the largest excluded term, over the main terms A_win + |B|
        lr = float(((kG[~maskw].max(axis=0) - Mx) / math.log(10) + math.log10(nexcl) - np.log10(A_win + np.abs(B))).max())
        outside[str(k)] = {'log10_bound_of_excluded_terms_over_main_terms_max_over_grid': '%.1f' % lr, 'grid_points_where_the_sign_differs_between_window_and_all_zeros': sign_diff}
    hbound = lambda u: math.log(u) - u + 1
    y_lo, y_hi = ys[0], ys[-1]
    hb, ha = hbound(40.0 ** 2 / y_lo), hbound(85.0 ** 2 / y_hi)
    out['zeros_outside_the_window'] = {'per_order': outside, 'h_at_edges': {'below_40_largest': '%.5f' % hb, 'above_85_largest': '%.5f' % ha},
                                       'excluded_terms_at_most': 'about N e^(k max(h)) with max(h) = %.5f, so below 10^(minus %.0f) times the number of excluded zeros for every k >= 3349' % (max(hb, ha), -max(hb, ha) * 3349 / math.log(10))}
    log('zeros outside (40, 85): largest h at the window edges %.5f (below) and %.5f (above), so each excluded term is at most e^(%.5f k), below 1e-%.0f for k >= 3349' % (hb, ha, max(hb, ha), -max(hb, ha) * 3349 / math.log(10)))
    for k, d in outside.items():
        log('   k = %s: log10 of the bound (count times largest excluded term)/(main terms), largest over the grid: %s; grid points where the sign of F differs between the window and all 1700 zeros: %d' % (k, d['log10_bound_of_excluded_terms_over_main_terms_max_over_grid'], d['grid_points_where_the_sign_differs_between_window_and_all_zeros']))
    check('the sign of F on the grid is the same with the window zeros and with all 1700 zeros at every order tested', all(d['grid_points_where_the_sign_differs_between_window_and_all_zeros'] == 0 for d in outside.values()))

    # mpmath confirmation of the sign at the first failing order and one below
    mpc = {}
    for tag, k, ygrid in (('K0_minus_1', K0 - 1, ys[amin[i1]]), ('K0', K0, ys[amin[i0]])):
        tot, F = mp_F(k, ygrid, win_strs, GP, DP)
        Fd = float(model.F(float(k))[int(np.argmin(np.abs(ys - ygrid)))])
        mpc[tag] = {'k': k, 'y': '%.9f' % ygrid, 'T_k_sign': 'negative' if tot < 0 else 'positive', 'T_k(k/y)': mp.nstr(tot, 10), 'F_mpmath_60_digits': mp.nstr(F, 12), 'F_float64': '%.12e' % Fd,
                    'relative_difference': '%.2e' % float(abs(F - Fd) / abs(F))}
        log('mpmath at k = %d, y = %.6f: T_k(k/y) = %s (%s), F = %s against float64 %.10e' % (k, ygrid, mp.nstr(tot, 8), mpc[tag]['T_k_sign'], mp.nstr(F, 10), Fd))
    out['mpmath_confirmation'] = mpc
    check('mpmath (60 digits) agrees with float64 on the sign and to 1e-6 relative at the minima of F for K_0 minus 1 and K_0', all(float(d['relative_difference']) < 1e-6 for d in mpc.values()) and mpc['K0']['T_k_sign'] == 'negative' and mpc['K0_minus_1']['T_k_sign'] == 'positive')

    # (5) further pairs
    pairs = [(60.5, 0.45), (60.5, 0.30), (60.5, 0.15), (30.5, 0.45)] if quick else [
        (60.5, 0.45), (30.5, 0.45), (120.5, 0.45), (250.5, 0.45), (60.5, 0.30), (60.5, 0.20), (60.5, 0.15), (60.5, 0.10), (60.5, 0.05), (60.5, 0.02),
        (120.5, 0.25), (120.5, 0.10), (120.5, 0.05), (250.5, 0.15)]
    rows = []
    for (gp, dp) in pairs:
        wlo, whi = (40.0, 85.0) if (gp, dp) == (60.5, 0.45) else (0.66 * gp, 1.4 * gp)
        gw = gam[(gam > wlo) & (gam < whi)]
        ystar = gp * gp - dp * dp
        yg = ystar * (1 + np.linspace(-0.06, 0.06, 3001))
        S = math.log1p(2 * dp * dp / ystar)
        sec1p = 2 * dp * dp / ystar
        mdl = Model(gw, gp, dp, yg)
        K, anom = find_K0(mdl, S)
        mall = Model(gam, gp, dp, yg)
        def fails_all(kk):
            A_, B_ = mall.parts(float(kk))
            return bool((A_ + B_).min() < 0)
        Kall = K if (fails_all(K) and not fails_all(K - 1)) else None
        mwide_ = Model(gw, gp, dp, ystar * (1 + np.linspace(-0.3, 0.3, 15001)))

        def fails_wide(kk):
            A_, B_ = mwide_.parts(float(kk))
            return bool((A_ + B_).min() < 0)
        wide_ok = (not fails_wide(K - 1)) and fails_wide(K)
        Kw, anom_w = find_K0(mwide_, S)
        near = gw[np.argsort(np.abs(gw ** 2 - ystar))[:2]]
        rows.append({'gamma': gp, 'delta': dp, 'window': [round(wlo, 2), round(whi, 2)], 'window_zeros': int(len(gw)), 'y_star': '%.4f' % ystar, 'sec_theta_minus_1': '%.6e' % sec1p, 'K0': K,
                     'K0_times_sec_theta_minus_1': '%.4f' % (K * sec1p), 'K0_over_gamma2_over_2delta2': '%.4f' % (K / (gp * gp / (2 * dp * dp))), 'gamma2_over_2delta2': '%.1f' % (gp * gp / (2 * dp * dp)),
                     'K0_times_sec_theta_minus_1_over_delta': '%.4f' % (K * sec1p / dp), 'K0_times_delta_over_y_star': '%.4f' % (K * dp / ystar),
                     'K0_with_all_1700_zeros_agrees': Kall is not None, 'K0_unchanged_on_the_grid_y_star_(1_pm_0.3)_15001_points': bool(wide_ok), 'K0_on_that_wider_grid': Kw,
                     'K0_wide_times_sec_theta_minus_1': '%.4f' % (Kw * sec1p), 'anomalies_near_K0_wide_grid': anom_w, 'anomalies_near_K0': anom, 'nearest_ordinates_to_gamma': ['%.4f' % v for v in sorted(near)]})
        log('pair (%.1f, %.2f): window %d zeros, K_0 = %d, K_0 (sec theta minus 1) = %.4f, K_0/(gamma^2/(2 delta^2)) = %.4f, gamma^2/(2 delta^2) = %.0f, K_0 (sec theta minus 1)/delta = %.3f, same K_0 with all zeros: %s, same K_0 on the wider grid of plus or minus 30 percent: %s (K_0 there %d, times sec theta minus 1 %.4f), anomalies: %s' % (
            gp, dp, len(gw), K, K * sec1p, K / (gp * gp / (2 * dp * dp)), gp * gp / (2 * dp * dp), K * sec1p / dp, Kall is not None, wide_ok, Kw, Kw * sec1p, anom or 'none'))
    out['further_pairs'] = rows
    ratios = [float(r['K0_times_sec_theta_minus_1']) for r in rows]
    slopes = {}
    for g0 in sorted(set(r['gamma'] for r in rows)):
        sub = [r for r in rows if r['gamma'] == g0]
        if len(sub) >= 3:
            slopes['%g' % g0] = '%.3f' % np.polyfit(np.log([r['delta'] for r in sub]), np.log([float(r['K0_times_sec_theta_minus_1']) for r in sub]), 1)[0]
    ratios_w = [float(r['K0_wide_times_sec_theta_minus_1']) for r in rows]
    slopes_w = {}
    for g0 in sorted(set(r['gamma'] for r in rows)):
        sub = [r for r in rows if r['gamma'] == g0]
        if len(sub) >= 3:
            slopes_w['%g' % g0] = '%.3f' % np.polyfit(np.log([r['delta'] for r in sub]), np.log([float(r['K0_wide_times_sec_theta_minus_1']) for r in sub]), 1)[0]
    out['heuristic_summary_wide_grid'] = {'K0_times_sec_theta_minus_1_min': '%.4f' % min(ratios_w), 'max': '%.4f' % max(ratios_w), 'loglog_slope_against_delta_at_fixed_gamma': slopes_w}
    log('on the wider grid (plus or minus 30 percent) K_0 (sec theta minus 1) lies between %.4f and %.4f, log log slopes against delta %s' % (min(ratios_w), max(ratios_w), slopes_w))
    out['heuristic_summary'] = {'K0_times_sec_theta_minus_1_min': '%.4f' % min(ratios), 'max': '%.4f' % max(ratios), 'loglog_slope_of_K0_sec_minus_1_against_delta_at_fixed_gamma': slopes,
                                'statement': 'K_0 (sec theta minus 1) lies between %.2f and %.2f over %d pairs; against delta at fixed gamma its log log slope is %s (slope 0 would be K_0 of order gamma^2/delta^2, slope 1 is K_0 of order gamma^2/delta)' % (min(ratios), max(ratios), len(ratios), slopes)}
    log(out['heuristic_summary']['statement'])

    # (6) Figure yo:fig:W and Table yo:tab:W (zero side)
    Tg = np.linspace(0.0, 60.0, 241 if quick else 1201)
    curves = [(x, W_curve(gam, Tg, x)) for x in (0.02, 0.1, 0.5)]
    zero_marks = [g for g in gam if g <= 60.0]
    os.makedirs(os.path.dirname(os.path.abspath(fig_path)), exist_ok=True)
    draw_W(fig_path, Tg, curves, zero_marks)
    Ts = np.arange(0.0, 60.0001, 0.1)
    Wd = {'T': ['%.1f' % t for t in Ts]}
    for x, _ in curves:
        Wd['W_x=%g' % x] = ['%.10g' % v for v in W_curve(gam, Ts, x)]
    Wd['zeros_marked'] = len(zero_marks)
    peaks = {'W(gamma_n, 0.5) for the %d zeros up to 60' % len(zero_marks): ['%.6f' % v for v in W_curve(gam, np.array(zero_marks), 0.5)]}
    Wd['peaks_at_x_0.5'] = peaks
    Wd['W(0, x) = 2 Theta(x)'] = {'%g' % x: '%.10g' % float(W_curve(gam, np.array([0.0]), x)[0]) for x in (0.02, 0.1, 0.5)}
    Wd['range_over_[0,60]'] = {'x=%g' % x: ['%.6f' % w.min(), '%.6f' % w.max()] for x, w in curves}
    out['W'] = Wd
    log('figure written to', fig_path, '; W(0, x) = 2 Theta(x) at x = 0.02, 0.1, 0.5:', Wd['W(0, x) = 2 Theta(x)'], '; ranges over [0, 60]:', Wd['range_over_[0,60]'])
    log('W(gamma_n, 0.5) for the zeros up to 60:', ' '.join(list(peaks.values())[0]))
    mp.mp.dps = 30
    tab = []
    for (x, T) in ((0.01, 0.0), (0.01, 14.134725), (0.05, 17.0), (0.2, 14.134725), (0.2, 17.0), (0.5, 21.02204), (0.5, 30.0)):
        tot = mp.mpf(0)
        xm, Tm = mp.mpf(str(x)), mp.mpf(str(T))
        for t in zero_strs:
            g = mp.mpf(t)
            tot += mp.exp(-xm * (g - Tm) ** 2) + mp.exp(-xm * (g + Tm) ** 2)
        tab.append({'x': str(x), 'T': str(T), 'zero_side': mp.nstr(tot, 15)})
    out['W_table_zero_side'] = tab
    log('Table yo:tab:W zero side to 15 digits:', '; '.join('(%s, %s): %s' % (r['x'], r['T'], r['zero_side']) for r in tab))

    out['claims'] = CLAIMS
    out['checks'] = CHECKS
    out['all_checks_pass'] = all(CHECKS.values())
    out['elapsed_s'] = round(time.time() - T0, 1)
    log('all checks pass:', out['all_checks_pass'], '; claims that differ from the printed digits:', [c['claim'] for c in CLAIMS if not c['agrees_to_printed_digits']] or 'none', '; elapsed %.1f s' % (time.time() - T0))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(out, fh, indent=1)
    with open(os.path.splitext(out_path)[0] + '.log', 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    log('wrote', out_path, 'and', os.path.splitext(out_path)[0] + '.log')


if __name__ == '__main__':
    main()
