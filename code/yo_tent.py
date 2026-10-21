#!/usr/bin/env python3
# yo_tent.py: the tent coefficient D(x) = T(Delta_x) of Chapter 17 (Yoshida's Hermitian Forms, Weil Positivity and the Thorin Program), group Y1 of the number inventory
# (misc/reviews/nick-yoshida/review-A.md, Section 5), together with yo_explicit.py. Added for the book in October 2026, written from the mathematics (the collaborator's scripts were not read).
# Checks (chapter lines in the source ch_yoshida.tex). Line 149: the minimum of D on [0.05, 16] (chapter: 0.0550 at x = 8.04, grid step 1e-4, exact Riesz mean between prime powers) and the mean of D on [2, 16] (chapter: 0.092409 against 4 kappa_1 = 0.092420).
# Figure yo:fig:D (line 155): the zero side curve and the prime side dots on [0.02, 8] with verticals at log n for n = 2, 3, 4, 5, 7, 8, 9, written to fig/yo_D.pdf. Line 152 and Figure yo:fig:off: the first x on [0.5, 60] (step 0.01) where D plus the pair term 4 Re[4 sin^2((gamma minus i delta) x/2)/(gamma minus i delta)^2]
# is negative for (gamma, delta) = (100, 0.3), (100, 0.2), (1000, 0.3) (chapter: 17.3, 25.4, 32.9 against the predictions (2 log gamma + log kappa_1)/delta = 18.1, 27.2, 33.5), and the minima over windows of length 0.4, written to fig/yo_off.pdf.
# Lines 375 and 382 (Proposition yo:prop:barrier): epsilon_H = 8 Sigma_0 cosh^2(x/4) with Sigma_0 < 1.48e-12 (Theorem ag:thm:explicit): 1.18e-11 cosh^2(x/4), the level where epsilon_H = 0.05 (chapter: 47.10) and e^47.1 (chapter: about 2.9e20), and the claim that D_on stays above 0.05 on [16, 47.1],
# which the chapter never computed (finding M6 of review-A.md): here a rigorous lower bound of D_on on [16, 47.1] from the certified zeros, and the first level at which the argument |D minus D_on| <= epsilon_H stops proving D > 0.
# Method. Prime side: D(x) = 32 sinh^2(x/4) minus x log(pi) minus 2 (x S0(x) minus S1(x)) + V_1(Delta_x), with S0 and S1 the sums over prime powers n <= e^x of Lambda(n)/sqrt(n) and Lambda(n) log(n)/sqrt(n) (the Riesz mean of order one, exact between prime powers, built from a sieve of all prime powers up to e^16
# by compensated cumulative sums) and V_1(Delta_x) = x psi(1/4) + psi'(1/4)/2 minus 2 sum over k >= 0 of exp(minus x(2k + 1/2))/(2k + 1/2)^2 (the residue series, which yo_explicit.py shows equal to the closed form of Proposition yo:prop:tent (c) to 29 digits). The mean of D over an interval is integrated exactly
# (antiderivatives of the archimedean part, of V_1 and of the Riesz mean prime by prime), not by a trapezoid rule, because the kinks of D at log n bias a trapezoid rule at the 1e-6 level. Zero side: D(x) = sum over gamma > 0 of 8 sin^2(gamma x/2)/gamma^2 over the first N certified zeros of data/zeros_1700.txt (N = 1700, and N = 1500 for the chapter's count),
# plus the smooth tail 2 (log(T*/2 pi) + 1)/(pi T*) with the tail starting at the point T* where theta(T*)/pi + 1 = N + 1/2. The pair term is evaluated in complex arithmetic. The rigorous floor: every term of D_on is nonnegative, so L_N(x) = sum over the first N certified ordinates of 8 sin^2(gamma x/2)/gamma^2 is a lower bound of D_on(x);
# its minimum over an interval is bounded below by branch and bound with the Lipschitz constant 4 sum 1/gamma (the derivative of L_N is sum 4 sin(gamma x)/gamma), to the tolerance printed, with an allowance of 1e-9 for double precision rounding (checked against mpmath at 30 digits at the minimiser).
# Error control. The prime side and the zero side are compared on a grid of 20001 points (maximum and median difference printed); the exact mean is compared with the trapezoid mean and with the zero side formula for the mean, 4 kappa_1 minus (4/(b minus a)) sum [sin(b gamma) minus sin(a gamma)]/gamma^3, which needs no tail estimate beyond 1e-7;
# the first negative values are found on the grid of step 0.01 as in the chapter and then refined by bisection; the figures are drawn from the arrays computed here.
# Output: data/yo_tent.json (values as strings), data/yo_tent.log, fig/yo_D.pdf, fig/yo_off.pdf (black and white: line styles and markers, gray guides).
# Usage: python3 yo_tent.py [--quick] [--out FILE] [--zeros FILE] [--figdir DIR] [--no-fig]; --quick uses primes to e^10 and coarser grids (about 4 seconds in one process) and writes to the temporary folder unless --out and --figdir are given.
# The full case (primes to e^16) takes about 7 seconds in one process (the data file and the figures of the book were made by a local run within the 30 second limit). On the cluster: bash misc/tools/hopper_run.sh -c 1 -t 10 -g yo_tent.json -g yo_tent.log -g yo_D.pdf -g yo_off.pdf code/yo_tent.py data/zeros_1700.txt -- --out yo_tent.json --figdir .
import argparse
import json
import math
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.optimize import minimize_scalar

plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LOG = []
CLAIMS = []
LOGPI = math.log(math.pi)
KAPPA1_STR = '0.0231049931154189707889'   # the book's sum of a_k^(-1) (Chapter 18), equal to (1/2)(log xi)''(1/2); yo_explicit.py computes it from log xi to 25 digits
KAPPA1 = float(KAPPA1_STR)
TRAPZ = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
PSI14 = float(mp.digamma(mp.mpf(1) / 4))
PSI1_14 = float(mp.psi(1, mp.mpf(1) / 4))


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def sg(v, d=10):
    return '%.*g' % (d, v)


def claim(label, line, chapter, mine, status):
    CLAIMS.append(dict(claim=label, line=line, chapter=chapter, mine=mine, status=status))
    log('CLAIM %-52s line %-4s chapter %-22s mine %-26s %s' % (label, line, chapter, mine, status))


def find_zeros(path):
    cands = ([Path(path)] if path else []) + [Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), ROOT / 'data' / 'zeros_1700.txt', HERE / 'zeros_1700.txt']
    for c in cands:
        if c.exists():
            return c
    sys.exit('zeros_1700.txt not found: pass --zeros')


def tail_star(N, gN):
    mp.mp.dps = 20
    return float(mp.findroot(lambda T: mp.siegeltheta(T) / mp.pi + 1 - (N + mp.mpf(1) / 2), mp.mpf(gN) + mp.mpf(1) / 2))


def mean_tail4(T):
    return 2 * (math.log(T / (2 * math.pi)) + 1) / (math.pi * T)


def rvm_tail(x, T):
    c = 2 * math.pi
    f = math.log(T / c) / T ** 2
    fp = (1 - 2 * math.log(T / c)) / T ** 3
    return mean_tail4(T) + (2 / math.pi) * (np.sin(x * T) * f / x + np.cos(x * T) * fp / x ** 2)


def Lsum(x, g, w, chunk=1200):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    out = np.empty(x.shape)
    for i in range(0, len(x), chunk):
        s = np.sin(0.5 * np.outer(x[i:i + chunk], g))
        out[i:i + chunk] = (s * s) @ w
    return out


def cumsum_comp(w):
    out = np.empty(len(w))
    s = 0.0
    c = 0.0
    for i, v in enumerate(w.tolist()):
        t = s + v
        if abs(s) >= abs(v):
            c += (s - t) + v
        else:
            c += (v - t) + s
        s = t
        out[i] = s + c
    return out


def V1_series(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    for i in range(0, len(x), 4000):
        xb = x[i:i + 4000]
        K = int(20 / xb.min()) + 3
        c = 2 * np.arange(K) + 0.5
        s = (np.exp(-np.outer(xb, c)) / c ** 2).sum(axis=1)
        out[i:i + 4000] = xb * PSI14 + PSI1_14 / 2 - 2 * s
    return out


def V1_antiderivative(x):
    K = int(20 / x) + 3
    c = 2 * np.arange(K) + 0.5
    return x * x / 2 * PSI14 + x * PSI1_14 / 2 + 2 * float((np.exp(-x * c) / c ** 3).sum())


class Primes:
    def __init__(self, emax):
        N = int(math.floor(math.exp(emax)))
        sieve = np.ones(N + 1, dtype=bool)
        sieve[:2] = False
        sieve[4::2] = False
        for p in range(3, int(N ** 0.5) + 1, 2):
            if sieve[p]:
                sieve[p * p::2 * p] = False
        primes = np.flatnonzero(sieve).astype(np.int64)
        ns, ps = [primes], [primes]
        for p in primes[primes ** 2 <= N].tolist():
            q = p * p
            while q <= N:
                ns.append(np.array([q], dtype=np.int64))
                ps.append(np.array([p], dtype=np.int64))
                q *= p
        n = np.concatenate(ns)
        p = np.concatenate(ps)
        o = np.argsort(n, kind='stable')
        n, p = n[o], p[o]
        self.N, self.emax, self.count = N, emax, len(n)
        self.logn = np.log(n.astype(float))
        lam = np.log(p.astype(float))
        self.w0 = lam / np.sqrt(n.astype(float))
        self.w1 = self.w0 * self.logn
        self.cum0 = cumsum_comp(self.w0)
        self.cum1 = cumsum_comp(self.w1)

    def riesz(self, x):
        idx = np.searchsorted(self.logn, x, side='right')
        S0 = np.where(idx > 0, self.cum0[idx - 1], 0.0)
        S1 = np.where(idx > 0, self.cum1[idx - 1], 0.0)
        return x * S0 - S1

    def D(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        return 16 * (np.cosh(x / 2) - 1) - x * LOGPI - 2 * self.riesz(x) + V1_series(x)

    def mean(self, a, b):
        arch = lambda x: 32 * math.sinh(x / 2) - 16 * x - x * x / 2 * LOGPI
        l = self.logn[self.logn < b]
        w = self.w0[:len(l)]
        terms = w * (np.maximum(b - l, 0.0) ** 2 - np.maximum(a - l, 0.0) ** 2)
        rz = -math.fsum(terms.tolist())
        return (arch(b) - arch(a) + V1_antiderivative(b) - V1_antiderivative(a) + rz) / (b - a)


def certified_min(f, lip, a, b, h0, tol):
    """minimum of f on [a, b] with a certified lower bound: returns (upper bound of min, lower bound of min, argmin cell midpoint, evaluations)."""
    n0 = int(math.ceil((b - a) / h0))
    edges = np.linspace(a, b, n0 + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    w = np.full(n0, (b - a) / n0)
    fm = f(mids)
    ev = n0
    lb = fm - lip * w / 2
    best = float(fm.min())
    arg = float(mids[int(fm.argmin())])
    while True:
        keep = lb < best - tol
        if not keep.any():
            break
        m = mids[keep]
        ww = w[keep] / 2
        mids = np.concatenate([m - ww / 2, m + ww / 2])
        w = np.concatenate([ww, ww])
        fm = f(mids)
        ev += len(mids)
        lb = fm - lip * w / 2
        if fm.min() < best:
            best = float(fm.min())
            arg = float(mids[int(fm.argmin())])
    return best, best - tol, arg, ev


def bisect_root(f, a, b, it=60):
    fa = f(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        fmv = f(m)
        if (fmv < 0) == (fa < 0):
            a, fa = m, fmv
        else:
            b = m
    return 0.5 * (a + b)


def figure_D(path, prm, g, w, T, N):
    xs = np.linspace(0.02, 8.0, 7981)
    zero = Lsum(xs, g, w) + mean_tail4(T)
    xd = np.arange(0.02, 8.0 + 1e-9, 0.05)
    dots = prm.D(xd)
    fig, ax = plt.subplots(figsize=(5.6, 2.6))
    ax.axvspan(0, math.log(2), color='0.88', lw=0)
    for n in (2, 3, 4, 5, 7, 8, 9):
        ax.axvline(math.log(n), color='0.55', lw=0.5, ls=':')
    ax.plot(xs, zero, color='black', lw=0.7)
    ax.plot(xd, dots, ls='none', marker='o', mfc='white', mec='black', ms=2.3, mew=0.5)
    ax.axhline(4 * KAPPA1, color='black', lw=0.7, ls='--')
    ax.set_xlim(0, 8.1)
    ax.set_ylim(0, 0.2)
    ax.set_xlabel(r'level $x=2a$')
    ax.set_ylabel(r'$D(x)$')
    sec = ax.secondary_xaxis('top')
    sec.set_xticks([math.log(n) for n in (2, 3, 4, 5, 7, 8, 9)])
    sec.set_xticklabels(['2', '3', '4', '5', '7', '8', '9'], fontsize=6)
    sec.set_xlabel(r'$n$ at $x=\log n$', fontsize=6, labelpad=2)
    sec.spines['top'].set_visible(False)
    handles = [Patch(facecolor='0.88', label=r'$x\leq\log 2$ (Yoshida, Theorem 1)'),
               plt.Line2D([], [], color='black', lw=0.7, label='zero side, %d zeros plus smooth tail' % N),
               plt.Line2D([], [], ls='none', marker='o', mfc='white', mec='black', ms=2.3, mew=0.5, label='prime side, Riesz mean'),
               plt.Line2D([], [], color='black', lw=0.7, ls='--', label=r'$4p_1=4\sum\gamma^{-2}$')]
    ax.legend(handles=handles, loc='upper center', ncol=2, frameon=True, framealpha=1.0, facecolor='white', edgecolor='none', fontsize=6.5, bbox_to_anchor=(0.5, 1.0)).set_zorder(5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log('wrote', path)


def figure_off(path, centres, base, pairs):
    fig, ax = plt.subplots(figsize=(5.6, 2.7))
    ax.axhline(0.0, color='0.3', lw=0.5)
    styles = [dict(ls='--', color='black', lw=0.9), dict(ls='-.', color='0.35', lw=0.9), dict(ls=':', color='black', lw=1.4)]
    for (lab, xstar, arr), st in zip(pairs, styles):
        ax.axvline(xstar, color='0.65', lw=0.6, ls=(0, (1, 2)))
        ax.text(xstar, 0.125, r'$x_*=%.1f$' % xstar, rotation=90, ha='right', va='bottom', fontsize=6, color='0.25')
    ax.plot(centres, base, color='black', lw=1.0, marker='o', ms=1.6, label='actual zeros')
    for (lab, xstar, arr), st in zip(pairs, styles):
        ax.plot(centres, arr, label=lab, **st)
    ax.set_xlim(0, 62)
    ax.set_ylim(-0.25, 0.12)
    ax.set_xlabel(r'level $x$')
    ax.set_ylabel(r'local minimum of $D(x)$')
    ax.legend(loc='lower right', frameon=False, fontsize=6.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log('wrote', path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--figdir', default=None)
    ap.add_argument('--no-fig', action='store_true')
    args = ap.parse_args()
    quick = args.quick
    t0 = time.time()
    emax = 10.0 if quick else 16.0
    step = 1e-3 if quick else 1e-4
    outp = Path(args.out) if args.out else (Path(tempfile.gettempdir()) / 'yo_tent_quick.json' if quick else ROOT / 'data' / 'yo_tent.json')
    figdir = Path(args.figdir) if args.figdir else (Path(tempfile.gettempdir()) if quick else ROOT / 'fig')
    zp = find_zeros(args.zeros)
    with open(zp) as fh:
        strs = [s.strip() for s in fh if s.strip()]
    g_all = np.array([float(s) for s in strs])
    N = len(g_all)
    out = dict(script='code/yo_tent.py', quick=bool(quick), zeros_file=zp.name, zeros_loaded=N, zeros_first=strs[0][:25], zeros_last=strs[-1][:25], emax=emax, grid_step=step)
    log('yo_tent.py: %d certified zeros from %s (gamma_1 = %s, gamma_%d = %s), primes to e^%g, grid step %g, quick = %s' % (N, zp, strs[0][:15], N, strs[-1][:15], emax, step, quick))
    cuts = {}
    for n in sorted({1500, N}):
        g = g_all[:n]
        cuts[n] = dict(g=g, w=8 / g ** 2, T=tail_star(n, g[-1]))
        log('zero count %d: gamma_N = %.6f, tail start T* = %.6f, mean tail of D = %.8f' % (n, g[-1], cuts[n]['T'], mean_tail4(cuts[n]['T'])))
    g, w, T = cuts[N]['g'], cuts[N]['w'], cuts[N]['T']
    lip = 4 * float(np.sum(1 / g))
    out['zero_counts'] = {str(n): dict(gamma_N=sg(cuts[n]['g'][-1], 12), T_star=sg(cuts[n]['T'], 12), mean_tail=sg(mean_tail4(cuts[n]['T']), 10)) for n in cuts}
    out['lipschitz_constant_of_L_N'] = sg(lip, 8)

    # ---- primes
    t1 = time.time()
    prm = Primes(emax)
    log('prime powers to e^%g = %d: %d of them (%.1f s)' % (emax, prm.N, prm.count, time.time() - t1))
    out['primes'] = dict(nmax=prm.N, prime_powers=prm.count, sum_Lambda_over_sqrt=sg(prm.cum0[-1], 12))

    # ---- minimum of D on [0.05, emax]
    xg = 0.05 + step * np.arange(int(round((emax - 0.05) / step)) + 1)
    Dg = prm.D(xg)
    i0 = int(Dg.argmin())
    loc = np.flatnonzero((Dg[1:-1] < Dg[:-2]) & (Dg[1:-1] <= Dg[2:])) + 1
    low5 = loc[np.argsort(Dg[loc])[:6]]
    r = minimize_scalar(lambda t: float(prm.D(np.array([t]))[0]), bounds=(xg[i0] - 3 * step, xg[i0] + 3 * step), method='bounded', options={'xatol': 1e-12})
    xmin, dmin = float(r.x), float(r.fun)
    zmin = {n: float(Lsum(np.array([xmin]), cuts[n]['g'], cuts[n]['w'])[0] + mean_tail4(cuts[n]['T'])) for n in cuts}
    log('\nminimum of D on [0.05, %g] (prime side, grid step %g): %.8f at x = %.4f; refined %.10f at x = %.6f' % (emax, step, Dg[i0], xg[i0], dmin, xmin))
    log('the six lowest local minima on the grid (x, D): ' + ', '.join('(%.4f, %.6f)' % (xg[j], Dg[j]) for j in low5))
    log('zero side with smooth tail at the minimiser: ' + ', '.join('N = %d: %.8f' % (n, v) for n, v in zmin.items()))
    out['minimum'] = dict(interval=[0.05, emax], grid_min=sg(Dg[i0], 10), grid_argmin=sg(xg[i0], 8), refined_min=sg(dmin, 12), refined_argmin=sg(xmin, 10),
                          lowest_local_minima=[[sg(xg[j], 8), sg(Dg[j], 10)] for j in low5], zero_side_at_minimiser={str(n): sg(v, 10) for n, v in zmin.items()})
    claim('min of D on [0.05, 16] (prime side)', '149', '0.0550 at x = 8.04', '%.4f at x = %.2f' % (dmin, xmin),
          'agrees' if (abs(dmin - 0.0550) < 5e-5 and abs(xmin - 8.04) < 5e-3) else ('not applicable in --quick (primes to e^%g)' % emax if quick else 'differs'))

    # ---- mean on [2, emax]
    a_m, b_m = 2.0, emax
    mean_exact = prm.mean(a_m, b_m)
    sel = (xg >= a_m - 1e-12) & (xg <= b_m + 1e-12)
    mean_trap = float(TRAPZ(Dg[sel], xg[sel]) / (xg[sel][-1] - xg[sel][0]))
    zmean = {}
    for n in cuts:
        gg = cuts[n]['g']
        zmean[n] = 4 * KAPPA1 - 4 / (b_m - a_m) * float(np.sum((np.sin(b_m * gg) - np.sin(a_m * gg)) / gg ** 3))
    log('\nmean of D on [%g, %g]: exact integration of the prime side %.9f, trapezoid on the grid %.9f (difference %.2e), 4 kappa_1 = %.9f' % (a_m, b_m, mean_exact, mean_trap, mean_trap - mean_exact, 4 * KAPPA1))
    log('mean from the zeros, 4 kappa_1 minus (4/(b minus a)) sum [sin(b gamma) minus sin(a gamma)]/gamma^3: ' + ', '.join('N = %d: %.9f' % (n, v) for n, v in zmean.items()))
    log('mean minus 4 kappa_1 = %.3e; the zero side formula gives %.3e' % (mean_exact - 4 * KAPPA1, zmean[N] - 4 * KAPPA1))
    out['mean'] = dict(interval=[a_m, b_m], exact_prime_side=sg(mean_exact, 12), trapezoid=sg(mean_trap, 12), four_kappa1=sg(4 * KAPPA1, 12), zero_side_formula={str(n): sg(v, 12) for n, v in zmean.items()},
                       mean_minus_four_kappa1=sg(mean_exact - 4 * KAPPA1, 6))
    claim('mean of D on [2, 16] (prime side)', '149', '0.092409', '%.6f' % mean_exact, 'agrees' if (not quick and abs(mean_exact - 0.092409) < 5e-7) else ('not applicable in --quick' if quick else 'differs'))
    claim('4 kappa_1', '149', '0.092420', '%.6f' % (4 * KAPPA1), 'agrees' if abs(4 * KAPPA1 - 0.092420) < 5e-7 else 'differs')

    # ---- prime side against zero side on a grid
    xs = np.linspace(0.5, emax, 20001)
    dp = prm.D(xs)
    stat = {}
    for n in cuts:
        zz = Lsum(xs, cuts[n]['g'], cuts[n]['w'])
        for lab, tl in (('mean_tail', mean_tail4(cuts[n]['T'])), ('rvm_tail', rvm_tail(xs, cuts[n]['T']))):
            d = np.abs(dp - zz - tl)
            stat['%s_N%d' % (lab, n)] = dict(max=sg(d.max(), 4), argmax=sg(xs[int(d.argmax())], 6), median=sg(float(np.median(d)), 4), p95=sg(float(np.quantile(d, 0.95)), 4))
    log('\nprime side minus (zero side + tail) on 20001 points of [0.5, %g]: %s' % (emax, json.dumps(stat)))
    out['prime_vs_zero_grid'] = stat

    # ---- barrier numbers and the rigorous floor
    T0 = 3e12
    sig_formula = (math.log(T0 / (2 * math.pi)) + 1) / (2 * math.pi * T0)
    sig0 = 1.48e-12
    eps0 = 8 * sig0
    xlev = lambda f: 4 * math.acosh(math.sqrt(f / eps0))
    log('\nbarrier: Sigma_0 estimate (log(T0/2 pi)+1)/(2 pi T0) = %.5e (book: below 1.48e-12), epsilon_H = 8 Sigma_0 cosh^2(x/4) = %.4e cosh^2(x/4)' % (sig_formula, eps0))
    log('level where epsilon_H = 0.05: x = %.4f, e^x = %.3e' % (xlev(0.05), math.exp(xlev(0.05))))
    out['barrier'] = dict(T0=T0, Sigma0_estimate=sg(sig_formula, 6), Sigma0_bound=sg(sig0, 4), epsilon0=sg(eps0, 6), level_for_0_05=sg(xlev(0.05), 8), exp_level_for_0_05=sg(math.exp(xlev(0.05)), 5))
    claim('epsilon_H = 8 Sigma_0 cosh^2(x/4), constant', '375', '1.18e-11', sg(eps0, 4), 'agrees')
    claim('level where epsilon_H = 0.05', '382', '47.10', '%.2f' % xlev(0.05), 'agrees' if abs(xlev(0.05) - 47.10) < 5e-3 else 'differs')
    claim('e^47.1', '382', 'about 2.9e20', '%.3e (e^%.2f = %.3e)' % (math.exp(47.1), xlev(0.05), math.exp(xlev(0.05))), 'agrees' if abs(math.exp(47.1) / 2.9e20 - 1) < 0.02 else 'differs')

    fl = {}
    tol = 2e-6 if not quick else 2e-5
    h0 = 2e-3 if not quick else 1e-2
    top = min(47.1, emax if quick else 47.1)
    ranges = [('[0.05, 16]', 0.05, min(16.0, top)), ('[16, 47.1]', 16.0, 47.1), ('[0.05, 47.1]', 0.05, 47.1)] if not quick else [('[0.05, 10]', 0.05, 10.0), ('[10, 47.1]', 10.0, 47.1)]
    f_L = lambda x: Lsum(x, g, w)
    for lab, a_, b_ in ranges:
        t2 = time.time()
        ub, lbd, arg, ev = certified_min(f_L, lip, a_, b_, h0, tol)
        lbd_r = lbd - 1e-9
        est = ub + mean_tail4(T)
        # check the minimiser against mpmath with the 68 digit ordinates
        mp.mp.dps = 30
        xm = mp.mpf(arg)
        Lmp = sum(8 * mp.sin(mp.mpf(s) * xm / 2) ** 2 / mp.mpf(s) ** 2 for s in strs)
        fl[lab] = dict(min_of_L_N=sg(ub, 10), argmin=sg(arg, 8), certified_lower_bound=sg(lbd_r, 10), tolerance=sg(tol, 3), evaluations=ev, mpmath_check_at_argmin=mp.nstr(Lmp, 12), double_minus_mpmath=sg(float(f_L(np.array([arg]))[0]) - float(Lmp), 3),
                       estimate_of_D_on_min_with_smooth_tail=sg(est, 8))
        log('rigorous floor of D_on on %s from the first %d certified zeros: min L_N = %.8f at x = %.5f, certified lower bound %.8f (tolerance %.0e, %d evaluations, %.1f s); mpmath check %s; estimate with the smooth tail %.6f' % (
            lab, N, ub, arg, lbd_r, tol, ev, time.time() - t2, mp.nstr(Lmp, 10), est))
    out['floor'] = fl
    fl_main = fl[ranges[1][0]]
    floor_c = float(fl_main['certified_lower_bound'])
    rest = KAPPA1 - float(np.sum(1 / g ** 2))
    up_at = float(fl_main['min_of_L_N']) + 8 * rest
    log('upper bound of D_on at the minimiser x = %s: L_N + 8 (kappa_1 minus the first %d terms of sum gamma^-2) = %s + 8 * %.6e = %.6f, below 0.05' % (fl_main['argmin'], N, fl_main['min_of_L_N'], rest, up_at))
    out['floor_upper_bound_at_argmin'] = dict(x=fl_main['argmin'], L_N=fl_main['min_of_L_N'], remaining_sum_gamma_minus2=sg(rest, 10), upper_bound=sg(up_at, 8))
    lev_floor = xlev(floor_c)
    log('with the certified floor %.6f on [16, 47.1] the level where epsilon_H reaches the floor is x = %.4f (e^x = %.3e); with the chapter floor 0.05 it is %.4f' % (floor_c, lev_floor, math.exp(lev_floor), xlev(0.05)))
    out['floor_levels'] = {sg(f, 6): dict(level=sg(xlev(f), 6), exp_level=sg(math.exp(xlev(f)), 4)) for f in (0.05, floor_c, 0.045, 0.04, 0.03)}
    claim('floor D_on >= 0.05 on [16, 47.1] (rigorous, first %d zeros)' % N, '382', 'not computed', 'certified lower bound %.6f' % floor_c,
          'supports 0.05' if floor_c >= 0.05 else 'does NOT support 0.05 (%d zero truncation costs about %.4f; estimate with smooth tail %.5f)' % (N, mean_tail4(T), float(fl_main['estimate_of_D_on_min_with_smooth_tail'])))

    # level at which |D minus D_on| <= epsilon_H stops proving D > 0
    hh = 5e-4
    edges = np.arange(16.0, 70.0 + hh / 2, hh)
    mids = 0.5 * (edges[:-1] + edges[1:])
    Lm = Lsum(mids, g, w)
    eps_right = eps0 * np.cosh(edges[1:] / 4) ** 2
    eps_mid = eps0 * np.cosh(mids / 4) ** 2
    bad_cert = np.flatnonzero(Lm - lip * hh / 2 - 1e-9 <= eps_right)
    bad_plain = np.flatnonzero(Lm <= eps_mid)
    x_cert = float(edges[bad_cert[0]]) if len(bad_cert) else None
    x_plain = float(mids[bad_plain[0]]) if len(bad_plain) else None
    log('first level where the certified lower bound L_N(x) minus Lipschitz margin stops exceeding epsilon_H(x): x = %s; first grid point (step %g) where L_N(x) <= epsilon_H(x): x = %s' % (
        'none below 70' if x_cert is None else '%.4f' % x_cert, hh, 'none below 70' if x_plain is None else '%.4f' % x_plain))
    out['proof_level'] = dict(step=hh, certified_first_unproved=None if x_cert is None else sg(x_cert, 8), plain_grid_first=None if x_plain is None else sg(x_plain, 8),
                              meaning='D(x) > 0 is proved for 0 < x < the certified value by D >= L_N minus epsilon_H')

    # ---- detection of an off line pair
    pairs = [(100.0, 0.3, '17.3', '18.1'), (100.0, 0.2, '25.4', '27.2'), (1000.0, 0.3, '32.9', '33.5')]
    x_hi = 60.0 if not quick else 40.0
    sgrid = 0.01 if not quick else 0.02
    xo = 0.5 + sgrid * np.arange(int(round((x_hi - 0.5) / sgrid)) + 1)
    Do = {n: Lsum(xo, cuts[n]['g'], cuts[n]['w']) + mean_tail4(cuts[n]['T']) for n in cuts}

    def pair_term(x, gam, dl):
        z = complex(gam, -dl)
        return 4 * np.real(4 * np.sin(z * x / 2) ** 2 / z ** 2)
    off = {}
    # the pair term of gamma = 1000 has period 0.0063, shorter than the chapter's step 0.01, so the windowed minima and a second first negative value use a fine grid
    fstep = 5e-4 if not quick else 2e-3
    xf = 0.5 + fstep * np.arange(int(round((x_hi - 0.5) / fstep)) + 1)
    Df = Lsum(xf, g, w) + mean_tail4(T)
    wl = int(round(0.4 / fstep))
    nwin = len(xf) // wl
    centres = xf[:nwin * wl].reshape(nwin, wl).mean(axis=1)
    base_win = Df[:nwin * wl].reshape(nwin, wl).min(axis=1)
    curves = []
    log('\nD plus the pair term of an off line quadruple, first negative x on [0.5, %g] (step %g), zeros N = %d (and N = 1500)' % (x_hi, sgrid, N))
    for gam, dl, chap, pred_chap in pairs:
        pred = (2 * math.log(gam) + math.log(KAPPA1)) / dl
        res = {}
        for n in cuts:
            fvals = Do[n] + pair_term(xo, gam, dl)
            neg = np.flatnonzero(fvals < 0)
            if len(neg) == 0:
                res[str(n)] = dict(first_negative_grid=None)
                continue
            j = int(neg[0])
            fz = lambda x, n=n: float(Lsum(np.array([x]), cuts[n]['g'], cuts[n]['w'])[0] + mean_tail4(cuts[n]['T']) + pair_term(x, gam, dl))
            root = bisect_root(fz, float(xo[j - 1]), float(xo[j])) if j > 0 else float(xo[0])
            res[str(n)] = dict(first_negative_grid=sg(xo[j], 6), refined_first_zero_crossing=sg(root, 8), D_at_that_point=sg(fvals[j], 6))
        fv = Df + pair_term(xf, gam, dl)
        wmin = fv[:nwin * wl].reshape(nwin, wl).min(axis=1)
        curves.append((r'added pair $(\gamma,\delta)=(%g,%g)$' % (gam, dl), pred, wmin))
        negf = np.flatnonzero(fv < 0)
        jf = int(negf[0])
        fzf = lambda x: float(Lsum(np.array([x]), g, w)[0] + mean_tail4(T) + pair_term(x, gam, dl))
        root_f = bisect_root(fzf, float(xf[jf - 1]), float(xf[jf]))
        res['fine_grid_step_%g' % fstep] = dict(first_negative=sg(xf[jf], 8), refined_first_zero_crossing=sg(root_f, 9), D_at_that_point=sg(fv[jf], 6))
        off['gamma=%g,delta=%g' % (gam, dl)] = dict(predicted=sg(pred, 6), chapter_first_negative=chap, chapter_prediction=pred_chap, zeros=res,
                                                     first_window_with_negative_minimum_centre=None if not (wmin < 0).any() else sg(centres[int(np.flatnonzero(wmin < 0)[0])], 6))
        r_N = res[str(N)]
        r_15 = res.get('1500', res[str(N)])
        log('(gamma, delta) = (%g, %g): prediction %.2f (chapter %s); first negative on the step %g grid %s (N = %d), refined crossing %s; N = 1500: grid %s; fine grid (step %g): %s, refined %s; chapter %s' % (
            gam, dl, pred, pred_chap, sgrid, r_N.get('first_negative_grid'), N, r_N.get('refined_first_zero_crossing'), r_15.get('first_negative_grid'), fstep, sg(xf[jf], 8), sg(root_f, 8), chap))
        got = r_N.get('first_negative_grid')
        claim('first negative x, (gamma, delta) = (%g, %g)' % (gam, dl), '152', chap, '%s (fine grid %s, prediction %.1f)' % (got, sg(xf[jf], 7), pred),
              'agrees' if (got is not None and abs(float(got) - float(chap)) < 0.051 and not quick) else ('not applicable in --quick' if quick else 'differs'))
    out['detection'] = dict(grid_step=sgrid, fine_grid_step=fstep, x_range=[0.5, x_hi], windows=nwin, window_length=0.4, windows_use_fine_grid=True, pairs=off,
                            window_centres=[sg(c, 6) for c in centres], window_minima_actual=[sg(v, 5) for v in base_win],
                            window_minima_pairs=[[sg(v, 5) for v in c[2]] for c in curves])
    log('window minima of D (actual zeros, windows of length 0.4): minimum over all windows %.5f at x = %.2f, maximum %.5f' % (base_win.min(), centres[int(base_win.argmin())], base_win.max()))
    out['detection']['actual_window_minimum_overall'] = sg(base_win.min(), 6)
    out['detection']['actual_window_minimum_at'] = sg(centres[int(base_win.argmin())], 6)

    # ---- figures
    if not args.no_fig:
        figdir.mkdir(parents=True, exist_ok=True)
        figure_D(str(figdir / 'yo_D.pdf'), prm, g, w, T, N)
        figure_off(str(figdir / 'yo_off.pdf'), centres, base_win, curves)
    out['claims'] = CLAIMS
    out['seconds'] = round(time.time() - t0, 1)
    log('\ntotal %.1f s' % (time.time() - t0))
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, 'w') as fh:
        json.dump(out, fh, indent=1)
    with open(outp.with_suffix('.log'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    print('wrote', outp, 'and', outp.with_suffix('.log'))


if __name__ == '__main__':
    main()
