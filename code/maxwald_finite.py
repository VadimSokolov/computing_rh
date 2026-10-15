# maxwald_finite.py: the finite maxima stages P_N and G_N, G_N^sym of Section sec:maxwald (ch/maxwald3.tex), subsection
# "Finite maxima" and "What the finite stages look like": Theorem mw:thm:finite, Corollary mw:cor:hurwitz, Tables mw:tab:GN and
# mw:tab:sym, Figure mw:fig:hurwitz (fig/maxwald_hurwitz.pdf), and the rows "low coefficients of P_N" and "zero of G_64 near rho_1"
# of the table of numerical checks. Added for the book in October 2026 from the checks of the authors' chapter,
# misc/reviews/nick-maxwald/ (review-3 items M2 and M5: the observed rate between the heuristic and the proved one, 72 terms at
# N = 8). With P_N(q) = prod_{n <= N} (1 - q^n) prod_{m <= N} (1 - q^(2m-1)) = sum p_m q^m and D_N(w) = sum_{m >= 1} p_m m^(-w),
#   G_N(s) = (s - 1) E M_N^(s/2) / (2 (1 - 2^(1-s))) = A(s) D_N(s/2) / (1 - 2^(1-s)),  A(s) = -(s - 1)/2 Gamma(1 + s/2) pi^(-s/2),
# G_N^sym(s) = (G_N(s) + G_N(1 - s))/2, and G_N -> xi(s) = s (s - 1)/2 pi^(-s/2) Gamma(s/2) zeta(s). The script computes
#  polys:      for N in 1, 2, 3, 4, 8, 16, 32, 48, 64, 96, 128 the exact integer coefficients p_m by multiplying out the product:
#              the degree (against N(3N+1)/2), the number of nonzero terms with m >= 1 (the terms of F_N), max |p_m|, the sum
#              S = sum |p_m| m^(-1/4), the leading coefficient, sum p_m (= 0), and whether p_m = 2 (-1)^n at m = n^2 and 0
#              otherwise for m <= N (the low coefficients of theta_4); and the first j >= 1 with sum p_m m^j != 0, which is 2N
#              with value (2N)!^2/2^N (review-3 M4: F_N has a pole at s = -4N, so it is not entire);
#  chapter_values: the nine zeros of G_N printed in Table mw:tab:GN (N = 16, 32, 64) polished by the secant method from the printed
#              values, to see that each is a zero of G_N whatever "nearest rho_k" means at that N;
#  G_zeros:    the zero of G_N nearest rho_k = 1/2 + i gamma_k, k = 1, ..., 5, for a list of N up to 64 (96 and 128 with --big), by
#              the secant method on D_N(s/2) from a set of starting points, in mpmath with 30 + log10 S digits (D_N is a sum of
#              terms of size S that adds up to size 1, so double precision is not enough from N = 48 on); a zero farther than 3
#              from rho_k is not reported; the spurious zero of G_16 near 1.055 + 17.984 i;
#  sym:        the sign changes of Re G_N(1/2 + it) on [10, 44] (grid 0.002, double precision, N <= 64), the real zeros
#              polished in mpmath, and their distance to gamma_1, ..., gamma_8 (nearest zero, and the k-th zero in order);
#  sup_error:  max over s in {0.3, 0.5, 0.7} x i{0, 2, ..., 10} of |G_N - xi| and |G_N^sym - xi| (absolute and relative) and the
#              exponent c in e^(-c pi sqrt N) between consecutive N, against the heuristic c = 1 and the proved c = 1/4;
#  figure:     |Re rho_N - 1/2| of the zeros of G_N nearest rho_1, ..., rho_5 against N, with the three rates.
# Run from a folder holding zeros_1700.txt: python3 maxwald_finite.py [--small | --big]; writes maxwald_finite.json and
# fig/maxwald_hurwitz.pdf there (the figure needs matplotlib; without it only the JSON is written, and
# python3 maxwald_finite.py --figure redraws it from the JSON). --small stops at N = 16 for a quick test.
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import mpmath as mp
from scipy.special import loggamma

LOGPI = math.log(math.pi)
NS_POLY = [1, 2, 3, 4, 8, 16, 32, 48, 64, 96, 128]
NS_ZERO = [4, 6, 8, 12, 16, 24, 32, 48, 64]
NS_SYM = [8, 16, 24, 32, 48, 64]
STARTS_FULL = [(0, 0), (0.3, 0), (0.8, 0), (1.4, 0), (-0.3, 0), (0.3, 0.7), (0.3, -0.7), (0.8, 0.7), (0.8, -0.7),
               (0.3, 1.5), (0.3, -1.5), (1.0, 1.5), (1.0, -1.5)]
STARTS_NEAR = [(0, 0), (0.05, 0), (-0.05, 0), (0.4, 0)]


def pn_coeffs(N):
    deg = N * (3 * N + 1) // 2
    p = [0] * (deg + 1)
    p[0] = 1
    cur = 0
    for a in list(range(1, N + 1)) + list(range(1, 2 * N, 2)):
        for m in range(cur + a, a - 1, -1):
            p[m] -= p[m - a]
        cur += a
    assert cur == deg
    return p


class Stage:
    def __init__(self, N):
        self.N = N
        p = pn_coeffs(N)
        idx = [m for m in range(1, len(p)) if p[m] != 0]
        self.m = np.array(idx, dtype=np.int64)
        self.cint = [p[m] for m in idx]
        self.c = np.array([float(x) for x in self.cint])
        self.logm = np.log(self.m.astype(float))
        self.S14 = float(np.sum(np.abs(self.c) * self.m ** -0.25))
        self.dps = 30 + int(math.ceil(math.log10(max(self.S14, 10.0))))
        self._L = None

    def D_np(self, w):
        w = np.atleast_1d(np.asarray(w, dtype=complex))
        out = np.empty(w.shape, dtype=complex)
        B = max(1, int(2_000_000 // len(self.c)))
        for i in range(0, len(w), B):
            out[i:i + B] = np.exp(-np.multiply.outer(w[i:i + B], self.logm)) @ self.c
        return out

    def D_mp(self, w):
        mp.mp.dps = self.dps
        if self._L is None:
            self._L = [mp.log(int(m)) for m in self.m]
        w = mp.mpmathify(w)
        acc = mp.mpc(0)
        for c, L in zip(self.cint, self._L):
            acc += c * mp.exp(-w * L)
        return acc


_ST = {}


def stage(N):
    if N not in _ST:
        _ST[N] = Stage(N)
    return _ST[N]


def logA(s):
    s = np.asarray(s, dtype=complex)
    return np.log(-(s - 1) / 2) + loggamma(1 + s / 2) - (s / 2) * LOGPI


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def G_mp(st, s):
    mp.mp.dps = st.dps
    s = mp.mpmathify(s)
    A = -(s - 1) / 2 * mp.gamma(1 + s / 2) * mp.pi ** (-s / 2)
    return A * st.D_mp(s / 2) / (1 - mp.power(2, 1 - s))


def gam_list(n=10):
    for p in ('zeros_1700.txt', 'data/zeros_1700.txt', '../data/zeros_1700.txt'):
        if os.path.exists(p):
            mp.mp.dps = 60
            return [mp.mpf(x) for x in open(p).read().split()[:n]]
    raise SystemExit('zeros_1700.txt not found')


def secant(f, x0, x1, tol=1e-25, maxit=80, rad=None):
    f0, f1 = f(x0), f(x1)
    for _ in range(maxit):
        if f1 == f0:
            return None
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if rad is not None and abs(x2 - rad[0]) > rad[1]:
            return None
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(x1 - x0) < tol * max(1, abs(x1)):
            return x1
    return None


def poly_stats(N):
    p = pn_coeffs(N)
    nz = [m for m in range(1, len(p)) if p[m]]
    mx = max(abs(p[m]) for m in nz)
    low = all(p[m] == (2 * (-1) ** math.isqrt(m) if math.isqrt(m) ** 2 == m else 0) for m in range(1, N + 1))
    # review-3 M4: D_N(-j) = sum p_m m^j = (q d/dq)^j P_N at q = 1 vanishes for 1 <= j <= 2N - 1 (P_N has a zero of order 2N at
    # q = 1) and equals (2N)!^2/2^N at j = 2N, so the pole of Gamma(1 + w) at w = -2N, that is s = -4N, survives in F_N
    mom_j, mom_v = None, None
    for j in range(1, 2 * N + 1):
        v = sum(p[m] * m ** j for m in nz)
        if v != 0:
            mom_j, mom_v = j, v
            break
    return {'first_nonvanishing_moment_j_of_sum_p_m_m_pow_j': mom_j, 'it_equals_2N': mom_j == 2 * N,
            'moment_at_2N_equals_factorial_2N_squared_over_2_pow_N': mom_v == math.factorial(2 * N) ** 2 // 2 ** N,
            'degree': len(p) - 1, 'degree_formula_N_3N_plus_1_over_2': N * (3 * N + 1) // 2,
            'nonzero_terms_m_ge_1': len(nz), 'max_abs_p_m': mx, 'argmax_m': [m for m in nz if abs(p[m]) == mx][0],
            'sum_abs_p_m_m_pow_minus_quarter': float(sum(abs(p[m]) * m ** -0.25 for m in nz)),
            'p_0': p[0], 'leading_coefficient': p[-1], 'sum_p_m': sum(p),
            'low_coefficients_equal_theta4_for_m_le_N': bool(low)}


def zero_task(a):
    N, k, rho_re, rho_im, big = a
    st = stage(N)
    mp.mp.dps = st.dps
    rho = mp.mpc(rho_re, rho_im)
    starts = STARTS_FULL if N <= 32 else STARTS_NEAR
    roots = []
    for ds, dt in starts:
        x0 = rho + mp.mpc(ds, dt)
        z = secant(lambda s: st.D_mp(s / 2), x0, x0 + mp.mpc(1e-3, 1e-3), rad=(rho, 25))
        if z is not None and abs(z - rho) < 3 and not any(abs(z - r) < 1e-12 for r in roots):
            roots.append(z)
    if not roots:
        return (N, k, None)
    z = min(roots, key=lambda r: abs(r - rho))
    return (N, k, {'re': float(z.real), 'im': float(z.imag), 're_minus_half': float(z.real - mp.mpf(1) / 2),
                   'abs_re_minus_half': float(abs(z.real - mp.mpf(1) / 2)), 'distance_to_rho_k': float(abs(z - rho)),
                   'zeros_found_within_3': len(roots), 'dps': st.dps})


def sym_task(a):
    N, gam = a
    st = stage(N)
    mp.mp.dps = st.dps
    g = [mp.mpf(x) for x in gam]

    def f(t):
        s = mp.mpc(0.5, t)
        A = -(s - 1) / 2 * mp.gamma(1 + s / 2) * mp.pi ** (-s / 2)
        ph = mp.arg(A) - mp.arg(1 - mp.power(2, 1 - s))
        return (mp.expj(ph) * st.D_mp(s / 2)).real
    out = {'N': N}
    if N <= 64:
        ts = 10 + 0.002 * np.arange(17001)
        s = 0.5 + 1j * ts
        ph = np.exp(1j * (logA(s).imag - np.angle(1 - 2 ** (1 - s))))
        lv = (ph * st.D_np(s / 2)).real
        chg = np.nonzero(np.sign(lv[1:]) != np.sign(lv[:-1]))[0]
        t_est = ts[chg] + 0.002 * lv[chg] / (lv[chg] - lv[chg + 1])
        out['sign_changes_on_10_44'] = int(len(chg))
        roots = []
        for t0 in t_est:
            r = secant(f, mp.mpf(float(t0)), mp.mpf(float(t0)) + mp.mpf(1e-7))
            roots.append(r if r is not None else mp.mpf(float(t0)))
    else:
        out['sign_changes_on_10_44'] = None
        roots = []
        for k in (1, 2, 3, 5, 8):
            r = secant(f, g[k - 1], g[k - 1] + mp.mpf(1e-3))
            if r is not None:
                roots.append(r)
    out['real_zeros'] = [float(r) for r in roots]
    out['err_nearest'] = {str(k): float(min(abs(r - g[k - 1]) for r in roots)) for k in range(1, 9)} if roots else {}
    srt = sorted(roots)
    out['err_kth_in_order'] = {str(k): float(abs(srt[k - 1] - g[k - 1])) for k in range(1, 9) if len(srt) >= k} if N <= 64 else {}
    return out


def sup_task(N):
    st = stage(N)
    mp.mp.dps = st.dps
    pts = [(sg, t) for sg in (0.3, 0.5, 0.7) for t in (0, 2, 4, 6, 8, 10)]
    G = {pt: G_mp(st, mp.mpc(pt[0], pt[1])) for pt in pts}
    e = er = es = esr = 0
    for sg, t in pts:
        x = xi(mp.mpc(sg, t))
        gs = (G[(sg, t)] + mp.conj(G[(round(1 - sg, 10), t)])) / 2
        e = max(e, abs(G[(sg, t)] - x))
        er = max(er, abs(G[(sg, t)] - x) / abs(x))
        es = max(es, abs(gs - x))
        esr = max(esr, abs(gs - x) / abs(x))
    return {'N': N, 'sup_abs_G_minus_xi': float(e), 'sup_rel_G_minus_xi': float(er),
            'sup_abs_Gsym_minus_xi': float(es), 'sup_rel_Gsym_minus_xi': float(esr), 'dps': st.dps}


CHAPTER_G = [(16, 0.51577, 14.09915), (16, 0.95703, 20.91901), (16, 1.56809, 34.68128), (32, 0.49948, 14.13483),
             (32, 0.48159, 21.02790), (32, 1.02797, 32.22350), (64, 0.50000052, 14.1347256),
             (64, 0.50003418, 21.0220191), (64, 0.51254, 32.92764)]


def chapter_task(a):
    """Polish a zero of G_N listed in Table tab:GN of the chapter (secant start at the printed value)."""
    N, re, im = a
    st = stage(N)
    mp.mp.dps = st.dps
    x0 = mp.mpc(re, im)
    z = secant(lambda s: st.D_mp(s / 2), x0, x0 + mp.mpc(1e-4, 1e-4), rad=(x0, 1))
    if z is None:
        return {'N': N, 'chapter': [re, im], 'ours': None}
    return {'N': N, 'chapter': [re, im], 'ours': [float(z.real), float(z.imag)],
            'abs_diff_to_chapter': float(abs(z - x0)), 'abs_re_minus_half': float(abs(z.real - mp.mpf(1) / 2))}


def spurious_task(_):
    st = stage(16)
    mp.mp.dps = st.dps
    x0 = mp.mpc(1.05, 17.98)
    z = secant(lambda s: st.D_mp(s / 2), x0, x0 + mp.mpc(1e-3, 1e-3), rad=(x0, 5))
    return {'re': float(z.real), 'im': float(z.imag)} if z is not None else None


def make_figure(res):
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
    ls = ['-', '--', ':', '-.', (0, (6, 1.5, 1, 1.5, 1, 1.5))]
    gr = ['0.0', '0.2', '0.35', '0.5', '0.62']
    mk = ['o', 's', '^', 'D', 'v']
    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    Z = res['G_zeros']
    Ns = sorted({int(n) for n in Z})
    d1 = {}
    for k in range(1, 6):
        xs = [n for n in Ns if Z[str(n)].get(str(k))]
        ys = [Z[str(n)][str(k)]['abs_re_minus_half'] for n in xs]
        ax.plot(xs, ys, ls=ls[k - 1], color=gr[k - 1], marker=mk[k - 1], ms=3.5, lw=1.2, label=f'zero {k}')
        if k == 1:
            d1 = dict(zip(xs, ys))
    if d1:
        N0 = max(d1)
        grid = np.linspace(min(d1), max(d1), 100)
        for c, st_, lab in ((1.0, (0, (1, 2.5)), r'$e^{-\pi\sqrt{N}}$, heuristic'),
                            (0.75, (0, (4, 2, 1, 2)), r'$e^{-0.75\pi\sqrt{N}}$, observed'),
                            (0.25, (0, (3, 3)), r'$e^{-\pi\sqrt{N}/4}$, proved')):
            ax.plot(grid, d1[N0] * np.exp(-c * math.pi * (np.sqrt(grid) - math.sqrt(N0))), ls=st_, color='k', lw=0.8,
                    alpha=0.6, label=lab)
    ax.set_xscale('log', base=2)
    ax.set_yscale('log')
    tk = [n for n in Ns if n in (4, 8, 16, 32, 64, 96, 128)]
    ax.set_xticks(tk)
    ax.set_xticklabels([str(n) for n in tk])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel(r'$N$')
    ax.set_ylabel(r'$|\mathrm{Re}\,\rho_N-1/2|$')
    ax.legend(fontsize=7, frameon=False, ncol=2, loc='upper right')
    fig.tight_layout()
    fig.savefig('fig/maxwald_hurwitz.pdf')
    print('wrote fig/maxwald_hurwitz.pdf')


def main():
    a = sys.argv[1:]
    if a and a[0] == '--figure':
        make_figure(json.load(open('maxwald_finite.json')))
        return
    small, big = '--small' in a, '--big' in a
    t0 = time.time()
    gam = gam_list(10)
    gams = [mp.nstr(x, 30) for x in gam]
    res = {'working_precision': 'mpmath, dps = 30 + ceil(log10 S) per N (S = sum |p_m| m^(-1/4)); grids and the pair of '
           'gamma factors in double precision (numpy, scipy loggamma); gamma_k from zeros_1700.txt (70 digits)',
           'gamma_1_to_10': gams}
    res['polys'] = {str(N): poly_stats(N) for N in ([n for n in NS_POLY if n <= 16] if small else NS_POLY)}
    nz = [n for n in NS_ZERO if n <= (16 if small else 64)] + ([96, 128] if big else [])
    ns = [n for n in NS_SYM if n <= (16 if small else 64)] + ([96, 128] if big else [])
    tasks_z = [(n, k, float(0.5), float(gam[k - 1]), big) for n in nz for k in range(1, 6)]
    tasks_z.sort(key=lambda x: -x[0])
    tasks_s = [(n, gams) for n in ns]
    tasks_u = list(ns)
    procs = 1 if small else min(8, os.cpu_count() or 1)
    with Pool(procs) as pool:
        r_sp = pool.map_async(spurious_task, [0])
        r_z = pool.map_async(zero_task, tasks_z, chunksize=1)
        r_s = pool.map_async(sym_task, tasks_s, chunksize=1)
        r_u = pool.map_async(sup_task, tasks_u, chunksize=1)
        r_c = pool.map_async(chapter_task, [c for c in CHAPTER_G if c[0] <= (16 if small else 64)], chunksize=1)
        res['G_zeros_polished_from_the_chapters_printed_values'] = r_c.get()
        res['G16_spurious_zero'] = r_sp.get()[0]
        Z = {}
        for N, k, d in r_z.get():
            Z.setdefault(str(N), {})[str(k)] = d
        res['G_zeros'] = Z
        res['sym'] = {str(d['N']): d for d in r_s.get()}
        res['sup_error'] = {str(d['N']): d for d in r_u.get()}
    ks = sorted(int(n) for n in res['sup_error'])
    rate = []
    for n1, n2 in zip(ks[:-1], ks[1:]):
        e1, e2 = res['sup_error'][str(n1)], res['sup_error'][str(n2)]
        den = math.pi * (math.sqrt(n2) - math.sqrt(n1))
        rate.append({'N1': n1, 'N2': n2,
                     'c_sup_abs_G': math.log(e1['sup_abs_G_minus_xi'] / e2['sup_abs_G_minus_xi']) / den,
                     'c_sup_abs_Gsym': math.log(e1['sup_abs_Gsym_minus_xi'] / e2['sup_abs_Gsym_minus_xi']) / den})
    sy = res['sym']
    for n1, n2 in zip(ns[:-1], ns[1:]):
        e1, e2 = sy[str(n1)]['err_nearest'].get('1'), sy[str(n2)]['err_nearest'].get('1')
        if e1 and e2:
            for r in rate:
                if r['N1'] == n1 and r['N2'] == n2:
                    r['c_error_at_gamma_1'] = math.log(e1 / e2) / (math.pi * (math.sqrt(n2) - math.sqrt(n1)))
    res['rate_exponent_c_in_exp_minus_c_pi_sqrt_N'] = rate
    if '32' in sy and '64' in sy:
        e32, e64 = sy['32']['err_nearest']['1'], sy['64']['err_nearest']['1']
        res['gamma_1_error_ratio_N32_to_N64'] = {'observed': e32 / e64, 'heuristic_exp_pi_8_minus_sqrt32': math.exp(math.pi * (8 - math.sqrt(32))),
                                                 'proved_fourth_root': math.exp(math.pi * (8 - math.sqrt(32)) / 4)}
    res['seconds'] = round(time.time() - t0, 1)
    json.dump(res, open('maxwald_finite.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('G_zeros', 'sym')}, indent=1)[:5000])
    make_figure(res)


if __name__ == '__main__':
    main()
