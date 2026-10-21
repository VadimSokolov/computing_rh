# yo_damped.py: Table yo:tab:lam, the floor row, Figure yo:fig:lam and the check of Remark yo:rem (alpha = 1, a = 0.8) of Chapter 17 (Yoshida's Hermitian Forms, Weil Positivity and the Thorin Program), Section yo:sec:damped (damped forms and sharp coercivity).
# Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were neither read nor run). Uses the shared module code/yo_matrix.py (the prime side matrix of the damped form on the modes chi_n(u) = (2a)^(-1/2) e^(i pi n u/a), |n| <= N).
# Parts (each writes its key into the JSON, which is updated in place, so partial runs add up). ladder: the least eigenvalue lambda_N(a, alpha) of the matrix on |n| <= N for a in {0.35, 0.8, 1.5, 2.5, 4}, alpha in {0.5, 0.55, 0.6, 0.75, 1}, double precision, N in the ladder 12, 24, 48, 96, 192, 384 and the sizes max(12, 12a), max(24, 24a) (the sizes the chapter's script takes), ceil(12a), ceil(24a) (the sizes the chapter's caption states); the rate p of the algebraic convergence lambda_N - lambda_inf ~ N^(-p) and the converged value lambda_inf from Richardson extrapolation on 96, 192, 384, with the uncertainty from 48, 96, 192; the comparison with the printed table; the size actually needed for 5e-8. half: alpha = 1/2 in extended precision (mpmath, the digits raised until two precisions agree to 8 digits), the sizes as above, and the decay in a against exp(-c L e^L) with L = 2a. floors: b_alpha = xi'(alpha)/xi(alpha) at 30 digits, the row 2 b_alpha, and the bounds of Theorem yo:thm:coercive (iv) with kappa_1 = (log xi)''(1/2)/2 and kappa_2 from the first 1700 zeros. line243: at alpha = 1, a = 0.8 the entry (0, 0) of the matrix against the symbol integral (2/pi) int_0^oo (2/a) sin^2(a tau)/tau^2 v_1(tau) dtau, v_1 = Re xi'/xi(1 + i tau), computed from the first 1700 zeros with a smooth tail on [0, T], plus the smooth part of v_1 beyond T (digamma), plus a model of the zeta'/zeta part beyond T from its Dirichlet series (the prime 5, with log 5 - 2a = 0.0094, makes the tail converge like 1/(0.0094 T^2)). cross: the conventions and cross checks for the shared module (Zhu's table zh:tab:hp and the deconvolution chapter's numbers at y = 0, informational). figure: fig/yo_lam.pdf from the JSON (left: the excess lambda_min/(2 b_alpha) - 1 over the floor for alpha = 0.55, 0.6, 0.75, 1; right: the least eigenvalue at alpha = 1/2 on the sizes of the table and the double precision floor; black and white).
# Usage (repository root): timeout 60 python3 -I code/yo_damped.py --part ladder|half|floors|line243|cross|figure [--quick] [--out FILE] [--procs P]; --quick restricts half to a in {0.35, 0.8} with N <= 48 (a smoke test of about 30 s). The full run is for the Hopper cluster:
# bash misc/tools/hopper_run.sh -c 8 -t 240 -g yo_damped.json -g yo_damped.log code/yo_damped.py code/yo_matrix.py data/zeros_1700.txt -- --part all --procs 8 --out yo_damped.json
# Output: data/yo_damped.json and data/yo_damped.log (numbers as strings: double precision values to 12 decimals, extended precision values to 25 digits).
import argparse
import json
import math
import os
import platform
import socket
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mpmath
import numpy as np

import yo_matrix as ym

A_TABLE = [0.35, 0.8, 1.5, 2.5, 4.0]
ALPHAS = [0.5, 0.55, 0.6, 0.75, 1.0]
PRINTED = {0.35: [1.20e-3, 0.0066203, 0.0120264, 0.0281575, 0.0547629],
           0.8: [None, 0.0047366, 0.0094674, 0.0236499, 0.0472537],
           1.5: [None, 0.0046491, 0.0092978, 0.0232406, 0.0464627],
           2.5: [None, 0.0046305, 0.0092609, 0.0231498, 0.0462848],
           4.0: [None, 0.0046246, 0.0092491, 0.0231206, 0.0462270]}
PRINTED_BOUND = {0.8: '1e-15', 1.5: '1e-13', 2.5: '1e-13', 4.0: '1e-12'}
FLOOR_PRINTED = {0.55: 0.0046210, 0.6: 0.0092418, 0.75: 0.0231027, 1.0: 0.0461914}


class Log:
    def __init__(self, path):
        self.fh = open(path, 'a')

    def __call__(self, msg=''):
        print(msg, flush=True)
        self.fh.write(msg + '\n')
        self.fh.flush()


def default_path(stem):
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    if os.path.basename(here) == 'code' and os.path.isdir(os.path.join(root, 'data')):
        return os.path.join(root, 'data', stem)
    return os.path.join(os.getcwd(), stem)


def find_zeros():
    here = os.path.dirname(os.path.abspath(__file__))
    for p in (os.path.join(os.path.dirname(here), 'data', 'zeros_1700.txt'), os.path.join(here, 'zeros_1700.txt'), 'zeros_1700.txt'):
        if os.path.exists(p):
            return np.array([float(s.split()[0]) for s in open(p) if s.strip()])
    return None


def load(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save(path, res):
    with open(path, 'w') as fh:
        json.dump(res, fh, indent=1)


def stamp(res, part, seconds, mode):
    res.setdefault('meta', {})
    res['meta'][part] = dict(seconds=round(seconds, 1), host=socket.gethostname(), mode=mode, date=time.strftime('%Y-%m-%d %H:%M'), python=platform.python_version(),
                             numpy=np.__version__, mpmath=mpmath.__version__)


def fmt(x, n=12):
    return '%.*f' % (n, x)


# ---------------------------------------------------------------------------------------------------------------------------------
def sizes_for(a):
    return sorted({max(12, int(math.ceil(12 * a))), max(24, int(math.ceil(24 * a))), int(math.ceil(12 * a)), int(math.ceil(24 * a)), 12, 24, 48, 96, 192, 384})


def richardson(ns, vals):
    """lambda_N - lambda_inf ~ c N^(-p) from the last three sizes (ratio 2): the rate p and the extrapolated value (None if not monotone)."""
    (n1, n2, n3), (v1, v2, v3) = ns[-3:], vals[-3:]
    d1, d2 = v1 - v2, v2 - v3
    if d1 > 0 and d2 > 0 and d1 > d2:
        p = math.log(d1 / d2) / math.log(2.0)
        return p, v3 - d2 / (2 ** p - 1)
    return None, None


def part_ladder(res, log, grid=False):
    t0 = time.time()
    out = res.setdefault('ladder', {})
    for a in A_TABLE:
        ns = sizes_for(a)
        out[str(a)] = {}
        for al in ALPHAS:
            vals = [float(ym.least_eigenvalue(a, al, basis='fourier', size=N)) for N in ns]
            ent = dict(N=ns, lam=[fmt(v) for v in vals])
            if al > 0.5:
                v = dict(zip(ns, vals))
                p1, ext1 = richardson([48, 96, 192], [v[48], v[96], v[192]])
                p2, ext2 = richardson([96, 192, 384], [v[96], v[192], v[384]])
                if ext2 is not None:
                    ent.update(rate=fmt(p2, 3), converged=fmt(ext2), uncertainty=('%.1e' % abs(ext2 - ext1)) if ext1 is not None else 'n/a')
                    need = [n for n, x in zip(ns, vals) if abs(x - ext2) < 5e-8]
                    ent['N_needed_5e-8'] = need[0] if need else '>%d' % ns[-1]
            out[str(a)][str(al)] = ent
        log('ladder a=%g: sizes %s done (%.1f s)' % (a, ns, time.time() - t0))
    # the comparison with the printed table
    cmp = res.setdefault('table_comparison', {})
    for a in A_TABLE:
        ns = out[str(a)]['0.5']['N']
        n1, n2 = max(12, int(math.ceil(12 * a))), max(24, int(math.ceil(24 * a)))
        n1c, n2c = int(math.ceil(12 * a)), int(math.ceil(24 * a))
        cmp[str(a)] = {}
        for k, al in enumerate(ALPHAS):
            e = out[str(a)][str(al)]
            lam = lambda N: float(e['lam'][ns.index(N)])
            row = dict(N_script_small=n1, N_script_big=n2, lam_N_script_small=fmt(lam(n1)), lam_N_script_big=fmt(lam(n2)),
                       lam_N_caption_small=fmt(lam(n1c)), lam_N_caption_big=fmt(lam(n2c)), N_caption=[n1c, n2c], printed=PRINTED[a][k])
            if al > 0.5:
                row['converged'] = e.get('converged')
                row['printed_minus_script_big'] = '%.2e' % (PRINTED[a][k] - lam(n2))
                row['printed_minus_caption_big'] = '%.2e' % (PRINTED[a][k] - lam(n2c))
                row['printed_minus_converged'] = ('%.2e' % (PRINTED[a][k] - float(e['converged']))) if 'converged' in e else 'n/a'
                row['small_to_big_change_script'] = '%.2e' % (lam(n1) - lam(n2))
            cmp[str(a)][str(al)] = row
    log('table (alpha > 1/2): printed | lambda at N = max(24, 24a) | lambda at N = ceil(24a) | converged (Richardson) | size needed for 5e-8')
    for a in A_TABLE:
        for al in ALPHAS[1:]:
            e, r = out[str(a)][str(al)], cmp[str(a)][str(al)]
            log('  a=%-4g alpha=%-4g printed %.7f | N=%-3d %s | N=%-3d %s | conv %s (rate %s, +-%s) | need N=%s' % (
                a, al, r['printed'], r['N_script_big'], r['lam_N_script_big'][:10], r['N_caption'][1], r['lam_N_caption_big'][:10], e.get('converged', 'n/a')[:10],
                e.get('rate'), e.get('uncertainty'), e.get('N_needed_5e-8')))
    log('alpha = 1/2 in double precision (noise, not values): ' + '; '.join('a=%g: %s' % (a, ' '.join('%.1e' % float(v) for v in out[str(a)]['0.5']['lam'][:3])) for a in A_TABLE))
    if grid:
        g = res.setdefault('figure_grid', {})
        for a in [0.2, 0.35, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]:
            N = max(96, int(math.ceil(24 * a)))
            g[str(a)] = {str(al): fmt(float(ym.least_eigenvalue(a, al, basis='fourier', size=N))) for al in ALPHAS[1:]}
            g[str(a)]['N'] = N
    stamp(res, 'ladder', time.time() - t0, 'double precision')


def part_legendre(res, log):
    """Independent basis: the Legendre matrix (M = 24, 48, 96) at alpha = 0.6 and 1, with its own Richardson extrapolation, against the Fourier converged value."""
    t0 = time.time()
    out = res.setdefault('legendre_check', {})
    for a in A_TABLE:
        out[str(a)] = {}
        for al in (0.6, 1.0):
            vals = [float(ym.least_eigenvalue(a, al, basis='legendre', size=M)) for M in (24, 48, 96)]
            p, ext = richardson([24, 48, 96], vals)
            fe = res.get('ladder', {}).get(str(a), {}).get(str(al), {}).get('converged')
            out[str(a)][str(al)] = dict(M=[24, 48, 96], lam=[fmt(v) for v in vals], rate=fmt(p, 3) if p else 'n/a', converged=fmt(ext) if ext else 'n/a', fourier_converged=fe,
                                        difference=('%.1e' % (ext - float(fe))) if (ext and fe) else 'n/a')
            log('legendre a=%g alpha=%g: M=24,48,96: %s, rate %s, converged %s against Fourier %s' % (a, al, ' '.join(fmt(v, 9) for v in vals), out[str(a)][str(al)]['rate'],
                                                                                           out[str(a)][str(al)]['converged'][:11], (fe or 'n/a')[:11]))
    stamp(res, 'legendre', time.time() - t0, 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
def half_case(args):
    a, N, dps_list = args
    t0 = time.time()
    prev, rel = None, None
    for dps in dps_list:
        lam = ym.least_eigenvalue(str(a), 0.5, basis='fourier', size=N, dps=dps)
        if prev is not None:
            rel = abs(lam - prev) / abs(prev)
            if rel < mpmath.mpf(10) ** -8:
                break
        prev = lam
    return dict(a=a, N=N, dps=dps, lam=mpmath.nstr(lam, 25), rel_change_from_previous_dps=mpmath.nstr(rel, 3) if rel is not None else 'n/a',
                seconds=round(time.time() - t0, 1))


def part_half(res, log, quick, procs, path, a_arg=None, max_n=None, min_n=None):
    t0 = time.time()
    out = res.setdefault('half', {})
    if a_arg:
        a_list = [float(x) for x in a_arg.split(',')]
    else:
        a_list = [0.35, 0.8] if quick else [0.25, 0.35, 0.5, 0.65, 0.8, 1.0, 1.5, 2.5, 4.0]
    if quick and max_n is None:
        max_n = 48
    cases = []
    for a in a_list:
        ns = {max(12, int(math.ceil(12 * a))), max(24, int(math.ceil(24 * a))), int(math.ceil(12 * a)), int(math.ceil(24 * a))}
        if a < 0.5:
            ns |= {4, 5, 8, 9}
        if a <= 1:
            ns |= {48, 96}
        if a <= 0.5:
            ns |= {192}
        for N in sorted(ns):
            if (max_n is None or N <= max_n) and (min_n is None or N >= min_n):
                cases.append((a, N, [50, 80, 130, 220, 380, 640, 1100]))
    cases.sort(key=lambda c: -(c[1] ** 3))

    def done(r):
        out['%g_%d' % (r['a'], r['N'])] = r
        log('half a=%g N=%d dps=%d lambda=%s (%.1f s)' % (r['a'], r['N'], r['dps'], mpmath.nstr(mpmath.mpf(r['lam']), 10), r['seconds']))
        save(path, res)

    if procs > 1:
        import multiprocessing
        with multiprocessing.Pool(procs) as pool:
            for r in pool.imap_unordered(half_case, cases):
                done(r)
    else:
        for c in cases:
            done(half_case(c))
    # decay against exp(-c L e^L), L = 2a: c = -ln(lambda)/(L e^L) on the largest size computed for each a
    dec = res.setdefault('half_decay', {})
    for a in sorted({v['a'] for v in out.values()}):
        big = max((v for v in out.values() if abs(v['a'] - a) < 1e-12), key=lambda v: v['N'])
        L = 2 * a
        c = -mpmath.log(mpmath.mpf(big['lam'])) / (L * mpmath.e ** L)
        dec[str(a)] = dict(N_largest=big['N'], lam=big['lam'], L=L, c=mpmath.nstr(c, 6), log10_lam=mpmath.nstr(mpmath.log10(mpmath.mpf(big['lam'])), 8))
        log('decay a=%g L=%g: lambda_N(N=%d) = %s, c = -ln(lambda)/(L e^L) = %s' % (a, L, big['N'], mpmath.nstr(mpmath.mpf(big['lam']), 8), dec[str(a)]['c']))
    stamp(res, 'half', time.time() - t0, 'local' if (quick or a_arg) else 'full')


# ---------------------------------------------------------------------------------------------------------------------------------
def part_floors(res, log, gam):
    t0 = time.time()
    fl = res.setdefault('floors', {})
    with mpmath.workdps(40):
        s = mpmath.mpf(1) / 2
        logxi2 = (-1 / s ** 2 - 1 / (s - 1) ** 2 + mpmath.psi(1, s / 2) / 4
                  + (mpmath.zeta(s, derivative=2) * mpmath.zeta(s) - mpmath.zeta(s, derivative=1) ** 2) / mpmath.zeta(s) ** 2)
        kappa1 = logxi2 / 2
        fl['kappa1'] = mpmath.nstr(kappa1, 15)
        fl['V_xi_log_xi_second_derivative_at_half'] = mpmath.nstr(logxi2, 15)
        k2 = float((gam ** -4.0).sum()) if gam is not None else None
        fl['kappa2_from_1700_zeros'] = repr(k2)
        log('kappa_1 = (log xi)\'\'(1/2)/2 = %s (printed 0.0231049931); kappa_2 from 1700 zeros = %s' % (fl['kappa1'], fl['kappa2_from_1700_zeros']))
        for al in (0.55, 0.6, 0.75, 1.0):
            b = ym.b_alpha(al, 30)
            c = mpmath.mpf(al) - mpmath.mpf(1) / 2
            row = dict(b_alpha=mpmath.nstr(b, 20), two_b=mpmath.nstr(2 * b, 20), printed=FLOOR_PRINTED[al], agrees_to_printed_digits=bool(abs(2 * b - FLOOR_PRINTED[al]) < 6e-8),
                       upper_4k1c=mpmath.nstr(4 * kappa1 * c, 12))
            if k2 is not None:
                row['lower_4k1c_minus_4k2c3'] = mpmath.nstr(4 * kappa1 * c - 4 * k2 * c ** 3, 12)
                row['bounds_hold'] = bool(4 * kappa1 * c - 4 * k2 * c ** 3 <= 2 * b <= 4 * kappa1 * c)
            fl[str(al)] = row
            log('floor alpha=%g: 2 b = %s (printed %.7f), bounds %s .. %s' % (al, row['two_b'][:12], FLOOR_PRINTED[al], row.get('lower_4k1c_minus_4k2c3'), row['upper_4k1c']))
    stamp(res, 'floors', time.time() - t0, 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
def zero_sum_symbol(alpha, tau, gam):
    """Re xi'/xi(alpha + i tau) = sum over zeros (Lorentzians of width c = alpha - 1/2) for the first K zeros, plus the smooth tail with the density
    theta'(g)/pi and the boundary term f(g_K) (theta(g_K)/pi + 1 - K) from N(g) = theta(g)/pi + 1 + S(g)."""
    from scipy.special import digamma, loggamma
    c = alpha - 0.5
    K = len(gam)
    f = lambda t, g: c / (c * c + (t - g) ** 2) + c / (c * c + (t + g) ** 2)
    tau = np.asarray(tau, dtype=float)
    out = np.zeros_like(tau)
    for i in range(0, len(tau), 2000):
        t = tau[i:i + 2000, None]
        out[i:i + 2000] = f(t, gam[None, :]).sum(axis=1)
    gK = gam[-1]
    theta = lambda g: (loggamma(0.25 + 0.5j * g).imag - 0.5 * g * math.log(math.pi))
    xs, ws = np.polynomial.legendre.leggauss(32)
    edges = gK * 2.0 ** np.arange(0, 24)
    pts = np.concatenate([(e0 + e1) / 2 + (e1 - e0) / 2 * xs for e0, e1 in zip(edges[:-1], edges[1:])])
    wts = np.concatenate([(e1 - e0) / 2 * ws for e0, e1 in zip(edges[:-1], edges[1:])])
    dens = (digamma(0.25 + 0.5j * pts).real - math.log(math.pi)) / (2 * math.pi)
    for i in range(0, len(tau), 2000):
        t = tau[i:i + 2000, None]
        out[i:i + 2000] += (f(t, pts[None, :]) * (dens * wts)[None, :]).sum(axis=1)
    big = edges[-1]
    out += 2 * c / (2 * math.pi) * (math.log(big / (2 * math.pi)) + 1) / big
    out -= f(tau, gK) * (K - theta(gK) / math.pi - 1)
    return out


def smooth_symbol_1(tau):
    """The part of Re xi'/xi(1 + i tau) other than zeta'/zeta: 1/(1 + tau^2) - log(pi)/2 + Re psi((1 + i tau)/2)/2."""
    from scipy.special import digamma
    return 1 / (1 + tau ** 2) - 0.5 * math.log(math.pi) + 0.5 * digamma(0.5 + 0.5j * tau).real


def part_line243(res, log, gam, quick):
    from scipy.special import sici
    t0 = time.time()
    a = 0.8
    ent = res.setdefault('line243', {})
    m_d = float(ym.build_matrix(a, 1.0, basis='fourier', size=0)[0, 0])
    m_e = ym.build_matrix('0.8', '1', basis='fourier', size=0, dps=40)[0, 0]
    ent['prime_side_entry_double'] = fmt(m_d)
    ent['prime_side_entry_40_digits'] = mpmath.nstr(m_e, 25)
    log('line 243: prime side entry (0, 0) at alpha = 1, a = 0.8: double %s, 40 digits %s (chapter 0.114058)' % (fmt(m_d, 10), mpmath.nstr(m_e, 18)))
    if gam is None:
        log('zeros_1700.txt not found: the symbol integral is skipped')
        return
    # the zero sum symbol against the exact one
    pts = [1.0, 5.0, 14.2, 20.0, 50.0, 100.0, 200.0, 399.7, 800.0, 1200.0]
    zs = zero_sum_symbol(1.0, np.array(pts), gam)
    ex = np.array([float(ym.symbol_exact(1, t, 20)) for t in pts])
    ent['symbol_zero_sum_vs_exact_max_diff'] = '%.2e' % np.abs(zs - ex).max()
    ent['symbol_check_points'] = pts
    log('  v_1 from 1700 zeros plus smooth tail against the exact xi\'/xi at %d points up to 1200: max difference %.2e' % (len(pts), np.abs(zs - ex).max()))
    w2 = lambda t: (2 / a) * np.sin(a * t) ** 2 / t ** 2
    xs, ws = np.polynomial.legendre.leggauss(10)

    def panels(lo, hi, width, xs, ws):
        edges = np.linspace(lo, hi, int(round((hi - lo) / width)) + 1)
        return (np.concatenate([(e0 + e1) / 2 + (e1 - e0) / 2 * xs for e0, e1 in zip(edges[:-1], edges[1:])]),
                np.concatenate([(e1 - e0) / 2 * ws for e0, e1 in zip(edges[:-1], edges[1:])]))

    t_max = 400.0 if quick else 1200.0
    cuts = [50.0, 100.0, 200.0, 400.0] + ([800.0, 1200.0] if not quick else [])
    # the integral of w2 v_1 over [0, 1200] in pieces (cumulative at the cuts)
    cum, lo, tot = {}, 0.0, 0.0
    for cut in cuts:
        tt, ww = panels(lo, cut, 0.5, xs, ws)
        tot += float((ww * w2(tt) * zero_sum_symbol(1.0, tt, gam)).sum())
        cum[cut] = tot * 2 / math.pi
        lo = cut
    # smooth tail beyond T: numerical to T2, then 1/(pi a) (log(T2/(2 pi)) + 1)/T2
    T2 = 40000.0
    xs6, ws6 = np.polynomial.legendre.leggauss(6)
    tail_s = {}
    for cut in cuts:
        tt, ww = panels(cut, T2, 1.0, xs6, ws6)
        tail_s[cut] = 2 / math.pi * float((ww * w2(tt) * smooth_symbol_1(tt)).sum()) + (math.log(T2 / (2 * math.pi)) + 1) / (math.pi * a * T2)
    # model of the zeta'/zeta part beyond T: -sum Lambda(n)/n cos(tau log n); J_n(T) = (4/(pi a)) (K(b)/2 - (K(2a + b) + K(2a - b))/4)
    qs, ps = ym.prime_powers(math.log(2000000.0))
    q, p = np.array(qs, dtype=float), np.array(ps, dtype=float)
    b, lam = np.log(q), np.log(p)

    def Kint(om, T):
        om = np.abs(om)
        return np.cos(om * T) / T - om * (math.pi / 2 - sici(om * T)[0])

    tail_z = {}
    for cut in cuts:
        J = 4 / (math.pi * a) * (Kint(b, cut) / 2 - (Kint(2 * a + b, cut) + Kint(2 * a - b, cut)) / 4)
        tail_z[cut] = float(-(lam / q * J).sum())
    rows = {}
    for cut in cuts:
        lead = (math.log(cut / (2 * math.pi)) + 1) / (math.pi * a * cut)
        rows[str(cut)] = dict(integral_to_T=fmt(cum[cut], 9), smooth_tail=fmt(tail_s[cut], 9), zeta_tail_model=fmt(tail_z[cut], 9),
                              leading_asymptotic_tail=fmt(lead, 9), integral_plus_leading_tail=fmt(cum[cut] + lead, 9),
                              smooth_tail_minus_leading='%.2e' % (tail_s[cut] - lead),
                              chapter_recipe_without_zeta_tail=fmt(cum[cut] + tail_s[cut], 9), total_with_zeta_tail=fmt(cum[cut] + tail_s[cut] + tail_z[cut], 9),
                              total_minus_prime_side='%.2e' % (cum[cut] + tail_s[cut] + tail_z[cut] - m_d),
                              recipe_minus_prime_side='%.2e' % (cum[cut] + tail_s[cut] - m_d))
        log('  T=%-6g integral to T %s + smooth tail %s = %s (prime side minus this: %s); zeta tail model %s; total %s (minus prime side %s)' % (
            cut, fmt(cum[cut], 7), fmt(tail_s[cut], 7), fmt(cum[cut] + tail_s[cut], 7), '%.2e' % (m_d - cum[cut] - tail_s[cut]), fmt(tail_z[cut], 7),
            fmt(cum[cut] + tail_s[cut] + tail_z[cut], 7), '%.1e' % (cum[cut] + tail_s[cut] + tail_z[cut] - m_d)))
    ent['by_cut'] = rows
    ent['note'] = ('v_1 on [0, T] from the first 1700 zeros plus the smooth tail of the zero density; the beyond T smooth part by digamma to 4e4 plus an analytic '
                   'remainder; the zeta tail model sums -Lambda(n)/n J_n(T) over n <= 2e6; log 5 - 2a = 0.0094 (e^1.6 = 4.953 < 5) is the near resonance')
    stamp(res, 'line243', time.time() - t0, 'local' if quick else 'full')


# ---------------------------------------------------------------------------------------------------------------------------------
def part_cross(res, log):
    t0 = time.time()
    ent = res.setdefault('cross', {})
    zhu = {}
    for L, K, even, odd in (('0.5', 8, '1.336e-6', '2.51e-4'), ('0.5', 14, '1.082e-6', '2.39e-4'), ('0.8', 8, '5.49e-15', '2.67e-12'), ('0.8', 14, '4.127e-17', '3.86e-14')):
        ev = ym.least_eigenvalue(L, 0.5, basis='cossin', size=K + 1, sector='even', dps=60)
        od = ym.least_eigenvalue(L, 0.5, basis='cossin', size=K, sector='odd', dps=60)
        zhu['L=%s modes=%d' % (L, K)] = dict(even=mpmath.nstr(ev, 8), odd=mpmath.nstr(od, 8), zhu_even=even, zhu_odd=odd, calls='cossin even size=%d, odd size=%d (cos modes 0..K, sin modes 1..K)' % (K + 1, K))
        log('Zhu table L=%s K=%d: even %s (printed %s), odd %s (printed %s)' % (L, K, mpmath.nstr(ev, 5), even, mpmath.nstr(od, 5), odd))
    ent['zhu_table_cossin'] = zhu
    dec = {}
    for a, M, printed in (('0.25', 24, '0.03336'), ('0.5', 24, '9.447e-7'), ('0.5', 32, '9.40e-7')):
        lam = ym.least_eigenvalue(a, y=0, basis='legendre', size=M, dps=30)
        dec['a=%s M=%d y=0' % (a, M)] = dict(value=mpmath.nstr(lam, 8), deconv_chapter=printed)
        log('deconvolution chapter y=0 Legendre a=%s M=%d: %s (printed %s)' % (a, M, mpmath.nstr(lam, 6), printed))
    ent['deconv_legendre_y0'] = dec
    stamp(res, 'cross', time.time() - t0, 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
def part_figure(res, log, dst):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.6, 2.5))
    g = res.get('figure_grid')
    styles = {0.55: ('-', 'o', '0.0'), 0.6: ('--', 's', '0.25'), 0.75: ('-.', '^', '0.4'), 1.0: (':', 'D', '0.5')}
    if g:
        a_s = sorted(float(k) for k in g)
        for al, (ls, mk, col) in styles.items():
            b = float(ym.b_alpha(al, 20))
            y = [float(g[str(a)][str(al)]) / (2 * b) - 1 for a in a_s]
            ax1.plot(a_s, y, ls=ls, marker=mk, ms=3, color=col, label=r'$\alpha=%g$' % al)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xticks([0.2, 0.5, 1, 2, 5])
        ax1.set_xticklabels(['0.2', '0.5', '1', '2', '5'])
        ax1.minorticks_off()
        ax1.set_xlabel(r'half length $a$')
        ax1.set_ylabel(r'$\lambda_{\min}/(2\xi^\prime(\alpha)/\xi(\alpha))-1$')
        ax1.legend(frameon=False, fontsize=6.5)
    half = res.get('half')
    if half:
        for key, ls, mk, fill, lab in (('big', '-', 's', '0.0', r'$N=\max(24,24a)$'), ('small', '--', 'o', 'white', r'$N=\max(12,12a)$')):
            pts = []
            for a in sorted({v['a'] for v in half.values()}):
                n = max(24, int(math.ceil(24 * a))) if key == 'big' else max(12, int(math.ceil(12 * a)))
                v = half.get('%g_%d' % (a, n))
                if v:
                    pts.append((a, float(mpmath.log10(mpmath.mpf(v['lam'])))))
            if pts:
                ax2.plot([p[0] for p in pts], [p[1] for p in pts], ls=ls, marker=mk, ms=3.5, color='black', markerfacecolor=fill, label=lab)
        dec = res.get('half_decay', {})
        cs = sorted(float(dec[k]['c']) for k in dec if float(k) <= 0.8 and float(k) >= 0.35)
        if cs:
            c = cs[len(cs) // 2]
            aa = np.linspace(0.3, max(float(k) for k in half and [v['a'] for v in half.values()]), 200)
            ax2.plot(aa, -c * 2 * aa * np.exp(2 * aa) / math.log(10), color='0.55', lw=0.8, ls=(0, (1, 1.5)), label=r'$\exp(-%.2g\,ae^{2a})$' % (2 * c))
        ax2.axhline(-16, color='0.5', lw=0.6)
        ax2.text(1.58, -14, 'double precision', ha='right', va='bottom', fontsize=6.5, color='0.35')
        ax2.set_xlim(0.15, 1.6)
        ax2.set_ylim(-70, 3)
        ax2.set_xlabel(r'half length $a$')
        ax2.set_ylabel(r'$\log_{10}\lambda_N$ at $\alpha=1/2$')
        ax2.legend(frameon=False, fontsize=6.5, loc='lower left')
    fig.tight_layout()
    fig.savefig(dst)
    plt.close(fig)
    log('wrote ' + dst)


# ---------------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Table yo:tab:lam and neighbours (see the header)')
    ap.add_argument('--part', default='all', help='ladder, half, floors, line243, cross, figure or all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--fig', default=None)
    ap.add_argument('--a', default=None, help='half: comma separated list of a')
    ap.add_argument('--maxN', type=int, default=None)
    ap.add_argument('--minN', type=int, default=None)
    ap.add_argument('--grid', action='store_true', help='add the dense a grid for the figure to the ladder part')
    args = ap.parse_args()
    out = args.out or default_path('yo_damped.json')
    log = Log(os.path.splitext(out)[0] + '.log')
    log('=== yo_damped.py --part %s%s on %s at %s' % (args.part, ' --quick' if args.quick else '', socket.gethostname(), time.strftime('%Y-%m-%d %H:%M')))
    res = load(out)
    gam = find_zeros()
    parts = ['ladder', 'legendre', 'half', 'floors', 'line243', 'cross', 'figure'] if args.part == 'all' else [args.part]
    for part in parts:
        if part == 'ladder':
            part_ladder(res, log, grid=args.grid or args.part == 'all')
        elif part == 'legendre':
            part_legendre(res, log)
        elif part == 'half':
            part_half(res, log, args.quick, args.procs, out, args.a, args.maxN, args.minN)
        elif part == 'floors':
            part_floors(res, log, gam)
        elif part == 'line243':
            part_line243(res, log, gam, args.quick)
        elif part == 'cross':
            part_cross(res, log)
        elif part == 'figure':
            here = os.path.dirname(os.path.abspath(__file__))
            dst = args.fig or (os.path.join(os.path.dirname(here), 'fig', 'yo_lam.pdf') if os.path.basename(here) == 'code' else 'yo_lam.pdf')
            part_figure(res, log, dst)
        save(out, res)
