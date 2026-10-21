# gl_zeta_moments.py: the numerical claims about the zeta function in Chapter 21 (Thorin Representations, Dirichlet L Functions and the
# Linnik Constant), Section gl:sec:numerics and the statements of Section gl:sec:critical that quote it: the moment table of equation
# gl:eq:moments, the anchored tails of equation gl:eq:anchored, the classical constant and the tilt by one quarter (Proposition gl:prop:tilt,
# equation gl:eq:classical-series), the tilt u0 = 1, the transfer function T(t) (Proposition gl:prop:transfer), the negative moments of
# Proposition gl:prop:negmom and the interpolation of Proposition gl:prop:interp. Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. H is the clock of xi, E exp(-u H) = xi(1/2)/xi(1/2 + sqrt u) = 1/G(u) with G(u) = sum_j [s^(2j)] xi(1/2 + s) u^j / xi(1/2);
# p_j = sum over the zeros gamma > 0 of gamma^(-2j); E H^k = k! h_k(p) with h_k the complete homogeneous symmetric polynomials (Newton).
# Function side. The Taylor coefficients of xi(1/2 + s) are the trapezoid sums on circles |s| = 3 and |s| = 2 (two radii, two node counts,
# their agreement is logged and is the error control); p_j and E H^k follow by series arithmetic. An independent determination of the
# coefficients of 1/G(u) is the trapezoid rule on |u| = 50 (inside the radius gamma_1^2 = 199.79 of convergence) applied to 1/G directly.
# Zero side. data/zeros_1700.txt (68 digits, the book's list). Tails: the Riemann von Mangoldt density (1/2 pi) log(t/2 pi) as the chapter
# states it, the smooth density theta'(t)/pi (closed form from the asymptotic series of theta', checked once against mp.quad of the
# digamma function) and the anchored tail of equation gl:eq:anchored, smooth tail minus S(T) T^(-2j) with S(T) = N(T) - theta(T)/pi - 1.
# Tails of sums of other functions of gamma^2 (the tilts, T(t)) use the exact series in gamma^(-2), so no quadrature enters the zero side.
# Negative moments. The two forms of equation gl:eq:negmom are the same integral after u = w^2, so their disagreement measures the
# quadrature at the singular endpoint; the trusted value is the substitution u = v^4 (smooth integrand). The entire continuation is
# M(s) = sum_k d_k/(s + k) + int_1^infty u^(s - 1)/G(u) du with d_k the Taylor coefficients of 1/G, evaluated at s = -n + delta (residue
# ratio) and on Re s = -2 (the bound |E H^(-s)| <= E H^(-Re s) holds because H is a positive variable with all moments finite).
# Output: data/gl_zeta_moments.json (numbers as strings with the digits that are trusted, parameters, counts) and the log gl_zeta_moments.log
# next to it. Nothing here proves anything: finite precision arithmetic on a list of 1700 ordinates.
# Usage: python3 gl_zeta_moments.py [--quick] [--out FILE] [--zeros FILE] [--dps D] [--line-dps D] [--taus LIST]; --quick runs the whole claim
# list at 50 digits in about 10 seconds in one process and writes to /tmp unless --out is given; the full case (100 digits, every cut N from 50
# to 1699, heights up to 80 on Re s = -2 at 200 digits, which costs 1.4 tau + 40 digits) is meant for the Hopper cluster (one core, up to an
# hour): bash misc/tools/hopper_run.sh -c 1 -t 180 -g gl_zeta_moments.json -g gl_zeta_moments.log code/gl_zeta_moments.py data/zeros_1700.txt
# -- --out gl_zeta_moments.json, then move the two files to data/.
import argparse
import json
import math
import statistics
import sys
import tempfile
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
T_START = time.time()
LOGFILE = None
HALF = mp.mpf(1) / 2
XIH = None   # xi(1/2), set in main
BK = None    # coefficients of the asymptotic series of theta'


def log(*a):
    line = '[%6.1fs] %s' % (time.time() - T_START, ' '.join(str(x) for x in a))
    print(line, flush=True)
    if LOGFILE is not None:
        LOGFILE.write(line + '\n')
        LOGFILE.flush()


def fs(x, d=25):
    return mp.nstr(x, d)


def rel(a, b):
    return abs(a - b) / abs(b)


def xi(s):
    """xi(s) = s (s - 1)/2 pi^(-s/2) Gamma(s/2) zeta(s), xi(0) = xi(1) = 1/2."""
    s = mp.mpmathify(s)
    if s == 1 or s == 0:
        return HALF
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def Fval(t):
    """xi(1/2)/xi(1/2 + t) for real t >= 0, the Laplace transform of H at u = t^2."""
    return XIH / xi(HALF + t)


def theta(t):
    return mp.im(mp.loggamma(mp.mpc(mp.mpf(1) / 4, t / 2))) - t / 2 * mp.log(mp.pi)


def dtheta(t):
    return mp.re(mp.digamma(mp.mpc(mp.mpf(1) / 4, t / 2))) / 2 - mp.log(mp.pi) / 2


def bcoef(K):
    """theta'(t) = (1/2) log(t/2 pi) - sum_k b_k t^(-2k), b_k = (1 - 2^(1-2k)) |B_2k| / (4k)."""
    return [None] + [(1 - mp.mpf(2) ** (1 - 2 * k)) * abs(mp.bernoulli(2 * k)) / (4 * k) for k in range(1, K + 1)]


def tail_rvm(T, j):
    """int_T^infty t^(-2j) (1/2 pi) log(t/2 pi) dt."""
    a = 2 * j - 1
    return T ** (-a) / a * (mp.log(T / (2 * mp.pi)) + mp.mpf(1) / a) / (2 * mp.pi)


def tail_smooth(T, j):
    """int_T^infty t^(-2j) theta'(t)/pi dt from the asymptotic series (T >= 20 here)."""
    a = 2 * j - 1
    s = T ** (-a) / a * (mp.log(T / (2 * mp.pi)) + mp.mpf(1) / a) / 2
    for k in range(1, len(BK)):
        s -= BK[k] * T ** (-a - 2 * k) / (a + 2 * k)
    return s / mp.pi


def tail_smooth_quad(T, j):
    with mp.workdps(30):
        return mp.quad(lambda t: t ** (-2 * j) * dtheta(t) / mp.pi, [T, 2 * T, 8 * T, 64 * T, mp.inf])


def tailp(kind, T, j, N):
    """Tail of sum gamma^(-2j) over gamma > T, N = number of ordinates below T."""
    if kind == 'none':
        return mp.mpf(0)
    if kind == 'rvm':
        return tail_rvm(T, j)
    s = tail_smooth(T, j)
    if kind == 'smooth':
        return s
    return s - (N - theta(T) / mp.pi - 1) * T ** (-2 * j)   # anchored


def tail_of_function(coefs, kind, T, N):
    """Tail of sum f(gamma) for f = sum_j coefs[j] gamma^(-2j)."""
    return sum(coefs[j] * tailp(kind, T, j, N) for j in range(1, len(coefs)))


class Zeros:
    def __init__(self, gam, J):
        self.gam = gam
        self.cum = [None]
        pw = [mp.mpf(1)] * len(gam)
        inv2 = [1 / (g * g) for g in gam]
        for j in range(1, J + 1):
            acc = mp.mpf(0)
            c = [mp.mpf(0)]
            for i in range(len(gam)):
                pw[i] = pw[i] * inv2[i]
                acc += pw[i]
                c.append(acc)
            self.cum.append(c)

    def mid(self, N):
        return (self.gam[N - 1] + self.gam[N]) / 2

    def p(self, j, N, kind, T=None):
        return self.cum[j][N] + tailp(kind, T, j, N)


def h_from_p(p, K):
    h = [mp.mpf(1)]
    for k in range(1, K + 1):
        h.append(sum(p[i] * h[k - i] for i in range(1, k + 1)) / k)
    return h


def moments_from_p(p, K):
    h = h_from_p(p, K)
    return [mp.factorial(k) * h[k] for k in range(K + 1)]


def circle_coeffs(vals, M, r, nmax):
    """Taylor coefficients 0..nmax from the values of a real analytic function on the upper half of |z| = r (nodes 2 pi m/M, m = 0..M/2)."""
    half = M // 2
    out = []
    for n in range(nmax + 1):
        acc = mp.re(vals[0]) + (-1) ** n * mp.re(vals[half])
        for m in range(1, half):
            acc += 2 * mp.re(vals[m] * mp.expjpi(-2 * mp.mpf(m * n) / M))
        out.append(acc / (M * r ** n))
    return out


def xi_taylor(J, r, M):
    """c_n = [s^n] xi(1/2 + s), n = 0..2J+1 (the odd ones vanish: that is the evenness check)."""
    r = mp.mpf(r)
    vals = [xi(HALF + r * mp.expjpi(2 * mp.mpf(m) / M)) for m in range(M // 2 + 1)]
    return circle_coeffs(vals, M, r, 2 * J + 1)


def series_inv(g, K):
    d = [mp.mpf(1)] + [mp.mpf(0)] * K
    for k in range(1, K + 1):
        d[k] = -sum(g[j] * d[k - j] for j in range(1, min(k, len(g) - 1) + 1))
    return d


def series_log(g, K):
    l = [mp.mpf(0)] * (K + 1)
    for k in range(1, K + 1):
        l[k] = g[k] - sum(mp.mpf(j) * l[j] * g[k - j] for j in range(1, k)) / k
    return l


def digits(x):
    return int(-mp.log10(x)) if x > 0 else 999


def main():
    global LOGFILE, XIH, BK
    ap = argparse.ArgumentParser(description='zeta function claims of Chapter 21 (moments, tails, tilts, transfer function, negative moments)')
    ap.add_argument('--quick', action='store_true', help='50 digits, thinned sweeps, writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/gl_zeta_moments.json; quick: a /tmp file)')
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    ap.add_argument('--dps', type=int, default=None)
    ap.add_argument('--line-dps', type=int, default=None, help='working precision of the evaluation on Re s = -2 (default: 200, quick: the working precision)')
    ap.add_argument('--taus', default=None, help='comma separated heights tau on Re s = -2 (default 0.5,...,80; quick 1,5,20)')
    a = ap.parse_args()
    quick = a.quick
    DPS = a.dps or (50 if quick else 100)
    out = Path(a.out) if a.out else (Path(tempfile.gettempdir()) / 'gl_zeta_moments_quick.json' if quick else HERE.parent / 'data' / 'gl_zeta_moments.json')
    LOGFILE = open(out.with_suffix('.log'), 'w')
    mp.mp.dps = DPS
    BK = bcoef(14)
    cands = [Path(a.zeros)] if a.zeros else []
    cands += [HERE.parent / 'data' / 'zeros_1700.txt', Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), HERE / 'zeros_1700.txt']
    zf = next((c for c in cands if c.exists()), None)
    if zf is None:
        sys.exit('zeros_1700.txt not found, pass --zeros')
    gam = [mp.mpf(x) for x in zf.read_text().split()]
    log('zeros', len(gam), 'from', zf, 'gamma_1', fs(gam[0], 20), 'gamma_800', fs(gam[799], 20), 'dps', DPS)
    R = dict(script='gl_zeta_moments.py', mode='quick' if quick else 'full')
    XIH = xi(HALF)
    R['xi_half'] = fs(XIH, 40)
    log('xi(1/2) =', fs(XIH, 30))
    JM = 10 if quick else 16
    Z = Zeros(gam, 8)

    # ---- function side
    M1, M2 = (48, 64) if quick else (96, 128)
    c1 = xi_taylor(JM, 3, M1)
    c2 = xi_taylor(JM, 2, M2)
    dev = max(rel(c1[2 * j], c2[2 * j]) for j in range(JM + 1))
    dev4 = max(rel(c1[2 * j], c2[2 * j]) for j in range(5))
    odd = max(abs(c1[2 * j + 1]) / abs(c1[0]) for j in range(JM + 1))
    log('Taylor coefficients of xi(1/2 + s): radii 3 and 2 agree to relative', fs(dev4, 3), 'for the coefficients of s^0..s^8 (the ones the claims use) and', fs(dev, 3),
        'up to s^%d (the high orders are limited by the working precision); odd coefficients below' % (2 * JM), fs(odd, 3), '; c_0 vs xi(1/2)', fs(rel(c1[0], XIH), 3))
    g = [c1[2 * j] / c1[0] for j in range(JM + 1)]
    lg = series_log(g, JM)
    pref = [None] + [(-1) ** (k + 1) * k * lg[k] for k in range(1, JM + 1)]
    dinv = series_inv(g, JM)
    EHref = [(-1) ** k * mp.factorial(k) * dinv[k] for k in range(JM + 1)]
    # independent determination: trapezoid on |u| = 50 for 1/G(u) = xi(1/2)/xi(1/2 + sqrt u)
    MU, NU = (96, 40) if quick else (256, 60)
    rho = mp.mpf(50)
    vu = [XIH / xi(HALF + mp.sqrt(rho * mp.expjpi(2 * mp.mpf(m) / MU))) for m in range(MU // 2 + 1)]
    dU = circle_coeffs(vu, MU, rho, NU)
    devU = max(rel(dU[k], dinv[k]) for k in range(JM + 1))
    devU4 = max(rel(dU[k], dinv[k]) for k in range(5))
    EHu = [(-1) ** k * mp.factorial(k) * dU[k] for k in range(NU + 1)]
    log('p_1 =', fs(pref[1], 30), ' p_2 =', fs(pref[2], 25), ' p_3 =', fs(pref[3], 25))
    log('u plane trapezoid (|u| = 50, M = %d) vs series of 1/G: max relative difference of d_0..d_4 = %s, of d_0..d_%d = %s' % (MU, fs(devU4, 3), JM, fs(devU, 3)))
    R['function_side'] = dict(p=[fs(pref[j], 40) for j in range(1, 9)], EH=[fs(EHref[k], 40) for k in range(0, 6)],
                              params=dict(radii=[3, 2], nodes=[M1, M2], u_plane=dict(rho=50, M=MU)), agree_radii_up_to_s8=fs(dev4, 3), agree_radii_all=fs(dev, 3), agree_uplane_d0_d4=fs(devU4, 3), agree_uplane_all=fs(devU, 3), odd_coeff=fs(odd, 3))
    log('E H^k function side k = 1..4:', [fs(EHref[k], 15) for k in range(1, 5)])

    # ---- the moment table (C:293 to 301)
    N = 800
    gN, mid = gam[N - 1], Z.mid(N)
    variants = [('rvm_from_gamma_N', 'rvm', gN), ('none', 'none', None), ('rvm_from_midpoint', 'rvm', mid),
                ('smooth_midpoint', 'smooth', mid), ('anchored_midpoint', 'anchored', mid)]
    tab = {}
    for name, kind, T in variants:
        p = [None] + [Z.p(j, N, kind, T) for j in range(1, 5)]
        eh = moments_from_p(p, 4)
        tab[name] = dict(zero_side=[fs(eh[k], 12) for k in range(1, 5)], rel_error=[fs(rel(eh[k], EHref[k]), 4) for k in range(1, 5)],
                         p=[fs(p[j], 15) for j in range(1, 4)])
        log('table N = 800,', name, ': zero side E H^k =', tab[name]['zero_side'], ' relative error', tab[name]['rel_error'])
    # sensitivity of p_1 to a missing zero and to an off line pair (what the agreement can and cannot detect, finding M11a)
    sens = {}
    for kk in (800, 400, 100):
        sens['relative_change_of_p1_if_gamma_%d_is_missing' % kk] = fs(1 / gam[kk - 1] ** 2 / pref[1], 4)
    dl = mp.mpf(1) / 4
    aoff = -(dl + mp.mpc(0, 1) * gam[799]) ** 2
    sens['relative_change_of_p1_if_a_pair_at_gamma_800_moves_off_the_line_by_0.25'] = fs((2 * mp.re(1 / aoff) - 2 / gam[799] ** 2) / pref[1], 4)
    log('sensitivity of p_1:', sens)
    R['table_N800'] = dict(T_gamma_N=fs(gN, 15), T_midpoint=fs(mid, 15), variants=tab, function_side=[fs(EHref[k], 12) for k in range(1, 5)], sensitivity=sens)

    # ---- anchored tails (C:303 to 310): N = 799, 200, 400 and a sweep
    check = {}
    for N, js in ((200, (1, 2, 3)), (799, (1, 2, 3))) if not quick else ((200, (1,)),):
        T = Z.mid(N)
        for j in js:
            q = tail_smooth_quad(T, j)
            check['N%d_j%d' % (N, j)] = fs(rel(tail_smooth(T, j), q), 3)
    log('smooth tail (asymptotic series) vs mp.quad of the digamma density, relative difference:', check)
    anc = {}
    for N in (799, 200, 400):
        T = Z.mid(N)
        row = dict(T=fs(T, 12), S_of_T=fs(N - theta(T) / mp.pi - 1, 6), smooth=[], anchored=[], gain=[])
        for j in (1, 2, 3):
            ref = pref[j]
            es = rel(Z.p(j, N, 'smooth', T), ref)
            ea = rel(Z.p(j, N, 'anchored', T), ref)
            row['smooth'].append(fs(es, 4))
            row['anchored'].append(fs(ea, 4))
            row['gain'].append(fs(es / ea, 4))
        anc['N%d' % N] = row
        log('anchored tails N = %d, T = %s, S(T) = %s: relative error of p_1,p_2,p_3 smooth %s anchored %s gain %s' % (N, row['T'], row['S_of_T'], row['smooth'], row['anchored'], row['gain']))
    R['anchored'] = anc
    step = 10 if quick else 1
    sweep = []
    for N in range(50, 1700, step):
        T = Z.mid(N)
        es = [rel(Z.p(j, N, 'smooth', T), pref[j]) for j in (1, 2, 3)]
        ea = [rel(Z.p(j, N, 'anchored', T), pref[j]) for j in (1, 2, 3)]
        sweep.append((N, T, es, ea, abs(N - theta(T) / mp.pi - 1)))
    summ = {}
    for ji, j in enumerate((1, 2, 3)):
        gains = [float(es[ji] / ea[ji]) for (_, _, es, ea, _) in sweep]
        gt = [float(es[ji] / ea[ji] / T) for (_, T, es, ea, _) in sweep]
        # the smooth error is |S(T)| T^(-2j)/p_j to leading order (the second term of equation gl:eq:anchored)
        pr = [float(es[ji] / (aS * T ** (-2 * j) / pref[j])) for (_, T, es, ea, aS) in sweep]
        dec = statistics.quantiles(gains, n=10)
        summ['p%d' % j] = dict(median_gain=fs(statistics.median(gains), 4), decile_1=fs(dec[0], 4), decile_9=fs(dec[-1], 4), min_gain=fs(min(gains), 4), max_gain=fs(max(gains), 4),
                               fraction_gain_below_10=fs(sum(1 for x in gains if x < 10) / len(gains), 3), median_gain_over_T=fs(statistics.median(gt), 4),
                               median_smooth_error_over_S_prediction=fs(statistics.median(pr), 5))
    summ['median_abs_S_at_cut'] = fs(statistics.median([float(x[4]) for x in sweep]), 4)
    R['anchored_sweep'] = dict(N_values=len(sweep), step=step, N_first=sweep[0][0], N_last=sweep[-1][0], summary=summ)
    log('anchored versus smooth over', len(sweep), 'cuts N =', sweep[0][0], 'to', sweep[-1][0], '(gain = smooth error / anchored error):')
    for key, val in summ.items():
        log('   ', key, val)

    # ---- classical constant and the tilt by 1/4 (C:318 to 341)
    cl = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
    ser = sum((-1) ** (j + 1) * mp.mpf(4) ** (1 - j) * pref[j] for j in range(1, JM + 1))
    cs = [None] + [(-1) ** (j + 1) * mp.mpf(4) ** (1 - j) for j in range(1, 9)]
    cp = [None] + [(-1) ** j * mp.mpf(4) ** (-j) / j for j in range(1, 9)]
    N = 800
    head_s = sum(1 / (x * x + mp.mpf(1) / 4) for x in gam[:N])
    head_p = -sum(mp.log(1 + 1 / (4 * x * x)) for x in gam[:N])
    cls = dict(constant=fs(cl, 30), two_xi_half=fs(2 * XIH, 30), p1=fs(pref[1], 30), p1_minus_constant=fs(pref[1] - cl, 8),
               series_value=fs(ser, 35), series_relative_error=fs(rel(ser, cl), 3), series_terms=JM, sum_N800={}, prod_N800={})
    for name, kind, T in (('rvm_from_gamma_N', 'rvm', gam[N - 1]), ('anchored_midpoint', 'anchored', Z.mid(N)), ('none', 'none', None)):
        cls['sum_N800'][name] = fs(head_s + tail_of_function(cs, kind, T, N), 15)
        cls['prod_N800'][name] = fs(mp.exp(head_p + tail_of_function(cp, kind, T, N)), 15)
    R['classical'] = cls
    log('classical constant', fs(cl, 20), '; 2 xi(1/2) =', fs(2 * XIH, 20), '; E H - constant =', fs(pref[1] - cl, 6), '; series of the p_j (%d terms) agrees to relative %s' % (JM, fs(rel(ser, cl), 3)))
    log('sum 1/(gamma^2 + 1/4), 800 zeros, by tail:', cls['sum_N800'], '; product gamma^2/(gamma^2 + 1/4):', cls['prod_N800'])

    # ---- the tilt u0 = 1 (C:199, C:731)
    tilt = dict(one_plus_gamma1_sq=fs(1 + gam[0] ** 2, 12), cases={})
    N = 400
    for s in (3, 50):
        s = mp.mpf(s)
        head = sum(mp.log((1 + x * x) / (1 + x * x + s)) for x in gam[:N])
        ct = [None] + [(-1) ** (j + 1) * (1 - (1 + s) ** j) / j for j in range(1, 9)]
        fn = xi(mp.mpf(3) / 2) / xi(HALF + mp.sqrt(1 + s))
        d = dict(function_side=fs(fn, 15))
        for name, kind, T in (('rvm_from_gamma_N', 'rvm', gam[N - 1]), ('anchored_midpoint', 'anchored', Z.mid(N)), ('none', 'none', None)):
            v = mp.exp(head + tail_of_function(ct, kind, T, N))
            d[name] = fs(v, 15)
            d[name + '_rel'] = fs(rel(v, fn), 4)
        tilt['cases']['s=%d' % int(s)] = d
        log('tilt u0 = 1, s = %d: function side %s; zero side (400 zeros) %s' % (int(s), d['function_side'], {k: (d[k], d[k + '_rel']) for k in ('rvm_from_gamma_N', 'anchored_midpoint', 'none')}))
    R['tilt_u0_1'] = tilt
    log('1 + gamma_1^2 =', fs(1 + gam[0] ** 2, 10))

    # ---- transfer function (C:377)
    tr = dict(values={}, small_t={}, p1=fs(pref[1], 20), minus_2p1=fs(-2 * pref[1], 20))
    Nn = 1699
    Tn = Z.mid(Nn)
    for t in (2, 5, 13):
        t = mp.mpf(t)
        z = xi(HALF + mp.mpc(0, t)) / xi(HALF + t)
        head = sum(mp.log(abs((x * x - t * t) / (x * x + t * t))) for x in gam[:Nn])
        sign = sum(1 for x in gam[:Nn] if x < t)
        ct = [None] + [(-2 * t ** (2 * j) / j if j % 2 == 1 else mp.mpf(0)) for j in range(1, 9)]
        zs = mp.exp(head + tail_of_function(ct, 'anchored', Tn, Nn))
        tr['values'][str(int(t))] = dict(T=fs(mp.re(z), 15), imag_over_real=fs(abs(mp.im(z) / mp.re(z)), 3), product_1699_zeros_anchored=fs(zs, 15), product_rel=fs(rel(zs, mp.re(z)), 4))
        log('T(%d) = %s (imaginary part relative %s); product over 1699 zeros with anchored tail %s (relative difference %s)' % (t, fs(mp.re(z), 12), fs(abs(mp.im(z) / mp.re(z)), 3), fs(zs, 12), fs(rel(zs, mp.re(z)), 3)))
    for t in ('0.1', '0.01', '0.001'):
        t = mp.mpf(t)
        z = xi(HALF + mp.mpc(0, t)) / xi(HALF + t)
        v = mp.log(mp.re(z)) / t ** 2
        tr['small_t'][str(t)] = fs(v, 15)
        log('log T(t)/t^2 at t =', fs(t, 3), '=', fs(v, 15), ' (limit -2 p_1 =', fs(-2 * pref[1], 15), ')')
    R['transfer'] = tr

    # ---- negative moments (C:405 to 411)
    V = 2
    while Fval(mp.mpf(V) ** 2) > mp.mpf(10) ** (-(DPS + 12)):
        V += 1
    nm = {}
    for s in (mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4):
        reg = 4 * mp.quad(lambda v: v ** (4 * s - 1) * Fval(v * v), [mp.mpf(i) for i in range(V + 1)]) / mp.gamma(s)
        with mp.workdps(30):
            uf = mp.quad(lambda u: u ** (s - 1) * Fval(mp.sqrt(u)), [0, 1, 10, 100, 1000, 10000]) / mp.gamma(s)
            wf = 2 * XIH * mp.quad(lambda w: w ** (2 * s - 1) / xi(HALF + w), [0, 1, 10, 100]) / mp.gamma(s)
        nm[fs(s, 3)] = dict(regularised=fs(reg, 40), u_form_30_digits=fs(uf, 20), w_form_30_digits=fs(wf, 20), u_form_vs_regularised=fs(rel(uf, reg), 3),
                            w_form_vs_regularised=fs(rel(wf, reg), 3), u_form_vs_w_form=fs(rel(uf, wf), 3))
        log('E H^(-s), s = %s: regularised %s; plain u form %s (relative %s); plain w form %s (relative %s); forms differ by %s'
            % (fs(s, 3), fs(reg, 22), fs(uf, 18), fs(rel(uf, reg), 3), fs(wf, 18), fs(rel(wf, reg), 3), fs(rel(uf, wf), 3)))
    R['negative_moments'] = dict(values=nm, upper_limit_v=V)

    # ---- residues at s = -n and the line Re s = -2 (C:428)
    def make_cont(dps_c):
        """M(s) = sum_k d_k/(s + k) + int_0^X e^(s x) F(e^x) dx at working precision dps_c (the coefficients d_k are recomputed there)."""
        global XIH
        mp.mp.dps = dps_c
        XIH = xi(HALF)
        Mu = 96 if quick else 256
        rho_c = mp.mpf(50)
        vu_c = [XIH / xi(HALF + mp.sqrt(rho_c * mp.expjpi(2 * mp.mpf(m) / Mu))) for m in range(Mu // 2 + 1)]
        d_c = circle_coeffs(vu_c, Mu, rho_c, NU)
        X = 10 if dps_c <= 100 else 12
        cache = {}

        def Fx(x):
            v = cache.get(x)
            if v is None:
                v = cache[x] = Fval(mp.exp(x / 2))
            return v

        seg = [mp.mpf(i) / 4 for i in range(4 * X + 1)]

        def contM(s):
            ser = sum(d_c[k] / (s + k) for k in range(len(d_c)))
            integ = mp.quad(lambda x: mp.exp(s * x) * Fx(x), seg, maxdegree=8)
            return ser + integ
        return contM

    contM = make_cont(DPS)
    s0 = mp.mpf(1) / 2
    chk = contM(s0) * mp.rgamma(s0)
    reg_half = R['negative_moments']['values']['0.5']['regularised']
    log('continuation formula at s = 1/2:', fs(chk, 20), 'vs regularised', reg_half[:22], '(relative difference', fs(rel(chk, mp.mpf(reg_half)), 3), ')')
    delta = mp.mpf(10) ** (-(DPS // 2))
    res = {}
    for n in range(5):
        sn = -n + delta
        v = contM(sn) * mp.rgamma(sn)
        # zero side moment with all 1699 ordinates and the anchored tail
        pz = [None] + [Z.p(j, 1699, 'anchored', Z.mid(1699)) for j in range(1, 5)]
        ehz = moments_from_p(pz, 4)
        res['n=%d' % n] = dict(residue_ratio=fs(v, 30), function_side=fs(EHref[n], 30), uplane_side=fs(EHu[n], 30), zero_side_1699_anchored=fs(ehz[n], 15),
                               ratio_vs_function=fs(rel(v, EHref[n]), 3), zero_vs_function=fs(rel(ehz[n], EHref[n]), 3))
        log('s = -%d + delta (delta = %s): E H^(-s) = %s; E H^%d function side %s; relative difference %s; zero side (1699, anchored) relative %s'
            % (n, fs(delta, 2), fs(v, 14), n, fs(EHref[n], 14), fs(rel(v, EHref[n]), 3), fs(rel(ehz[n], EHref[n]), 3)))
    taus = [float(x) for x in a.taus.split(',')] if a.taus else ([1, 5, 20] if quick else [0.5, 1, 2, 5, 10, 20, 40, 80])
    line = {}
    bound = EHref[2]
    # M(-2 + i tau) is of size exp(-pi tau) while its pieces are of size one: the digits needed are 1.4 tau + 40
    dps_line = a.line_dps or (DPS if quick else 200)
    if dps_line != DPS:
        contM = make_cont(dps_line)
        log('line Re s = -2 at working precision', dps_line, 'digits (cancellation of about 1.4 tau digits)')
    for tau in taus:
        s = mp.mpc(-2, tau)
        Mv = contM(s)
        v = Mv * mp.rgamma(s)
        line['tau=%g' % tau] = dict(abs_E_H_neg_s=fs(abs(v), 12), abs_M=fs(abs(Mv), 8), abs_one_over_gamma=fs(abs(mp.rgamma(s)), 8), below_bound=bool(abs(v) <= bound))
        log('Re s = -2, tau = %s: |E H^(-s)| = %s (bound E H^2 = %s); |M(s)| = %s; |1/Gamma(s)| = %s' % (tau, fs(abs(v), 8), fs(bound, 8), fs(abs(Mv), 6), fs(abs(mp.rgamma(s)), 6)))
    mp.mp.dps = DPS
    XIH = xi(HALF)
    tau33 = mp.findroot(lambda t: -mp.re(mp.loggamma(mp.mpc(-2, t))) - 33 * mp.log(10), 45)
    growth = {str(n): fs(EHu[n], 10) for n in (0, 1, 2, 3, 4, 10, 20, 30, 40) if n <= NU}
    ratios = {str(n): fs(EHu[n + 1] / EHu[n] * gam[0] ** 2 / (n + 1), 6) for n in (10, 20, 30, 39) if n + 1 <= NU}
    R['interpolation'] = dict(residues=res, delta=fs(delta, 3), line_re_s_minus_2=line, bound_E_H2=fs(bound, 12), tau_where_1_over_gamma_is_1e33=fs(tau33, 8),
                              moments_for_growth=growth, ratio_EH_np1_over_EH_n_times_gamma1_sq_over_np1=ratios)
    log('|1/Gamma(-2 + i tau)| = 1e33 at tau =', fs(tau33, 6), '; E H^n for n = 0,1,2,3,4,10,20,30,40:', growth)
    log('E H^(n+1)/E H^n * gamma_1^2/(n+1) for n = 10, 20, 30, 39 (tends to 1 if the growth is factorial over gamma_1^(2n)):', ratios)

    R['parameters'] = dict(dps=DPS, zeros=len(gam), zeros_file=str(zf), zero_digits=68, taylor_J=JM, quick=quick, seconds=fs(time.time() - T_START, 4))
    json.dump(R, open(out, 'w'), indent=1)
    log('wrote', out)


if __name__ == '__main__':
    main()
