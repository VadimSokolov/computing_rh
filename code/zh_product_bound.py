# zh_product_bound.py: Zhu's product test function in the Thorin reading, Section zh:sec:extremal of Chapter 18 (Zhu's compact windows, the Landau Widom
# horizon and the Thorin program), group Z1. It checks the product table of lines 97 to 99 (best bound, atoms killed, height of the last atom, 2 pi e^(2L)), line 103
# (4.0e-9 at L = 0.8 against the true 2e-17), adds N(T*) and 2 pi e^(2L+1) to it, and evaluates every bound unconditionally on the geometric side. Added for the book in
# October 2026, written from the mathematics (no script of the chapter's author was read or run).
# The test function is F_0 = B h, B(r) = product over the first K ordinates of (1 minus r^2/gamma_j^2), h(r) = (sin(L r/m)/(L r/m))^m with m even and m >= 2K + 2, so that F_0 is
# O(r^(0 minus 2)), f_0 is continuous and supported in [minus L, L]. Under RH, Q(f_0) = 2 sum over j > K of F_0(gamma_j)^2 and ||f_0||^2 = (1/pi) integral over (0, infinity) of F_0^2, and
# R(K, m) = Q(f_0)/||f_0||^2 is an upper bound for lambda*(L) = inf Q(f)/||f||^2 (zero side, under RH). The unconditional version is Zhu's formula (2) (compared with zhu.txt lines 112 to 127: the same
# normalisation as the chapter's lines 22 to 27): Q(f) = 2 F(i/2)^2 + (1/2 pi) integral over R of |F(t)|^2 Psi_L(t) dt, Psi_L(t) = Re psi(1/4 + i t/2) minus log pi minus P_L(t), P_L the prime comb over log n < 2L.
# Method, zero side. Log arithmetic in float64: log |F_0(gamma_j)| = sum over i <= K of log |1 minus gamma_j^2/gamma_i^2| + m log |sin(a gamma_j)/(a gamma_j)|, a = L/m, summed over the 1700
# ordinates of data/zeros_1700.txt (to T_0 = 2197.2587); the top optimum of every L is recomputed in mpmath at 50 digits. Tail beyond T_0. Upper bound: |F_0(t)|^2 <= C t^(0 minus p) for t >= gamma_K,
# p = 2m minus 4K, C = (product of gamma_j^(0 minus 4)) (m/L)^(2m), and sum over gamma > T_0 of C gamma^(0 minus p) <= integral over t > T_0 of (N^+(t) minus 1700) p C t^(0 minus p minus 1) dt, with N^+ the
# Riemann von Mangoldt bound as the book uses it (Hasanalizade, Shen and Wong 2022: |N(T) minus (T/2 pi) log(T/2 pi e) minus 7/8| <= 2.5167 + 1/(50 T) for e <= T <= 30610046000, and
# <= 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50 T) for all T >= e). Central estimate: the same integral against the smooth density log(t/2 pi)/(2 pi), with the exact polynomial B(t)^2 and the mean
# C(2m, m)/4^m of sin^(2m). Norm: exact, ||f_0||^2 = sum_i c_i (0 minus 1)^i p^(2i)(0), where B(r)^2 = sum_i c_i r^(2i) and p is the density of a sum of 2m independent uniform variables on [minus a, a]
# (an Irwin Hall density; p^(2i)(0) = (2a)^(0 minus 1 minus 2i) W(m, i), W an exact rational number), all terms positive; checked at the optimum by adaptive quadrature (scipy quad on each interval of length pi/a,
# relative 1e-13, the tail beyond R bounded by the envelope to 1e-12). Search: K in [0, 15], m in {2K + 2, 2K + 4, ..., 2K + 40}, L in {0.35, 0.5, 0.8, 1.0, 1.2}; when the optimum sits on the edge of
# the box the box is enlarged (K to 25, m to 2K + 80) and both optima are written. The selected optimum is the minimum of the explicit upper bound R_upper; the minimum of the central estimate is written too.
# Method, geometric side (option --geometric). The autocorrelation g = f_0 * f_0 reflected is a piecewise polynomial with knots at the multiples of 2a: g(x) = sum_i c_i (0 minus 1)^i p^(2i)(x), evaluated exactly
# in truncated powers. Comb term: sum over n < e^(2L) of 2 Lambda(n) n^(0 minus 1/2) g(log n). Pole term: 2 F_0(i/2)^2, F_0(i/2) = product of (1 + 1/(4 gamma_j^2)) times (sinh(a/2)/(a/2))^m. Archimedean term
# (1/2 pi) integral of |F|^2 (Re psi(1/4 + i t/2) minus log pi) dt = minus (gamma_E + log pi) g(0) + 2 integral over (0, 2L) of (g(0) e^(0 minus 2x) minus g(x) e^(0 minus x/2))/(1 minus e^(0 minus 2x)) dx minus g(0) log(1 minus e^(0 minus 4L)),
# from psi(s) = minus gamma_E + integral over (0, infinity) of (e^(0 minus t) minus e^(0 minus s t))/(1 minus e^(0 minus t)) dt, by Gauss Legendre on each of the m pieces.
# Error control. (a) Every geometric value is computed twice, at (dps, ngl) and at (dps + 40, ngl + 10) Gauss Legendre nodes per piece, and only the digits on which the two agree are written. (b) g(0) from the
# piecewise form is compared with the exact norm. (c) The geometric value must lie in the bracket of the zero side: 2 S/||f||^2 <= Q/||f||^2 <= 2 (S + tail_upper)/||f||^2 under RH (S the sum over the 1700
# ordinates), so the geometric formula, its normalisation, the pole, archimedean and comb terms and the zero data are tested together; the central tail estimate is printed beside the observed difference.
# Output: data/zh_product_bound.json (values as strings, the claims with their status), the log on standard output (data/zh_product_bound.log).
# Usage: python3 zh_product_bound.py [--quick] [--out FILE] [--zeros FILE] [--geometric] [--Ls L1,L2,..] [--geom-Ls L1,..] [--kmax K] [--mspan M] [--dps D] [--ngl N] [--no-extend];
# --quick searches L = 0.5 and 0.8 with K <= 8 and m <= 2K + 12 and checks the optimum on the geometric side at 50 digits (a few seconds), writing to /tmp unless --out is given; the full geometric
# case (five values of L, 100 digits, 40 nodes per piece) is meant for the Hopper cluster (bash misc/tools/hopper_run.sh, see the commands at the end of misc/reviews/nick-zhu/checks-Z1.md).
# The options --geom-cases best --dps 100 --ngl 12 give the local small case that wrote the stored data/zh_product_bound.json (16 s, the two runs agree to 15 digits or more); the cluster run with
# the defaults (--geom-cases all, 100 and 140 digits, 40 and 50 nodes per piece) adds the other geometric cases (K = 0 with m = 2, the optimum of the unenlarged box) and more digits.
import argparse
import json
import math
import sys
import tempfile
import time
import warnings
from math import comb, factorial
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import integrate

warnings.simplefilter('ignore', integrate.IntegrationWarning)
T0 = time.time()
X0 = 30610046000.0
LS = ['0.35', '0.5', '0.8', '1.0', '1.2']
# the chapter's table of lines 97 to 99: L -> (best bound, atoms killed, height of the last atom, 2 pi e^(2L))
CHAPTER = {'0.5': ('3.9e-4', 2, '21.0', '17.1'), '0.8': ('4.0e-9', 6, '37.6', '31.1'), '1.0': ('7.0e-14', 9, '48.0', '46.4')}
ZERO = {}
LJ = {}


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def sig(x, nd):
    return '%.*e' % (nd - 1, float(x))


def find_zeros(path):
    here = Path(__file__).resolve().parent
    for c in ([Path(path)] if path else []) + [Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), here / 'zeros_1700.txt',
                                                here.parent / 'data' / 'zeros_1700.txt']:
        if c.is_file():
            return c
    raise SystemExit('zeros_1700.txt not found, use --zeros')


def nplus(t):
    main = t / (2 * math.pi) * math.log(t / (2 * math.pi * math.e)) + 0.875
    err = 2.5167 if t <= X0 else 0.1095 * math.log(t) + 0.2042 * math.log(math.log(t)) + 3.0305
    return main + err + 1.0 / (50 * t)


def prime_powers(N):
    sv = bytearray([1]) * (N + 1)
    sv[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
    out = []
    for p in range(2, N + 1):
        if sv[p]:
            q = p
            while q <= N:
                out.append((q, p))
                q *= p
    out.sort()
    return out


def poly_B2(gam, K):
    """coefficients c_0..c_2K of B(r)^2 as a polynomial in s = r^2, B = product over j <= K of (1 minus s/gamma_j^2)"""
    P = [mp.mpf(1)]
    for j in range(K):
        q = 1 / gam[j] ** 2
        Q = [mp.mpf(0)] * (len(P) + 1)
        for i, v in enumerate(P):
            Q[i] += v
            Q[i + 1] -= v * q
        P = Q
    C = [mp.mpf(0)] * (2 * K + 1)
    for i, u in enumerate(P):
        for j, v in enumerate(P):
            C[i + j] += u * v
    return C


_W = {}


def wnum(m, i):
    """(numerator, denominator) of W(m, i) = (1/(2m minus 1 minus 2i)!) sum over r < m of (0 minus 1)^r C(2m, r) (m minus r)^(2m minus 1 minus 2i)"""
    v = _W.get((m, i))
    if v is None:
        e = 2 * m - 1 - 2 * i
        s = 0
        for r in range(m):
            t = comb(2 * m, r) * (m - r) ** e
            s += -t if r & 1 else t
        v = (s, factorial(e))
        _W[(m, i)] = v
    return v


def norm_exact(c, m, L):
    """||f_0||^2 = (1/pi) integral over (0, infinity) of F_0^2, exactly: sum_i c_i (0 minus 1)^i (2a)^(0 minus 1 minus 2i) W(m, i)"""
    a2 = 2 * L / m
    tot = mp.mpf(0)
    for i, ci in enumerate(c):
        num, den = wnum(m, i)
        tot += ci * (-1) ** i * mp.mpf(num) / den / a2 ** (1 + 2 * i)
    return tot


def log_arrays(gd, Kmax):
    out = [np.zeros(len(gd))]
    cur = np.zeros(len(gd))
    for K in range(1, Kmax + 1):
        with np.errstate(divide='ignore'):
            cur = cur + np.log(np.abs(1.0 - (gd / gd[K - 1]) ** 2))
        out.append(cur.copy())
    return out


def zsum(K, m, Lf, lb):
    gd = ZERO['gd']
    x = (Lf / m) * gd[K:]
    with np.errstate(divide='ignore'):
        lw = m * (np.log(np.abs(np.sin(x))) - np.log(x))
    return float(np.exp(2 * (lb[K:] + lw)).sum())


def zsum_mp(K, m, L, gam, dps=50):
    with mp.workdps(dps):
        a = L / m
        tot = mp.mpf(0)
        for j in range(K, len(gam)):
            g = gam[j]
            B = mp.mpf(1)
            for i in range(K):
                B *= 1 - g * g / gam[i] ** 2
            tot += (B * (mp.sin(a * g) / (a * g)) ** m) ** 2
        return tot


def log_J(p):
    """log of integral over t > T_0 of (N^+(t) minus 1700) p t^(0 minus p minus 1) dt"""
    v = LJ.get(p)
    if v is None:
        T0z, ng = ZERO['T0'], ZERO['n']
        vx = math.log(X0 / T0z)

        def f(u):
            return (nplus(T0z * math.exp(u)) - ng) * p * math.exp(-p * u)
        i1 = integrate.quad(f, 0, vx, epsabs=0, epsrel=1e-11, limit=200)[0]
        i2 = integrate.quad(f, vx, vx + 60, epsabs=0, epsrel=1e-8, limit=200)[0]
        v = -p * math.log(T0z) + math.log(i1 + i2)
        LJ[p] = v
    return v


def tail_central(K, m, Lf):
    gk = ZERO['gd'][:K]
    T0z = ZERO['T0']
    lavg = math.log(comb(2 * m, m)) - 2 * m * math.log(2)
    lm = 2 * m * math.log(m / Lf)

    def f(u):
        t = T0z * math.exp(u)
        lb = float(np.sum(np.log(np.abs(1.0 - (t / gk) ** 2)))) if K else 0.0
        nu = math.log(t / (2 * math.pi)) / (2 * math.pi)
        return math.exp(2 * lb + lm - 2 * m * math.log(t) + math.log(nu) + math.log(t) + lavg)
    return integrate.quad(f, 0, 40, epsabs=0, epsrel=1e-7, limit=100)[0]


def log_C(K, m, Lf):
    gk = ZERO['gd'][:K]
    return (-4 * float(np.sum(np.log(gk))) if K else 0.0) + 2 * m * math.log(m / Lf)


def norm_quad(K, m, Lf, norm_ex):
    """adaptive quadrature of (1/pi) integral of F_0^2 over (0, R), intervals of length pi/a, tail beyond R bounded by the envelope"""
    gk = ZERO['gd'][:K]
    a = Lf / m
    p = 2 * m - 4 * K
    target = 1e-12 * math.pi * norm_ex
    lR = (log_C(K, m, Lf) - math.log(p - 1) - math.log(target)) / (p - 1)
    R = min(max(math.exp(lR), 3 * (gk[-1] if K else 1.0)), 5e5)
    tailbound = math.exp(log_C(K, m, Lf) - (p - 1) * math.log(R) - math.log(p - 1))

    def f(r):
        lb = float(np.sum(np.log(np.abs(1.0 - (r / gk) ** 2)))) if K else 0.0
        s = abs(math.sin(a * r))
        if s == 0.0:
            return 0.0
        return math.exp(2 * (lb + m * (math.log(s) - math.log(a * r))))
    edges = np.arange(0.0, R + math.pi / a, math.pi / a)
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        tot += integrate.quad(f, lo, hi, epsabs=0, epsrel=1e-13, limit=100)[0]
    return tot / math.pi, tailbound / (math.pi * norm_ex), R


def search(Ls, Kmax, mspan, cK, lbs):
    gd = ZERO['gd']
    out = {}
    for Lstr in Ls:
        L = mp.mpf(Lstr)
        Lf = float(L)
        recs = []
        for K in range(Kmax + 1):
            for m in range(2 * K + 2, 2 * K + mspan + 1, 2):
                p = 2 * m - 4 * K
                nrm = float(norm_exact(cK[K], m, L))
                S = zsum(K, m, Lf, lbs[K])
                Tu = math.exp(log_C(K, m, Lf) + log_J(p))
                Tc = tail_central(K, m, Lf)
                recs.append({'K': K, 'm': m, 'S': S, 'norm': nrm, 'Tu': Tu, 'Tc': Tc, 'R_zero': 2 * S / nrm, 'R_central': 2 * (S + Tc) / nrm,
                             'R_upper': 2 * (S + Tu) / nrm})
        out[Lstr] = recs
    return out


def gl_nodes(n):
    xs, ws = [], []
    eps = mp.mpf(10) ** (-(mp.mp.dps - 3))
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(60):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < eps:
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def geom_Q(gam, K, m, Lstr, dps, ngl):
    """Q(f_0) on the geometric side, Zhu's formula (2), with ||f_0||^2 = g(0); returns mp numbers at precision dps"""
    with mp.workdps(dps):
        L = mp.mpf(Lstr)
        c = poly_B2(gam, K)
        a = L / m
        a2 = 2 * a
        g0 = norm_exact(c, m, L)
        nI = len(c)
        coef = [c[i] * (-1) ** i / a2 ** (1 + 2 * i) / mp.mpf(factorial(2 * m - 1 - 2 * i)) for i in range(nI)]
        bn = [(-1 if r & 1 else 1) * comb(2 * m, r) for r in range(2 * m + 1)]
        emin = 2 * m - 1 - 2 * (nI - 1)

        def g(x):
            y = m + x / a2
            kk = min(int(mp.floor(y)), 2 * m)
            T = [mp.mpf(0)] * nI
            for r in range(kk + 1):
                z = y - r
                if z <= 0:
                    continue
                zp = z ** emin
                z2 = z * z
                sc = bn[r]
                for i in range(nI - 1, -1, -1):
                    T[i] += sc * zp
                    zp *= z2
            return mp.fsum(coef[i] * T[i] for i in range(nI))
        gcheck = g(mp.mpf(0))
        gend = g(2 * L)
        X, W = gl_nodes(ngl)
        I = mp.mpf(0)
        for q in range(m):
            lo, hi = a2 * q, a2 * (q + 1)
            half, mid = (hi - lo) / 2, (hi + lo) / 2
            s = mp.mpf(0)
            for xk, wk in zip(X, W):
                x = mid + half * xk
                e2 = mp.exp(-2 * x)
                s += wk * (g0 * e2 - g(x) * mp.exp(-x / 2)) / (1 - e2)
            I += s * half
        FB = mp.mpf(1)
        for j in range(K):
            FB *= 1 + 1 / (4 * gam[j] ** 2)
        pole = 2 * (FB * (mp.sinh(a / 2) / (a / 2)) ** m) ** 2
        arch = -(mp.euler + mp.log(mp.pi)) * g0 + 2 * I - g0 * mp.log(1 - mp.exp(-4 * L))
        comb_sum = mp.mpf(0)
        for n, p in prime_powers(int(mp.floor(mp.exp(2 * L))) + 1):
            if mp.log(n) < 2 * L:
                comb_sum += 2 * mp.log(p) / mp.sqrt(n) * g(mp.log(n))
        Q = pole + arch - comb_sum
        return {'g0': g0, 'g_at_0': gcheck, 'g_at_2L': gend, 'pole': pole, 'arch': arch, 'comb': comb_sum, 'Q': Q, 'R': Q / g0}


def digits_agree(x, y):
    if x == y:
        return 999.0
    return float(-mp.log10(abs(x - y) / max(abs(x), abs(y))))


def rec_json(r, Lstr, T_star):
    gd = ZERO['gd']
    K = r['K']
    return {'K': K, 'm': r['m'], 'R_zero_sum_1700': sig(r['R_zero'], 8), 'R_central': sig(r['R_central'], 8), 'R_upper': sig(r['R_upper'], 8),
            'norm_exact': sig(r['norm'], 12), 'tail_upper_over_S': sig(r['Tu'] / r['S'], 3), 'tail_central_over_S': sig(r['Tc'] / r['S'], 3),
            'gamma_K': ('%.10g' % gd[K - 1]) if K else None, 'gamma_K_minus_T_star': ('%.6g' % (gd[K - 1] - T_star)) if K else None}


def main():
    ap = argparse.ArgumentParser(description="Zhu's product test function in the Thorin reading (Chapter 18, group Z1)")
    ap.add_argument('--quick', action='store_true', help='L = 0.5 and 0.8, K <= 8, m <= 2K + 12, geometric check at 50 digits, writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/zh_product_bound.json, or beside the script when there is no data folder)')
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    ap.add_argument('--geometric', action='store_true', help='evaluate Q(f_0) on the geometric side at the optima (always on with --quick)')
    ap.add_argument('--Ls', default=None, help='comma separated list of L (default 0.35,0.5,0.8,1.0,1.2)')
    ap.add_argument('--geom-Ls', default=None, help='the L for which the geometric side is evaluated (default: all of --Ls)')
    ap.add_argument('--kmax', type=int, default=None)
    ap.add_argument('--mspan', type=int, default=None)
    ap.add_argument('--dps', type=int, default=None)
    ap.add_argument('--ngl', type=int, default=None)
    ap.add_argument('--no-extend', action='store_true', help='do not enlarge the box when the optimum sits on its edge')
    ap.add_argument('--geom-cases', default='all', choices=['all', 'best'], help='all: the optima, K = 0 with m = 2 and the chapter K; best: the optimum only (saves time)')
    args = ap.parse_args()
    quick = args.quick
    Ls = args.Ls.split(',') if args.Ls else (['0.5', '0.8'] if quick else LS)
    Kmax = args.kmax if args.kmax is not None else (8 if quick else 15)
    mspan = args.mspan if args.mspan is not None else (12 if quick else 40)
    geometric = args.geometric or quick
    gLs = args.geom_Ls.split(',') if args.geom_Ls else Ls
    dps = args.dps or (50 if quick else 100)
    ngl = args.ngl or (20 if quick else 40)
    Kext = Kmax + 10
    mp.mp.dps = 70
    zpath = find_zeros(args.zeros)
    lines = [ln.strip() for ln in open(zpath) if ln.strip()]
    gam = [mp.mpf(s) for s in lines]
    ZERO['gd'] = np.array([float(s) for s in lines])
    ZERO['n'] = len(lines)
    ZERO['T0'] = float(ZERO['gd'][-1])
    gd = ZERO['gd']
    log('zeros', ZERO['n'], 'to height %.6f' % ZERO['T0'], 'file', zpath, 'Ls', Ls, 'box K <=', Kmax, 'm <= 2K +', mspan, 'geometric', geometric, 'dps', dps, 'ngl', ngl)
    res = {'script': 'zh_product_bound.py', 'mode': 'quick' if quick else 'full',
           'parameters': {'argv': sys.argv[1:], 'geometric_cases': args.geom_cases, 'Ls': Ls, 'K_max': Kmax, 'm_rule': 'even m from 2K + 2 to 2K + %d' % mspan, 'zeros': ZERO['n'], 'T0': lines[-1][:16], 'dps_search': 60,
                          'geometric': geometric, 'geometric_Ls': gLs if geometric else [], 'dps_geometric': dps, 'ngl': ngl, 'ngl_second_run': ngl + 10, 'dps_second_run': dps + 40,
                          'tail': 'N(T) bound of Hasanalizade, Shen and Wong 2022: 2.5167 + 1/(50T) for e <= T <= 30610046000, 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50T) beyond'},
           'results': []}
    res['parameters']['geometric_Ls_not_run'] = [x for x in Ls if not (geometric and x in gLs)]
    mp.mp.dps = 60
    cK = [poly_B2(gam, K) for K in range(Kext + 1)]
    lbs = log_arrays(gd, Kext)
    allrecs = search(Ls, Kmax, mspan, cK, lbs)
    claims = []
    summary = []
    for Lstr in Ls:
        L = mp.mpf(Lstr)
        Lf = float(L)
        T_star = float(2 * mp.pi * mp.exp(2 * L))
        T_lw = float(mp.e * 2 * mp.pi * mp.exp(2 * L))
        Nst = int(np.searchsorted(gd, T_star, side='right'))
        recs = allrecs[Lstr]
        bu = min(recs, key=lambda r: r['R_upper'])
        bc = min(recs, key=lambda r: r['R_central'])
        byK = []
        for K in range(Kmax + 1):
            rk = min((r for r in recs if r['K'] == K), key=lambda r: r['R_upper'])
            byK.append(rk)
        edge = {'K_on_edge': bu['K'] == Kmax, 'm_on_edge': bu['m'] == 2 * bu['K'] + mspan, 'central_K_on_edge': bc['K'] == Kmax, 'central_m_on_edge': bc['m'] == 2 * bc['K'] + mspan}
        log('L', Lstr, 'T* = %.4f, 2 pi e^(2L+1) = %.4f, N(T*) = %d, dimension L T*/pi = %.4f' % (T_star, T_lw, Nst, Lf * T_star / math.pi))
        for r in byK:
            log('   K %2d best m %3d  R_upper %.4e  R_central %.4e  R_zero(1700) %.4e  tail_upper/S %.2e' % (r['K'], r['m'], r['R_upper'], r['R_central'], r['R_zero'], r['Tu'] / r['S']))
        log('   best by R_upper: K %d m %d R_upper %.6e R_central %.6e; best by R_central: K %d m %d R_central %.6e; edge %s' % (bu['K'], bu['m'], bu['R_upper'], bu['R_central'], bc['K'], bc['m'], bc['R_central'], edge))
        entry = {'L': Lstr, 'T_star_2pi_e_2L': '%.8g' % T_star, 'T_LW_2pi_e_2Lplus1': '%.8g' % T_lw, 'N_T_star': Nst, 'dimension_L_T_star_over_pi': '%.8g' % (Lf * T_star / math.pi),
                 'best_upper': rec_json(bu, Lstr, T_star), 'best_central': rec_json(bc, Lstr, T_star), 'edge': edge,
                 'by_K': [{'K': r['K'], 'm': r['m'], 'R_upper': sig(r['R_upper'], 6), 'R_central': sig(r['R_central'], 6), 'R_zero_sum_1700': sig(r['R_zero'], 6)} for r in byK]}
        entry['R_upper_at_K_equal_N_T_star'] = sig(byK[Nst]['R_upper'], 6) if Nst <= Kmax else None
        entry['K_minus_N_T_star_best_upper'] = bu['K'] - Nst
        big, eu = [], None
        if (edge['K_on_edge'] or edge['m_on_edge'] or edge['central_K_on_edge'] or edge['central_m_on_edge']) and not args.no_extend:
            big = search([Lstr], Kext, mspan + 40, cK, lbs)[Lstr]
            eu = min(big, key=lambda r: r['R_upper'])
            ec = min(big, key=lambda r: r['R_central'])
            entry['extended_box'] = {'K_max': Kext, 'm_rule': 'even m from 2K + 2 to 2K + %d' % (mspan + 40), 'best_upper': rec_json(eu, Lstr, T_star), 'best_central': rec_json(ec, Lstr, T_star)}
            log('   box enlarged to K <= %d, m <= 2K + %d: best by R_upper K %d m %d R_upper %.6e; best by R_central K %d m %d R_central %.6e' % (Kext, mspan + 40, eu['K'], eu['m'], eu['R_upper'], ec['K'], ec['m'], ec['R_central']))
        # error control of the optimum: exact norm against adaptive quadrature, float64 sum against mpmath
        chk = {}
        for tag, r in (('upper', bu), ('central', bc)):
            if tag == 'central' and (bc['K'], bc['m']) == (bu['K'], bu['m']):
                continue
            nq, tb, R = norm_quad(r['K'], r['m'], Lf, r['norm'])
            Smp = zsum_mp(r['K'], r['m'], L, gam)
            chk[tag] = {'K': r['K'], 'm': r['m'], 'norm_exact': sig(r['norm'], 14), 'norm_quad': sig(nq, 14), 'norm_relative_difference': '%.2e' % abs(nq / r['norm'] - 1),
                        'quad_range_R': '%.6g' % R, 'quad_tail_bound_relative': '%.2e' % tb, 'S_float64': sig(r['S'], 14), 'S_mpmath_50_digits': mp.nstr(Smp, 20),
                        'S_relative_difference': '%.2e' % float(abs(Smp / mp.mpf(r['S']) - 1))}
            log('   check %s: norm exact %.12e quad %.12e (diff %.1e, R %.4g, envelope tail %.1e); S float64 %.12e mpmath %s (diff %.1e)' % (
                tag, r['norm'], nq, abs(nq / r['norm'] - 1), R, tb, r['S'], mp.nstr(Smp, 14), float(abs(Smp / mp.mpf(r['S']) - 1))))
        entry['error_control'] = chk
        # geometric side
        if geometric and Lstr in gLs:
            cases = [(bu['K'], bu['m'], 'best by R_upper'), (bc['K'], bc['m'], 'best by R_central'), (0, 2, 'K = 0, m = 2 (no atom killed)')]
            if Lstr in CHAPTER:
                kc = CHAPTER[Lstr][1]
                if kc <= Kmax:
                    rkc = byK[kc]
                    cases.append((kc, rkc['m'], 'the chapter K = %d, its best m' % kc))
            if eu is not None:
                cases.insert(0, (eu['K'], eu['m'], 'best by R_upper in the enlarged box'))
            if args.geom_cases == 'best':
                cases = cases[:1]
            seen, gl = set(), []
            for K, m, why in cases:
                if (K, m) in seen:
                    continue
                seen.add((K, m))
                r = next(x for x in recs + big if x['K'] == K and x['m'] == m)
                t1 = time.time()
                A1 = geom_Q(gam, K, m, Lstr, dps, ngl)
                A2 = geom_Q(gam, K, m, Lstr, dps + 40, ngl + 10)
                with mp.workdps(dps + 40):
                    dg = digits_agree(A1['R'], A2['R'])
                    Rg = A2['R']
                    Szero = zsum_mp(K, m, mp.mpf(Lstr), gam, 60)
                    nrm = A2['g0']
                    Rzero = 2 * Szero / nrm
                    lo, hi = Rzero, Rzero + 2 * r['Tu'] / float(nrm)
                    inside = bool(lo <= Rg <= hi)
                    delta = Rg - Rzero
                    ratio_c = float(delta / (2 * r['Tc'] / float(nrm))) if r['Tc'] > 0 else None
                    trusted = max(6, min(int(dg) - 1, 45))
                    gdict = {'K': K, 'm': m, 'why': why, 'R_geometric': mp.nstr(Rg, trusted), 'Q': mp.nstr(A2['Q'], trusted), 'digits_agreeing_between_the_two_runs': '%.1f' % dg,
                             'g0_exact_norm': mp.nstr(nrm, 30), 'g_at_0_minus_norm': mp.nstr(A2['g_at_0'] - nrm, 3), 'g_at_2L': mp.nstr(A2['g_at_2L'], 3),
                             'pole_term': mp.nstr(A2['pole'], 30), 'archimedean_term': mp.nstr(A2['arch'], 30), 'comb_term': mp.nstr(A2['comb'], 30),
                             'R_zero_sum_1700_mpmath': mp.nstr(Rzero, 20), 'R_zero_plus_tail_upper': mp.nstr(hi, 12), 'R_zero_plus_tail_central': sig(Rzero + 2 * r['Tc'] / float(nrm), 12),
                             'geometric_inside_zero_side_bracket': inside, 'geometric_minus_zero_over_central_tail': ('%.4f' % ratio_c) if ratio_c is not None else None,
                             'R_central_float64': sig(r['R_central'], 10), 'R_upper_float64': sig(r['R_upper'], 10), 'seconds': '%.1f' % (time.time() - t1)}
                gl.append(gdict)
                log('   geometric K %d m %d (%s): R = %s (two runs agree to %.1f digits); zero side 2S/norm = %s, plus central tail %s, plus upper tail %s; inside bracket %s; (R minus zero)/central tail %s' % (
                    K, m, why, mp.nstr(Rg, 14), dg, mp.nstr(Rzero, 14), gdict['R_zero_plus_tail_central'], gdict['R_zero_plus_tail_upper'], inside, gdict['geometric_minus_zero_over_central_tail']))
            entry['geometric'] = gl
        res['results'].append(entry)
        # the chapter's claims
        if Lstr in CHAPTER:
            bound, katoms, hlast, hprinted = CHAPTER[Lstr]
            geo_best = next((g for g in entry.get('geometric', []) if g['K'] == bu['K'] and g['m'] == bu['m']), None)
            claims.append({'line': '97 to 99', 'claim': 'best bound at L = %s' % Lstr, 'chapter': bound, 'mine': 'central %s, upper %s%s' % (sig(bc['R_central'], 3), sig(bu['R_upper'], 3),
                           (', geometric %s' % sig(mp.mpf(geo_best['R_geometric']), 3)) if geo_best else ''), 'status': 'agrees' if sig(bc['R_central'], 2) == sig(float(bound), 2) else 'differs'})
            claims.append({'line': '97 to 99', 'claim': 'atoms killed at L = %s' % Lstr, 'chapter': str(katoms), 'mine': 'best by R_central K = %d, by R_upper K = %d' % (bc['K'], bu['K']),
                           'status': 'agrees' if bc['K'] == katoms else 'differs'})
            claims.append({'line': '97 to 99', 'claim': 'height of the last atom killed at L = %s' % Lstr, 'chapter': hlast, 'mine': '%.4f (gamma_%d), %.4f for the R_upper optimum (gamma_%d)' % (
                gd[bc['K'] - 1], bc['K'], gd[bu['K'] - 1], bu['K']), 'status': 'agrees' if ('%.1f' % gd[bc['K'] - 1]) == hlast else 'differs'})
            claims.append({'line': '97 to 99', 'claim': '2 pi e^(2L) at L = %s' % Lstr, 'chapter': hprinted, 'mine': '%.2f' % T_star, 'status': 'agrees' if ('%.1f' % T_star) == hprinted else 'differs'})
            kn = bc['K'] - Nst
            claims.append({'line': '103', 'claim': 'the last atom killed sits at the horizon 2 pi e^(2L) "to within the spacing of the zeros", L = %s' % Lstr, 'chapter': 'to within the spacing',
                           'mine': 'K = %d, N(T*) = %d, so %d killed atoms above T*; gamma_K minus T* = %.2f against the gap %.2f below gamma_K and %.2f above; R_upper at K = N(T*): %s' % (
                               bc['K'], Nst, kn, gd[bc['K'] - 1] - T_star, gd[bc['K'] - 1] - gd[bc['K'] - 2], gd[bc['K']] - gd[bc['K'] - 1], entry['R_upper_at_K_equal_N_T_star']),
                           'status': 'agrees' if kn <= 1 else 'differs'})
            claims.append({'line': '103', 'claim': '"killing an atom above the horizon costs more growth than the window kernel can absorb", L = %s' % Lstr,
                           'chapter': 'no atom above 2 pi e^(2L) is worth killing', 'mine': 'the optimum kills %d atom(s) above T*; R_upper(K = N(T*)) / R_upper(best) = %.3g' % (
                               kn, byK[Nst]['R_upper'] / bu['R_upper']), 'status': 'differs' if kn >= 1 else 'agrees'})
        summary.append((Lstr, bc, bu))
    if '0.8' in Ls:
        e = next(x for x in res['results'] if x['L'] == '0.8')
        geo = next((g for g in e.get('geometric', []) if g['K'] == int(e['best_upper']['K']) and g['m'] == int(e['best_upper']['m'])), None)
        b = float(e['best_central']['R_central'])
        claims.append({'line': '103', 'claim': 'product bound at L = 0.8 against the true 2e-17', 'chapter': '4.0e-9 against 2e-17', 'mine': '%s against Zhu upper bound 2.27e-17 (ratio %.2e)' % (sig(b, 3), b / 2.27e-17),
                       'status': 'agrees' if sig(b, 2) == '4.0e-09' else 'differs'})
    res['claims'] = claims
    for c in claims:
        log('CLAIM %-14s %-8s %s | chapter %s | mine %s' % (c['line'], c['status'], c['claim'], c['chapter'], c['mine']))
    if args.out:
        out = Path(args.out)
    elif quick:
        out = Path(tempfile.gettempdir()) / 'zh_product_bound_quick.json'
    else:
        out = Path('data/zh_product_bound.json') if Path('data').is_dir() else Path(__file__).resolve().parent / 'zh_product_bound.json'
    out.write_text(json.dumps(res, indent=1, default=str))
    log('wrote', out)


if __name__ == '__main__':
    main()
