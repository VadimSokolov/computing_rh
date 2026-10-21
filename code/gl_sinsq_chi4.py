# gl_sinsq_chi4.py: the numerical check of the sine squared lemma for the character chi_4, the paragraph "Numerical check for chi_4"
# of Section gl:sec:sinsq (Proposition gl:prop:sinsq, equation gl:eq:sinsq, Theorem gl:thm:sinsq) of Chapter 21 (Thorin
# Representations, Dirichlet L Functions and the Linnik Constant). Added for the book in October 2026, written from the mathematics
# (the first author's own scripts were never read, copied or run).
# Objects. chi_4 is the odd character modulo 4 (a = 1), xi(s) = xi(s, chi_4) = (4/pi)^((s+1)/2) Gamma((s+1)/2) L(s, chi_4), and for
# alpha > 1 pi M_alpha(b) = Re (xi'/xi)(alpha + i b) = (1/2) log(4/pi) + (1/2) Re psi((alpha+1+ib)/2) - sum_n chi_4(n) Lambda(n)
# n^(-alpha) cos(b log n), the prime power sum running over all n (the chapter's gl:eq:Malpha); every L value and every digamma
# value is computed by mpmath at 30 digits.
# Statements checked, in the order of the chapter's paragraph. (1) The kernel identity gl:eq:sinsq at (w, x) = (0.7, 1.3) and (2,
# 0.4): both sides at 30 digits, the integral over b by exact integration over each half period plus the asymptotic expansion of the
# tail of the oscillatory integral (K = 40 and K = 60 half periods, which agree). (2) log xi(alpha+w)/xi(alpha) against the integral
# over b of log(1 + w^2/b^2) M_alpha(b), alpha = 2, w = 1 and 3, with M_alpha from the prime powers up to N = 1e3, 1e4, 1e5, 1e6
# (and 2e6): the b integral of each prime term is the kernel identity, int_0^inf log(1+w^2/b^2) cos(b x) db = pi (1 - e^(-wx))/x, so
# the integral of the truncated M_alpha is known term by term and the difference to log xi(alpha+w)/xi(alpha) is the tail of log
# L(alpha+w)/L(alpha) (written D_N below); the gamma factor part is integrated numerically at 30 digits and compared with its closed
# form; one genuine numerical b quadrature of the whole truncated M_alpha (Gauss Legendre on panels, N = nq) tests the interchange;
# the loss from cutting the b integral at B is computed. (3) M_alpha(b) at alpha = 2 and 3/2, b = 0, 6.02, 10, 25, 50: the prime
# side (cutoffs up to 2e6) against the direct value pi^(-1) Re (xi'/xi)(alpha+ib) and against the zero side (1/pi) sum over the
# zeros of both signs of c/(c^2+(b-gamma)^2), c = alpha - 1/2, from the zeros computed here below T = 300, with the tail beyond T.
# (4) The minimum of pi M_alpha on 0 <= b <= 150 for alpha = 1, 3/4, 3/5, 11/20 (grid step 0.1 as in the chapter, and a grid of step
# 0.005), the exact value at b = 0 (direct, and zero sum with tail), the bound (alpha - 1/2) V_f = 2 (alpha - 1/2) p_1, the smallest
# local minimum away from the origin. (5) The peaks of F_N(b) = - sum over n <= N of chi_4(n) Lambda(n) n^(-1/2) cos(b log n)
# cos^2((pi/2) log n / log N), N = 2e6, on 0 <= b <= 45, grid step 0.005, against the zeros of L(s, chi_4): peak finder = strict
# local maxima above half the largest value, each refined by a parabola through three grid points and by bounded maximisation of F_N
# itself.
# Zero side and its tail. The positive ordinates are found as sign changes of the real function xi(1/2+it) on a grid of step 1/32
# (the method of code/dbeta_zeros.py, here in mpmath at 20 digits) and refined by Brent's method (the imaginary part of xi relative
# to its modulus is recorded: it is rounding only). With N^+(t) the number of zeros of height at most t and M(t) = (t/2pi) log(4t/(2
# pi e)) + 1/8, the tail over the zeros above T of a decreasing function g(gamma) is, by Stieltjes integration by parts, int_T^inf g
# M' dt - g(T) (N^+(T) - M(T)) + int_T^inf (N^+ - M)(-g') dt. The last term is bounded by int B(t) |g'| dt, with B = (0.22737 l + 2
# log(1+l) - 0.5)/2, l = log(4(t+2)/(2 pi)), half of the Bennett, Martin, O'Bryant and Rechnitzer bound for |N(t, chi) - ((t/pi)
# log(qt/(2 pi e)) - chi(-1)/4)| (Math. Comp. 90 (2021), the bound used in Theorem thm:dbeta:minall(v)). So each tail comes with an
# enclosure R that is rigorous if the zero list is complete up to T (checked here only by comparison with the book's certified table
# below height 100 and with M; not certified) and every zero above T satisfies GRH (for alpha > 1 this is unnecessary).
# Error control. mpmath works at 30 digits (the quadrature of the gamma part is compared with its closed form, which gives its true
# error); prime sums are float64 pairwise sums of at most 1.5e5 terms (error about 1e-14), so differences of size 1e-9 and above are
# meaningful. mpmath's quadosc is NOT used: on this integrand it is wrong by 8e-5 (checked against the closed form).
# Output: data/gl_sinsq_chi4.json (numbers as strings, with the zero list, the cutoffs and all quadrature parameters) and
# data/gl_sinsq_chi4.log (the printed lines). Parts can be run separately (--parts) and accumulated into one JSON file (--merge);
# every part records its own cutoffs and zero source. The comparison with the chapter's numbers is in
# misc/reviews/nick-linnik/checks-L3.md.
# Usage: python3 gl_sinsq_chi4.py [--quick] [--out FILE] [--parts 1,2,3,4,5] [--nmax N] [--tmax T] [--zeros compute|book|FILE]
# [--bmax B] [--nq N] [--workers W] [--merge]. --quick uses N = 1e4, zeros below 45 and b up to 20 in part 4; it takes 4 seconds in
# one process and writes to /tmp unless --out is given. The full case (N = 2e6, zeros below 300) is meant for the Hopper cluster:
#   bash misc/tools/hopper_run.sh -c 8 -t 60 -g gl_sinsq_chi4.json -g gl_sinsq_chi4.log code/gl_sinsq_chi4.py data/dbeta_zeros.json \
#        -- --workers 8 --out gl_sinsq_chi4.json
# (the cluster has the same mpmath, numpy and scipy; python-flint is not needed; the input file is only used to compare the zeros
# below 100). --zeros book reads the book's certified ordinates below 100 (data/dbeta_zeros.json) and is the way to run the parts
# that need zeros locally (parts 3 and 4 are then partial: T = 100, b up to 75).
import argparse
import json
import math
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'data'
mp.mp.dps = 30
CHI = [0, 1, 0, -1]
T_START = time.time()
LOGF = None


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGF is not None:
        LOGF.write(s + '\n')
        LOGF.flush()


def sf(x, d=12):
    """String with d significant digits of a float or an mpmath number."""
    if isinstance(x, (mp.mpf, mp.mpc)):
        return mp.nstr(x, d)
    return '%.*g' % (d, float(x))


def elapsed():
    return '%.1fs' % (time.time() - T_START)


# ---------------------------------------------------------------------------------------------------------------------------------
# functions of chi_4 at 30 digits
# ---------------------------------------------------------------------------------------------------------------------------------
def L4(s):
    return mp.dirichlet(s, CHI)


def dlogL4(s):
    return mp.dirichlet(s, CHI, 1) / mp.dirichlet(s, CHI)


def logxi4(s):
    """log xi(s, chi_4) for real s > 0 (xi is positive on the real axis)."""
    return (s + 1) / 2 * mp.log(4 / mp.pi) + mp.loggamma((s + 1) / 2) + mp.log(L4(s))


def dlogxi4(s):
    return mp.log(4 / mp.pi) / 2 + mp.psi(0, (s + 1) / 2) / 2 + dlogL4(s)


def agamma(alpha, b):
    """The gamma factor part (1/2) log(4/pi) + (1/2) Re psi((alpha+1+ib)/2) of pi M_alpha(b)."""
    return mp.log(4 / mp.pi) / 2 + mp.re(mp.psi(0, (alpha + 1 + 1j * b) / 2)) / 2


# ---------------------------------------------------------------------------------------------------------------------------------
# prime powers
# ---------------------------------------------------------------------------------------------------------------------------------
def iroot(N, k):
    r = int(round(N ** (1.0 / k)))
    while r ** k > N:
        r -= 1
    while (r + 1) ** k <= N:
        r += 1
    return r


def prime_powers(N):
    """All prime powers p^k <= N (n_all, k_all, sorted) and the odd ones with the weights a = chi_4(n) log p, x = log n."""
    N = int(N)
    sv = np.ones(N + 1, dtype=bool)
    sv[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if sv[p]:
            sv[p * p::p] = False
    primes = np.nonzero(sv)[0].astype(np.int64)
    ns, ks, ls = [], [], []
    k = 1
    while 2 ** k <= N:
        pr = primes[primes <= iroot(N, k)]
        ns.append(pr ** k)
        ks.append(np.full(len(pr), k, dtype=np.int64))
        ls.append(np.log(pr.astype(np.float64)))
        k += 1
    n_all = np.concatenate(ns)
    k_all = np.concatenate(ks)
    lp_all = np.concatenate(ls)
    o = np.argsort(n_all, kind='stable')
    n_all, k_all, lp_all = n_all[o], k_all[o], lp_all[o]
    odd = (n_all % 2 == 1)
    n, lp = n_all[odd], lp_all[odd]
    a = np.where(n % 4 == 1, 1.0, -1.0) * lp
    return dict(N=N, n=n, x=np.log(n.astype(np.float64)), a=a, n_all=n_all, k_all=k_all, nprimes=int(len(primes)))


# ---------------------------------------------------------------------------------------------------------------------------------
# parallel helper and the zeros of L(s, chi_4)
# ---------------------------------------------------------------------------------------------------------------------------------
def pmap(fn, items, workers):
    if workers <= 1:
        return [fn(i) for i in items]
    with get_context('fork').Pool(workers) as p:
        return p.map(fn, items, chunksize=1)


def zline(t):
    """xi(1/2 + it, chi_4) divided by the modulus of its gamma factor, and the ratio |Im|/modulus (rounding only)."""
    with mp.workdps(20):
        s = mp.mpc(mp.mpf(1) / 2, t)
        z = (s + 1) / 2
        gam = mp.exp(z * mp.log(4 / mp.pi) + mp.loggamma(z))
        val = gam * mp.dirichlet(s, CHI)
        ab = abs(gam)
        return float(val.real / ab), float(abs(val.imag) / ab)


def _zl_chunk(ts):
    return [zline(t) for t in ts]


def _refine(ab):
    a, b = ab
    return brentq(lambda t: zline(t)[0], a, b, xtol=1e-13, rtol=4 * np.finfo(float).eps, maxiter=200)


def compute_zeros(tmax, workers, h):
    n = int(tmax / h)
    ts = h * np.arange(1, n + 1)
    chunks = [list(ts[i:i + 40]) for i in range(0, n, 40)]
    vals = []
    imr = 0.0
    for r in pmap(_zl_chunk, chunks, workers):
        for v, im in r:
            vals.append(v)
            imr = max(imr, im)
    vals = np.array(vals)
    sg = np.sign(vals)
    br = [(ts[i], ts[i + 1]) for i in range(n - 1) if sg[i] * sg[i + 1] < 0]
    roots = np.array(sorted(pmap(_refine, br, workers)))
    return roots, imr


def book_zeros():
    f = DATA / 'dbeta_zeros.json'
    if not f.exists():
        f = HERE / 'dbeta_zeros.json'        # on the cluster the input file is copied beside the script
    if not f.exists():
        return None
    d = json.load(open(f))
    return np.array([float(s) for s in d['chi4']['ordinates']])


def Mcount(t, q=4):
    return t / (2 * np.pi) * np.log(q * t / (2 * np.pi * np.e)) + 1.0 / 8


def Mprime(t, q=4):
    return np.log(q * t / (2 * np.pi)) / (2 * np.pi)


def Bbound(t, q=4):
    l = np.log(q * (t + 2) / (2 * np.pi))
    return (0.22737 * l + 2 * np.log1p(l) - 0.5) / 2


def gfun(u, c):
    return c / (c * c + u * u)


def mgfun(u, c):
    """-g'(u)."""
    return 2 * c * u / (c * c + u * u) ** 2


def zero_part(bs, c, zeros):
    """Partial sum over the listed zeros of both signs of the Poisson kernel c/(c^2+(b-gamma)^2) (units of pi M), in chunks of b."""
    bs = np.atleast_1d(np.asarray(bs, dtype=float))
    out = np.empty(len(bs))
    for i in range(0, len(bs), 2000):
        b = bs[i:i + 2000, None]
        out[i:i + 2000] = (gfun(b - zeros[None, :], c) + gfun(b + zeros[None, :], c)).sum(axis=1)
    return out


def zero_side(bs, c, zeros, T, with_R=True):
    """Poisson smoothing at height c of the zeros of both signs: partial sum, smooth tail, boundary term, enclosure R (all in units of pi M)."""
    bs = np.atleast_1d(np.asarray(bs, dtype=float))
    part = zero_part(bs, c, zeros)
    main, _ = quad_vec(lambda t: Mprime(t) * (gfun(t - bs, c) + gfun(t + bs, c)), T, np.inf, epsabs=1e-15, epsrel=1e-11, norm='max')
    if with_R:
        R, _ = quad_vec(lambda t: Bbound(t) * (np.abs(mgfun(t - bs, c)) + np.abs(mgfun(t + bs, c))), T, np.inf, epsabs=1e-15, epsrel=1e-9,
                        norm='max')
    else:
        R = np.zeros(len(bs))
    ET = len(zeros) - Mcount(T)
    bnd = -ET * (gfun(T - bs, c) + gfun(T + bs, c))
    return dict(part=part, main=np.atleast_1d(main), bnd=bnd, R=np.atleast_1d(R))


# ---------------------------------------------------------------------------------------------------------------------------------
# part 1: the kernel identity
# ---------------------------------------------------------------------------------------------------------------------------------
def hder(b, w, j):
    """j-th derivative in b of h(b) = log(1 + w^2/b^2) = log(b+iw) + log(b-iw) - 2 log b."""
    if j == 0:
        return mp.log(1 + w ** 2 / b ** 2)
    z = mp.mpc(0, w)
    return (-1) ** (j - 1) * mp.factorial(j - 1) * mp.re((b + z) ** (-j) + (b - z) ** (-j) - 2 * b ** (-j))


def osc_int(w, x, K):
    """int_0^inf log(1+w^2/b^2) cos(b x) db: quadrature over K half periods plus the asymptotic expansion of the tail."""
    h = lambda b: mp.log(1 + w ** 2 / b ** 2)
    pts = [k * mp.pi / x for k in range(K + 1)]
    s = mp.mpf(0)
    for a, b in zip(pts[:-1], pts[1:]):
        s += mp.quad(lambda t: h(t) * mp.cos(t * x), [a, b])
    P = pts[-1]
    tail = mp.mpc(0)
    for j in range(0, 90):
        tail += (-1) ** (j + 1) * hder(P, w, j) / (1j * x) ** (j + 1)
    return s + mp.re(mp.exp(1j * P * x) * tail)


def plain_int(w, P=100):
    """int_0^inf log(1+w^2/b^2) db (= pi w): quadrature on [0, P] plus the series of the tail."""
    s = mp.quad(lambda b: mp.log(1 + w ** 2 / b ** 2), [0, 1, 10, P])
    t = sum((-1) ** (m + 1) * w ** (2 * m) / (m * (2 * m - 1) * mp.mpf(P) ** (2 * m - 1)) for m in range(1, 60))
    return s + t


def part1():
    log('== part 1: the kernel identity gl:eq:sinsq at 30 digits ==', elapsed())
    rows = []
    for ws, xs in (('0.7', '1.3'), ('2', '0.4')):
        w, x = mp.mpf(ws), mp.mpf(xs)
        lhs = mp.exp(-w * x) - 1 + w * x
        I1 = plain_int(w)
        out = {'w': ws, 'x': xs, 'lhs': mp.nstr(lhs, 25)}
        for K in (40, 60):
            I2 = osc_int(w, x, K)
            rhs = (x / mp.pi) * (I1 - I2)
            out['rhs_K%d' % K] = mp.nstr(rhs, 25)
            out['diff_K%d' % K] = mp.nstr(lhs - rhs, 4)
        out['int_log_db_minus_pi_w'] = mp.nstr(I1 - mp.pi * w, 4)
        rows.append(out)
        log('w=%s x=%s lhs=%s rhs=%s |lhs-rhs|=%s (K=40), %s (K=60)' % (ws, xs, out['lhs'], out['rhs_K60'], out['diff_K40'], out['diff_K60']))
    return {'rows': rows, 'method': 'int_0^inf log(1+w^2/b^2)(1-cos bx) db = I1 - I2, I1 by quadrature on [0,100] plus series, I2 by quadrature over K half '
            'periods plus asymptotic tail; mpmath dps 30; mp.quadosc NOT used (error 8e-5 on this integrand)'}


# ---------------------------------------------------------------------------------------------------------------------------------
# part 2: log xi(alpha+w)/xi(alpha) against the b integral of log(1+w^2/b^2) M_alpha
# ---------------------------------------------------------------------------------------------------------------------------------
def gl_prime_integral(pp, alpha, w, nq, B, nodes=12):
    """Numerical int_0^inf log(1+w^2/b^2) sum_{n<=nq} u_n cos(b log n) db, u_n = -chi(n) Lambda(n) n^(-alpha)/pi (Gauss Legendre panels)."""
    n, x, a = pp['n'], pp['x'], pp['a']
    k = int(np.searchsorted(n, nq, side='right'))
    xs = x[:k]
    u = -a[:k] * np.exp(-alpha * xs) / np.pi
    gx, gw = np.polynomial.legendre.leggauss(nodes)
    eps = 1e-10
    geo = eps * 10.0 ** (0.5 * np.arange(0, 21))
    uni = np.arange(1.0, B + 1e-9, 0.25)
    edges = np.concatenate([geo[:-1], uni])
    mid = (edges[:-1] + edges[1:]) / 2
    half = (edges[1:] - edges[:-1]) / 2
    bn = (mid[:, None] + half[:, None] * gx[None, :]).ravel()
    wn = (half[:, None] * gw[None, :]).ravel()
    total = 0.0
    step = max(1000, int(4e6 // max(k, 1)))
    for i in range(0, len(bn), step):
        bb = bn[i:i + step]
        cs = np.cos(np.outer(bb, xs)) @ u
        total += float(np.sum(wn[i:i + step] * np.log1p(w * w / (bb * bb)) * cs))
    first = eps * math.log1p(w * w / (eps * eps)) + 2 * w * math.atan(eps / w)
    total += float(u.sum()) * first
    hB = math.log1p(w * w / (B * B))
    dhB = -2 * w * w / (B * (B * B + w * w))
    tail = float(np.sum(u * (-hB * np.sin(B * xs) / xs - dhB * np.cos(B * xs) / (xs * xs))))
    return total + tail


def part2(pp, cuts, nq, Bq):
    log('== part 2: log xi(alpha+w)/xi(alpha) against int log(1+w^2/b^2) M_alpha db, alpha = 2 ==', elapsed())
    alpha = mp.mpf(2)
    af = 2.0
    n, x, a = pp['n'], pp['x'], pp['a']
    n_all, k_all = pp['n_all'], pp['k_all']
    logzeta = mp.log(mp.zeta(alpha))
    out = {'alpha': '2', 'cuts': [int(c) for c in cuts], 'rows': [],
           'method': 'D_N = (b integral of M^N against log(1+w^2/b^2), each prime term integrated exactly by the kernel identity) minus log xi(alpha+w)/xi(alpha); '
                     'gamma part: mp.quad (tanh sinh, intervals 0,1,10,100,1000,inf, 30 digits) against its closed form; absolute bound = log zeta(alpha) minus the sum of '
                     'p^(-k alpha)/k over p^k <= N; b cut loss: mp.quad of the gamma part integrand from B to infinity (15 digits); genuine quadrature of the truncated M: '
                     'Gauss Legendre, 12 nodes per panel, geometric panels (factor sqrt 10) from 1e-10 to 1, panels of width 0.25 up to B, analytic first panel, '
                     'tail beyond B by two terms of the asymptotic expansion of the oscillatory integral'}
    for ws in (1, 3):
        w = mp.mpf(ws)
        wf = float(ws)
        ref = logxi4(alpha + w) - logxi4(alpha)
        Gex = w / 2 * mp.log(4 / mp.pi) + mp.loggamma((alpha + 1 + w) / 2) - mp.loggamma((alpha + 1) / 2)
        logLr = mp.log(L4(alpha + w) / L4(alpha))
        gint = lambda b: mp.log(1 + w ** 2 / b ** 2) * (mp.log(4 / mp.pi) + mp.re(mp.psi(0, (alpha + 1 + 1j * b) / 2))) / (2 * mp.pi)
        Gnum, Gerr = mp.quad(gint, [0, 1, 10, 100, 1000, mp.inf], error=True)
        row = {'w': ws, 'log_xi_ratio': mp.nstr(ref, 25), 'gamma_part_exact': mp.nstr(Gex, 22), 'gamma_part_numeric': mp.nstr(Gnum, 22),
               'gamma_part_diff': mp.nstr(Gnum - Gex, 4), 'gamma_part_quad_error_estimate': mp.nstr(Gerr, 3), 'cuts': {}}
        log('w=%d: log xi(alpha+w)/xi(alpha) = %s; gamma part numeric minus closed form = %s (quadrature error estimate %s)' %
            (ws, mp.nstr(ref, 20), row['gamma_part_diff'], row['gamma_part_quad_error_estimate']))
        for N in cuts:
            i = int(np.searchsorted(n, N, side='right'))
            PN = float(np.sum(-a[:i] * np.exp(-af * x[:i]) * (-np.expm1(-wf * x[:i])) / x[:i]))
            DN = mp.mpf(PN) - logLr
            j = int(np.searchsorted(n_all, N, side='right'))
            Sall = float(np.sum(np.exp(-af * np.log(n_all[:j].astype(np.float64))) / k_all[:j]))
            absb = float(logzeta) - Sall
            row['cuts'][str(int(N))] = {'D_N_integral_minus_log_ratio': sf(DN, 6), 'abs_tail_bound_sum_Lambda_over_log_n_n_minus_alpha': sf(absb, 4)}
            log('   N=%8d  D_N = (integral of M^N) - log ratio = %s   absolute value bound on the truncation = %s' % (N, sf(DN, 6), sf(absb, 4)))
        # chapter claims for comparison
        row['chapter_claim_at_1e5'] = {1: '-1.16e-5', 3: '-3.3e-5'}[ws]
        # loss from cutting the b integral at B (gamma part only; the prime part loss is oscillatory and smaller)
        row['b_cut_loss_gamma_part'] = {}
        with mp.workdps(15):
            for Bc in (100, 1000, 10000, 100000, 1000000):
                v = mp.quad(gint, [Bc, 10 * Bc, mp.inf])
                row['b_cut_loss_gamma_part'][str(Bc)] = mp.nstr(v, 5)
        log('   loss from cutting the b integral at B = 1e2, 1e3, 1e4, 1e5, 1e6 (gamma part): ' +
            ', '.join(row['b_cut_loss_gamma_part'][str(b)] for b in (100, 1000, 10000, 100000, 1000000)))
        # one genuine numerical quadrature of the whole truncated M_alpha
        t0 = time.time()
        pnum = gl_prime_integral(pp, af, wf, nq, Bq)
        i = int(np.searchsorted(n, nq, side='right'))
        pcl = float(np.sum(-a[:i] * np.exp(-af * x[:i]) * (-np.expm1(-wf * x[:i])) / x[:i]))
        row['numeric_b_quadrature'] = {'nq': int(nq), 'B': Bq, 'prime_part_numeric': sf(pnum, 15), 'prime_part_closed_form': sf(pcl, 15),
                                       'diff': sf(pnum - pcl, 4), 'total_numeric_minus_log_ratio': mp.nstr(Gnum + mp.mpf(pnum) - ref, 6),
                                       'total_closed_minus_log_ratio_at_nq': mp.nstr(Gex + mp.mpf(pcl) - ref, 6)}
        log('   numerical b quadrature at N=%d, b up to %g (%.1fs): prime part %s against closed form %s, difference %s' %
            (nq, Bq, time.time() - t0, sf(pnum, 12), sf(pcl, 12), sf(pnum - pcl, 4)))
        out['rows'].append(row)
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
# part 3: M_alpha(b) from the primes, directly, and from the zeros
# ---------------------------------------------------------------------------------------------------------------------------------
def prime_sum(pp, alpha, b, N):
    n, x, a = pp['n'], pp['x'], pp['a']
    i = int(np.searchsorted(n, N, side='right'))
    xs = x[:i]
    return float(np.sum(a[:i] * np.exp(-alpha * xs) * np.cos(b * xs)))


def part3(pp, cuts, zeros, T, zsrc):
    log('== part 3: prime side, direct value and zero side of M_alpha(b), chi_4 ==', elapsed())
    Nplus = len(zeros)
    out = {'cuts': [int(c) for c in cuts], 'T': T, 'zero_source': zsrc, 'count_zeros_below_T': Nplus, 'M_of_T': sf(Mcount(T), 8),
           'E_T': sf(Nplus - Mcount(T), 6), 'rows': [],
           'method': 'differences in M units (divide pi M by pi); prime side = (gamma part - float64 pairwise sum over prime powers up to N)/pi; direct = mpmath 30 digits '
                     '(digamma and L\'/L from mp.dirichlet with derivative 1); zero side = partial sum over the listed zeros of both signs, + smooth tail int_T^inf M\'(t)(g(t-b)+g(t+b))dt, '
                     '+ boundary term -(N(T)-M(T))(g(T-b)+g(T+b)), g(u) = c/(c^2+u^2); enclosure R = int_T^inf B(t)(|g\'(t-b)|+|g\'(t+b)|)dt (Bennett Martin O\'Bryant Rechnitzer)'}
    bs = [0.0, 6.02, 10.0, 25.0, 50.0]
    for alpha_s, af in (('2', 2.0), ('3/2', 1.5)):
        alpha = mp.mpf(af)
        c = af - 0.5
        ok = [b for b in bs if b <= T - 30]
        zs = zero_side(ok, c, zeros, T)
        for ib, b in enumerate(ok):
            bm = mp.mpf(b)
            direct = mp.re(dlogxi4(alpha + 1j * bm))      # pi M_alpha(b)
            ag = float(agamma(alpha, bm))
            dm = float(direct)
            row = {'alpha': alpha_s, 'b': sf(b, 6), 'pi_M_direct': mp.nstr(direct, 22), 'M_direct': mp.nstr(direct / mp.pi, 22), 'prime_side': {}}
            for N in cuts:
                ps = (ag - prime_sum(pp, af, b, N)) / math.pi
                row['prime_side'][str(int(N))] = sf(ps - dm / math.pi, 4)
            zA = zs['part'][ib]
            zB = zA + zs['main'][ib]
            zC = zB + zs['bnd'][ib]
            row['zero_side_minus_direct_in_M'] = {'partial_sum_only': sf((zA - dm) / math.pi, 4), 'with_smooth_tail': sf((zB - dm) / math.pi, 4),
                                                  'with_tail_and_boundary_term': sf((zC - dm) / math.pi, 4), 'enclosure_half_width_R': sf(zs['R'][ib] / math.pi, 4)}
            row['zero_side_pi_M'] = sf(zC, 15)
            out['rows'].append(row)
            nmx = max(cuts)
            log('alpha=%-4s b=%-6s M_direct=%s  prime side minus direct (M units) at N=%s;  zero side minus direct: partial %s, +tail %s, +boundary %s (R=%s)' %
                (alpha_s, sf(b, 6), mp.nstr(direct / mp.pi, 12), ', '.join('%d: %s' % (N, row['prime_side'][str(int(N))]) for N in cuts),
                 row['zero_side_minus_direct_in_M']['partial_sum_only'], row['zero_side_minus_direct_in_M']['with_smooth_tail'],
                 row['zero_side_minus_direct_in_M']['with_tail_and_boundary_term'], row['zero_side_minus_direct_in_M']['enclosure_half_width_R']))
        # scaling of the prime side error with N at b = 0 and the absolute bound
        out.setdefault('abs_tail_bound_formula', {})[alpha_s] = {str(int(N)): sf(N ** (1 - af) / (af - 1), 3) for N in cuts}
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
# part 4: the minimum of pi M_alpha on 0 <= b <= bmax
# ---------------------------------------------------------------------------------------------------------------------------------
def part4(pp, cuts, zeros, T, bmax, zsrc):
    log('== part 4: the minima of pi M_alpha on [0, %g] for alpha = 1, 3/4, 3/5, 11/20 ==' % bmax, elapsed())
    Nplus = len(zeros)
    claims = {'1': '0.077784', '3/4': '0.038978', '3/5': '0.015601', '11/20': '0.007801'}
    fracs = {'1': mp.mpf(1), '3/4': mp.mpf(3) / 4, '3/5': mp.mpf(3) / 5, '11/20': mp.mpf(11) / 20}
    out = {'bmax': bmax, 'T': T, 'zero_source': zsrc, 'count_zeros_below_T': Nplus, 'rows': [],
           'method': 'pi M_alpha(0) direct = Re (xi\'/xi)(alpha) at 30 digits; grids of step 0.1 and 0.005 on [0, bmax] by the zero side (partial sum + smooth tail + boundary term), '
                     'local minima of the fine grid with b > 0.5 refined by bounded minimisation (xatol 1e-9) on the partial sum, the tail added for the five lowest; '
                     'second derivative at 0 = 2 sum over the listed zeros of phi\'\'(gamma), phi(u) = c/(c^2+u^2)'}
    bg = np.round(np.arange(0, int(round(bmax / 0.1)) + 1) * 0.1, 10)
    bf = np.round(np.arange(0, int(round(bmax / 0.005)) + 1) * 0.005, 10)
    # p_1 = sum over gamma > 0 of gamma^-2: the listed zeros plus the tail by the same method (g = t^-2), enclosure int B 2 t^-3 dt
    p1_part = float(np.sum(zeros ** -2.0))
    p1_main = float(quad_vec(lambda t: np.array([Mprime(t) / t ** 2]), T, np.inf, epsabs=1e-16, epsrel=1e-11)[0][0])
    p1_bnd = -(Nplus - Mcount(T)) / T ** 2
    p1_R = float(quad_vec(lambda t: np.array([Bbound(t) * 2 / t ** 3]), T, np.inf, epsabs=1e-16, epsrel=1e-9)[0][0])
    p1 = p1_part + p1_main + p1_bnd
    out['p1'] = {'value': sf(p1, 12), 'enclosure_R': sf(p1_R, 3), 'listed_zeros_only': sf(p1_part, 12), 'V_f_book': '0.1560296499', '2_p1': sf(2 * p1, 12)}
    log('p_1 = sum gamma^-2 = %s (listed zeros %s, tail %s, enclosure %s); 2 p_1 = %s against V_f = 0.1560296499 of the book' %
        (sf(p1, 12), sf(p1_part, 12), sf(p1_main + p1_bnd, 5), sf(p1_R, 3), sf(2 * p1, 12)))
    for key, alpha in fracs.items():
        af = float(alpha)
        c = af - 0.5
        t0 = time.time()
        exact = mp.re(dlogxi4(alpha))                    # pi M_alpha(0)
        z0 = zero_side([0.0], c, zeros, T)
        zsum0 = float(z0['part'][0] + z0['main'][0] + z0['bnd'][0])
        Rb0 = float(z0['R'][0])
        zg = zero_side(bg, c, zeros, T, with_R=False)
        vg = zg['part'] + zg['main'] + zg['bnd']
        zf = zero_side(bf, c, zeros, T, with_R=False)
        vf = zf['part'] + zf['main'] + zf['bnd']
        ig = int(np.argmin(vg))
        jf = int(np.argmin(vf))
        kg = int(np.nonzero(bg >= 1.0)[0][0])
        kf = int(np.nonzero(bf >= 1.0)[0][0])
        ig1 = kg + int(np.argmin(vg[kg:]))
        jf1 = kf + int(np.argmin(vf[kf:]))
        # the local minima of the fine grid away from the origin, refined on the partial sum; the tail is added for the five lowest
        loc = [j for j in range(1, len(vf) - 1) if vf[j] < vf[j - 1] and vf[j] <= vf[j + 1] and bf[j] > 0.5]
        ref = []
        for j in loc:
            r = minimize_scalar(lambda bb: float(zero_part([bb], c, zeros)[0]), bounds=(bf[j] - 0.005, bf[j] + 0.005), method='bounded',
                                options={'xatol': 1e-9})
            ref.append((float(r.x), float(r.fun)))
        ref.sort(key=lambda r: r[1])
        cand = []
        for (bb, vv) in ref[:5]:
            zz = zero_side([bb], c, zeros, T, with_R=False)
            cand.append((bb, vv + float(zz['main'][0] + zz['bnd'][0])))
        cand.sort(key=lambda r: r[1])
        d2 = 2 * float(np.sum(2 * c * (3 * zeros ** 2 - c * c) / (c * c + zeros ** 2) ** 3))     # pi M'' (0), listed zeros only
        ag0 = float(agamma(alpha, mp.mpf(0)))
        ps0 = {str(int(N)): sf(ag0 - prime_sum(pp, af, 0.0, N) - float(exact), 4) for N in cuts}    # plain prime side at b = 0 minus direct (pi M units)
        row = {'alpha': key, 'c': sf(c, 6), 'pi_M_at_0_direct': mp.nstr(exact, 22), 'pi_M_at_0_zero_sum': sf(zsum0, 12), 'zero_sum_enclosure_R': sf(Rb0, 3),
               'zero_sum_minus_direct': sf(zsum0 - float(exact), 4), 'chapter_value': claims[key],
               'two_c_p1': sf(2 * c * p1, 12), 'bound_(alpha-1/2)V_f_book': sf(c * 0.1560296499, 10), 'direct_over_2c_p1': sf(float(exact) / (2 * c * p1), 8),
               'grid_step_0.1': {'min': sf(vg[ig], 10), 'argmin': sf(bg[ig], 6), 'min_over_b_ge_1': sf(vg[ig1], 10), 'argmin_b_ge_1': sf(bg[ig1], 6)},
               'grid_step_0.005': {'min': sf(vf[jf], 10), 'argmin': sf(bf[jf], 6), 'min_over_b_ge_1': sf(vf[jf1], 10), 'argmin_b_ge_1': sf(bf[jf1], 6)},
               'smallest_local_minimum_away_from_origin': ({'b': sf(cand[0][0], 8), 'pi_M': sf(cand[0][1], 10), 'ratio_to_value_at_0': sf(cand[0][1] / float(exact), 6)}
                                                           if cand else None),
               'n_local_minima_b_gt_0.5': len(loc), 'second_derivative_at_0_listed_zeros': sf(d2, 8), 'plain_prime_side_at_0_minus_direct_pi_M': ps0}
        out['rows'].append(row)
        sm = row['smallest_local_minimum_away_from_origin']
        log('alpha=%-5s direct pi M(0)=%s  zero sum=%s (R=%s, minus direct %s)  chapter %s;  grid 0.1: min %s at b=%s (b>=1: %s at %s);  grid 0.005: min at b=%s '
            '(b>=1: %s at %s);  smallest local min away from 0: %s at b=%s (ratio %s);  2c p1=%s;  pi M2(0)=%s;  %.1fs' %
            (key, mp.nstr(exact, 14), sf(zsum0, 12), sf(Rb0, 2), sf(zsum0 - float(exact), 3), claims[key], sf(vg[ig], 8), sf(bg[ig], 4), sf(vg[ig1], 6), sf(bg[ig1], 5),
             sf(bf[jf], 4), sf(vf[jf1], 6), sf(bf[jf1], 6), sm['pi_M'] if sm else '-', sm['b'] if sm else '-', sm['ratio_to_value_at_0'] if sm else '-',
             sf(2 * c * p1, 10), sf(d2, 6), time.time() - t0))
        log('          plain prime side at b = 0 minus direct (pi M units) by cutoff: ' + ', '.join('%s: %s' % (k, v) for k, v in ps0.items()))
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
# part 5: the peaks of the smoothed prime sum against the zeros
# ---------------------------------------------------------------------------------------------------------------------------------
_G5 = {}


def _f5_chunk(rng):
    i0, i1, h = rng
    x, wg = _G5['x'], _G5['w']
    out = np.empty(i1 - i0)
    buf = np.empty_like(x)
    for k, i in enumerate(range(i0, i1)):
        np.multiply(x, i * h, out=buf)
        np.cos(buf, out=buf)
        out[k] = float(np.dot(wg, buf))
    return out


def part5(pp, N, zeros, workers):
    log('== part 5: the peaks of F_N(b) = - sum chi Lambda n^(-1/2) cos(b log n) cos^2(pi log n/(2 log N)), N = %d ==' % N, elapsed())
    n, x, a = pp['n'], pp['x'], pp['a']
    i = int(np.searchsorted(n, N, side='right'))
    xs = x[:i]
    L = math.log(N)
    wg = -a[:i] * np.exp(-0.5 * xs) * np.cos(np.pi / 2 * xs / L) ** 2
    _G5['x'], _G5['w'] = xs, wg
    h = 0.005
    bmax = 45.0
    J = int(round(bmax / h))
    chunks = [(s, min(s + 300, J + 1), h) for s in range(0, J + 1, 300)]
    F = np.concatenate(pmap(_f5_chunk, chunks, workers))
    bg = h * np.arange(J + 1)
    Fs = lambda b: float(np.dot(wg, np.cos(b * xs)))
    thr = 0.5 * F.max()
    js = [j for j in range(1, J) if F[j] > F[j - 1] and F[j] >= F[j + 1] and F[j] > thr and bg[j] > 1.0]
    zs = zeros[zeros < bmax]
    rows = []
    peaks_exact = []
    for j in js:
        y0, y1, y2 = F[j - 1], F[j], F[j + 1]
        dl = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
        bpar = bg[j] + dl * h
        r = minimize_scalar(lambda b: -Fs(b), bounds=(bg[j] - h, bg[j] + h), method='bounded', options={'xatol': 1e-10})
        peaks_exact.append((bg[j], bpar, float(r.x), -float(r.fun)))
    # FWHM of the first and last peaks: width at half of (peak - nearby trough)
    def fwhm(j):
        lo, hi = max(0, j - 200), min(J, j + 200)
        tr = F[lo:hi + 1].min()
        half = (F[j] + tr) / 2
        l = j
        while l > 0 and F[l] > half:
            l -= 1
        r_ = j
        while r_ < J and F[r_] > half:
            r_ += 1
        return (r_ - l) * h
    pk = np.array([p[2] for p in peaks_exact])
    maxg = maxp = maxe = 0.0
    for z in zs:
        k = int(np.argmin(np.abs(pk - z)))
        dg, dp, de = peaks_exact[k][0] - z, peaks_exact[k][1] - z, peaks_exact[k][2] - z
        maxg, maxp, maxe = max(maxg, abs(dg)), max(maxp, abs(dp)), max(maxe, abs(de))
        rows.append({'zero': sf(z, 11), 'peak_grid': sf(peaks_exact[k][0], 8), 'peak_parabola': sf(peaks_exact[k][1], 9), 'peak_exact': sf(peaks_exact[k][2], 10),
                     'dev_grid': sf(dg, 3), 'dev_parabola': sf(dp, 3), 'dev_exact': sf(de, 3), 'height': sf(peaks_exact[k][3], 8)})
    heights = [p[3] for p in peaks_exact]
    n40 = int((zeros < 40).sum())
    out = {'method': 'F_N on the grid b = 0.005 j, j = 0..9000 (direct cosines, float64 dot products); peaks = strict local maxima of the grid values above half the largest value with b > 1; '
                     'parabola through three grid points; exact = bounded maximisation of F_N (xatol 1e-10); each zero matched to the nearest exact peak; FWHM at half height above the lowest value within 1 of the peak',
           'N': int(N), 'terms': int(i), 'grid_step': h, 'bmax': bmax, 'peak_threshold': sf(thr, 6), 'number_of_peaks': len(js), 'zeros_below_45': int(len(zs)),
           'zeros_below_40': n40, 'first_zero': sf(zeros[0], 9), 'last_zero_below_40': sf(zeros[zeros < 40][-1], 9),
           'max_dev_grid': sf(maxg, 4), 'max_dev_parabola': sf(maxp, 4), 'max_dev_exact': sf(maxe, 4),
           'height_min': sf(min(heights), 8), 'height_max': sf(max(heights), 8), 'fwhm_first_peak': sf(fwhm(js[0]), 4), 'fwhm_last_peak': sf(fwhm(js[-1]), 4),
           'L_over_2': sf(L / 2, 8), 'rows': rows}
    log('N=%d, %d terms, grid step %g on [0, %g]: %d peaks above %s for %d zeros below 45 (%d below 40, from %s to %s)' %
        (N, i, h, bmax, len(js), sf(thr, 5), len(zs), n40, out['first_zero'], out['last_zero_below_40']))
    for r in rows:
        log('  zero %-13s grid peak %-9s parabola %-12s exact %-13s dev: grid %-9s parabola %-9s exact %-9s height %s' %
            (r['zero'], r['peak_grid'], r['peak_parabola'], r['peak_exact'], r['dev_grid'], r['dev_parabola'], r['dev_exact'], r['height']))
    log('largest deviation: grid peak %s, parabola %s, exact maximum %s; peak heights %s to %s (L/2 = %s); FWHM of the first and last peak %s, %s; %s' %
        (out['max_dev_grid'], out['max_dev_parabola'], out['max_dev_exact'], out['height_min'], out['height_max'], out['L_over_2'],
         out['fwhm_first_peak'], out['fwhm_last_peak'], elapsed()))
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
def sanity():
    one = mp.mpf(1)
    half = mp.mpf(1) / 2
    s = mp.mpc('0.3', '4.1')
    return {'dlogxi4(1)': mp.nstr(mp.re(dlogxi4(one)), 18), 'book (ch16)': '0.0777839899618', 'xi4(1/2)': mp.nstr(mp.exp(logxi4(half)), 15),
            'book xi(1/2)': '0.98071361405771', 'xi4(1)': mp.nstr(mp.exp(logxi4(one)), 15), 'L4(1/2)': mp.nstr(L4(half), 10),
            'functional_equation_residual_at_0.3+4.1i': mp.nstr(abs(mp.exp(((s + 1) / 2) * mp.log(4 / mp.pi) + mp.loggamma((s + 1) / 2)) * L4(s) -
                                                                  mp.exp(((2 - s) / 2) * mp.log(4 / mp.pi) + mp.loggamma((2 - s) / 2)) * L4(1 - s)), 3)}


def main():
    global LOGF
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--parts', default='1,2,3,4,5')
    ap.add_argument('--nmax', type=int, default=None)
    ap.add_argument('--tmax', type=float, default=None)
    ap.add_argument('--zeros', default='compute')
    ap.add_argument('--bmax', type=float, default=None)
    ap.add_argument('--nq', type=int, default=None)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--merge', action='store_true')
    args = ap.parse_args()
    quick = args.quick
    nmax = args.nmax or (10000 if quick else 2000000)
    tmax = args.tmax or (45.0 if quick else 300.0)
    bmax4 = args.bmax or (20.0 if quick else 150.0)
    nq = args.nq or (100 if quick else 1000)
    Bq = 2000.0 if quick else 20000.0
    parts = set(int(p) for p in args.parts.split(','))
    out = Path(args.out) if args.out else (Path('/tmp/gl_sinsq_chi4_quick.json') if quick else DATA / 'gl_sinsq_chi4.json')
    LOGF = open(out.with_suffix('.log'), 'a' if args.merge else 'w')
    log('gl_sinsq_chi4.py', ' '.join(sys.argv[1:]), '| mode', 'quick' if quick else 'full', '| nmax', nmax, '| tmax', tmax, '| zeros', args.zeros,
        '| workers', args.workers, '| mpmath', mp.__version__, 'numpy', np.__version__)
    cuts = sorted(set([int(c) for c in (1e3, 1e4, 1e5, 1e6, 2e6) if c <= nmax] + [nmax]))
    res = {}
    res['sanity'] = sanity()
    log('sanity:', json.dumps(res['sanity']))
    t0 = time.time()
    pp = prime_powers(nmax)
    log('prime powers up to %d: %d primes, %d prime powers, %d odd (%.1fs)' % (nmax, pp['nprimes'], len(pp['n_all']), len(pp['n']), time.time() - t0))
    zeros = T = zsrc = None
    if parts & {3, 4, 5}:
        book = book_zeros()
        if args.zeros == 'book':
            zeros, T, zsrc = book, 100.0, 'data/dbeta_zeros.json (chi4 ordinates below 100, certified count, 10 digits)'
        elif args.zeros != 'compute':
            zeros = np.sort(np.array(json.load(open(args.zeros)) if args.zeros.endswith('.json') else np.loadtxt(args.zeros), dtype=float))
            T, zsrc = tmax, args.zeros
        else:
            t0 = time.time()
            tz = tmax if parts & {3, 4} else 50.0
            zeros, imr = compute_zeros(tz, args.workers, 1 / 16 if quick else 1 / 32)
            T, zsrc = tz, 'computed here: sign changes of xi(1/2+it) on a grid of step %s, Brent refinement to 1e-13' % ('1/16' if quick else '1/32')
            e = np.arange(1, len(zeros) + 1) - 0.5 - Mcount(zeros)
            log('zeros below %g: %d found in %.1fs; main term M(T) = %.4f, N(T) - M(T) = %.4f; max |k - 1/2 - M(gamma_k)| = %.4f (<= B = %.3f); max |Im xi|/|xi| on the grid = %.2e' %
                (T, len(zeros), time.time() - t0, Mcount(T), len(zeros) - Mcount(T), np.abs(e).max(), Bbound(T), imr))
            res['zero_checks'] = {'T': T, 'count': int(len(zeros)), 'M(T)': sf(Mcount(T), 8), 'max_abs_k_minus_half_minus_M': sf(np.abs(e).max(), 5),
                                  'max_imag_over_modulus_on_grid': sf(imr, 3), 'first': sf(zeros[0], 12), 'count_below_40': int((zeros < 40).sum())}
            if book is not None:
                m = min(len(book), int((zeros < 100).sum()))
                res['zero_checks']['book_comparison'] = {'count_below_100': int((zeros < 100).sum()), 'book_count': int(len(book)),
                                                         'max_abs_diff': sf(np.abs(zeros[:m] - book[:m]).max(), 3)}
                log('comparison with the book table below 100: counts %d and %d, largest difference %.2e' % ((zeros < 100).sum(), len(book), np.abs(zeros[:m] - book[:m]).max()))
        res['zeros'] = [sf(z, 14) for z in zeros]
        res['zero_source'] = zsrc
        log('zero list: %d zeros, T = %g, source: %s' % (len(zeros), T, zsrc))
    if 1 in parts:
        res['part1'] = part1()
    if 2 in parts:
        res['part2'] = part2(pp, cuts, nq, Bq)
    if 3 in parts:
        res['part3'] = part3(pp, cuts, zeros, T, zsrc)
    if 4 in parts:
        res['part4'] = part4(pp, cuts, zeros, T, min(bmax4, T - 25.0 if T < 120 else T - 50.0), zsrc)
    if 5 in parts:
        res['part5'] = part5(pp, nmax, zeros, args.workers)
    res['run'] = {'argv': [a.replace(str(HERE.parent) + '/', '') for a in sys.argv[1:]], 'quick': quick, 'nmax': nmax, 'tmax': tmax, 'workers': args.workers, 'parts': sorted(parts), 'seconds': round(time.time() - T_START, 1)}
    old = {}
    if args.merge and out.exists():
        old = json.load(open(out))
        runs = old.get('runs', [])
    else:
        runs = []
    runs.append(res.pop('run'))
    old.update(res)
    old['runs'] = runs
    json.dump(old, open(out, 'w'), indent=1)
    log('wrote', out, elapsed())


if __name__ == '__main__':
    main()
