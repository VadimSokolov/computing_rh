# dc_weil.py: the Weil side of the chapter "Can the Deconvolution Be Done?" (label prefix dc:), group D1 of the number checks:
# the Gaussian check of the explicit formula (the sentence after the proof of dc:thm:weil), the constants of Section dc:sec:numerics
# (2 + gamma_E - log 4 pi, m(y), V_xi, p_1, the margin 0.09 y, the column pi^2 w_y''(0)/4 of Table dc:tab:rate), the identity
# dc:eq:Qw without zeros for the cosine test function, and the first identity of dc:thm:weil for a triangle Fourier transform.
# Added for the book in October 2026, written from the mathematics (no script of the chapter's author was read or run).
# Conventions. fhat(t) = int f(x) e^{itx} dx. For an even g, W(g) = ghat(i/2) + ghat(-i/2) - (gamma_E + log pi) g(0)
# + int_0^inf 2 (e^{-2x} g(0) - e^{-x/2} g(x))/(1 - e^{-2x}) dx - 2 sum_n Lambda(n) n^{-1/2} g(log n), which is the sum over the
# zeros of ghat(gamma_rho). Q_y(h) = W(e^{-y|x|} (h * htilde)). xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s), so that
# xi'/xi = 1/s + 1/(s-1) - (log pi)/2 + psi(s/2)/2 + zeta'/zeta. w_y(theta) = Re xi'/xi(1/2 + y + i theta), v_y = w_y/pi, m(y) = 2 xi'/xi(1/2 + y).
# Method. W1: g = exp(-(x/c)^2), c = 1/2 and 1, ghat(t) = c sqrt(pi) e^{-c^2 t^2/4}. The prime side is the formula above with the
# closed form tail g(0) (-log(1 - e^{-2X})) of the g(0) part beyond X (X = 12 or more), mp.quad on unit intervals, prime powers while
# g(log n) > 10^-(dps+20); the zero side is 2 sum_{gamma>0} ghat(gamma) over data/zeros_1700.txt. The prime side with the tail dropped
# at X = 12 is also written (the convention of a truncated script). W2: every constant at 50 digits (30 in --quick) from mpmath; the
# second derivative of xi'/xi at 1/2 + y three ways (the explicit formula in zeta', zeta'', zeta''', a Cauchy circle of radius 1/4,
# and at s = 1 the Stieltjes constants); V_xi also by mp.diff; p_1 from the 1700 zeros plus the tail int_{gamma_N}^inf t^-2 dNbar
# - S(gamma_N+)/gamma_N^2 with Nbar = theta/pi + 1 (Riemann von Mangoldt), also from 850 zeros. W3: h(u) = cos(pi u/2a) on [-a, a],
# a = 1; (h*h)(x) = (2a - |x|) cos(kx)/2 + sin(k|x|)/(2k), k = pi/2a; |hhat|^2 = 4 k^2 cos^2(a theta)/(k^2 - theta^2)^2. The prime
# side is W of e^{-y|x|}(h*h) (finitely many primes, n <= e^{2a}); the other side is (1/pi) int |hhat|^2 w_y dtheta by the scheme below.
# W4: phihat = triangle (1 - |x|/2)_+, g = phihat/2pi, phi(theta) = sin^2(theta)/(pi theta^2); P_{y0}^{-1} phi has Fourier transform
# phihat e^{y0|x|}, so (P^{-1} phi)(theta) = pi^-1 Re[(e^{2z} - 1)/(2 z^2) - 1/z], z = y0 + i theta; D = int v_{y0} (P^{-1} phi) dtheta
# at y0 = 1/2 and 1, against W(g); also the zero sum 2 sum_{gamma>0} phi(gamma) with the same kind of tail as for p_1.
# The theta integral (W3, W4). Write w = A + Z, A the gamma factor and pole terms (digamma), Z = Re zeta'/zeta; the weight Phi = Phi_s + Phi_o
# splits into a smooth and an oscillating part (frequency 2a, resp. 2) beyond T_s. Then (pi/2) L = int_0^{T_s} A Phi + int_{T_s}^inf A Phi_s
# (mp.quad, closed tail) + int_{T_s}^{T_1} chi A Phi_o + int_0^{T_1} chi Z Phi, chi = 1 up to T_0 and a smooth step down to 0 at T_1, all on
# Gauss Legendre panels of half width 0.8 (y) (the distance to the nearest pole of xi'/xi is y), L = (2/pi)(that sum) + L_N. Z has no
# constant Fourier component, so the taper loses only terms of size (omega Delta)^-4. Dirichlet terms: for sigma > 1 and in the sense
# of the Dirichlet series also at sigma = 1, Z_N = -sum_{n<=N} Lambda(n) n^-sigma cos(theta log n) integrates against Phi in closed form,
# int Phi Z_N = -sum Lambda(n) n^-sigma G(log n), G the Fourier transform of Phi (hybrid mode, N = 100 or more: it removes the slow
# components of Z, for which the taper needs a long window); N = 0 is the pure brute force that uses no prime at all.
# Error control. (a) the Gauss Legendre panels on [0, T_s] against panels of half the width; (b) the closed forms of hhat, h*h and of
# P^{-1} phi against direct quadrature at three points; (c) mp.quad error estimates of the tail integral; (d) every theta integral
# for several windows and tapers (smoother step, raised cosine) and for N = 0 and N in Nex; the spread is the error estimate; (e) the
# zeros and the prime side at 20 more digits than the printed result; (f) the gamma factor of xi'/xi tested by the limit at s = 1.
# Nothing here is a proof: the theta integrals are floating point (double precision accumulation of mp values at 20 digits).
# Output: data/dc_weil.json (numbers as strings, parameters, per part; a part written again replaces only its own entry) and the log
# beside it (same name, .log). Usage: python3 dc_weil.py [--quick] [--out FILE] [--zeros FILE] [--workers W] [--parts W1,W2,W3,W4] [--nzeros N]
# [--w3y 0.5,0.2] [--w4y 0.5,1]; --quick uses small sizes (about 3 s for W1 and W2, 16 s for W3 and 14 s for W4 in one
# process, so run the parts one at a time to stay under 30 s each) and writes to /tmp unless --out is given; the full sizes are meant for the Hopper cluster:
#   bash misc/tools/hopper_run.sh -c 16 -t 120 -g dc_weil.json -g dc_weil.log code/dc_weil.py data/zeros_1700.txt -- --workers 16 --out dc_weil.json
import argparse
import json
import math
import os
import sys
import tempfile
import time
from multiprocessing import get_context
from pathlib import Path

import mpmath as mp
import numpy as np

T0 = time.time()
LOGF = None

SIZES = {
    'quick': dict(dps=30, nzeros=300, circle=32,
                  W3=dict(Ts=20.0, T0w=40.0, T1=100.0, nodes=14, Nex=[100, 300]),
                  W4=dict(Ts=20.0, T0w=40.0, T1=150.0, nodes=14, Nex=[100, 300])),
    'full': dict(dps=50, nzeros=1700, circle=64,
                 W3=dict(Ts=40.0, T0w=200.0, T1=600.0, nodes=18, Nex=[100, 1000]),
                 W4=dict(Ts=40.0, T0w=500.0, T1=1500.0, nodes=18, Nex=[100, 1000])),
}
PRINTED_LAMBDA_A1 = {'0.5': 0.04684, '0.2': 0.01875}   # Table dc:tab:ladder at a = 1 (variational upper bounds)


def log(*a):
    line = '[%6.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(line, flush=True)
    if LOGF is not None:
        LOGF.write(line + '\n')
        LOGF.flush()


def ns(x, n=30):
    return mp.nstr(x, n)


def find_zeros(path):
    here = Path(__file__).resolve().parent
    for c in ([Path(path)] if path else []) + [Path('data/zeros_1700.txt'), Path('zeros_1700.txt'),
                                               here / 'zeros_1700.txt', here.parent / 'data' / 'zeros_1700.txt']:
        if c.is_file():
            return c
    sys.exit('zeros_1700.txt not found: pass --zeros FILE')


def load_zeros(path, count):
    out = []
    with mp.workdps(80):
        for line in open(path):
            t = line.split()
            if t:
                out.append(mp.mpf(t[-1]))
            if len(out) >= count:
                break
    return out


def prime_powers(N):
    """pairs (n, p) with n = p^k <= N, so that Lambda(n) = log p"""
    N = int(N)
    if N < 2:
        return []
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    out = []
    for p in np.nonzero(sieve)[0]:
        p = int(p)
        q = p
        while q <= N:
            out.append((q, p))
            q *= p
    out.sort()
    return out


def weil_terms(g, Lq, dps, poles=None, Nmax=None):
    """the terms of W(g) for an even g that vanishes (to working precision) outside [-Lq, Lq]; the pole term by quadrature unless given"""
    with mp.workdps(dps + 20):
        Lq = mp.mpf(Lq)
        g0 = g(mp.mpf(0))
        m = max(1, int(mp.ceil(Lq)))
        pts = [Lq * j / m for j in range(m + 1)]
        if poles is None:
            poles = 4 * mp.quad(lambda x: g(x) * mp.cosh(x / 2), pts)
        const = -(mp.euler + mp.log(mp.pi)) * g0
        arch = mp.quad(lambda x: 2 * (mp.exp(-2 * x) * g0 - mp.exp(-x / 2) * g(x)) / (-mp.expm1(-2 * x)), pts)
        tail = g0 * (-mp.log1p(-mp.exp(-2 * Lq)))
        if Nmax is None:
            Nmax = int(mp.floor(mp.exp(Lq)))
        pr = mp.mpf(0)
        for n, p in prime_powers(Nmax):
            pr += mp.log(p) / mp.sqrt(n) * g(mp.log(n))
        primes = -2 * pr
        total = poles + const + arch + tail + primes
        terms = dict(poles=poles, const=const, arch=arch, tail=tail, primes=primes, g0=g0)
    return terms, total


# ---------------------------------------------------------------- xi'/xi and its derivatives
def xi_logder(s):
    s = mp.mpmathify(s)
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s)


def xi_logder_1():
    """xi'/xi(1) = 1 + gamma_E/2 - (log 4 pi)/2, the limit of 1/(s-1) + zeta'/zeta(s) being gamma_E"""
    return 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2


def xi_logder_d1(s):
    L = mp.zeta(s, 1, 1) / mp.zeta(s)
    Lp = mp.zeta(s, 1, 2) / mp.zeta(s) - L ** 2
    return -1 / s ** 2 - 1 / (s - 1) ** 2 + mp.psi(1, s / 2) / 4 + Lp


def xi_logder_d2(s):
    z0, z1, z2, z3 = [mp.zeta(s, 1, k) for k in range(4)]
    L = z1 / z0
    Lp = z2 / z0 - L ** 2
    Lpp = z3 / z0 - (z2 / z0) * L - 2 * L * Lp
    return 2 / s ** 3 + 2 / (s - 1) ** 3 + mp.psi(2, s / 2) / 8 + Lpp


def xi_logder_d2_at_1():
    """(xi'/xi)''(1) with G(s) = 1/(s-1) + zeta'/zeta, G''(1) = 2 g0^3 + 6 g0 g1 + 3 g2 (Stieltjes constants)"""
    g0, g1, g2 = mp.stieltjes(0), mp.stieltjes(1), mp.stieltjes(2)
    return 2 + mp.psi(2, mp.mpf(1) / 2) / 8 + 2 * g0 ** 3 + 6 * g0 * g1 + 3 * g2


def cauchy_derivs(F, s0, r, N, nmax):
    """F^(n)(s0), n = 0..nmax, by the trapezoid rule on the circle |s - s0| = r (N points)"""
    vals = [F(s0 + r * mp.expjpi(mp.mpf(2 * j) / N)) for j in range(N)]
    out = []
    for n in range(nmax + 1):
        acc = sum(v * mp.expjpi(-mp.mpf(2 * j * n) / N) for j, v in enumerate(vals))
        out.append(mp.re(acc) / (N * mp.mpf(r) ** n) * mp.factorial(n))
    return out


def theta_prime(t):
    return mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * t / 2)) / 2 - mp.log(mp.pi) / 2


def p1_from_zeros(zeros, nz):
    """sum_{k<=nz} gamma_k^-2 and the tail int_{gamma_N}^inf t^-2 dNbar(t) - gamma_N^-2 S(gamma_N+), Nbar = theta/pi + 1"""
    gam = zeros[:nz]
    S = mp.fsum(1 / g ** 2 for g in gam)
    gN = gam[-1]
    jump = nz - (mp.siegeltheta(gN) / mp.pi + 1)
    smooth = mp.quad(lambda t: theta_prime(t) / (mp.pi * t ** 2), [gN, 2 * gN, 8 * gN, 64 * gN, 1024 * gN, mp.inf])
    return S, smooth - jump / gN ** 2


def m_of_y(ystr):
    s = mp.mpf('0.5') + mp.mpf(ystr)
    return 2 * (xi_logder_1() if s == 1 else xi_logder(s))


# ---------------------------------------------------------------- W1
def run_W1(zeros, dps):
    res = {'dps': dps, 'definition': 'g(x) = exp(-(x/c)^2); prime side = W(g) by the formula of the header, zero side = 2 sum_{gamma>0} ghat(gamma)'}
    for cs in ('0.5', '1'):
        with mp.workdps(dps + 20):
            c = mp.mpf(cs)
            g = lambda x: mp.exp(-(x / c) ** 2)
            L = c * mp.sqrt((dps + 20) * mp.log(10))
            X = max(mp.mpf(12), L)
            poles = 2 * c * mp.sqrt(mp.pi) * mp.exp(c ** 2 / 16)
            Nmax = int(mp.floor(mp.exp(L)))
        terms, W = weil_terms(g, X, dps, poles=poles, Nmax=Nmax)
        with mp.workdps(dps + 20):
            zero_side = 2 * mp.fsum(c * mp.sqrt(mp.pi) * mp.exp(-c ** 2 * gm ** 2 / 4) for gm in zeros)
            tail12 = -mp.log1p(-mp.exp(-24))
            diff = W - zero_side
            e = {'c': cs, 'prime_side_full': ns(W, 30), 'zero_side': ns(zero_side, 30), 'difference': ns(diff, 6),
                 'zeros_with_term_above_1e-60': sum(1 for gm in zeros if 2 * c * mp.sqrt(mp.pi) * mp.exp(-c ** 2 * gm ** 2 / 4) > mp.mpf(10) ** -60),
                 'prime_side_truncated_at_12_no_tail': ns(W - terms['g0'] * tail12, 20),
                 'omitted_tail_minus_log_1_minus_exp_minus_24': ns(tail12, 12),
                 'terms': {k: ns(v, 30) for k, v in terms.items()},
                 'X': ns(X, 6), 'prime_powers_up_to': Nmax, 'zeros_used': len(zeros)}
            res['c=' + cs] = e
        log('W1 c=%s: prime side %s, zero side %s, difference %s; truncated at 12 without tail %s'
            % (cs, ns(W, 18), ns(zero_side, 18), ns(diff, 4), e['prime_side_truncated_at_12_no_tail'][:12]))
    return res


# ---------------------------------------------------------------- W2
def run_W2(zeros, dps, nz, circ):
    out = {'dps': dps}
    with mp.workdps(dps + 10):
        C0 = 2 + mp.euler - mp.log(4 * mp.pi)
        out['two_plus_gammaE_minus_log4pi'] = ns(C0, 25)
        F1 = xi_logder_1()
        out['check_m_half_closed_form_minus_constant'] = ns(2 * F1 - C0, 5)
        eps = mp.mpf(10) ** -20
        out['check_limit_at_1_F(1+1e-20)-F1'] = ns(xi_logder(1 + eps) - F1, 5)
        ys = ['0.5', '0.375', '0.2', '0.1', '0.05']
        # V_xi three ways
        half = mp.mpf('0.5')
        v_exp = xi_logder_d1(half)
        v_circ = cauchy_derivs(xi_logder, half, mp.mpf('0.25'), circ, 1)[1]
        v_diff = mp.diff(xi_logder, half)
        out['V_xi'] = ns(v_exp, 30)
        out['V_xi_circle_minus_explicit'] = ns(v_circ - v_exp, 5)
        out['V_xi_mpdiff_minus_explicit'] = ns(v_diff - v_exp, 5)
        out['V_xi_over_2'] = ns(v_exp / 2, 30)
        out['F_at_half_should_be_zero'] = ns(xi_logder(half), 5)
        # m(y) and the margin
        mt = {}
        for ys_ in ys:
            m = m_of_y(ys_)
            y = mp.mpf(ys_)
            mt[ys_] = dict(m=ns(m, 25), m_over_y=ns(m / y, 12), two_Vxi_y=ns(2 * v_exp * y, 25),
                           two_Vxi_y_minus_m=ns(2 * v_exp * y - m, 10), printed_digits_m=ns(m, 4))
        out['m_of_y'] = mt
        out['two_Vxi'] = ns(2 * v_exp, 12)
        # second derivative of xi'/xi and the constant column of Table 1.2
        cons = {}
        for ys_ in ('0.5', '0.2', '0.1'):
            s = mp.mpf('0.5') + mp.mpf(ys_)
            d2c = cauchy_derivs(xi_logder, s, mp.mpf('0.25'), circ, 2)[2]
            d2e = xi_logder_d2_at_1() if s == 1 else xi_logder_d2(s)
            wpp = -d2c
            cons[ys_] = dict(Fpp_circle=ns(d2c, 25), Fpp_explicit=ns(d2e, 25), circle_minus_explicit=ns(d2c - d2e, 5),
                             w_pp_0=ns(wpp, 20), pi2_wpp_over_4=ns(mp.pi ** 2 * wpp / 4, 12))
        out['table_rate_constant_column'] = cons
        # p_1 from the zeros
        S1, T1_ = p1_from_zeros(zeros, nz)
        S2, T2_ = p1_from_zeros(zeros, nz // 2)
        out['p1_zeros_used'] = nz
        out['p1_partial_sum'] = ns(S1, 20)
        out['p1_tail'] = ns(T1_, 10)
        out['p1'] = ns(S1 + T1_, 15)
        out['p1_minus_Vxi_over_2'] = ns(S1 + T1_ - v_exp / 2, 5)
        out['p1_from_half_the_zeros_minus_Vxi_over_2'] = ns(S2 + T2_ - v_exp / 2, 5)
        out['p1_partial_sum_only_minus_Vxi_over_2'] = ns(S1 - v_exp / 2, 5)
    log('W2: 2+gE-log4pi = %s, V_xi = %s (circle - explicit %s, mpdiff - explicit %s), V_xi/2 = %s'
        % (out['two_plus_gammaE_minus_log4pi'][:12], out['V_xi'][:14], out['V_xi_circle_minus_explicit'],
           out['V_xi_mpdiff_minus_explicit'], out['V_xi_over_2'][:14]))
    for k, v in out['m_of_y'].items():
        log('W2: y = %s: m = %s, m/y = %s, 2 V_xi y = %s, 2 V_xi y - m = %s' % (k, v['m'][:12], v['m_over_y'][:10], v['two_Vxi_y'][:12], v['two_Vxi_y_minus_m']))
    for k, v in out['table_rate_constant_column'].items():
        log('W2: y = %s: (xi\'/xi)\'\' = %s (circle - explicit %s), pi^2 w\'\'(0)/4 = %s' % (k, v['Fpp_circle'][:14], v['circle_minus_explicit'], v['pi2_wpp_over_4'][:10]))
    log('W2: p_1 = %s from %d zeros (partial %s + tail %s), minus V_xi/2: %s; from %d zeros: %s'
        % (out['p1'], nz, out['p1_partial_sum'][:12], out['p1_tail'], out['p1_minus_Vxi_over_2'], nz // 2, out['p1_from_half_the_zeros_minus_Vxi_over_2']))
    return out


# ---------------------------------------------------------------- the theta integrals
def panels(a, b, width, n):
    npan = max(1, int(math.ceil((b - a) / width - 1e-9)))
    edges = np.linspace(a, b, npan + 1)
    x, w = np.polynomial.legendre.leggauss(n)
    mid = 0.5 * (edges[:-1] + edges[1:])
    rad = 0.5 * (edges[1:] - edges[:-1])
    return (mid[:, None] + rad[:, None] * x[None, :]).ravel(), (rad[:, None] * w[None, :]).ravel()


def _init(dps):
    mp.mp.dps = dps


def _eval_chunk(args):
    sig, ths = args
    sg = mp.mpf(sig)
    out = np.empty((len(ths), 2))
    lpi2 = mp.log(mp.pi) / 2
    for i, t in enumerate(ths):
        s = mp.mpc(sg, mp.mpf(float(t)))
        arch = 1 / s + 1 / (s - 1) - lpi2 + mp.digamma(s / 2) / 2
        zl = mp.zeta(s, 1, 1) / mp.zeta(s)
        out[i, 0] = float(arch.real)
        out[i, 1] = float(zl.real)
    return out


def eval_many(sig, ths, workers, dps, label, chunk=48):
    chunks = [(sig, ths[i:i + chunk]) for i in range(0, len(ths), chunk)]
    if workers <= 1 or len(chunks) < 2:
        _init(dps)
        res = [_eval_chunk(c) for c in chunks]
    else:
        res = []
        try:
            ctx = get_context('fork')
            with ctx.Pool(workers, initializer=_init, initargs=(dps,)) as pool:
                for k, r in enumerate(pool.imap(_eval_chunk, chunks, chunksize=1)):
                    res.append(r)
                    if (k + 1) % max(1, len(chunks) // 10) == 0:
                        log('  %s: %d of %d chunks' % (label, k + 1, len(chunks)))
        except Exception as exc:
            log('  %s: the pool failed (%r), continuing in one process' % (label, exc))
            _init(dps)
            res = [_eval_chunk(c) for c in chunks]
    return np.concatenate(res)


def taper(theta, T0w, T1, kind):
    u = np.clip((theta - T0w) / (T1 - T0w), 0.0, 1.0)
    s = u ** 4 * (35 - 84 * u + 70 * u ** 2 - 20 * u ** 3) if kind == 'smoother' else 0.5 * (1 - np.cos(np.pi * u))
    return 1.0 - s


_PP = {}


def ZN(theta, N, sigma):
    if N <= 0:
        return np.zeros_like(theta)
    if (N, sigma) not in _PP:
        pp = prime_powers(N)
        _PP[(N, sigma)] = (np.array([math.log(n) for n, p in pp]), np.array([math.log(p) * n ** (-sigma) for n, p in pp]))
    logs, coef = _PP[(N, sigma)]
    out = np.empty_like(theta)
    for i in range(0, len(theta), 4000):
        out[i:i + 4000] = -(np.cos(np.outer(theta[i:i + 4000], logs)) @ coef)
    return out


def LN_closed(N, sigma_str, G, dps):
    with mp.workdps(dps):
        sg = mp.mpf(sigma_str)
        s = mp.mpf(0)
        for n, p in prime_powers(N):
            s += mp.log(p) * mp.power(n, -sg) * G(mp.log(n))
        return -s / mp.pi


def theta_side(label, y, W, Sz, workers, dps, headline_hybrid):
    """(1/pi) int_R w_y(theta) Phi(theta) dtheta, sigma = 1/2 + y; W holds the weight Phi, Phi_s, Phi_o (numpy), Phi_s (mp) and G (mp)"""
    Ts, T0w, T1, nn = Sz['Ts'], Sz['T0w'], Sz['T1'], Sz['nodes']
    with mp.workdps(30):
        sigma_str = ns(mp.mpf('0.5') + mp.mpf(y), 20)
    sigma = float(sigma_str)
    hfar = 0.8 * float(y)
    far = 2 * hfar
    near = min(0.5, far)

    def region0(ref):
        a1, w1 = panels(0.0, 4.0, near / ref, nn)
        a2, w2 = panels(4.0, Ts, far / ref, nn)
        return np.concatenate([a1, a2]), np.concatenate([w1, w2])

    th0, wt0 = region0(1)
    th0r, wt0r = region0(2)
    th1, wt1 = panels(Ts, T1, far, nn)
    log('%s: sigma = %s, panels half width %.3f, %d + %d nodes (T_s %g, window %g to %g)' % (label, sigma_str, hfar, len(th0), len(th1), Ts, T0w, T1))
    AZ0 = eval_many(sigma_str, th0, workers, dps, label + ' region 0')
    AZ0r = eval_many(sigma_str, th0r, workers, dps, label + ' region 0 refined')
    AZ1 = eval_many(sigma_str, th1, workers, dps, label + ' region 1')
    phi0, phi0r = W['phi'](th0), W['phi'](th0r)
    D = dict(th0=th0, wt0=wt0, AZ0=AZ0, th1=th1, wt1=wt1, AZ1=AZ1, phi0=phi0, sigma=sigma,
             phis1=W['phis'](th1), phio1=W['phio'](th1))
    base0 = math.fsum(wt0 * (AZ0[:, 0] + AZ0[:, 1]) * phi0)
    ref0 = math.fsum(wt0r * (AZ0r[:, 0] + AZ0r[:, 1]) * phi0r)
    # the smooth tail of the gamma factor part
    with mp.workdps(20):
        sg = mp.mpf(sigma_str)

        def Aint(t):
            s = mp.mpc(sg, t)
            return mp.re(1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2) * W['phis_mp'](t)
        pts = [mp.mpf(Ts) * f for f in (1, 2, 8, 32, 256, 4096, 65536)] + [mp.inf]
        I2, e2 = mp.quad(Aint, pts, error=True)
        I2 = float(I2)
        e2 = float(e2)

    def assemble(N, win, kind):
        zc0 = AZ0[:, 1] - ZN(th0, N, sigma)
        zc1 = AZ1[:, 1] - ZN(th1, N, sigma)
        chi = taper(th1, win[0], win[1], kind)
        I1 = math.fsum(wt0 * AZ0[:, 0] * phi0)
        I4a = math.fsum(wt0 * zc0 * phi0)
        I3 = math.fsum(wt1 * chi * AZ1[:, 0] * D['phio1'])
        I4b = math.fsum(wt1 * chi * zc1 * (D['phis1'] + D['phio1']))
        return 2 / math.pi * (I1 + I3 + I4a + I4b + I2)

    wins = [(T0w, T1), (T0w, 0.75 * T1), (max(Ts, 0.5 * T0w), T1)]
    LNc = {0: 0.0}
    for N in Sz['Nex']:
        LNc[N] = float(LN_closed(N, sigma_str, W['G'], 30))
    vals = []
    for N in [0] + Sz['Nex']:
        for k, win in enumerate(wins):
            vals.append(dict(N=N, window=[win[0], win[1]], kind='smoother', L=assemble(N, win, 'smoother') + LNc[N]))
        vals.append(dict(N=N, window=[wins[0][0], wins[0][1]], kind='cosine', L=assemble(N, wins[0], 'cosine') + LNc[N]))
    Nh = Sz['Nex'][0] if headline_hybrid else 0
    head = [v['L'] for v in vals if v['N'] == Nh and v['kind'] == 'smoother']
    L_head = [v['L'] for v in vals if v['N'] == Nh and v['kind'] == 'smoother'][0]
    res = dict(sigma=sigma_str, T_s=Ts, T0=T0w, T1=T1, nodes=nn, panel_half_width=hfar, headline_N=Nh,
               headline=repr(L_head), headline_spread_over_windows=repr(max(head) - min(head)),
               panel_check_region0_refined_minus_base=repr(ref0 - base0), smooth_tail_integral=repr(I2),
               smooth_tail_quad_error_estimate=repr(e2), values=[dict(N=v['N'], window=v['window'], kind=v['kind'], L=repr(v['L'])) for v in vals])
    log('%s: panel check %.2e, smooth tail %.10g (quad error %.1e), hybrid N=%d L = %.12f, spread over windows %.2e'
        % (label, ref0 - base0, I2, e2, Nh, L_head, max(head) - min(head)))
    for v in vals:
        if v['kind'] == 'smoother' and v['window'] == [wins[0][0], wins[0][1]]:
            log('    N = %4d, window [%g, %g]: L = %.13f' % (v['N'], v['window'][0], v['window'][1], v['L']))
    return res, L_head


# ---------------------------------------------------------------- W3
def hh_conv_mp(a):
    k = mp.pi / (2 * a)

    def f(x):
        x = abs(x)
        if x >= 2 * a:
            return mp.mpf(0)
        return (2 * a - x) * mp.cos(k * x) / 2 + mp.sin(k * x) / (2 * k)
    return f


def run_W3(Sz, ys, workers, dps, dps_eval):
    res = {'a': 1, 'h': 'cos(pi u/(2a)) on [-a, a]', 'size': Sz,
           'note_y_0.2': 'at y = 0.2 the identity presupposes a zero free half plane Re s > 0.7, so this checks the normalisation, not the hypothesis'}
    a = 1.0
    k = math.pi / (2 * a)
    with mp.workdps(dps + 20):
        aa = mp.mpf(1)
        kk = mp.pi / (2 * aa)
        conv = hh_conv_mp(aa)
        chk = []
        for th in ('0.7', '3.3', '25.1'):
            t = mp.mpf(th)
            num = mp.quad(lambda u: mp.cos(kk * u) * mp.cos(t * u), [-aa, 0, aa])
            chk.append(abs(num - 2 * kk * aa * mp.sinc(aa * (kk - t)) / (kk + t)))
        for x in ('0.3', '1.2', '1.9'):
            xx = mp.mpf(x)
            num = mp.quad(lambda u: mp.cos(kk * u) * mp.cos(kk * (xx - u)), [max(-aa, xx - aa), min(aa, xx + aa)])
            chk.append(abs(num - conv(xx)))
        res['closed_form_checks_max_abs_error'] = ns(max(chk), 5)
    log('W3: closed form checks of hhat and h*h against quadrature, max error', res['closed_form_checks_max_abs_error'])

    def phi(t):
        return (2 * k * a * np.sinc(a * (k - t) / np.pi) / (k + t)) ** 2

    def phis(t):
        return 2 * k ** 2 / (k ** 2 - t ** 2) ** 2

    def phio(t):
        return 2 * k ** 2 * np.cos(2 * a * t) / (k ** 2 - t ** 2) ** 2

    def phis_mp(t):
        kk_ = mp.pi / 2
        return 2 * kk_ ** 2 / (kk_ ** 2 - t ** 2) ** 2

    def G(u):
        return 2 * mp.pi * conv(u)

    W = dict(phi=phi, phis=phis, phio=phio, phis_mp=phis_mp, G=G)
    for ystr in ys:
        with mp.workdps(dps + 20):
            y = mp.mpf(ystr)
            g = lambda x: mp.exp(-y * abs(x)) * conv(x)
        terms, Q = weil_terms(g, 2, dps)
        with mp.workdps(dps + 20):
            m = m_of_y(ystr)
            qa = Q / aa
            e = dict(y=ystr, Q_prime_side=ns(Q, 25), terms={kk_: ns(v, 25) for kk_, v in terms.items()},
                     rayleigh_Q_over_norm2=ns(qa, 15), m_y=ns(m, 15), rayleigh_minus_m=ns(qa - m, 12),
                     a2_times_excess=ns(qa - m, 12), printed_lambda_1_y=PRINTED_LAMBDA_A1.get(ystr),
                     rayleigh_minus_printed_lambda=ns(qa - mp.mpf(str(PRINTED_LAMBDA_A1[ystr])), 6) if ystr in PRINTED_LAMBDA_A1 else None)
        log('W3 y=%s: prime side Q_y(h) = %s, Q/||h||^2 = %s, m(y) = %s, Q/||h||^2 - m = %s' % (ystr, ns(Q, 18), ns(qa, 12), ns(m, 12), ns(qa - m, 8)))
        t, L = theta_side('W3 y=' + ystr, ystr, W, Sz, workers, dps_eval, headline_hybrid=Sz['headline_hybrid'])
        e['theta_side'] = t
        e['theta_side_minus_prime_side'] = repr(L - float(Q))
        e['pure_brute_force_minus_prime_side'] = repr(float(t['values'][0]['L']) - float(Q))
        log('W3 y=%s: theta side %.13f minus prime side %.13f = %.3e; pure brute force (N = 0, first window) minus prime side %s'
            % (ystr, L, float(Q), L - float(Q), e['pure_brute_force_minus_prime_side']))
        res['y=' + ystr] = e
    return res


# ---------------------------------------------------------------- W4
def zero_sum_triangle(zeros, nz):
    """2 sum_{gamma>0} phi(gamma), phi = sin^2(t)/(pi t^2), with the tail by Nbar = theta/pi + 1 and two terms of the expansion of the oscillating integral"""
    gam = zeros[:nz]
    S = mp.fsum(mp.sin(g) ** 2 / g ** 2 for g in gam)
    gN = gam[-1]
    jump = nz - (mp.siegeltheta(gN) / mp.pi + 1)
    fo = lambda t: theta_prime(t) / (2 * mp.pi * t ** 2)
    smooth = mp.quad(fo, [gN, 2 * gN, 8 * gN, 64 * gN, 1024 * gN, mp.inf])
    osc = -(-fo(gN) * mp.sin(2 * gN) / 2 - mp.diff(fo, gN) * mp.cos(2 * gN) / 4)
    tail = smooth + osc - mp.sin(gN) ** 2 / gN ** 2 * jump
    return 2 / mp.pi * S, 2 / mp.pi * tail


def run_W4(Sz, y0s, zeros, nz, workers, dps, dps_eval):
    res = {'phihat': 'triangle (1 - |x|/2) on [-2, 2]', 'g': 'phihat/(2 pi)', 'phi': 'sin^2(theta)/(pi theta^2)', 'size': Sz}

    def g(x):
        x = abs(x)
        return (1 - x / 2) / (2 * mp.pi) if x < 2 else mp.mpf(0)
    terms, Wg = weil_terms(g, 2, dps)
    res['W_g'] = ns(Wg, 25)
    res['W_g_terms'] = {k: ns(v, 25) for k, v in terms.items()}
    log('W4: W(g) = %s (poles %s, constant %s, arch %s + tail %s, primes %s)' % (ns(Wg, 18), ns(terms['poles'], 10), ns(terms['const'], 10), ns(terms['arch'], 10), ns(terms['tail'], 6), ns(terms['primes'], 10)))
    with mp.workdps(dps):
        zs_all, zs_tail = zero_sum_triangle(zeros, nz)
        zs_h_all, zs_h_tail = zero_sum_triangle(zeros, nz // 2)
        res['zero_sum'] = dict(zeros_used=nz, partial=ns(zs_all, 15), tail=ns(zs_tail, 6), value=ns(zs_all + zs_tail, 15),
                               minus_W_g=ns(zs_all + zs_tail - Wg, 4), partial_only_minus_W_g=ns(zs_all - Wg, 4),
                               value_from_half_the_zeros_minus_W_g=ns(zs_h_all + zs_h_tail - Wg, 4))
    log('W4: zero sum 2 sum phi(gamma) over %d zeros: partial %s + tail %s = %s, minus W(g) %s (partial only %s; from %d zeros %s)'
        % (nz, res['zero_sum']['partial'][:12], res['zero_sum']['tail'], res['zero_sum']['value'][:14], res['zero_sum']['minus_W_g'], res['zero_sum']['partial_only_minus_W_g'], nz // 2, res['zero_sum']['value_from_half_the_zeros_minus_W_g']))
    out_vals = {}
    for ystr in y0s:
        y0 = float(ystr)
        with mp.workdps(30):
            y0m = mp.mpf(ystr)
            chk = []
            for th in ('0.7', '3.3', '25.1'):
                t = mp.mpf(th)
                num = mp.quad(lambda x: (1 - x / 2) * mp.exp(y0m * x) * mp.cos(t * x), [0, 2]) / mp.pi
                z = mp.mpc(y0m, t)
                chk.append(abs(num - mp.re((mp.exp(2 * z) - 1) / (2 * z ** 2) - 1 / z) / mp.pi))
            cf = ns(max(chk), 5)

        def phi(t, y0=y0):
            z = y0 + 1j * t
            return (((np.exp(2 * z) - 1) / (2 * z ** 2) - 1 / z)).real / math.pi

        def phis(t, y0=y0):
            z = y0 + 1j * t
            return (-(1 / (2 * z ** 2) + 1 / z)).real / math.pi

        def phio(t, y0=y0):
            z = y0 + 1j * t
            return (np.exp(2 * z) / (2 * z ** 2)).real / math.pi

        def phis_mp(t, y0m=None):
            z = mp.mpc(mp.mpf(ystr), t)
            return mp.re(-(1 / (2 * z ** 2) + 1 / z)) / mp.pi

        def G(u, ystr=ystr):
            return (1 - u / 2) * mp.exp(mp.mpf(ystr) * u) if 0 <= u <= 2 else mp.mpf(0)

        W = dict(phi=phi, phis=phis, phio=phio, phis_mp=phis_mp, G=G)
        log('W4 y0=%s: closed form of P^{-1} phi against quadrature, max error %s' % (ystr, cf))
        t, L = theta_side('W4 y0=' + ystr, ystr, W, Sz, workers, dps_eval, headline_hybrid=Sz['headline_hybrid'])
        pure = float(t['values'][0]['L'])
        e = dict(y0=ystr, closed_form_check=cf, theta_side=t, theta_side_minus_W_g=repr(L - float(Wg)),
                 pure_brute_force_minus_W_g=repr(pure - float(Wg)), theta_side_minus_zero_sum=repr(L - float(zs_all + zs_tail)))
        log('W4 y0=%s: D = %.13f, W(g) = %.13f, difference %.3e; pure brute force (N = 0, first window) minus W(g) %.3e'
            % (ystr, L, float(Wg), L - float(Wg), pure - float(Wg)))
        out_vals[ystr] = L
        res['y0=' + ystr] = e
    if len(out_vals) == 2:
        ks = sorted(out_vals)
        res['D_at_first_minus_D_at_second'] = repr(out_vals[ks[0]] - out_vals[ks[1]])
        log('W4: D(y0=%s) - D(y0=%s) = %.3e' % (ks[0], ks[1], out_vals[ks[0]] - out_vals[ks[1]]))
    return res


# ---------------------------------------------------------------- main
def save_json(path, part, payload, meta):
    data = {}
    if os.path.exists(path):
        try:
            data = json.load(open(path))
        except Exception:
            data = {}
    data.setdefault('meta', {}).update(meta)
    data[part] = payload
    with open(path, 'w') as f:
        json.dump(data, f, indent=1)


def main():
    global LOGF
    ap = argparse.ArgumentParser(description='Weil side checks of the deconvolution chapter (group D1)')
    ap.add_argument('--quick', action='store_true', help='small sizes, local, writes to /tmp unless --out is given')
    ap.add_argument('--out')
    ap.add_argument('--zeros')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--parts', default='W1,W2,W3,W4')
    ap.add_argument('--nzeros', type=int, help='number of zeros used (default 300 with --quick, else 1700)')
    ap.add_argument('--w3y', help='comma separated list of y for W3 (default 0.5,0.2; with --quick 0.5)')
    ap.add_argument('--w4y', help='comma separated list of y0 for W4 (default 0.5,1; with --quick 1)')
    args = ap.parse_args()
    mode = 'quick' if args.quick else 'full'
    S = SIZES[mode]
    if args.out:
        out = Path(args.out)
    elif args.quick:
        out = Path(tempfile.gettempdir()) / 'dc_weil_quick.json'
    else:
        d = Path(__file__).resolve().parent.parent / 'data'
        out = (d if d.is_dir() else Path('.')) / 'dc_weil.json'
    LOGF = open(out.with_suffix('.log'), 'a')
    log('dc_weil.py', ' '.join(sys.argv[1:]), '| mode', mode, '| mpmath', mp.__version__, '| numpy', np.__version__, '| out', out)
    parts = [p.strip() for p in args.parts.split(',') if p.strip()]
    nz = args.nzeros or S['nzeros']
    zeros = None
    if any(p in parts for p in ('W1', 'W2', 'W4')):
        zpath = find_zeros(args.zeros)
        zeros = load_zeros(zpath, nz)
        log('zeros: %d from %s, first %s, last %s' % (len(zeros), zpath, ns(zeros[0], 20), ns(zeros[-1], 20)))
    meta = {'script': 'dc_weil.py', 'last_run_mode': mode}
    dps = S['dps']
    if 'W1' in parts:
        save_json(out, 'W1', dict(run_W1(zeros, dps), size=mode), meta)
    if 'W2' in parts:
        save_json(out, 'W2', dict(run_W2(zeros, dps, nz, S['circle']), size=mode), meta)
    if 'W3' in parts:
        ys = (args.w3y or ('0.5' if args.quick else '0.5,0.2')).split(',')
        Sz = dict(S['W3'], headline_hybrid=args.quick)
        save_json(out, 'W3', dict(run_W3(Sz, ys, args.workers, dps, 20), size=mode), meta)
    if 'W4' in parts:
        ys = (args.w4y or ('1' if args.quick else '0.5,1')).split(',')
        Sz = dict(S['W4'], headline_hybrid=args.quick)
        save_json(out, 'W4', dict(run_W4(Sz, ys, zeros, nz, args.workers, dps, 20), size=mode), meta)
    log('done; written', out)


if __name__ == '__main__':
    main()
