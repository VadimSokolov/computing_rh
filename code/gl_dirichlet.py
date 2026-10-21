# gl_dirichlet.py: the Dirichlet L function claims of Chapter 21 (Thorin Representations, Dirichlet L Functions and the Linnik
# Constant, label prefix gl:), group L2a of the check of that chapter: lines C:541 to C:584 (complex characters, the Dedekind route),
# C:586 to C:620 (the Levy decomposition of xi_F with its signs) and C:701 to C:731 (numerical verification for the characters
# modulo 4 and 5). Added for the book in October 2026, written from the mathematics (the chapter author's scripts were not read).
# Conventions. chi_4 odd (a = 1), chi_5 the Legendre symbol modulo 5 (even, a = 0), psi the complex character modulo 5 with
# psi(2) = i (odd, since psi(4) = -1). Lambda(s, chi) = (q/pi)^((s+a)/2) Gamma((s+a)/2) L(s, chi) = xi(s, chi), no factor s(s-1) for
# a nonprincipal character (this is the chapter's normalisation, C:513 to C:540); xi(s) is Riemann's xi, xi_F is the chapter's
# eq. gl:eq:xiF. All characters are explicit tables of exact Gaussian integers; L(s, chi) = q^(-s) sum_b chi(b) zeta(s, b/q) (Hurwitz).
# Parts (--parts, default all).
#  zeros     C:705. Certified count of the zeros of L(s, chi) with 0 < Im s < T by the argument principle with interval arithmetic
#            (python-flint balls): N = (theta(T) + Delta arg L along 1/2 -> 2 -> 2+iT -> 1/2+iT)/pi, the half contour of the rectangle
#            [-1, 2] x [0, T] (the functional equation doubles it). Delta arg L is a sum of principal arguments over boxes on each of
#            which the ball enclosure of L excludes 0 and lies in a disc of radius half its centre modulus (so the argument is
#            continuous); on Re s = 2 it is the principal argument of L(2+iT)/L(2) because |L - 1| <= zeta(2) - 1 < 1; near s = 1 a
#            Cauchy estimate on the circle |s - 1| = 0.3 replaces the Hurwitz enclosure. The certificate is an integer inside an
#            enclosure of width below 1e-6. The same T is scanned for sign changes of the real function xi(1/2 + it, chi) on a
#            grid (certified signs); if their number equals N every zero is on the line and simple. Ordinates are refined by
#            Illinois iteration at 110 bits and enclosed by two certified signs. The main term is compared in the chapter's form
#            (T/2pi) log(qT/2pi e) and in the exact form theta(T)/pi (which differs by chi(-1)/8 with a sign), S(T) = N - theta/pi.
#  real      C:707. The grid of step 0.005, L(1/2, chi), and a certificate (ball enclosures, bisection) that L(sigma, chi) > 0 on
#            [1/2, 2]; with xi(1 - s) = xi(s) this gives xi > 0 on [0, 1], so no real zero in the critical strip and L > 0 on (0, 1).
#  classical C:707, C:729. (xi'/xi)(1, chi) by the formula (gl:eq:chi-classical) and by the series sum (-1)^(j+1) p_j / 4^(j-1) from
#            the Taylor coefficients at 1/2, xi(1/2, chi_4), xi(1, chi_4), E exp(-H/4) = xi(1/2)/xi(1) against the zero product,
#            both sides of the sentence "p_1 - (xi'/xi)(1) = p_2/4 - p_3/16 + ...".
#  table     C:710 to C:721. p_j (j = 1, 2, 3) from the Taylor coefficients of log xi(1/2 + s, chi) (Cauchy integral on a circle
#            of radius 2, 96 points, 60 digits, radii 2 and 3 compared (half integer radii hit the poles of the gamma factor, which
#            cancel against the trivial zeros of L); p_j = (-1)^(j+1) j b_(2j)) against sum_(gamma < T) gamma^(-2j)
#            plus the tail integral from T; relative errors; the residual times T^(2j) against S(T) (the scaling sentence C:721).
#  transfer  C:723. T_chi(t) = xi(1/2 + it)/xi(1/2 + t) at t = 2, 5, the Hadamard product at s = 3 (with and without tail).
#  negmom    C:727. E H_chi^(-1/2) = (2 xi(1/2, chi)/sqrt(pi)) int_0^inf dw/xi(1/2 + w, chi), checked on zeta (6.733272312756635)
#            and against the zero clock (2/sqrt(pi)) int_0^inf prod_gamma (1 + w^2/gamma^2)^(-1) dw.
#  complex   C:541, C:556. xi(1/2 + s, psi)/xi(1/2 - s, psi) and the product xi_(psi psibar) at s = 0.7, with the root number.
#  dedekind  C:584. xi_F/(xi xi(., chi)) at s = 0.3, 2, 3.7 for Q(i), Q(sqrt 5), zeta_F built independently: Euler product over
#            the prime ideals (splitting by Dedekind Kummer: roots of x^2 + 1 and x^2 - x - 1 modulo p, by the Euler criterion for
#            large p) with the prime number theorem tail for s = 2, 3.7, and the approximate functional equation of zeta_F
#            (Mellin Barnes with the gamma factor of F, residue from the class number formula) for s = 0.3, validated at s = 2 and
#            3.7 against the Euler product and on zeta itself; additivity of E H_F for Q(i) from xi_F built by that equation.
#  levy      C:615. The two sides of gl:eq:levyF at (alpha, s) = (3, 1) for Q, (3, 3/2) for Q(i), (5/2, 1) for Q(sqrt 5): the log ratio
#            of xi_F from the L function values against the gamma integrals plus the prime ideal sum up to N (Lambda_F(n) from
#            the splitting), with and without the density one tail correction.
# Error control. Hurwitz L against flint's dirichlet_l (Conrey characters) at test points; AFE against zeta; reality of the
# rotated functions is a rounding check only (it follows from the functional equation).
# Output: --out (default gl_dirichlet.json in the current folder, /tmp/gl_dirichlet_quick.json with --quick) with the log beside it
# (same stem, .log); numbers are strings. A partial run (--parts) merges into an existing output file.
# Usage: python3 gl_dirichlet.py [--quick] [--out FILE] [--parts zeros,real,...] [--T 300] [--h 0.03125] [--P 10000000] [--N 2000000]
# --quick runs T = 50 and small sums in about 13 s locally. The full case (T = 300, prime cutoff 2e6, ideal cutoff 1e7) took 15.8 s for
# --parts zeros,real,classical and 11.5 s for --parts table,transfer,negmom,complex,dedekind,levy on a laptop (one process each); as one job
# on the Hopper cluster:
#   bash misc/tools/hopper_run.sh -c 1 -t 90 -g gl_dirichlet.json -g gl_dirichlet.log code/gl_dirichlet.py -- --out gl_dirichlet.json
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext

import mpmath as mp
import numpy as np
from flint import acb, arb, ctx, dirichlet_char

LOGF = None


def log(*a):
    msg = ' '.join(str(x) for x in a)
    print(msg, flush=True)
    if LOGF is not None:
        LOGF.write(msg + '\n')
        LOGF.flush()


class prec:
    """Context manager for the python-flint working precision (bits)."""

    def __init__(self, p):
        self.p = p

    def __enter__(self):
        self.old = ctx.prec
        ctx.prec = self.p

    def __exit__(self, *e):
        ctx.prec = self.old


# ---------------------------------------------------------------------------------------------------------------- characters
class Char:
    def __init__(self, name, q, tab, conrey):
        self.name = name
        self.q = q
        self.tab = [(b, v[0], v[1]) for b, v in sorted(tab.items())]
        m1 = dict(tab)[q - 1]
        self.a = 0 if m1 == (1, 0) else 1
        self.conrey = conrey
        self.real = all(vi == 0 for _, _, vi in self.tab)


CH = {
    'chi4': Char('chi4', 4, {1: (1, 0), 3: (-1, 0)}, (4, 3)),
    'chi5': Char('chi5', 5, {1: (1, 0), 2: (-1, 0), 3: (-1, 0), 4: (1, 0)}, (5, 4)),
    'psi': Char('psi', 5, {1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}, (5, 2)),
    'psibar': Char('psibar', 5, {1: (1, 0), 2: (0, -1), 3: (0, 1), 4: (-1, 0)}, (5, 3)),
}


def L_fl(ch, s):
    tot = acb(0)
    for b, vr, vi in ch.tab:
        z = s.zeta(acb(b) / ch.q)
        tot = tot + (z * vr if vi == 0 else z * acb(vr, vi))
    return tot * (acb(ch.q) ** (-s))


def G_fl(ch, s):
    w = (s + ch.a) / 2
    return (acb(ch.q) / arb.pi()) ** w * w.gamma()


def logG_fl(ch, s):
    w = (s + ch.a) / 2
    return w * (acb(ch.q) / arb.pi()).log() + w.lgamma()


def L_mp(ch, s):
    tot = mp.mpc(0)
    for b, vr, vi in ch.tab:
        tot += mp.mpc(vr, vi) * mp.zeta(s, mp.mpf(b) / ch.q)
    return tot * mp.mpf(ch.q) ** (-s)


def G_mp(ch, s):
    w = (s + ch.a) / 2
    return (mp.mpf(ch.q) / mp.pi) ** w * mp.gamma(w)


def xi_mp(ch, s):
    return G_mp(ch, s) * L_mp(ch, s)


def xiR_mp(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def a2mp(x, digits=34):
    return mp.mpf(x.mid().str(digits, radius=False))


def nst(x, d=25):
    return mp.nstr(x, d)


def agree(claimed, mine):
    """Compare a printed value (string) with mine (mpf): returns 'agrees', 'agrees within 1 unit of the last digit' or 'differs'."""
    getcontext().prec = 80
    c = Decimal(claimed)
    e = c.as_tuple().exponent
    q = Decimal(1).scaleb(e)
    m = Decimal(mp.nstr(mine, 50))
    if m.quantize(q) == c:
        return 'agrees'
    if abs(m - c) < q:
        return 'agrees within 1 unit of the last digit'
    return 'differs'


def claim(label, claimed, mine, digits=12):
    st = agree(claimed, mine)
    log('  claim %-46s chapter %-16s mine %-22s %s' % (label, claimed, nst(mine, digits), st))
    return {'chapter': claimed, 'mine': nst(mine, 20), 'status': st}


# ----------------------------------------------------------------------------------------------- near s = 1 enclosure of L
class LEnc:
    """Rigorous enclosures of L(s, chi) on balls, with the pole of the Hurwitz terms at s = 1 removed by a Cauchy estimate."""

    R0 = 0.08
    RC = 0.3

    def __init__(self, ch):
        self.ch = ch
        self.L1 = acb(1).dirichlet_l(dirichlet_char(*ch.conrey))
        M = 0.0
        n = 128
        for k in range(n):
            th = 2 * math.pi * k / n
            E = L_fl(ch, acb(arb(1 + self.RC * math.cos(th), 0.02), arb(self.RC * math.sin(th), 0.02)))
            M = max(M, float(abs(E).upper()))
        self.M = M * (1 + 1e-9)
        self.C = self.M * self.RC / (self.RC - self.R0) ** 2

    def __call__(self, s):
        m = s.mid()
        dx = abs(float(m.real) - 1.0)
        dy = abs(float(m.imag))
        rr = float(s.real.rad())
        ri = float(s.imag.rad())
        if rr == 0 and ri == 0:
            if dx == 0 and dy == 0:
                return self.L1
            return L_fl(self.ch, s)
        dmax = math.hypot(dx + rr, dy + ri)
        if dmax <= self.R0:
            # Lipschitz bound |L'| <= C on |z - 1| <= R0 (Cauchy on |z - 1| = RC): L(z) lies within C |z - m| of L(m), m the centre
            r = math.hypot(rr, ri) * (1 + 1e-9)
            Lm = self.L1 if (dx == 0 and dy == 0) else L_fl(self.ch, m)
            return Lm + acb(arb(0, self.C * r), arb(0, self.C * r))
        if dx <= rr and dy <= ri:
            return None
        return L_fl(self.ch, s)


def good(E):
    if E is None:
        return False
    lo = abs(E).lower()
    if not lo > 0:
        return False
    rho = float(E.rad())
    mid = float(abs(E.mid()))
    return rho < 0.5 * mid


def chain_arg(f, a, b, maxdepth=34):
    """Certified change of arg f along the segment a -> b (python complex numbers with exact binary parts)."""
    cache = {}

    def fp(z):
        key = (z.real, z.imag)
        v = cache.get(key)
        if v is None:
            v = f(acb(z.real, z.imag))
            cache[key] = v
        return v

    def rec(a, b, depth):
        m = (a + b) / 2
        hx = abs(b.real - a.real) / 2
        hy = abs(b.imag - a.imag) / 2
        E = f(acb(arb(m.real, hx), arb(m.imag, hy)))
        if good(E):
            return (fp(b) / fp(a)).arg()
        if depth >= maxdepth:
            raise RuntimeError('chain_arg: no certified box near %r, %r' % (a, b))
        return rec(a, m, depth + 1) + rec(m, b, depth + 1)

    return rec(a, b, 0)


def sgn(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def count_positive(ch, T, lenc):
    """Certified number of zeros of L(s, ch) with 0 < Im s < T (all in 0 < Re s < 1): returns (N, enclosure of N, theta)."""
    theta = logG_fl(ch, acb(0.5, T)).imag
    P1 = chain_arg(lenc, complex(0.5, 0.0), complex(2.0, 0.0))
    P2 = (L_fl(ch, acb(2.0, T)) / L_fl(ch, acb(2.0, 0.0))).arg()
    P3 = chain_arg(lenc, complex(2.0, T), complex(0.5, T))
    tot = (theta + P1 + P2 + P3) / arb.pi()
    n = int(round(float(tot)))
    if not (tot - n).abs_upper() < 1e-3:
        raise RuntimeError('count_positive: enclosure %s is not within 1e-3 of an integer' % tot.str(12))
    return n, tot, theta


def certify_positive_real(lenc, lo, hi, pieces, depth=0):
    mid = (lo + hi) / 2
    rad = (hi - lo) / 2
    E = lenc(acb(arb(mid, rad), arb(0)))
    if E is not None and E.real > 0:
        pieces.append(float(E.real.lower()))
        return True
    if depth >= 40:
        return False
    return certify_positive_real(lenc, lo, mid, pieces, depth + 1) and certify_positive_real(lenc, mid, hi, pieces, depth + 1)


def make_fx(ch, rot=None):
    def fx(t):
        s = acb(arb(0.5), t)
        v = G_fl(ch, s) * L_fl(ch, s)
        if rot is not None:
            v = v * rot
        return v.real

    return fx


def scan_brackets(ch, T, h, rot=None):
    n = int(round(T / h))
    fx = make_fx(ch, rot)
    out = []
    pt = 0.0
    ps = sgn(fx(arb(0.0)))
    for k in range(1, n + 1):
        t = k * h
        sg = sgn(fx(arb(t)))
        if sg == 0:
            for p in (128, 256):
                with prec(p):
                    sg = sgn(fx(arb(t)))
                if sg != 0:
                    break
        if sg != ps:
            out.append((pt, t))
        pt, ps = t, sg
    return out


def refine_zero(fx, a, b, tol=1e-22, maxit=80):
    a, b = arb(a), arb(b)
    fa, fb = fx(a), fx(b)
    sa = sgn(fa)
    side = 0
    for _ in range(maxit):
        c = ((a * fb - b * fa) / (fb - fa)).mid()
        fc = fx(c)
        sc = sgn(fc)
        if sc == 0:
            a = b = c
            break
        if sc == sa:
            a, fa = c, fc
            if side == -1:
                fb = fb / 2
            side = -1
        else:
            b, fb = c, fc
            if side == 1:
                fa = fa / 2
            side = 1
        if abs(float(b - a)) < tol:
            break
    c = ((a + b) / 2).mid()
    d = max(abs(float(b - a)), 1e-25)
    for _ in range(40):
        lo, hi = c - d, c + d
        s1, s2 = sgn(fx(lo)), sgn(fx(hi))
        if s1 != 0 and s2 != 0 and s1 != s2:
            return c, lo, hi
        d *= 4
    raise RuntimeError('refine_zero: no certified bracket near %s' % c.str(20))


# ---------------------------------------------------------------------------------------------------------- series helpers
def taylor_coeffs(fun, r, N, kmax):
    om = [mp.expjpi(mp.mpf(2 * m) / N) for m in range(N)]
    vals = [fun(r * w) for w in om]
    c = []
    for k in range(kmax + 1):
        acc = mp.mpc(0)
        for m in range(N):
            acc += vals[m] * om[(-m * k) % N]
        c.append(acc / (N * mp.mpf(r) ** k))
    return c


def log_coeffs(c):
    K = len(c) - 1
    b = [mp.mpc(0)] * (K + 1)
    for k in range(1, K + 1):
        acc = c[k]
        for m in range(1, k):
            acc -= mp.mpf(m) / k * b[m] * c[k - m]
        b[k] = acc / c[0]
    return b


def pj_function(fun, jmax, r, N):
    """p_j = (-1)^(j+1) j b_(2j) for fun(s) = xi(1/2 + s), by a Cauchy integral on |s| = r."""
    c = taylor_coeffs(fun, r, N, 2 * jmax)
    b = log_coeffs(c)
    odd = max(abs(b[k]) for k in range(1, 2 * jmax + 1, 2))
    imag = max(abs(b[k].imag) for k in range(2, 2 * jmax + 1, 2))
    p = [(-1) ** (j + 1) * j * b[2 * j].real for j in range(1, jmax + 1)]
    return p, odd, imag, c[0]


def tailj(j, T, q):
    """int_T^inf t^(-2j) (2 pi)^(-1) log(q t/2 pi) dt."""
    T = mp.mpf(T)
    n = 2 * j
    return (T ** (1 - n) * (mp.log(q * T / (2 * mp.pi)) / (n - 1) + mp.mpf(1) / (n - 1) ** 2)) / (2 * mp.pi)


def get_taylor(S, name, args):
    tay = S.setdefault('taylor', {})
    key = '%s_%s_%s' % (name, args.taylor_r, args.taylor_N)
    if key in tay:
        return [mp.mpf(x) for x in tay[key]['p']]
    mp.mp.dps = 60
    ch = CH[name]
    jmax = 12
    p, odd, imag, c0 = pj_function(lambda s: xi_mp(ch, mp.mpf(1) / 2 + s), jmax, args.taylor_r, args.taylor_N)
    p2, _, _, _ = pj_function(lambda s: xi_mp(ch, mp.mpf(1) / 2 + s), 6, args.taylor_r + 1, args.taylor_N + 32)
    dig = [int(-mp.log10(abs(p[j] - p2[j]) / abs(p[j]) + mp.mpf(10) ** -59)) for j in range(6)]
    tay[key] = {'p': [nst(x, 40) for x in p], 'max_abs_odd_log_coefficient': nst(odd, 3), 'max_abs_imag_even': nst(imag, 3),
                'xi_half': nst(c0.real, 30), 'agreeing_digits_two_radii_j1_to_6': dig}
    log('  Taylor %s: p_1 = %s, p_2 = %s, p_3 = %s (two radii agree to %s digits)' % (name, nst(p[0], 20), nst(p[1], 15), nst(p[2], 15), dig))
    return p


# --------------------------------------------------------------------------------------------------------------------- parts
def need_zeros(S, args):
    z = S.get('zeros', {})
    if all(n in z and z[n].get('T') == args.T for n in ('chi4', 'chi5')):
        return
    part_zeros(S, args)


def zero_list(S, name):
    return [mp.mpf(x) for x in S['zeros'][name]['ordinates']]


def part_zeros(S, args):
    log('== zeros (C:705): certified count and ordinates, T = %s' % args.T)
    out = S.setdefault('zeros', {})
    for name in ('chi4', 'chi5'):
        ch = CH[name]
        t0 = time.time()
        T = args.T
        h = args.h
        with prec(64):
            lenc = LEnc(ch)
            for attempt in range(4):
                try:
                    N, Nenc, theta = count_positive(ch, T, lenc)
                    break
                except RuntimeError as e:
                    log('   ', e, ': moving T by 1/32')
                    T += 1.0 / 32
            br = scan_brackets(ch, T, h)
            if len(br) != N:
                log('   sign changes %d differ from the certified count %d: refining the grid by 4' % (len(br), N))
                br = scan_brackets(ch, T, h / 4)
        fx = make_fx(ch)
        with prec(110):
            roots = [refine_zero(fx, a, b) for a, b in br]
        gam = [a2mp(c) for c, _, _ in roots]
        wid = max(float(hi - lo) for _, lo, hi in roots)
        mp.mp.dps = 40
        q = ch.q
        main = T / (2 * mp.pi) * mp.log(q * T / (2 * mp.pi * mp.e))
        thp = a2mp(theta) / mp.pi
        S_T = N - thp
        log('  %s: certified N(%s) = %d (enclosure %s), sign changes %d, all simple and on the line: %s' %
            (name, T, N, Nenc.str(12), len(br), len(br) == N))
        log('  %s: first ordinates %s, %s, %s ; max enclosure width %.2e ; main term (T/2pi)log(qT/2pi e) = %s ; theta/pi = %s ; S(T) = N - theta/pi = %s ; %.1f s' %
            (name, nst(gam[0], 12), nst(gam[1], 12), nst(gam[2], 12), wid, nst(main, 8), nst(thp, 8), nst(S_T, 6), time.time() - t0))
        out[name] = {'T': args.T, 'T_used': T, 'count': N, 'count_enclosure': Nenc.str(15), 'sign_changes': len(br),
                     'certificate': 'argument principle, flint balls, half contour 1/2 -> 2 -> 2+iT -> 1/2+iT, plus certified sign changes on a grid of step %s' % h,
                     'main_term_chapter_form': nst(main, 12), 'theta_over_pi': nst(thp, 12), 'S_T': nst(S_T, 8),
                     'chapter_main_term_minus_count': nst(main - N, 6), 'chi_minus_1': (-1) ** ch.a,
                     'ordinates': [nst(x, 28) for x in gam], 'max_enclosure_width': '%.3e' % wid,
                     'enclosures_first_five': [[lo.str(30, radius=False), hi.str(30, radius=False)] for _, lo, hi in roots[:5]]}
        if name == 'chi4':
            claim('chi_4 first ordinate', '6.0209489', gam[0], 9)
            if args.T == 300:
                claim('chi_4 count', '203', mp.mpf(N), 3)
                claim('chi_4 main term', '203.03', main, 8)
        else:
            claim('chi_5 first ordinate', '6.6484533', gam[0], 9)
            if args.T == 300:
                claim('chi_5 count', '214', mp.mpf(N), 3)
                claim('chi_5 main term', '213.68', main, 8)


def part_real(S, args):
    log('== real axis (C:707)')
    out = S.setdefault('real_axis', {})
    for name in ('chi4', 'chi5'):
        ch = CH[name]
        with prec(64):
            vals = []
            for k in range(1, 200):
                v = L_fl(ch, acb(k) / 200).real
                vals.append(float(v.mid()))
            imin = min(range(len(vals)), key=lambda i: vals[i])
            ihalf = [i for i in range(len(vals)) if i >= 99]
            imin_half = min(ihalf, key=lambda i: vals[i])
            lenc = LEnc(ch)
            pieces = []
            ok = certify_positive_real(lenc, 0.5, 2.0, pieces)
        with prec(110):
            Lhalf = L_fl(ch, acb(1) / 2).real
        mp.mp.dps = 30
        Lp0 = mp.diff(lambda s: L_mp(ch, s).real, 0, 1) if name == 'chi5' else None
        log('  %s: grid min of L over (0,1), step 0.005: %.10g at sigma = %.3f ; grid min over [1/2, 1): %.10g at %.3f ; L(1/2) = %s' %
            (name, vals[imin], 0.005 * (imin + 1), vals[imin_half], 0.005 * (imin_half + 1), Lhalf.str(20, radius=False)))
        log('  %s: L(sigma) > 0 certified on [1/2, 2] by %d ball pieces (smallest lower bound %.4g): %s' % (name, len(pieces), min(pieces) if pieces else 0, ok))
        if Lp0 is not None:
            log('  %s: L(0) = 0 (trivial zero, even character), L\'(0) = %s, so L(0.005) is about %s' % (name, nst(Lp0, 12), nst(Lp0 * 0.005, 8)))
        out[name] = {'grid_min_0_1': '%.10g' % vals[imin], 'grid_min_at_sigma': '%.3f' % (0.005 * (imin + 1)),
                     'grid_min_half_1': '%.10g' % vals[imin_half], 'grid_min_half_1_at': '%.3f' % (0.005 * (imin_half + 1)),
                     'L_half': Lhalf.str(22, radius=False), 'certified_positive_on_[1/2,2]': bool(ok), 'pieces': len(pieces),
                     'min_certified_lower_bound': '%.6g' % (min(pieces) if pieces else 0),
                     'Lprime_0': nst(Lp0, 15) if Lp0 is not None else None,
                     'conclusion': 'xi(sigma, chi) = G(sigma) L(sigma) is positive on [1/2, 1] and by xi(1 - s) = xi(s) on [0, 1]: no real zero in the critical strip'}
        if name == 'chi4':
            claim('L(1/2, chi_4)', '0.66769', mp.mpf(Lhalf.mid().str(25, radius=False)), 8)
            claim('min L chi_4 grid', '0.502', mp.mpf(vals[imin]), 8)
        else:
            claim('L(1/2, chi_5)', '0.23175', mp.mpf(Lhalf.mid().str(25, radius=False)), 8)
            claim('min L chi_5 grid', '0.0024', mp.mpf(vals[imin]), 8)


def part_classical(S, args):
    log('== classical identity and the values at 1 (C:707, C:729)')
    need_zeros(S, args)
    out = S.setdefault('classical', {})
    mp.mp.dps = 50
    for name in ('chi4', 'chi5'):
        ch = CH[name]
        p = get_taylor(S, name, args)
        mp.mp.dps = 50
        N = 64
        r = mp.mpf(1) / 2
        om = [mp.expjpi(mp.mpf(2 * m) / N) for m in range(N)]
        vals = [L_mp(ch, 1 + r * w) for w in om]
        L1 = sum(vals) / N
        L1p = sum(vals[m] * om[(-m) % N] for m in range(N)) / (N * r)
        L1 = L1.real
        L1p = L1p.real
        a = ch.a
        form = mp.log(ch.q / mp.pi) / 2 + mp.psi(0, mp.mpf(1 + a) / 2) / 2 + L1p / L1
        series = sum((-1) ** (j + 1) * p[j - 1] / mp.mpf(4) ** (j - 1) for j in range(1, len(p) + 1))
        xi_half = xi_mp(ch, mp.mpf(1) / 2).real
        xi_one = G_mp(ch, mp.mpf(1)) * L1
        T = args.T
        zs = zero_list(S, name)
        prod = mp.mpf(1)
        for g in zs:
            prod /= (1 + 1 / (4 * g * g))
        tailv = mp.quad(lambda t: mp.log(1 + 1 / (4 * t * t)) * mp.log(ch.q * t / (2 * mp.pi)) / (2 * mp.pi), [S['zeros'][name]['T_used'], 10 * S['zeros'][name]['T_used'], mp.inf])
        zprod = prod * mp.exp(-tailv)
        diff = p[0] - form
        rhs3 = p[1] / 4 - p[2] / 16
        rhs_all = sum((-1) ** j * p[j - 1] / mp.mpf(4) ** (j - 1) for j in range(2, len(p) + 1))
        log('  %s: (xi\'/xi)(1) formula = %s ; series sum (-1)^(j+1) p_j/4^(j-1) = %s ; agree to %d digits' %
            (name, nst(form, 20), nst(series, 20), int(-mp.log10(abs(form - series) + mp.mpf(10) ** -49))))
        log('  %s: xi(1/2) = %s, xi(1) = %s, E exp(-H/4) = xi(1/2)/xi(1) = %s, zero product with tail = %s (relative difference %s)' %
            (name, nst(xi_half, 15), nst(xi_one, 15), nst(xi_half / xi_one, 15), nst(zprod, 15), nst(abs(zprod / (xi_half / xi_one) - 1), 3)))
        log('  %s: p_1 - (xi\'/xi)(1) = %s ; p_2/4 - p_3/16 = %s (difference %s) ; p_2/4 - p_3/16 + ... to j = %d = %s (difference %s)' %
            (name, nst(diff, 12), nst(rhs3, 12), nst(diff - rhs3, 3), len(p), nst(rhs_all, 12), nst(diff - rhs_all, 3)))
        out[name] = {'xi_logder_at_1_formula': nst(form, 30), 'xi_logder_at_1_series': nst(series, 30),
                     'xi_half': nst(xi_half, 30), 'xi_one': nst(xi_one, 30), 'E_exp_minus_H_over_4': nst(xi_half / xi_one, 30),
                     'zero_product_with_tail': nst(zprod, 20), 'p1_minus_logder': nst(diff, 20),
                     'p2_over_4_minus_p3_over_16': nst(rhs3, 20), 'difference_three_terms': nst(diff - rhs3, 6),
                     'full_series_difference': nst(diff - rhs_all, 6), 'L_1': nst(L1, 30), 'Lprime_1': nst(L1p, 30)}
        if name == 'chi4':
            claim('(xi\'/xi)(1, chi_4) chapter', '0.0777840', form, 12)
            claim('(xi\'/xi)(1, chi_4) book', '0.0777839899618', form, 15)
            claim('xi(1/2, chi_4) chapter', '0.9807136', xi_half, 12)
            claim('xi(1/2, chi_4) book', '0.98071361405771', xi_half, 15)
            claim('xi(1, chi_4)', '1', xi_one, 15)
            claim('E exp(-H_chi4/4)', '0.9807136', xi_half / xi_one, 12)
        else:
            claim('(xi\'/xi)(1, chi_5)', '0.0782785', form, 12)


def part_table(S, args):
    log('== power sums table (C:710 to C:721)')
    need_zeros(S, args)
    out = S.setdefault('table', {})
    chap = {'chi4': {1: ('0.0780237', '0.0780148', '1.1e-4'), 2: ('9.288459e-4', '9.288458e-4', '1.0e-7'), 3: ('2.216962e-5', '2.216962e-5', '2.3e-11')},
            'chi5': {1: ('0.0784643', '0.0784571', '9.1e-5'), 2: ('7.178002e-4', '7.178001e-4', '1.1e-7'), 3: ('1.316152e-5', '1.316152e-5', '1.7e-11')}}
    for name in ('chi4', 'chi5'):
        ch = CH[name]
        p = get_taylor(S, name, args)
        mp.mp.dps = 40
        zs = zero_list(S, name)
        T = mp.mpf(args.T)
        ST = mp.mpf(S['zeros'][name]['S_T'])
        rows = {}
        for j in (1, 2, 3):
            zside = sum(g ** (-2 * j) for g in zs) + tailj(j, T, ch.q)
            fside = p[j - 1]
            rel = abs(zside - fside) / fside
            resid = zside - fside
            scaled = resid * T ** (2 * j)
            log('  %s j=%d: zero side %s, function side %s, relative error %s ; residual*T^(2j) = %s (S(T) = %s)' %
                (name, j, nst(zside, 10), nst(fside, 10), nst(rel, 2), nst(scaled, 4), nst(ST, 4)))
            # the same with the cut at the ordinate of the last zero below T (the tail integral from there)
            glast = zs[-1]
            zside2 = sum(g ** (-2 * j) for g in zs) + tailj(j, glast, ch.q)
            rel2 = abs(zside2 - fside) / fside
            log('        cut at the last zero %s: zero side %s, relative error %s' % (nst(glast, 12), nst(zside2, 10), nst(rel2, 2)))
            rows[j] = {'zero_side': nst(zside, 20), 'function_side': nst(fside, 25), 'relative_error': nst(rel, 6),
                       'residual': nst(resid, 8), 'residual_times_T_2j': nst(scaled, 8), 'S_T': nst(ST, 8),
                       'last_zero': nst(glast, 20), 'zero_side_cut_at_last_zero': nst(zside2, 20), 'relative_error_cut_at_last_zero': nst(rel2, 6)}
        out[name] = rows
        if args.T == 300:
            for j in (1, 2, 3):
                cz, cf, cr = chap[name][j]
                claim('%s j=%d zero side (cut T = 300)' % (name, j), cz, mp.mpf(rows[j]['zero_side']), 10)
                claim('%s j=%d zero side (cut at last zero)' % (name, j), cz, mp.mpf(rows[j]['zero_side_cut_at_last_zero']), 10)
                claim('%s j=%d function side' % (name, j), cf, mp.mpf(rows[j]['function_side']), 10)
                claim('%s j=%d relative error (cut T = 300)' % (name, j), cr, mp.mpf(rows[j]['relative_error']), 3)
                claim('%s j=%d relative error (cut at last zero)' % (name, j), cr, mp.mpf(rows[j]['relative_error_cut_at_last_zero']), 3)


def part_transfer(S, args):
    log('== transfer function and Hadamard product (C:723)')
    need_zeros(S, args)
    out = S.setdefault('transfer', {})
    mp.mp.dps = 40
    for name in ('chi4', 'chi5'):
        ch = CH[name]
        row = {}
        for t in (2, 5):
            v = xi_mp(ch, mp.mpf(1) / 2 + 1j * mp.mpf(t)) / xi_mp(ch, mp.mpf(1) / 2 + mp.mpf(t))
            row['T_%d' % t] = nst(v.real, 15)
            row['T_%d_imag' % t] = nst(v.imag, 3)
            log('  %s: T_chi(%d) = %s (imaginary part %s: rounding only, reality follows from the functional equation)' % (name, t, nst(v.real, 12), nst(v.imag, 2)))
        zs = zero_list(S, name)
        Tu = mp.mpf(S['zeros'][name]['T_used'])
        s = mp.mpf(3)
        prod = mp.mpf(1)
        for g in zs:
            prod *= 1 + s * s / (g * g)
        tailv = mp.quad(lambda t: mp.log(1 + s * s / (t * t)) * mp.log(ch.q * t / (2 * mp.pi)) / (2 * mp.pi), [Tu, 10 * Tu, mp.inf])
        fun = xi_mp(ch, mp.mpf(1) / 2 + s) / xi_mp(ch, mp.mpf(1) / 2)
        fun = fun.real
        glast = zs[-1]
        tailv2 = mp.quad(lambda t: mp.log(1 + s * s / (t * t)) * mp.log(ch.q * t / (2 * mp.pi)) / (2 * mp.pi), [glast, 10 * glast, mp.inf])
        row['hadamard_product_tail_from_last_zero'] = nst(prod * mp.exp(tailv2), 15)
        row['last_zero'] = nst(glast, 20)
        row['hadamard_product_no_tail'] = nst(prod, 15)
        row['hadamard_product_with_tail'] = nst(prod * mp.exp(tailv), 15)
        log('  %s: Hadamard product with the tail started at the last zero %s: %s (relative difference to xi(7/2)/xi(1/2): %s)' %
            (name, nst(glast, 10), nst(prod * mp.exp(tailv2), 12), nst(abs(prod * mp.exp(tailv2) / fun - 1), 3)))
        row['xi_7_2_over_xi_1_2'] = nst(fun, 15)
        row['relative_difference_with_tail'] = nst(abs(prod * mp.exp(tailv) / fun - 1), 4)
        row['relative_difference_no_tail'] = nst(abs(prod / fun - 1), 4)
        log('  %s: Hadamard product at s = 3: without tail %s, with tail %s, xi(7/2)/xi(1/2) = %s ; relative difference with tail %s' %
            (name, nst(prod, 12), nst(prod * mp.exp(tailv), 12), nst(fun, 12), row['relative_difference_with_tail']))
        out[name] = row
        if name == 'chi4':
            claim('T_chi4(2)', '0.53522', mp.mpf(row['T_2']), 8)
            claim('T_chi4(5)', '0.014573', mp.mpf(row['T_5']), 8)
            if args.T == 300:
                claim('chi_4 product s=3 (tail from T = 300) vs 1.952622677', '1.952622677', prod * mp.exp(tailv), 12)
                claim('chi_4 product s=3 (tail from last zero) vs 1.952622677', '1.952622677', prod * mp.exp(tailv2), 12)
            claim('chi_4 xi(7/2)/xi(1/2) vs 1.952467549', '1.952467549', fun, 12)
        else:
            if args.T == 300:
                claim('chi_5 product s=3 (tail from T = 300) vs 1.973707919', '1.973707919', prod * mp.exp(tailv), 12)
                claim('chi_5 product s=3 (tail from last zero) vs 1.973707919', '1.973707919', prod * mp.exp(tailv2), 12)
            claim('chi_5 xi(7/2)/xi(1/2) vs 1.973580921', '1.973580921', fun, 12)


def gl_nodes(a, b, n=40):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def part_negmom(S, args):
    log('== negative moment E H^(-1/2) (C:727)')
    need_zeros(S, args)
    out = S.setdefault('negmom', {})
    mp.mp.dps = 30
    pts = [0, 0.3, 1.1, 3, 6, 10, 15, 25, 40, 70]

    def eneg(xi_real, xi_half):
        return 2 * xi_half / mp.sqrt(mp.pi) * mp.quad(lambda w: 1 / xi_real(w), pts)

    half = mp.mpf(1) / 2
    ez = eneg(lambda w: xiR_mp(half + w), xiR_mp(half))
    log('  check on zeta: E H^(-1/2) = %s (chapter C:407: 6.733272312756635)' % nst(ez, 18))
    out['zeta_check'] = nst(ez, 20)
    claim('zeta E H^(-1/2)', '6.733272312756635', ez, 17)
    for name in ('chi4', 'chi5'):
        ch = CH[name]
        xh = xi_mp(ch, half).real
        e1 = eneg(lambda w: xi_mp(ch, half + w).real, xh)
        zs = np.array([float(g) for g in zero_list(S, name)])
        g2 = zs * zs
        Tu = float(S['zeros'][name]['T_used'])
        t1 = float(tailj(1, Tu, ch.q))
        t2 = float(tailj(2, Tu, ch.q))
        tot = 0.0
        edges = [0, 1, 2, 3, 4.5, 6, 8, 10, 13, 16, 20, 25, 32, 40, 55, 80]
        for lo, hi in zip(edges[:-1], edges[1:]):
            x, w = gl_nodes(lo, hi)
            f = np.array([np.exp(-np.sum(np.log1p(xx * xx / g2)) - xx * xx * t1 + xx ** 4 / 2 * t2) for xx in x])
            tot += float(np.sum(w * f))
        e2 = 2 / math.sqrt(math.pi) * tot
        log('  %s: E H^(-1/2) function side (Mellin integral on the real axis) = %s ; zero clock (%d zeros and tail) = %s ; relative difference %.2e' %
            (name, nst(e1, 15), len(zs), '%.10f' % e2, abs(e2 / float(e1) - 1)))
        out[name] = {'function_side': nst(e1, 20), 'zero_clock': '%.12f' % e2, 'relative_difference': '%.3e' % abs(e2 / float(e1) - 1)}
        if name == 'chi4':
            claim('E H_chi4^(-1/2)', '3.7574677', e1, 12)
        else:
            claim('E H_chi5^(-1/2)', '3.7097072', e1, 12)


def part_complex(S, args):
    log('== complex character modulo 5 (C:541, C:556)')
    out = S.setdefault('complex', {})
    mp.mp.dps = 40
    ps, pb = CH['psi'], CH['psibar']
    s = mp.mpf('0.7')
    half = mp.mpf(1) / 2
    r = xi_mp(ps, half + s) / xi_mp(ps, half - s)
    rb = (xi_mp(ps, half + s) * xi_mp(pb, half + s)) / (xi_mp(ps, half - s) * xi_mp(pb, half - s))
    tau = sum(mp.mpc(vr, vi) * mp.expjpi(2 * mp.mpf(b) / ps.q) for b, vr, vi in ps.tab)
    W = tau / (mp.mpc(0, 1) ** ps.a * mp.sqrt(ps.q))
    xv = xi_mp(ps, half + s)
    pred = xv / mp.conj(xv) / W
    log('  psi(2) = i, a = %d, root number W = %s (modulus %s)' % (ps.a, nst(W, 12), nst(abs(W), 12)))
    log('  xi(1/2+s, psi)/xi(1/2-s, psi) at s = 0.7 = %s, modulus %s ; W^(-1) xi/conj(xi) = %s' % (nst(r, 12), nst(abs(r), 15), nst(pred, 12)))
    log('  same ratio for xi_(psi psibar) = %s (an identity: xi_(psi psibar)(s) = xi_(psi psibar)(1 - s))' % nst(rb, 15))
    out['ratio_psi_s_0.7'] = {'real': nst(r.real, 15), 'imag': nst(r.imag, 15), 'modulus': nst(abs(r), 20)}
    out['root_number'] = {'real': nst(W.real, 15), 'imag': nst(W.imag, 15)}
    out['ratio_xi_psi_psibar_s_0.7'] = {'real': nst(rb.real, 20), 'imag': nst(rb.imag, 6)}
    claim('ratio psi real part', '0.99222', r.real, 8)
    claim('ratio psi imaginary part', '-0.12452', r.imag, 8)
    claim('ratio xi_(psi psibar)', '1', rb.real, 15)


# ---------------------------------------------------------------------------------------------------- number fields
FIELDS = {
    'Q': dict(r1=1, r2=0, d=1, poly=None),
    'Qi': dict(r1=0, r2=1, d=-4, poly=(0, 1)),     # x^2 + b x + c with (b, c) = (0, 1)
    'Q5': dict(r1=2, r2=0, d=5, poly=(-1, -1)),    # x^2 - x - 1
}


def primes_upto(n):
    sv = np.ones(n + 1, dtype=bool)
    sv[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return np.nonzero(sv)[0].astype(np.int64)


def roots_mod(poly, p):
    b, c = poly
    return sum(1 for x in range(p) if (x * x + b * x + c) % p == 0)


def classify(fd, ps):
    """0 ramified, 1 split, 2 inert, for the rational primes ps (Dedekind Kummer, rings of integers Z[x]/(f), index 1)."""
    if fd['poly'] is None:
        return np.zeros(len(ps), dtype=np.int8)
    b, c = fd['poly']
    disc = b * b - 4 * c
    typ = np.zeros(len(ps), dtype=np.int8)
    big = ps > 200
    pb = ps[big]
    if pb.size:
        e = (pb - 1) // 2
        base = np.mod(disc, pb)
        res = np.ones_like(pb)
        sq = base.copy()
        for bit in range(int(e.max()).bit_length()):
            mask = ((e >> bit) & 1).astype(bool)
            res[mask] = (res[mask] * sq[mask]) % pb[mask]
            sq = (sq * sq) % pb
        tb = np.where(res == 1, 1, 2).astype(np.int8)
        typ[big] = tb
    for i in np.nonzero(~big)[0]:
        n = roots_mod(fd['poly'], int(ps[i]))
        typ[i] = {0: 2, 1: 0, 2: 1}[n]
    # cross check of the Euler criterion against brute force on the primes up to 2000
    for p in ps[(ps > 5) & (ps < 2000)][:80]:
        i = int(np.searchsorted(ps, p))
        e1 = (int(p) - 1) // 2
        v = pow(disc % int(p), e1, int(p))
        assert (v == 1) == (roots_mod(fd['poly'], int(p)) == 2)
    return typ


def euler_logzeta(fd, s, ps, typ):
    p = ps.astype(np.float64)
    a = np.log1p(-p ** (-s))
    b = np.log1p(-p ** (-2 * s))
    tot = np.where(typ == 0, -a, np.where(typ == 1, -2 * a, -b))
    return float(np.sum(tot[::-1]))


def GammaR(s):
    return mp.pi ** (-s / 2) * mp.gamma(s / 2)


def GammaC(s):
    return 2 * (2 * mp.pi) ** (-s) * mp.gamma(s)


def xiF_from_zeta(fd, s, zF):
    return s * (s - 1) * mp.mpf(abs(fd['d'])) ** (s / 2) * GammaR(s) ** fd['r1'] * GammaC(s) ** fd['r2'] * zF


def ideal_counts(fd, nmax):
    """a_n = number of ideals of norm n, multiplicative, from the splitting type of the primes dividing n."""
    pr = [int(p) for p in primes_upto(nmax)]
    typ = classify(fd, np.array(pr, dtype=np.int64))
    ty = dict(zip(pr, [int(t) for t in typ]))
    a = [0, 1]
    for n in range(2, nmax + 1):
        m, v = n, 1
        for p in pr:
            if p * p > m:
                break
            if m % p == 0:
                k = 0
                while m % p == 0:
                    m //= p
                    k += 1
                v *= 1 if ty[p] == 0 else (k + 1 if ty[p] == 1 else (1 if k % 2 == 0 else 0))
        if m > 1:
            v *= 1 if ty[m] == 0 else (2 if ty[m] == 1 else 0)
        a.append(v)
    return a


class R2Kernel:
    """Panels of a 24 point Gauss Legendre rule for I(s, x_n) = int_(x_n)^inf t^(s-1) K_0(2t) dt, x_n = n/A, shared by all s."""

    def __init__(self, A, nmax):
        gx, gw = np.polynomial.legendre.leggauss(24)
        self.xs = [mp.mpf(n) / A for n in range(1, nmax + 1)]
        edges = self.xs + [self.xs[-1] + 2 * k for k in range(1, 21)]
        self.panels = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            t = [(hi - lo) / 2 * mp.mpf(float(x)) + (hi + lo) / 2 for x in gx]
            wk = [(hi - lo) / 2 * mp.mpf(float(w)) * mp.besselk(0, 2 * ti) for w, ti in zip(gw, t)]
            self.panels.append((t, wk))

    def I(self, s, n):
        tot = mp.mpf(0)
        for t, wk in self.panels[n - 1:]:
            tot += mp.fsum(wi * ti ** (s - 1) for ti, wi in zip(t, wk))
        return tot


def afe_V(kind, s, n, A, kern=None):
    """V_s(x) = (1/2 pi i) int_(c) gamma(s+u) x^(-s-u) du/u at x = n/A: the incomplete gamma (or Bessel) form for the gamma factor of F."""
    x = mp.mpf(n) / A
    if kind == 'C':
        return x ** (-s) * mp.gammainc(s, x, mp.inf)
    if kind == 'R2':
        return 4 * x ** (-s) * kern.I(s, n)
    if kind == 'R1':
        return x ** (-s) * mp.gammainc(s / 2, x * x, mp.inf)
    raise ValueError(kind)


def afe_Lambda(s, kind, A, a, R1, kern=None):
    """Lambda(s) = A^s gamma(s) zeta_F(s) from the coefficients a_n and the functional equation Lambda(s) = Lambda(1 - s)."""
    tot = mp.mpc(0)
    for n in range(1, len(a)):
        if a[n]:
            tot += a[n] * (afe_V(kind, s, n, A, kern) + afe_V(kind, 1 - s, n, A, kern))
    return tot - R1 * (1 / (1 - s) + 1 / s)


def afe_setup(fname):
    if fname == 'Q':
        return 'R1', 1 / mp.sqrt(mp.pi), mp.mpf(1)
    if fname == 'Qi':
        return 'C', 1 / mp.pi, mp.mpf(1) / 4
    eps = (1 + mp.sqrt(5)) / 2
    return 'R2', mp.sqrt(5) / mp.pi, 2 * mp.log(eps)


def part_dedekind(S, args):
    log('== Dedekind zeta functions (C:584)')
    out = S.setdefault('dedekind', {})
    mp.mp.dps = 30
    t0 = time.time()
    # self check of the approximate functional equation on zeta itself
    kind, A, R1 = afe_setup('Q')
    ones = [0] + [1] * args.afe_n
    s0 = mp.mpf('0.3')
    afeQ = afe_Lambda(s0, kind, A, ones, R1)
    exact = mp.pi ** (-s0 / 2) * mp.gamma(s0 / 2) * mp.zeta(s0)
    log('  AFE self check on zeta at s = 0.3: %s against %s (difference %s)' % (nst(afeQ.real, 15), nst(exact, 15), nst(abs(afeQ - exact), 3)))
    out['afe_selfcheck_zeta_0.3'] = nst(abs(afeQ - exact), 4)
    ps = primes_upto(args.P)
    for fname, ch in (('Qi', CH['chi4']), ('Q5', CH['chi5'])):
        fd = FIELDS[fname]
        typ = classify(fd, ps)
        kind, A, R1 = afe_setup(fname)
        a = ideal_counts(fd, args.afe_n)
        kern = R2Kernel(A, args.afe_n) if kind == 'R2' else None
        row = {}
        for s in (mp.mpf(2), mp.mpf('3.7'), mp.mpf('0.3')):
            lhs_prod = xiR_mp(s) * xi_mp(ch, s)
            sf = float(s)
            if s > 1:
                lz = euler_logzeta(fd, sf, ps, typ)
                tail = mp.e1((s - 1) * mp.log(args.P))
                zF = mp.exp(lz + tail)
                zF_raw = mp.exp(lz)
                xiF = xiF_from_zeta(fd, s, zF)
                xiF_raw = xiF_from_zeta(fd, s, zF_raw)
                lam = afe_Lambda(s, kind, A, a, R1, kern)
                cF = mp.mpf(abs(fd['d'])) ** (s / 2) * GammaR(s) ** fd['r1'] * GammaC(s) ** fd['r2'] / (A ** s * (mp.gamma(s) if kind == 'C' else mp.gamma(s / 2) ** 2))
                xiF_afe = s * (s - 1) * cF * lam
                row['s=%s' % sf] = {'ratio_euler_with_tail': nst((xiF / lhs_prod).real, 15), 'ratio_euler_no_tail': nst((xiF_raw / lhs_prod).real, 15),
                                    'ratio_afe': nst((xiF_afe / lhs_prod).real, 15), 'zeta_F_euler_minus_zeta_L_relative': nst(abs(zF / (mp.zeta(s) * L_mp(ch, s)) - 1), 3),
                                    'afe_vs_euler_relative': nst(abs(xiF_afe / xiF - 1), 3)}
                log('  %s s=%s: Euler product (P = %d, tail by the prime number theorem) ratio %s, without tail %s ; AFE ratio %s ; AFE against Euler %s' %
                    (fname, sf, args.P, nst((xiF / lhs_prod).real, 12), nst((xiF_raw / lhs_prod).real, 12), nst((xiF_afe / lhs_prod).real, 12), nst(abs(xiF_afe / xiF - 1), 3)))
            else:
                lam = afe_Lambda(s, kind, A, a, R1, kern)
                cF = mp.mpf(abs(fd['d'])) ** (s / 2) * GammaR(s) ** fd['r1'] * GammaC(s) ** fd['r2'] / (A ** s * (mp.gamma(s) if kind == 'C' else mp.gamma(s / 2) ** 2))
                xiF_afe = s * (s - 1) * cF * lam
                row['s=%s' % sf] = {'ratio_afe': nst((xiF_afe / lhs_prod).real, 15), 'imag': nst((xiF_afe / lhs_prod).imag, 3), 'completion_constant': nst(cF.real, 15)}
                log('  %s s=%s: AFE built from the ideal counts a_n (n <= %d) and the gamma factor of F alone: ratio %s ; constant between xi_F and s(s-1)Lambda = %s' %
                    (fname, sf, args.afe_n, nst((xiF_afe / lhs_prod).real, 12), nst(cF.real, 12)))
        out[fname] = row
    claim('Q(i) ratio s=2', '1', mp.mpf(out['Qi']['s=2.0']['ratio_afe']), 12)
    claim('Q(sqrt5) ratio s=2', '2', mp.mpf(out['Q5']['s=2.0']['ratio_afe']), 12)
    claim('Q(i) ratio s=0.3', '1', mp.mpf(out['Qi']['s=0.3']['ratio_afe']), 12)
    claim('Q(sqrt5) ratio s=0.3', '2', mp.mpf(out['Q5']['s=0.3']['ratio_afe']), 12)
    # additivity of E H_F for Q(i): p_1 of xi_F built from the AFE, against p_1(zeta) + p_1(chi_4)
    kind, A, R1 = afe_setup('Qi')
    fd = FIELDS['Qi']
    a = ideal_counts(fd, args.afe_n)
    mp.mp.dps = 25

    def xiF_afe_fun(s):
        lam = afe_Lambda(s, kind, A, a, R1)
        return 2 * s * (s - 1) * lam

    pF, _, _, _ = pj_function(lambda s: xiF_afe_fun(mp.mpf(1) / 2 + s), 2, mp.mpf(2), args.afe_circle)
    pz, _, _, _ = pj_function(lambda s: xiR_mp(mp.mpf(1) / 2 + s), 2, mp.mpf(2), 48)
    mp.mp.dps = 50
    p4 = get_taylor(S, 'chi4', args)
    mp.mp.dps = 25
    log('  E H_F for Q(i) from xi_F built by the AFE: p_1(F) = %s ; p_1(zeta) = %s ; p_1(chi_4) = %s ; sum = %s' % (nst(pF[0], 12), nst(pz[0], 12), nst(p4[0], 12), nst(pz[0] + p4[0], 12)))
    out['EH_Qi'] = {'p1_F_afe': nst(pF[0], 15), 'p1_zeta': nst(pz[0], 15), 'p1_chi4': nst(p4[0], 15), 'sum': nst(pz[0] + p4[0], 15), 'p2_F_afe': nst(pF[1], 12), 'p2_zeta_plus_chi4': nst(pz[1] + p4[1], 12)}
    claim('E H_zeta', '0.0231050', pz[0], 10)
    claim('E H_chi4', '0.0780148', p4[0], 10)
    claim('E H_F = sum', '0.1011198', pF[0], 10)
    log('  (%.1f s)' % (time.time() - t0))


def phi_mp(x, s):
    z = s * x
    if z < mp.mpf('0.01'):
        return sum((-z) ** k / mp.factorial(k) for k in range(2, 14))
    return mp.expm1(-z) + z


def levy_prime_sum(fname, alpha, s, N):
    fd = FIELDS[fname]
    ps = primes_upto(N)
    typ = classify(fd, ps)
    p = ps.astype(np.float64)
    lp = np.log(p)
    tot = 0.0
    if fname == 'Q':
        groups = [(np.ones(len(ps), dtype=bool), 1, 1.0)]
    else:
        groups = [(typ == 0, 1, 1.0), (typ == 1, 1, 2.0), (typ == 2, 2, 1.0)]
    for mask, e0, c in groups:
        pm, lm = p[mask], lp[mask]
        k = 1
        while True:
            n = pm ** (e0 * k)
            keep = n <= N
            if not keep.any():
                break
            x = e0 * k * lm[keep]
            nn = n[keep]
            z = s * x
            phi = np.where(z < 0.01, z * z / 2 - z ** 3 / 6 + z ** 4 / 24, np.expm1(-z) + z)
            tot += float(np.sum((c / k) * nn ** (-alpha) * phi))
            k += 1
    return tot


def part_levy(S, args):
    log('== Levy decomposition (C:615)')
    out = S.setdefault('levy', {})
    mp.mp.dps = 30
    t0 = time.time()
    cases = [('Q', 3, 1, 'xi'), ('Qi', 3, mp.mpf(3) / 2, 'chi4'), ('Q5', mp.mpf(5) / 2, 1, 'chi5')]
    chap = {'Q': '5.05e-10', 'Qi': '7.83e-10', 'Q5': '1.17e-7'}
    half = mp.mpf(1) / 2

    def ld_xi(s):
        return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.psi(0, s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s)

    def ld_L(ch, s):
        Lv = L_mp(ch, s)
        Lp = -mp.log(ch.q) * Lv + mp.mpf(ch.q) ** (-s) * sum(vr * mp.zeta(s, mp.mpf(b) / ch.q, 1) for b, vr, vi in ch.tab)
        return mp.log(ch.q / mp.pi) / 2 + mp.psi(0, (s + ch.a) / 2) / 2 + Lp / Lv

    for fname, alpha, s, kind in cases:
        alpha = mp.mpf(alpha)
        s = mp.mpf(s)
        fd = FIELDS[fname]
        if fname == 'Q':
            lhs = mp.log(xiR_mp(alpha + s) / xiR_mp(alpha)) - s * ld_xi(alpha)
        else:
            ch = CH[kind]
            lhs = mp.log(xiR_mp(alpha + s) * xi_mp(ch, alpha + s) / (xiR_mp(alpha) * xi_mp(ch, alpha))) - s * (ld_xi(alpha) + ld_L(ch, alpha))
        lhs = lhs.real if hasattr(lhs, 'real') else lhs

        def integ(f):
            return mp.quad(f, [0, mp.mpf(1) / 10, 1, 5, 20, 60, 120])

        def rr(x):
            return phi_mp(x, s) * mp.exp(-alpha * x) / (-mp.expm1(-2 * x)) / x

        def cc(x):
            return phi_mp(x, s) * mp.exp(-alpha * x) / (-mp.expm1(-x)) / x

        def ll(x):
            return -phi_mp(x, s) * mp.exp(-alpha * x) * (1 + mp.exp(x)) / x

        IR = integ(rr)
        IC = integ(cc)
        IL = integ(ll)
        binetR = mp.log(GammaR(alpha + s) / GammaR(alpha)) - s * (-mp.log(mp.pi) / 2 + mp.psi(0, alpha / 2) / 2)
        binetC = mp.log(GammaC(alpha + s) / GammaC(alpha)) - s * (-mp.log(2 * mp.pi) + mp.psi(0, alpha))
        gam = fd['r1'] * IR + fd['r2'] * IC + IL
        row = {'alpha': nst(alpha, 5), 's': nst(s, 5), 'lhs': nst(lhs, 20), 'gamma_and_polynomial_integrals': nst(gam, 20),
               'binet_R_check': nst(IR - binetR, 3), 'binet_C_check': nst(IC - binetC, 3) if fd['r2'] else None, 'by_N': {}}
        log('  %s (alpha, s) = (%s, %s): LHS %s ; gamma and polynomial integrals %s ; Binet checks %s, %s' %
            (fname, nst(alpha, 4), nst(s, 4), nst(lhs, 16), nst(gam, 16), nst(IR - binetR, 3), nst(IC - binetC, 3) if fd['r2'] else 'n/a'))
        for NN in sorted({30000, args.N}):
            sN = levy_prime_sum(fname, float(alpha), float(s), NN)
            tail = mp.quad(lambda x: x ** (-alpha) * phi_mp(mp.log(x), s) / mp.log(x), [NN, 10 * NN, mp.inf])
            rhs = gam + mp.mpf(sN)
            rhs_t = rhs + tail
            log('     N = %d: prime ideal sum %s ; difference LHS - RHS %s ; density one tail %s ; difference with tail %s' %
                (NN, nst(mp.mpf(sN), 16), nst(lhs - rhs, 4), nst(tail, 4), nst(lhs - rhs_t, 4)))
            row['by_N'][str(NN)] = {'prime_ideal_sum': nst(mp.mpf(sN), 20), 'difference': nst(lhs - rhs, 6), 'tail_density_one': nst(tail, 8),
                                    'difference_with_tail': nst(lhs - rhs_t, 6)}
            if NN == 30000:
                claim('%s difference at N = 30000 (chapter value, N not stated)' % fname, chap[fname], abs(lhs - rhs), 4)
        out[fname] = row
    log('  (%.1f s)' % (time.time() - t0))


ORDER = ['zeros', 'real', 'classical', 'table', 'transfer', 'negmom', 'complex', 'dedekind', 'levy']
PARTS = {'zeros': part_zeros, 'real': part_real, 'classical': part_classical, 'table': part_table, 'transfer': part_transfer,
         'negmom': part_negmom, 'complex': part_complex, 'dedekind': part_dedekind, 'levy': part_levy}


def check_tables():
    with prec(64):
        for name, ch in CH.items():
            fc = dirichlet_char(*ch.conrey)
            for b, vr, vi in ch.tab:
                assert float(abs(fc(b) - acb(vr, vi))) < 1e-9, (name, b)
            s = acb(0.5, 17.25)
            dl = s.dirichlet_l(fc)
            hz = L_fl(ch, s)
            assert float(abs(dl - hz)) < 1e-12, (name, dl, hz)
    log('character tables agree with flint (Conrey labels) and Hurwitz L agrees with dirichlet_l at 1/2 + 17.25 i')


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--parts', default=','.join(ORDER))
    ap.add_argument('--T', type=float, default=None)
    ap.add_argument('--h', type=float, default=None)
    ap.add_argument('--P', type=int, default=None)
    ap.add_argument('--N', type=int, default=None)
    ap.add_argument('--afe-n', dest='afe_n', type=int, default=None)
    ap.add_argument('--afe-circle', dest='afe_circle', type=int, default=None)
    ap.add_argument('--taylor-N', dest='taylor_N', type=int, default=None)
    ap.add_argument('--taylor-r', dest='taylor_r', type=float, default=2.0)
    args = ap.parse_args()
    q = args.quick
    args.T = args.T if args.T is not None else (50.0 if q else 300.0)
    args.h = args.h if args.h is not None else (1 / 16 if q else 1 / 32)
    args.P = args.P if args.P is not None else (200000 if q else 10000000)
    args.N = args.N if args.N is not None else (30000 if q else 2000000)
    args.afe_n = args.afe_n if args.afe_n is not None else (12 if q else 30)
    args.afe_circle = args.afe_circle if args.afe_circle is not None else (24 if q else 48)
    args.taylor_N = args.taylor_N if args.taylor_N is not None else (64 if q else 96)
    out = args.out or ('/tmp/gl_dirichlet_quick.json' if q else 'gl_dirichlet.json')
    LOGF = open(os.path.splitext(out)[0] + '.log', 'a')
    log('# gl_dirichlet.py %s' % ' '.join(sys.argv[1:]))
    log('# started %s ; python-flint arb, mpmath %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), mp.__version__))
    S = {}
    if os.path.exists(out):
        try:
            S = json.load(open(out))
        except Exception:
            S = {}
    S['script'] = 'gl_dirichlet.py'
    S['mode'] = 'quick' if q else 'full'
    S['args'] = {k: v for k, v in vars(args).items()}
    check_tables()
    todo = [p for p in ORDER if p in args.parts.split(',')]
    for name in todo:
        t0 = time.time()
        PARTS[name](S, args)
        S.setdefault('parts_done', {})[name] = {'T': args.T, 'seconds': round(time.time() - t0, 1), 'mode': S['mode']}
        json.dump(S, open(out, 'w'), indent=1)
    log('# done, wrote %s ; total %.1f s' % (out, time.time() - START))


START = time.time()
if __name__ == '__main__':
    main()
