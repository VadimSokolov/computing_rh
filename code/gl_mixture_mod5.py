# gl_mixture_mod5.py: the mixture of clocks modulo 5 of Chapter 21 (Thorin Representations, Dirichlet L Functions and the Linnik
# Constant, label prefix gl:), group L2b of the check of that chapter: the last paragraph of Section "The family: sum and mixture"
# (C:753, Proposition gl:prop:mixture applied to the family of characters modulo 5). Added for the book in October 2026, written from
# the mathematics (the chapter author's scripts were not read).
# Objects. chi_5 the Legendre symbol modulo 5 (even), psi the complex character modulo 5 with psi(2) = i (so psi(3) = -i, psi(4) = -1,
# odd, a = 1), psibar its conjugate. xi(s, chi) = (q/pi)^((s+a)/2) Gamma((s+a)/2) L(s, chi). For y real
#   A(y) = xi(1/2, chi_5)/xi(1/2 + iy, chi_5),   B(y) = xi_(psi psibar)(1/2)/xi_(psi psibar)(1/2 + iy),
#   xi_(psi psibar)(s) = xi(s, psi) xi(s, psibar), M_lambda = lambda A + (1 - lambda) B,
# the Laplace transforms of the two clocks on the negative axis u = -y^2 (poles at the ordinates of the zeros).
# Parts (--parts, default all).
#  zeros   (1) The zeros of xi(1/2 + it, chi_5) for 0 < t < 12 and of L(s, psi) and L(s, psibar) for |Im s| < 12. A certified count of
#          the zeros of L with 0 < Im s < T by the argument principle with interval arithmetic (python-flint balls; N = (theta(T) + Delta
#          arg L along 1/2 -> 2 -> 2+iT -> 1/2+iT)/pi, the half contour of [-1, 2] x [0, T], valid for every primitive character because
#          xi(1 - conj s) = conj xi(s) up to the unimodular root number); for psi the zeros with negative ordinates are those of psibar
#          with the sign changed, counted as zeros of L(s, psibar) with 0 < Im s < T. The sign changes of the real function
#          Re(W^(-1/2) xi(1/2 + it, chi)) (W the root number) on the grid are counted; their number equals the certified count, so all
#          the zeros are on the line and simple. Ordinates refined by Illinois iteration at 128 bits, enclosed by two certified signs.
#          The ordinates of psibar are checked to be the negatives of those of psi. The atoms of the two clocks in the variable y
#          (|gamma| for each ordinate of psi, both signs; the ordinates of chi_5) and the counts of atoms below the two ends of the
#          intervals of Proposition gl:prop:mixture (b) (three below 6.648, five below 9.443, odd) are recorded.
#  mixture (2) For lambda in {0.01, 0.1, 0.5, 0.9, 0.99} and the two intervals (6.6484533, 8.4572292), (9.4429311, 9.8314444) between a
#          pole of A and a pole of B: every zero of M_lambda. With D_A = xi(1/2 + iy, chi_5) and D_B = xi_(psi psibar)(1/2 + iy) (both real)
#          the zeros of M_lambda in an interval free of poles are the zeros of the entire function N = lambda c_A D_B + (1 - lambda) c_B D_A
#          (c_A = D_A(0), c_B = D_B(0)), which stays of size one at the poles. N is evaluated at 128 bits on a nonuniform grid (geometric
#          offsets 1e-9 to 0.5 times the length from both ends, so a zero at 3e-5 from a pole is resolved, plus a uniform grid); every
#          sign change of N is a certified bracket (two evaluations with certified opposite signs, N continuous), refined by bisection
#          to width 1e-22; the bracket is checked to contain no pole (balls over the bracket exclude 0 for D_A and D_B). Uniqueness:
#          the argument principle for N on the rectangle [pole, pole] x [-eta, eta] (eta = 0.25, 0.1, 0.03, first that succeeds) gives
#          the number of zeros of N inside; equality with the number of real sign changes shows there is no other zero (a nonreal
#          pair would add two). The values M_lambda at the distances 1e-9, 1e-7, 1e-5, 2.4e-5 from each of the four poles, on
#          both sides, are recorded (the claim "values of size at least 5e5").
#  central (3) The central values L(1/2, chi_5), L(1/2, psi), L(1/2, psibar) at 128 bits, by Hurwitz zeta and by flint's dirichlet_l (they
#          agree to 1e-30); the enclosure of |L| excludes 0, which is the hypothesis L(1/2, chi) != 0 of the clocks of Propositions gl:prop:Lreal
#          and gl:prop:Lcomplex for the members modulo 5 (finding M5 of misc/reviews/nick-linnik/final-X21.md). Under 0.1 s.
# Error control. Hurwitz L against flint's dirichlet_l (Conrey characters); the reality of the rotated functions and of D_B on the line is a
# rounding check only (it follows from the functional equation); the zero counts agree with the sign changes.
# Output: --out (default gl_mixture_mod5.json in the current folder, /tmp/gl_mixture_quick.json with --quick) with the log beside it (same
# stem, .log); numbers are strings. A partial run (--parts, --lambdas) merges into an existing output file.
# Usage: python3 gl_mixture_mod5.py [--quick] [--out FILE] [--parts zeros,central,mixture] [--lambdas 0.01,0.1,0.5,0.9,0.99] [--no-rect]
# --quick (lambda = 1/2, no rectangles) runs in about 1 s locally; the full case (all five lambdas, with the rectangle certificates) took
# 8.9 s on a laptop (one process); as one job on the Hopper cluster:
#   bash misc/tools/hopper_run.sh -c 1 -t 60 -g gl_mixture_mod5.json -g gl_mixture_mod5.log code/gl_mixture_mod5.py -- --out gl_mixture_mod5.json
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext

import mpmath as mp
from flint import acb, arb, ctx, dirichlet_char

LOGF = None
START = time.time()


def log(*a):
    msg = ' '.join(str(x) for x in a)
    print(msg, flush=True)
    if LOGF is not None:
        LOGF.write(msg + '\n')
        LOGF.flush()


class prec:
    def __init__(self, p):
        self.p = p

    def __enter__(self):
        self.old = ctx.prec
        ctx.prec = self.p

    def __exit__(self, *e):
        ctx.prec = self.old


class Char:
    def __init__(self, name, q, tab, conrey):
        self.name = name
        self.q = q
        self.tab = [(b, v[0], v[1]) for b, v in sorted(tab.items())]
        m1 = dict(tab)[q - 1]
        self.a = 0 if m1 == (1, 0) else 1
        self.conrey = conrey


CHI5 = Char('chi5', 5, {1: (1, 0), 2: (-1, 0), 3: (-1, 0), 4: (1, 0)}, (5, 4))
PSI = Char('psi', 5, {1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}, (5, 2))
PSIB = Char('psibar', 5, {1: (1, 0), 2: (0, -1), 3: (0, 1), 4: (-1, 0)}, (5, 3))


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


def xi_line(ch, y):
    """xi(1/2 + i y, ch) for y an acb (real or complex)."""
    s = acb(0.5) + acb(0, 1) * y
    return G_fl(ch, s) * L_fl(ch, s)


def a2mp(x, digits=34):
    return mp.mpf(x.mid().str(digits, radius=False))


def nst(x, d=25):
    return mp.nstr(x, d)


def agree(claimed, mine):
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


def claim(label, claimed, mine, digits=10):
    st = agree(claimed, mine)
    log('  claim %-52s chapter %-12s mine %-18s %s' % (label, claimed, nst(mine, digits), st))
    return {'chapter': claimed, 'mine': nst(mine, 16), 'status': st}


class LEnc:
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
    return float(E.rad()) < 0.5 * float(abs(E.mid()))


def chain_arg(f, a, b, maxdepth=34):
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
    theta = logG_fl(ch, acb(0.5, T)).imag
    P1 = chain_arg(lenc, complex(0.5, 0.0), complex(2.0, 0.0))
    P2 = (L_fl(ch, acb(2.0, T)) / L_fl(ch, acb(2.0, 0.0))).arg()
    P3 = chain_arg(lenc, complex(2.0, T), complex(0.5, T))
    tot = (theta + P1 + P2 + P3) / arb.pi()
    n = int(round(float(tot)))
    if not (tot - n).abs_upper() < 1e-3:
        raise RuntimeError('count_positive: enclosure %s is not within 1e-3 of an integer' % tot.str(12))
    return n, tot


def root_rotation(ch):
    """W^(-1/2) for a primitive character: xi(s, chi) = W xi(1 - s, chibar), W = tau(chi)/(i^a sqrt q)."""
    tau = acb(0)
    for b, vr, vi in ch.tab:
        tau = tau + acb(vr, vi) * (2 * arb.pi() * acb(0, 1) * b / ch.q).exp()
    W = tau / (acb(0, 1) ** ch.a * acb(ch.q).sqrt())
    return W, W ** acb(-0.5)


class Line:
    """The real function t -> Re(rot * xi(1/2 + i t, ch)) and the largest imaginary part seen (a rounding check)."""

    def __init__(self, ch, rot=None):
        self.ch, self.rot = ch, rot
        self.imax = 0.0

    def __call__(self, t):
        v = xi_line(self.ch, acb(t))
        if self.rot is not None:
            v = v * self.rot
        self.imax = max(self.imax, float(abs(v.imag).upper()))
        return v.real


def scan_brackets(fx, t0, t1, h):
    n = int(round((t1 - t0) / h))
    out = []
    pt = t0
    ps = sgn(fx(arb(t0)))
    for k in range(1, n + 1):
        t = t0 + k * h
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


def refine_zero(fx, a, b, tol=1e-24, maxit=80):
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
    d = max(abs(float(b - a)), 1e-28)
    for _ in range(40):
        lo, hi = c - d, c + d
        s1, s2 = sgn(fx(lo)), sgn(fx(hi))
        if s1 != 0 and s2 != 0 and s1 != s2:
            return c, lo, hi
        d *= 4
    raise RuntimeError('refine_zero: no certified bracket near %s' % c.str(20))


def check_tables():
    with prec(64):
        for ch in (CHI5, PSI, PSIB):
            fc = dirichlet_char(*ch.conrey)
            for b, vr, vi in ch.tab:
                assert float(abs(fc(b) - acb(vr, vi))) < 1e-9, (ch.name, b)
            s = acb(0.5, 9.25)
            assert float(abs(s.dirichlet_l(fc) - L_fl(ch, s))) < 1e-12
    log('character tables agree with flint (Conrey labels 4, 2, 3 modulo 5) and Hurwitz L agrees with dirichlet_l at 1/2 + 9.25 i')


# --------------------------------------------------------------------------------------------------------------------- zeros
def part_zeros(S, args):
    T = args.T
    log('== zeros below %s of xi(1/2 + it, chi_5) and of L(s, psi), L(s, psibar) (C:753)' % T)
    out = S.setdefault('zeros', {})
    h = args.h
    res = {}
    t0 = time.time()
    with prec(64):
        l5, lp, lb = LEnc(CHI5), LEnc(PSI), LEnc(PSIB)
        N5, N5e = count_positive(CHI5, T, l5)
        Np, Npe = count_positive(PSI, T, lp)
        Nb, Nbe = count_positive(PSIB, T, lb)
        Wp, rotp = root_rotation(PSI)
        Wb, rotb = root_rotation(PSIB)
        f5, fp_, fb_ = Line(CHI5), Line(PSI, rotp), Line(PSIB, rotb)
        br5 = scan_brackets(f5, 0.0, T, h)
        brp = scan_brackets(fp_, -T, T, h)
        brb = scan_brackets(fb_, -T, T, h)
    log('  certified counts (0 < Im s < %s): chi_5 %d, psi %d, psibar %d (enclosures %s, %s, %s)' % (T, N5, Np, Nb, N5e.str(10), Npe.str(10), Nbe.str(10)))
    log('  sign changes of the real functions: chi_5 on (0, T): %d ; psi on (-T, T): %d ; psibar on (-T, T): %d' % (len(br5), len(brp), len(brb)))
    log('  zeros of psi with negative ordinate = zeros of psibar with positive ordinate: counts %d = %d ; with positive ordinate %d' % (
        sum(1 for a, b in brp if b <= 0), Nb, Np))
    with prec(128):
        z5 = [refine_zero(f5, a, b) for a, b in br5]
        zp = [refine_zero(fp_, a, b) for a, b in brp]
        zb = [refine_zero(fb_, a, b) for a, b in brb]
    mp.mp.dps = 40
    g5 = [a2mp(c) for c, _, _ in z5]
    gp = [a2mp(c) for c, _, _ in zp]
    gb = [a2mp(c) for c, _, _ in zb]
    neg = max(abs(x - y) for x, y in zip(sorted(gp), sorted([-v for v in gb]))) if len(gp) == len(gb) else None
    log('  chi_5 ordinates: %s' % ', '.join(nst(x, 11) for x in g5))
    log('  psi ordinates (signed, increasing): %s' % ', '.join(nst(x, 11) for x in gp))
    log('  psibar ordinates (signed): %s' % ', '.join(nst(x, 11) for x in gb))
    log('  psibar ordinates are the negatives of the psi ordinates: largest difference %s' % (nst(neg, 3) if neg is not None else 'counts differ'))
    log('  root numbers W(psi) = %s, W(psibar) = %s ; largest imaginary part of the rotated functions: %.2e, %.2e, %.2e (rounding only)' % (
        Wp.str(12, radius=False), Wb.str(12, radius=False), f5.imax, fp_.imax, fb_.imax))
    atoms_pair = sorted(abs(x) for x in gp)
    atoms_chi5 = sorted(g5)
    union = sorted([(x, 'chi5') for x in atoms_chi5] + [(x, 'pair') for x in atoms_pair])
    below_a = sum(1 for x, _ in union if x <= g5[0] + mp.mpf(10) ** -20)
    pair_next = [x for x, w in union if w == 'pair' and x > g5[0]][0]
    p1 = [x for x in atoms_pair if x > g5[0]][0]
    p2 = [x for x in atoms_pair if x < g5[1]][-1]
    below_b = sum(1 for x, _ in union if x <= p2 + mp.mpf(10) ** -20)
    log('  atoms in y (union of the two measures): %s' % ', '.join('%s (%s)' % (nst(x, 9), w) for x, w in union))
    log('  atoms in (0, %s]: %d (odd) and the next atom %s belongs to the pair; atoms in (0, %s]: %d (odd) and the next atom %s belongs to chi_5' % (
        nst(g5[0], 9), below_a, nst(pair_next, 9), nst(p2, 9), below_b, nst(g5[1], 9)))
    out.update({'T': T, 'h': h, 'certified_counts': {'chi5': N5, 'psi_positive': Np, 'psibar_positive': Nb},
                'sign_changes': {'chi5': len(br5), 'psi': len(brp), 'psibar': len(brb)},
                'chi5_ordinates': [nst(x, 28) for x in g5], 'psi_ordinates_signed': [nst(x, 28) for x in gp],
                'psibar_ordinates_signed': [nst(x, 28) for x in gb], 'psibar_minus_psi_negatives_max_difference': nst(neg, 3) if neg is not None else None,
                'enclosure_half_widths_max': '%.2e' % max(float(hi - lo) / 2 for _, lo, hi in z5 + zp + zb),
                'atoms': [[nst(x, 20), w] for x, w in union], 'atoms_up_to_first_chi5': below_a, 'atoms_up_to_pair_end_of_second_interval': below_b,
                'root_number_psi': Wp.str(15, radius=False), 'root_number_psibar': Wb.str(15, radius=False),
                'certificate': 'argument principle, flint balls, half contour 1/2 -> 2 -> 2+iT -> 1/2+iT; sign changes of the rotated real function on a grid of step %s' % h})
    claim('chi_5 first ordinate', '6.6484533', g5[0], 9)
    claim('chi_5 second ordinate', '9.8314444', g5[1], 9)
    claim('chi_5 third ordinate', '11.958846', g5[2], 9)
    gps = sorted(gp, key=lambda x: abs(x))
    for k, cl in enumerate(('-4.1329037', '6.1835782', '8.4572292', '-9.4429311')):
        claim('psi ordinate %d by modulus (signed)' % (k + 1), cl, gps[k], 9)
    claim('psi ordinate 5 by modulus (chapter prints no sign)', '11.282896', abs(gps[4]), 9)
    log('  sign of the fifth ordinate of psi by modulus: %s (the chapter prints 11.282896 without a sign)' % ('negative' if gps[4] < 0 else 'positive'))
    out['psi_ordinates_by_modulus'] = [nst(x, 20) for x in gps]
    log('  (%.1f s)' % (time.time() - t0))


# ------------------------------------------------------------------------------------------------------------------ mixture
class Mix:
    def __init__(self):
        with prec(256):
            self.c5 = xi_line(CHI5, acb(0)).real
            self.cB = (xi_line(PSI, acb(0)) * xi_line(PSIB, acb(0))).real
        self.imax = 0.0

    def pair(self, y):
        """(D_A(y), D_B(y)) real balls at the working precision, y an arb."""
        x5 = xi_line(CHI5, acb(y))
        xp = xi_line(PSI, acb(y))
        xb = xi_line(PSIB, acb(y))
        db = xp * xb
        self.imax = max(self.imax, float(abs(x5.imag).upper()), float(abs(db.imag).upper()))
        return x5.real, db.real

    def N(self, lam, da, db):
        return lam * self.c5 * db + (1 - lam) * self.cB * da

    def M(self, lam, da, db):
        return lam * self.c5 / da + (1 - lam) * self.cB / db

    def Ncomplex(self, lam):
        def f(y):
            x5 = xi_line(CHI5, y)
            db = xi_line(PSI, y) * xi_line(PSIB, y)
            return lam * self.c5 * db + (1 - lam) * self.cB * x5
        return f


def offsets(length, n_geo=32, lo=1e-9):
    r = (0.5 / lo) ** (1.0 / (n_geo - 1))
    return [length * lo * r ** k for k in range(n_geo)]


def part_central(S, args):
    log('== central values L(1/2, chi) for chi_5, psi, psibar (the nonvanishing hypothesis of the clocks modulo 5)')
    out = S.setdefault('central', {})
    with prec(128):
        half = acb(1) / 2
        for ch in (CHI5, PSI, PSIB):
            v = L_fl(ch, half)
            w = half.dirichlet_l(dirichlet_char(*ch.conrey))
            assert float(abs(v - w)) < 1e-30, ch.name
            out[ch.name] = {'re': v.real.str(20, radius=False), 'im': v.imag.str(20, radius=False), 'abs': abs(v).str(20, radius=False),
                            'nonzero_certified': bool(abs(v) > 0), 'hurwitz_minus_dirichlet_l': '%.1e' % float(abs(v - w))}
            log('  L(1/2, %s) = %s + %s i, |L| = %s, enclosure excludes 0: %s' % (
                ch.name, out[ch.name]['re'], out[ch.name]['im'], out[ch.name]['abs'], out[ch.name]['nonzero_certified']))


def part_mixture(S, args):
    log('== zeros of lambda A + (1 - lambda) B (C:753)')
    need = S.get('zeros')
    if not need or 'chi5_ordinates' not in need:
        part_zeros(S, args)
    Z = S['zeros']
    mp.mp.dps = 40
    g5 = [mp.mpf(x) for x in Z['chi5_ordinates']]
    gp = [mp.mpf(x) for x in Z['psi_ordinates_signed']]
    pair = sorted(abs(x) for x in gp)
    intervals = {'I1': (g5[0], [x for x in pair if x > g5[0]][0]),
                 'I2': ([x for x in pair if x < g5[1]][-1], g5[1])}
    log('  intervals: I1 = (%s, %s) from a pole of A to a pole of B ; I2 = (%s, %s) from a pole of B to a pole of A' % (
        nst(intervals['I1'][0], 10), nst(intervals['I1'][1], 10), nst(intervals['I2'][0], 10), nst(intervals['I2'][1], 10)))
    out = S.setdefault('mixture', {})
    out['intervals'] = {k: [nst(v[0], 25), nst(v[1], 25)] for k, v in intervals.items()}
    mix = Mix()
    lams = [(round(float(x) * 100), 100) for x in args.lambdas.split(',')]
    claims_zero = {('0.5', 'I1'): '6.6786579', ('0.5', 'I2'): '9.8285930', ('0.1', 'I1'): '6.6516953', ('0.9', 'I1'): '6.9811236', ('0.9', 'I2'): '9.8069278'}
    for key, (pl_mp, pr_mp) in intervals.items():
        pl = arb(nst(pl_mp, 40))
        pr = arb(nst(pr_mp, 40))
        with prec(128):
            length = float(pr - pl)
            off = offsets(length)
            nuni = args.nuni if key == 'I1' else max(args.nuni // 4, 40)
            pts = [pl + o for o in off] + [pl + length * (k / nuni) for k in range(1, nuni)] + [pr - o for o in reversed(off)]
            pts = sorted({p.mid().str(36, radius=False): p.mid() for p in pts}.items(), key=lambda kv: float(kv[1]))
            pts = [p for _, p in pts]
            vals = [mix.pair(p) for p in pts]
        log('  %s: grid of %d points between the poles, offsets from %.1e to %.2f of the length from both ends' % (key, len(pts), off[0] / length, off[-1] / length))
        for num, den in lams:
            lam = arb(num) / den
            lname = ('%g' % (num / den))
            row = {'zeros': [], 'sign_changes': 0}
            with prec(128):
                sg = [sgn(mix.N(lam, da, db)) for da, db in vals]
                cuts = [(i, i + 1) for i in range(len(pts) - 1) if sg[i] != sg[i + 1] and sg[i] != 0 and sg[i + 1] != 0]
                zeros = []
                for i, j in cuts:
                    a, b = pts[i], pts[j]
                    sa = sg[i]
                    while float(b - a) > 1e-22:
                        m = ((a + b) / 2).mid()
                        sm = sgn(mix.N(lam, *mix.pair(m)))
                        if sm == 0:
                            a = b = m
                            break
                        if sm == sa:
                            a = m
                        else:
                            b = m
                    ball = (a + (b - a) / 2)
                    ballr = arb(ball.mid(), float(b - a) / 2 + 1e-30)
                    dA, dB = mix.pair(ballr)
                    nopole = (not (dA.lower() <= 0 <= dA.upper())) and (not (dB.lower() <= 0 <= dB.upper()))
                    sA, sB = sgn(dA), sgn(dB)
                    zeros.append({'y': a2mp(ball), 'lo': a, 'hi': b, 'no_pole_in_bracket': bool(sA != 0 and sB != 0)})
            row['sign_changes'] = len(zeros)
            row['zeros'] = [{'y': nst(z['y'], 24), 'bracket_lo': z['lo'].str(26, radius=False), 'bracket_hi': z['hi'].str(26, radius=False),
                             'no_pole_in_bracket': z['no_pole_in_bracket'],
                             'distance_to_left_pole': nst(z['y'] - pl_mp, 6), 'distance_to_right_pole': nst(pr_mp - z['y'], 6)} for z in zeros]
            # rectangle certificate
            rect = None
            if args.rect:
                for eta in (0.25, 0.1, 0.03):
                    try:
                        with prec(64):
                            f = mix.Ncomplex(lam)
                            xl, xr = float(pl_mp), float(pr_mp)
                            tot = (chain_arg(f, complex(xl, -eta), complex(xr, -eta)) + chain_arg(f, complex(xr, -eta), complex(xr, eta)) +
                                   chain_arg(f, complex(xr, eta), complex(xl, eta)) + chain_arg(f, complex(xl, eta), complex(xl, -eta))) / (2 * arb.pi())
                        n = int(round(float(tot)))
                        if (tot - n).abs_upper() < 1e-3:
                            rect = {'eta': eta, 'count': n, 'enclosure': tot.str(10), 'unique': (n == len(zeros))}
                            break
                    except RuntimeError as e:
                        log('    rectangle eta = %s failed: %s' % (eta, e))
            row['rectangle'] = rect
            log('  lambda = %s, %s: %d zero(s): %s ; poles in brackets excluded: %s ; rectangle %s' % (
                lname, key, len(zeros), ', '.join(nst(z['y'], 14) for z in zeros), all(z['no_pole_in_bracket'] for z in zeros),
                ('count %d, eta %s' % (rect['count'], rect['eta'])) if rect else 'not run'))
            # values near the four poles, on both sides
            poles = {'left_pole': pl, 'right_pole': pr}
            near = {}
            with prec(128):
                for pname, pole in poles.items():
                    for d in (1e-9, 1e-7, 1e-5, 2.4e-5):
                        for side, sign in (('inside', 1 if pname == 'left_pole' else -1), ('outside', -1 if pname == 'left_pole' else 1)):
                            y = (pole + sign * arb(d)).mid()
                            da, db = mix.pair(y)
                            near['%s_%s_%.1e' % (pname, side, d)] = mix.M(lam, da, db).str(8, radius=False)
            row['values_near_poles'] = near
            out.setdefault('lambda=' + lname, {})[key] = row
            for z in zeros:
                ck = claims_zero.get((lname, key))
                if ck is not None and len(zeros) >= 1:
                    claim('lambda=%s %s zero' % (lname, key), ck, z['y'], 9)
        if key == 'I2':
            pass
    # summary of the sizes at the poles
    vals_by_delta = {}
    for lk, d in out.items():
        if not lk.startswith('lambda='):
            continue
        for ik, row in d.items():
            for name, v in row['values_near_poles'].items():
                delta = name.rsplit('_', 1)[1]
                vals_by_delta.setdefault(delta, []).append((abs(float(v)), lk, ik, name))
    summ = {}
    for delta, lst in sorted(vals_by_delta.items()):
        mn = min(lst)
        mx = max(lst)
        summ[delta] = {'min_abs': '%.4g' % mn[0], 'min_at': '%s %s %s' % mn[1:], 'max_abs': '%.4g' % mx[0]}
        log('  |M| at distance %s from a pole (all lambda, both intervals, both sides): smallest %.4g at %s, largest %.4g' % (delta, mn[0], ' '.join(mn[1:3]) + ' ' + mn[3], mx[0]))
    out['sizes_at_poles'] = summ
    out['reality_check_max_imaginary_part'] = '%.2e' % mix.imax


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--parts', default='zeros,central,mixture')
    ap.add_argument('--lambdas', default=None)
    ap.add_argument('--no-rect', dest='rect', action='store_false')
    ap.add_argument('--T', type=float, default=12.0)
    ap.add_argument('--h', type=float, default=None)
    ap.add_argument('--nuni', type=int, default=None)
    args = ap.parse_args()
    q = args.quick
    args.h = args.h if args.h is not None else (1 / 16 if q else 1 / 64)
    args.nuni = args.nuni if args.nuni is not None else (60 if q else 400)
    args.lambdas = args.lambdas or ('0.5' if q else '0.01,0.1,0.5,0.9,0.99')
    if q:
        args.rect = False
    out = args.out or ('/tmp/gl_mixture_quick.json' if q else 'gl_mixture_mod5.json')
    LOGF = open(os.path.splitext(out)[0] + '.log', 'a')
    log('# gl_mixture_mod5.py %s' % ' '.join(sys.argv[1:]))
    log('# started %s ; python-flint arb, mpmath %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), mp.__version__))
    S = {}
    if os.path.exists(out):
        try:
            S = json.load(open(out))
        except Exception:
            S = {}
    S['script'] = 'gl_mixture_mod5.py'
    S['mode'] = 'quick' if q else 'full'
    S['args'] = {k: v for k, v in vars(args).items()}
    check_tables()
    for name, fn in (('zeros', part_zeros), ('central', part_central), ('mixture', part_mixture)):
        if name in args.parts.split(','):
            t0 = time.time()
            fn(S, args)
            S.setdefault('parts_done', {})[name] = {'seconds': round(time.time() - t0, 1), 'mode': S['mode']}
            json.dump(S, open(out, 'w'), indent=1)
    log('# done, wrote %s ; total %.1f s' % (out, time.time() - START))


if __name__ == '__main__':
    main()
