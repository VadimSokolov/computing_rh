#!/usr/bin/env python3
# yo_explicit.py: the explicit formula numbers of Chapter 17 (Yoshida's Hermitian Forms, Weil Positivity and the Thorin Program), group Y1 of the number inventory
# (misc/reviews/nick-yoshida/review-A.md, Section 5), together with yo_tent.py. Added for the book in October 2026, written from the mathematics (the collaborator's scripts were not read).
# Checks (chapter lines in the source ch_yoshida.tex). Line 56: kappa_1 = sum over gamma > 0 of gamma^(minus 2) = 0.0231049931 as half of (log xi)'' at 1/2 (Taylor coefficients of log xi at 1/2 from the polygamma and zeta derivatives), and kappa_2 = sum gamma^(minus 4)
# (used at line 266) from the fourth derivative, both against the certified zeros with a smooth tail and against the book's sums of a_k^(minus 1) and a_k^(minus 2) (Chapter 18, p_1 and p_2, a_k = gamma_k^2). Line 80: Theta(0.003) from the prime side of (yo:EF) at F = h_x and from
# the zeros, and V_1(h_0.003) by both forms of (yo:V1) with Pf(1/|y|)(h_x) = h_x(0)(log 4x minus gamma_E). Line 139: V_1(Delta_1) by the closed form of Proposition yo:prop:tent (c) and by the residue series (derived here by closing the contour in the t integral of (yo:V1)).
# Table yo:tab:D: prime side by Proposition yo:prop:tent (c), zero side with the first 1500 and the first 1700 certified zeros, the mean tail 2(log(T/2 pi) + 1)/(pi T), the tail integrated against the Riemann von Mangoldt density, and Landau's correction 2 Lambda(n)/(pi sqrt(n) T) at x = log n;
# V_1(Delta_x) both ways. Lines 149 and 402: <chi_0, chi_0> = D(log 2)/log 2. Line 199: the entries (0,0), (0,1), (1,2), (2,minus 1), (3,3) of the prime side matrix at a = 1 (undamped form, alpha = 1/2), computed from F = chi_n * conj chi_m by quadrature on (yo:EF) in the variable y,
# against Theorem yo:thm:lowner (i) from the zeros; the difference should be the tail (2/a) sum over gamma > T of gamma^(minus 2). Table yo:tab:W: prime side by Proposition yo:prop:W (ii), zero side over both signs, and eq:ch12:Garith of the book at t = x/2 coded literally. Line 454: xi'/xi(1).
# The inequality of line 266, 4 kappa_1 c minus 4 kappa_2 c^3 <= 2 xi'/xi(alpha) <= 4 kappa_1 c with c = alpha minus 1/2, is tested at alpha = 0.55, 0.6, 0.75, 1.
# Method. mpmath at 30 digits (20 with --quick). The conventions are those of the chapter: gamma_rho = (rho minus 1/2)/i, Fhat(t) = integral of F(y) exp(i t y) dy, h_x(y) = (4 pi x)^(minus 1/2) exp(minus y^2/4x), Delta_x = (x minus |y|)_+. The two forms of V_1 are
# (1/2 pi) times the integral of Fhat(t) Re psi(1/4 + i t/2) over the real line, and the real variable form with the finite part Pf(1/|y|)(F) = lim of [integral over |y| > eps of F/|y| plus 2 F(0) log eps] and g(y) = e^(y/2)/(2 sinh y) minus 1/(2y). The residue series is
# V_1(Delta_x) = x psi(1/4) + psi'(1/4)/2 minus 2 sum over k >= 0 of exp(minus x(2k + 1/2))/(2k + 1/2)^2. Zeros: data/zeros_1700.txt, the book's certified ordinates (68 digits, one per line, all on the line); the count used is printed with every table. The zero tails use the smooth count
# N(T) = theta(T)/pi + 1 with the tail starting at the point T* where N(T*) = N + 1/2 (and, for comparison with the chapter, at gamma_N itself).
# Error control. Every identity is computed twice (two forms of V_1, prime side against zero side, closed form against residue series), and the difference is printed. Quadrature is Gauss Legendre on pieces on which the integrand is analytic, with extra working digits where g(y) has a removable singularity.
# The summation over zeros is exact for the zeros used; what is left is the tail, which is estimated and never hidden: the tables show the truncated sum, the sum with each tail and the residual against the prime side.
# Output: data/yo_explicit.json (values as strings with the digits trusted, the parameters, the zero counts, the precision, the chapter values and a status per claim) and data/yo_explicit.log (the printed report).
# Usage: python3 yo_explicit.py [--quick] [--out FILE] [--zeros FILE] [--dps D]; --quick runs a reduced case (three rows of Table yo:tab:D, two entries of the matrix, two rows of Table yo:tab:W, dps 20) in about 2 seconds in one process and writes to the temporary folder unless --out is given.
# The full case is the default and takes about 16 seconds in one process (the data file of the book was made by a local run within the 30 second limit). On the cluster: bash misc/tools/hopper_run.sh -c 1 -t 10 -g yo_explicit.json -g yo_explicit.log code/yo_explicit.py data/zeros_1700.txt -- --out yo_explicit.json
import argparse
import json
import sys
import tempfile
import time
from decimal import Decimal
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LOG = []
CLAIMS = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def fs(v, d=20):
    return mp.nstr(v, d)


def digits_agree(a, b):
    if a == b:
        return 99
    r = abs(a - b) / max(abs(b), mp.mpf(10) ** -30)
    return int(mp.floor(-mp.log10(r)))


def claim(label, line, chapter, mine, kind='printed'):
    """chapter is the printed string; status says whether mine rounds to it."""
    d = Decimal(chapter)
    dec = -d.as_tuple().exponent
    diff = mine - mp.mpf(chapter)
    unit = mp.mpf(10) ** (-dec)
    if abs(diff) <= unit * mp.mpf('0.5000001'):
        status, short = 'agrees to the printed digits', 'agrees'
    elif abs(diff) <= unit:
        status, short = 'agrees to within one unit of the last printed digit (truncated, not rounded)', 'AGREES TO ONE UNIT (truncated)'
    else:
        status, short = 'differs', 'DIFFERS'
    CLAIMS.append(dict(claim=label, line=line, chapter=chapter, mine=fs(mine, 17), difference=fs(diff, 4), status=status))
    log('CLAIM %-44s line %-4s chapter %-20s mine %-24s diff %-11s %s' % (label, line, chapter, fs(mine, 17), fs(diff, 3), short))
    return status == 'agrees to the printed digits'


def within(val, printed):
    """True if val rounds to the printed decimal string."""
    dec = -Decimal(printed).as_tuple().exponent
    return bool(abs(val - mp.mpf(printed)) <= mp.mpf('0.5000001') * mp.mpf(10) ** (-dec))


# ---------------------------------------------------------------- data
def find_zeros(path):
    cands = ([Path(path)] if path else []) + [Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), ROOT / 'data' / 'zeros_1700.txt', HERE / 'zeros_1700.txt']
    for c in cands:
        if c.exists():
            return c
    sys.exit('zeros_1700.txt not found: pass --zeros')


def prime_powers(nmax):
    nmax = int(nmax)
    s = bytearray([1]) * (nmax + 1)
    s[0] = 0
    if nmax >= 1:
        s[1] = 0
    for i in range(2, int(nmax ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, nmax + 1, i)))
    out = []
    for p in range(2, nmax + 1):
        if s[p]:
            q = p
            while q <= nmax:
                out.append((q, p))
                q *= p
    out.sort()
    return out


def tail1(T):
    return (mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T)


def tail2(T):
    return (3 * mp.log(T / (2 * mp.pi)) + 1) / (18 * mp.pi * T ** 3)


def tail_start(N, gN):
    return mp.findroot(lambda T: mp.siegeltheta(T) / mp.pi + 1 - (N + mp.mpf(1) / 2), gN + mp.mpf(1) / 2)


def D_rvm_tail(x, T):
    c = 2 * mp.pi
    f = mp.log(T / c) / T ** 2
    fp = (1 - 2 * mp.log(T / c)) / T ** 3
    return 4 * tail1(T) + (2 / mp.pi) * (mp.sin(x * T) * f / x + mp.cos(x * T) * fp / x ** 2)


# ---------------------------------------------------------------- pieces of the explicit formula
BREAKS = [mp.mpf(v) for v in ('0', '0.5', '1', '2', '3', '4', '6', '8', '12', '16', '24', '32', '48', '64', '96', '128', '192', '256', '384', '512')]


def pieces(a, b):
    pts = [a, b]
    for g in BREAKS:
        for s in (g, -g):
            if a < s < b:
                pts.append(s)
    return sorted(set(pts))


def psi_re(t):
    return mp.re(mp.digamma(mp.mpc(mp.mpf(1) / 4, t / 2)))


def qd(f, a, b):
    return mp.quad(f, pieces(a, b), method='gauss-legendre', maxdegree=9)


def gfun(y):
    with mp.workdps(mp.mp.dps + 30):
        y = mp.mpf(y)
        if y == 0:
            return mp.mpf(1) / 4
        return mp.exp(y / 2) / (2 * mp.sinh(y)) - 1 / (2 * y)


def lamsum(nmax_log, fn):
    """sum over prime powers k >= 2 with log k <= nmax_log of Lambda(k)/sqrt(k) * fn(k, log k)."""
    s = mp.mpf(0)
    for k, p in prime_powers(int(mp.exp(nmax_log)) + 1):
        lk = mp.log(k)
        if lk <= nmax_log:
            s += mp.log(p) / mp.sqrt(k) * fn(k, lk)
    return s


def gauss_cut(x, margin=5):
    return mp.sqrt(4 * x * (mp.mp.dps + margin) * mp.log(10))


def V1_gauss_t(x):
    """(1/2 pi) int e^{-t^2 x} Re psi(1/4 + i t/2) dt, over the real line."""
    L = mp.sqrt((mp.mp.dps + 5) * mp.log(10) / x)
    return 2 * qd(lambda t: mp.exp(-t * t * x) * psi_re(t), mp.mpf(0), L) / (2 * mp.pi)


def V1_gauss_y(x):
    """real variable form for F = h_x."""
    h0 = 1 / mp.sqrt(4 * mp.pi * x)
    Y = gauss_cut(x)
    pf = h0 * (mp.log(4 * x) - mp.euler)
    gint = 2 * qd(lambda y: gfun(y) * h0 * mp.exp(-y * y / (4 * x)), mp.mpf(0), Y)
    return -pf / 2 - (mp.euler + mp.log(2)) * h0 - gint, pf, gint


def V1_tent_closed(x):
    x = mp.mpf(x)
    return x - x * mp.log(x) - (mp.euler + mp.log(2)) * x - 2 * mp.quad(lambda y: gfun(y) * (x - y), [0, x])


def V1_tent_series(x):
    x = mp.mpf(x)
    s = mp.mpf(0)
    k = 0
    eps = mp.mpf(10) ** (-(mp.mp.dps + 8))
    while True:
        c = 2 * k + mp.mpf(1) / 2
        term = mp.exp(-x * c) / c ** 2
        s += term
        if term < eps and k > 3:
            break
        k += 1
    return x * mp.digamma(mp.mpf(1) / 4) + mp.psi(1, mp.mpf(1) / 4) / 2 - 2 * s


def tent_prime_side(x):
    x = mp.mpf(x)
    pr = lamsum(x, lambda k, lk: (x - lk) if lk < x else mp.mpf(0))
    v1 = V1_tent_closed(x)
    return 32 * mp.sinh(x / 4) ** 2 - x * mp.log(mp.pi) - 2 * pr + v1, v1


def zero_D(x, zs):
    return sum(8 * mp.sin(g * x / 2) ** 2 / g ** 2 for g in zs)


# ---------------------------------------------------------------- sections
def sec_kappa(zs, out):
    log('\n== kappa_1, kappa_2 (lines 56 and 266) ==')
    s0 = mp.mpf(1) / 2
    q14 = mp.mpf(1) / 4
    c = [-mp.zeta(s0, 1, k) / mp.factorial(k) for k in range(5)]
    l = [mp.log(c[0])]
    for k in range(1, 5):
        l.append(c[k] / c[0] - sum(j * l[j] * c[k - j] for j in range(1, k)) / (k * c[0]))
    a1 = -mp.log(mp.pi) / 2 + mp.digamma(q14) / 2 + l[1]
    a3 = mp.psi(2, q14) / 48 + l[3]
    a2 = -4 + mp.psi(1, q14) / 8 + l[2]
    a4 = -8 + mp.psi(3, q14) / 384 + l[4]
    k1, k2 = a2, -2 * a4
    N = len(zs)
    T = tail_start(N, zs[-1])
    z1 = sum(1 / g ** 2 for g in zs)
    z2 = sum(1 / g ** 4 for g in zs)
    log('Taylor coefficients of log xi(1/2 + z): a_1 = %s, a_3 = %s (both must vanish by the functional equation)' % (fs(a1, 3), fs(a3, 3)))
    log('kappa_1 = a_2 = %s' % fs(k1, 25))
    log('kappa_2 = minus 2 a_4 = %s' % fs(k2, 25))
    log('zeros: N = %d, gamma_N = %s, tail start T* = %s' % (N, fs(zs[-1], 12), fs(T, 12)))
    log('sum_{n<=N} gamma^-2 = %s, plus smooth tail %s = %s, difference to kappa_1 = %s' % (fs(z1, 15), fs(tail1(T), 8), fs(z1 + tail1(T), 15), fs(z1 + tail1(T) - k1, 3)))
    log('sum_{n<=N} gamma^-4 = %s, plus smooth tail %s = %s, difference to kappa_2 = %s' % (fs(z2, 15), fs(tail2(T), 8), fs(z2 + tail2(T), 15), fs(z2 + tail2(T) - k2, 3)))
    book1, book2 = mp.mpf('0.0231049931154189707889'), mp.mpf('0.0000371725992852696862')
    log('book (Chapter 18, p_1 and p_2): %s and %s; differences %s and %s' % (fs(book1, 22), fs(book2, 22), fs(k1 - book1, 3), fs(k2 - book2, 3)))
    claim('kappa_1 = (1/2)(log xi)\'\'(1/2)', '56', '0.0231049931', k1)
    out['kappa'] = dict(kappa1=fs(k1, 28), kappa2=fs(k2, 28), a1_check=fs(a1, 5), a3_check=fs(a3, 5), zeros=N, zero_sum1=fs(z1, 20), zero_sum2=fs(z2, 20),
                        tail_start=fs(T, 15), tail1=fs(tail1(T), 12), tail2=fs(tail2(T), 12), zero_plus_tail1=fs(z1 + tail1(T), 20), zero_plus_tail2=fs(z2 + tail2(T), 20),
                        book_p1='0.0231049931154189707889', book_p2='0.0000371725992852696862', diff_book1=fs(k1 - book1, 4), diff_book2=fs(k2 - book2, 4))
    return k1, k2


def dlogxi(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s)


def sec_xi1(k1, k2, out):
    log('\n== xi\'/xi(1) (line 454) and the inequality of line 266 ==')
    closed = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2

    def lx(s):
        return mp.log(mp.mpf(1) / 2 * s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s))
    h = mp.mpf(10) ** -4
    num = (-lx(1 + 2 * h) + 8 * lx(1 + h) - 8 * lx(1 - h) + lx(1 - 2 * h)) / (12 * h)
    log('xi\'/xi(1) closed form %s; numerical derivative of log xi (step 1e-4, fourth order) %s; difference %s' % (fs(closed, 20), fs(num, 20), fs(num - closed, 3)))
    claim('xi\'/xi(1) = 1 + gamma_E/2 minus (1/2) log 4 pi', '454', '0.0230957', closed)
    rows = []
    for al in ('0.55', '0.6', '0.75', '1'):
        c = mp.mpf(al) - mp.mpf(1) / 2
        v = closed if al == '1' else dlogxi(mp.mpf(al))
        lo, hi = 4 * k1 * c - 4 * k2 * c ** 3, 4 * k1 * c
        ok = lo <= 2 * v <= hi
        log('alpha = %-4s 2 xi\'/xi = %s, lower bound %s, upper bound %s, inside: %s' % (al, fs(2 * v, 10), fs(lo, 10), fs(hi, 10), ok))
        rows.append(dict(alpha=al, two_xi_prime_over_xi=fs(2 * v, 15), lower=fs(lo, 15), upper=fs(hi, 15), inside=bool(ok)))
    out['xi_prime_1'] = dict(closed=fs(closed, 25), numerical=fs(num, 25), inequality_line_266=rows)


def sec_theta(zs, out):
    log('\n== Theta(0.003) and V_1(h_x) (line 80) ==')
    x = mp.mpf('0.003')
    h0 = 1 / mp.sqrt(4 * mp.pi * x)
    Y = gauss_cut(x)
    v1t = V1_gauss_t(x)
    v1y, pf, gint = V1_gauss_y(x)
    log('V_1(h_x) from the t integral: %s' % fs(v1t, 24))
    log('V_1(h_x) from the real variable form: %s (Pf = %s, integral of g h_x = %s)' % (fs(v1y, 24), fs(pf, 15), fs(gint, 15)))
    log('difference of the two forms: %s' % fs(v1t - v1y, 3))
    claim('V_1(h_0.003), t integral form', '80', '6.08466798944601', v1t)
    claim('V_1(h_0.003), real variable form', '80', '6.08466798944601', v1y)
    pole = 2 * mp.exp(x / 4)
    pr = 2 * lamsum(Y, lambda k, lk: h0 * mp.exp(-lk * lk / (4 * x)))
    arch = -mp.log(mp.pi) * h0
    Tfun = pole + arch - pr + v1t
    theta_p = Tfun / 2
    theta_z = sum(mp.exp(-g * g * x) for g in zs)
    log('prime side of (yo:EF): pole %s, minus log(pi) h_x(0) = %s, prime sum %s, V_1 %s; T(h_x) = %s' % (fs(pole, 15), fs(arch, 15), fs(pr, 5), fs(v1t, 15), fs(Tfun, 25)))
    log('Theta(0.003) prime side %s, zero side (%d zeros) %s, difference %s' % (fs(theta_p, 25), len(zs), fs(theta_z, 25), fs(theta_p - theta_z, 3)))
    claim('Theta(0.003), prime side', '80', '1.09522012145284', theta_p)
    claim('Theta(0.003), zero side', '80', '1.09522012145284', theta_z)
    out['theta'] = dict(x='0.003', V1_t_form=fs(v1t, 28), V1_y_form=fs(v1y, 28), V1_difference=fs(v1t - v1y, 4), Pf=fs(pf, 20), g_integral=fs(gint, 20), prime_term=fs(pr, 8),
                        Theta_prime_side=fs(theta_p, 28), Theta_zero_side=fs(theta_z, 28), difference=fs(theta_p - theta_z, 4), zeros_used=len(zs))


def sec_tent_v1(out):
    log('\n== V_1(Delta_1) (line 139) ==')
    a = V1_tent_closed(1)
    b = V1_tent_series(1)
    log('closed form of Proposition yo:prop:tent (c): %s' % fs(a, 28))
    log('residue series:                              %s' % fs(b, 28))
    log('difference %s (agreement to %d digits)' % (fs(a - b, 3), digits_agree(a, b)))
    claim('V_1(Delta_1), closed form', '139', '-0.50847597187659', a)
    claim('V_1(Delta_1), residue series', '139', '-0.50847597187659', b)
    out['v1_tent_1'] = dict(closed=fs(a, 30), residue_series=fs(b, 30), difference=fs(a - b, 4), digits_agree=digits_agree(a, b))


ROWS_D = [('0.2', '0.2', '0.1089996', '0.1089981', '0.25787887'), ('0.5', '0.5', '0.0803252', '0.0803270', '0.15008051'), ('log2', 'log2', '0.1271112', '0.1269547', '-0.04998528'),
          ('1.0', '1', '0.0880146', '0.0880147', '-0.50847597'), ('1.5', '1.5', '0.0898189', '0.0898186', '-1.52909231'), ('log5', 'log5', '0.1343264', '0.1340964', '-1.78866457'),
          ('2.0', '2', '0.1066822', '0.1066850', '-2.80144646'), ('3.0', '3', '0.0794186', '0.0794176', '-5.86891443'), ('4.0', '4', '0.0702819', '0.0702855', '-9.39384635'),
          ('6.0', '6', '0.0962464', '0.0962437', '-17.16435327')]


def sec_table_D(zs, out, quick):
    log('\n== Table yo:tab:D and <chi_0, chi_0> (lines 149, 402, 410 to 419) ==')
    rows = ROWS_D if not quick else [ROWS_D[0], ROWS_D[2], ROWS_D[6]]
    cuts = {1500: zs[:1500], len(zs): zs}
    info = {}
    for N, z in cuts.items():
        T = tail_start(N, z[-1])
        info[N] = dict(T=T, Tg=z[-1], mean=4 * tail1(T), mean_g=4 * tail1(z[-1]))
        log('zero count %d: gamma_N = %s, T* = %s, mean tail 2(log(T/2pi)+1)/(pi T) = %s (at T*) and %s (at gamma_N)' % (N, fs(z[-1], 10), fs(T, 10), fs(info[N]['mean'], 8), fs(info[N]['mean_g'], 8)))
    out_rows = []
    log('%-6s %-12s %-12s | %-12s %-12s %-12s %-12s | residuals prime minus (zeros + mean tail at T*) | V1 closed | V1 series' % ('x', 'prime', 'chap prime', 'zeros1500', 'z1500+mean', 'z1500+rvm', 'z1700+mean'))
    for lab, xs, cp, cz, cv in rows:
        x = mp.log(int(lab[3:])) if lab.startswith('log') else mp.mpf(xs)
        pside, v1c = tent_prime_side(x)
        v1s = V1_tent_series(x)
        r = dict(x=lab, x_value=fs(x, 15), prime_side=fs(pside, 18), V1_closed=fs(v1c, 18), V1_series=fs(v1s, 18), V1_difference=fs(v1c - v1s, 4), chapter=dict(prime=cp, zero=cz, V1=cv))
        r['prime_agrees_chapter'] = bool(abs(pside - mp.mpf(cp)) <= mp.mpf('0.5000001') * mp.mpf(10) ** (-(-Decimal(cp).as_tuple().exponent)))
        r['V1_agrees_chapter'] = bool(abs(v1c - mp.mpf(cv)) <= mp.mpf('0.5000001') * mp.mpf(10) ** (-(-Decimal(cv).as_tuple().exponent)))
        landau = {}
        for N, z in cuts.items():
            zsum = zero_D(x, z)
            T = info[N]['T']
            mean = info[N]['mean']
            rvm = D_rvm_tail(x, T)
            corner = 0
            if lab.startswith('log'):
                n = int(lab[3:])
                p = [pp for kk, pp in prime_powers(n) if kk == n][0]
                corner = 2 * mp.log(p) / (mp.pi * mp.sqrt(n) * T)
            landau[N] = corner
            r['zeros_%d' % N] = dict(truncated=fs(zsum, 15), plus_mean_tail=fs(zsum + mean, 15), plus_mean_tail_at_gamma_N=fs(zsum + info[N]['mean_g'], 15), plus_rvm_tail=fs(zsum + rvm, 15),
                                     landau_correction=fs(corner, 8), plus_rvm_and_landau=fs(zsum + rvm + corner, 15), residual_mean=fs(pside - zsum - mean, 4),
                                     residual_mean_at_gamma_N=fs(pside - zsum - info[N]['mean_g'], 4), residual_rvm=fs(pside - zsum - rvm, 4), residual_rvm_landau=fs(pside - zsum - rvm - corner, 4))
        z15, z17 = r['zeros_1500'], r['zeros_%d' % len(zs)]
        rep = mp.mpf(z15['plus_mean_tail_at_gamma_N'])
        r['zero_side_chapter_reproduced'] = dict(first_1500_plus_mean_tail_at_gamma_1500=fs(rep, 10), matches_chapter_zero_column=within(rep, cz), chapter_zero=cz,
                                                 first_1500_plus_mean_tail_at_Tstar=fs(mp.mpf(z15['plus_mean_tail']), 10), matches_with_Tstar=within(mp.mpf(z15['plus_mean_tail']), cz),
                                                 first_1700_plus_mean_tail_at_gamma_N=fs(mp.mpf(z17['plus_mean_tail_at_gamma_N']), 10), matches_with_1700=within(mp.mpf(z17['plus_mean_tail_at_gamma_N']), cz))
        r['chapter_prime_minus_zero'] = fs(mp.mpf(cp) - mp.mpf(cz), 3)
        r['mine_prime_minus_chapter_zero'] = fs(pside - mp.mpf(cz), 3)
        log('%-6s prime %-12s (chapter %s) | 1500 zeros %-11s +mean tail at gamma_N %-11s (chapter zero %s, reproduced: %s) | +mean at T* %-11s +rvm %-11s | z1700+mean(gamma_N) %-11s' % (
            lab, fs(pside, 10), cp, fs(mp.mpf(z15['truncated']), 8), fs(rep, 8), cz, within(rep, cz), fs(mp.mpf(z15['plus_mean_tail']), 8), fs(mp.mpf(z15['plus_rvm_tail']), 8), fs(mp.mpf(z17['plus_mean_tail_at_gamma_N']), 8)))
        log('       prime minus chapter zero column %s (chapter prime minus zero %s); residuals prime minus (zeros + tail): 1500 mean(T*) %s, 1500 rvm %s, 1500 rvm+Landau %s (Landau %s), 1700 rvm %s, 1700 rvm+Landau %s; V1 closed %s series %s (difference %s)' % (
            fs(pside - mp.mpf(cz), 3), r['chapter_prime_minus_zero'], z15['residual_mean'], z15['residual_rvm'], z15['residual_rvm_landau'], z15['landau_correction'], z17['residual_rvm'], z17['residual_rvm_landau'], fs(v1c, 11), fs(v1s, 11), fs(v1c - v1s, 3)))
        out_rows.append(r)
    out['table_D'] = dict(zero_counts=[1500, len(zs)], tails={str(N): dict(T_star=fs(info[N]['T'], 12), gamma_N=fs(info[N]['Tg'], 12), mean_tail=fs(info[N]['mean'], 12), mean_tail_at_gamma_N=fs(info[N]['mean_g'], 12)) for N in cuts}, rows=out_rows)
    # chi_0
    x2 = mp.log(2)
    p2, _ = tent_prime_side(x2)
    c0 = p2 / x2
    log('D(log 2) = %s, <chi_0, chi_0> = D(log 2)/log 2 = %s' % (fs(p2, 15), fs(c0, 15)))
    claim('<chi_0, chi_0> at a = (1/2) log 2', '149', '0.18338', c0)
    out['chi0'] = dict(D_log2=fs(p2, 20), chi0_chi0=fs(c0, 20))


def Fnm(n, m, a, y):
    ph = mp.expj(mp.pi * m * y / a)
    if n == m:
        return ph * (2 * a - y) / (2 * a)
    b = mp.pi * (n - m) / a
    return ph * (mp.expj(b * a) - mp.expj(b * (y - a))) / (2 * a * mp.mpc(0, 1) * b)


def matrix_prime(n, m, a):
    a = mp.mpf(a)
    L = 2 * a
    F0 = 1 if n == m else 0
    wp = mp.mp.dps + 12

    def reF(y):
        with mp.workdps(wp):
            return mp.re(Fnm(n, m, a, y))
    pole = 4 * mp.quad(lambda y: reF(y) * mp.cosh(y / 2), [0, L])
    prm = lamsum(L, lambda k, lk: 2 * reF(lk))
    pf = mp.quad(lambda y: (2 * reF(y) - 2 * F0) / y, [0, L]) + 2 * F0 * mp.log(L)
    gF = mp.quad(lambda y: gfun(y) * 2 * reF(y), [0, L])
    v1 = -pf / 2 - (mp.euler + mp.log(2)) * F0 - gF
    return pole - mp.log(mp.pi) * F0 - prm + v1, dict(pole=pole, prime_sum=prm, V1=v1)


def matrix_zero(n, m, a, zs):
    a = mp.mpf(a)
    nn, nm = mp.pi * n / a, mp.pi * m / a
    s = mp.mpf(0)
    for g in zs:
        s += mp.sin(a * g) ** 2 * (1 / ((g + nn) * (g + nm)) + 1 / ((-g + nn) * (-g + nm)))
    return (-1) ** (n + m) * (2 / a) * s


def sec_matrix(zs, out, quick):
    log('\n== entries of the prime side matrix at a = 1, alpha = 1/2 against Theorem yo:thm:lowner (i) (line 199) ==')
    a = mp.mpf(1)
    entries = [(0, 0), (0, 1), (1, 2), (2, -1), (3, 3)]
    if quick:
        entries = entries[:2]
    res = []
    for N in (1500, len(zs)):
        T = tail_start(N, zs[N - 1])
        log('zero count %d, T* = %s, predicted tail magnitude (2/a) sum_{gamma > T*} gamma^-2 = %s (at gamma_N: %s)' % (N, fs(T, 10), fs(2 / a * tail1(T), 8), fs(2 / a * tail1(zs[N - 1]), 8)))
    for (n, m) in entries:
        pv, parts = matrix_prime(n, m, a)
        r = dict(n=n, m=m, prime_side=fs(pv, 18), pole=fs(parts['pole'], 15), prime_sum=fs(parts['prime_sum'], 15), V1=fs(parts['V1'], 15))
        for N in (1500, len(zs)):
            zv = matrix_zero(n, m, a, zs[:N])
            T = tail_start(N, zs[N - 1])
            pred = (-1) ** (n + m) * (2 / a) * tail1(T)
            r['zeros_%d' % N] = dict(zero_side=fs(zv, 15), difference=fs(pv - zv, 6), predicted_tail=fs(pred, 6), residual=fs(pv - zv - pred, 3))
        z15, z17 = r['zeros_1500'], r['zeros_%d' % len(zs)]
        log('(n,m) = (%d,%d): prime side %s; zero side N=1500 %s, difference %s (predicted %s, residual %s); N=%d %s, difference %s (predicted %s, residual %s)' % (
            n, m, fs(pv, 12), z15['zero_side'][:12], z15['difference'], z15['predicted_tail'], z15['residual'], len(zs), z17['zero_side'][:12], z17['difference'], z17['predicted_tail'], z17['residual']))
        res.append(r)
    # the (0,0) entry must be D(1)/... : <chi_0,chi_0> at a = 1 is D(2)/2
    d2, _ = tent_prime_side(2)
    log('check: entry (0,0) at a = 1 equals D(2)/2 = %s' % fs(d2 / 2, 15))
    out['matrix'] = dict(a='1', alpha='1/2', entries=res, check_00_equals_D2_over_2=fs(d2 / 2, 18))


ROWS_W = [('0.01', '0', '0.299396022507559', '0.299396022507559'), ('0.01', '14.134725', '2.0336517', '2.0336504'), ('0.05', '17', '1.149233992297324', '1.149233992297324'),
          ('0.2', '14.134725', '1.000075829748193', '1.000075829748193'), ('0.2', '17', '0.232948895021655', '0.232948895021655'),
          ('0.5', '21.022040', '1.000350786528200', '1.000350786528200'), ('0.5', '30', '0.927167615181088', '0.927167615181088')]


def W_arch_integral(T, x, lo=None, hi=None):
    L = mp.sqrt((mp.mp.dps + 5) * mp.log(10) / x)
    a, b = (T - L, T + L) if lo is None else (lo, hi)
    return qd(lambda t: mp.exp(-(t - T) ** 2 * x) * psi_re(t), a, b)


def W_prime(T, x):
    T, x = mp.mpf(T), mp.mpf(x)
    polar = 2 * mp.exp(-(T * T - mp.mpf(1) / 4) * x) * mp.cos(T * x)
    lg = -mp.log(mp.pi) / mp.sqrt(4 * mp.pi * x)
    arch = W_arch_integral(T, x) / (2 * mp.pi)
    Y = gauss_cut(x)
    pr = lamsum(Y, lambda k, lk: mp.cos(T * lk) * mp.exp(-lk * lk / (4 * x))) / mp.sqrt(mp.pi * x)
    return polar + lg + arch - pr, dict(polar=polar, lg=lg, arch=arch, prime=pr)


def W_garith(u, t):
    """eq:ch12:Garith literally: G(u,t) + G(-u,t) = 2 e^{t/2 - 2 t u^2} cos(2 t u) + A(u,t) - P(u,t)."""
    u, t = mp.mpf(u), mp.mpf(t)
    first = 2 * mp.exp(t / 2 - 2 * t * u * u) * mp.cos(2 * t * u)
    L = mp.sqrt((mp.mp.dps + 5) * mp.log(10) / (2 * t))

    def m(r):
        return psi_re(r) - mp.log(mp.pi)
    A = (qd(lambda r: mp.exp(-2 * t * (r - u) ** 2) * m(r), u - L, u + L) + qd(lambda r: mp.exp(-2 * t * (r + u) ** 2) * m(r), -u - L, -u + L)) / (4 * mp.pi)
    Y = mp.sqrt(8 * t * (mp.mp.dps + 5) * mp.log(10))
    P = lamsum(Y, lambda k, lk: mp.exp(-lk * lk / (8 * t)) * mp.cos(u * lk)) / mp.sqrt(2 * mp.pi * t)
    return first + A - P


def W_zero(T, x, zs):
    T, x = mp.mpf(T), mp.mpf(x)
    s = mp.mpf(0)
    lim = (mp.mp.dps + 6) * mp.log(10)
    for g in zs:
        if (g - T) ** 2 * x > lim + 50 and g > T:
            break
        s += mp.exp(-(g - T) ** 2 * x) + mp.exp(-(g + T) ** 2 * x)
    return s


def sec_table_W(zs, out, quick):
    log('\n== Table yo:tab:W (lines 432 to 438) and eq:ch12:Garith at t = x/2 ==')
    rows = ROWS_W if not quick else [ROWS_W[3], ROWS_W[6]]
    res = []
    for xs, Ts, cp, cz in rows:
        pv, parts = W_prime(Ts, xs)
        zv = W_zero(Ts, xs, zs)
        gv = W_garith(Ts, mp.mpf(xs) / 2)
        r = dict(x=xs, T=Ts, prime_side=fs(pv, 20), zero_side=fs(zv, 20), difference=fs(pv - zv, 4), garith=fs(gv, 20), garith_minus_prime=fs(gv - pv, 4), digits_agree=digits_agree(pv, zv), chapter=dict(prime=cp, zero=cz))
        for key, val in (('prime', pv), ('zero', zv)):
            cv = r['chapter'][key]
            dec = -Decimal(cv).as_tuple().exponent
            r['chapter_%s_agrees' % key] = bool(abs(val - mp.mpf(cv)) <= mp.mpf('0.5000001') * mp.mpf(10) ** (-dec))
            r['chapter_%s_diff' % key] = fs(val - mp.mpf(cv), 4)
        res.append(r)
        log('x = %-5s T = %-10s prime %s zero %s | difference %s (%d digits) | Garith minus prime %s | chapter prime %s (mine minus it %s), chapter zero %s (mine minus it %s)' % (
            xs, Ts, fs(pv, 17), fs(zv, 17), fs(pv - zv, 3), digits_agree(pv, zv), fs(gv - pv, 3), cp, r['chapter_prime_diff'], cz, r['chapter_zero_diff']))
    # the seven digit row: T = 14.134725 is gamma_1 rounded; the effect of the rounding and of cutting the t integral
    if not quick:
        xs, Ts = '0.01', '14.134725'
        T, x = mp.mpf(Ts), mp.mpf(xs)
        T1 = zs[0]
        shift = W_zero(T1, xs, zs) - W_zero(T, xs, zs)
        pv = mp.mpf([r for r in res if r['x'] == xs and r['T'] == Ts][0]['prime_side'])
        cut = {}
        for Lc in (30, 35, 40, 50, 60):
            full = W_arch_integral(T, x) / (2 * mp.pi)
            part = W_arch_integral(T, x, lo=-mp.mpf(Lc), hi=mp.mpf(Lc)) / (2 * mp.pi)
            cut[str(Lc)] = fs(part - full, 4)
        log('row x = 0.01, T = 14.134725: zero side at T = gamma_1 minus zero side at T = 14.134725 is %s; the archimedean integral cut to |t| <= L changes the prime side by %s' % (fs(shift, 3), cut))
        trap = {}
        full = W_arch_integral(T, x) / (2 * mp.pi)
        for hs in ('1', '0.5', '0.4', '0.3', '0.25', '0.2'):
            hq = mp.mpf(hs)
            K = int(110 / hq)
            trap[hs] = fs(hq * sum(mp.exp(-(k * hq - T) ** 2 * x) * psi_re(k * hq) for k in range(-K, K + 1)) / (2 * mp.pi) - full, 4)
        log('the archimedean integral of this row by the trapezoid rule of step h on [-110, 110]: error of the prime side %s (the chapter row is high by 1.29e-6)' % trap)
        out['table_W_seven_digit_row'] = dict(shift_from_rounding_T=fs(shift, 4), effect_of_cutting_t_integral_at_L=cut, trapezoid_rule_error_by_step=trap,
                                              chapter_prime_minus_exact='1.29e-6', verdict='the chapter prime side is in error; prime and zero sides agree to 30 digits at 2.03365040977979')
    out['table_W'] = dict(zero_count=len(zs), rows=res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--dps', type=int, default=None)
    args = ap.parse_args()
    dps = args.dps or (20 if args.quick else 30)
    mp.mp.dps = dps
    t0 = time.time()
    zp = find_zeros(args.zeros)
    with open(zp) as fh:
        zs = [mp.mpf(s.strip()) for s in fh if s.strip()]
    outp = Path(args.out) if args.out else (Path(tempfile.gettempdir()) / 'yo_explicit_quick.json' if args.quick else ROOT / 'data' / 'yo_explicit.json')
    log('yo_explicit.py: dps %d, %d certified zeros from %s (first %s, last %s), quick = %s' % (dps, len(zs), zp, fs(zs[0], 15), fs(zs[-1], 15), args.quick))
    out = dict(script='code/yo_explicit.py', dps=dps, quick=bool(args.quick), zeros_file=str(zp.name), zeros_loaded=len(zs), zeros_first=fs(zs[0], 20), zeros_last=fs(zs[-1], 20))
    k1, k2 = sec_kappa(zs, out)
    sec_xi1(k1, k2, out)
    sec_theta(zs, out)
    sec_tent_v1(out)
    sec_table_D(zs, out, args.quick)
    sec_matrix(zs, out, args.quick)
    sec_table_W(zs, out, args.quick)
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
