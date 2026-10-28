# code/xg_defect.py (October 2026, from misc/reviews/nick-xiggc/: review-A1.md, the row thm:sin2, item (d)8 and error E7): the rate of the
# compound Poisson defect of Theorem xg:thm:defect (Section sec:xg:defect, the end of Section 5.1), in ball arithmetic. The bound is written from
# the mathematics, not from the verify scripts of the authors' chapter. Notation of Section 5.1: a zero rho = 1/2 + delta + i gamma, gamma > 0, has
# the exponent a = -(rho - 1/2)^2 = gamma^2 - delta^2 - 2 i delta gamma; an off line quadruple 1/2 +- delta +- i gamma, 0 < delta < 1/2, gives the
# pair {a, abar} with A = Re a = gamma^2 - delta^2 and B = |Im a| = 2 delta gamma, and the pair contributes 2 e^{-A t} cos(B t) = 2 e^{-A t} -
# 4 e^{-A t} sin^2(B t / 2) to the heat trace W. So W = (Laplace transform of U_+) - R with R(t) = 4 sum_pairs e^{-A t} sin^2(B t / 2) >= 0, and
#   lambda = int_0^oo R(t) dt / t = sum_pairs log(1 + B^2/A^2)
# by Tonelli and int_0^oo e^{-A t} sin^2(B t / 2) dt / t = (1/4) log(1 + B^2/A^2) (the real part of Frullani's integral); lambda is the total
# Levy mass of the defect Y, which is therefore compound Poisson, and Pr(Y = 0) = e^{-lambda}.
# (1) The two integrals of the proof, the Frullani type integral and -log E e^{-uY} = int_0^oo (1 - e^{-ut}) R(t) dt / t for a single pair, which
# must equal log[(1 + u/A)^2 / ((1 + u/a)(1 + u/abar))], are checked by quadrature (mpmath, 40 digits) for three pairs; a check of the
# identities, not a certificate.
# (2) The bound. The zeros of height at most H = 3e12 lie on the line (Platt and Trudgian 2021), so every pair has gamma > H. For fixed gamma,
# B^2/A^2 = 4 delta^2 gamma^2/(gamma^2 - delta^2)^2 increases with delta, so for delta < 1/2 it is below gamma^2/(gamma^2 - 1/4)^2 = C_gamma
# gamma^(-2) with C_gamma = (1 - 1/(4 gamma^2))^(-2), which decreases in gamma; hence B^2/A^2 <= C gamma^(-2) with C = C_H, and log(1 + x) <= x.
# A pair uses two zeros 1/2 +- delta + i gamma of the same ordinate, so lambda <= (C/2) S_H with S_H = sum over all zeros with gamma > H of
# gamma^(-2), with multiplicity (the worst case: every zero above H off the line at distance almost 1/2). Partial summation gives
# S_H = 2 int_H^oo (N(t) - N(H)) t^(-3) dt. With N_0(T) = (T/2 pi) log(T/(2 pi e)) and E(T) = c1 log T + c2 log log T + c3, (c1, c2, c3) =
# (0.1038, 0.2573, 9.3675), HSW22 (Hasanalizade, Shen, Wong, arXiv:2107.06506, Corollary 1.2) states |N(T) - N_0(T)| <= E(T) for T >= e, N(T) =
# #{0 < gamma <= T}. So S_H <= 2 int_H^oo (N_0(t) - N_0(H) + E(t) + E(H)) t^(-3) dt = M + R with M = (log(H/2 pi) + 1)/(2 pi H) exactly and
# R = 2 int_H^oo (E(t) + E(H)) t^(-3) dt <= (2 E(H) + c1/2 + c2/(2 log H))/H^2, using int_H^oo t^(-3) dt = 1/(2 H^2), int_H^oo t^(-3) log t dt =
# (2 log H + 1)/(4 H^2) and, by concavity of log, log log t <= log log H + (log t - log H)/log H, so int_H^oo t^(-3) log log t dt <= log log H/(2 H^2)
# + 1/(4 H^2 log H). Both closed forms are compared with quadrature (mpmath, 40 digits); the certificate is the ball inequality lambda <= (C/2)(M + R)
# < 1e-12. The same sum is bounded in the proof of Proposition mw:prop:approx (1.4798e-12), which is compared.
# (3) A consistency check of the partial summation with the book's zeros: for the heights h in a short list, the exact sum of gamma^(-2) over
# h < gamma <= gamma_8000 (zeros_8000.txt, 8148.19) plus M(gamma_8000) differs from M(h) by at most R(h) + R(gamma_8000), since S_h and S_gamma8000 lie
# within R of M at their heights; this would fail by a large factor if the formula for M had a wrong constant. Not part of the certificate.
# Run on Hopper from a folder that may hold zeros_8000.txt: python3 xg_defect.py (seconds, one core); writes xg_defect.json there.
import json
import os
import random
import time

import mpmath as mp
from flint import arb, ctx

PREC = 256
ctx.prec = PREC
mp.mp.dps = 40
t0 = time.time()


def say(*a):
    print(*a, flush=True)


def up(x):
    """A float that is at least the upper end of the ball x (the factor 1 + 1e-12 covers the rounding of float)."""
    return float((x.upper() * arb('1.000000000001')).upper())


def tomp(x):
    """The upper end of the ball x as an mpmath number (45 digits)."""
    return mp.mpf(x.upper().str(45, radius=False))


# ---------------------------------------------------------------- 1. the identities of the proof, by quadrature
checks = []
for A, B, u in [('3', '2', '0.7'), ('199.7904548', '1.3', '5'), ('0.9', '4.5', '0.2')]:
    A, B, u = mp.mpf(A), mp.mpf(B), mp.mpf(u)
    cut = [0, 1 / A, 10 / A, 100 / A, mp.inf]
    frul = mp.quad(lambda t: mp.exp(-A * t) * mp.sin(B * t / 2) ** 2 / t, cut)
    frul_closed = mp.log(1 + B ** 2 / A ** 2) / 4
    laplace = mp.quad(lambda t: (1 - mp.exp(-u * t)) * 4 * mp.exp(-A * t) * mp.sin(B * t / 2) ** 2 / t, cut)
    a = mp.mpc(A, -B)
    laplace_closed = mp.log((1 + u / A) ** 2 / ((1 + u / a) * (1 + u / mp.conj(a)))).real
    e1, e2 = abs(frul - frul_closed), abs(laplace - laplace_closed)
    assert e1 < mp.mpf(10) ** -25 and e2 < mp.mpf(10) ** -25, (A, B, u, e1, e2)
    checks.append(dict(A=float(A), B=float(B), u=float(u), frullani_error=float(e1), laplace_error=float(e2)))
say('identities of the proof by quadrature: Frullani type integral and the transform of the defect, errors below 1e-25 for %d pairs' % len(checks))

# ---------------------------------------------------------------- 2. the bound S_H <= M + R and lambda <= (C/2) S_H
c1, c2, c3 = arb('0.1038'), arb('0.2573'), arb('9.3675')
pi = arb.pi()


def E_(T):
    L = T.log()
    return c1 * L + c2 * L.log() + c3


def M_(T):
    return ((T / (2 * pi)).log() + 1) / (2 * pi * T)


def R_(T):
    L = T.log()
    return (2 * E_(T) + c1 / 2 + c2 / (2 * L)) / T ** 2


H = arb(3) * arb(10) ** 12
assert H > arb(1).exp()
EH, M, R = E_(H), M_(H), R_(H)
S_up = M + R
C = 1 / (1 - 1 / (4 * H ** 2)) ** 2
lam_up = C / 2 * S_up
threshold = 1 / arb(10) ** 12
assert C < arb('1.01')
assert lam_up < threshold, 'the rate bound does not reach 1e-12'
say('H = 3e12: E(H) = %.6f, M = %.8e, R = %.3e, S_H <= %.8e' % (float(EH), up(M), up(R), up(S_up)))
say('C - 1 = %.3e, lambda <= (C/2) S_H = %.8e; threshold 1e-12, ratio %.5f, margin %.4e' % (
    up(C - 1), up(lam_up), up(lam_up / threshold), float((threshold - lam_up).lower())))
say('Pr(Y = 0) = exp(-lambda) >= 1 - lambda >= %.15f' % float((1 - lam_up).lower()))
book_eps = arb('1.4798e-12')
assert S_up < book_eps, 'larger than the bound for eps_H in the proof of Proposition mw:prop:approx'
chapter_lam = arb('1.01') / 2 * S_up
assert chapter_lam < arb('0.75e-12')
say('the same sum is below the book value 1.4798e-12 of Proposition mw:prop:approx; with the chapter factor 1.01 the rate is %.4e < 0.75e-12' % up(chapter_lam))

# the closed forms against quadrature, and the pair bound B^2/A^2 <= C gamma^(-2) on a random sample
mH = mp.mpf(3) * mp.mpf(10) ** 12
N0 = lambda t: t / (2 * mp.pi) * mp.log(t / (2 * mp.pi * mp.e))
Ef = lambda t: mp.mpf('0.1038') * mp.log(t) + mp.mpf('0.2573') * mp.log(mp.log(t)) + mp.mpf('9.3675')
cut = [1, 10, 10 ** 3, 10 ** 6, 10 ** 12, mp.inf]
M_num = 2 / mH ** 2 * mp.quad(lambda s: (N0(mH * s) - N0(mH)) / s ** 3, cut)
R_num = 2 / mH ** 2 * mp.quad(lambda s: (Ef(mH * s) + Ef(mH)) / s ** 3, cut)
assert abs(M_num - tomp(M)) < mp.mpf(10) ** -25 * M_num, (M_num, M)
assert R_num <= tomp(R), (R_num, R)
say('quadrature: M = %s (closed form %.10e), R = %.4e <= the closed form bound %.4e' % (mp.nstr(M_num, 12), float(M), float(R_num), up(R)))
random.seed(1)
worst = mp.mpf(0)
for _ in range(2000):
    g = mH * (1 + mp.mpf(random.random()) * 10 ** random.randint(0, 6))
    d = mp.mpf(random.random()) * mp.mpf('0.5')
    ratio = (4 * d ** 2 * g ** 2 / (g ** 2 - d ** 2) ** 2) * g ** 2
    worst = max(worst, ratio)
assert worst <= tomp(C)
say('pair bound B^2/A^2 <= C gamma^(-2): largest sampled ratio %s <= C, C - 1 = %.3e' % (mp.nstr(worst, 20), up(C - 1)))

# ---------------------------------------------------------------- 3. consistency of M with the book's zeros
zfile = next((p for p in ('zeros_8000.txt', 'data/zeros_8000.txt', '../data/zeros_8000.txt') if os.path.exists(p)), None)
zero_rows = None
if zfile:
    gam = [mp.mpf(s) for s in open(zfile).read().split()]
    G = gam[-1]
    M_fun = lambda T: (mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T)
    zero_rows = []
    for h in (50, 100, 300, 1000, 3000):
        P = sum((g ** -2 for g in gam if g > h), mp.mpf(0))
        dev = abs(P + M_fun(G) - M_fun(h))
        bound = tomp(R_(arb(h))) + tomp(R_(arb(str(G))))
        assert dev <= bound, (h, dev, bound)
        zero_rows.append(dict(h=h, deviation=float(dev), bound=float(bound), ratio=float(dev / bound)))
        say('zeros_8000.txt: h = %5d  |sum_{h<gamma<=%.2f} gamma^-2 + M(gamma_8000) - M(h)| = %.3e <= R(h) + R(gamma_8000) = %.3e  (ratio %.3f)' % (h, float(G), dev, bound, dev / bound))
else:
    say('zeros_8000.txt not found: the consistency check of M with the zeros was skipped')

out = dict(
    script='xg_defect.py', working_precision_bits=PREC, H=3e12,
    N_bound='HSW22 Corollary 1.2: |N(T) - (T/2pi) log(T/(2 pi e))| <= 0.1038 log T + 0.2573 log log T + 9.3675, T >= e',
    E_H=up(EH), tail_sum_main_term_M=up(M), tail_sum_error_term_R=up(R), tail_sum_upper=up(S_up),
    pair_constant_C_minus_1=up(C - 1), lambda_upper=up(lam_up), threshold=1e-12, lambda_over_threshold=up(lam_up / threshold),
    margin_to_threshold=float((threshold - lam_up).lower()), prob_Y_zero_lower=float((1 - lam_up).lower()),
    book_value_eps_H_mw_prop_approx=1.4798e-12, tail_sum_below_book_value=bool(S_up < book_eps),
    chapter_factor_rate=up(chapter_lam), chapter_rate_bound=0.75e-12,
    quadrature_M=float(M_num), quadrature_R=float(R_num), pair_bound_largest_sampled_ratio=float(worst),
    identity_checks=checks, zeros_consistency=zero_rows, seconds=round(time.time() - t0, 1))
json.dump(out, open('xg_defect.json', 'w'), indent=1)
say('wrote xg_defect.json in %.1f s' % (time.time() - t0))
