# maxwald_weak.py: a rigorous bound for the Laplace distance between the reciprocal Xi law and its explicit GGC approximant, for
# Section sec:maxwald (Proposition mw:prop:approx, the chapter's Theorem thm:weak (iv), where the chapter claims 3.2e-10).
# Added for the book in October 2026 from the checks of the authors' chapter, misc/reviews/nick-maxwald/ (review-2 MINOR C:362 and
# the recomputation list; the bound is written from the mathematics, not from the chapter's verify scripts). Notation: rho runs over
# the zeros with gamma = Im rho > 0, delta = Re rho - 1/2 (|delta| < 1/2), a_rho = -(rho - 1/2)^2 = gamma^2 - delta^2 - 2 i delta gamma,
# so Re a_rho = gamma^2 - delta^2 > 0 and |a_rho| = gamma^2 + delta^2 >= gamma^2; then 1/G(lambda) = xi(1/2)/xi(1/2 + sqrt(lambda)) =
# prod_rho (1 + lambda/a_rho)^(-1) (Hadamard, genus 0 in lambda), and 1/G_H is the product over gamma <= H, all on the line
# (Platt and Trudgian 2021, H = 3e12), where a_rho = gamma^2. For lambda >= 0 each factor F = (1 + lambda/a)^(-1) has |1 + lambda/a| >=
# 1 + lambda Re(a)/|a|^2 >= 1, hence |F| <= 1 and |F - 1| = |F| lambda/|a| <= lambda/gamma^2, and |prod F - 1| <= sum |F - 1|. So
#   D_H = sup_{lambda >= 0} |1/G - 1/G_H| <= sup_lambda lambda E_H/G_H(lambda) <= E_H S_K,   E_H = sum_{gamma > H} gamma^(-2),
#   S_K = sup_{lambda >= 0} lambda/G_K(lambda), G_K = prod_{k <= K} (1 + lambda/gamma_k^2), any K <= N(H) (more factors only decrease
# 1/G_H, and any K distinct zeros on the line below H will do, so no count of zeros is needed); K = 1 gives gamma_1^2 E_H, the chapter's bound.
# (1) E_H. Let N_0(T) = (T/2 pi) log(T/(2 pi e)) and E(T) = c1 log T + c2 log log T + c3 with (c1, c2, c3) = (0.1038, 0.2573, 9.3675); HSW22
# (Hasanalizade, Shen, Wong, arXiv:2107.06506, Corollary 1.2) states |N(T) - N_0(T)| <= E(T) for T >= e, N(T) = #{0 < gamma <= T}.
# Partial summation, E_H = 2 int_H^oo (N(t) - N(H)) t^(-3) dt <= 2 int_H^oo (N_0(t) - N_0(H) + E(t) + E(H)) t^(-3) dt, gives
#   E_H <= M + R,   M = (log(H/2 pi) + 1)/(2 pi H),   R = (2 E(H) + c1/2 + c2/(2 log H))/H^2
# (int_H^oo t^(-3) log log t dt <= log log H/(2 H^2) + 1/(4 H^2 log H) by concavity of log). M is the chapter's 1.4797e-12; R is below 3e-24.
# (2) S_K. The logarithmic derivative of lambda/G_K is -phi(lambda)/lambda with phi(lambda) = sum_k lambda/(lambda + gamma_k^2) - 1, which
# increases from -1 to K - 1, so lambda/G_K increases up to the unique root lambda* of phi and decreases after it. The script finds the
# root in floating point, certifies phi(lo) < 0 < phi(hi) for lo = lambda*(1 - 1e-9), hi = lambda*(1 + 1e-9) in ball arithmetic (python-flint,
# Arb, 128 bits), and bounds the supremum by hi/G_K(lo) >= lambda*/G_K(lambda*). The zeros are the balls of acb.zeta_zeros (Arb isolates each
# zero rigorously; the script checks that every real part is exactly 1/2 and that the ordinates increase), so no table is trusted; the
# book's table zeros_1700.txt is compared if it sits in the folder. Outside the certificate the script prints floating point values.
# Run on Hopper (a minute, one core), from a folder that may hold zeros_1700.txt: python3 maxwald_weak.py [KMAX]; writes maxwald_weak.json
# there (default KMAX = 1700; the bound for each K in a short list is reported).
import json
import math
import os
import sys
import time

import numpy as np
from flint import arb, acb, ctx

PREC = 128
ctx.prec = PREC
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1700
K_LIST = sorted(set([k for k in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 1700) if k <= KMAX] + [KMAX]))
CHAPTER_BOUND = 3.2e-10
CHAPTER_EPS = 1.5e-12
t0 = time.time()


def up(x):
    """A float that is at least the upper end of the ball x (the factor 1 + 1e-12 covers the rounding of float)."""
    return float((x.upper() * arb('1.000000000001')).upper())


def say(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- 0. the two lemmas on a factor, tested on random a with Re a > 0
rng = np.random.default_rng(1)
worst = 0.0
for _ in range(20000):
    a = complex(rng.uniform(0.01, 50) * 10 ** rng.uniform(-2, 6), rng.uniform(-1, 1) * 10 ** rng.uniform(-2, 6))
    lam = 10 ** rng.uniform(-6, 8)
    F = 1 / (1 + lam / a)
    assert abs(F) <= 1 + 1e-12 and abs(F - 1) <= lam / abs(a) + 1e-14, (a, lam)   # absolute slack: F - 1 cancels in double precision
    if lam / abs(a) > 1e-6:
        worst = max(worst, abs(F - 1) / (lam / abs(a)))
say('factor lemma: |F| <= 1 and |F - 1| <= lambda/|a| on 20000 random (a, lambda), worst ratio %.6f' % worst)

# ---------------------------------------------------------------- 1. the zeros (Arb)
Z = acb.zeta_zeros(1, KMAX)
say('zeros: %d balls in %.1f s, last ordinate %s' % (len(Z), time.time() - t0, Z[-1].imag.str(15, radius=False)))
half = arb(1) / 2
assert all(z.real == half for z in Z), 'a real part is not exactly 1/2'
gam = [z.imag for z in Z]
assert all(gam[k] < gam[k + 1] for k in range(KMAX - 1)), 'ordinates not increasing: the balls overlap'
H = arb(3) * arb(10) ** 12
assert all(g < H for g in gam), 'a zero above H'
assert gam[0] > 14 and gam[0] < 14.2
maxrad = max(float(g.rad()) for g in gam)
say('largest ball radius of an ordinate: %.3g' % maxrad)
table_ok = None
if os.path.exists('zeros_1700.txt'):
    tab = [arb(l.strip()) for l in open('zeros_1700.txt') if l.strip()]
    table_ok = all(gam[k].overlaps(tab[k]) for k in range(min(KMAX, len(tab))))
    say('book table zeros_1700.txt agrees with the Arb balls for the first %d zeros: %s' % (min(KMAX, len(tab)), table_ok))
    assert table_ok

# ---------------------------------------------------------------- 2. E_H from HSW22
c1, c2, c3 = arb('0.1038'), arb('0.2573'), arb('9.3675')
assert H > arb(1).exp()
L = H.log()
LL = L.log()
E_H_err = c1 * L + c2 * LL + c3                                  # E(H)
M = ((H / (2 * arb.pi())).log() + 1) / (2 * arb.pi() * H)         # the main term, the chapter's epsilon_H
R = (2 * E_H_err + c1 / 2 + c2 / (2 * L)) / H ** 2                # the contribution of the error term of N(T)
EH = M + R
g1sq = gam[0] * gam[0]
crude = g1sq * EH
say('E(H) = %.6f, M = %.6e, R = %.3e, E_H <= %.8e (chapter: %.1e)' % (float(E_H_err), up(M), up(R), up(EH), CHAPTER_EPS))
say('crude bound gamma_1^2 E_H <= %.6e against the chapter %.1e: %s' % (up(crude), CHAPTER_BOUND, 'holds' if up(crude) < CHAPTER_BOUND else 'FAILS'))

# ---------------------------------------------------------------- 3. S_K and the bound E_H S_K
rows = []
for K in K_LIST:
    g2 = [g * g for g in gam[:K]]
    if K == 1:
        S_plus = g2[0]                      # lambda/(1 + lambda/g2) increases to g2 and never reaches it
        row = dict(K=1, lambda_star=None, lambda_lo=None, lambda_hi=None, S_upper=up(S_plus), S_float=float(S_plus.mid()))
    else:
        g2f = [float(q.mid()) for q in g2]
        lo, hi = 1e-9, 1e12
        for _ in range(300):
            mid = math.sqrt(lo * hi)
            if math.fsum(mid / (mid + q) for q in g2f) < 1:
                lo = mid
            else:
                hi = mid
        ls = math.sqrt(lo * hi)
        a_lo, a_hi = arb(ls * (1 - 1e-9)), arb(ls * (1 + 1e-9))        # exact balls of two doubles
        phi_lo = sum(a_lo / (a_lo + q) for q in g2) - 1
        phi_hi = sum(a_hi / (a_hi + q) for q in g2) - 1
        assert phi_lo < 0 and phi_hi > 0, 'bracket of the root not certified'
        prod = arb(1)
        for q in g2:
            prod = prod * (1 + a_lo / q)
        S_plus = a_hi / prod
        pf = math.fsum(math.log1p(ls / q) for q in g2f)
        row = dict(K=K, lambda_star=ls, lambda_lo=float(a_lo), lambda_hi=float(a_hi), S_upper=up(S_plus), S_float=ls * math.exp(-pf))
    row['bound_D_upper'] = up(EH * S_plus)
    row['bound_over_chapter'] = row['bound_D_upper'] / CHAPTER_BOUND
    rows.append(row)
    say('K = %4d: lambda* = %s, sup lambda/G_K <= %.8g, D_H <= %.6e (%.4f of 3.2e-10)' % (
        K, 'inf' if row['lambda_star'] is None else '%.6f' % row['lambda_star'], row['S_upper'], row['bound_D_upper'], row['bound_over_chapter']))
best = rows[-1]

# ---------------------------------------------------------------- 4. information outside the certificate
# If RH also holds above H, the distance is D = sup_lambda P_H (1 - Q) with P_H = prod_{gamma <= H} (1 + lambda/gamma^2)^(-1) and
# Q = prod_{gamma > H} (1 + lambda/gamma^2)^(-1) = exp(-lambda E_H (1 + o(1))) for lambda << H^2. The model below takes the first K
# zeros exactly, the zeros K < k <= N(H) through exp(-lambda A) with A = sum_{K < k, gamma <= H} gamma^(-2) from the smooth density
# (log(T/2 pi) + 1)/(2 pi T) (the fluctuation S(T) moves A by about 1e-7), and E_H by its main term M. It shows how far the certified
# bound is from the true distance, not a bound itself.
Hf = 3e12
gK = float(gam[KMAX - 1].mid())
g2f = np.array([float(g.mid()) ** 2 for g in gam])
smooth = lambda T: (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T)
A_tail = smooth(gK) - smooth(Hf)
lams = best['lambda_star'] * np.exp(np.linspace(-1.5, 1.5, 3001)) if best['lambda_star'] else np.array([1e6])
logP = -np.log1p(lams[:, None] / g2f[None, :]).sum(axis=1) - lams * A_tail
model = np.exp(logP) * (-np.expm1(-lams * smooth(Hf)))
info = dict(note='floating point, not part of the certificate',
            RH_model_distance=float(model.max()), RH_model_argmax_lambda=float(lams[model.argmax()]),
            tail_sum_A_between_gamma_K_and_H=A_tail, bound_over_RH_model=best['bound_D_upper'] / float(model.max()),
            chapter_bound_over_RH_model=CHAPTER_BOUND / float(model.max()))
say('model of the true distance if RH holds above H too (floating point): %.4e at lambda = %.3f; the certified bound is %.4f times it, the chapter bound %.2f times it' % (
    info['RH_model_distance'], info['RH_model_argmax_lambda'], info['bound_over_RH_model'], info['chapter_bound_over_RH_model']))

out = dict(
    script='maxwald_weak.py', working_precision_bits=PREC, zeros_from='acb.zeta_zeros (Arb), K up to %d, real part exactly 1/2' % KMAX,
    max_ordinate_ball_radius=maxrad, table_zeros_1700_agrees=table_ok,
    H=3e12, N_bound='HSW22 Corollary 1.2: |N(T) - (T/2pi) log(T/(2 pi e))| <= 0.1038 log T + 0.2573 log log T + 9.3675, T >= e',
    E_H_main_term_M=up(M), E_H_error_term_R=up(R), E_H_upper=up(EH), E_H_chapter_value=CHAPTER_EPS,
    gamma1_squared=float(g1sq.mid()), crude_bound_gamma1sq_EH_upper=up(crude), chapter_bound=CHAPTER_BOUND,
    crude_bound_below_chapter_bound=bool(up(crude) < CHAPTER_BOUND),
    by_K=rows, best_K=best['K'], best_bound_D_upper=best['bound_D_upper'],
    best_bound_below_chapter_bound=bool(best['bound_D_upper'] < CHAPTER_BOUND),
    information=info, seconds=round(time.time() - t0, 1))
name = os.environ.get('MAXWALD_OUT', 'maxwald_weak.json')
json.dump(out, open(name, 'w'), indent=1)
say('wrote %s in %.1f s' % (name, time.time() - t0))
