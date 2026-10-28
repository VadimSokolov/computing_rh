# code/xg_kernel.py (October 2026, from misc/reviews/nick-xiggc/: review-A1.md, the rows prop:harm, ex:harm, ex:g2, thm:notid and thm:nogamma, items (d)4 to (d)6): the numbers
# of Section sec:xg:kernel (the end of Section 7): the clock as a harmonic mean of Brownian exit time transforms (Proposition xg:prop:harm), the two examples
# xg:ex:harm and xg:ex:g2, and the two bounds behind Theorem xg:thm:nogamma. Written from the mathematics, not from the authors' scripts.
# Setting. Phi(x) = 2 sum_n (2 pi^2 n^4 e^{9x/2} - 3 pi n^2 e^{5x/2}) exp(-pi n^2 e^{2x}) is the Polya kernel of eq:ch2:Phi (even; every term is positive for x >= 0),
# p = Phi/xi(1/2) is the density of the Polya variable X, and G(s) = E cosh(sqrt(s) X) = xi(1/2 + sqrt s)/xi(1/2) (eq:ch2:mgf); the clock T has E exp(-sT) = 1/G(s).
# With L_x(s) = 1/cosh(x sqrt s) = E exp(-x^2 s C/2), C the exit time of Brownian motion from (-1,1), 1/G(s) = 1/E[1/L_X(s)], and by Jensen's inequality 1/G(s) <= E L_X(s).
# (1) Quadrature on [0, 4] with the terms n <= 12 of Phi (the n = 13 term of the integral below, and the integrand at x = 4, are printed as the truncation evidence):
#   check: the total mass of p, and E cosh(wX) against xi(1/2 + w)/xi(1/2) at w = 1, 3, 10;
#   number 1: 1/G(100) = E exp(-100 T), from the quadrature and from xi;  number 2: E L_X(100) = E 1/cosh(10X), by quadrature.
#   Both are repeated in the variable u = x/2 of the authors (U = |X|/2 has the density 8 Phi_c(u)/xi(1/2) on (0, oo), with Phi_c(u) = sum_n (2 pi^2 n^4 e^{9u}
#   - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}) = Phi(2u)/2, and 2U sqrt(s) = |X| sqrt(s)), on [0, 2].
# (2) Example xg:ex:harm: the harmonic mean of (1+s)^-2 and (1+2s)^-2, its poles, its inverse Laplace transform 2 exp(-3t/5) sin(t/5) (exact, sympy; and the
#   Laplace transform of this function by quadrature at s = 1, 2) and the sign of the function (negative exactly on (5 pi, 10 pi) in its first period).
# (3) Example xg:ex:g2: the zeros of the numerator (1 + lam s)^2 + (1 + s)^2 of E (1 + sY)^-2, Y = 1 or lam with probability 1/2: exact for symbolic lam, and the
#   transform itself at the zeros for lam = 2, 1/2, 3, 10 (the poles -1 and -1/lam are real).
# (4) Theorem xg:thm:nogamma: the bound 1/G(s) <= c_2 exp(-x_0 sqrt s) with c_2 = 2/P(|X| >= x_0), checked at x_0 = 1, 2 and s = 1 to 1e6 (it follows from
#   G(s) >= cosh(x_0 sqrt s) P(|X| >= x_0)).
# (5) gamma_1 and 2 gamma_1 from zeros_1700.txt, and Xi(gamma_1) = 2 int_0^oo Phi(x) cos(gamma_1 x) dx = 0 (the characteristic function of X vanishes at gamma_1, that
#   of X/2 at 2 gamma_1 = 28.27).
# Run from a folder that holds zeros_1700.txt: python3 xg_kernel.py (seconds, one core); writes xg_kernel.json there.
import json
import os
import time

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
NMAX = 12
t0 = time.time()


def say(*a):
    print(*a, flush=True)


def find(name):
    for p in (name, 'data/' + name, '../data/' + name):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)


def S(x, k=15):
    return mp.nstr(x, k)


Z = [mp.mpf(s) for s in open(find('zeros_1700.txt')).read().split()]
assert len(Z) == 1700

HALF = mp.mpf(1) / 2


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def term(n, x):                         # the n-th term of the Polya kernel Phi(x) (book normalisation), x >= 0
    return 2 * (2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * x / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * x / 2)) * mp.exp(-mp.pi * n ** 2 * mp.exp(2 * x))


def Phi(x):
    return mp.fsum(term(n, x) for n in range(1, NMAX + 1))


def Phic(u):                            # the authors' kernel in the variable u = x/2: Phi_c(u) = Phi(2u)/2
    return mp.fsum((2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u) - 3 * mp.pi * n ** 2 * mp.exp(5 * u)) * mp.exp(-mp.pi * n ** 2 * mp.exp(4 * u)) for n in range(1, NMAX + 1))


XI0 = xi(HALF)
PANELS = [mp.mpf(v) for v in ('0', '0.25', '0.5', '1', '2', '4')]
PANELS_U = [v / 2 for v in PANELS]


def mean(h):                            # E h(X) for an even h, X ~ p, by quadrature on [0, 4]
    return 2 * mp.quad(lambda x: Phi(x) * h(x), PANELS) / XI0


def mean_u(h):                          # E h(2U), the same in the variable u of the authors, on [0, 2]
    return mp.quad(lambda u: 8 * Phic(u) / XI0 * h(2 * u), PANELS_U)


say('Section sec:xg:kernel: the clock as a harmonic mean of Brownian exit time transforms, two examples, and the Gamma mixture bound (%d digits, %d terms of Phi)' % (mp.mp.dps, NMAX))
say('xi(1/2) = %s' % S(XI0, 20))

# ---------------------------------------------------------------- 1. the harmonic mean representation
mass = mean(lambda x: 1)
mass_u = mean_u(lambda x: 1)
say('(1) total mass of p on R: %s (quadrature in x), %s (in u of the authors); 1 - mass = %s' % (S(mass, 20), S(mass_u, 20), S(1 - mass, 3)))
checks = []
for w in (1, 3, 10):
    q = mean(lambda x, w=w: mp.cosh(w * x))
    ex = xi(HALF + w) / XI0
    checks.append(dict(w=w, quadrature=S(q, 25), closed_form=S(ex, 25), relative_difference=float(abs(q / ex - 1))))
    say('    E cosh(%d X): quadrature %s, xi(1/2+w)/xi(1/2) %s, relative difference %s' % (w, S(q, 20), S(ex, 20), S(abs(q / ex - 1), 3)))

s0 = 100
r0 = mp.sqrt(s0)
G_q = mean(lambda x: mp.cosh(r0 * x))
G_xi = xi(HALF + r0) / XI0
L_q = mean(lambda x: 1 / mp.cosh(r0 * x))
G_u = mean_u(lambda x: mp.cosh(r0 * x))
L_u = mean_u(lambda x: 1 / mp.cosh(r0 * x))
tail13 = 2 * mp.quad(lambda x: term(13, x) * mp.cosh(r0 * x), PANELS) / XI0
end4 = 2 * Phi(4) * mp.cosh(r0 * 4) / XI0
say('    s = 100: G(100) = E cosh(10 X): quadrature %s, xi(21/2)/xi(1/2) %s, u of the authors %s' % (S(G_q, 20), S(G_xi, 20), S(G_u, 20)))
say('    number 1: E exp(-100 T) = 1/G(100) = %s (quadrature), %s (xi), %s (u of the authors)' % (S(1 / G_q, 15), S(1 / G_xi, 15), S(1 / G_u, 15)))
say('    number 2: E L_X(100) = E 1/cosh(10 X) = %s (quadrature), %s (u of the authors); Jensen: %s <= %s: %s' % (S(L_q, 15), S(L_u, 15), S(1 / G_q, 6), S(L_q, 6), 1 / G_q <= L_q))
say('    truncation: the n = 13 term of the integral of G(100) is %s, the integrand at x = 4 is %s (both relative to G(100) = %s)' % (S(tail13, 3), S(end4, 3), S(G_q, 8)))
say('    harmonic mean over arithmetic mean at s = 100: %s' % S(1 / G_q / L_q, 6))
harmonic = dict(total_mass_x=S(mass, 30), total_mass_u=S(mass_u, 30), checks=checks, s=s0,
                G100_quadrature=S(G_q, 25), G100_xi=S(G_xi, 25), G100_u=S(G_u, 25),
                clock_transform_100=S(1 / G_q, 25), clock_transform_100_float=float(1 / G_q),
                LX_100=S(L_q, 25), LX_100_u=S(L_u, 25), LX_100_float=float(L_q),
                jensen_holds=bool(1 / G_q <= L_q), tail_n13=S(tail13, 5), integrand_at_x4=S(end4, 5))

# ---------------------------------------------------------------- 2. Example xg:ex:harm
z = sp.Symbol('z')
t = sp.Symbol('t', positive=True)
f1 = (1 + z) ** -2
f2 = (1 + 2 * z) ** -2
Hm = sp.simplify(2 / (1 / f1 + 1 / f2))
assert sp.simplify(Hm - 1 / (1 + 3 * z + sp.Rational(5, 2) * z ** 2)) == 0
poles = sp.solve(sp.Eq(5 * z ** 2 + 6 * z + 2, 0), z)
assert set(poles) == {sp.Rational(-3, 5) + sp.I / 5, sp.Rational(-3, 5) - sp.I / 5}
inverse = 2 * sp.exp(-3 * t / 5) * sp.sin(t / 5)
assert sp.simplify(sp.laplace_transform(inverse, t, z, noconds=True) - Hm) == 0
say('(2) Example xg:ex:harm: harmonic mean = %s, poles %s, inverse transform %s (its Laplace transform equals the harmonic mean exactly)' % (Hm, poles, inverse))
finv = lambda tt: 2 * mp.exp(-mp.mpf(3) / 5 * tt) * mp.sin(tt / 5)
lap = []
for sv in (1, 2):
    lq = mp.quad(lambda tt: finv(tt) * mp.exp(-sv * tt), [0, 10 * mp.pi, 30 * mp.pi, mp.inf])
    hv = 1 / (1 + 3 * mp.mpf(sv) + mp.mpf(5) / 2 * sv ** 2)
    lap.append(dict(s=sv, quadrature=S(lq, 20), harmonic_mean=S(hv, 20), difference=float(abs(lq - hv))))
    say('    Laplace transform of the inverse at s = %d: quadrature %s, harmonic mean %s' % (sv, S(lq, 18), S(hv, 18)))
tmin = 5 * (mp.pi + mp.atan(mp.mpf(1) / 3))
say('    f(5 pi) = %s, f(10 pi) = %s; f(2.5 pi) = %s > 0, f(7.5 pi) = %s < 0, f(12.5 pi) = %s > 0; the minimum on (5 pi, 10 pi) is %s at t = %s' % (
    S(finv(5 * mp.pi), 3), S(finv(10 * mp.pi), 3), S(finv(2.5 * mp.pi), 5), S(finv(7.5 * mp.pi), 5), S(finv(12.5 * mp.pi), 5), S(finv(tmin), 6), S(tmin, 6)))
assert finv(2.5 * mp.pi) > 0 and finv(7.5 * mp.pi) < 0 and finv(12.5 * mp.pi) > 0 and finv(tmin) < 0 and 5 * mp.pi < tmin < 10 * mp.pi
ex_harm = dict(harmonic_mean=str(Hm), poles=[str(p) for p in poles], inverse=str(inverse), laplace_checks=lap,
               negative_on=['5 pi', '10 pi'], minimum=S(finv(tmin), 8), minimum_at=S(tmin, 8))

# ---------------------------------------------------------------- 3. Example xg:ex:g2
lam = sp.Symbol('lambda', real=True)
num = sp.expand((1 + lam * z) ** 2 + (1 + z) ** 2)
rp = (-(lam + 1) + sp.I * (lam - 1)) / (lam ** 2 + 1)
rm = (-(lam + 1) - sp.I * (lam - 1)) / (lam ** 2 + 1)
assert sp.simplify(num.subs(z, rp)) == 0 and sp.simplify(num.subs(z, rm)) == 0
say('(3) Example xg:ex:g2: numerator %s vanishes at s = (-(lam+1) +- i (lam-1))/(lam^2+1), whose imaginary part +-(lam-1)/(lam^2+1) is not zero for real lam != 1' % num)
g2rows = []
for lv in (2, mp.mpf(1) / 2, 3, 10):
    lv = mp.mpmathify(lv)
    r = mp.mpc(-(lv + 1) / (lv ** 2 + 1), (lv - 1) / (lv ** 2 + 1))
    tv = ((1 + lv * r) ** -2 + (1 + r) ** -2) / 2
    g2rows.append(dict(lam=S(lv, 6), zero=S(r, 12), transform_at_zero=float(abs(tv)), poles=[-1.0, float(-1 / lv)]))
    say('    lam = %s: zero %s, |E(1+sY)^-2| there %s, poles -1 and %s' % (S(lv, 6), S(r, 12), S(abs(tv), 3), S(-1 / lv, 6)))
ex_g2 = dict(numerator=str(num), zeros='(-(lam+1) +- i (lam-1))/(lam^2+1)', rows=g2rows)

# ---------------------------------------------------------------- 4. Theorem xg:thm:nogamma
rows4 = []
for x0 in (1, 2):
    P = 2 * mp.quad(lambda x: Phi(x) / XI0, [mp.mpf(x0)] + [q for q in PANELS if q > x0])
    c2 = 2 / P
    ratios = []
    for sv in (1, 10, 100, 10 ** 4, 10 ** 6):
        Ginv = XI0 / xi(HALF + mp.sqrt(sv))
        ratios.append(Ginv / (c2 * mp.exp(-x0 * mp.sqrt(sv))))
    rows4.append(dict(x0=x0, P_abs_X_ge_x0=S(P, 15), c2=S(c2, 10), ratios=[S(q, 4) for q in ratios], max_ratio=float(max(ratios))))
    say('(4) x_0 = %d: P(|X| >= x_0) = %s, c_2 = 2/P = %s; (1/G(s))/(c_2 exp(-x_0 sqrt s)) at s = 1, 10, 100, 1e4, 1e6: %s; all <= 1: %s' % (
        x0, S(P, 10), S(c2, 8), [S(q, 3) for q in ratios], all(q <= 1 for q in ratios)))
    assert all(q <= 1 for q in ratios)

# ---------------------------------------------------------------- 5. the first zero
g1 = Z[0]
Xi_g1 = 2 * mp.quad(lambda x: Phi(x) * mp.cos(g1 * x), [mp.mpf(v) for v in ('0', '0.125', '0.25', '0.375', '0.5', '0.75', '1', '1.5', '2', '3', '4')])
say('(5) gamma_1 = %s, 2 gamma_1 = %s; Xi(gamma_1) = 2 int_0^4 Phi(x) cos(gamma_1 x) dx = %s (Xi(0) = %s)' % (S(g1, 20), S(2 * g1, 20), S(Xi_g1, 3), S(XI0, 10)))
assert abs(Xi_g1) < mp.mpf(10) ** -30

json.dump(dict(script='xg_kernel.py', zeros='zeros_1700.txt', digits=mp.mp.dps, terms_of_Phi=NMAX, harmonic_mean=harmonic, example_harm=ex_harm, example_g2=ex_g2,
               nogamma=rows4, gamma_1=S(g1, 25), two_gamma_1=S(2 * g1, 25), Xi_at_gamma_1=S(Xi_g1, 5), seconds=round(time.time() - t0, 1)),
          open('xg_kernel.json', 'w'), indent=1)
say('wrote xg_kernel.json in %.1f s' % (time.time() - t0))
