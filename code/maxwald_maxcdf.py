# maxwald_maxcdf.py: the maximum M of independent exponentials with rates pi n and pi (2m-1), whose distribution function is theta_4(t),
# and the Khinchin factorisation of its tilted law, for Section sec:maxwald (ch/maxwald.tex): mw:sec:max with mw:prop:sum (Theorem
# thm:max of the authors' chapter), mw:fig:maxcdf, mw:sec:khinchin and mw:thm:khinchin. Added for the book in October 2026 from the checks
# of the authors' chapter (misc/reviews/nick-maxwald/), recomputed from the mathematics (mpmath at 30 digits; nothing is taken from the
# chapter's scripts). With theta_4(t) = sum_n (-1)^n e^{-pi n^2 t} and xi(s) = s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)/2:
#  1. (i)  theta_4(t) = prod_{n>=1}(1 - e^{-pi n t}) prod_{m>=1}(1 - e^{-pi(2m-1)t}) on a grid of t (2000 factors of each kind), the largest
#          difference; theta_4 by its series (t >= 1) or by Poisson summation (t < 1).
#  2. (ii) E e^{-lambda M} = lambda int_0^oo e^{-lambda t} theta_4(t) dt (quadrature) against sqrt(pi lambda)/sinh sqrt(pi lambda) and the
#          Euler product prod_n (1 + lambda/(pi n^2))^{-1}.
#  3. (iii) E M^{s/2} = (s/2) int_0^oo t^{s/2-1} (1 - theta_4(t)) dt (quadrature) against 2 (1 - 2^{1-s}) xi(s)/(s-1), with log 2 at s = 1,
#          at real and complex s; E M = pi/6, E M^{1/2} = log 2, E M^2 = 7 pi^2/180; the Kolmogorov value theta_4(2/pi) = 0.7300003.
#  4. (iv) P(M_1 + M_2 <= u) = int_0^u theta_4(u-t) theta_4'(t) dt against P(Y^2 <= u) = sum_n (1 - 2 pi n^2 u) e^{-pi n^2 u}; the Polya
#          kernel Phi(u) = sum_{n>=1}(2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}) against e^{3u} f_Y(e^{2u})/4 (f_Y by numerical
#          differentiation of the distribution function), Phi even, and int_0^oo Phi(u) cos(zu) du = Xi(z/2)/8, Xi(t) = xi(1/2 + it).
#  5. Theorem thm:khinchin: for sigma > 1, a = 2^{1-sigma}, the moments E[(M~ 4^{-K})^q] with M~ the m^{sigma/2}-tilted M, taken from the
#          quadrature of 2. and 3. (no zeta function), against E[(G/(pi N^2))^q] = pi^{-q} Gamma(1+sigma/2+q) zeta(sigma+2q)/(Gamma(1+sigma/2)
#          zeta(sigma)); and the same ratio of Mellin transforms at complex q = it/2 (a characteristic function check).
#  6. fig/maxwald_maxcdf.pdf: theta_4(t) = P(M <= t) with the products truncated after one and after three factors of each kind.
# Run from a folder holding this script: python3 maxwald_maxcdf.py [--quick]; writes maxwald_maxcdf.json to the folder and the figure to
# fig/ (created if absent). --quick is the laptop test (20 digits, a few points). On Hopper (about 5 minutes, one core):
# bash misc/tools/hopper_run.sh -c 1 -t 60 -g maxwald_maxcdf.json -g fig/maxwald_maxcdf.pdf code/maxwald_maxcdf.py
import json
import os
import sys
import time

import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QUICK = '--quick' in sys.argv
DPS = 20 if QUICK else 30
mp.mp.dps = DPS
pi = mp.pi
HALF = mp.mpf(1) / 2
TINY = mp.mpf(10) ** (-(DPS + 10))
T0 = time.time()
out = {'working_precision_dps': DPS, 'quick': QUICK, 'library': 'mpmath ' + mp.__version__}


def S(x, n=None):
    return mp.nstr(x, n or DPS)


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


# ---------------------------------------------------------------- theta_4, its derivative and survival function
def theta4(t):
    """theta_4(t) = sum_{n in Z} (-1)^n e^{-pi n^2 t}: direct for t >= 1, Poisson summation t^{-1/2} sum_n e^{-pi (n+1/2)^2/t} for t < 1."""
    t = mp.mpf(t)
    if t >= 1:
        s = mp.mpf(1)
        n = 1
        while True:
            term = mp.exp(-pi * n * n * t)
            s += 2 * (-1) ** n * term
            if term < TINY:
                return s
            n += 1
    if pi / (4 * t) > (DPS + 20) * mp.log(10):
        return mp.mpf(0)
    s = mp.mpf(0)
    n = 0
    while True:
        term = mp.exp(-pi * (n + HALF) ** 2 / t)
        s += 2 * term
        if term < TINY * s:
            return s / mp.sqrt(t)
        n += 1


def dtheta4(t):
    """The density of M, theta_4'(t)."""
    t = mp.mpf(t)
    if t >= 1:
        s = mp.mpf(0)
        n = 1
        while True:
            term = n * n * mp.exp(-pi * n * n * t)
            s += (-1) ** n * term
            if term < TINY:
                return -2 * pi * s
            n += 1
    if pi / (4 * t) > (DPS + 20) * mp.log(10):
        return mp.mpf(0)
    s = mp.mpf(0)
    d = mp.mpf(0)
    n = 0
    while True:
        e = mp.exp(-pi * (n + HALF) ** 2 / t)
        s += 2 * e
        d += 2 * e * pi * (n + HALF) ** 2 / t ** 2
        if e < TINY * s:
            break
        n += 1
    return d / mp.sqrt(t) - s / (2 * t ** mp.mpf(1.5))


def sf4(t):
    """P(M > t) = 1 - theta_4(t) without cancellation for t >= 1/4."""
    t = mp.mpf(t)
    if t >= mp.mpf(1) / 4:
        s = mp.mpf(0)
        n = 1
        while True:
            term = 2 * (-1) ** (n + 1) * mp.exp(-pi * n * n * t)
            s += term
            if abs(term) < TINY * abs(s):
                return s
            n += 1
    return 1 - theta4(t)


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def prod_max(t, N):
    t = mp.mpf(t)
    p = mp.mpf(1)
    for n in range(1, N + 1):
        p *= 1 - mp.exp(-pi * n * t)
    for m in range(1, N + 1):
        p *= 1 - mp.exp(-pi * (2 * m - 1) * t)
    return p


# ---------------------------------------------------------------- 1. (i)
NPROD = 300 if QUICK else 2000
tgrid = ['0.05', '0.1', '0.25', '0.5', '1', '2'] if QUICK else ['0.02', '0.05', '0.1', '0.25', '0.5', '1', '2', '4']
rows = []
worst = mp.mpf(0)
for ts in tgrid:
    a, b = theta4(ts), prod_max(ts, NPROD)
    rel = abs(a - b) / a
    worst = max(worst, rel)
    rows.append({'t': ts, 'theta4': S(a), 'product': S(b), 'rel_diff': mp.nstr(rel, 3)})
dcheck = max(abs(dtheta4(t) - mp.diff(theta4, mp.mpf(t))) / abs(dtheta4(t)) for t in ['0.3', '2'])
out['thm_max_i'] = {'factors_of_each_kind': NPROD, 'rows': rows, 'max_rel_diff': mp.nstr(worst, 3),
                    'dtheta4_analytic_vs_numerical_rel_diff': mp.nstr(dcheck, 3)}
log('(i) max relative difference between theta_4 and the product', mp.nstr(worst, 3))

# ---------------------------------------------------------------- 2. (ii)
def euler_product(lam):
    """prod_{n>=1} (1 + lambda/(pi n^2))^{-1}: 100 factors, and the tail sum_{n>100} log(1 + a/n^2) = sum_j (-1)^{j+1} a^j zeta(2j, 101)/j, a = lambda/pi."""
    a = mp.mpf(lam) / pi
    N = 100
    s = mp.fsum(mp.log(1 + a / (n * n)) for n in range(1, N + 1))
    tail = mp.fsum((-1) ** (j + 1) * a ** j * mp.zeta(2 * j, N + 1) / j for j in range(1, 16))
    return mp.exp(-(s + tail))


rows = []
for ls in (['1'] if QUICK else ['0.1', '1', '10', '100']):
    lam = mp.mpf(ls)
    q = mp.quad(lambda t: lam * mp.exp(-lam * t) * theta4(t), [0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 30, mp.inf])
    c = mp.sqrt(pi * lam) / mp.sinh(mp.sqrt(pi * lam))
    e = euler_product(lam)
    rows.append({'lambda': ls, 'quadrature': S(q), 'closed_form': S(c), 'euler_product': S(e),
                 'rel_diff_quadrature': mp.nstr(abs(q - c) / c, 3), 'rel_diff_product': mp.nstr(abs(e - c) / c, 3)})
out['thm_max_ii'] = rows
log('(ii) done', [r['rel_diff_quadrature'] for r in rows])

# ---------------------------------------------------------------- 3. (iii)
def mellin_closed(s):
    """E M^{s/2} = 2 (1 - 2^{1-s}) xi(s)/(s-1), continued at s = 1 by log 2."""
    s = mp.mpmathify(s)
    if abs(s - 1) < mp.mpf(10) ** (-DPS + 4):
        return mp.log(2)
    return 2 * (1 - 2 ** (1 - s)) * xi(s) / (s - 1)


def mellin_quad(s):
    """E M^{s/2} = (s/2) int_0^oo t^{s/2-1} P(M > t) dt for Re s > 0."""
    s = mp.mpmathify(s)
    return (s / 2) * mp.quad(lambda t: t ** (s / 2 - 1) * sf4(t), [0, 0.1, 0.5, 2, 8, 30, mp.inf])


rows = []
slist = ['1', '2', '3+2j'] if QUICK else ['0.5', '1', '2', '3', '4', '6', '3+2j', '5-4j']
for ss in slist:
    c, q = mellin_closed(ss), mellin_quad(ss)
    rows.append({'s': ss, 'closed_form': S(c), 'quadrature': S(q), 'rel_diff': mp.nstr(abs(c - q) / abs(c), 3)})
out['thm_max_iii'] = rows
EM, ESQ, EM2 = mellin_closed(2), mellin_closed(1), mellin_closed(4)
out['moments'] = {'E_M': S(EM), 'pi_over_6': S(pi / 6), 'E_sqrt_M': S(ESQ), 'log_2': S(mp.log(2)), 'E_M2': S(EM2),
                  '7_pi2_over_180': S(7 * pi ** 2 / 180)}
kol = theta4(2 / pi)
kol_series = 1 + 2 * mp.nsum(lambda k: (-1) ** int(k) * mp.exp(-2 * k * k), [1, mp.inf])
out['kolmogorov'] = {'theta4_2_over_pi': S(kol), 'series_sum_(-1)^k e^{-2k^2}': S(kol_series),
                     'meaning': 'P(sup|b| <= 1) for a Brownian bridge b, since sup|b| = sqrt(pi M/2) in law'}
log('(iii) done', [r['rel_diff'] for r in rows])

# ---------------------------------------------------------------- 4. (iv)
def cdf_Y2(u):
    """P(Y^2 <= u) = sum_{n in Z} (1 - 2 pi n^2 u) e^{-pi n^2 u} (Biane, Pitman, Yor)."""
    u = mp.mpf(u)
    s = mp.mpf(1)
    n = 1
    while True:
        term = (1 - 2 * pi * n * n * u) * mp.exp(-pi * n * n * u)
        s += 2 * term
        if abs(term) < TINY and n > 3:
            return s
        n += 1


def cdf_sum2(u):
    u = mp.mpf(u)
    return mp.quad(lambda t: theta4(u - t) * dtheta4(t), [0, u / 4, u / 2, 3 * u / 4, u])


rows = []
for us in (['0.5', '2'] if QUICK else ['0.2', '0.5', '1', '2', '4']):
    a, b = cdf_sum2(us), cdf_Y2(us)
    rows.append({'u': us, 'P(M1+M2<=u) quadrature': S(a), 'P(Y^2<=u) series': S(b), 'abs_diff': mp.nstr(abs(a - b), 3)})
out['thm_max_iv_sum'] = rows
mom4 = 2 * mellin_closed(4) + 2 * mellin_closed(2) ** 2
out['thm_max_iv_moments'] = {'E(M1+M2)': S(2 * EM), '2 xi(2)': S(2 * xi(2)), 'E(M1+M2)^2': S(mom4), '2 xi(4)': S(2 * xi(4))}


def Phi(u):
    u = mp.mpf(u)
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    s = mp.mpf(0)
    n = 1
    while True:
        term = (2 * pi ** 2 * n ** 4 * e9 - 3 * pi * n ** 2 * e5) * mp.exp(-pi * n * n * e4)
        s += term
        if abs(term) < TINY and n > 3:
            return s
        n += 1


def fY_numeric(y):
    return mp.diff(lambda v: cdf_Y2(v * v), mp.mpf(y))


rows = []
for us in (['0', '0.25'] if QUICK else ['-0.5', '-0.25', '0', '0.1', '0.25', '0.5']):
    u = mp.mpf(us)
    a = Phi(u)
    b = mp.exp(3 * u) * fY_numeric(mp.exp(2 * u)) / 4
    rows.append({'u': us, 'Phi(u)': S(a), 'e^{3u} f_Y(e^{2u})/4': S(b), 'rel_diff': mp.nstr(abs(a - b) / abs(a), 3)})
out['thm_max_iv_kernel'] = rows
out['kernel_even'] = [{'u': us, 'Phi(u)': S(Phi(us)), 'Phi(-u)': S(Phi('-' + us)),
                       'rel_diff': mp.nstr(abs(Phi(us) - Phi('-' + us)) / Phi(us), 3)} for us in (['0.25'] if QUICK else ['0.1', '0.25', '0.5'])]
rows = []
for zs in (['10'] if QUICK else ['0', '1', '10']):
    z = mp.mpf(zs)
    a = mp.quad(lambda u: Phi(u) * mp.cos(z * u), [0, 0.25, 0.5, 0.75, 1, 1.5])
    b = mp.re(xi(HALF + 1j * z / 2)) / 8
    rows.append({'z': zs, 'int_0^oo Phi(u) cos(zu) du': S(a), 'Xi(z/2)/8': S(b), 'rel_diff': mp.nstr(abs(a - b) / abs(b), 3)})
out['kernel_fourier'] = rows
out['kernel_mass'] = {'Xi(0)/4': S(xi(HALF) / 4), 'int_R Phi = twice int_0^oo Phi': S(2 * mp.quad(Phi, [0, 0.25, 0.5, 1, 1.5]))}
log('(iv) done')

# ---------------------------------------------------------------- 5. Theorem thm:khinchin
rows = []
for sg in (['2'] if QUICK else ['1.5', '2', '3.5']):
    sigma = mp.mpf(sg)
    a = 2 ** (1 - sigma)
    base = mellin_quad(sigma)
    for qs in (['1', '2'] if QUICK else ['0.5', '1', '2', '3']):
        q = mp.mpf(qs)
        lhs = mellin_quad(sigma + 2 * q) / base * (1 - a) / (1 - a * 4 ** (-q))
        rhs = pi ** (-q) * mp.gamma(1 + sigma / 2 + q) / mp.gamma(1 + sigma / 2) * mp.zeta(sigma + 2 * q) / mp.zeta(sigma)
        rows.append({'sigma': sg, 'q': qs, 'E[(M~ 4^-K)^q] (quadrature)': S(lhs), 'E[(G/(pi N^2))^q]': S(rhs),
                     'rel_diff': mp.nstr(abs(lhs - rhs) / abs(rhs), 3)})
out['thm_khinchin_real_moments'] = rows
rows = []
for sg in (['2'] if QUICK else ['1.5', '2', '3.5']):
    sigma = mp.mpf(sg)
    a = 2 ** (1 - sigma)
    for tt in (['5'] if QUICK else ['1', '5', '20']):
        t = mp.mpf(tt)
        lhs = mellin_closed(sigma + 1j * t) / mellin_closed(sigma) * (1 - a) / (1 - a * mp.exp(-1j * t * mp.log(2)))
        rhs = pi ** (-1j * t / 2) * mp.gamma(1 + sigma / 2 + 1j * t / 2) / mp.gamma(1 + sigma / 2) * mp.zeta(sigma + 1j * t) / mp.zeta(sigma)
        rows.append({'sigma': sg, 't': tt, 'lhs': S(lhs, 20), 'rhs': S(rhs, 20), 'abs_diff': mp.nstr(abs(lhs - rhs), 3)})
out['thm_khinchin_characteristic_function'] = rows
log('Theorem thm:khinchin done')

# ---------------------------------------------------------------- 6. figure
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
os.makedirs('fig', exist_ok=True)
tt = np.linspace(0.005, 1.5, 60 if QUICK else 300)
exact = np.array([float(theta4(t)) for t in tt])


def trunc(N):
    p = np.ones_like(tt)
    for n in range(1, N + 1):
        p *= 1 - np.exp(-np.pi * n * tt)
    for m in range(1, N + 1):
        p *= 1 - np.exp(-np.pi * (2 * m - 1) * tt)
    return p


fig, ax = plt.subplots(figsize=(4.4, 2.8))
ax.plot(tt, exact, color='black', lw=1.8, label=r'$\theta_4(t)=\mathbb{P}(M\leq t)$')
ax.plot(tt, trunc(3), '--', color='#555555', lw=1.3, label='product truncated at 3')
ax.plot(tt, trunc(1), ':', color='#888888', lw=1.8, label='product truncated at 1')
ax.set_xlabel(r'$t$')
ax.legend(frameon=False, fontsize=8, loc='lower right')
fig.tight_layout()
fig.savefig('fig/maxwald_maxcdf.pdf')
out['figure'] = {'file': 'fig/maxwald_maxcdf.pdf', 't_range': [0.005, 1.5],
                 'max_gap_theta4_minus_product_truncated_at_3': float(np.max(np.abs(exact - trunc(3)))),
                 'max_gap_theta4_minus_product_truncated_at_1': float(np.max(np.abs(exact - trunc(1))))}
out['seconds'] = round(time.time() - T0, 1)
with open('maxwald_maxcdf.json', 'w') as f:
    json.dump(out, f, indent=1)
log('done; wrote maxwald_maxcdf.json and fig/maxwald_maxcdf.pdf')
