# maxwald_basepoint.py: the numbers of the basepoint family of Section sec:maxwald (Proposition mw:prop:lambda, the comparison of the
# line, zero and prime forms of the Thorin density, and the curves of Figure mw:fig:basepoint). Added for the book in October 2026 from
# the checks of the authors' chapter, misc/reviews/nick-maxwald/ (review-2 M2 and M7 and the recomputation list); written from the
# mathematics, not from the chapter's verify scripts. For the basepoint law R_alpha(lambda) = xi(alpha)/xi(alpha + sqrt(lambda)) the
# Thorin density is u_alpha(lambda) = Re (xi'/xi)(alpha + i sqrt(lambda))/(2 pi sqrt(lambda)) (line form), which equals under RH
# sum_{gamma > 0} term_gamma(lambda), term_gamma = delta (lambda + gamma^2 + delta^2)/(pi sqrt(lambda) |(delta - i gamma)^2 + lambda|^2),
# delta = alpha - 1/2 (zero form), and for alpha > 1 the 2 sin^2 prime form
#   2 pi sqrt(lambda) u = xi'/xi(alpha) + sum_{n >= 2} Lambda(n) n^(-alpha) 2 sin^2(sqrt(lambda) log(n)/2)
#                         + sum_{k >= 0} lambda/((alpha + 2k)((alpha + 2k)^2 + lambda)) - lambda/(alpha (alpha^2 + lambda)) - lambda/((alpha - 1)((alpha - 1)^2 + lambda)).
# 1. One zero's term in the lambda chart (review-2 M2), at alpha in {0.55, 0.7, 1, 1.5, 2} for gamma_1 and at alpha = 0.55 for gamma_5:
#    term_gamma = (lambda + gamma^2 + delta^2)/(2 gamma sqrt(lambda)) times the Cauchy density centred at gamma^2 - delta^2 with half
#    width 2 delta gamma (identity checked at eight values of lambda, and against the Stieltjes inversion (1/pi) Im (g'/g)(-lambda - i0) of
#    the pair factor g = 1 + a sqrt(z) + b z); its mass over (0, oo) (quadrature in sqrt(lambda)); its full width at half height at its
#    true maximum and at half of its value at the Cauchy centre, against 4 delta gamma; the factor at lambda = 1, at the centre and at
#    4 gamma^2; the mass of the Cauchy density restricted to (0, oo).
# 2. The three forms at alpha in {2, 1.5} and lambda in {100, 250, 500, 1000} (mpmath, 30 digits). The line form is the reference. The
#    prime form is summed to N = NMAX/10 and N = NMAX (NMAX = 2e6, the chapter used 2e5), raw and with the tail sum_{n > N} Lambda(n)
#    n^(-alpha) 2 sin^2(.) added back as (a) the smooth tail N^(1 - alpha)/(alpha - 1), (b) the same with its oscillating part
#    N^(1 - alpha)/(alpha - 1) - Re N^(1 - alpha + i sqrt(lambda))/(alpha - 1 - i sqrt(lambda)), (c) (b) corrected by the Chebyshev remainder
#    psi(N) - N: - (psi(N) - N) N^(-alpha) (1 - cos(sqrt(lambda) log N)); for the gamma sum the digamma closed form (1/2) Re[psi((alpha + i sqrt(lambda))/2)
#    - psi(alpha/2)] is used (checked against the series once). The zero form is summed over the first 300 and the first 1700 zeros
#    (the book's table zeros_1700.txt), raw, with the chapter's tail (smooth density log(gamma/2 pi)/(2 pi) and the large gamma form of the
#    term), and with our tail: the integral of the term against the exact smooth density theta'(x)/pi, theta the Riemann-Siegel theta,
#    from the last zero gamma_K, minus term(gamma_K)(K - theta(gamma_K)/pi - 1), the boundary term of the partial summation against
#    S(T) = N(T) - theta(T)/pi - 1. For every variant the script reports the absolute and relative gap to the line form and the digits.
# 3. The curves of the figure: u_alpha by the line form for alpha = 1.5, 1, 0.7, 0.55 on lambda in [20, 1750], on a grid that is refined
#    around gamma_k^2, with the squares gamma_k^2 < 1750 (the maxwald_basepoint_fig.py script draws fig/maxwald_basepoint.pdf from them).
# 4. Two checks of the remark on the three forms: the minimum of u_alpha on [1, 2500] for alpha = 0.55, 0.6, 0.75, 1, 1.5, and the
#    coefficient (1 - phi_alpha(lambda))/sqrt(lambda) at lambda = 1e-6 for alpha = 0.6 and 1.2 against xi'/xi(alpha) (added 14 October 2026).
# Run on Hopper (8 cores, about a minute), from a folder containing zeros_1700.txt: python3 maxwald_basepoint.py [quick|chart]; writes
# maxwald_basepoint.json there ('quick' is a smoke test at small sizes: NMAX = 2e5 and NMAX/10 = 2e4 for the prime sums, one process,
# a coarse grid for the figure; 'chart' stops after part 1).
import json
import math
import os
import sys
import time
import multiprocessing

import numpy as np
import mpmath as mp

QUICK = 'quick' in sys.argv[1:]
DPS = 30
mp.mp.dps = DPS
PI = mp.pi
NMAX = 2 * 10 ** (5 if QUICK else 6)
NSMALL = NMAX // 10
KS = (300, 1700)
NPROC = 1 if QUICK else int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1))
t0 = time.time()
out = dict(script='maxwald_basepoint.py', mp_dps=DPS, NMAX=NMAX, NMAX_small=NSMALL, zero_counts=list(KS), quick=QUICK,
           prime_sums='float64 (math.fsum of numpy terms), digamma closed form for the gamma sum', processes=NPROC)


def say(*a):
    print(*a, flush=True)


def digits(rel):
    return 30.0 if rel == 0 else float(-mp.log10(rel))


gam = [mp.mpf(l.strip()) for l in open('zeros_1700.txt') if l.strip()]
assert len(gam) == 1700 and abs(gam[0] - mp.mpf('14.134725141734693790457251983562')) < mp.mpf('1e-28')
HALF = mp.mpf(1) / 2


def term(lam, alpha, g):
    """The term of the zero 1/2 + i g in the Thorin density of the basepoint law (RH: alpha - beta = alpha - 1/2)."""
    d = alpha - HALF
    return d * (d * d + g * g + lam) / (PI * mp.sqrt(lam) * abs((d - 1j * g) ** 2 + lam) ** 2)


def prefactor(lam, alpha, g):
    d = alpha - HALF
    return (lam + g * g + d * d) / (2 * g * mp.sqrt(lam))


def cauchy(lam, alpha, g):
    d = alpha - HALF
    return (2 * d * g) / (PI * ((lam - (g * g - d * d)) ** 2 + (2 * d * g) ** 2))


def stieltjes_form(lam, alpha, g):
    """(1/pi) Im (g'/g)(-lam - i0) for the pair factor g(z) = 1 + a sqrt(z) + b z, a = 2 (alpha - 1/2)/|c|^2, b = 1/|c|^2, |c|^2 = delta^2 + gamma^2."""
    d = alpha - HALF
    a, b = 2 * d / (d * d + g * g), 1 / (d * d + g * g)
    r = -1j * mp.sqrt(lam)                      # the principal root of -lam - i0
    return mp.im((a / (2 * r) + b) / (1 + a * r - b * lam)) / PI


def bisect(f, lo, hi, n=110):
    flo = f(lo)
    assert flo * f(hi) < 0, 'no sign change'
    for _ in range(n):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def argmax(f, a, b, n=150):
    for _ in range(n):
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        if f(m1) < f(m2):
            a = m1
        else:
            b = m2
    return (a + b) / 2


# ---------------------------------------------------------------- 1. one zero's term in the lambda chart
def chart_checks(alpha_s, k):
    alpha, g = mp.mpf(alpha_s), gam[k - 1]
    d = alpha - HALF
    c, w = g * g - d * d, 2 * d * g
    tests = [mp.mpf(1), c / 2, c - w, c, c + w, 2 * c, 4 * g * g, 10 * g * g]
    f = lambda x: term(x, alpha, g)
    ident = max(abs(f(x) - prefactor(x, alpha, g) * cauchy(x, alpha, g)) / f(x) for x in tests)
    stj = max(abs(f(x) - stieltjes_form(x, alpha, g)) / f(x) for x in tests)
    pts = sorted(set([mp.sqrt(c)] + [mp.sqrt(c + j * w) for j in (-125, -25, -5, -1, 1, 5, 25, 125) if c + j * w > 0]))
    mass = mp.quad(lambda s: 2 * s * f(s * s), [mp.mpf(0)] + pts + [mp.inf])
    xm = argmax(f, c - w, c + w)
    fm = f(xm)
    lo_b, hi_b = max(xm - 3 * w, xm / 4), xm + 3 * w
    fwhm = bisect(lambda x: f(x) - fm / 2, xm, hi_b) - bisect(lambda x: f(x) - fm / 2, lo_b, xm)
    f0 = f(c)
    fwhm0 = bisect(lambda x: f(x) - f0 / 2, c, c + 3 * w) - bisect(lambda x: f(x) - f0 / 2, max(c - 3 * w, c / 4), c)
    return dict(alpha=alpha_s, k=k, gamma=mp.nstr(g, 12), centre=mp.nstr(c, 12), half_width=mp.nstr(w, 12),
                identity_max_relative_error=float(ident), stieltjes_form_max_relative_error=float(stj),
                mass=mp.nstr(mass, 15), mass_minus_one=float(mass - 1),
                argmax=mp.nstr(xm, 12), fwhm_at_true_maximum=mp.nstr(fwhm, 12), fwhm_at_half_of_value_at_centre=mp.nstr(fwhm0, 12),
                fwhm_lorentzian_4_delta_gamma=mp.nstr(4 * d * g, 12), fwhm_relative_excess=float(fwhm / (4 * d * g) - 1),
                prefactor_at_lambda_1=mp.nstr(prefactor(mp.mpf(1), alpha, g), 8), prefactor_at_centre=mp.nstr(prefactor(c, alpha, g), 10),
                one_plus_delta2_over_2gamma2=mp.nstr(1 + d * d / (2 * g * g), 10), prefactor_at_4_gamma2=mp.nstr(prefactor(4 * g * g, alpha, g), 8),
                cauchy_mass_on_0_oo=mp.nstr(HALF + mp.atan(c / w) / PI, 8))


say('== 1. one zero in the lambda chart')
cases = [('0.55', 1), ('0.7', 1), ('1.0', 1), ('1.5', 1), ('2.0', 1), ('0.55', 5)]
out['lambda_chart_term'] = {}
for a_s, k in cases:
    r = chart_checks(a_s, k)
    out['lambda_chart_term']['alpha=%s,gamma_index=%d' % (a_s, k)] = r
    say('alpha=%s k=%d: identity %.1e, Stieltjes form %.1e, mass-1 %.1e, FWHM %s (true max) %s (half of centre value) against 4 d g = %s' % (
        a_s, k, r['identity_max_relative_error'], r['stieltjes_form_max_relative_error'], r['mass_minus_one'],
        r['fwhm_at_true_maximum'][:11], r['fwhm_at_half_of_value_at_centre'][:11], r['fwhm_lorentzian_4_delta_gamma'][:11]))

if 'chart' in sys.argv[1:]:        # part 1 only (a test of the checks of the lambda chart)
    name = os.environ.get('MAXWALD_OUT', 'maxwald_basepoint.json')
    json.dump(out, open(name, 'w'), indent=1)
    say('part 1 only: wrote %s in %.1f s' % (name, time.time() - t0))
    sys.exit(0)


# ---------------------------------------------------------------- 2. the three forms
def ld(s):  # xi'/xi
    return 1 / s + 1 / (s - 1) - mp.log(PI) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)


def u_line(alpha, lam):
    return mp.re(ld(alpha + 1j * mp.sqrt(lam))) / (2 * PI * mp.sqrt(lam))


def prime_powers(N):
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    primes = np.nonzero(sieve)[0].astype(np.int64)
    ns, ls = [primes], [np.log(primes.astype(float))]
    for p in primes[primes <= int(N ** 0.5)]:
        q, lp, qs = int(p) * int(p), math.log(int(p)), []
        while q <= N:
            qs.append(q)
            q *= int(p)
        ns.append(np.array(qs, dtype=np.int64))
        ls.append(np.full(len(qs), lp))
    n, L = np.concatenate(ns), np.concatenate(ls)
    o = np.argsort(n, kind='stable')
    return n[o].astype(np.float64), L[o]


nn, LL = prime_powers(NMAX)
say('prime powers up to %d: %d, psi(NMAX) - NMAX = %.1f (%.1f s)' % (NMAX, len(nn), math.fsum(LL) - NMAX, time.time() - t0))
psiN = {N: math.fsum(LL[:np.searchsorted(nn, N, side='right')]) for N in (NSMALL, NMAX)}


def prime_variants(alpha_s, lam, N):
    alpha, t = float(alpha_s), math.sqrt(lam)
    m = np.searchsorted(nn, N, side='right')
    n = nn[:m]
    P = math.fsum(LL[:m] * n ** (-alpha) * 2 * np.sin(0.5 * t * np.log(n)) ** 2)
    a, T, Nm = mp.mpf(alpha_s), mp.sqrt(lam), mp.mpf(N)
    tail1 = Nm ** (1 - a) / (a - 1)
    tail2 = tail1 - mp.re(Nm ** (1 - a + 1j * T) / (a - 1 - 1j * T))
    tail3 = tail2 - (mp.mpf(psiN[N]) - Nm) * Nm ** (-a) * (1 - mp.cos(T * mp.log(Nm)))
    return mp.mpf(P), dict(raw=mp.mpf(0), smooth=tail1, smooth_osc=tail2, smooth_osc_psi=tail3)


def prime_u(alpha_s, lam, P, tail):
    a, x = mp.mpf(alpha_s), mp.mpf(lam)
    T = mp.sqrt(x)
    gam_sum = mp.re(mp.digamma((a + 1j * T) / 2) - mp.digamma(a / 2)) / 2
    pole = x / (a * (a * a + x)) + x / ((a - 1) * ((a - 1) ** 2 + x))
    return (ld(a) + P + tail + gam_sum - pole) / (2 * PI * T)


def zero_partial(alpha_s, lam, K):
    a, x = mp.mpf(alpha_s), mp.mpf(lam)
    return mp.fsum(term(x, a, g) for g in gam[:K])


def zero_tail_chapter(alpha_s, lam, K):
    d, x = mp.mpf(alpha_s) - HALF, mp.mpf(lam)
    return mp.quad(lambda g: mp.log(g / (2 * PI)) / (2 * PI) * d * (g * g + x) / (PI * mp.sqrt(x) * (g * g - x) ** 2), [gam[K - 1], mp.inf])


def zero_tail_ours(alpha_s, lam, K):
    a, x, G = mp.mpf(alpha_s), mp.mpf(lam), gam[K - 1]
    dNs = lambda y: (mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * y)) - mp.log(PI)) / (2 * PI)
    I = mp.quad(lambda y: term(x, a, y) * dNs(y), [G, 2 * G, 10 * G, 100 * G, 10 ** 4 * G, mp.inf])
    S = K - (mp.siegeltheta(G) / PI + 1)
    return I - term(x, a, G) * S


def gap_record(ref, val, scale):
    gap = ref - val
    rel = abs(gap) / abs(ref)
    return dict(u=mp.nstr(val, 15), gap=float(gap), gap_times_2pi_sqrt_lambda=float(gap * scale), relative_gap=float(rel), digits=digits(rel))


say('== 2. the three forms')
chk = mp.nsum(lambda k: mp.mpf(100) / ((2 + 2 * k) * ((2 + 2 * k) ** 2 + 100)), [0, mp.inf])
out['gamma_sum_series_vs_digamma_alpha2_lambda100'] = float(chk - mp.re(mp.digamma((2 + 10j) / 2) - mp.digamma(mp.mpf(1))) / 2)
ALPHAS_B = ['2.0', '1.5']
LAMS_B = [100, 250, 500, 1000]
out['three_forms'] = {}
for a_s in ALPHAS_B:
    for lam in LAMS_B:
        x, T = mp.mpf(lam), mp.sqrt(lam)
        scale = 2 * PI * T
        ref = u_line(mp.mpf(a_s), x)
        rec = dict(u_line=mp.nstr(ref, 20), prime={}, zero={})
        for N in (NSMALL, NMAX):
            P, tails = prime_variants(a_s, lam, N)
            rec['prime']['N=%d' % N] = {name: gap_record(ref, prime_u(a_s, lam, P, tl), scale) for name, tl in tails.items()}
            rec['prime']['N=%d' % N]['smooth_tail_value_in_2pi_sqrt_lambda_u'] = mp.nstr(tails['smooth'], 10)
        for K in KS:
            part = zero_partial(a_s, lam, K)
            rec['zero']['K=%d' % K] = dict(raw=gap_record(ref, part, scale),
                                           chapter_tail=gap_record(ref, part + zero_tail_chapter(a_s, lam, K), scale),
                                           our_tail=gap_record(ref, part + zero_tail_ours(a_s, lam, K), scale))
        out['three_forms']['alpha=%s,lambda=%d' % (a_s, lam)] = rec
        p5, p6 = rec['prime']['N=%d' % NSMALL], rec['prime']['N=%d' % NMAX]
        say('alpha=%s lambda=%4d u_line=%s | prime N=%.0e raw %.1f dig, smooth %.1f, +osc %.1f, +psi %.1f | N=%.0e raw %.1f, smooth %.1f, +osc %.1f, +psi %.1f | zeros(%d) raw %.1f, ch.tail %.1f, our tail %.1f (%.1fs)' % (
            a_s, lam, mp.nstr(ref, 10), NSMALL, p5['raw']['digits'], p5['smooth']['digits'], p5['smooth_osc']['digits'], p5['smooth_osc_psi']['digits'],
            NMAX, p6['raw']['digits'], p6['smooth']['digits'], p6['smooth_osc']['digits'], p6['smooth_osc_psi']['digits'], KS[-1],
            rec['zero']['K=%d' % KS[-1]]['raw']['digits'], rec['zero']['K=%d' % KS[-1]]['chapter_tail']['digits'],
            rec['zero']['K=%d' % KS[-1]]['our_tail']['digits'], time.time() - t0))

# ---------------------------------------------------------------- 3. the curves of the figure
FIG_ALPHAS = ['1.5', '1.0', '0.7', '0.55']
LO, HI = 20.0, 1750.0
base = np.linspace(LO, HI, 60 if QUICK else 1731)
refine = [np.linspace(float(g * g) - 8 * 0.1 * float(g), float(g * g) + 8 * 0.1 * float(g), 9 if QUICK else 161) for g in gam if g * g < HI]
grid = np.unique(np.concatenate([base] + refine))
grid = grid[(grid >= LO) & (grid <= HI)]


def work(task):
    a_s, lam = task
    with mp.workdps(20):
        return float(u_line(mp.mpf(a_s), mp.mpf(lam)))


tasks = [(a_s, float(lam)) for a_s in FIG_ALPHAS for lam in grid]
if NPROC == 1:
    vals = [work(tk) for tk in tasks]
else:
    with multiprocessing.get_context('fork').Pool(NPROC) as pool:
        vals = pool.map(work, tasks, chunksize=40)
vals = np.array(vals).reshape(len(FIG_ALPHAS), len(grid))
out['figure'] = dict(alphas=[float(a) for a in FIG_ALPHAS], lambda_grid=[float('%.9g' % v) for v in grid],
                     u={a: [float('%.9g' % v) for v in vals[i]] for i, a in enumerate(FIG_ALPHAS)},
                     gamma_squared=[float('%.12g' % (g * g)) for g in gam if g * g < HI], working_dps=20, points=int(len(grid)))
mw = max(float(abs(mp.mpf(vals[2][len(grid) // 3]) - u_line(mp.mpf('0.7'), mp.mpf(float(grid[len(grid) // 3]))))), 0.0)
out['figure']['spot_check_alpha_0.7_abs_difference_at_30_digits'] = mw
say('figure curves: %d points x %d alphas in %.1f s' % (len(grid), len(FIG_ALPHAS), time.time() - t0))


# ---------------------------------------------------------------- 4. two checks of the remark on the three forms
# (a) The minimum of u_alpha on [1, 2500] for alpha = 0.55, 0.6, 0.75, 1, 1.5, by the line form at 20 digits on a grid of step 1/4
#     (step 25 with quick), refined by golden section search within one step of the smallest grid value. On this window every zero is
#     known to lie on the line, so a positive minimum checks the computation and does not test RH.
# (b) The coefficient (1 - phi_alpha(lambda))/sqrt(lambda) at lambda = 1e-6, phi_alpha(lambda) = xi(alpha)/xi(alpha + sqrt(lambda)),
#     at 30 digits for alpha = 0.6 and 1.2, against its limit xi'/xi(alpha); the offset is eps (xi''/(2 xi) - (xi'/xi)^2) + O(eps^2),
#     eps = sqrt(lambda), and its coefficient is printed too.
MIN_ALPHAS = ['0.55', '0.6', '0.75', '1.0', '1.5']
MIN_LO, MIN_HI, MIN_STEP = 1.0, 2500.0, (25.0 if QUICK else 0.25)
min_grid = np.arange(MIN_LO, MIN_HI + MIN_STEP / 2, MIN_STEP)


def work4(task):
    a_s, lam = task
    with mp.workdps(20):
        return float(u_line(mp.mpf(a_s), mp.mpf(lam)))


def golden_min(a_s, lo, hi, n=40):
    with mp.workdps(20):
        f = lambda x: u_line(mp.mpf(a_s), x)
        lo, hi = mp.mpf(lo), mp.mpf(hi)
        r = (mp.sqrt(5) - 1) / 2
        x1, x2 = hi - r * (hi - lo), lo + r * (hi - lo)
        f1, f2 = f(x1), f(x2)
        for _ in range(n):
            if f1 < f2:
                hi, x2, f2 = x2, x1, f1
                x1 = hi - r * (hi - lo)
                f1 = f(x1)
            else:
                lo, x1, f1 = x1, x2, f2
                x2 = lo + r * (hi - lo)
                f2 = f(x2)
        x = (lo + hi) / 2
        return float(x), float(f(x))


tasks4 = [(a_s, float(lam)) for a_s in MIN_ALPHAS for lam in min_grid]
if NPROC == 1:
    vals4 = [work4(tk) for tk in tasks4]
else:
    with multiprocessing.get_context('fork').Pool(NPROC) as pool:
        vals4 = pool.map(work4, tasks4, chunksize=40)
vals4 = np.array(vals4).reshape(len(MIN_ALPHAS), len(min_grid))
out['minimum_of_u_alpha'] = {'window': [MIN_LO, MIN_HI], 'grid_step': MIN_STEP, 'grid_points': int(len(min_grid)), 'working_dps': 20}
for i, a_s in enumerate(MIN_ALPHAS):
    j = int(np.argmin(vals4[i]))
    lo, hi = max(MIN_LO, float(min_grid[j]) - MIN_STEP), min(MIN_HI, float(min_grid[j]) + MIN_STEP)
    xr, vr = golden_min(a_s, lo, hi)
    out['minimum_of_u_alpha']['alpha=' + a_s] = {'grid_min': float('%.9g' % vals4[i][j]), 'grid_argmin': float(min_grid[j]),
                                                 'refined_min': float('%.9g' % vr), 'refined_argmin': float('%.9g' % xr),
                                                 'positive': bool(min(vals4[i].min(), vr) > 0)}
say('minimum of u_alpha on [1, 2500]: %s' % {a: (r['refined_min'], r['refined_argmin']) for a, r in out['minimum_of_u_alpha'].items()
                                             if a.startswith('alpha=')})


def xi_real(s):
    return s * (s - 1) / 2 * mp.power(PI, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


out['small_lambda_coefficient'] = {'lambda': '1e-6', 'working_dps': 30}
with mp.workdps(30):
    lam = mp.mpf('1e-6')
    eps = mp.sqrt(lam)
    for a_s in ['0.6', '1.2']:
        a = mp.mpf(a_s)
        coef = (1 - xi_real(a) / xi_real(a + eps)) / eps
        l1 = ld(a)
        l2 = mp.diff(ld, a) + l1 ** 2          # xi''/xi
        c2 = l2 / 2 - l1 ** 2
        out['small_lambda_coefficient']['alpha=' + a_s] = {
            'one_minus_phi_over_sqrt_lambda': mp.nstr(coef, 12), 'xi_prime_over_xi': mp.nstr(l1, 12),
            'difference': mp.nstr(coef - l1, 6), 'next_taylor_coefficient': mp.nstr(c2, 6),
            'difference_minus_first_order_prediction': mp.nstr(coef - l1 - eps * c2, 3)}
say('small lambda coefficient: %s' % out['small_lambda_coefficient'])

out['seconds'] = round(time.time() - t0, 1)
name = os.environ.get('MAXWALD_OUT', 'maxwald_basepoint.json')
json.dump(out, open(name, 'w'), indent=1)
say('wrote %s in %.1f s' % (name, time.time() - t0))
