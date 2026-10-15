# maxwald_FT.py: the residue expansion of the distribution function of the reciprocal Xi law T (Laplace transform R(lambda) = xi(1/2)/xi(1/2 +
# sqrt(lambda))), the failure of the maximum representation, and Gnedenko's norming, for Section sec:maxwald (ch/maxwald2.tex): mw:prop:nomax
# with mw:fig:FT and mw:prop:gnedenko (Propositions prop:nomax and prop:gnedenko of the authors' chapter). Added for the book in October 2026
# from the checks of the authors' chapter (misc/reviews/nick-maxwald/), recomputed from the mathematics (mpmath at 30 digits, 60 to 100 where
# said; the zeros are the first of data/zeros_1700.txt, 70 digits, read from the run folder or ../data; nothing is taken from the chapter's scripts).
# With Xi(t) = xi(1/2 + it) and gamma_k the ordinates of the zeros, a_k = gamma_k^2:
#  1. c_k = 2 Xi(0)/(gamma_k Xi'(gamma_k)), the residue of R(lambda) e^{lambda t}/lambda at -a_k, with Xi'(gamma) = i xi'(1/2 + i gamma) =
#     i (s(s-1) pi^{-s/2} Gamma(s/2)/2) zeta'(s) at the zero s (zeta(s) = 0); c_1 to c_4 at 30 digits, a check by numerical differentiation of
#     Xi, the alternation of the signs, and the growth log|c_k| against pi gamma_k/4 (150 zeros).
#  2. F_T(t) = 1 + sum_{k<=150} c_k e^{-a_k t} on a printed grid of [0.008, 0.06] (step 0.0005): monotone, density f_T = -sum c_k a_k e^{-a_k t} > 0,
#     the size of the last term; the figure fig/maxwald_FT.pdf.
#  3. Proposition prop:nomax: a_1, a_2, a_3 and 3 a_1, the zeros with a_k < 3 a_1; the three atoms nu = (-c_1, c_1^2/2, -c_2) of -log F_T below
#     3 a_1 and their signs; the remainder -log F_T + c_1 e^{-a_1 t} - c_1^2 e^{-2 a_1 t}/2 + c_2 e^{-a_2 t} times e^{3 a_1 t} against its next
#     terms (-c_1^3/3 - c_3 e^{-(a_3 - 3 a_1) t} + c_1 c_2 e^{-(a_2 - 2 a_1) t}); for theta_4 the coefficient -2 + c_1^4/4 = 2; a witness that
#     -log F_T is not completely monotone: the signs of (-1)^n d^n/dt^n (-log F_T) at t_0 from the Taylor series of log F_T at t_0 (exact recurrence),
#     for n up to 1.4 a_2 t_0 (negative near n = a_2 t_0, where the atom of c_2 > 0 dominates).
#  4. Proposition prop:gnedenko: (i) 1 - theta_4(t) = 2 e^{-pi t} + O(e^{-4 pi t}) and the Gumbel limit of P(pi max M_j - log 2n <= x) for
#     n = 10^2, 10^4, 10^6; (ii) P(T > t)/(|c_1| e^{-a_1 t}) and the Gumbel limit of P(a_1 max T_j - log(|c_1| n) <= x) for n up to 10^8, from F_T.
#  5. The zero free inversion: the Talbot inversion (mpmath) of (1 - R(lambda))/lambda, which uses xi but no zero, at t = 0.04, 0.05, 0.06, 0.07, 0.08;
#     the local log slopes on the four windows against a_1, the intercept |c_1| at t = 0.08, and the comparison with the residue series.
# Run from a folder holding this script and zeros_1700.txt: python3 maxwald_FT.py [--quick]; writes maxwald_FT.json and fig/maxwald_FT.pdf
# (created if absent). --quick is the laptop test (30 zeros, one witness). On Hopper (about 15 minutes, one core):
# bash misc/tools/hopper_run.sh -c 2 -t 90 -g maxwald_FT.json -g fig/maxwald_FT.pdf code/maxwald_FT.py data/zeros_1700.txt
import json
import os
import sys
import time

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QUICK = '--quick' in sys.argv
DPS = 30
mp.mp.dps = DPS
NZ = 30 if QUICK else 150
HALF = mp.mpf(1) / 2
pi = mp.pi
T0 = time.time()
out = {'working_precision_dps': DPS, 'quick': QUICK, 'zeros_used': NZ, 'library': 'mpmath ' + mp.__version__}


def S(x, n=None):
    return mp.nstr(x, n or DPS)


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def load_zeros(n, dps=DPS):
    here = os.path.dirname(os.path.abspath(__file__))
    for p in ('zeros_1700.txt', os.path.join(here, 'zeros_1700.txt'), os.path.join(here, '..', 'data', 'zeros_1700.txt')):
        if os.path.exists(p):
            with mp.workdps(dps):
                z = [mp.mpf(l.strip()) for l in open(p) if l.strip()]
            if len(z) >= n:
                return z[:n]
    with mp.workdps(dps):
        return [mp.zetazero(k).imag for k in range(1, n + 1)]


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def pref(s):
    return s * (s - 1) / 2 * pi ** (-s / 2) * mp.gamma(s / 2)


def c_of(g, xi0):
    """(c, Im Xi'(g)) with c = 2 Xi(0)/(g Xi'(g)), Xi'(g) = i xi'(1/2 + i g) and xi'(rho) = pref(rho) zeta'(rho) at a zero."""
    rho = HALF + 1j * g
    dXi = 1j * pref(rho) * mp.zeta(rho, derivative=1)
    return 2 * xi0 / (g * mp.re(dXi)), mp.im(dXi)


XI0 = xi(HALF)
gam = load_zeros(NZ)
a = [g * g for g in gam]
cs, imparts = [], []
for g in gam:
    c, im = c_of(g, XI0)
    cs.append(c)
    imparts.append(im)
log('c_k for', NZ, 'zeros; c_1 =', S(cs[0], 12))

# ---------------------------------------------------------------- 1. the residues
out['Xi(0)'] = S(XI0)
rows = []
for k in range(8):
    rows.append({'k': k + 1, 'gamma_k': S(gam[k], 25), 'a_k = gamma_k^2': S(a[k]), 'c_k': S(cs[k]), 'abs(Im Xi\'(gamma_k)) (should be 0)': mp.nstr(abs(imparts[k]), 3)})
out['residues_first_eight'] = rows
chk = []
for k in range(4):
    dnum = mp.diff(lambda z: mp.re(xi(HALF + 1j * z)), gam[k])
    cnum = 2 * XI0 / (gam[k] * dnum)
    chk.append({'k': k + 1, 'c_k by numerical differentiation of Xi': S(cnum), 'rel_diff': mp.nstr(abs(cnum - cs[k]) / abs(cs[k]), 3)})
out['residues_numerical_derivative_check'] = chk
CHC = {'c_1': (-50.871, 0.001, cs[0]), 'c_2': (2664.55, 0.01, cs[1]), 'c_3': (-31419.8, 0.1, cs[2]), 'c_4': (1.355e6, 1e3, cs[3])}
out['chapter_comparison_residues'] = {k: {'chapter': v[0], 'unit': v[1], 'ours': float(v[2]), 'agree': bool(abs(float(v[2]) - v[0]) <= 0.5 * v[1] * (1 + 1e-9))}
                                      for k, v in CHC.items()}
out['all_residues'] = {'gamma_k': [S(g, 25) for g in gam], 'c_k': [S(c, 25) for c in cs]}
out['signs_alternate_c_k_c_k+1_negative'] = bool(all(cs[k] * cs[k + 1] < 0 for k in range(NZ - 1)))
out['max_abs_Im_Xi_prime'] = mp.nstr(max(abs(x) for x in imparts), 3)
lc = np.array([float(mp.log(abs(c))) for c in cs])
gg = np.array([float(g) for g in gam])
dev = lc - np.pi * gg / 4
sl = np.polyfit(gg[10:], lc[10:], 1)
out['growth'] = {'log|c_k| - pi gamma_k/4: min, max over k': [float(dev.min()), float(dev.max())], 'least squares slope of log|c_k| on gamma_k (k > 10)': float(sl[0]),
                 'pi/4': float(np.pi / 4), 'log|c_k| - pi gamma_k/4 at k = 10, 50, last': [float(dev[9]), float(dev[min(49, NZ - 1)]), float(dev[-1])]}

# ---------------------------------------------------------------- 2. F_T on a printed grid
def FT(t):
    t = mp.mpf(t)
    return 1 + mp.fsum(c * mp.exp(-ak * t) for c, ak in zip(cs, a))


def fT(t):
    t = mp.mpf(t)
    return -mp.fsum(c * ak * mp.exp(-ak * t) for c, ak in zip(cs, a))


def sfT(t):
    t = mp.mpf(t)
    return -mp.fsum(c * mp.exp(-ak * t) for c, ak in zip(cs, a))


step = mp.mpf('0.002') if QUICK else mp.mpf('0.0005')
tgrid = [mp.mpf('0.008') + j * step for j in range(int(round((mp.mpf('0.06') - mp.mpf('0.008')) / step)) + 1)]
Fv = [FT(t) for t in tgrid]
fv = [fT(t) for t in tgrid]
out['F_T_grid'] = {'t': [float(t) for t in tgrid], 'F_T': [S(x, 20) for x in Fv], 'f_T': [S(x, 12) for x in fv]}
out['F_T_monotone_on_grid'] = bool(all(Fv[j + 1] > Fv[j] for j in range(len(Fv) - 1)))
out['f_T_positive_on_grid'] = bool(all(x > 0 for x in fv))
out['F_T_range_on_grid'] = [S(min(Fv), 10), S(max(Fv), 15)]
out['smallest_increment_of_F_T_on_grid'] = S(min(Fv[j + 1] - Fv[j] for j in range(len(Fv) - 1)), 10)
out['F_T(0.008)'] = S(Fv[0], 15)
out['log10_of_last_term_at_t_0.008'] = float(mp.log10(abs(cs[-1] * mp.exp(-a[-1] * mp.mpf('0.008')))))
tt12 = [mp.mpf('0.008') + j * (mp.mpf('0.052') / 119) for j in range(120)]
v12 = [FT(t) for t in tt12]
out['F_T_on_linspace(0.008,0.06,120)'] = {'monotone': bool(all(v12[j + 1] > v12[j] for j in range(119))), 'min': S(min(v12), 10), 'max': S(max(v12), 15)}
log('F_T monotone on grid', out['F_T_monotone_on_grid'], 'density positive', out['f_T_positive_on_grid'])

plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
os.makedirs('fig', exist_ok=True)
tt = np.linspace(0.008, 0.06, 60 if QUICK else 200)
fig, ax = plt.subplots(figsize=(4.4, 2.8))
ax.plot(tt, [float(FT(t)) for t in tt], color='black', lw=1.8)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$F_T(t)$')
fig.tight_layout()
fig.savefig('fig/maxwald_FT.pdf')
out['figure'] = {'file': 'fig/maxwald_FT.pdf', 't_range': [0.008, 0.06], 'zeros': NZ}

# ---------------------------------------------------------------- 3. Proposition prop:nomax
nom = {'a_1': S(a[0]), 'a_2': S(a[1]), 'a_3': S(a[2]), '3 a_1': S(3 * a[0]), '2 a_1': S(2 * a[0]), 'a_2/a_1': S(a[1] / a[0]),
       'number of k with a_k < 3 a_1': int(sum(1 for x in a if x < 3 * a[0])), '2 a_1 < a_2 < 3 a_1': bool(2 * a[0] < a[1] < 3 * a[0]),
       'c_2 > 0': bool(cs[1] > 0)}
nom['atoms_nu_below_3a_1'] = {'-c_1 at a_1': S(-cs[0], 15), 'c_1^2/2 at 2 a_1': S(cs[0] ** 2 / 2, 15), '-c_2 at a_2': S(-cs[1], 15),
                              'is_positive (expected True, True, False)': [bool(-cs[0] > 0), bool(cs[0] ** 2 / 2 > 0), bool(-cs[1] > 0)]}
nom['theta_4 coefficient of e^{-4 pi t}: -c_2 + c_1^4/4 with c_n = 2(-1)^n'] = S(-mp.mpf(2) + mp.mpf(-2) ** 4 / 4)
CHN = {'3 a_1 (L496)': (599.4, 0.1, float(3 * a[0])), 'a_1 (L496)': (199.79, 0.01, float(a[0])), 'a_2 (L496)': (441.93, 0.01, float(a[1])),
       'a_2/a_1 (L503)': (2.2119, 0.0001, float(a[1] / a[0]))}
nom['chapter_comparison'] = {k: {'chapter': v[0], 'unit': v[1], 'ours': v[2], 'agree': bool(abs(v[2] - v[0]) <= 0.5 * v[1] * (1 + 1e-9))} for k, v in CHN.items()}
rows = []
with mp.workdps(80):
    g80 = load_zeros(12, 80)
    x80 = mp.mpf(1) / 2
    xi80 = xi(x80)
    c80 = [c_of(g, xi80)[0] for g in g80]
    a80 = [g * g for g in g80]
    for ts in ['0.06', '0.08', '0.1', '0.15']:
        t = mp.mpf(ts)
        F = 1 + mp.fsum(c * mp.exp(-ak * t) for c, ak in zip(c80, a80))
        gfun = -mp.log(F)
        main3 = -c80[0] * mp.exp(-a80[0] * t) + c80[0] ** 2 / 2 * mp.exp(-2 * a80[0] * t) - c80[1] * mp.exp(-a80[1] * t)
        scaled = (gfun - main3) * mp.exp(3 * a80[0] * t)
        pred = -c80[0] ** 3 / 3 - c80[2] * mp.exp(-(a80[2] - 3 * a80[0]) * t) + c80[0] * c80[1] * mp.exp(-(a80[1] - 2 * a80[0]) * t)
        rows.append({'t': ts, 'remainder times e^{3 a_1 t}': mp.nstr(scaled, 12), 'next three terms': mp.nstr(pred, 12),
                     'rel_diff': mp.nstr(abs(scaled - pred) / abs(pred), 3), '-c_1^3/3': mp.nstr(-c80[0] ** 3 / 3, 12)})
nom['expansion_remainder'] = rows


def cm_witness(t0, kmax, J, dps):
    """Signs of s_n = (-1)^n d^n/dt^n (-log F_T) at t_0 for n = 1..kmax, from the Taylor series of log F_T(t_0 + h) (exact recurrence)."""
    with mp.workdps(dps):
        xi_d = xi(mp.mpf(1) / 2)
        gz = load_zeros(J, dps)
        ck = [c_of(g, xi_d)[0] for g in gz]
        ak = [g * g for g in gz]
        t0m = mp.mpf(t0)
        w = [c * mp.exp(-x * t0m) for c, x in zip(ck, ak)]
        p = [1 + mp.fsum(w)]
        term = list(w)
        for j in range(1, kmax + 1):
            for k in range(J):
                term[k] = term[k] * (-ak[k]) / j
            p.append(mp.fsum(term))
        ell = [mp.mpf(0)] * (kmax + 1)
        for n in range(1, kmax + 1):
            acc = n * p[n]
            for m in range(1, n):
                acc -= m * ell[m] * p[n - m]
            ell[n] = acc / (n * p[0])
        # g = -log F, g^{(n)}(t_0) = -n! ell_n, so sign(s_n) = -(-1)^n sign(ell_n)
        sg = [0] + [-(-1) ** n * (1 if ell[n] > 0 else -1) for n in range(1, kmax + 1)]
        return sg, float(ak[0]), float(ak[1])


wit = []
plan = [(1, 40), (1, 70)] if QUICK else [(0.5, 60), (0.5, 100), (1, 60), (1, 100), (2, 80), (2, 120)]
sg_prev = {}
for t0, dps in plan:
    t1 = time.time()
    a2f = float(a[1])
    kmax = int(1.4 * a2f * t0) + 5
    sg, a1f, a2f = cm_witness(t0, kmax, 12 if QUICK else 20, dps)
    neg = [n for n in range(1, kmax + 1) if sg[n] < 0]
    runs = []
    for n in range(1, kmax + 1):
        if not runs or runs[-1][0] != sg[n]:
            runs.append([sg[n], n, n])
        else:
            runs[-1][2] = n
    d = {'t0': t0, 'dps': dps, 'orders_tested': kmax, 'orders_n_with_negative_s_n': len(neg),
         'first_negative_order': neg[0] if neg else None, 'last_negative_order': neg[-1] if neg else None,
         'order_a_1 t0': int(round(a1f * t0)), 'sign_at_order_a_1 t0': sg[max(1, int(round(a1f * t0)))],
         'order_a_2 t0': int(round(a2f * t0)), 'sign_at_order_a_2 t0': sg[int(round(a2f * t0))],
         'orders_1_to_10_all_positive': bool(all(sg[n] > 0 for n in range(1, 11))), 'sign_runs_(sign, from_order, to_order)': runs,
         'seconds': round(time.time() - t1, 1)}
    if t0 in sg_prev:
        d['sign_vector_equal_to_the_run_at_lower_precision'] = bool(sg_prev[t0] == sg)
    sg_prev[t0] = sg
    wit.append(d)
    log('witness t0 =', t0, 'dps', dps, ':', len(neg), 'negative orders, first', d['first_negative_order'], 'last', d['last_negative_order'])
nom['not_completely_monotone_witness'] = wit
out['prop_nomax'] = nom

# ---------------------------------------------------------------- 4. Proposition prop:gnedenko
def theta4(t):
    t = mp.mpf(t)
    if t >= mp.mpf('0.3'):
        s = mp.mpf(1)
        n = 1
        while True:
            term = mp.exp(-pi * n * n * t)
            s += 2 * (-1) ** n * term
            if term < mp.mpf(10) ** (-DPS - 8):
                return s
            n += 1
    s = mp.mpf(0)
    n = 0
    while True:
        term = mp.exp(-pi * (n + HALF) ** 2 / t)
        s += 2 * term
        if term < mp.mpf(10) ** (-DPS - 8) * s:
            return s / mp.sqrt(t)
        n += 1


gum = {'(i) M: 1 - theta_4(t) over 2 e^{-pi t}': [{'t': ts, 'ratio': S((1 - theta4(ts)) / (2 * mp.exp(-pi * mp.mpf(ts))), 15), '1 - e^{-3 pi t}': S(1 - mp.exp(-3 * pi * mp.mpf(ts)), 15)}
                                                  for ts in ['0.25', '0.5', '1', '2']]}
xs = [mp.mpf(-3) + mp.mpf(11) * j / 44 for j in range(45)]
rows = []
for n in [10 ** 2, 10 ** 4, 10 ** 6]:
    dm = max(abs(theta4((x + mp.log(2 * n)) / pi) ** n - mp.exp(-mp.exp(-x))) for x in xs)
    rows.append({'n': n, 'sup_x_in_[-3,8]_|P(pi max M_j - log 2n <= x) - exp(-e^{-x})|': mp.nstr(dm, 5)})
gum['(i) M: Gumbel limit'] = rows
rows = []
for ts in ['0.02', '0.04', '0.08', '0.2']:
    t = mp.mpf(ts)
    rows.append({'t': ts, 'P(T>t) over |c_1| e^{-a_1 t}': S(sfT(t) / (-cs[0] * mp.exp(-a[0] * t)), 12)})
gum['(ii) T: tail over the leading term'] = rows
rows = []
for n in [10 ** 2, 10 ** 4, 10 ** 6, 10 ** 8]:
    dm = mp.mpf(0)
    for x in xs:
        t = (x + mp.log(-cs[0] * n)) / a[0]
        dm = max(dm, abs(FT(t) ** n - mp.exp(-mp.exp(-x))))
    rows.append({'n': n, 'sup_x_in_[-3,8]_|P(a_1 max T_j - log(|c_1| n) <= x) - exp(-e^{-x})|': mp.nstr(dm, 5)})
gum['(ii) T: Gumbel limit'] = rows
out['prop_gnedenko'] = gum
log('Gumbel checks done')

# ---------------------------------------------------------------- 5. zero free inversion
def R(lam):
    return XI0 / xi(HALF + mp.sqrt(lam))


def tail_hat(t):
    return mp.invertlaplace(lambda lam: (1 - R(lam)) / lam, mp.mpf(t), method='talbot')


tl = ['0.04', '0.05', '0.06', '0.07', '0.08']
ph = [tail_hat(t) for t in tl]
ps = [sfT(t) for t in tl]
rows = []
for j, t in enumerate(tl):
    rows.append({'t': t, 'P(T>t) by Talbot inversion (no zeros)': S(ph[j], 20), 'P(T>t) by the residue series': S(ps[j], 20), 'rel_diff': mp.nstr(abs(ph[j] - ps[j]) / ps[j], 3)})
sl_hat = [-(mp.log(ph[j + 1]) - mp.log(ph[j])) / (mp.mpf(tl[j + 1]) - mp.mpf(tl[j])) for j in range(4)]
sl_ser = [-(mp.log(ps[j + 1]) - mp.log(ps[j])) / (mp.mpf(tl[j + 1]) - mp.mpf(tl[j])) for j in range(4)]
c1_hat_true = ph[-1] * mp.exp(a[0] * mp.mpf('0.08'))
c1_hat_free = ph[-1] * mp.exp(sl_hat[-1] * mp.mpf('0.08'))
inv = {'rows': rows, 'windows': ['[0.04,0.05]', '[0.05,0.06]', '[0.06,0.07]', '[0.07,0.08]'], 'local_log_slopes_zero_free': [S(x, 12) for x in sl_hat],
       'local_log_slopes_residue_series': [S(x, 12) for x in sl_ser], 'a_1 = gamma_1^2': S(a[0], 12),
       'intercept |c_1| at t = 0.08 with the true a_1': S(c1_hat_true, 10), 'intercept |c_1| at t = 0.08 with the zero free slope': S(c1_hat_free, 10),
       'residue value |c_1|': S(-cs[0], 10), 'max_rel_diff_Talbot_vs_series': mp.nstr(max(abs(ph[j] - ps[j]) / ps[j] for j in range(5)), 3)}
CHI = [('slope [0.04,0.05] (L528)', 199.496, 0.001, sl_hat[0]), ('slope [0.05,0.06] (L528)', 199.764, 0.001, sl_hat[1]),
       ('slope [0.06,0.07] (L528)', 199.788, 0.001, sl_hat[2]), ('slope [0.07,0.08] (L528)', 199.790, 0.001, sl_hat[3]),
       ('gamma_1^2 (L528)', 199.7905, 0.0001, a[0]), ('intercept |c_1| at t = 0.08 with the true a_1 (L528)', 50.87103, 0.00001, c1_hat_true),
       ('residue value |c_1| (L528)', 50.87104, 0.00001, -cs[0])]
inv['chapter_comparison'] = {k: {'chapter': v, 'unit': u, 'ours': float(o), 'agree': bool(abs(float(o) - v) <= 0.5 * u * (1 + 1e-9))} for k, v, u, o in CHI}
out['zero_free_inversion'] = inv
out['seconds'] = round(time.time() - T0, 1)
with open('maxwald_FT.json', 'w') as f:
    json.dump(out, f, indent=1)
log('done; wrote maxwald_FT.json and fig/maxwald_FT.pdf')
