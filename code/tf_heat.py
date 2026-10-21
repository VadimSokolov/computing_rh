# tf_heat.py: the numbers of Section tf:sec:heat ("The heat flow: a stable step that cannot carry real zeros") of the chapter "Three Families of Generalised Gamma Convolutions and the Deconvolution Obstruction" (lines 271 to 273, 283, 299, 301, 327 and 25 of the chapter source). Added for the book in October 2026 and written from the mathematics (the chapter author's script three_families.py was neither read nor run).
# Definitions. Kernel and flow in the normalisation of Rodgers and Tao, as in Section sec:ch21:flow of the book (the chapter only refers to that section): Phi_dBN(u) = sum_{n>=1} (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u)) exp(-pi n^2 e^(4u)), H_t(z) = int_0^inf e^(t u^2) Phi_dBN(u) cos(z u) du, H_0(z) = Xi(z/2)/8, Xi(theta) = xi(1/2 + i theta), xi(s) = s (s - 1) pi^(-s/2) Gamma(s/2) zeta(s)/2. All terms of the series for Phi_dBN are positive for u >= 0.
# Checks. (1) H_0(z) from the integral against Xi(z/2)/8 at z = 0, 5, 20, 28.27 (and 10): the integral is a Gauss Legendre rule on [0, 1.5] (200 nodes, nodes and weights polished by Newton's method at the working precision; a second rule with 160 nodes gives the quadrature error), Xi(z/2)/8 comes from mpmath's zeta and gamma, a different route; both at 40 digits, and the number of agreeing digits is reported. (2) The Gaussian average of Theorem tf:thm:heat(b): E H_{0.2}(z + sqrt(0.4) N) against H_0(z) at z = 0, 10, 28.27 (5 and 20 as well), with Gauss Legendre in u (200 nodes) and Gauss Hermite in N (80 nodes), first in double precision (numpy, the recipe of line 327) and then at 40 digits, which tells the rounding error from the quadrature error; relative and absolute errors; an extra case t = s = 1/2. z = 28.27 lies 5.5e-4 from the zero 2 gamma_1 of H_0, where |H_0| is about 5e-8, so the relative error of the double precision value is dominated by rounding in values of H_{0.2} of size 0.1, and the absolute error is the meaningful number. The double precision value is compared with H_0 at the double nearest to z. (3) Phi_dBN(1.5) and a proved bound for every u >= 1.5 (claimed below 1e-200 in line 327). (4) sympy: E (z + sqrt(2 s) N)^2 - c = z^2 - c + 2 s, whose zeros are real exactly when s <= c/2 (line 283). (5) The arithmetic of lines 301 and 25: 2 eps*^2 at eps* = 1/2 and at 3/8 = 7/8 - 1/2, 9/32 against 0.2, 1/2 + sqrt(0.1).
# Error control. The mpmath integrals use fixed rules whose error is read from the second rule (160 nodes against 200, Hermite 40 against 80) and from the agreement with Xi(z/2)/8; the sums of the weights are checked against 2 (Legendre on [-1, 1]) and sqrt(pi) (Hermite). Nothing here is a proof; the bound of (3) is the only analytic statement, and it is proved in the comment at tail_bound: for u >= 3/2, Phi_dBN(u) <= 4 pi^2 exp(9u - pi e^(4u)), which decreases because 9 - 4 pi e^(4u) < 0, so the value at 3/2 bounds all of them.
# Output: data/tf_heat.json (values as strings, with the parameters and the precision); the log is data/tf_heat.log. Usage: python3 tf_heat.py [--quick] [--out FILE] [--dps D]; --quick uses 25 digits, 60 and 40 Legendre nodes, 30 and 16 Hermite nodes and z = 0, 28.27 only (a few seconds in one process) and writes to /tmp unless --out is given; the full case takes about 2 s in one process (run locally on 17 October 2026).
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np

T0 = time.time()


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def sd(x, n=30):
    return mp.nstr(x, n)


def phi_dbn(u):
    """Phi_dBN(u), a sum of positive terms, to the working precision."""
    u = mp.mpf(u)
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    stop = mp.mpf(10) ** (-(mp.mp.dps + 6))
    tot = mp.mpf(0)
    n = 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        tot += term
        if term < stop * tot:
            return tot
        n += 1


def phi_double(u):
    u = np.asarray(u, dtype=float)
    e4, e5, e9 = np.exp(4 * u), np.exp(5 * u), np.exp(9 * u)
    tot = np.zeros_like(u)
    for n in range(1, 13):
        tot += (2 * np.pi ** 2 * n ** 4 * e9 - 3 * np.pi * n ** 2 * e5) * np.exp(-np.pi * n * n * e4)
    return tot


def xi_func(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def gl_half(n):
    """Positive nodes and weights of the n point Gauss Legendre rule on [-1, 1] (n even), polished at the working precision."""
    x0 = np.polynomial.legendre.leggauss(n)[0][n // 2:]
    tol = mp.mpf(10) ** (-mp.mp.dps - 4)
    xs, ws = [], []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(12):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < tol:
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def gh_half(n):
    """Positive nodes and weights of the n point Gauss Hermite rule for the weight exp(-x^2) (n even), polished."""
    x0 = np.polynomial.hermite.hermgauss(n)[0][n // 2:]
    tol = mp.mpf(10) ** (-mp.mp.dps - 4)
    xs, ws = [], []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(12):
            h0, h1 = mp.mpf(1), 2 * x
            for k in range(1, n):
                h0, h1 = h1, 2 * x * h1 - 2 * k * h0
            dx = h1 / (2 * n * h0)
            x -= dx
            if abs(dx) < tol:
                break
        h0, h1 = mp.mpf(1), 2 * x
        for k in range(1, n):
            h0, h1 = h1, 2 * x * h1 - 2 * k * h0
        xs.append(x)
        ws.append(mp.mpf(2) ** (n - 1) * mp.factorial(n) * mp.sqrt(mp.pi) / (n * n * h0 * h0))
    return xs, ws


def make_rule(n, U):
    xs, ws = gl_half(n)
    u, W = [], []
    for x, w in zip(xs, ws):
        u += [U / 2 * (1 + x), U / 2 * (1 - x)]
        W += [U / 2 * w, U / 2 * w]
    phi = [phi_dbn(v) for v in u]
    return u, W, phi


def h0_rule(z, rule):
    u, W, phi = rule
    return mp.fsum(Wk * pk * mp.cos(z * uk) for uk, Wk, pk in zip(u, W, phi))


def avg_mp(zs_list, t, s, rule, gh):
    """E H_t(z + sqrt(2 s) N) = (1/sqrt(pi)) sum_i w_i H_t(z + 2 sqrt(s) x_i), the nodes paired as cos((z+a)u) + cos((z-a)u) = 2 cos(zu) cos(au)."""
    u, W, phi = rule
    taus, wts = gh
    A = [Wk * mp.exp(t * uk * uk) * pk for uk, Wk, pk in zip(u, W, phi)]
    a = [2 * mp.sqrt(s) * tau for tau in taus]
    C = [[mp.cos(ai * uk) for uk in u] for ai in a]
    out = {}
    for zs in zs_list:
        z = mp.mpf(zs)
        B = [Ak * mp.cos(z * uk) for Ak, uk in zip(A, u)]
        tot = mp.mpf(0)
        for wi, Ci in zip(wts, C):
            tot += wi * mp.fsum(Bk * Cik for Bk, Cik in zip(B, Ci))
        out[zs] = 2 * tot / mp.sqrt(mp.pi)
    return out


def double_nodes_numpy(nu, nh, U=1.5):
    """The nodes and weights numpy returns (its documentation says leggauss and hermgauss are only tested up to degree 100)."""
    xl, wl = np.polynomial.legendre.leggauss(nu)
    th, wh = np.polynomial.hermite.hermgauss(nh)
    return 0.5 * U * (xl + 1.0), 0.5 * U * wl, th, wh


def double_nodes_polished(rule, gh):
    """The nodes and weights polished at the working precision (40 digits) and then rounded to double."""
    u = np.array([float(v) for v in rule[0]])
    W = np.array([float(v) for v in rule[1]])
    tp = np.array([float(v) for v in gh[0]])
    wp = np.array([float(v) for v in gh[1]])
    return u, W, np.concatenate((-tp[::-1], tp)), np.concatenate((wp[::-1], wp))


def avg_double(z, t, s, nodes):
    u, W, th, wh = nodes
    A = W * np.exp(t * u * u) * phi_double(u)
    x = z + math.sqrt(2.0 * s) * math.sqrt(2.0) * th
    H = np.cos(np.outer(x, u)) @ A
    return float(np.dot(wh, H) / math.sqrt(math.pi))


def digits_agree(a, b):
    d = abs(a - b)
    if d == 0:
        return 'all'
    return mp.nstr(-mp.log10(d / abs(b)), 4)


def tail_bound():
    """For u >= 3/2, x = e^(4u) >= e^6 > 400 and n^2 - 1 >= 3 (n - 1), so sum_n n^4 exp(-pi n^2 x) <= 2 exp(-pi x); the terms of Phi_dBN are at most 2 pi^2 n^4 e^(9u) exp(-pi n^2 x), hence Phi_dBN(u) <= 4 pi^2 exp(9u - pi e^(4u)) =: B(u); B'/B = 9 - 4 pi e^(4u) < 0 for u >= 3/2, so B(u) <= B(3/2)."""
    u0 = mp.mpf('1.5')
    B = 4 * mp.pi ** 2 * mp.exp(9 * u0 - mp.pi * mp.exp(4 * u0))
    phi15 = phi_dbn(u0)
    # the tail of the integral for t = 0.2: the exponent 0.2 u^2 + 9 u - pi e^(4u) also decreases beyond 3/2
    slope = mp.mpf('0.4') * u0 + 9 - 4 * mp.pi * mp.exp(4 * u0)
    out = {
        'Phi_dBN(1.5)': sd(phi15, 12), 'log10_Phi_dBN(1.5)': sd(mp.log10(phi15), 8),
        'bound_B(1.5)': sd(B, 12), 'log10_bound_B(1.5)': sd(mp.log10(B), 8),
        'B_decreasing_beyond_1.5 (9 - 4 pi e^6 < 0)': bool(9 - 4 * mp.pi * mp.exp(4 * u0) < 0),
        'exponent_slope_for_t_0.2_at_1.5': sd(slope, 8),
        'Phi_dBN(2)': sd(phi_dbn(2), 8), 'log10_Phi_dBN(2)': sd(mp.log10(phi_dbn(2)), 8),
        'Phi_dBN(3)': sd(phi_dbn(3), 8), 'log10_Phi_dBN(3)': sd(mp.log10(phi_dbn(3)), 8),
        'below_1e-200_for_all_u_ge_1.5': bool(B < mp.mpf(10) ** -200),
    }
    return out


def sympy_part():
    import sympy as sp
    from sympy.stats import Normal, E
    z, c = sp.symbols('z c', real=True)
    s = sp.symbols('s', positive=True)
    N = Normal('N', 0, 1)
    avg = sp.simplify(E((z + sp.sqrt(2 * s) * N) ** 2 - c))
    target = z ** 2 - c + 2 * s
    ok = sp.simplify(avg - target) == 0
    disc = sp.simplify(sp.discriminant(sp.expand(avg), z))
    thr = sp.solve(sp.Eq(disc, 0), s)
    return {'E[(z + sqrt(2 s) N)^2 - c]': str(avg), 'equals z^2 - c + 2 s': bool(ok),
            'discriminant_in_z': str(disc), 'discriminant_zero_at_s': [str(v) for v in thr],
            'real_zeros_exactly_when': 's <= c/2' if [str(v) for v in thr] == ['c/2'] else 'check'}


def arithmetic_part():
    half = Fraction(1, 2)
    e38 = Fraction(7, 8) - half
    y = mp.sqrt(mp.mpf('0.1'))
    return {'eps_star_from_7_over_8': str(e38), '2 eps*^2 at eps*=1/2': str(2 * half ** 2),
            '2 eps*^2 at eps*=3/8': str(2 * e38 ** 2), '9/32 as decimal': str(float(Fraction(9, 32))),
            '9/32 greater than 0.2': bool(Fraction(9, 32) > Fraction(1, 5)),
            'y with 2 y^2 = 0.2 (sqrt 0.1)': sd(y, 12), '1/2 + sqrt(0.1)': sd(half + y, 12),
            'rounds_to_0.816': bool(abs(float(half + y) - 0.816) < 5e-4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--dps', type=int, default=None)
    args = ap.parse_args()
    quick = args.quick
    dps = args.dps or (25 if quick else 40)
    mp.mp.dps = dps
    out = args.out or ('/tmp/tf_heat_quick.json' if quick else 'data/tf_heat.json')
    nu, nu2, nh, nh2 = (60, 40, 30, 16) if quick else (200, 160, 80, 40)
    zs_norm = ['0', '28.27'] if quick else ['0', '5', '20', '28.27']
    zs_avg = ['0', '28.27'] if quick else ['0', '10', '28.27']
    zs_all = sorted(set(zs_norm + zs_avg), key=float)
    cases = [('0.2', '0.2')] if quick else [('0.2', '0.2'), ('0.5', '0.5')]
    U = mp.mpf('1.5')
    log('mode', 'quick' if quick else 'full', 'dps', dps, 'Legendre', nu, nu2, 'Hermite', nh, nh2)

    # nodes and weights
    ruleA, ruleB = make_rule(nu, U), make_rule(nu2, U)
    wsum_err = max(abs(mp.fsum(ruleA[1]) - U), abs(mp.fsum(ruleB[1]) - U))
    ghA, ghB = gh_half(nh), gh_half(nh2)
    gh_err = max(abs(2 * mp.fsum(ghA[1]) - mp.sqrt(mp.pi)), abs(2 * mp.fsum(ghB[1]) - mp.sqrt(mp.pi)))
    gh_m2 = abs(2 * mp.fsum(w * x * x for x, w in zip(*ghA)) - mp.sqrt(mp.pi) / 2)
    log('node checks: Legendre weight sum error %s, Hermite weight sum error %s, second moment error %s' % (sd(wsum_err, 3), sd(gh_err, 3), sd(gh_m2, 3)))
    nodes_np = double_nodes_numpy(nu, nh)
    nodes_acc = double_nodes_polished(ruleA, ghA)
    idx = np.argsort(nodes_acc[0])
    node_diff = {'legendre_max_abs_node_diff': '%.3g' % np.max(np.abs(nodes_np[0] - nodes_acc[0][idx])),
                 'legendre_max_rel_weight_diff': '%.3g' % np.max(np.abs(nodes_np[1] / nodes_acc[1][idx] - 1)),
                 'hermite_max_abs_node_diff': '%.3g' % np.max(np.abs(nodes_np[2] - nodes_acc[2])),
                 'hermite_max_rel_weight_diff': '%.3g' % np.max(np.abs(nodes_np[3] / nodes_acc[3] - 1))}
    log('numpy nodes against the polished ones:', node_diff)

    # (1) normalisation
    norm = []
    ref = {}
    for zs in zs_all:
        z = mp.mpf(zs)
        ha, hb = h0_rule(z, ruleA), h0_rule(z, ruleB)
        xv = xi_func(mp.mpc(mp.mpf(1) / 2, z / 2)) / 8
        xr = xv.real
        ref[zs] = ha
        rec = {'z': zs, 'H0_integral': sd(ha, dps - 3), 'Xi_over_8': sd(xr, dps - 3), 'abs_diff': sd(abs(ha - xr), 4),
               'rel_diff': sd(abs(ha - xr) / abs(xr), 4), 'digits_agreeing': digits_agree(ha, xr),
               'rule_B_abs_diff_from_A': sd(abs(ha - hb), 4), 'imag_part_of_Xi_over_8': sd(abs(xv.imag), 3)}
        norm.append(rec)
        log('H0(%s) = %s  Xi(z/2)/8 = %s  agree %s digits  (second rule differs by %s)' % (zs, sd(ha, 20), sd(xr, 20), rec['digits_agreeing'], rec['rule_B_abs_diff_from_A']))

    # (2) the Gaussian average
    avg = []
    for ts, ss in cases:
        t, s = mp.mpf(ts), mp.mpf(ss)
        vA = avg_mp(zs_avg, t, s, ruleA, ghA)
        vB = avg_mp(zs_avg, t, s, ruleB, ghB)
        for zs in zs_avg:
            r = ref[zs]
            zd = float(zs)
            r_d = h0_rule(mp.mpf(zd), ruleA)
            rec = {'t': ts, 's': ss, 'sigma': sd(mp.sqrt(2 * s), 12), 'z': zs, 'H0_reference': sd(r, 25),
                   'mp_average': sd(vA[zs], dps - 3), 'mp_abs_err': sd(abs(vA[zs] - r), 4), 'mp_rel_err': sd(abs(vA[zs] - r) / abs(r), 4),
                   'mp_abs_err_smaller_rules': sd(abs(vB[zs] - r), 4), 'double_reference_at_nearest_double_z': sd(r_d, 20)}
            for label, nodes in (('double_numpy_nodes', nodes_np), ('double_accurate_nodes', nodes_acc)):
                vd = avg_double(zd, float(ts), float(ss), nodes)
                e_d = abs(mp.mpf(vd) - r_d)
                rec[label] = {'average': repr(vd), 'abs_err': sd(e_d, 4), 'rel_err': sd(e_d / abs(r_d), 4)}
            avg.append(rec)
            log('t=s=%s z=%s: double (numpy nodes) abs %s rel %s | double (accurate nodes) abs %s rel %s | mp (%d x %d) abs %s rel %s | mp (%d x %d) abs %s' % (ts, zs, rec['double_numpy_nodes']['abs_err'], rec['double_numpy_nodes']['rel_err'], rec['double_accurate_nodes']['abs_err'], rec['double_accurate_nodes']['rel_err'], nu, nh, rec['mp_abs_err'], rec['mp_rel_err'], nu2, nh2, rec['mp_abs_err_smaller_rules']))

    # (3) tail, (4) sympy, (5) arithmetic, (6) the zero near 28.27
    tail = tail_bound()
    log('tail:', tail['log10_Phi_dBN(1.5)'], 'log10 Phi(1.5);', tail['log10_bound_B(1.5)'], 'log10 bound; below 1e-200 for all u >= 1.5:', tail['below_1e-200_for_all_u_ge_1.5'])
    sym = sympy_part()
    log('sympy:', sym['E[(z + sqrt(2 s) N)^2 - c]'], 'equals z^2 - c + 2 s:', sym['equals z^2 - c + 2 s'], 'real zeros exactly when', sym['real_zeros_exactly_when'])
    ar = arithmetic_part()
    log('arithmetic:', ar)
    g1 = mp.zetazero(1).imag
    zero = {'gamma_1': sd(g1, 25), '2 gamma_1': sd(2 * g1, 25), 'distance_of_28.27': sd(abs(mp.mpf('28.27') - 2 * g1), 6),
            'within_6e-4': bool(abs(mp.mpf('28.27') - 2 * g1) < mp.mpf('6e-4')), 'H0(28.27)': sd(ref['28.27'], 12)}
    log('zero:', zero)

    res = {'script': 'tf_heat.py', 'mode': 'quick' if quick else 'full',
           'parameters': {'dps': dps, 'U': '1.5', 'legendre_nodes': nu, 'legendre_nodes_second_rule': nu2, 'hermite_nodes': nh,
                          'hermite_nodes_second_rule': nh2, 'double_precision_legendre': nu, 'double_precision_hermite': nh,
                          'node_checks': {'legendre_weight_sum_err': sd(wsum_err, 3), 'hermite_weight_sum_err': sd(gh_err, 3), 'hermite_second_moment_err': sd(gh_m2, 3)},
                          'numpy_nodes_against_polished_nodes': node_diff},
           'definition': 'Phi_dBN(u) = sum_{n>=1} (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u)) exp(-pi n^2 e^(4u)); H_t(z) = int_0^inf e^(t u^2) Phi_dBN(u) cos(z u) du; Xi(theta) = xi(1/2 + i theta), xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s); Rodgers and Tao normalisation, Section sec:ch21:flow of the book',
           'normalisation': norm, 'average': avg, 'tail': tail, 'sympy': sym, 'arithmetic': ar, 'zero_near_28.27': zero,
           'seconds': round(time.time() - T0, 1)}
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as fh:
        json.dump(res, fh, indent=1)
    log('wrote', out)


if __name__ == '__main__':
    main()
