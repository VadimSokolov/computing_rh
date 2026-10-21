# tf_primes.py: Table tf:tab:mellin of the chapter "Three Families of Generalised Gamma Convolutions and the Deconvolution Obstruction" (lines 141 to 148, and the statements of lines 12, 135 and 154 that rest on it) and Proposition tf:prop:mellin itself (equation tf:eq:mellin). Added for the book in October 2026, written from the mathematics (the chapter author's script three_families.py was neither read nor run).
# Definitions. chi_4 is the character modulo 4 (chi_4(n) = 0 for n even, +1 for n = 1 mod 4, -1 for n = 3 mod 4), Lambda the von Mangoldt function, s = alpha + i b. The error is E(N, s) = sum_{n <= 40 N} chi_4(n) Lambda(n) n^(-s) e^(-n/N) + L'/L(s, chi_4); the cut at 40 N drops a tail below e^(-40) (4e-18). For alpha in {1.5, 1.2, 1.0, 0.9, 0.8, 0.7, 0.6}, N in {1e3, 1e4, 1e5, 2.5e5}, b in {10, 25} the table prints the mean over the two b of |E| and the least squares slope of log10 of that mean on log10 N, against the prediction 1/2 - alpha.
# Method. (a) The prime powers up to 40 N_max = 1e7 are sieved with numpy (664579 primes, no powers of 2 are kept because chi_4 vanishes on them); each sum is a double precision sum of n^(-s) e^(-n/N) over the prime powers, added with math.fsum, so its error comes from the phases b log n (about 1e-13 in absolute terms; four of the sums at N = 1e3 are repeated in mpmath at 25 digits and the largest difference is printed). (b) L'/L(s, chi_4) is computed by mpmath at 25 digits from dirichlet(s, chi, derivative = 1) and dirichlet(s, chi, 0). (c) The zeros of L(s, chi_4) up to height T (62 in the full case) are the sign changes of the real Hardy function Z(t) = e^(i theta(t)) L(1/2 + i t, chi_4), theta(t) = arg of (4/pi)^((s+1)/2) Gamma((s+1)/2) at s = 1/2 + i t (this product times L is the completed L function, which is real on the line because chi_4 is real and odd with root number 1), found on a grid of step 0.1 at 15 digits and refined by the Illinois method at 25 digits. (d) Proposition tf:prop:mellin is tested: the right side of equation tf:eq:mellin, minus the sum over the zeros rho = 1/2 +- i gamma (all those found), minus the sum over the trivial zeros -1, -3, ..., -15 of Gamma(rho - s) N^(rho - s), plus the sum over m = 1, ..., 9 of (-1)^m/m! (-L'/L)(s - m, chi_4) N^(-m), is compared with E(N, s) for all 56 cases; the zeros above height T contribute less than e^(-pi (T - 25)/2) (below 1e-24 in the full case) and the terms after m = 9 less than 1e-26. (e) The claims of line 135 and 154: the largest deviation of a slope from 1/2 - alpha over alpha = 1.2, ..., 0.6 (the chapter says 0.005), the number of decades in which the error at alpha = 0.6 falls by a factor 10 (the chapter says 10), and, at alpha = 1.5, the sizes of the zero sum and of the O(1/N) part (the sum over m) separately.
# Error control. The comparison with the explicit formula is an end to end test of the sieve, of L'/L, of the zeros and of the formula; a zero missed by the scan or an error in a sign would show up as a difference of the size of the missing term, and the observed largest difference is printed and stored (it measures the double precision noise of the sieve sums). The Hardy function is checked to be real (the largest imaginary part over the scan is printed). Nothing is proved here, and nothing is said about GRH: the zeros found are the real zeros of Z up to height T, and the agreement with the explicit formula shows that no zero within about 25 of the heights used was missed.
# Output: data/tf_primes.json (numbers as strings; E as real and imaginary parts, the table, the slopes, the zeros, the explicit formula comparison); the log is data/tf_primes.log. Usage: python3 tf_primes.py [--quick] [--out FILE]; --quick keeps alpha in {1.0, 0.6}, N in {1e3, 1e4}, sieves to 4e5 and uses zeros up to height 40 (about 2 s in one process) and writes to /tmp unless --out is given; the full case takes about 6 s on one core (run locally within the 30 s limit on 17 October 2026) and also fits a single Hopper core (bash misc/tools/hopper_run.sh -c 1 -t 30 -g tf_primes.json code/tf_primes.py -- --out tf_primes.json).
import argparse
import json
import math
import os
import sys
import time

import mpmath as mp
import numpy as np

T0 = time.time()
CHI = [0, 1, 0, -1]
ALPHAS = ['1.5', '1.2', '1.0', '0.9', '0.8', '0.7', '0.6']
NS = ['1000', '10000', '100000', '250000']
BS = ['10', '25']
CUT = 40
MMAX = 9
NTRIV = 8


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def sd(x, n=12):
    return mp.nstr(x, n)


def sieve_prime_powers(X):
    sv = np.ones(X + 1, dtype=bool)
    sv[:2] = False
    sv[4::2] = False
    for p in range(3, math.isqrt(X) + 1, 2):
        if sv[p]:
            sv[p * p::2 * p] = False
    primes = np.nonzero(sv)[0].astype(np.int64)
    ns, bases = [primes], [primes]
    for p in primes[primes <= math.isqrt(X)]:
        q = int(p) * int(p)
        qs = []
        while q <= X:
            qs.append(q)
            q *= int(p)
        if qs:
            ns.append(np.array(qs, dtype=np.int64))
            bases.append(np.full(len(qs), int(p), dtype=np.int64))
    n = np.concatenate(ns)
    base = np.concatenate(bases)
    keep = (n % 2 == 1)
    n, base = n[keep], base[keep]
    order = np.argsort(n, kind='stable')
    return n[order], base[order], len(primes)


def LL(s):
    return mp.dirichlet(s, CHI, 1) / mp.dirichlet(s, CHI, 0)


def hardy(t):
    s = mp.mpc(mp.mpf(1) / 2, t)
    pre = (mp.mpf(4) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2)
    return pre * mp.dirichlet(s, CHI, 0) / abs(pre)


def find_zeros(T, step='0.1'):
    step = mp.mpf(step)
    ts, vals, maxim = [], [], mp.mpf(0)
    with mp.workdps(15):
        for k in range(1, int(T / step) + 1):
            t = k * step
            v = hardy(t)
            ts.append(t)
            vals.append(v.real)
            maxim = max(maxim, abs(v.imag))
    roots = []
    with mp.workdps(25):
        for i in range(len(ts) - 1):
            if vals[i] * vals[i + 1] < 0:
                roots.append(mp.findroot(lambda x: hardy(x).real, (ts[i], ts[i + 1]), solver='illinois'))
    return roots, maxim


def explicit_rhs(a, b, N, zeros, LLs):
    """Right side of equation tf:eq:mellin minus the first term -L'/L(s): -sum_rho Gamma(rho - s) N^(rho - s) + sum_m (-1)^m/m! (-L'/L)(s - m) N^(-m)."""
    s = mp.mpc(mp.mpf(a), mp.mpf(b))
    Nm = mp.mpf(N)
    logN = mp.log(Nm)
    Znt = mp.mpc(0)
    for g in zeros:
        for sgn in (1, -1):
            rho = mp.mpc(mp.mpf(1) / 2, sgn * g)
            Znt += mp.gamma(rho - s) * mp.exp((rho - s) * logN)
    Ztr = mp.mpc(0)
    for j in range(1, NTRIV + 1):
        rho = mp.mpf(-(2 * j - 1))
        Ztr += mp.gamma(rho - s) * mp.exp((rho - s) * logN)
    M = mp.mpc(0)
    for m in range(1, MMAX + 1):
        M += (-1) ** m / mp.factorial(m) * (-LLs[m]) * Nm ** (-m)
    return -Znt - Ztr + M, Znt, Ztr, M


def mp_prime_sum(a, b, N, n, base):
    s = mp.mpc(mp.mpf(a), mp.mpf(b))
    Nm = mp.mpf(N)
    m = int(np.searchsorted(n, CUT * int(N), side='right'))
    tot = mp.mpc(0)
    for k in range(m):
        nn = int(n[k])
        chi = 1 if nn % 4 == 1 else -1
        tot += chi * mp.log(int(base[k])) * mp.exp(-s * mp.log(nn) - nn / Nm)
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    quick = args.quick
    alphas = ['1.0', '0.6'] if quick else ALPHAS
    Ns = ['1000', '10000'] if quick else NS
    T = 40 if quick else 62
    mp.mp.dps = 25
    out = args.out or ('/tmp/tf_primes_quick.json' if quick else 'data/tf_primes.json')
    X = CUT * max(int(N) for N in Ns)
    log('mode', 'quick' if quick else 'full', 'alphas', alphas, 'N', Ns, 'sieve to', X, 'zeros up to', T)

    n, base, nprimes = sieve_prime_powers(X)
    lam = np.log(base.astype(np.float64))
    chi = np.where(n % 4 == 1, 1.0, -1.0)
    w = chi * lam
    logn = np.log(n.astype(np.float64))
    log('sieve: %d primes (%d odd prime powers kept) up to %d' % (nprimes, len(n), X))

    sums = {}
    for a in alphas:
        for b in BS:
            wz = w * np.exp(-float(a) * logn) * np.exp(-1j * float(b) * logn)
            for N in Ns:
                m = int(np.searchsorted(n, CUT * int(N), side='right'))
                t = wz[:m] * np.exp(-n[:m].astype(np.float64) / float(N))
                sums[(a, b, N)] = (math.fsum(t.real), math.fsum(t.imag))
    log('prime sums done:', len(sums))

    # L'/L(s - m) for m = 0..MMAX
    LLs = {}
    for a in alphas:
        for b in BS:
            LLs[(a, b)] = [LL(mp.mpc(mp.mpf(a) - m, mp.mpf(b))) for m in range(MMAX + 1)]
    log('L\'/L values done:', len(LLs) * (MMAX + 1))

    # the check of the double precision sums in mpmath at N = 1e3
    chk = []
    for a in ([alphas[0], alphas[-1]]):
        for b in BS:
            ex = mp_prime_sum(a, b, '1000', n, base)
            sd_ = sums[(a, b, '1000')]
            chk.append(abs(ex - mp.mpc(*sd_)))
    sieve_check = max(chk)
    log('largest |double sum - mpmath sum| at N = 1000 over %d sums: %s' % (len(chk), sd(sieve_check, 3)))

    zeros, maxim = find_zeros(T)
    gaps = [zeros[i + 1] - zeros[i] for i in range(len(zeros) - 1)]
    log('zeros of L(s, chi_4) up to %s: %d; first five %s; smallest gap %s; max |Im Z| over the scan %s' % (T, len(zeros), [sd(z, 12) for z in zeros[:5]], sd(min(gaps), 4), sd(maxim, 3)))

    # table, explicit formula comparison
    table, comp = [], []
    mean = {}
    worst = mp.mpf(0)
    worst_rel = mp.mpf(0)
    for a in alphas:
        for N in Ns:
            Es, absE = [], []
            for b in BS:
                E = mp.mpc(*sums[(a, b, N)]) + LLs[(a, b)][0]
                rhs, Znt, Ztr, M = explicit_rhs(a, b, N, zeros, LLs[(a, b)])
                diff = abs(E - rhs)
                worst = max(worst, diff)
                worst_rel = max(worst_rel, diff / abs(E))
                Es.append(E)
                absE.append(abs(E))
                comp.append({'alpha': a, 'N': N, 'b': b, 'E_re': sd(E.real, 12), 'E_im': sd(E.imag, 12), 'abs_E': sd(abs(E), 8),
                             'rhs_re': sd(rhs.real, 12), 'rhs_im': sd(rhs.imag, 12), 'abs_diff': sd(diff, 4),
                             'abs_zero_sum': sd(abs(Znt), 6), 'abs_trivial_zero_sum': sd(abs(Ztr), 6), 'abs_m_sum': sd(abs(M), 6)})
            mean[(a, N)] = (absE[0] + absE[1]) / 2
            table.append({'alpha': a, 'N': N, 'abs_E_b10': sd(absE[0], 8), 'abs_E_b25': sd(absE[1], 8), 'mean_abs_E': sd(mean[(a, N)], 8)})
    log('explicit formula against E: largest |E - rhs| over %d cases %s, largest relative %s' % (len(comp), sd(worst, 4), sd(worst_rel, 4)))

    # slopes
    x = np.log10([float(N) for N in Ns])
    slopes = []
    dev = {}
    for a in alphas:
        y = np.log10([float(mean[(a, N)]) for N in Ns])
        sl = float(np.polyfit(x, y, 1)[0])
        pred = 0.5 - float(a)
        dev[a] = abs(sl - pred)
        slopes.append({'alpha': a, 'slope': '%.6f' % sl, 'prediction': '%.2f' % pred, 'abs_deviation': '%.6f' % abs(sl - pred)})

    # the printed table, in the form of the chapter
    log('alpha | mean |E| at N = ' + ', '.join(Ns) + ' | slope | 1/2 - alpha')
    for a, sl in zip(alphas, slopes):
        log('%5s | %s | %s | %s' % (a, '  '.join('%.3e' % float(mean[(a, N)]) for N in Ns), sl['slope'], sl['prediction']))
    claims = {}
    inrange = [a for a in alphas if 0.6 <= float(a) <= 1.2]
    if inrange:
        claims['max_slope_deviation_alpha_0.6_to_1.2'] = '%.6f' % max(dev[a] for a in inrange)
        claims['alpha_of_max_deviation'] = max(inrange, key=lambda a: dev[a])
    if '0.6' in dev:
        sl06 = [s for s in slopes if s['alpha'] == '0.6'][0]['slope']
        claims['alpha_0.6_decades_for_a_factor_10'] = '%.4f' % (1 / abs(float(sl06)))
    if '1.5' in alphas:
        claims['alpha_1.5_zero_sum_and_m_sum'] = [{'N': c['N'], 'b': c['b'], 'abs_zero_sum': c['abs_zero_sum'], 'abs_m_sum': c['abs_m_sum'], 'abs_E': c['abs_E']} for c in comp if c['alpha'] == '1.5']
    log('claims:', {k: v for k, v in claims.items() if k != 'alpha_1.5_zero_sum_and_m_sum'})
    if '1.5' in alphas:
        for c in claims['alpha_1.5_zero_sum_and_m_sum']:
            log('alpha 1.5 N %s b %s: |zero sum| %s |m sum| %s |E| %s' % (c['N'], c['b'], c['abs_zero_sum'], c['abs_m_sum'], c['abs_E']))

    res = {'script': 'tf_primes.py', 'mode': 'quick' if quick else 'full',
           'parameters': {'alphas': alphas, 'N': Ns, 'b': BS, 'cut_factor': CUT, 'sieve_limit': X, 'primes': nprimes, 'mp_dps': 25,
                          'zero_height_T': T, 'm_terms': MMAX, 'trivial_zeros_used': NTRIV,
                          'double_sum_vs_mpmath_largest_abs_diff_N1000': sd(sieve_check, 4)},
           'zeros_of_L_chi4': {'count_up_to_T': len(zeros), 'ordinates': [sd(z, 22) for z in zeros], 'smallest_gap': sd(min(gaps), 6),
                               'max_abs_imag_of_hardy_over_scan': sd(maxim, 3), 'scan_step': '0.1'},
           'table': table, 'slopes': slopes, 'claims': claims,
           'explicit_formula': {'max_abs_diff': sd(worst, 6), 'max_rel_diff': sd(worst_rel, 6), 'cases': comp},
           'seconds': round(time.time() - T0, 1)}
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as fh:
        json.dump(res, fh, indent=1)
    log('wrote', out)


if __name__ == '__main__':
    main()
