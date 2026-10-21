# zh_comb.py: the comb barrier, the horizon heights and the zero data of Chapter 18 (Zhu's compact windows, the Landau Widom horizon and
# the Thorin program), group Z1. Added for the book in October 2026, written from the mathematics (no script of the chapter's author was read or run).
# What it checks, in the order of the chapter. (1) Table zh:tab:comb: the comb constant A_L(alpha), the sum over log n < 2L of 2 Lambda(n) n^(0 minus alpha),
# and the barrier height T_1 = 2 pi exp(A_L) at L = 0.35, 0.8, 1.6, 2, 3, 4, 5 and alpha = 1/2, 3/4, 1 (42 values), with its caption: A_0.8(1/2) = 2.94197
# and T_1 at L = 1.6 and 2 (quoted from Zhu), minus 2 zeta'/zeta(alpha) at alpha = 1.1, 1.5, 2. (2) The three asymptotic forms of the display after the
# heading of Section zh:sec:comb, 4 e^L, 8 e^(L/2) and 4 L minus 2 gamma_E, against the exact values (ratios). (3) The heights 2 pi e^(2L) of the horizon at
# the half lengths of the product table, the window table and Section zh:sec:horizon, the heights 2 pi e^(2L+1) of the definition T_LW, the Landau Widom
# dimension L T*/pi against the count N(T*) of ordinates up to T* = 2 pi e^(2L), and the half length (1/2) log(3e12/2 pi) of the Platt and Trudgian verification.
# (4) The ordinates quoted at lines 53, 97 to 99 and 160. (5) kappa_1 = sum over gamma > 0 of gamma^(0 minus 2), the crossovers lambda/(4 kappa_1) of Section
# zh:sec:damp and the sandwich 4 kappa_1 c minus 4 kappa_2 c^3 <= 2 xi'/xi(1/2 + c) <= 4 kappa_1 c behind the display. (6) Columns 1 and 2 of Table zh:tab:det (the
# midpoints of consecutive ordinates and (1/2) log(gamma/2 pi)) and the offsets and the amplification quoted in the text after the table.
# Method. Lambda(n) from an exact sieve of Eratosthenes to exp(2 L_max); every sum in mpmath at 40 digits (at most 2600 positive terms, so no cancellation); the
# test log n < 2L is made in mpmath arithmetic (the nearest integer to exp(2L) among the table's L is 2 at L = 0.35, where exp(0.7) = 2.0138). Minus 2 zeta'/zeta
# from mpmath's zeta with derivative, checked against the partial sums of the Dirichlet series to 2e6 with the main term N^(1 minus alpha)/(alpha minus 1) as tail
# (three or four digits, which is all that route gives). kappa_1 from the identity sum over gamma > 0 of 1/(gamma^2 + 1/4) = (2 + gamma_E minus log 4 pi)/2 plus the
# rapidly convergent sum of 1/(4 gamma^2 (gamma^2 + 1/4)) over the 1700 ordinates of data/zeros_1700.txt and an explicit Riemann von Mangoldt bound for the rest (Hasanalizade, Shen
# and Wong 2022: |N(T) minus T/(2 pi) log(T/(2 pi e)) minus 7/8| <= 2.5167 + 1/(50 T) for e <= T <= 30610046000, and 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50 T) beyond);
# independently as the limit of xi'/xi(1/2 + c)/(2c) at c = 1e-10 (xi'/xi from the Gamma and zeta factors at 50 digits). Everything about the zeros is under RH.
# Error control. The table values are compared with the chapter's printed digits (half up rounding), and the number of terms of each sum is recorded.
# Output: data/zh_comb.json (values as strings, the claims with their status), the log on standard output (data/zh_comb.log).
# Usage: python3 zh_comb.py [--quick] [--out FILE] [--zeros FILE]; --quick uses 20 digits and L up to 3 (about 1 s) and writes to /tmp unless --out is given.
import argparse
import json
import math
import sys
import tempfile
import time
from pathlib import Path

import mpmath as mp
import numpy as np

T0 = time.time()
CLAIMS = []
X0 = 30610046000.0

LS = ['0.35', '0.8', '1.6', '2.0', '3.0', '4.0', '5.0']
ALPHAS = ['0.5', '0.75', '1']
# the table of lines 177 to 183 as printed: L -> [(A_L, log10 T_1) for alpha = 1/2, 0.75, 1]
PRINTED = {
    '0.35': [('0.980', '1.22'), ('0.824', '1.16'), ('0.693', '1.10')],
    '0.8': [('2.942', '2.08'), ('2.278', '1.79'), ('1.772', '1.57')],
    '1.6': [('14.323', '7.02'), ('8.497', '4.49'), ('5.222', '3.07')],
    '2.0': [('24.383', '11.39'), ('12.566', '6.26'), ('6.874', '3.78')],
    '3.0': [('75.190', '33.45'), ('26.647', '12.37'), ('10.856', '5.51')],
    '4.0': [('213.548', '93.54'), ('49.930', '22.48'), ('14.855', '7.25')],
    '5.0': [('588.024', '256.17'), ('88.183', '39.10'), ('18.844', '8.98')]}
# Table zh:tab:det as printed: gamma -> smallest L for delta = 0.01, 0.05, 0.2, 0.45
DET = {'51.37': ('1.051', ['1.20', '1.15', '1.15', '1.10']), '199.64': ('1.729', ['1.85', '1.80', '1.75', '1.75']),
       '1000.57': ('2.535', ['2.65', '2.65', '2.60', '2.50']), '1999.99': ('2.882', ['3.05', '3.00', '2.95', '2.45'])}


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def ms(x, n=25):
    return mp.nstr(x, n)


def rnd(x, nd):
    q = mp.floor(x * mp.mpf(10) ** nd + mp.mpf(1) / 2) / mp.mpf(10) ** nd
    return '%.*f' % (nd, float(q))


def sig(x, nd):
    return '%.*e' % (nd - 1, float(x))


def claim(line, text, chapter, mine, ok, note=''):
    CLAIMS.append({'line': line, 'claim': text, 'chapter': chapter, 'mine': mine, 'status': 'agrees' if ok else 'differs', 'note': note})


def find_zeros(path):
    here = Path(__file__).resolve().parent
    for c in ([Path(path)] if path else []) + [Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), here / 'zeros_1700.txt',
                                                here.parent / 'data' / 'zeros_1700.txt']:
        if c.is_file():
            return c
    raise SystemExit('zeros_1700.txt not found, use --zeros')


def prime_powers(N):
    sv = bytearray([1]) * (N + 1)
    sv[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
    out = []
    for p in range(2, N + 1):
        if sv[p]:
            q = p
            while q <= N:
                out.append((q, p))
                q *= p
    out.sort()
    return out


def nplus(t):
    """upper bound for the number N(t) of ordinates in (0, t] (HSW22, as used in the book): main term, 7/8 and the explicit error"""
    main = t / (2 * math.pi) * math.log(t / (2 * math.pi * math.e)) + 0.875
    err = 2.5167 if t <= X0 else 0.1095 * math.log(t) + 0.2042 * math.log(math.log(t)) + 3.0305
    return main + err + 1.0 / (50 * t)


def mtz(alpha):
    a = mp.mpf(alpha)
    return -2 * mp.zeta(a, derivative=1) / mp.zeta(a)


def xixi(s):
    """xi'/xi(s) from xi(s) = s (s minus 1) pi^(0 minus s/2) Gamma(s/2) zeta(s)/2"""
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)


def dirichlet_check(alpha, N):
    """float64 partial sum of sum Lambda(n) n^(0 minus alpha) to N plus the main term N^(1 minus alpha)/(alpha minus 1), times 2"""
    sv = np.ones(N + 1, dtype=bool)
    sv[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    tot = 0.0
    for p in np.nonzero(sv)[0]:
        q = int(p)
        lp = math.log(q)
        while q <= N:
            tot += lp * q ** (-alpha)
            q *= int(p)
    return 2 * (tot + N ** (1 - alpha) / (alpha - 1))


def main():
    ap = argparse.ArgumentParser(description='comb barrier, horizon heights and zero data of Chapter 18 (group Z1)')
    ap.add_argument('--quick', action='store_true', help='20 digits, L up to 3, writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/zh_comb.json, or beside the script when there is no data folder; quick: /tmp)')
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    args = ap.parse_args()
    mp.mp.dps = 20 if args.quick else 40
    zpath = find_zeros(args.zeros)
    lines = [ln.strip() for ln in open(zpath) if ln.strip()]
    gam = [mp.mpf(s) for s in lines]
    zd = np.array([float(s) for s in lines])
    ng = len(gam)
    T0z = float(zd[-1])
    log('zeros', ng, 'to height', '%.6f' % T0z, 'file', zpath, 'dps', mp.mp.dps)
    res = {'script': 'zh_comb.py', 'mode': 'quick' if args.quick else 'full',
           'parameters': {'argv': sys.argv[1:], 'dps': mp.mp.dps, 'zeros': ng, 'zero_file': str(zpath), 'highest_ordinate': lines[-1][:20]}}

    # (1) Table zh:tab:comb
    Ls = [x for x in LS if (float(x) <= 3 or not args.quick)]
    Lmax = mp.mpf(max(Ls, key=float))
    N = int(mp.floor(mp.exp(2 * Lmax))) + 1
    pw = prime_powers(N)
    logn = [mp.log(n) for n, p in pw]
    wt = [2 * mp.log(p) for n, p in pw]
    terms = {a: [wt[i] * mp.power(pw[i][0], -mp.mpf(a)) for i in range(len(pw))] for a in ALPHAS}

    def comb_A(L, a):
        cnt = 0
        while cnt < len(pw) and logn[cnt] < 2 * L:
            cnt += 1
        return mp.fsum(terms[a][:cnt]), cnt

    ln10 = mp.log(10)
    rows = []
    nagree = 0
    A = {}
    for Ls_ in Ls:
        L = mp.mpf(Ls_)
        for k, a in enumerate(ALPHAS):
            val, cnt = comb_A(L, a)
            A[(Ls_, a)] = val
            logT = mp.log10(2 * mp.pi) + val / ln10
            pa, pt = PRINTED[Ls_][k]
            oka, okt = rnd(val, 3) == pa, rnd(logT, 2) == pt
            nagree += int(oka) + int(okt)
            rows.append({'L': Ls_, 'alpha': a, 'prime_powers_summed': cnt, 'A_L': ms(val), 'log10_T1': ms(logT), 'T1': ms(2 * mp.pi * mp.exp(val), 12),
                         'printed_A': pa, 'printed_log10_T1': pt, 'mine_A_rounded': rnd(val, 3), 'mine_log10_T1_rounded': rnd(logT, 2),
                         'A_agrees': oka, 'log10_T1_agrees': okt})
            if not (oka and okt):
                log('TABLE DIFFERS at L', Ls_, 'alpha', a, 'mine', rnd(val, 3), rnd(logT, 2), 'printed', pa, pt)
    res['table'] = rows
    res['table_summary'] = {'entries_compared': 2 * len(rows), 'entries_agreeing': nagree}
    log('Table zh:tab:comb: %d of %d printed entries agree to the printed digits (%d prime powers to e^(2L_max))' % (nagree, 2 * len(rows), len(pw)))
    claim('177 to 183', 'Table zh:tab:comb, %d values (A_L and log10 T_1)' % (2 * len(rows)), 'as printed', '%d of %d agree' % (nagree, 2 * len(rows)),
          nagree == 2 * len(rows))

    # caption
    cap = {}
    a08 = A[('0.8', '0.5')]
    cap['A_0.8(1/2)'] = ms(a08)
    claim('186', 'A_0.8(1/2) (Zhu: 2.94197)', '2.94197', rnd(a08, 5), rnd(a08, 5) == '2.94197')
    if '1.6' in Ls:
        T16 = 2 * mp.pi * mp.exp(A[('1.6', '0.5')])
        cap['T1(1.6)'] = ms(T16, 12)
        claim('186, 33', 'T_1 at L = 1.6 (about 1e7)', '1e7', sig(T16, 3), sig(T16, 1) == '1e+07')
    if '2.0' in Ls:
        T20 = 2 * mp.pi * mp.exp(A[('2.0', '0.5')])
        cap['T1(2.0)'] = ms(T20, 12)
        claim('186, 33', 'T_1 at L = 2 (2.4e11)', '2.4e11', sig(T20, 3), sig(T20, 2) == '2.4e+11')
    T08 = 2 * mp.pi * mp.exp(a08)
    cap['T1(0.8)'] = ms(T08, 12)
    chk = []
    for al, pr in (('1.1', '18.88'), ('1.5', '3.01'), ('2', '1.14')):
        v = mtz(al)
        d = dirichlet_check(float(al), 200000 if args.quick else 2000000)
        chk.append({'alpha': al, 'minus_2_zeta_prime_over_zeta': ms(v), 'printed': pr, 'dirichlet_partial_sum_check': '%.5g' % d,
                    'check_relative_difference': '%.2g' % abs(d / float(v) - 1)})
        claim('186', 'minus 2 zeta\'/zeta(%s)' % al, pr, rnd(v, 2), rnd(v, 2) == pr, 'Dirichlet series check %.5g' % d)
    cap['minus_2_zeta_over_zeta'] = chk
    res['caption'] = cap
    log('caption: A_0.8(1/2) =', ms(a08, 12), 'T_1(0.8) =', ms(T08, 8), [c['minus_2_zeta_prime_over_zeta'][:8] for c in chk])

    # (2) asymptotic forms
    asy = []
    for Ls_ in Ls:
        L = mp.mpf(Ls_)
        forms = [('1/2', 4 * mp.exp(L), A[(Ls_, '0.5')]), ('3/4', 8 * mp.exp(L / 2), A[(Ls_, '0.75')]), ('1', 4 * L - 2 * mp.euler, A[(Ls_, '1')])]
        asy.append({'L': Ls_, 'forms': [{'alpha': al, 'asymptotic': ms(f, 12), 'exact': ms(ex, 12), 'ratio_exact_over_asymptotic': ms(ex / f, 8)} for al, f, ex in forms]})
        log('asymptotics L', Ls_, ['%s: %s vs %s ratio %s' % (al, ms(ex, 8), ms(f, 8), ms(ex / f, 6)) for al, f, ex in forms])
    res['asymptotics'] = asy

    # (3) heights
    pt_half = mp.log(mp.mpf(3) * mp.mpf(10) ** 12 / (2 * mp.pi)) / 2
    hlist = [('(1/2) log 2 (Yoshida)', mp.log(2) / 2), ('0.347', mp.mpf('0.347')), ('0.35', mp.mpf('0.35')), ('0.5', mp.mpf('0.5')), ('0.8', mp.mpf('0.8')),
             ('1.0', mp.mpf(1)), ('1.2', mp.mpf('1.2')), ('1.6', mp.mpf('1.6')), ('2.0', mp.mpf(2)), ('13.4', mp.mpf('13.4')), ('Platt Trudgian 13.446', pt_half)]
    heights = []
    for name, L in hlist:
        Ts = 2 * mp.pi * mp.exp(2 * L)
        Tlw = mp.e * Ts
        Nst = int(np.searchsorted(zd, float(Ts), side='right')) if float(Ts) < T0z else None
        nmain = Ts / (2 * mp.pi) * mp.log(Ts / (2 * mp.pi * mp.e)) if Ts > 2 * mp.pi * mp.e else mp.mpf(0)
        Nlw = int(np.searchsorted(zd, float(Tlw), side='right')) if float(Tlw) < T0z else None
        heights.append({'L': name, 'T_star_2pi_e_2L': ms(Ts, 12), 'T_LW_2pi_e_2Lplus1': ms(Tlw, 12), 'N_T_star': Nst, 'N_T_LW': Nlw,
                        'N_main_term_T_star': ms(nmain, 8), 'dimension_L_T_star_over_pi': ms(L * Ts / mp.pi, 8),
                        'dimension_L_T_LW_over_pi': ms(L * Tlw / mp.pi, 8)})
        log('height L', name, 'T*', ms(Ts, 8), '2 pi e^(2L+1)', ms(Tlw, 8), 'N(T*)', Nst, 'N(T_LW)', Nlw, 'dim L T*/pi', ms(L * Ts / mp.pi, 6), 'dim L T_LW/pi', ms(L * Tlw / mp.pi, 6))
    res['heights'] = heights
    res['platt_trudgian'] = {'half_length': ms(pt_half, 12), 'level_2_log_3e12': ms(2 * mp.log(mp.mpf(3) * mp.mpf(10) ** 12), 12),
                             'T_star_at_13.4': ms(2 * mp.pi * mp.exp(mp.mpf('26.8')), 8), 'T_star_at_2L_level': ms(2 * mp.pi * mp.exp(2 * pt_half), 8)}
    Tof = {h['L']: mp.mpf(h['T_star_2pi_e_2L']) for h in heights}
    claim('152', 'height 2 pi e^(2L) at the Yoshida window, L = 0.347', '12.6', rnd(Tof['0.347'], 1), rnd(Tof['0.347'], 1) == '12.6')
    claim('97', 'height 2 pi e^(2L) at L = 0.5', '17.1', rnd(Tof['0.5'], 1), rnd(Tof['0.5'], 1) == '17.1')
    claim('98, 153', 'height 2 pi e^(2L) at L = 0.8', '31.1 (31)', rnd(Tof['0.8'], 1), rnd(Tof['0.8'], 1) == '31.1')
    claim('99', 'height 2 pi e^(2L) at L = 1.0', '46.4', rnd(Tof['1.0'], 1), rnd(Tof['1.0'], 1) == '46.4')
    claim('154, 160, 216', 'height 2 pi e^(2L) at L = 1.6 (about 154, "near 150")', '154', rnd(Tof['1.6'], 0), rnd(Tof['1.6'], 0) == '154')
    claim('155, 18', 'Platt and Trudgian half length (1/2) log(3e12/2 pi) (about 13.4)', '13.4', rnd(pt_half, 1), rnd(pt_half, 1) == '13.4')
    claim('155', 'height 2 pi e^(2L) at L = 13.4 (3e12)', '3e12', sig(Tof['13.4'], 2), sig(Tof['13.4'], 1) == '3e+12',
          'exact at 13.446; at the rounded half length 13.4 the height is 2.7e12')
    nst = {h['L']: h['N_T_star'] for h in heights}
    claim('18, 160', 'Zhu certificate L = 0.8 is worth about the first four zeros: N(T*)', '4', str(nst['0.8']), nst['0.8'] == 4)
    log('Platt Trudgian half length', ms(pt_half, 10), 'level', ms(2 * mp.log(mp.mpf(3) * mp.mpf(10) ** 12), 8))

    # (4) ordinates
    quoted = [(1, '14.13', 2, '53'), (2, '21.0', 1, '97'), (4, '30.42', 2, '160'), (5, '32.94', 2, '160'), (6, '37.6', 1, '98'), (9, '48.0', 1, '99')]
    ords = {}
    for j, pr, nd, ln in quoted:
        ords[str(j)] = ms(gam[j - 1], 15)
        claim(ln, 'ordinate gamma_%d' % j, pr, rnd(gam[j - 1], nd), rnd(gam[j - 1], nd) == pr)
    for j in (8,):
        ords[str(j)] = ms(gam[j - 1], 15)
    res['ordinates'] = ords
    log('ordinates', {k: v[:9] for k, v in ords.items()})

    # (5) kappa_1, crossovers, sandwich
    S2 = mp.fsum(1 / (g * g * (g * g + mp.mpf(1) / 4)) for g in gam) / 4
    S4 = mp.fsum(1 / g ** 4 for g in gam)

    def tail4(coef):
        return mp.quad(lambda t: (nplus(float(t)) - ng) * 4 * t ** -5, [T0z, 1e4, 1e8, mp.inf]) * coef
    t4 = tail4(1)
    A1 = (2 + mp.euler - mp.log(4 * mp.pi)) / 2
    kappa1 = A1 + S2 + t4 / 4
    kappa2 = S4 + t4
    with mp.workdps(60):
        c0 = mp.mpf(10) ** -10
        kappa1_xi = xixi(mp.mpf(1) / 2 + c0) / (2 * c0)
    log('kappa_1 =', ms(kappa1, 14), '(identity route, tail bound', mp.nstr(t4 / 4, 3), ') and', mp.nstr(kappa1_xi, 14), '(xi\'/xi at 1/2 + 1e-10); kappa_2 =', ms(kappa2, 10))
    claim('yo:56, 208', 'kappa_1 = sum over gamma > 0 of gamma^(0 minus 2)', '0.0231049931', rnd(kappa1, 10), rnd(kappa1, 10) == '0.0231049931',
          'identity route and the xi\'/xi limit agree to %s' % mp.nstr(abs(kappa1 - kappa1_xi), 2))
    xs = []
    ok_sand = True
    for c in ('0.01', '0.1', '0.25', '0.4'):
        cm = mp.mpf(c)
        with mp.workdps(40):
            val = 2 * xixi(mp.mpf(1) / 2 + cm)
        # kappa_1 from the xi'/xi limit (exact to 1e-23), since the remainder bound of the identity route is an upper bound only
        lo, hi = 4 * kappa1_xi * cm - 4 * kappa2 * cm ** 3, 4 * kappa1_xi * cm
        okc = (lo <= val <= hi)
        ok_sand = ok_sand and okc
        xs.append({'c': c, 'alpha': ms(mp.mpf(1) / 2 + cm), 'two_xi_prime_over_xi': ms(val, 15), 'lower_4k1c_minus_4k2c3': ms(lo, 15), 'upper_4k1c': ms(hi, 15), 'inside': bool(okc)})
    log('sandwich of M11 holds at c = 0.01, 0.1, 0.25, 0.4:', ok_sand)
    claim('208', '2 xi\'/xi(alpha) = 4 kappa_1 (alpha minus 1/2) + O(c^3), sandwich 4k1c minus 4k2c^3 <= . <= 4k1c', 'as printed', 'holds at c = 0.01, 0.1, 0.25, 0.4' if ok_sand else 'fails', ok_sand)
    cross = []
    for Lv, lam, who in (('0.5', '1.05e-6', 'Zhu upper bound'), ('0.5', '9.743e-7', 'chapter 26 mode value'), ('0.8', '2.27e-17', 'Zhu upper bound'),
                         ('0.8', '2.144e-17', 'chapter 26 mode value'), ('0.8', '8.9e-18', 'Zhu certified lower bound')):
        x = mp.mpf(lam) / (4 * kappa1)
        cross.append({'L': Lv, 'lambda': lam, 'source': who, 'lambda_over_4kappa1': mp.nstr(x, 6), 'two_digits': sig(x, 2)})
        log('crossover L', Lv, lam, who, '->', mp.nstr(x, 5))
    res['crossovers'] = cross
    byw = {(c['L'], c['source']): c for c in cross}
    claim('212', 'crossover at L = 0.5 from Zhu upper bound 1.05e-6', '1.1e-5', byw[('0.5', 'Zhu upper bound')]['two_digits'], byw[('0.5', 'Zhu upper bound')]['two_digits'] == '1.1e-05')
    claim('212', 'crossover at L = 0.5 from the chapter value 9.743e-7', '1.1e-5', byw[('0.5', 'chapter 26 mode value')]['two_digits'],
          byw[('0.5', 'chapter 26 mode value')]['two_digits'] == '1.1e-05', 'rounds to 1.1e-5 only because 1.054e-5 is rounded up')
    claim('212', 'crossover at L = 0.8 from Zhu upper bound 2.27e-17', '2.5e-16', byw[('0.8', 'Zhu upper bound')]['two_digits'], byw[('0.8', 'Zhu upper bound')]['two_digits'] == '2.5e-16')
    claim('212', 'crossover at L = 0.8 from the chapter value 2.144e-17', '2.5e-16', byw[('0.8', 'chapter 26 mode value')]['two_digits'],
          byw[('0.8', 'chapter 26 mode value')]['two_digits'] == '2.5e-16', 'the printed 2.5e-16 uses Zhu upper bound, not the chapter value')
    res['kappa'] = {'kappa_1': ms(kappa1, 14), 'kappa_1_xi_limit': mp.nstr(kappa1_xi, 14), 'kappa_2': ms(kappa2, 12), 'tail_bound_kappa_1': mp.nstr(t4 / 4, 3),
                    'identity_sum_gamma_over_gamma2_plus_quarter': ms(A1, 14), 'sandwich': xs}

    # (6) Table zh:tab:det columns 1 and 2, offsets
    mids = [(zd[j] + zd[j + 1]) / 2 for j in range(ng - 1)]
    det = []
    for H, (pcol, pl) in DET.items():
        jm = min(range(ng - 1), key=lambda j: abs(mids[j] - float(H)))
        mid = (gam[jm] + gam[jm + 1]) / 2
        half = mp.log(mid / (2 * mp.pi)) / 2
        hv = float(H)
        jc = max(j for j in range(ng - 1) if zd[j] <= hv)
        det.append({'target': H, 'pair_index_1based': [jm + 1, jm + 2], 'gamma_left': ms(gam[jm], 12), 'gamma_right': ms(gam[jm + 1], 12), 'midpoint': ms(mid, 12),
                    'half_log_over_2pi': ms(half, 10), 'gap': ms(gam[jm + 1] - gam[jm], 6), 'gap_containing_target_starts_at_index': jc + 1,
                    'mean_spacing_2pi_over_log': ms(2 * mp.pi / mp.log(mid / (2 * mp.pi)), 6),
                    'offsets_L_det_minus_half_log': [ms(mp.mpf(x) - half, 5) for x in pl]})
        claim('127 to 130', 'Table zh:tab:det midpoint (column 1) at target %s' % H, H, rnd(mid, 2), rnd(mid, 2) == H, 'pair %d,%d' % (jm + 1, jm + 2))
        claim('127 to 130', 'Table zh:tab:det (1/2) log(gamma/2 pi) (column 2) at %s' % H, pcol, rnd(half, 3), rnd(half, 3) == pcol)
        log('det pair', jm + 1, jm + 2, 'mid', ms(mid, 10), 'half log', ms(half, 8), 'gap', ms(gam[jm + 1] - gam[jm], 5), 'offsets', [ms(mp.mpf(x) - half, 4) for x in pl])
    res['det_columns_1_2'] = det
    small = []
    for d in det:
        for k in range(3):
            small.append(mp.mpf(d['offsets_L_det_minus_half_log'][k]))
    lo_off, hi_off = min(small), max(small)
    big = mp.mpf(det[3]['offsets_L_det_minus_half_log'][3])
    claim('144', 'offsets of the 12 thresholds with delta <= 0.2 above (1/2) log(gamma/2 pi), "between 0.02 and 0.17"', '0.02 to 0.17', '%s to %s' % (rnd(lo_off, 3), rnd(hi_off, 3)),
          rnd(lo_off, 2) == '0.02' and rnd(hi_off, 2) == '0.17')
    claim('144', 'large displacement seen earlier "by up to 0.43 at height 2000" (delta = 0.45)', '0.43', rnd(-big, 3), rnd(-big, 2) == '0.43')
    gap2000 = mp.mpf(det[3]['gap'])
    claim('144', 'on line neighbours "spaced 0.9 apart" at height 2000', '0.9', rnd(gap2000, 3), rnd(gap2000, 1) == '0.9', 'the gap of the pair; the mean spacing is %s' % det[3]['mean_spacing_2pi_over_log'])
    amp = mp.exp(2 * mp.mpf('2.45') * mp.mpf('0.45'))
    claim('144', 'amplification e^(2 L delta) at L = 2.45, delta = 0.45 (about 9.1)', 'about e^(2 L delta)', rnd(amp, 2), True, 'an upper bound of Phragmen Lindelof type, not a measured gain')
    res['det_text_numbers'] = {'offset_range_delta_le_0.2': [ms(lo_off, 5), ms(hi_off, 5)], 'offset_delta_0.45_at_2000': ms(big, 5), 'e_2L_delta': ms(amp, 6)}

    res['claims'] = CLAIMS
    bad = [c for c in CLAIMS if c['status'] != 'agrees']
    for c in CLAIMS:
        log('CLAIM %-14s %-8s %s | chapter %s | mine %s %s' % (c['line'], c['status'], c['claim'], c['chapter'], c['mine'], c['note']))
    log('claims: %d, agreeing %d, differing %d' % (len(CLAIMS), len(CLAIMS) - len(bad), len(bad)))
    if args.out:
        out = Path(args.out)
    elif args.quick:
        out = Path(tempfile.gettempdir()) / 'zh_comb_quick.json'
    else:
        out = Path('data/zh_comb.json') if Path('data').is_dir() else Path(__file__).resolve().parent / 'zh_comb.json'
    out.write_text(json.dumps(res, indent=1))
    log('wrote', out)


if __name__ == '__main__':
    main()
