# polya_qrh.py: Table pol:tab:qrh, Lemma pol:lem:window, Theorem pol:thm:qrhcomp and the threshold alpha_* of Section pol:sec:qrh of
# Chapter 16 (Polya and the zeros of entire functions). Added for the book in October 2026, written from the mathematics (the
# chapter author's own scripts were not read, run or copied).
# What it checks. (1) The column B_alpha = (xi'/xi)(alpha) of the table (alpha = 0.9, 0.925, 0.95, 1, 1.5, 2, 5, 7.8; B_1 = 1 + gamma/2 - (log pi)/2 - log 2),
# mpmath at 30 digits. (2) Lemma window: nu(r, d) = (theta(r + d) - theta(r - d))/pi - s(r + d) - s(r - d), r > d + e, in two variants: T, the
# chapter's s(t) = 0.112 log t + 0.278 log log t + 2.510 (Trudgian, which the book rejects as incompletely proved), and H, the bound of Hasanalizade,
# Shen and Wong used by the book: s(t) = 2.5167 for e <= t <= X_0 = 30 610 046 000 and s(t) = 0.1095 log t + 0.2042 log log t + 3.0305 for t > X_0
# (variant H has two downward jumps in nu, at r = X_0 - d and r = X_0 + d, which the grids contain). theta is the Riemann Siegel theta function, evaluated
# in double precision through the exact difference theta(r + d) - theta(r - d) = (1/2)(r log((r + d)/(r - d)) + d log((r^2 - d^2)/(4 pi^2))) - d + c(r + d) - c(r - d)
# with the Stirling series c(t) = sum_{k = 1..7} (1 - 2^(1 - 2k)) |B_2k| / (4k(2k - 1)) t^-(2k-1) (valid for t >= 30, tested against mpmath siegeltheta at
# 50 digits, including r = 1e12), so that no digits are lost at height 1e12 to 1e15. A sanity test counts the actual zeros of data/zeros_1700.txt in windows
# [r - d, r + d] up to height 2300 against nu (the lower bound must not exceed the count). (3) Theorem qrhcomp: c_min = alpha - 7/8, c_max = alpha - 1/8
# (off line) or c_min = c_max = alpha - 1/2 (on line, only below the verified height H = 3e12), f(c, d) = c/(c^2 + d^2), m(d) = min(f(c_min, d), f(c_max, d)),
# k = floor(B_alpha/m(d)) + 1. T_0(alpha, d) is the supremum of the set of heights r (r >= d + 30, r <= 1e15 off line, r <= H - d on line) at which the
# requirement fails; it is found on a grid of 1500 points per decade (plus the jump points of variant H) by a suffix minimum of nu and refined by bisection to
# machine precision; the requirement holds for every r above T_0 on the grid. THREE definitions of the requirement are computed and written, because the chapter's
# numbers need the second: cont (the theorem as stated: nu(r, d) m(d) > B_alpha, strict), kzero (the window contains at least k zeros: nu(r, d) >= k), int (the
# count of zeros is an integer at least nu, so it is at least k exactly when nu > k - 1; this is the sharpest valid form). T_0(alpha) is the minimum over the
# half width d on the grid d = 0.5, 1, ..., 40 (ties go to the smaller d). (4) Verification, not assumption, of what the theorem needs for all r >= T_0: on a grid of
# 20000 points per decade from T_0 to 1e15 (or H - d) the requirement holds at every point (including one sided limits at the jumps of variant H) and d nu/dr > 0 at
# every point; beyond 1e15 the derivative is at least d/(pi r) - 2 s'(r - d), positive when d > 2 pi (r/(r - d)) (c1 + c2/log(r - d)), which is checked at r = 1e15
# (the right side decreases in r). The value of nu at the reported T_0 is recomputed with mpmath at 50 digits. (5) alpha_*: the smallest alpha whose off line T_0 is
# at most H = 3e12, by a scan of 400 values of alpha in [0.8751, 0.96] (monotonicity of the predicate is tested) and bisection to 1e-12.
# Error control. Variant T reproduces (or not) the printed numbers, which validates the implementation of nu; the window lemma is tested against the real zeros;
# the mpmath recomputation of nu at T_0 is compared with the level. Nothing here is a proof: the grids are finite, and the Stirling series error is below 1e-15 at t >= 30.
# Output: data/polya_qrh.json (numbers as strings; the three definitions and both variants, the pasteable table in the chapter's form), the log is the standard output
# (data/polya_qrh.log). Usage: python3 polya_qrh.py [--quick] [--out FILE] [--workers W] [--zeros FILE]; --quick uses coarser grids (about 20 s, writes to /tmp unless
# --out is given); the full case is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8 -g polya_qrh.json code/polya_qrh.py data/zeros_1700.txt -- --out polya_qrh.json
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import mpmath as mp
import numpy as np

T0 = time.time()
PI = math.pi
X0 = 30610046000.0
H = 3e12
KS = 7
CK = [(1 - 2.0 ** (1 - 2 * k)) * abs(float(mp.bernoulli(2 * k))) / (4 * k * (2 * k - 1)) for k in range(1, KS + 1)]
ALPHAS = [0.9, 0.925, 0.95, 1.0, 1.5, 2.0, 5.0, 7.8]
DEFS = ['cont', 'kzero', 'int']
EXTRA = [0.8885, 0.8887, 0.8888, 0.8889, 0.889]   # off line rows around the chapter's alpha_* = 0.8888
CHAPTER = {
    0.9: ('0.01848', ('5.2e8', 4, 12), ('270', 14.5, 10)),
    0.925: ('0.01963', ('1.4e6', 4.5, 8), None),
    0.95: ('0.02079', ('1.4e5', 6, 10), ('270', 14.5, 10)),
    1.0: ('0.02310', ('1.8e4', 6.5, 8), ('270', 14.5, 10)),
    1.5: ('0.04614', ('721', 11, 9), None),
    2.0: ('0.06907', ('464', 12, 9), ('270', 14.5, 10)),
    5.0: ('0.20166', ('351', 12, 8), ('299', 12.5, 8)),
    7.8: ('0.31322', ('420', 13, 10), ('365', 16.5, 14)),
}


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def default_out(name):
    return os.path.join('data', name) if os.path.isdir('data') else name


def find_zeros(path):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in ([path] if path else []) + ['data/zeros_1700.txt', 'zeros_1700.txt', os.path.join(here, 'zeros_1700.txt'), os.path.join(here, '..', 'data', 'zeros_1700.txt')]:
        if os.path.exists(c):
            return c
    raise FileNotFoundError('zeros_1700.txt')


# ---------------------------------------------------------------- xi'/xi
def xi_logder(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)


def B_val(a):
    a = mp.mpf(a)
    if a == 1:
        return 1 + mp.euler / 2 - mp.log(mp.pi) / 2 - mp.log(2)
    return mp.re(xi_logder(a))


# ---------------------------------------------------------------- the window count (double precision, vectorised)
def c_theta(t):
    return sum(CK[k - 1] / t ** (2 * k - 1) for k in range(1, KS + 1))


def dc_theta(t):
    return -sum((2 * k - 1) * CK[k - 1] / t ** (2 * k) for k in range(1, KS + 1))


def theta_diff(r, d):
    r = np.asarray(r, dtype=np.float64)
    L = 2.0 * np.arctanh(d / r)
    M = np.log((r - d) * (r + d)) - 2.0 * math.log(2 * PI)
    return 0.5 * (r * L + d * M) - d + c_theta(r + d) - c_theta(r - d)


def dtheta_diff(r, d):
    r = np.asarray(r, dtype=np.float64)
    return np.arctanh(d / r) + dc_theta(r + d) - dc_theta(r - d)


def s_fun(var, t):
    t = np.asarray(t, dtype=np.float64)
    if var == 'N':
        return np.zeros_like(t)
    lt = np.log(t)
    if var == 'T':
        return 0.112 * lt + 0.278 * np.log(lt) + 2.510
    return np.where(t <= X0, 2.5167, 0.1095 * lt + 0.2042 * np.log(lt) + 3.0305)


def ds_fun(var, t):
    t = np.asarray(t, dtype=np.float64)
    if var == 'N':
        return np.zeros_like(t)
    lt = np.log(t)
    if var == 'T':
        return (0.112 + 0.278 / lt) / t
    return np.where(t <= X0, 0.0, (0.1095 + 0.2042 / lt) / t)


def nu_fun(var, r, d):
    r = np.asarray(r, dtype=np.float64)
    return theta_diff(r, d) / PI - s_fun(var, r + d) - s_fun(var, r - d)


def dnu_fun(var, r, d):
    r = np.asarray(r, dtype=np.float64)
    return dtheta_diff(r, d) / PI - ds_fun(var, r + d) - ds_fun(var, r - d)


def nu_mp(var, r, d):
    mp.mp.dps = 50
    r, d = mp.mpf(r), mp.mpf(d)

    def s(t):
        lt = mp.log(t)
        if var == 'T':
            return mp.mpf('0.112') * lt + mp.mpf('0.278') * mp.log(lt) + mp.mpf('2.510')
        if t <= X0:
            return mp.mpf('2.5167')
        return mp.mpf('0.1095') * lt + mp.mpf('0.2042') * mp.log(lt) + mp.mpf('3.0305')
    return (mp.siegeltheta(r + d) - mp.siegeltheta(r - d)) / mp.pi - s(r + d) - s(r - d)


def f_lor(c, d):
    return c / (c * c + d * d)


def m_of_d(alpha, d, online):
    if online:
        c = alpha - 0.5
        return f_lor(c, d)
    return min(f_lor(alpha - 0.875, d), f_lor(alpha - 0.125, d))


# ---------------------------------------------------------------- grids and T_0
def build_store(var, DS, per_decade):
    glob = 10.0 ** np.linspace(math.log10(31.0), 15.0, int(per_decade * (15 - math.log10(31.0))) + 1)
    store = {}
    for d in DS:
        r = glob[glob >= d + 30.0]
        if var == 'H':
            ex = [b + o for b in (X0 - d, X0 + d) for o in (-1.0, -0.01, 0.0, 0.01, 1.0)]
            r = np.unique(np.concatenate([r, np.array(ex)]))
        nu = nu_fun(var, r, d)
        sm_off = np.minimum.accumulate(nu[::-1])[::-1]
        sm_on = np.minimum.accumulate(np.where(r <= H - d, nu, np.inf)[::-1])[::-1]
        store[d] = (r, sm_off, sm_on)
    return store


def level_of(defname, B, m):
    k = math.floor(B / m) + 1
    if defname == 'cont':
        return B / m, True, k
    if defname == 'kzero':
        return float(k), False, k
    return float(k - 1), True, k


def t0_search(var, d, lev, strict, online, store, refine=True):
    r, sm_off, sm_on = store[d]
    sm = sm_on if online else sm_off
    i0 = int(np.searchsorted(sm, lev, side='right' if strict else 'left'))
    if i0 >= len(sm):
        return None
    if online and r[i0] > H - d:
        return None
    if i0 == 0:
        return float(r[0]), True
    lo, hi = float(r[i0 - 1]), float(r[i0])
    if not refine:
        return hi, False
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        v = float(nu_fun(var, mid, d))
        if (v > lev) if strict else (v >= lev):
            hi = mid
        else:
            lo = mid
    return hi, False


def best_over_d(var, alpha, B, online, defname, store, DS):
    cands = []
    for d in DS:
        m = m_of_d(alpha, d, online)
        lev, strict, k = level_of(defname, B, m)
        res = t0_search(var, d, lev, strict, online, store, refine=False)
        if res is not None:
            cands.append((res[0], d, lev, strict, k, m, res[1]))
    if not cands:
        return None
    cands.sort(key=lambda x: x[0])
    best = None
    for c in cands[:4]:
        res = t0_search(var, c[1], c[2], c[3], online, store, refine=True)
        t0v = res[0]
        if best is None or t0v < best['T0'] * (1 - 1e-12):
            best = {'T0': t0v, 'd': c[1], 'k': c[4], 'level': c[2], 'strict': c[3], 'm': c[5], 'B_over_m': B / c[5], 'at_grid_start': res[1]}
    return best


def fmt_T0(x):
    if x is None:
        return None
    if x < 1000:
        return '%d' % int(round(x))
    mant, ex = ('%.1e' % x).split('e')
    return '%s\\times10^{%d}' % (mant, int(ex))


def fmt_chapter_style(x):
    if x is None:
        return None
    if x < 1000:
        return '%d' % int(round(x))
    mant, ex = ('%.1e' % x).split('e')
    return '%se%d' % (mant, int(ex))


def verify_case(args):
    var, online, d, lev, strict, T0v, npd = args
    Rmax = H - d if online else 1e15
    start = T0v * (1 + 1e-9)
    if start >= Rmax:
        return {'trivial': True}
    ndec = math.log10(Rmax / start)
    rr = start * 10.0 ** np.linspace(0.0, ndec, max(int(npd * ndec), 2000))
    if var == 'H':
        ex = np.array([b + o for b in (X0 - d, X0 + d) for o in (-1.0, -0.01, 0.0, 0.01, 1.0)])
        ex = ex[(ex > start) & (ex < Rmax)]
        rr = np.sort(np.concatenate([rr, ex]))
    nu = nu_fun(var, rr, d)
    ok = (nu > lev) if strict else (nu >= lev)
    dnu = dnu_fun(var, rr, d)
    out = {'points': int(len(rr)), 'all_satisfy_requirement': bool(ok.all()), 'min_nu_minus_level': '%.6g' % float((nu - lev).min()),
           'min_dnu_dr': '%.6g' % float(dnu.min()), 'dnu_positive_everywhere': bool((dnu > 0).all())}
    if not online:
        c1, c2 = (0.112, 0.278) if var == 'T' else (0.1095, 0.2042)
        R = Rmax
        out['tail_certificate_beyond_1e15'] = bool(d / PI > 2.0 * (R / (R - d)) * (c1 + c2 / math.log(R - d)))
    out['verified'] = bool(out['all_satisfy_requirement'] and out['dnu_positive_everywhere'] and out.get('tail_certificate_beyond_1e15', True))
    return out


def tmin_for(var, alpha, B, defname, store, DS):
    b = best_over_d(var, alpha, B, False, defname, store, DS)
    return (b['T0'], b) if b else (math.inf, None)


def alpha_star(var, defname, store, DS, npts, iters):
    def pred(a):
        B = float(B_val(a))
        t, b = tmin_for(var, a, B, defname, store, DS)
        return t <= H, t, b
    grid = [0.8751 + (0.96 - 0.8751) * i / (npts - 1) for i in range(npts)]
    flags = [pred(a)[0] for a in grid]
    first_true = flags.index(True) if True in flags else None
    monotone = first_true is not None and all(flags[first_true:]) and not any(flags[:first_true])
    if first_true is None or first_true == 0:
        return {'found': False}
    lo, hi = grid[first_true - 1], grid[first_true]
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if pred(mid)[0]:
            hi = mid
        else:
            lo = mid
    ok, t, b = pred(hi)
    up = math.ceil(hi * 1e4 - 1e-9) / 1e4
    ok4, t4, b4 = pred(up)
    return {'found': True, 'monotone_predicate_on_scan': monotone, 'alpha_star': '%.10f' % hi, 'false_side': '%.10f' % lo,
            'T0_at_alpha_star': '%.6g' % t, 'd': b['d'], 'k': b['k'], 'B_over_m': '%.6f' % b['B_over_m'],
            'rounded_up_to_4_decimals': {'alpha': '%.4f' % up, 'T0': '%.6g' % t4, 'd': b4['d'], 'k': b4['k'], 'B_over_m': '%.6f' % b4['B_over_m'], 'T0_le_H': bool(ok4)}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--workers', type=int)
    ap.add_argument('--zeros')
    args = ap.parse_args()
    quick = args.quick
    workers = 1 if quick else (args.workers or min(8, os.cpu_count() or 1))
    out_path = args.out or ('/tmp/polya_qrh_quick.json' if quick else default_out('polya_qrh.json'))
    if workers > 1:
        pool = get_context('fork').Pool(workers)

        def pmap(f, tasks):
            return pool.map(f, tasks, chunksize=1)
    else:
        def pmap(f, tasks):
            return [f(t) for t in tasks]
    zeros = [float(x) for x in open(find_zeros(args.zeros)).read().split()]
    DS = [0.5 * i for i in range(1, 81)] if not quick else [0.5 * i for i in range(1, 41)]
    per_decade = 1500 if not quick else 300
    npd_verify = 20000 if not quick else 1500
    result = {'script': 'polya_qrh.py', 'quick': quick, 'workers': workers,
              'definitions': {
                  'cont': 'nu(r, d) m(d) > B_alpha strict, as in the statement of Theorem qrhcomp; T_0 = sup of the failing heights',
                  'kzero': 'nu(r, d) >= k with k = floor(B_alpha/m(d)) + 1 (the window holds at least k zeros); T_0 = sup of the failing heights',
                  'int': 'nu(r, d) > k - 1 (the zero count is an integer at least nu, so at least k); T_0 = sup of the failing heights'},
              'grid': {'d': 'step 0.5 from 0.5 to %g' % DS[-1], 'r_points_per_decade': per_decade, 'r_range': 'from d + 30 to 1e15 (off line), to H - d (on line)',
                       'verification_points_per_decade': npd_verify}}

    # ---------------- self tests of the window count
    mp.mp.dps = 50
    tests = []
    for r, d in [(300.0, 14.5), (1000.0, 6.5), (1e6, 4.0), (1e12, 3.5), (1.45e12, 3.5)]:
        exact = (mp.siegeltheta(mp.mpf(r) + d) - mp.siegeltheta(mp.mpf(r) - d))
        ours = float(theta_diff(r, d))
        tests.append({'r': r, 'd': d, 'exact': mp.nstr(exact, 17), 'ours': '%.17g' % ours, 'abs_error': '%.3g' % abs(float(exact) - ours)})
    result['theta_difference_self_test'] = tests
    log('theta difference self test, largest abs error %s' % max(t['abs_error'] for t in tests))
    result['window_vs_real_zeros'] = {}
    zz = np.array(zeros)
    for d in ([3.5, 14.5] if quick else [3.5, 4.0, 6.0, 6.5, 11.0, 12.0, 13.0, 14.5, 16.5]):
        r = np.arange(d + 30.0, zz[-1] - d, 0.05)
        cnt = np.searchsorted(zz, r + d, side='right') - np.searchsorted(zz, r - d, side='left')
        nT, nH, nN = nu_fun('T', r, d), nu_fun('H', r, d), nu_fun('N', r, d)
        result['window_vs_real_zeros'][str(d)] = {
            'heights': '[%g, %g]' % (r[0], r[-1]), 'min_count': int(cnt.min()), 'min_count_minus_nu_T': '%.4f' % float((cnt - nT).min()),
            'min_count_minus_nu_H': '%.4f' % float((cnt - nH).min()), 'mean_count_minus_theta_difference_over_pi': '%.4f' % float((cnt - nN).mean()),
            'max_abs_count_minus_theta_difference_over_pi': '%.4f' % float(np.abs(cnt - nN).max())}
    log('window against the real zeros (d: min count minus nu_T, nu_H):', {d: (v['min_count_minus_nu_T'], v['min_count_minus_nu_H']) for d, v in result['window_vs_real_zeros'].items()})

    # ---------------- B column
    mp.mp.dps = 30
    Bs = {a: B_val(a) for a in ALPHAS + [0.8888]}
    mp.mp.dps = 80
    b1_lim = mp.re(xi_logder(mp.mpf(1) + mp.mpf('1e-30')))
    mp.mp.dps = 30
    result['B_column'] = {'%g' % a: {'B': mp.nstr(Bs[a], 15), 'printed_5_decimals': '%.5f' % float(Bs[a]), 'chapter': CHAPTER[a][0] if a in CHAPTER else None,
                                     'agrees': ('%.5f' % float(Bs[a]) == CHAPTER[a][0]) if a in CHAPTER else None} for a in Bs}
    result['B_column']['B_1_limit_check'] = mp.nstr(b1_lim, 25)
    extra = {}
    for a in ['0.9', '0.95']:
        extra[a] = '%.5f' % float(B_val(a))
    result['B_column']['line_672_values_at_0.9_0.95_1'] = {'0.9': extra['0.9'], '0.95': extra['0.95'], '1': '%.5f' % float(Bs[1.0]), 'chapter': ['0.01848', '0.02079', '0.02310']}
    log('B column:', {k: v['printed_5_decimals'] for k, v in result['B_column'].items() if isinstance(v, dict) and 'printed_5_decimals' in v})

    # ---------------- tables
    result['table'] = {}
    result['pasteable_table'] = {}
    result['alpha_star'] = {}
    result['verification'] = {}
    result['nu_at_T0_mpmath_check'] = {}
    vtasks, vkeys = [], []
    for var in ['T', 'H']:
        log('variant %s: building the grids' % var)
        store = build_store(var, DS, per_decade)
        result['table'][var] = {}
        result['pasteable_table'][var] = {}
        for dn in DEFS:
            rows, lines = {}, []
            for a in ALPHAS:
                B = float(Bs[a])
                row = {'B': '%.5f' % B}
                for online in (False, True):
                    b = best_over_d(var, a, B, online, dn, store, DS)
                    key = 'on_line' if online else 'off_line'
                    if b is None:
                        row[key] = None
                        continue
                    row[key] = {'T0': '%.6g' % b['T0'], 'd': b['d'], 'k': b['k'], 'm': '%.8g' % b['m'], 'B_over_m': '%.5f' % b['B_over_m'], 'at_grid_start': b['at_grid_start']}
                    vtasks.append((var, online, b['d'], b['level'], b['strict'], b['T0'], npd_verify))
                    vkeys.append((var, dn, a, key))
                    nuv = float(nu_mp(var, b['T0'], b['d']))
                    result['nu_at_T0_mpmath_check'].setdefault(var, {}).setdefault(dn, {})['%g_%s' % (a, key)] = {'nu_mp': '%.9f' % nuv, 'level': '%.9f' % b['level'], 'minus': '%.3g' % (nuv - b['level'])}
                ch = CHAPTER.get(a)
                if ch is not None and var == 'T':
                    for key, idx in (('off_line', 1), ('on_line', 2)):
                        if ch[idx] is not None and row.get(key) is not None:
                            t_ch, d_ch, k_ch = ch[idx]
                            row[key]['chapter'] = {'T0': t_ch, 'd': d_ch, 'k': k_ch}
                            row[key]['T0_rounds_to_chapter'] = fmt_chapter_style(float(row[key]['T0'])) == t_ch
                            row[key]['d_k_agree'] = (row[key]['d'] == d_ch and row[key]['k'] == k_ch)
                if ch is not None and var == 'T' and dn == 'kzero':
                    for key, idx, onl in (('off_line', 1, False), ('on_line', 2, True)):
                        if ch[idx] is None or row.get(key) is None:
                            continue
                        t_ch, d_ch, k_ch = ch[idx]
                        m_ch = m_of_d(a, float(d_ch), onl)
                        lev_ch, strict_ch, k_own = level_of('kzero', B, m_ch)
                        res_ch = t0_search('T', float(d_ch), lev_ch, strict_ch, onl, store)
                        row[key]['at_chapter_d'] = {'d': d_ch, 'k_from_formula': k_own, 'k_chapter': k_ch, 'k_agrees': bool(k_own == k_ch),
                                                    'T0': ('%.6g' % res_ch[0]) if res_ch else None, 'T0_chapter': t_ch,
                                                    'chapter_over_ours': ('%.4f' % (float(t_ch) / res_ch[0])) if res_ch else None}
                rows['%g' % a] = row
                off, on = row.get('off_line'), row.get('on_line')
                lines.append('$%.3f$ & $%s$ & $%s$ & $(%g,%d)$ & %s & %s\\\\' % (
                    a, row['B'], fmt_T0(float(off['T0'])) if off else '--', off['d'] if off else 0, off['k'] if off else 0,
                    ('$%s$' % fmt_T0(float(on['T0']))) if on else '', ('$(%g,%d)$' % (on['d'], on['k'])) if on else ''))
            for a in EXTRA:
                mp.mp.dps = 30
                B = float(B_val(a))
                b = best_over_d(var, a, B, False, dn, store, DS)
                if b is None:
                    rows['%g' % a] = {'B': '%.5f' % B, 'off_line': None}
                    continue
                rows['%g' % a] = {'B': '%.5f' % B, 'off_line': {'T0': '%.6g' % b['T0'], 'd': b['d'], 'k': b['k'], 'm': '%.8g' % b['m'], 'B_over_m': '%.5f' % b['B_over_m'], 'at_grid_start': b['at_grid_start']}}
                vtasks.append((var, False, b['d'], b['level'], b['strict'], b['T0'], npd_verify))
                vkeys.append((var, dn, a, 'off_line'))
            result['table'][var][dn] = rows
            result['pasteable_table'][var][dn] = lines
            if var == 'T':
                nmatch = sum(1 for a in ALPHAS for key in ('off_line', 'on_line')
                             if rows['%g' % a].get(key) and rows['%g' % a][key].get('T0_rounds_to_chapter'))
                ndk = sum(1 for a in ALPHAS for key in ('off_line', 'on_line')
                          if rows['%g' % a].get(key) and rows['%g' % a][key].get('d_k_agree'))
                nall = sum(1 for a in ALPHAS for key in ('off_line', 'on_line') if CHAPTER[a][1 if key == 'off_line' else 2] is not None)
                log('variant T, definition %s: T_0 rounds to the chapter in %d of %d cases, (d, k) agree in %d of %d' % (dn, nmatch, nall, ndk, nall))
                if dn == 'kzero':
                    own = {'%g_%s' % (a, key): rows['%g' % a][key]['at_chapter_d'] for a in ALPHAS for key in ('off_line', 'on_line')
                           if rows['%g' % a].get(key) and 'at_chapter_d' in rows['%g' % a][key]}
                    log('variant T, kzero, at the chapter\'s own d: k from the formula agrees in %d of %d; chapter T_0 over ours: %s' % (
                        sum(1 for v in own.values() if v['k_agrees']), len(own), {k: v['chapter_over_ours'] for k, v in own.items()}))
                    result['table']['T']['kzero_at_chapter_d_summary'] = own
                result['table'][var][dn]['_summary'] = {'T0_match': nmatch, 'dk_match': ndk, 'cases': nall}
            else:
                log('variant H, definition %s: off line T_0 %s' % (dn, {a: rows['%g' % a]['off_line']['T0'] if rows['%g' % a]['off_line'] else None for a in ALPHAS}))
        log('variant %s: alpha_* scans' % var)
        result['alpha_star'][var] = {}
        for dn in DEFS:
            als = alpha_star(var, dn, store, DS, 60 if quick else 400, 25 if quick else 55)
            result['alpha_star'][var][dn] = als
            if als.get('found'):
                log('alpha_* (%s, %s) = %s: T_0 %s, d %s, k %s; rounded up %s' % (var, dn, als['alpha_star'], als['T0_at_alpha_star'], als['d'], als['k'], als['rounded_up_to_4_decimals']))
            else:
                log('alpha_* (%s, %s) not found in [0.8751, 0.96]' % (var, dn))
        del store
    log('verifying monotonicity and the requirement beyond T_0 for %d reported cases' % len(vtasks))
    vres = pmap(verify_case, vtasks)
    for key, v in zip(vkeys, vres):
        result['verification'].setdefault(key[0], {}).setdefault(key[1], {})['%g_%s' % (key[2], key[3])] = v
    allv = all(v.get('verified', True) for v in vres)
    result['verification_all_verified'] = bool(allv)
    log('verification of all reported T_0: %s (%d cases)' % (allv, len(vres)))
    result['seconds'] = round(time.time() - T0, 1)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=1)
    log('written', out_path)


if __name__ == '__main__':
    main()
