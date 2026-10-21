# yo_gaps.py: the gap lemma of Chapter 17 (Yoshida's Hermitian Forms, Weil Positivity and the Thorin Program), Lemma yo:lem:gap
# ("consecutive ordinates of the zeros of zeta differ by less than 8.3"), and its use in the proof of Theorem yo:thm:coercive (ii)
# (the margin 0.057233 minus 0.010010 minus 0.046210 at c = 1/2), lines 276 to 295 of the chapter source. Added for the book in
# October 2026, written from the mathematics (the chapter's own gaps.py was never read or run).
# What is checked, in order. (A) The ordinates of data/zeros_1700.txt (68 digits, to height 2197.26) as balls: every gap, the first
# gap gamma_2 minus gamma_1, the largest gap above gamma_2, and the summary of the certified 1800 zeros (to height 2304.3) that the
# book's v1min.py wrote to riesz_program/data/v1min.json (the individual ordinates 1701 to 1800 are not stored anywhere in the
# repository, so those are quoted, not recomputed); also the sawtooth S(T) = N(T) minus theta(T)/pi minus 1 at both sides of the 1700
# ordinates (its extreme values, to compare with the bounds below). (B) The lemma in the chapter's form, with the Trudgian bound
# 0.112 log T + 0.278 log log T + 2.510: the two sides h theta'(T)/pi and 2 B(T + h) at T = 2524 and h = 8.3 (exact theta' from the
# digamma function), the threshold in T, the same with the published extra term 0.2/T and with the bounds taken at each endpoint,
# the height at which the inequality fails (below the threshold, in particular at gamma_1500 = 1980.91 and at the end 2304.3 of the
# certified zeros), and the step h that Trudgian's bound would need there. (C) The version with the explicit bounds of Hasanalizade,
# Shen and Wong (HSW22, quoted in the proof of Proposition bp:prop:v1min): |S(T)| <= 2.5167 for 0 <= T <= X0 = 30 610 046 000 and
# |S(T)| <= 0.1095 log T + 0.2042 log log T + 3.0305 for T >= e; the threshold in T for h = 8.3, 6.8874 and 5.36, and the smallest
# step h that the constant bound gives on [2304.3, X0]. (D) Above X0: with D_h(T) = h theta'(T)/pi minus 2 B_2(T + h), B_2 the second
# HSW22 bound, the lemma holds for every T >= X0 minus h if D_h is positive at X0 minus h and increasing. The script gives the
# value at X0 (the review quotes 29.4 against 12.6 for h = 8.3), the derivative D_h' in closed form with the exact theta'', the
# explicit lower bound for it, and a geometric grid check of D_h and D_h' in ball arithmetic as far as 10^3000 (the chapter's grid
# stopped at 10^30); a grid proves nothing between its points, the monotonicity does. (E) The conclusion: every gap is at most
# gamma_2 minus gamma_1 = 6.8873144970..., hence below 6.8874, and the gaps above gamma_2 are at most 5.4141. (F) The algebra of the
# proof of Theorem yo:thm:coercive (ii): the equivalence P_gamma(tau) >= P_gamma(0) iff tau^2 <= 3 gamma^2 minus c^2 (exact rational
# spot checks), kappa_1 = (1/2)(log xi)''(1/2) from its closed form (mpmath, 40 digits) against the zeros and a Riemann von Mangoldt
# tail, and the margins at c = 1/2 with the distance 4.15 = 8.3/2 of the chapter and 3.4437 = 6.8874/2 of the corrected lemma.
# Method and error control. Everything in (A) to (E) is python-flint ball arithmetic at 320 bits (192 with --quick); the ordinates
# are balls with the radius 10^(minus the number of decimals) of the file; a comparison is accepted only if it holds for the whole
# ball. theta' = (1/2) Re psi(1/4 + iT/2) minus (1/2) log(pi), theta'' = minus (1/4) Im psi'(1/4 + iT/2) (positive for every T > 0,
# since psi'(z) is the sum of (z + n)^(minus 2)); both are checked against a difference quotient of theta (from log Gamma) and against
# mpmath's siegeltheta. The thresholds are found by bisection on exact binary rationals with the sign certified by the balls, and
# the monotonicity of the left minus right side (needed to pass from a threshold to every larger T) is shown by the derivative.
# Nothing here is a proof of the HSW22 or Trudgian bounds themselves: they are quoted from the book (basepoint.tex) and the chapter.
# Output: data/yo_gaps.json (balls as strings with their radius, parameters, the claims table with the chapter's printed digits
# against the computed ones) and data/yo_gaps.log (the log below, written next to the JSON).
# Usage: python3 yo_gaps.py [--quick] [--out FILE] [--zeros FILE] [--v1min FILE]; --quick uses 300 zeros, 192 bits and a short grid
# (about 5 s) and writes to /tmp unless --out is given; the full case takes under 30 s in one process.
import argparse
import json
import os
import random
import sys
import time
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from fractions import Fraction
from pathlib import Path

import flint
import mpmath as mp
from flint import acb, arb, ctx

T0 = time.time()
LOG = []
CLAIMS = []
CHECKS = {}
S1 = None   # |S(T)| <= 2.5167 for T <= X0 (HSW22), set in main
X0 = None


def log(*a):
    line = '[%6.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(line, flush=True)
    LOG.append(line)


def s(x, d=20):
    return x.str(d)


def sf(x, d=10):
    return x.str(d, radius=False)


def rnd(x, nd):
    getcontext().prec = 80
    m = Decimal(x.mid().str(50, radius=False))
    return str(m.quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_EVEN))


def claim(key, line, printed, value, nd):
    got = rnd(value, nd)
    ok = Decimal(got) == Decimal(printed)
    getcontext().prec = 80
    trunc = str(Decimal(value.mid().str(50, radius=False)).quantize(Decimal(1).scaleb(-nd), rounding='ROUND_DOWN'))
    note = '' if ok else ('the printed value is the truncation, not the rounding' if Decimal(trunc) == Decimal(printed) else '')
    CLAIMS.append({'claim': key, 'chapter_line': line, 'chapter': printed, 'computed': got, 'computed_ball': s(value, 25), 'agrees_to_printed_digits': ok, 'note': note})
    log('claim %-52s line %-4s chapter %-10s computed %-10s %s %s' % (key, line, printed, got, 'agrees' if ok else 'DIFFERS', note))


def check(name, cond):
    CHECKS[name] = bool(cond)
    log(('ok   ' if cond else 'FAIL ') + name)
    return cond


def find_file(arg, rels):
    here = Path(__file__).resolve().parent
    cands = ([Path(arg)] if arg else []) + [here.parent / r for r in rels] + [here / Path(r).name for r in rels] + [Path(Path(r).name) for r in rels]
    for c in cands:
        if c.is_file():
            return c
    sys.exit('input file not found: ' + ', '.join(rels))


def load_zeros(path, nmax):
    zs = []
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            d = len(t[0].split('.')[1])
            zs.append(arb(t[0] + ' +/- 1e-%d' % d))
            if len(zs) >= nmax:
                break
    return zs


# the Riemann Siegel theta function and its derivatives, exact (digamma, polygamma, log Gamma), T an arb
def theta(T):
    z = acb(arb(1) / 4, T / 2)
    return z.lgamma().imag - T / 2 * arb.pi().log()


def thetap(T):
    z = acb(arb(1) / 4, T / 2)
    return (z.digamma().real - arb.pi().log()) / 2


def thetapp(T):
    z = acb(arb(1) / 4, T / 2)
    return -(z.polygamma(1).imag) / 4


# bounds for S(T): Trudgian as printed in the chapter (published=False) or with the 0.2/T of the published theorem; HSW22 second bound
def b_tr(t, published=False):
    t = arb(t)
    L = t.log()
    v = arb('0.112') * L + arb('0.278') * L.log() + arb('2.510')
    return v + arb('0.2') / t if published else v


def b_tr_p(t, published=False):
    t = arb(t)
    L = t.log()
    v = arb('0.112') / t + arb('0.278') / (t * L)
    return v - arb('0.2') / (t * t) if published else v


def b2(t):
    t = arb(t)
    L = t.log()
    return arb('0.1095') * L + arb('0.2042') * L.log() + arb('3.0305')


def b2_p(t):
    t = arb(t)
    L = t.log()
    return arb('0.1095') / t + arb('0.2042') / (t * L)


def tr_sides(T, h, variant):
    T, h = arb(T), arb(h)
    lhs = h * thetap(T) / arb.pi()
    if variant == 'chapter':        # 2 max |S| with the bound taken at T + h, as the chapter writes it
        rhs = 2 * b_tr(T + h)
    elif variant == 'endpoints':    # |S(T)| + |S(T + h)| with each bound at its own point
        rhs = b_tr(T) + b_tr(T + h)
    else:                           # 'published': as chapter, with the extra 0.2/T of Trudgian's theorem
        rhs = 2 * b_tr(T + h, True)
    return lhs, rhs


def tr_diff(T, h, variant):
    lhs, rhs = tr_sides(T, h, variant)
    return lhs - rhs


def hsw_diff(T, h):
    return arb(h) * thetap(arb(T)) / arb.pi() - 2 * S1


def d_h(T, h):
    T, h = arb(T), arb(h)
    return h * thetap(T) / arb.pi() - 2 * b2(T + h)


def d_h_p(T, h):
    T, h = arb(T), arb(h)
    return h * thetapp(T) / arb.pi() - 2 * b2_p(T + h)


def bisect_true(pred, lo, hi, iters=200):
    # pred is certified true at hi and not certified at lo; returns a bracket (lo, hi) of floats, pred(hi) certified
    assert pred(hi) and not pred(lo), (lo, hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return lo, hi


def grid(T_from, T_to, ratio):
    T = arb(T_from)
    r = arb(ratio)
    end = arb(T_to)
    while True:
        yield T
        if T > end:
            return
        T = T * r


def main():
    global S1, X0
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--zeros')
    ap.add_argument('--v1min')
    args = ap.parse_args()
    quick = args.quick
    ctx.prec = 192 if quick else 320
    out_path = args.out or ('/tmp/yo_gaps_quick.json' if quick else str(Path(__file__).resolve().parent.parent / 'data' / 'yo_gaps.json'))
    zpath = find_file(args.zeros, ['data/zeros_1700.txt'])
    vpath = find_file(args.v1min, ['riesz_program/data/v1min.json'])
    S1 = arb('2.5167')
    X0 = arb(30610046000)
    pi = arb.pi()
    out = {'python_flint': flint.__version__, 'precision_bits': ctx.prec, 'quick': quick, 'zeros_file': 'data/zeros_1700.txt',
           'HSW22_bounds': {'S_bound_up_to_X0': '2.5167', 'X0': '30610046000', 'second_bound': '0.1095 log T + 0.2042 log log T + 3.0305 for T >= e'},
           'Trudgian_bound_as_in_chapter': '0.112 log T + 0.278 log log T + 2.510 (the published theorem has + 0.2/T)'}

    # (A) the ordinates
    zs = load_zeros(zpath, 300 if quick else 1700)
    n = len(zs)
    gaps = [zs[i + 1] - zs[i] for i in range(n - 1)]
    check('A: the %d ordinates are strictly increasing (certified)' % n, all(g > 0 for g in gaps))
    order = sorted(range(n - 1), key=lambda i: -float(gaps[i].mid()))
    i_max, i_later = order[0], max(range(1, n - 1), key=lambda i: float(gaps[i].mid()))
    check('A: the largest gap is the first one, gamma_2 minus gamma_1 (certified against every other gap)', i_max == 0 and all(gaps[0] > gaps[i] for i in range(1, n - 1)))
    check('A: the largest later gap is gamma_4 minus gamma_3 (certified against every gap with index at least 2)', i_later == 2 and all(gaps[2] > gaps[i] for i in range(3, n - 1)) and gaps[2] > gaps[1])
    out['ordinates'] = {'count': n, 'gamma_1': s(zs[0], 30), 'gamma_2': s(zs[1], 30), 'gamma_3': s(zs[2], 30), 'gamma_4': s(zs[3], 30),
                        'gamma_1500': s(zs[1499], 20) if n >= 1500 else None, 'gamma_last': s(zs[-1], 20),
                        'first_gap': s(gaps[0], 30), 'largest_later_gap': s(gaps[i_later], 30), 'largest_later_gap_between_indices': [i_later + 1, i_later + 2],
                        'largest_gaps': [{'n': i + 1, 'gap': s(gaps[i], 15)} for i in order[:8]], 'smallest_gap': s(min(gaps, key=lambda g: float(g.mid())), 10)}
    log('A: gamma_1 = %s, gamma_2 = %s, first gap %s, largest later gap %s (between ordinates %d and %d)' % (sf(zs[0], 12), sf(zs[1], 12), sf(gaps[0], 16), sf(gaps[i_later], 16), i_later + 1, i_later + 2))
    log('A: largest gaps (n, gap):', ', '.join('%d: %s' % (i + 1, sf(gaps[i], 8)) for i in order[:8]))
    v = json.load(open(vpath))['zeros']
    out['certified_1800_from_v1min_json'] = {k: v[k] for k in ('zeros_used', 'largest_ordinate', 'max_gap_gamma2_minus_gamma1', 'largest_later_gap', 'largest_later_gap_between')}
    vg1, vg2 = arb(v['max_gap_gamma2_minus_gamma1'].strip('[]')), arb(v['largest_later_gap'].strip('[]'))
    if not quick:
        check('A: v1min.json first gap overlaps the first gap from the 1700 ordinates (same value)', abs(vg1 - gaps[0]).upper() < arb('1e-25'))
        check('A: v1min.json largest later gap overlaps the largest later gap of the 1700 ordinates', abs(vg2 - gaps[i_later]).upper() < arb('1e-12'))
    T_cert = arb(v['largest_ordinate'].strip('[]'))
    log('A: certified summary: %s zeros, largest ordinate %s, max gap %s, largest later gap %s' % (v['zeros_used'], sf(T_cert, 14), v['max_gap_gamma2_minus_gamma1'][:26], v['largest_later_gap'][:26]))
    sp = [zs[i] for i in range(n)]
    Smax, Smin = None, None
    for i in range(n):
        base = theta(sp[i]) / pi + 1
        for N_here in (i, i + 1):              # N just below and just above the ordinate (simple zeros)
            Sv = arb(N_here) - base
            f = float(Sv.mid())
            Smax = f if Smax is None or f > Smax else Smax
            Smin = f if Smin is None or f < Smin else Smin
    out['sawtooth_S_on_the_ordinates'] = {'max': '%.6f' % Smax, 'min': '%.6f' % Smin, 'range': 'both sides of the first %d ordinates' % n}
    log('A: S(T) = N(T) minus theta(T)/pi minus 1 at both sides of the %d ordinates lies in [%.4f, %.4f]' % (n, Smin, Smax))

    # exact theta' and theta'' against independent computations
    Tt = arb(2524)
    eps = arb('1e-12')
    dq = (theta(Tt + eps) - theta(Tt - eps)) / (2 * eps)
    check('theta prime (digamma) equals the difference quotient of theta (log Gamma) at T = 2524 to 1e-9', abs(dq - thetap(Tt)).upper() < arb('1e-9'))
    dq2 = (thetap(Tt + eps) - thetap(Tt - eps)) / (2 * eps)
    check('theta second derivative (polygamma) equals the difference quotient of theta prime at T = 2524 to 1e-9', abs(dq2 - thetapp(Tt)).upper() < arb('1e-9'))
    try:
        mp.mp.dps = 30
        d1 = mp.siegeltheta(mp.mpf(2524), derivative=1)
        d2 = mp.siegeltheta(mp.mpf(2524), derivative=2)
        check('theta prime and theta second derivative agree with mpmath siegeltheta at T = 2524 to 1e-20',
              abs(arb(str(d1)) - thetap(Tt)).upper() < arb('1e-20') and abs(arb(str(d2)) - thetapp(Tt)).upper() < arb('1e-20'))
    except Exception as e:
        log('mpmath siegeltheta check skipped:', type(e).__name__, e)
    check('theta second derivative is positive at T = 284, 2524, X0 and 1e300 (it is positive for every T > 0)', all(thetapp(arb(t)) > 0 for t in (284, 2524, 30610046000)) and thetapp(arb(10) ** 300) > 0)

    check('theta second derivative >= 1/(2T) minus 1/T^2 at T = 284, 2524, X0, 1e30 (the bound used in (D); a ball of 320 bits cannot see the slack 1/T^2 beyond 1e45)',
          all(thetapp(arb(t)) >= 1 / (2 * arb(t)) - 1 / arb(t) ** 2 for t in (284, 2524, 30610046000, 10 ** 30)))

    # (B) Trudgian form of the lemma
    h83 = arb('8.3')
    tr = {}
    for variant in ('chapter', 'endpoints', 'published'):
        lhs, rhs = tr_sides(2524, h83, variant)
        lo, hi = bisect_true(lambda T: tr_diff(T, h83, variant) > 0, 100.0, 1.0e5)
        tr[variant] = {'lhs_at_2524': s(lhs, 12), 'rhs_at_2524': s(rhs, 12), 'margin_at_2524': s(lhs - rhs, 6), 'threshold_bracket': [repr(lo), repr(hi)]}
        log('B: variant %-9s at T = 2524, h = 8.3: left %s, right %s, margin %s; threshold in (%.6f, %.6f]' % (variant, sf(lhs, 8), sf(rhs, 8), sf(lhs - rhs, 3), lo, hi))
    claim('Trudgian left side 8.3 theta prime(2524)/pi', '291', '7.92027', tr_sides(2524, h83, 'chapter')[0], 5)
    claim('Trudgian right side 2 B(2524 + 8.3), no 0.2/T', '291', '7.92018', tr_sides(2524, h83, 'chapter')[1], 5)
    thr_ch = float(tr['chapter']['threshold_bracket'][1])
    log('B: chapter variant threshold %.4f, so the lemma needs T >= %d; published variant (with 0.2/T) threshold %.4f' % (thr_ch, int(thr_ch) + 1, float(tr['published']['threshold_bracket'][1])))
    CLAIMS.append({'claim': 'Trudgian threshold T >= 2524 (h = 8.3)', 'chapter_line': '291', 'chapter': '2524', 'computed': str(int(thr_ch) + 1), 'computed_ball': tr['chapter']['threshold_bracket'][1],
                   'agrees_to_printed_digits': int(thr_ch) + 1 == 2524})
    log('claim %-52s line %-4s chapter %-10s computed %-10s %s' % ('Trudgian threshold T >= 2524 (h = 8.3)', '291', '2524', str(int(thr_ch) + 1), 'agrees' if int(thr_ch) + 1 == 2524 else 'DIFFERS'))
    # the left minus right side increases with T: derivative h theta''/pi minus 2 B'(T + h) > 0 on a grid from the threshold to 1e300
    fp_ok, npts = True, 0
    for T in grid(2524, arb(10) ** 30, 1.5 if quick else 1.02):
        d = h83 * thetapp(T) / pi - 2 * b_tr_p(T + h83)
        fp_ok = fp_ok and d > 0
        npts += 1
    check('B: derivative of (left minus right) is positive on the geometric grid from 2524 to 1e30 (%d points); beyond, theta\'\' >= 1/(2T) minus 1/T^2 and 2 T B\'(T + h) <= 0.224 + 0.556/log T give it' % npts, fp_ok)
    gr_ok, npts = True, 0
    for T in grid(2524, arb(10) ** 30, 1.5 if quick else 1.01):
        gr_ok = gr_ok and tr_diff(T, h83, 'chapter') > 0
        npts += 1
    check('B: Trudgian inequality (chapter variant) holds at every point of the geometric grid from 2524 to 1e30 (%d points)' % npts, gr_ok)
    fail = {}
    for name, Tv in (('gamma_1500', zs[1499] if n >= 1500 else None), ('gamma_1700', zs[1699] if n >= 1700 else None), ('end_of_certified_1800', T_cert), ('T_2000', arb(2000)), ('T_2500', arb(2500))):
        if Tv is None:
            continue
        lhs, rhs = tr_sides(Tv, h83, 'chapter')
        lo_h, hi_h = bisect_true(lambda hh: tr_diff(Tv, hh, 'chapter') > 0, 1.0, 60.0)
        fail[name] = {'T': sf(Tv, 10), 'left_h_8.3': sf(lhs, 7), 'right_h_8.3': sf(rhs, 7), 'inequality_holds': bool(lhs > rhs), 'step_needed_chapter_variant': '%.4f' % hi_h,
                      'left_h_8.38': sf(arb('8.38') * thetap(Tv) / pi, 6), 'right_h_8.38': sf(2 * b_tr(Tv + arb('8.38')), 6)}
        log('B: at %s (T = %s): h = 8.3 gives %s against %s (%s); the step Trudgian needs here is %.4f' % (name, sf(Tv, 8), sf(lhs, 6), sf(rhs, 6), 'holds' if lhs > rhs else 'FAILS', hi_h))
    out['Trudgian'] = {'variants_at_T_2524_h_8.3': tr, 'failure_range': fail, 'threshold_chapter_variant': tr['chapter']['threshold_bracket'],
                       'note': 'the inequality fails for every T below the threshold (left minus right is increasing), so with the certified zeros to 2304.3 the range from 2304.3 to the threshold is uncovered'}
    check('B: the Trudgian inequality fails (certified false or undecided) at the end 2304.3 of the certified zeros', not (tr_diff(T_cert, h83, 'chapter') > 0))

    # (C) HSW22
    hsw = {}
    for hname in ('8.3', '6.8874', '5.36'):
        hh = arb(hname)
        lo, hi = bisect_true(lambda T: hsw_diff(T, hh) > 0, 7.0, 1.0e5)
        hsw[hname] = {'threshold_bracket': [repr(lo), repr(hi)]}
        log('C: HSW22 constant bound, h = %s: h theta prime(T)/pi > 2 x 2.5167 for T >= %.5f (bracket top %r)' % (hname, hi, hi))
    out['HSW22_thresholds'] = hsw
    claim('HSW22 threshold in T for h = 8.3 (round up)', 'review M10', '284', arb(float(hsw['8.3']['threshold_bracket'][1])).ceil(), 0)
    h_min_cert = 2 * S1 * pi / thetap(T_cert)
    out['HSW22_step_on_[T_cert, X0 minus h]'] = {'T_cert': sf(T_cert, 14), 'h_min': s(h_min_cert, 12)}
    log('C: smallest step the constant bound gives for T >= %s: h = 2 x 2.5167 x pi / theta prime(T) = %s' % (sf(T_cert, 10), sf(h_min_cert, 10)))
    CLAIMS.append({'claim': 'every gap below 5.36 on [2304.3, X0]', 'chapter_line': 'review M10', 'chapter': '5.36 (an upper bound)', 'computed': sf(h_min_cert, 8),
                   'computed_ball': s(h_min_cert, 15), 'agrees_to_printed_digits': bool(h_min_cert < arb('5.36'))})
    log('claim %-52s line %-4s chapter %-10s computed %-10s %s' % ('every gap below 5.36 on [2304.3, X0]', 'M10', '5.36', sf(h_min_cert, 8), 'agrees (computed below 5.36)' if h_min_cert < arb('5.36') else 'DIFFERS'))
    check('C: the HSW22 step on [2304.3, X0 minus h] is below 5.36', h_min_cert < arb('5.36'))

    # (D) above X0
    dd = {}
    for hname in ('8.3', '6.8874', '5.36'):
        hh = arb(hname)
        lhs = hh * thetap(X0) / pi
        rhs = 2 * b2(X0 + hh)
        Tm = X0 - hh
        dmin = d_h(Tm, hh)
        # explicit lower bound of D_h' for T >= X0 minus h, from theta'' >= 1/(2T) minus 1/T^2 (integral comparison for the sum over n)
        # and 2 T B_2'(T + h) <= 2 x 0.1095 + 2 x 0.2042/log(T): both worst at T = X0 minus h
        lb = hh / (2 * pi) - hh / (pi * Tm) - 2 * arb('0.1095') - 2 * arb('0.2042') / Tm.log()
        pp = d_h_p(Tm, hh)
        dd[hname] = {'left_at_X0': s(lhs, 8), 'right_at_X0': s(rhs, 8), 'D_at_X0_minus_h': s(dmin, 8), 'D_prime_at_X0_minus_h_times_T': s(pp * Tm, 8), 'coefficient_lower_bound_for_T_D_prime': s(lb, 8)}
        log('D: h = %s: at X0 left %s against right %s; D_h(X0 minus h) = %s; T D_h\'(X0 minus h) = %s; explicit lower bound for T D_h\' = %s' % (hname, sf(lhs, 7), sf(rhs, 7), sf(dmin, 6), sf(pp * Tm, 6), sf(lb, 6)))
        check('D: h = %s: D_h(X0 minus h) > 0 and the explicit lower bound of T D_h\' is positive' % hname, dmin > 0 and lb > 0 and pp > 0)
    out['above_X0'] = dd
    claim('left side at X0, h = 8.3 (review: 29.4)', 'review M10', '29.4', arb('8.3') * thetap(X0) / pi, 1)
    claim('right side at X0, h = 8.3 (review: 12.6)', 'review M10', '12.6', 2 * b2(X0 + arb('8.3')), 1)
    gtop = 30 if quick else 3000
    grid_info = {}
    for hname in ('8.3', '6.8874'):
        hh = arb(hname)
        okD, okP, npts, minD = True, True, 0, None
        t_start = X0 - hh
        lim_p = arb(10) ** 45          # D_h' is evaluated with the exact polygamma only up to here (it is slow beyond; the explicit bound covers the rest)
        segs = [(t_start, arb(10) ** 30, 1.5 if quick else 1.01), (arb(10) ** 30, arb(10) ** min(gtop, 300), 1.5 if quick else 1.1), (arb(10) ** 300, arb(10) ** gtop, 2.0)]
        for (a, b, r) in segs:
            if not (b > a):
                continue
            for T in grid(a, b, r):
                dv = d_h(T, hh)
                okD = okD and dv > 0
                if hname == '6.8874' and T < lim_p:
                    okP = okP and d_h_p(T, hh) > 0
                npts += 1
                if minD is None or float(dv.mid()) < minD:
                    minD = float(dv.mid())
        grid_info[hname] = {'points': npts, 'from': 'X0 minus h', 'to': '1e%d' % gtop, 'D_positive_at_every_point': bool(okD), 'D_prime_positive_at_every_point': bool(okP) if hname == '6.8874' else None, 'min_D_over_grid': '%.6f' % minD}
        check('D: h = %s: D_h > 0 at all %d grid points from X0 to 1e%d%s' % (hname, npts, gtop, ' and D_h\' > 0 at those below 1e45' if hname == '6.8874' else ''), okD and (okP if hname == '6.8874' else True))
    out['above_X0_grid'] = grid_info
    out['above_X0_argument'] = ('For T >= X0 minus h put D_h(T) = h theta\'(T)/pi minus 2 B_2(T + h), B_2 the second HSW22 bound (increasing, at least 2.5167 there), so that N(T + h) minus N(T) >= '
                                'h theta\'(T)/pi minus |S(T)| minus |S(T + h)| >= D_h(T) by convexity of theta. D_h\'(T) = h theta\'\'(T)/pi minus 2 B_2\'(T + h) with theta\'\'(T) >= 1/(2T) minus 1/T^2 '
                                'and 2 T B_2\'(T + h) <= 0.219 + 0.4084/log T, so T D_h\'(T) >= h/(2 pi) minus h/(pi T) minus 0.219 minus 0.4084/log T, positive for h >= 1.5; hence D_h(T) >= D_h(X0 minus h) > 0.')

    # (E) conclusion
    hstar = arb('6.8874')
    lo6, hi6 = float(hsw['6.8874']['threshold_bracket'][0]), float(hsw['6.8874']['threshold_bracket'][1])
    nz = sum(1 for z in zs if z < arb(hi6))
    check('E: the first gap gamma_2 minus gamma_1 is below 6.8874 (certified)', gaps[0] < hstar)
    check('E: the HSW22 lemma with h = 6.8874 starts at T = %.2f, below the end 2304.3 of the certified zeros, so the two ranges overlap' % hi6, arb(hi6) < T_cert)
    check('E: the largest gap among the ordinates above gamma_2 (5.414...) is below 5.42', gaps[i_later] < arb('5.42'))
    out['conclusion'] = {'every_gap_at_most': s(gaps[0], 20), 'every_gap_below': '6.8874', 'gaps_above_gamma_2_at_most': s(max(gaps[i_later], h_min_cert), 15),
                         'HSW22_lemma_starts_at_T': '%.4f' % hi6, 'ordinates_below_that_height_in_the_file': nz,
                         'proof': 'gaps between ordinates below the end of the certified zeros: computed; above T = %.2f: N(T + h) minus N(T) > 0 for h = 6.8874 by (C) and (D)' % hi6}
    log('E: every gap is at most %s (< 6.8874); above gamma_2 at most max(%s, %s) = %s; the HSW22 argument with h = 6.8874 starts at T = %.2f (%d ordinates of the file lie below it)' % (sf(gaps[0], 12), sf(gaps[i_later], 10), sf(h_min_cert, 8), sf(max(gaps[i_later], h_min_cert), 10), hi6, nz))

    # (F) the algebra and the margins of the proof of Theorem yo:thm:coercive (ii)
    rng = random.Random(20261017)
    bad = 0
    nsamp = 400 if quick else 4000
    for _ in range(nsamp):
        g = Fraction(rng.randint(141, 40000), 10)
        c = Fraction(rng.randint(1, 500), 1000)
        t = Fraction(rng.randint(0, 80000), 10)
        P = lambda tt: c / (c * c + (tt - g) ** 2) + c / (c * c + (tt + g) ** 2)
        geq = P(t) >= P(Fraction(0))
        if geq != (t * t <= 3 * g * g - c * c):
            bad += 1
    tb = Fraction(3) * Fraction(212, 10) ** 2 - Fraction(1, 4)       # tau^2 = 3 gamma^2 minus c^2 exactly: equality
    gg, cc = Fraction(212, 10), Fraction(1, 2)
    Pf = lambda tt: cc / (cc * cc + (tt - gg) ** 2) + cc / (cc * cc + (tt + gg) ** 2)
    eq_ok = False
    # the boundary tau is irrational in general, so test it on P as a function of v = tau^2: P = 2c(A+v)/((A+v)^2 minus 4 gamma^2 v)
    A = cc * cc + gg * gg
    v = 3 * gg * gg - cc * cc
    eq_ok = (2 * cc * (A + v) / ((A + v) ** 2 - 4 * gg * gg * v)) == 2 * cc / A
    check('F: P_gamma(tau) >= P_gamma(0) iff tau^2 <= 3 gamma^2 minus c^2, exact rational arithmetic, %d random samples, and equality at the boundary' % nsamp, bad == 0 and eq_ok)
    mp.mp.dps = 40
    z0 = mp.mpf(1) / 2
    zeta0, zeta1, zeta2 = mp.zeta(z0), mp.zeta(z0, 1, 1), mp.zeta(z0, 1, 2)
    two_k1 = -1 / z0 ** 2 - 1 / (z0 - 1) ** 2 + mp.psi(1, z0 / 2) / 4 + (zeta2 / zeta0 - (zeta1 / zeta0) ** 2)
    kappa1 = arb(mp.nstr(two_k1 / 2, 38))
    kappa1_ball = arb(mp.nstr(two_k1 / 2, 30) + ' +/- 1e-30')
    sum_inv = arb(0)
    for z in zs:
        sum_inv += 1 / (z * z)
    Tm = zs[-1] + (zs[-1] - zs[-2]) / 2
    tail = ((Tm / (2 * pi)).log() + 1) / (2 * pi * Tm)
    k1_zeros = sum_inv + tail
    out['kappa_1'] = {'closed_form_mpmath_40_digits': mp.nstr(two_k1 / 2, 25), 'from_zeros_plus_RvM_tail': s(k1_zeros, 10), 'partial_sum_of_zeros': s(sum_inv, 12), 'tail_estimate': s(tail, 8)}
    log('F: kappa_1 = (1/2)(log xi)\'\'(1/2) = %s from the closed form; sum over the %d zeros %s plus the tail %s = %s' % (mp.nstr(two_k1 / 2, 15), n, sf(sum_inv, 10), sf(tail, 6), sf(k1_zeros, 10)))
    check('F: kappa_1 from the zeros and the tail agrees with the closed form to %s' % ('3e-5' if quick else '3e-6'), abs(k1_zeros - kappa1_ball).upper() < arb('3e-5' if quick else '3e-6'))
    claim('kappa_1 = sum 1/gamma^2 (line 56)', '56', '0.0231049931', kappa1_ball, 10)
    claim('2 kappa_1 = L(0) (line 298)', '298', '0.046210', 2 * kappa1_ball, 6)
    c = arb(1) / 2
    mar = {}
    for name, hg in (('chapter_8.3', arb('8.3')), ('corrected_6.8874', arb('6.8874'))):
        d = hg / 2
        lor = 1 / (c * c + d * d)
        g1 = 2 / (zs[0] * zs[0])
        two_k = 2 * kappa1_ball
        margin = lor - g1 - two_k
        mar[name] = {'distance': s(d, 8), 'one_over_c2_plus_d2': s(lor, 10), 'two_over_gamma1_squared': s(g1, 10), 'two_kappa1': s(two_k, 10), 'margin_chapter_form': s(margin, 10),
                     'margin_without_the_2_over_gamma1_squared_term': s(lor - two_k, 10)}
    out['coercive_margins_c_half'] = mar
    claim('1/(c^2 + 4.15^2) at c = 1/2', '278', '0.057233', 1 / (c * c + arb('4.15') ** 2), 6)
    claim('1/(c^2 + 3.4437^2) at c = 1/2', 'review M10', '0.082583', 1 / (c * c + arb('3.4437') ** 2), 6)
    claim('2/gamma_1^2', '278', '0.010010', 2 / (zs[0] * zs[0]), 6)
    claim('margin 0.057233 minus 0.010010 minus 0.046210 (unrounded inputs)', '278', '0.001013', 1 / (c * c + arb('4.15') ** 2) - 2 / (zs[0] * zs[0]) - 2 * kappa1_ball, 6)
    claim('corrected margin (unrounded inputs)', 'review M10', '0.026362', 1 / (c * c + arb('3.4437') ** 2) - 2 / (zs[0] * zs[0]) - 2 * kappa1_ball, 6)
    m_ch = 1 / (c * c + arb('4.15') ** 2) - 2 / (zs[0] * zs[0]) - 2 * kappa1_ball
    m_co = 1 / (c * c + arb('3.4437') ** 2) - 2 / (zs[0] * zs[0]) - 2 * kappa1_ball
    log('F: at c = 1/2 the chapter margin is %s and the corrected margin %s; without the term 2/gamma_1^2 (which repeats P_gamma_n(0), already inside the sum of all P_gamma(0) = xi\'/xi(alpha) <= 2 c kappa_1) they are %s and %s' % (
        sf(m_ch, 6), sf(m_co, 6), sf(1 / (c * c + arb('4.15') ** 2) - 2 * kappa1_ball, 6), sf(1 / (c * c + arb('3.4437') ** 2) - 2 * kappa1_ball, 6)))
    check('F: 4.15 = 8.3/2 and 3.4437 = 6.8874/2 exactly; the largest distance from tau >= gamma_1 to the nearest ordinate is half of the largest gap, %s' % sf(gaps[0] / 2, 12),
          abs(arb('8.3') / 2 - arb('4.15')).upper() < arb('1e-50') and abs(arb('6.8874') / 2 - arb('3.4437')).upper() < arb('1e-50') and gaps[0] / 2 < arb('3.4437'))

    claim('gamma_2 minus gamma_1 (line 291: 6.887; review: 6.8873)', '291', '6.8873', gaps[0], 4)
    claim('largest gap above gamma_2, gamma_4 minus gamma_3 (review: 5.414)', 'review M10', '5.414', gaps[i_later], 3)
    if n >= 1500:
        claim('gamma_1500 (review M9)', 'review M9', '1980.91', zs[1499], 2)
    log('claims that differ from the printed digits:', [c_['claim'] for c_ in CLAIMS if not c_['agrees_to_printed_digits']] or 'none')

    out['claims'] = CLAIMS
    out['checks'] = CHECKS
    out['all_checks_pass'] = all(CHECKS.values())
    out['elapsed_s'] = round(time.time() - T0, 1)
    log('all checks pass:', out['all_checks_pass'], '; elapsed %.1f s' % (time.time() - T0))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(out, fh, indent=1)
    with open(os.path.splitext(out_path)[0] + '.log', 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    log('wrote', out_path, 'and', os.path.splitext(out_path)[0] + '.log')


if __name__ == '__main__':
    main()
