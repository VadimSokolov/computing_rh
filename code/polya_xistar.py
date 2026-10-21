# polya_xistar.py: Polya's function Xi*, Chapter 16 (Polya and the Zeros of Entire Functions), Section pol:sec:xistar (Proposition
# pol:prop:bessel, Proposition pol:prop:gram and the Numerics paragraph that follows them) and the right series of Figure
# pol:fig:one. It checks: the quadrature of the kernel against the Bessel form, Xi*(t) = 2 int_0^oo Phi*(u) cos(tu) du =
# 8 pi^2 Re K_{9/4 + it/2}(2 pi), at t = 0 and t = 10 (and at larger t); the number of zeros of Xi* in (0,100] (29, as for Xi, whose
# ordinates are in data/zeros_1700.txt) and the first eight; theta(t_k*)/pi for k = 1..5, 20, 25; the 138 zeros up to height 300 and
# t*_138; the deviation (theta(t_k*)/pi - (k - 3/2)) t_k* against -(2 pi^2 - 3)/pi for all zeros, with the remainder of the expansion
# of Proposition pol:prop:gram; the phase correction t (arg K_nu(2 pi) - theta(t) - pi) against 2 pi^2 - 3 at t = 100, 200, 400, 800
# and further (to 6400); the counting formula N*(T) = floor(theta(T)/pi + 3/2); the ratio Phi*/Phi at u = 0.5, 1, 1.5. Added for the
# book in October 2026, written from the mathematics (the chapter author's own scripts were neither read nor run).
# Method. Everything is mpmath. The kernel integral is done by tanh sinh quadrature on subintervals of [0, 3] (at most 1.5 cycles of
# cos(tu) each), with the working precision raised to 0.69 t + 55 digits because Xi*(t) is of size exp(-pi t/4) while the integrand is
# of size one. The Bessel form uses mpmath's besselk of complex order at the same precision (its hypergeometric series loses 2 pi t/4
# nats to cancellation, which the precision covers, and the result is compared at two precisions). The zeros are found by scanning
# the normalised function Re K / |K| (so a sign change is a zero of Xi*) on a grid of step 0.05 up to height 320, refined by Brent's
# method in double precision and polished by the secant method at 0.69 t + 80 digits (the zeros are then good to about 30 digits).
# theta is mpmath's siegeltheta. The phase of K is the principal argument of K reduced modulo 2 pi after subtracting theta + pi (the
# correction is small, so no branch can be lost). The analytic coefficient of the next term, (phase - theta - pi) = C/t + c3/t^3 + ...
# with C = 2 pi^2 - 3 and c3 = 21/4 - (25/2) pi^2 + 4 pi^4 (from arg s + arg(s+2) = pi - 3/t + (21/4)/t^3 and from the series
# 1 + pi^2/(1-nu) + pi^4/(2(1-nu)(2-nu)) + ..., whose logarithm is pi^2/(1-nu) - pi^4/(2 (-i t/2)^3) + O(t^-4)), is derived by hand and
# compared with (t (phase - theta - pi) - C) t^2 at the sampled heights. Error control. (a) The kernel and Bessel values are written with
# their difference. (b) K_nu(2 pi) from the Bessel function is compared with the integral int_0^oo exp(-2 pi cosh v) cosh(nu v) dv
# at t = 100 and 200 (real and imaginary parts). (c) Besselk is evaluated at two precisions at t = 800 and 6400, and the phase is formed at 80
# digits once K is known (theta is then also taken at 80 digits). (d) The
# zero count is compared with the formula floor(theta(T)/pi + 3/2) away from the zeros. Nothing here is a proof: the zero count
# comes from sign changes on a grid, the real rootedness of Xi* is Polya's theorem.
# Output: data/polya_xistar.json (values as strings) and data/polya_xistar.log, both next to --out (the log has the same name with
# the extension .log).
# Usage: python3 polya_xistar.py [--quick] [--out FILE] [--workers W] [--zeros FILE]; --quick scans to height 60 (about 4 s on one
# core) and writes to /tmp unless --out is given; the full case is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8
# -t 90 -g polya_xistar.json -g polya_xistar.log code/polya_xistar.py data/zeros_1700.txt -- --out polya_xistar.json
import argparse
import atexit
import json
import math
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
T0 = time.time()
LOGF = None


def log(*a):
    line = '[%7.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(line, flush=True)
    if LOGF is not None:
        LOGF.write(line + '\n')
        LOGF.flush()


def fs(x, n=12):
    return format(float(x), '.%dg' % n)


def ms(x, n=15):
    return mp.nstr(x, n)


def dps_for(t):
    return int(0.69 * float(t) + 55)


# ---------------------------------------------------------------- the functions
def phistar(u):
    return 8 * mp.pi ** 2 * mp.cosh(9 * u / 2) * mp.exp(-2 * mp.pi * mp.cosh(2 * u))


def phi(u):
    u = mp.mpf(u)
    s = mp.mpf(0)
    for n in range(1, 60):
        term = (2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * u / 2)) * mp.exp(-mp.pi * n ** 2 * mp.exp(2 * u))
        s += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 10) * abs(s):
            break
    return 2 * s


def K_bessel(t):
    return mp.besselk(mp.mpc(mp.mpf(9) / 4, mp.mpf(t) / 2), 2 * mp.pi)


def xistar_kernel(t):
    t = mp.mpf(t)
    U = mp.mpf(3)
    m = max(2, int(float(t) * 3 / (2 * math.pi) / 1.5) + 2)
    pts = [U * i / m for i in range(m + 1)]
    return 2 * mp.quad(lambda u: phistar(u) * mp.cos(t * u), pts)


def K_kernel(t):
    t = mp.mpf(t)
    nu = mp.mpc(mp.mpf(9) / 4, t / 2)
    V = mp.mpf('4.5')
    m = max(2, int(float(t) * 4.5 / (4 * math.pi) / 1.5) + 2)
    pts = [V * i / m for i in range(m + 1)]
    return mp.quad(lambda v: mp.exp(-2 * mp.pi * mp.cosh(v)) * mp.cosh(nu * v), pts)


# ---------------------------------------------------------------- workers
def w_scan(ts):
    mp.mp.dps = dps_for(max(ts))
    out = []
    for t in ts:
        K = K_bessel(t)
        out.append(float(mp.re(K) / abs(K)))
    return out


def w_zero(task):
    k, t0 = task
    mp.mp.dps = dps_for(t0) + 25
    f = lambda x: mp.re(K_bessel(x))
    x0 = mp.mpf(t0)
    f0 = f(x0)
    x1 = x0 + mp.mpf(10) ** -10
    f1 = f(x1)
    for _ in range(12):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
        f1 = f(x1)
        if abs(x1 - x0) < mp.mpf(10) ** -32:
            break
    t = x1
    th = mp.siegeltheta(t) / mp.pi
    dev = th - (k - mp.mpf(3) / 2)
    C = 2 * mp.pi ** 2 - 3
    D = dev * t
    r = (dev + C / (mp.pi * t)) * t ** 2
    Tsw = mp.findroot(lambda T: mp.siegeltheta(T) - (k - mp.mpf(3) / 2) * mp.pi, t + mp.mpf('0.01'))
    return dict(k=k, t=ms(t, 30), theta_over_pi=ms(th, 20), dev=ms(dev, 15), D=ms(D, 15), r=ms(r, 12), window=ms(Tsw - t, 8),
                residual_at_root=ms(abs(f(t)) / abs(K_bessel(t)), 3))


def w_phase(t):
    try:
        return w_phase_inner(t)
    except Exception as e:
        return dict(t=str(t), error=repr(e))


def w_phase_inner(t):
    mp.mp.dps = dps_for(t) + 25
    t = mp.mpf(t)
    K = K_bessel(t)
    K2 = None
    if t in (800, 6400):
        mp.mp.dps = dps_for(t) + 95
        K2 = K_bessel(t)
    mp.mp.dps = 80
    d = mp.arg(K) - mp.siegeltheta(t) - mp.pi
    d = (d + mp.pi) % (2 * mp.pi) - mp.pi
    P = t * d
    C = 2 * mp.pi ** 2 - 3
    out = dict(t=ms(t, 8), P=ms(P, 15), P_minus_C=ms(P - C, 10), E=ms((P - C) * t ** 2, 12))
    if K2 is not None:
        out['besselk_relative_difference_two_precisions'] = ms(abs(K - K2) / abs(K2), 3)
    return out


def w_check(task):
    kind, t = task
    mp.mp.dps = dps_for(t) + 15
    if kind == 'xistar':
        kern = xistar_kernel(t)
        bes = 8 * mp.pi ** 2 * mp.re(K_bessel(t))
        return dict(kind=kind, t=ms(t, 6), kernel=ms(kern, 25), bessel=ms(bes, 25), difference=ms(kern - bes, 3),
                    relative_difference=ms(abs(kern - bes) / abs(bes), 3))
    kk = K_kernel(t)
    kb = K_bessel(t)
    return dict(kind=kind, t=ms(t, 6), K_kernel=ms(kk, 20), K_bessel=ms(kb, 20), relative_difference=ms(abs(kk - kb) / abs(kb), 3))


def w_task(task):
    name, arg = task
    return globals()[name](arg)


_POOL = None


def pmap(name, args, workers):
    global _POOL
    tasks = [(name, a) for a in args]
    if workers <= 1 or len(tasks) == 1:
        return [w_task(t) for t in tasks]
    if _POOL is None:
        try:
            _POOL = get_context('fork').Pool(workers)
            atexit.register(_POOL.terminate)
        except Exception as e:
            log('the worker pool failed (%r): running in one process' % (e,))
            _POOL = False
    if not _POOL:
        return [w_task(t) for t in tasks]
    return _POOL.map(w_task, tasks, chunksize=1)


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--zeros', default=None)
    a = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    out_path = Path(a.out) if a.out else (Path('/tmp/polya_xistar_quick.json') if a.quick else root / 'data' / 'polya_xistar.json')
    LOGF = open(out_path.with_suffix('.log'), 'w')
    workers = 1 if a.quick else a.workers
    zf = None
    for cand in ([Path(a.zeros)] if a.zeros else []) + [Path('zeros_1700.txt'), root / 'data' / 'zeros_1700.txt']:
        if cand.exists():
            zf = cand
            break
    ordinates = np.array([float(x) for x in zf.read_text().split()])
    if a.quick:
        P = dict(tmax=60.0, step=0.1, phase=[50.0, 100.0], checks=[('xistar', 0), ('xistar', 10), ('xistar', 30)], count_T=[30.0, 50.0])
    else:
        P = dict(tmax=320.0, step=0.05, phase=[100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0, 6400.0],
                 checks=[('xistar', 0), ('xistar', 10), ('xistar', 50), ('xistar', 100), ('xistar', 200), ('K', 100), ('K', 200)],
                 count_T=[50.0, 100.0, 200.0, 300.0])
    R = {'meta': {'quick': a.quick, 'workers': workers, 'mpmath': mp.__version__, 'parameters': {k: str(v) for k, v in P.items()}}}
    C = 2 * mp.pi ** 2 - 3
    mp.mp.dps = 30
    R['constants'] = {'C_2pi2_minus_3': ms(C, 15), 'minus_C_over_pi': ms(-C / mp.pi, 15),
                      'c3_analytic': ms(mp.mpf(21) / 4 - mp.mpf(25) / 2 * mp.pi ** 2 + 4 * mp.pi ** 4, 12)}
    log('workers', workers, 'C = 2 pi^2 - 3 =', R['constants']['C_2pi2_minus_3'], '-C/pi =', R['constants']['minus_C_over_pi'], 'c3 (analytic) =',
        R['constants']['c3_analytic'], 'N(100) =', int(np.sum(ordinates <= 100)), 'N(300) =', int(np.sum(ordinates <= 300)))

    # ---- kernel against Bessel
    R['kernel_vs_bessel'] = pmap('w_check', P['checks'], workers)
    for r in R['kernel_vs_bessel']:
        log('check', r)
    log('chapter: Xi*(0) = 0.10502899694, Xi*(10) = -0.00345296144711')

    # ---- scan and zeros
    n = int(round(P['tmax'] / P['step']))
    ts = [(j + 1) * P['step'] for j in range(n)]
    chunks = [ts[i:i + 40] for i in range(0, len(ts), 40)]
    g = np.concatenate([np.array(x) for x in pmap('w_scan', chunks, workers)])
    tarr = np.array(ts)
    s = np.sign(g)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    log('scan: %d points, %d sign changes up to %s; first value Re K/|K| %s (Xi*(0) > 0 so the first sign change is a zero)' % (len(ts), len(idx), fs(tarr[-1]), fs(g[0], 4)))
    roots_d = []
    for k in idx:
        mp.mp.dps = dps_for(tarr[k + 1])
        def gg(t):
            mp.mp.dps = dps_for(t)
            K = K_bessel(t)
            return float(mp.re(K) / abs(K))
        roots_d.append(brentq(gg, tarr[k], tarr[k + 1], xtol=1e-14, rtol=1e-15, maxiter=100))
    zs = pmap('w_zero', [(k + 1, r) for k, r in enumerate(roots_d)], workers)
    t_all = [mp.mpf(z['t']) for z in zs]
    nle = lambda T: int(sum(1 for t in t_all if t <= T))
    R['zeros'] = {'count_le_100': nle(100), 'count_le_300': nle(300), 'count_total_scanned': len(zs), 'scan_height': fs(tarr[-1]),
                  'N_Xi_le_100': int(np.sum(ordinates <= 100)), 'first8': [z['t'][:12] for z in zs[:8]],
                  'first8_Xi_ordinates': [fs(x, 10) for x in ordinates[:8]], 't_138': zs[137]['t'] if len(zs) > 137 else None,
                  't_139': zs[138]['t'] if len(zs) > 138 else None,
                  'smallest_gap': fs(min(float(t_all[i + 1] - t_all[i]) for i in range(len(t_all) - 1)), 6), 'all': [z['t'] for z in zs]}
    log('zeros: count <= 100:', nle(100), '(Xi:', int(np.sum(ordinates <= 100)), ')  count <= 300:', nle(300), ' t*_138 =', zs[min(137, len(zs) - 1)]['t'][:14],
        ' t*_139 =', zs[min(138, len(zs) - 1)]['t'][:14] if len(zs) > 138 else 'beyond the scan')
    log('first eight zeros of Xi*:', [z['t'][:10] for z in zs[:8]])
    log('first eight ordinates of Xi:', [fs(x, 9) for x in ordinates[:8]])
    log('chapter: 8.9928 19.0654 24.4748 28.9750 33.0069 36.7367 40.2491 43.5947; 14.1347 21.0220 25.0109 30.4249 32.9351 37.5862 40.9187 43.3271; t*_138 = 299.627')
    sel = [1, 2, 3, 4, 5, 20, 25]
    R['theta_over_pi'] = {str(k): zs[k - 1]['theta_over_pi'] for k in sel if k <= len(zs)}
    log('theta(t_k*)/pi at k = 1..5, 20, 25:', {k: v[:9] for k, v in R['theta_over_pi'].items()}, '(chapter -1.042 0.209 1.277 2.313 3.336, 18.431, 23.440)')

    # ---- the expansion of Proposition gram for all zeros
    D = [float(mp.mpf(z['D'])) for z in zs]
    r = [float(mp.mpf(z['r'])) for z in zs]
    tt = [float(z) for z in t_all]
    last = min(137, len(zs) - 1)
    R['gram'] = {'predicted_deviation_minus_C_over_pi': ms(-C / mp.pi, 12), 'D_first': fs(D[0], 10), 'D_at_last_zero_below_300': zs[last]['D'][:12],
                 'last_zero': zs[last]['t'][:14], 'max_abs_remainder_scaled_by_t2': fs(max(abs(x) for x in r), 8),
                 'remainder_scaled_by_t2_last': fs(r[last], 8), 'remainder_scaled_by_t3_last': fs(r[last] * tt[last], 8),
                 'max_abs_D_minus_predicted_over_zeros_ge_100': fs(max((abs(D[i] - float(-C / mp.pi)) for i in range(len(zs)) if tt[i] >= 100), default=0.0), 8),
                 'analytic_remainder_coefficient_minus_c3_over_pi': ms(-(mp.mpf(21) / 4 - mp.mpf(25) / 2 * mp.pi ** 2 + 4 * mp.pi ** 4) / mp.pi, 10)}
    log('Proposition gram: D_k = (theta/pi - (k - 3/2)) t_k at the last zero below 300:', zs[last]['D'][:12], 'against -C/pi =', R['gram']['predicted_deviation_minus_C_over_pi'],
        '(chapter -5.329 against -5.328)')
    log('   remainder r_k = (theta/pi - (k - 3/2) + C/(pi t_k)) t_k^2: max |r_k| =', R['gram']['max_abs_remainder_scaled_by_t2'], 'last', R['gram']['remainder_scaled_by_t2_last'],
        'r_k t_k at the last zero', R['gram']['remainder_scaled_by_t3_last'], '; analytic limit of r_k t_k is -c3/pi =', R['gram']['analytic_remainder_coefficient_minus_c3_over_pi'])
    log('   largest |D_k + C/pi| over the zeros with t_k >= 100:', R['gram']['max_abs_D_minus_predicted_over_zeros_ge_100'])
    # counting formula N*(T) = floor(theta(T)/pi + 3/2): fails in windows after each zero
    mp.mp.dps = 30
    cnt = {}
    for T in P['count_T']:
        mT = mp.mpf(T) + mp.mpf('0.0123')
        cnt[fs(T)] = {'T_tested': ms(mT, 8), 'N_star': nle(mT), 'floor_formula': int(mp.floor(mp.siegeltheta(mT) / mp.pi + mp.mpf(3) / 2))}
    wins = [float(mp.mpf(z['window'])) for z in zs]
    sel_w = [(tt[i], wins[i]) for i in range(len(zs)) if 100 <= tt[i] <= 300]
    R['counting_formula'] = {'tests': cnt, 'window_after_zero_last_below_300': zs[last]['window'],
                             'window_after_zero_1': zs[0]['window'],
                             'fraction_of_heights_100_to_300_where_floor_formula_is_off_by_one': fs(sum(w for _, w in sel_w) / 200.0, 5) if sel_w else None,
                             'predicted_window_C_over_t_thetaprime_last': ms(C / (mp.mpf(zs[last]['t']) * (mp.log(mp.mpf(zs[last]['t']) / (2 * mp.pi)) / 2)), 6)}
    mism, inwin = [], 0
    for k in range(1, len(zs)):
        Tm = (t_all[k - 1] + t_all[k]) / 2
        if int(mp.floor(mp.siegeltheta(Tm) / mp.pi + mp.mpf(3) / 2)) != k:
            mism.append(k)
    for k in range(1, len(zs) + 1):
        Tw = t_all[k - 1] + mp.mpf(zs[k - 1]['window']) / 2
        inwin += int(int(mp.floor(mp.siegeltheta(Tw) / mp.pi + mp.mpf(3) / 2)) == k - 1)
    R['counting_formula'].update({'midpoints_between_zeros_k_and_k_plus_1_where_floor_formula_is_not_k': mism,
                                  'zeros_k_whose_window_midpoint_gives_k_minus_1': inwin, 'zeros_checked': len(zs)})
    log('counting formula: N*(T) against floor(theta(T)/pi + 3/2) at', cnt)
    log('   at the midpoints between consecutive zeros the formula fails for k =', mism, '; just after zero k (the window midpoint) it gives k - 1 for', inwin, 'of', len(zs), 'zeros')
    log('   the floor formula undercounts by one in a window after each zero: width', zs[last]['window'], 'after the last zero (predicted C/(t theta\')', R['counting_formula']['predicted_window_C_over_t_thetaprime_last'], ')',
        'and', zs[0]['window'], 'after t_1; fraction of heights in [100,300] affected', R['counting_formula']['fraction_of_heights_100_to_300_where_floor_formula_is_off_by_one'])

    # ---- phase correction at large t
    ph = pmap('w_phase', sorted(P['phase'], reverse=True), workers)
    ph.sort(key=lambda d: float(d['t']))
    R['phase_correction'] = ph
    for d in ph:
        log('phase correction t =', d['t'], ' t (phase - theta - pi) =', d.get('P'), ' minus C =', d.get('P_minus_C'), ' (that minus C) t^2 =', d.get('E'),
            d.get('besselk_relative_difference_two_precisions', ''), d.get('error', ''))
    log('chapter: 16.767 16.746 16.741 16.7396 at t = 100 200 400 800, against 2 pi^2 - 3 = 16.7392; analytic c3 =', R['constants']['c3_analytic'])

    # ---- Phi*/Phi
    mp.mp.dps = 30
    R['ratio_phistar_over_phi'] = {fs(u, 3): ms(phistar(mp.mpf(u)) / phi(u), 12) for u in (0.5, 1.0, 1.5)}
    R['ratio_phistar_over_phi'].update({fs(u, 3): ms(phistar(mp.mpf(u)) / phi(u), 12) for u in (0.0, 2.0, 3.0)})
    log('Phi*/Phi:', R['ratio_phistar_over_phi'], '(chapter 0.386 0.699 0.876 at u = 0.5 1 1.5)')

    # ---- Figure 1 right series
    R['figure1_right'] = {'k': list(range(1, len(zs) + 1)), 't': [z['t'][:14] for z in zs], 'theta_over_pi_minus_k_minus_3_over_2': [z['dev'][:14] for z in zs],
                          'asymptote_constant_minus_C_over_pi': ms(-C / mp.pi, 12), 'zeros_up_to_300': nle(300)}
    log('figure 1 right: %d values' % len(zs))
    R['meta']['seconds'] = fs(time.time() - T0, 4)
    out_path.write_text(json.dumps(R, indent=1))
    log('wrote', out_path)


if __name__ == '__main__':
    main()
