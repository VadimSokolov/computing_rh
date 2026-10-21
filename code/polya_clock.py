# polya_clock.py: the clock of the base point law and the even and odd parts of Polya's family, Chapter 16 (Polya and the Zeros of
# Entire Functions). It checks the numerical statements of Section pol:sec:cosh (Table pol:tab:clock: zeros of Re xi(alpha + it) on
# (0,100], the clock U_alpha(100), its rounding, the largest distance of the clock from Z + 1/2 at the computed zeros; the first
# three atoms at alpha = 1 and 0.52; the strict increase of the clock), of Theorem pol:thm:exact (0 < R_alpha < 1, the defects
# log((1 + R_alpha(s))/2) at s = 1 and 100, their linearity in a = alpha - 1/2 against -a (xi'/xi)(1/2 + sqrt s)), of Section
# pol:sec:progress (the sign changes of Re and Im xi(alpha + it); the argument principle counts of the even part A_alpha and the odd
# part B_alpha in the rectangles of Proposition pol:prop:local and Theorem pol:thm:accounting, for |Re w| <= 1 and for |Re w| <= 2;
# the zeros of xi'(1/2 + w)) and the left series of Figure pol:fig:one (the clock U_{3/4} on [0,60] and its half integer crossings).
# Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were neither read nor run).
# Method. xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s) is evaluated with python-flint (ball arithmetic at 100 bits, the midpoints are
# kept) or, when python-flint is missing or fails its self test against mpmath, with mpmath at 20 digits (the local quick run). The
# clock U_alpha(tau) = (1/pi) arg xi(alpha + i tau) is unwrapped on a grid: the principal argument of the ratio of consecutive unit
# vectors is summed (the step is safe when every increment stays below pi, the largest is recorded). Table pol:tab:clock is built
# in the way the text describes (grid of step 0.01, zeros of Re xi by linear interpolation of Re xi between the grid points where
# it changes sign, the clock at such a zero by linear interpolation of the grid clock) and again with the finer method (grid of
# step 0.0025, each zero refined by Brent's method on Re xi / |xi| to 1e-14, the clock at the refined zero from the unwrapped value
# at the grid point below plus the principal argument of the short ratio); a third variant evaluates the clock exactly at the
# linear interpolation zeros, which separates the error of the zero from the error of the interpolated clock. As an independent
# check of U_alpha(T), the clock is also integrated from U' = (1/pi) Re xi'/xi, with xi'/xi = 1/s + 1/(s-1) - (1/2) log pi +
# (1/2) psi(s/2) + zeta'/zeta(s), by Gauss Legendre on panels no longer than 5 (alpha - 1/2) (the Lorentzian width). The argument
# principle counts the zeros of A_alpha(w) = (xi(alpha + w) + xi(alpha - w))/2 and B_alpha(w) = (xi(alpha + w) - xi(alpha - w))/(2i)
# inside the rectangle |Re w| <= X, y0 <= Im w <= T as (1/2 pi) times the change of the argument along its boundary (four sides, no
# symmetry used), the sides sampled with steps 0.02 (vertical) and 0.01 (horizontal) and every pair of consecutive samples whose
# argument jumps by more than 0.5 rad, or whose modulus changes by more than a factor e, bisected until none does. The same count is
# done for g(w) = xi'(1/2 + w) in |Re w| <= 0.6, 1 <= Im w <= 100, next to the sign changes of the real function
# t -> -i xi'(1/2 + it) on (1,100]. The defects of Theorem pol:thm:exact use mpmath at 30 digits on the real axis (xi(x) = xi(1 - x)).
# Error control. (a) The self test compares flint with mpmath for xi and xi'/xi at four points. (b) The largest argument increment of
# every unwrapped grid is recorded (it stays far below pi). (c) U(T) from unwrapping is compared with U(T) from the integral of
# U'. (d) Every count of an argument principle contour is reported with its distance from an integer, the number of samples, the
# smallest modulus on the boundary (the boundary keeps away from zeros) and the number of unresolved jumps. (e) The refined zeros
# are checked to sit at half integer levels of the clock to about 1e-12. Nothing here is a proof: the counts are floating point.
# Output: data/polya_clock.json (values as strings) and data/polya_clock.log, both next to --out (the log has the same name with the
# extension .log).
# Usage: python3 polya_clock.py [--quick] [--out FILE] [--workers W] [--zeros FILE]; --quick uses two base points, height 15 and
# coarse grids (about 3 s on one core with mpmath) and writes to /tmp unless --out is given; the full case is meant for the Hopper
# cluster: bash misc/tools/hopper_run.sh -c 8 -t 60 -g polya_clock.json -g polya_clock.log code/polya_clock.py data/zeros_1700.txt
# -- --out polya_clock.json
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

try:
    from flint import acb, acb_series
    from flint import ctx as fctx
    fctx.prec = 100
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False

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


# ---------------------------------------------------------------- the evaluation of xi and xi'/xi
def _xi_mp(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def _xil_mp(s):
    z = mp.zeta(s)
    xi = s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * z
    L = 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / z
    return xi, L


def _xi_fl(s):
    return s * (s - 1) / 2 * acb.pi() ** (-s / 2) * (s / 2).gamma() * s.zeta()


def _xil_fl(s):
    ser = acb_series([s, 1], prec=2).zeta()
    try:
        c = ser.coeffs()
    except Exception:
        c = [ser[0], ser[1]]
    z, dz = c[0], c[1]
    xi = s * (s - 1) / 2 * acb.pi() ** (-s / 2) * (s / 2).gamma() * z
    L = 1 / s + 1 / (s - 1) - acb.pi().log() / 2 + (s / 2).digamma() / 2 + dz / z
    return xi, L


def XI(s):
    """xi(s) as a Python complex; Re s < 1/2 is reflected by xi(s) = xi(1 - s), and xi(1) = 1/2 (zeta has its pole there)."""
    s = complex(s)
    if s.real < 0.5:
        s = 1 - s
    if abs(s - 1) < 1e-9:
        return 0.5 + 0j
    if HAVE_FLINT:
        return complex(_xi_fl(acb(s.real, s.imag)))
    return complex(_xi_mp(mp.mpc(s.real, s.imag)))


def XIL(s):
    """(xi(s), xi'/xi(s)) as Python complexes; the reflection gives xi'/xi(s) = - xi'/xi(1 - s)."""
    s = complex(s)
    if s.real < 0.5:
        x, L = XIL(1 - s)
        return x, -L
    if HAVE_FLINT:
        a, b = _xil_fl(acb(s.real, s.imag))
    else:
        a, b = _xil_mp(mp.mpc(s.real, s.imag))
    return complex(a), complex(b)


def selftest():
    global HAVE_FLINT
    mp.mp.dps = 25
    pts = [complex(0.75, 31.4), complex(0.5, 7.0), complex(-0.1, 55.5), complex(3.0, 0.0)]
    worst = 0.0
    if HAVE_FLINT:
        try:
            for s in pts:
                b = complex(_xi_mp(mp.mpc(s)))
                b2, l2 = _xil_mp(mp.mpc(s))
                a = XI(s)
                a2, la = XIL(s)
                worst = max(worst, abs(a - b) / abs(b), abs(la - complex(l2)) / max(1.0, abs(complex(l2))),
                            abs(a2 - complex(b2)) / abs(complex(b2)))
            if not worst < 1e-10:
                raise ValueError('flint disagrees with mpmath: %g' % worst)
        except Exception as e:
            log('python-flint failed its self test (%r): using mpmath' % (e,))
            HAVE_FLINT = False
    mp.mp.dps = 20
    log('backend', 'python-flint 100 bits' if HAVE_FLINT else 'mpmath 20 digits', 'self test worst relative difference', fs(worst, 3))


def w_eval(task):
    kind, alpha, pts = task
    out = []
    if kind == 'xi':
        for w in pts:
            out.append(XI(w))
    elif kind == 'A':
        for w in pts:
            out.append((XI(alpha + w) + XI(alpha - w)) / 2)
    elif kind == 'B':
        for w in pts:
            out.append((XI(alpha + w) - XI(alpha - w)) / 2j)
    elif kind == 'dxi':
        for w in pts:
            x, L = XIL(0.5 + w)
            out.append(x * L)
    elif kind == 'L':
        for w in pts:
            out.append(XIL(w)[1])
    else:
        raise ValueError(kind)
    return out


_POOL = None


def get_pool(workers):
    global _POOL
    if workers <= 1:
        return None
    if _POOL is None:
        try:
            _POOL = get_context('fork').Pool(workers)
            atexit.register(_POOL.terminate)
        except Exception as e:
            log('the worker pool failed (%r): running in one process' % (e,))
            _POOL = False
    return _POOL or None


def pmap(kind, alpha, pts, workers, chunk=250):
    pts = [complex(p) for p in pts]
    if not pts:
        return np.zeros(0, complex)
    tasks = [(kind, alpha, pts[i:i + chunk]) for i in range(0, len(pts), chunk)]
    pool = get_pool(workers) if len(tasks) > 1 else None
    res = pool.map(w_eval, tasks, chunksize=1) if pool else [w_eval(t) for t in tasks]
    return np.concatenate([np.asarray(r, complex) for r in res])


# ---------------------------------------------------------------- the clock
def unwrap(z):
    """Clock (1/pi) arg z from unit vectors; z[0] is real and positive. Returns the clock and the largest increment in radians."""
    w = z / np.abs(z)
    d = np.angle(w[1:] * np.conj(w[:-1]))
    return np.concatenate([[0.0], np.cumsum(d)]) / np.pi, float(np.abs(d).max())


def dist_half(U):
    return np.abs((np.asarray(U, float) % 1.0) - 0.5)


def sign_changes(f):
    s = np.sign(f)
    return np.nonzero(s[:-1] * s[1:] < 0)[0]


def lin_zeros(t, f):
    idx = sign_changes(f)
    t0, t1, f0, f1 = t[idx], t[idx + 1], f[idx], f[idx + 1]
    return idx, t0 - f0 * (t1 - t0) / (f1 - f0)


def refine(alpha, a, b, part):
    def g(t):
        z = XI(complex(alpha, t))
        return (z.real if part == 'Re' else z.imag) / abs(z)
    return brentq(g, a, b, xtol=1e-14, rtol=1e-15, maxiter=200)


def clock_exact(alpha, tau, taus, z, U, h):
    """Clock at tau from the unwrapped grid value below it and the principal argument of the short ratio."""
    k = min(int(tau / h), len(taus) - 2)
    w1 = XI(complex(alpha, tau))
    w0 = z[k]
    return U[k] + np.angle((w1 / abs(w1)) * np.conj(w0 / abs(w0))) / np.pi


def clock_by_derivative(alpha, T, workers, lcap, nodes):
    c = alpha - 0.5
    L = min(lcap, 5 * c)
    npan = int(math.ceil(T / L))
    edges = np.linspace(0.0, T, npan + 1)
    x, w = np.polynomial.legendre.leggauss(nodes)
    pts, wts = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        pts.extend(complex(alpha, 0.5 * (b - a) * xi + 0.5 * (a + b)) for xi in x)
        wts.extend(0.5 * (b - a) * wi for wi in w)
    v = pmap('L', 0.0, pts, workers).real / np.pi
    return float(np.dot(v, np.array(wts))), npan


def table_alpha(alpha, T, hf, hg, workers, ordinates, with_im=False, lcap=0.5, nodes=40, deriv=True):
    nf = int(round(T / hf))
    taus = np.arange(nf + 1) * hf
    z = pmap('xi', 0.0, [complex(alpha, t) for t in taus], workers)
    Uf, stepf = unwrap(z)
    sub = int(round(hg / hf))
    assert abs(sub * hf - hg) < 1e-12
    zc, tc = z[::sub], taus[::sub]
    Uc, stepc = unwrap(zc)
    # variant 1: the text's method, grid hg, linear interpolation of Re xi and of the clock
    ic, tstar = lin_zeros(tc, zc.real)
    Ustar = np.interp(tstar, tc, Uc)
    d1 = dist_half(Ustar)
    # variant 2: zeros as variant 1, clock evaluated exactly there
    Uex = np.array([clock_exact(alpha, t, taus, z, Uf, hf) for t in tstar])
    d2 = dist_half(Uex)
    # variant 3: zeros refined by root finding on the fine grid, clock exact
    iff = sign_changes(z.real)
    roots = np.array([refine(alpha, taus[k], taus[k + 1], 'Re') for k in iff])
    Ur = np.array([clock_exact(alpha, t, taus, z, Uf, hf) for t in roots])
    d3 = dist_half(Ur)
    levels = np.rint(Ur + 0.5) - 0.5
    rec = {
        'alpha': fs(alpha), 'T': fs(T),
        'zeros_grid': int(len(ic)), 'zeros_fine_grid': int(len(iff)), 'zeros_refined': int(len(roots)),
        'U_T_grid': fs(Uc[-1]), 'U_T_fine': fs(Uf[-1]), 'round_U_T': int(math.floor(Uf[-1] + 0.5)), 'floor_U_T': int(math.floor(Uf[-1])),
        'maxdist_grid_linear': fs(d1.max(), 5), 'maxdist_grid_zero_exact_clock': fs(d2.max(), 5), 'maxdist_refined': fs(d3.max(), 5),
        'min_increment_grid': fs(np.diff(Uc).min(), 5), 'min_increment_fine': fs(np.diff(Uf).min(), 5),
        'strictly_increasing_grid': bool(np.diff(Uc).min() > 0), 'strictly_increasing_fine': bool(np.diff(Uf).min() > 0),
        'max_unwrap_step_rad_grid': fs(stepc, 5), 'max_unwrap_step_rad_fine': fs(stepf, 5),
        'levels_are_k_minus_half_in_order': bool(np.array_equal(levels, np.arange(len(roots)) + 0.5)),
        'atoms_first3': [fs(t, 10) for t in roots[:3]],
        'atoms_first3_grid_linear': [fs(t, 10) for t in tstar[:3]],
        'ordinates_first3': [fs(t, 10) for t in ordinates[:3]],
        'last_atom': fs(roots[-1], 10) if len(roots) else None,
        'N_T': int(np.sum(ordinates <= T)),
    }
    if deriv:
        Ud, npan = clock_by_derivative(alpha, T, workers, lcap, nodes)
        rec['U_T_by_derivative'] = fs(Ud)
        rec['U_T_derivative_minus_unwrapped'] = fs(Ud - Uf[-1], 3)
        rec['derivative_panels'] = npan
    if with_im:
        ii = sign_changes(z.imag[1:])
        rec['sign_changes_Re'] = int(len(iff))
        rec['sign_changes_Im'] = int(len(ii))
        rec['Im_sign_at_first_grid_point_positive'] = bool(z.imag[1] > 0)
        rec['atoms_B_first3'] = [fs(refine(alpha, taus[1 + k], taus[2 + k], 'Im'), 10) for k in ii[:3]]
    return rec, dict(taus=taus, z=z, Uf=Uf, roots=roots, tstar=tstar)


# ---------------------------------------------------------------- the argument principle on rectangles
def trace(kind, alpha, verts, steps, workers, tol=0.5, maxit=14):
    pts = []
    for i, a in enumerate(verts):
        b = verts[(i + 1) % len(verts)]
        n = max(1, int(math.ceil(abs(b - a) / steps[i])))
        pts.extend(a + (b - a) * k / n for k in range(n))
    pts.append(verts[0])
    P = np.array(pts, complex)
    F = pmap(kind, alpha, P, workers)
    bad = np.zeros(0, int)
    for rounds in range(maxit + 1):
        r = F[1:] / F[:-1]
        bad = np.nonzero((np.abs(np.angle(r)) > tol) | (np.abs(np.log(np.abs(r))) > 1.0))[0]
        if len(bad) == 0 or rounds == maxit:
            break
        mids = 0.5 * (P[bad] + P[bad + 1])
        Fm = pmap(kind, alpha, mids, workers)
        P = np.insert(P, bad + 1, mids)
        F = np.insert(F, bad + 1, Fm)
    d = np.angle(F[1:] / F[:-1])
    return dict(count=float(d.sum() / (2 * math.pi)), npts=int(len(P)), min_abs=float(np.abs(F).min()),
                max_jump=float(np.abs(d).max()), unresolved=int(len(bad)))


def rect_count(kind, alpha, X, y0, T, sv, sh, workers):
    verts = [complex(-X, y0), complex(X, y0), complex(X, T), complex(-X, T)]
    r = trace(kind, alpha, verts, [sh, sv, sh, sv], workers)
    return {'kind': kind, 'alpha': fs(alpha), 'X': fs(X), 'y0': fs(y0), 'T': fs(T), 'count': '%.6f' % r['count'],
            'distance_from_integer': fs(abs(r['count'] - round(r['count'])), 3), 'samples': r['npts'],
            'min_modulus_on_boundary': fs(r['min_abs'], 4), 'max_arg_jump': fs(r['max_jump'], 4), 'unresolved_jumps': r['unresolved']}


# ---------------------------------------------------------------- Theorem exact on the real axis
def xi_real(x):
    x = mp.mpf(x)
    if x < 0.5:
        x = 1 - x
    if abs(x - 1) < mp.mpf(10) ** -25:
        return mp.mpf(1) / 2
    return x * (x - 1) / 2 * mp.pi ** (-x / 2) * mp.gamma(x / 2) * mp.zeta(x)


def L_real(r):
    r = mp.mpf(r)
    return 1 / r + 1 / (r - 1) - mp.log(mp.pi) / 2 + mp.digamma(r / 2) / 2 + mp.zeta(r, derivative=1) / mp.zeta(r)


def theorem_exact(alphas):
    mp.mp.dps = 30
    out = {}
    for s in (1, 100):
        Ls = L_real(mp.mpf(1) / 2 + mp.sqrt(s))
        row = {'L_half_plus_sqrt_s': mp.nstr(Ls, 12), 'alphas': {}}
        for al in alphas:
            a = mp.mpf(al) - mp.mpf(1) / 2
            R = xi_real(mp.mpf(al) - mp.sqrt(s)) / xi_real(mp.mpf(al) + mp.sqrt(s))
            defect = mp.log((1 + R) / 2)
            pred = -a * Ls
            row['alphas'][fs(al)] = {'R': mp.nstr(R, 12), 'defect': mp.nstr(defect, 12), 'defect_over_a': mp.nstr(defect / a, 12),
                                     'prediction_minus_a_L': mp.nstr(pred, 12), 'defect_over_prediction': mp.nstr(defect / pred, 8),
                                     'relative_gap_to_prediction': mp.nstr((defect - pred) / pred, 4)}
        out['s=%d' % s] = row
    # 0 < R < 1 on a grid, and R at s = 0
    bad, worst_hi, worst_lo, n = [], mp.mpf(0), mp.mpf(1), 0
    for al in (1, 0.75, 0.6, 0.55, 0.52, 0.505, 0.501):
        for s in (1e-6, 1e-4, 0.01, 0.1, 0.5, 1, 2, 5, 10, 25, 50, 100, 400, 2500, 10000):
            R = xi_real(mp.mpf(al) - mp.sqrt(s)) / xi_real(mp.mpf(al) + mp.sqrt(s))
            n += 1
            worst_hi, worst_lo = max(worst_hi, R), min(worst_lo, R)
            if not 0 < R < 1:
                bad.append((al, s, mp.nstr(R, 8)))
    out['R_in_open_unit_interval'] = {'points': n, 'violations': bad, 'largest_R': mp.nstr(worst_hi, 12), 'smallest_R': mp.nstr(worst_lo, 6),
                                      'R_at_s_0': mp.nstr(xi_real(mp.mpf('0.75')) / xi_real(mp.mpf('0.75')), 5)}
    mp.mp.dps = 20
    return out


# ---------------------------------------------------------------- main
def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--zeros', default=None)
    a = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    out_path = Path(a.out) if a.out else (Path('/tmp/polya_clock_quick.json') if a.quick else root / 'data' / 'polya_clock.json')
    LOGF = open(out_path.with_suffix('.log'), 'w')
    workers = 1 if a.quick else a.workers
    selftest()
    zf = None
    for cand in ([Path(a.zeros)] if a.zeros else []) + [Path('zeros_1700.txt'), root / 'data' / 'zeros_1700.txt']:
        if cand.exists():
            zf = cand
            break
    ordinates = np.array([float(x) for x in zf.read_text().split()])
    log('workers', workers, 'zeros file', zf, 'ordinates', len(ordinates), 'N(100) =', int(np.sum(ordinates <= 100)),
        'N(300) =', int(np.sum(ordinates <= 300)), 'gamma_138 =', fs(ordinates[137], 10), 'gamma_139 =', fs(ordinates[138], 10))
    if a.quick:
        P = dict(T=15.0, hf=0.02, hg=0.04, alphas=[1.0, 0.75], along=0.75, Tlong=15.0, sv=0.4, sh=0.2, lcap=2.0, nodes=16,
                 xp_lo=1.0, xp_hi=15.0, xp_step=0.1, figT=15.0)
    else:
        P = dict(T=100.0, hf=0.0025, hg=0.01, alphas=[1.0, 0.75, 0.6, 0.52], along=0.55, Tlong=300.0, sv=0.02, sh=0.01, lcap=0.5,
                 nodes=40, xp_lo=1.0, xp_hi=100.0, xp_step=0.01, figT=60.0)
    R = {'meta': {'quick': a.quick, 'backend': 'python-flint 100 bits' if HAVE_FLINT else 'mpmath 20 digits', 'workers': workers,
                  'parameters': {k: (fs(v) if isinstance(v, float) else v) for k, v in P.items()}}}
    chapter = {'1': ('29', '29.1112', '29', '2.7e-5'), '0.75': ('29', '29.0548', '29', '5.6e-5'),
               '0.6': ('29', '29.0218', '29', '1.8e-4'), '0.52': ('29', '29.0044', '29', '1.6e-3')}

    # ---- Table pol:tab:clock, the atoms, the monotonicity, the figure series
    R['table'] = {}
    extras = {}
    for al in P['alphas']:
        rec, ex = table_alpha(al, P['T'], P['hf'], P['hg'], workers, ordinates, lcap=P['lcap'], nodes=P['nodes'])
        R['table'][fs(al)] = rec
        extras[al] = ex
        ch = chapter.get(fs(al), ('', '', '', ''))
        log('alpha', fs(al), 'zeros grid/fine/refined', rec['zeros_grid'], rec['zeros_fine_grid'], rec['zeros_refined'], '(chapter', ch[0] + ')',
            'U(T)', rec['U_T_grid'], rec['U_T_fine'], 'derivative route', rec['U_T_by_derivative'], '(chapter', ch[1] + ')', 'round',
            rec['round_U_T'], 'maxdist grid/exact clock/refined', rec['maxdist_grid_linear'], rec['maxdist_grid_zero_exact_clock'],
            rec['maxdist_refined'], '(chapter', ch[3] + ')', 'increasing', rec['strictly_increasing_grid'], rec['strictly_increasing_fine'],
            'min increment', rec['min_increment_grid'], 'max unwrap step', rec['max_unwrap_step_rad_fine'], 'levels ok',
            rec['levels_are_k_minus_half_in_order'])
        log('   atoms', rec['atoms_first3'], 'grid linear atoms', rec['atoms_first3_grid_linear'], 'ordinates', rec['ordinates_first3'])
    # ---- Figure pol:fig:one, left series (alpha = 3/4)
    ex = extras[0.75]
    step = int(round(0.02 / P['hf']))
    nfig = int(round(P['figT'] / 0.02))
    taus_f = ex['taus'][::step][:nfig + 1]
    U_f = ex['Uf'][::step][:nfig + 1]
    zr = ex['roots'][ex['roots'] <= P['figT']]
    R['figure1_left'] = {'alpha': '0.75', 'step': '0.02', 'tau': [fs(t, 6) for t in taus_f], 'U': [fs(u, 12) for u in U_f],
                         'zeros_of_Re_xi': [fs(t, 12) for t in zr], 'levels': [fs(k + 0.5, 3) for k in range(len(zr))],
                         'ordinates_below_height': [fs(t, 12) for t in ordinates[ordinates <= P['figT']]],
                         'N_at_height': int(np.sum(ordinates <= P['figT']))}
    log('figure 1 left: %d clock values, %d crossings, N(%s) = %d' % (len(U_f), len(zr), fs(P['figT']), int(np.sum(ordinates <= P['figT']))))

    # ---- the long run at alpha = 0.55 (height 300): sign changes of Re and Im, the clock
    al = P['along']
    rec, ex = table_alpha(al, P['Tlong'], P['hf'], P['hg'], workers, ordinates, with_im=True, lcap=P['lcap'], nodes=P['nodes'])
    R['long'] = rec
    log('long run alpha', fs(al), 'T', fs(P['Tlong']), 'sign changes Re/Im', rec['sign_changes_Re'], rec['sign_changes_Im'], 'U(T)', rec['U_T_fine'],
        'derivative route', rec['U_T_by_derivative'], 'round/floor', rec['round_U_T'], rec['floor_U_T'], 'N(T)', rec['N_T'],
        'increasing', rec['strictly_increasing_fine'], 'maxdist refined', rec['maxdist_refined'], 'first B zeros', rec['atoms_B_first3'])

    # ---- argument principle counts
    R['argument_principle'] = []
    jobs = []
    if a.quick:
        jobs = [('A', 1.0, 1.0, 0.0, P['T']), ('A', 0.75, 2.0, 0.0, P['T']), ('B', 0.75, 1.0, 0.5, P['T']),
                ('B', P['along'], 2.0, 0.5, P['Tlong'])]
    else:
        for al_ in P['alphas']:
            for X in (1.0, 2.0):
                jobs.append(('A', al_, X, 0.0, P['T']))
        for al_ in (0.6, 0.75):
            jobs.append(('B', al_, 1.0, 0.5, P['T']))
            jobs.append(('B', al_, 2.0, 0.5, P['T']))
        for X in (1.0, 2.0):
            jobs.append(('A', P['along'], X, 0.0, P['Tlong']))
            jobs.append(('B', P['along'], X, 0.5, P['Tlong']))
    for kind, al_, X, y0, T in jobs:
        rc = rect_count(kind, al_, X, y0, T, P['sv'], P['sh'], workers)
        R['argument_principle'].append(rc)
        log('argument principle', kind, 'alpha', rc['alpha'], '|Re w| <=', rc['X'], rc['y0'], '<= Im w <=', rc['T'], 'count', rc['count'], 'dist',
            rc['distance_from_integer'], 'samples', rc['samples'], 'min |f|', rc['min_modulus_on_boundary'], 'unresolved', rc['unresolved_jumps'])

    # ---- the zeros of xi'(1/2 + w)
    t = np.arange(int(round((P['xp_hi'] - P['xp_lo']) / P['xp_step'])) + 1) * P['xp_step'] + P['xp_lo']
    t_low = np.arange(1, 100) * 0.01
    tt = np.concatenate([t_low[t_low < P['xp_lo'] - 1e-9], t])
    v = -1j * pmap('dxi', 0.0, [1j * x for x in tt], workers)
    ic = sign_changes(v.real)
    def hnorm(x):
        xv, Lv = XIL(complex(0.5, x))
        q = -1j * xv * Lv
        return q.real / abs(q)
    roots = np.array([brentq(hnorm, tt[k], tt[k + 1], xtol=1e-14, rtol=1e-15) for k in ic])
    gaps = np.searchsorted(ordinates, roots)
    per_gap = {}
    for g_ in gaps:
        per_gap[int(g_)] = per_gap.get(int(g_), 0) + 1
    inside = roots[(roots > P['xp_lo']) & (roots <= P['xp_hi'])]
    xp_box = rect_count('dxi', 0.0, 0.6, P['xp_lo'], P['xp_hi'], P['sv'], P['sh'], workers)
    R['xi_prime'] = {'grid': [fs(tt[0]), fs(tt[-1]), fs(P['xp_step'])], 'sign_changes_of_minus_i_xi_prime_on_line_in_range': int(len(inside)),
                     'sign_changes_below_lower_end': int(np.sum(roots <= P['xp_lo'])), 'max_relative_imaginary_part_of_minus_i_xi_prime':
                     fs(np.max(np.abs(v.imag) / np.abs(v)), 3), 'zeros_first5': [fs(x, 12) for x in inside[:5]],
                     'zeros_last3': [fs(x, 12) for x in inside[-3:]], 'zeros_per_gap_of_ordinates_index_from_0': {str(k): c for k, c in sorted(per_gap.items())},
                     'gaps_between_consecutive_ordinates_with_exactly_one_zero': int(sum(1 for k, c in per_gap.items() if 1 <= k < int(np.sum(ordinates <= P['xp_hi'])) and c == 1)),
                     'zeros_before_first_ordinate_after_lower_end': int(np.sum((inside < ordinates[0]))),
                     'zeros_after_last_ordinate_below_top': int(np.sum(inside > ordinates[ordinates <= P['xp_hi']].max())),
                     'argument_principle_box': xp_box}
    log('xi prime: sign changes in (%s,%s]' % (fs(P['xp_lo']), fs(P['xp_hi'])), len(inside), 'below', int(np.sum(roots <= P['xp_lo'])),
        'max |Im|/|.|', R['xi_prime']['max_relative_imaginary_part_of_minus_i_xi_prime'], 'first', R['xi_prime']['zeros_first5'], 'last',
        R['xi_prime']['zeros_last3'])
    log('xi prime: zeros per gap of the ordinates (index 0 = before gamma_1)', R['xi_prime']['zeros_per_gap_of_ordinates_index_from_0'])
    log('xi prime: argument principle box |Re w| <= 0.6', xp_box['count'], 'dist', xp_box['distance_from_integer'], 'samples', xp_box['samples'],
        'min |g|', xp_box['min_modulus_on_boundary'], 'unresolved', xp_box['unresolved_jumps'])

    # ---- Theorem exact
    te = theorem_exact([1.0, 0.75, 0.6, 0.52])
    R['theorem_exact'] = te
    chd = {'s=1': {'1': '-0.02279', '0.75': '-0.01147', '0.6': '-0.0046', '0.52': '-0.00092'},
           's=100': {'1': '-0.1832', '0.75': '-0.0967', '0.6': '-0.0399', '0.52': '-0.0081'}}
    for sk in ('s=1', 's=100'):
        log('theorem exact', sk, 'L(1/2 + sqrt s) =', te[sk]['L_half_plus_sqrt_s'])
        for k, rec in te[sk]['alphas'].items():
            log('   alpha', k, 'R', rec['R'], 'defect', rec['defect'], '(chapter', chd[sk].get(k, '?') + ')', 'defect/a', rec['defect_over_a'],
                'prediction -a L', rec['prediction_minus_a_L'], 'ratio', rec['defect_over_prediction'])
    log('0 < R < 1 on', te['R_in_open_unit_interval']['points'], 'points: violations', te['R_in_open_unit_interval']['violations'], 'largest R',
        te['R_in_open_unit_interval']['largest_R'], 'smallest R', te['R_in_open_unit_interval']['smallest_R'], 'R at s = 0',
        te['R_in_open_unit_interval']['R_at_s_0'])

    R['meta']['seconds'] = fs(time.time() - T0, 4)
    out_path.write_text(json.dumps(R, indent=1))
    log('wrote', out_path)


if __name__ == '__main__':
    main()
