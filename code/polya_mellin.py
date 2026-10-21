# polya_mellin.py: Polya's Mellin criterion applied to the Riemann kernel, Chapter 16 (Polya and the Zeros of Entire Functions), Section
# pol:sec:mellin. With Phi(u) = 2 sum_n (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)) and M(z) = int_0^oo Phi(u) u^(2z-1)
# du it checks: M(k + 1/2) = (1/2) (2k)! c_k > 0 with xi(1/2 + w) = sum c_k w^(2k), and 2 M(1/2) = xi(1/2); the zeros of M in the right
# half plane by the argument principle (the chapter states 2 zeros in [0.3,6] x [0.2,4], 7 up to height 6, 22 up to height 12) and the first
# zero z = 0.36822 + 3.74341 i. Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were
# neither read nor run).
# Method. For Re z < 1/2 the integrand has the integrable singularity u^(2 Re z - 1) at u = 0 and oscillates like cos(2 Im z log u) there, so
# plain Gauss Legendre quadrature fails. Phi is even and analytic in the strip |Im u| < pi/4 (the series needs Re e^(2u) > 0), hence has a
# Taylor series sum_j p_j u^(2j) of radius pi/4, and the piece from 0 to delta = 0.25 is integrated termwise exactly:
# int_0^delta Phi(u) u^(2z-1) du = sum_j p_j delta^(2z+2j)/(2z+2j) (36 terms; p_j from the Cauchy integral on |u| = 1/2 with 256 points at 50 digits,
# aliasing error about 1e-50). The rest, delta to 2.1 (Phi(2.1) is below 1e-90), is Gauss Legendre with 48 nodes on each of seven panels.
# Everything is vectorised double precision, M' (for Newton's method) uses the same nodes with the factor 2 log u. The argument principle
# counts the zeros inside a rectangle as (1/2 pi) times the change of arg M along its boundary, sampled with step 0.01 and bisected wherever
# consecutive samples differ in argument by more than 0.5 rad or in modulus by more than a factor e. To find all zeros the region
# 0.05 <= Re z <= 8.05, 0.05 <= Im z <= 14.05 is cut into cells of side 1/2, the winding number of every cell is computed and Newton's method
# is started inside every cell with a nonzero winding number (the number of distinct zeros found must equal the winding number). The counts of
# the chapter's rectangles are then read from the zero list and recomputed by the argument principle on those rectangles.
# Error control. (a) M(z) is compared at several points with three other evaluations in mpmath: the substitution u = e^x (tanh sinh on unit
# subintervals of x from the left end where e^(2 Re z x) < 10^-(dps+6) to log 2.1), at 20 digits; for 0 < Re z < 1/2 the Fourier route
# M(z) = (1/pi) Gamma(2z) cos(pi z) int_0^oo Xi(t) t^(-2z) dt (derived from Phi(u) = (1/pi) int_0^oo Xi(t) cos(tu) dt, with the Taylor series of Xi on
# (0,1) and quadrature on (1,80), Xi from mpmath's zeta); and the exact values at z = k + 1/2 from the Taylor coefficients of xi at 1/2
# (Cauchy integral, 40 digits). (b) A second Gauss Legendre layout (delta = 0.35, 64 nodes, other panels) is compared on 400 points. (c) If a zero
# is found it is polished with mpmath (secant method on the substitution evaluation at 22 digits); the value of M at the chapter's first zero
# 0.36822 + 3.74341 i is recorded, and Newton's method on the plain Gauss Legendre sum (n nodes on [0,U], the singularity untreated) is started
# there to show where such a quadrature puts its roots. The same cell scan is run on the continuation of M to -3.95 <= Re z <= 0.05 (the formula
# above is the meromorphic continuation, with poles at z = 0, -1, -2, ...). (d) The conditioning kappa = (sum of the absolute
# terms)/|M| on every counted boundary is recorded (the relative accuracy of M is about 1e-16 kappa), together with the distance of every count
# from an integer and the number of unresolved jumps. Nothing here is a proof: the counts are floating point.
# Output: data/polya_mellin.json (values as strings) and data/polya_mellin.log, both next to --out (the log has the same name with the extension .log).
# Usage: python3 polya_mellin.py [--quick] [--out FILE] [--workers W]; --quick checks a few points and one rectangle (about 2 s on one core) and
# writes to /tmp unless --out is given; the full case is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8 -t 60 -g polya_mellin.json
# -g polya_mellin.log code/polya_mellin.py -- --out polya_mellin.json
import argparse
import atexit
import json
import math
import time
from multiprocessing import get_context
from pathlib import Path

import mpmath as mp
import numpy as np

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


def fz(z, n=12):
    return '%s%s%si' % (fs(z.real, n), '+' if z.imag >= 0 else '-', fs(abs(z.imag), n))


# ---------------------------------------------------------------- Phi and xi in mpmath
def phi_mp(u):
    u = mp.mpmathify(u)
    e2 = mp.exp(2 * u)
    s = 0
    for n in range(1, 120):
        term = (2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * u / 2)) * mp.exp(-mp.pi * n ** 2 * e2)
        s += term
        if n > 2 and abs(term) <= mp.mpf(10) ** (-mp.mp.dps - 8) * abs(s):
            break
    return 2 * s


def xi_mp(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def taylor_phi(J, r=0.5, N=256):
    mp.mp.dps = 50
    rr = mp.mpf(r)
    vals = [phi_mp(rr * mp.expj(2 * mp.pi * m / N)) for m in range(N)]
    out = []
    for j in range(J):
        k = 2 * j
        s = sum(vals[m] * mp.expj(-2 * mp.pi * m * k / N) for m in range(N)) / N
        out.append(float(mp.re(s) / rr ** k))
    return out


def xi_taylor(K, r=2, N=128):
    mp.mp.dps = 40
    vals = [xi_mp(mp.mpf(1) / 2 + r * mp.expj(2 * mp.pi * m / N)) for m in range(N)]
    out = []
    for k in range(K):
        n = 2 * k
        s = sum(vals[m] * mp.expj(-2 * mp.pi * m * n / N) for m in range(N)) / N
        out.append(mp.re(s) / mp.mpf(r) ** n)
    return out


def phi_np(u):
    u = np.asarray(u, float)
    s = np.zeros_like(u)
    for n in range(1, 9):
        s += (2 * np.pi ** 2 * n ** 4 * np.exp(4.5 * u) - 3 * np.pi * n ** 2 * np.exp(2.5 * u)) * np.exp(-np.pi * n ** 2 * np.exp(2 * u))
    return 2 * s


class Mellin:
    def __init__(self, delta=0.25, J=36, edges=(0.25, 0.35, 0.5, 0.7, 0.95, 1.3, 1.7, 2.1), nodes=48, P=None):
        self.delta = delta
        self.logd = math.log(delta)
        self.P = np.array(P if P is not None else taylor_phi(J), float)
        self.J = len(self.P)
        x, w = np.polynomial.legendre.leggauss(nodes)
        U, W = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            U.append(0.5 * (b - a) * x + 0.5 * (a + b))
            W.append(0.5 * (b - a) * w)
        self.U = np.concatenate(U)
        self.W = np.concatenate(W)
        self.base = phi_np(self.U) * self.W / self.U
        self.logU = np.log(self.U)
        self.k2 = 2.0 * np.arange(self.J)

    def M(self, z, deriv=False, cond=False, step=1500):
        z = np.atleast_1d(np.asarray(z, complex))
        out = np.empty(z.shape, complex)
        outd = np.empty(z.shape, complex)
        kap = np.empty(z.shape, float)
        for i in range(0, len(z), step):
            zz = z[i:i + step, None]
            E = np.exp(2 * zz * self.logU[None, :])
            s = 2 * zz + self.k2[None, :]
            Tt = self.P[None, :] * np.exp(s * self.logd) / s
            val = (E * self.base[None, :]).sum(1) + Tt.sum(1)
            out[i:i + step] = val
            if deriv:
                dT = self.P[None, :] * np.exp(s * self.logd) * (2 * self.logd / s - 2 / s ** 2)
                outd[i:i + step] = (E * (2 * self.logU[None, :]) * self.base[None, :]).sum(1) + dT.sum(1)
            if cond:
                Ea = np.exp(2 * zz.real * self.logU[None, :])
                kap[i:i + step] = ((Ea * np.abs(self.base)[None, :]).sum(1) + np.abs(Tt).sum(1)) / np.abs(val)
        res = [out]
        if deriv:
            res.append(outd)
        if cond:
            res.append(kap)
        return res[0] if len(res) == 1 else tuple(res)


# ---------------------------------------------------------------- independent evaluations (mpmath)
def M_logvar(z, dps):
    mp.mp.dps = dps
    z = mp.mpc(z)
    xhi = mp.log(mp.mpf('2.1'))
    X0 = int((dps + 6) * 2.303 / (2 * float(z.real))) + 3
    n = int(2 * (xhi + X0)) + 1
    edges = [-X0 + (xhi + X0) * i / n for i in range(n + 1)]
    return mp.quad(lambda x: phi_mp(mp.exp(x)) * mp.exp(2 * z * x), edges)


def M_fourier(z, ck, dps):
    mp.mp.dps = dps
    z = mp.mpc(z)
    ser = mp.mpf(0)
    for k, c in enumerate(ck):
        ser += (-1) ** k * c / (2 * k + 1 - 2 * z)
    pts = [1, 2, 3.5, 5, 7, 9.5, 12, 15, 18, 21.5, 25, 29, 33, 37, 41, 46, 51, 56, 62, 68, 74, 80]
    integ = mp.quad(lambda t: mp.re(xi_mp(mp.mpf(1) / 2 + 1j * t)) * t ** (-2 * z), pts)
    return mp.gamma(2 * z) * mp.cos(mp.pi * z) * (ser + integ) / mp.pi


def w_indep(task):
    kind, z, dps, ck = task
    t = time.time()
    v = M_logvar(z, dps) if kind == 'logvar' else M_fourier(z, ck, dps)
    return kind, z, complex(v), mp.nstr(v, 16), time.time() - t


# ---------------------------------------------------------------- the argument principle
def trace_M(Mf, verts, step=0.01, tol=0.5, maxit=14):
    pts = []
    for i, a in enumerate(verts):
        b = verts[(i + 1) % len(verts)]
        n = max(1, int(math.ceil(abs(b - a) / step)))
        pts.extend(a + (b - a) * k / n for k in range(n))
    pts.append(verts[0])
    P = np.array(pts, complex)
    F = Mf(P)
    bad = np.zeros(0, int)
    for rounds in range(maxit + 1):
        r = F[1:] / F[:-1]
        bad = np.nonzero((np.abs(np.angle(r)) > tol) | (np.abs(np.log(np.abs(r))) > 1.0))[0]
        if len(bad) == 0 or rounds == maxit or len(P) > 150000:
            break
        mids = 0.5 * (P[bad] + P[bad + 1])
        P = np.insert(P, bad + 1, mids)
        F = np.insert(F, bad + 1, Mf(mids))
    d = np.angle(F[1:] / F[:-1])
    return dict(count=float(d.sum() / (2 * math.pi)), npts=len(P), min_abs=float(np.abs(F).min()), unresolved=int(len(bad))), P


def rect_M(Mo, x0, x1, y0, y1, step=0.01, with_cond=True):
    verts = [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)]
    r, P = trace_M(Mo.M, verts, step)
    rec = {'rectangle': [fs(x0), fs(x1), fs(y0), fs(y1)], 'count': '%.6f' % r['count'], 'distance_from_integer': fs(abs(r['count'] - round(r['count'])), 3),
           'samples': r['npts'], 'min_modulus_on_boundary': fs(r['min_abs'], 4), 'unresolved_jumps': r['unresolved']}
    if with_cond:
        rec['max_condition_kappa'] = fs(Mo.M(P, cond=True)[1].max(), 4)
    return rec, r['count']


def newton(Mo, z0, iters=80, rmin=0.02):
    z = np.array(z0, complex)
    for _ in range(iters):
        f, df = Mo.M(z, deriv=True)
        step = f / df
        big = np.abs(step) > 0.3
        step[big] = step[big] / np.abs(step[big]) * 0.3
        z = z - step
        z = np.where(z.real < rmin, rmin + 1j * z.imag, z)
    f, df = Mo.M(z, deriv=True)
    return z, np.abs(f / df)


def scan_cells(Mo, xe, ye, rmin):
    """Winding number of every cell of the grid xe x ye; Newton's method inside cells with a nonzero winding number."""
    zeros, wsum, warn = [], 0, []
    for i in range(len(xe) - 1):
        for j in range(len(ye) - 1):
            rec, cnt = rect_M(Mo, xe[i], xe[i + 1], ye[j], ye[j + 1], 0.01, with_cond=False)
            w = int(round(cnt))
            if abs(cnt - w) > 0.05 or rec['unresolved_jumps']:
                warn.append((xe[i], ye[j], rec['count'], rec['unresolved_jumps']))
            if w <= 0:
                continue
            wsum += w
            starts = [complex(xe[i] + fx * 0.5, ye[j] + fy * 0.5) for fx in (0.2, 0.5, 0.8) for fy in (0.2, 0.5, 0.8)]
            zz, err = newton(Mo, starts, rmin=rmin)
            found = []
            for z, e in zip(zz, err):
                if e < 1e-9 and xe[i] - 1e-9 <= z.real <= xe[i + 1] + 1e-9 and ye[j] - 1e-9 <= z.imag <= ye[j + 1] + 1e-9:
                    if all(abs(z - q) > 1e-6 for q in found):
                        found.append(z)
            if len(found) != w:
                warn.append(('cell', xe[i], ye[j], 'winding', w, 'zeros found', len(found)))
            zeros.extend(found)
    zeros.sort(key=lambda z: (z.imag, z.real))
    return zeros, wsum, warn


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    out_path = Path(a.out) if a.out else (Path('/tmp/polya_mellin_quick.json') if a.quick else root / 'data' / 'polya_mellin.json')
    LOGF = open(out_path.with_suffix('.log'), 'w')
    workers = 1 if a.quick else a.workers
    R = {'meta': {'quick': a.quick, 'workers': workers, 'mpmath': mp.__version__}}
    Mo = Mellin()
    log('Taylor coefficients of Phi at 0 (p_0, p_1, p_2, p_10, p_30):', [fs(Mo.P[j], 10) for j in (0, 1, 2, 10, 30)], 'M(1/2) =', fs(Mo.M([0.5])[0].real, 15))

    # ---- (a) the values at z = k + 1/2
    ck = xi_taylor(30)
    K = 4 if a.quick else 7
    rows = []
    for k in range(K):
        val = Mo.M([k + 0.5])[0]
        exact = mp.factorial(2 * k) * ck[k] / 2
        rows.append({'k': k, 'M_quadrature': fs(val.real, 15), 'imag_part': fs(val.imag, 3), 'half_factorial_2k_c_k': mp.nstr(exact, 16),
                     'relative_difference': fs(abs(val.real - float(exact)) / float(exact), 3), 'positive': bool(val.real > 0)})
    R['values_at_k_plus_half'] = rows
    R['xi_half'] = mp.nstr(xi_mp(mp.mpf(1) / 2), 16)
    R['two_M_half_minus_xi_half'] = fs(2 * Mo.M([0.5])[0].real - float(xi_mp(mp.mpf(1) / 2)), 3)
    for r_ in rows:
        log('M(k + 1/2): k =', r_['k'], 'quadrature', r_['M_quadrature'], '(1/2)(2k)! c_k', r_['half_factorial_2k_c_k'], 'relative difference', r_['relative_difference'], 'positive', r_['positive'])
    log('2 M(1/2) - xi(1/2) =', R['two_M_half_minus_xi_half'], ' xi(1/2) =', R['xi_half'])

    # ---- (b) independent evaluations
    zs_chk = [complex(0.36822, 3.74341), complex(0.4, 1.0), complex(0.45, 7.5), complex(1.5, 2.0), complex(3.0, 8.0)]
    if a.quick:
        zs_chk = [complex(0.4, 1.0), complex(1.5, 2.0)]
    dps = 15 if a.quick else 22
    tasks = [('logvar', z, dps, None) for z in zs_chk] + ([] if a.quick else [('fourier', z, 20, ck) for z in zs_chk if 0 < z.real < 0.5])
    if workers > 1:
        pool = get_context('fork').Pool(workers)
        atexit.register(pool.terminate)
        res = pool.map(w_indep, tasks, chunksize=1)
    else:
        res = [w_indep(t) for t in tasks]
    ind = []
    for kind, z, v, vs, sec in res:
        m = Mo.M([z])[0]
        ind.append({'method': kind, 'z': fz(z, 8), 'M_gauss_legendre': fz(m, 14), 'M_independent': vs, 'relative_difference': fs(abs(m - v) / abs(v), 3), 'seconds': fs(sec, 3)})
        log('independent', kind, 'z =', fz(z, 8), 'Gauss Legendre', fz(m, 14), 'independent', vs[:44], 'relative difference', fs(abs(m - v) / abs(v), 3), '(%.0f s)' % sec)
    R['independent_evaluations'] = ind

    # ---- (c) consistency of two quadrature layouts
    Mo2 = Mellin(delta=0.35, edges=(0.35, 0.45, 0.6, 0.8, 1.05, 1.4, 1.8, 2.2), nodes=64, J=48)
    rng = np.random.default_rng(1)
    zr = rng.uniform(0.3, 6.0, 400) + 1j * rng.uniform(0.2, 12.0, 400)
    m1, m2 = Mo.M(zr), Mo2.M(zr)
    kap = Mo.M(zr, cond=True)[1]
    rel = np.abs(m1 - m2) / np.abs(m1)
    R['two_layouts'] = {'points': 400, 'max_relative_difference': fs(rel.max(), 3), 'median_relative_difference': fs(np.median(rel), 3), 'max_kappa': fs(kap.max(), 3)}
    log('two quadrature layouts on 400 random points in [0.3,6] x [0.2,12]: max relative difference', fs(rel.max(), 3), 'median', fs(np.median(rel), 3), 'max kappa', fs(kap.max(), 3))

    # ---- (d) the zeros: cells, Newton
    if a.quick:
        xe = [0.05 + 0.5 * i for i in range(3)]
        ye = [0.05 + 0.5 * j for j in range(9)]
    else:
        xe = [0.05 + 0.5 * i for i in range(17)]
        ye = [0.05 + 0.5 * j for j in range(29)]
    zeros, wsum, warn = scan_cells(Mo, xe, ye, 0.02)
    log('cell scan: %d cells, total winding %d, zeros found %d, warnings %s' % ((len(xe) - 1) * (len(ye) - 1), wsum, len(zeros), warn))
    big, _ = rect_M(Mo, xe[0], xe[-1], ye[0], ye[-1], 0.01, with_cond=False)
    log('the whole region', big['rectangle'], 'argument principle count', big['count'], 'distance', big['distance_from_integer'], 'unresolved', big['unresolved_jumps'])
    R['cell_scan'] = {'region': [fs(xe[0]), fs(xe[-1]), fs(ye[0]), fs(ye[-1])], 'cells': (len(xe) - 1) * (len(ye) - 1), 'sum_of_cell_windings': wsum,
                      'zeros_found': len(zeros), 'warnings': [str(w) for w in warn], 'whole_region_count': big['count']}
    R['zeros'] = [{'z': fz(z, 12), 're': fs(z.real, 12), 'im': fs(z.imag, 12)} for z in zeros]
    for z in zeros[:30]:
        log('   zero', fz(z, 12))
    log('zeros with Re z below 0.3:', [fz(z, 6) for z in zeros if z.real < 0.3], ' with Re z above 6:', [fz(z, 6) for z in zeros if z.real > 6])
    # counts in the chapter's rectangles from the list, and by the argument principle
    def cnt(x0, x1, y0, y1):
        return int(sum(1 for z in zeros if x0 < z.real < x1 and y0 < z.imag < y1))
    rects = [('[0.3,6] x [0.2,4]', 0.3, 6.0, 0.2, 4.0, 2), ('[0.3,6] x [0.2,6]', 0.3, 6.0, 0.2, 6.0, 7), ('[0.3,6] x [0.2,12]', 0.3, 6.0, 0.2, 12.0, 22),
             ('[0.3,6] x (0,6]', 0.3, 6.0, 0.0, 6.0, 7), ('[0.3,6] x (0,12]', 0.3, 6.0, 0.0, 12.0, 22), ('Re z > 0, Im z in (0,6]', 0.0, 99.0, 0.0, 6.0, 7),
             ('Re z > 0, Im z in (0,12]', 0.0, 99.0, 0.0, 12.0, 22)]
    if a.quick:
        rects = rects[:1]
    R['rectangles'] = []
    for name, x0, x1, y0, y1, claimed in rects:
        inlist = cnt(x0, x1, y0, y1)
        entry = {'name': name, 'chapter_claims': claimed, 'zeros_in_list': inlist}
        if x0 > 0 and y0 > 0 and y1 <= 14 and x1 <= 8:
            rec, c = rect_M(Mo, x0, x1, y0, y1, 0.01)
            entry['argument_principle'] = rec
        log('rectangle', name, 'chapter', claimed, 'zeros in list', inlist, 'argument principle', entry.get('argument_principle', {}).get('count'), entry.get('argument_principle', {}).get('max_condition_kappa'))
        R['rectangles'].append(entry)
    byheight = {fs(h): int(sum(1 for z in zeros if z.imag < h)) for h in (4, 6, 8, 10, 12, 14)}
    R['zeros_below_height_any_real_part'] = byheight
    log('zeros with 0.05 < Re z < 8.05 below each height:', byheight)

    # ---- (e) the chapter's first zero, and the first zero found (polished with mpmath)
    zc = complex(0.36822, 3.74341)
    mc = Mo.M([zc])[0]
    nb, kap_c = Mo.M([zc], cond=True)
    R['chapter_first_zero'] = {'z': fz(zc, 8), 'M_there': fz(mc, 12), 'abs_M_there': fs(abs(mc), 6), 'abs_M_at_1': fs(abs(Mo.M([1.0])[0]), 6),
                               'kappa': fs(kap_c[0], 4), 'nearest_zero_found': fz(min(zeros, key=lambda z: abs(z - zc)), 8) if zeros else None}
    log('the chapter first zero z =', fz(zc, 8), ': M(z) =', fz(mc, 12), ' |M| =', fs(abs(mc), 6), ' (|M(1)| =', fs(abs(Mo.M([1.0])[0]), 6), ') nearest zero found:', R['chapter_first_zero']['nearest_zero_found'])
    if zeros:
        z1 = zeros[0]
        R['first_zero_double'] = fz(z1, 14)
        log('first zero found (double precision Newton):', fz(z1, 14))
        if not a.quick:
            mp.mp.dps = 22
            zp = mp.findroot(lambda w: M_logvar(w, 22), mp.mpc(z1), solver='secant', tol=1e-17, maxsteps=14, verify=False)
            R['first_zero_mpmath'] = mp.nstr(zp, 18)
            log('first zero polished with mpmath at 22 digits:', mp.nstr(zp, 18), ' difference from the double precision value', mp.nstr(abs(zp - mp.mpc(z1)), 3))
    else:
        R['first_zero_double'] = None
        log('no zero found in the scanned region of the right half plane')

    # ---- (f) why a plain quadrature shows zeros there: Newton's method on the plain Gauss Legendre sum over [0, U] in u
    naive = []
    for U in (2.0, 3.0, 4.0):
        for n in ((100, 256, 500) if a.quick else (50, 100, 128, 200, 256, 400, 500, 512, 1000)):
            x, w = np.polynomial.legendre.leggauss(n)
            u = 0.5 * U * (x + 1)
            ww = 0.5 * U * w * phi_np(u)
            lu = np.log(u)
            f = lambda z: np.sum(ww * np.exp((2 * z - 1) * lu))
            z = zc
            try:
                for _ in range(40):
                    h = 1e-6
                    z = z - f(z) / ((f(z + h) - f(z - h)) / (2 * h))
                naive.append({'U': fs(U), 'nodes': n, 'root_near_chapter_zero': fz(z, 6), 'abs_plain_sum_there': fs(abs(f(z)), 3)})
            except Exception as e:
                naive.append({'U': fs(U), 'nodes': n, 'error': repr(e)})
    R['plain_quadrature_roots_started_at_chapter_zero'] = naive
    log('plain Gauss Legendre roots started at the chapter zero (U, nodes, root):', [(d['U'], d['nodes'], d.get('root_near_chapter_zero')) for d in naive])

    # ---- (g) beyond the chapter: the continuation of M to Re z < 0 (the formula above is the meromorphic continuation, poles at z = 0, -1, -2, ...)
    if not a.quick:
        xl = [-3.95 + 0.5 * i for i in range(9)]
        zl, wl, warnl = scan_cells(Mo, xl, ye, -4.5)
        kl = Mo.M(np.array([complex(x, y) for x in np.arange(-3.9, 0.0, 0.3) for y in np.arange(0.2, 14.0, 0.7)]), cond=True)[1]
        pts_l = np.array([complex(x, y) for x in np.arange(-3.9, 0.0, 0.3) for y in np.arange(0.2, 14.0, 0.7)])
        rel_l = np.abs(Mo.M(pts_l) - Mo2.M(pts_l)) / np.abs(Mo.M(pts_l))
        R['left_half_plane'] = {'region': [fs(xl[0]), fs(xl[-1]), fs(ye[0]), fs(ye[-1])], 'sum_of_cell_windings': wl, 'zeros_found': len(zl), 'warnings': [str(w) for w in warnl],
                                'zeros': [fz(z, 10) for z in zl], 'max_kappa_on_sample': fs(kl.max(), 3), 'max_relative_difference_of_two_layouts_on_sample': fs(rel_l.max(), 3)}
        log('continuation to the left half plane, -3.95 <= Re z <= 0.05, 0.05 <= Im z <= 14.05: total winding', wl, 'zeros found', len(zl), 'warnings', warnl,
            'max kappa', fs(kl.max(), 3), 'two layouts max relative difference', fs(rel_l.max(), 3))
        for z in zl[:40]:
            log('   zero', fz(z, 10))
    R['meta']['seconds'] = fs(time.time() - T0, 4)
    out_path.write_text(json.dumps(R, indent=1))
    log('wrote', out_path)


if __name__ == '__main__':
    main()
