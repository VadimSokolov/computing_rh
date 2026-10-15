# phase_gas.py: the zero gas of the de Bruijn Newman flow, Proposition cp:prop:gas (parts (i) and (ii)) in Section cp:sec:gas
# (The zero gas along the flow) of Chapter 14 (The Thorin phase diagram); the stop pair of part (i) is tied to the figure
# cp:fig:flow, whose drawing is made by phase_collide.py from the JSON written here. Added for the book in October 2026, written
# from the mathematics (the chapter's own script part4.py was read for comparison only, never run).
# Method. In the variable z = 2 gamma the zeros of H_t move, while simple, by dz_k/dt = 2 sum_{w != z_k} 1/(z_k - w), the sum over
# all other zeros w of both signs (Section sec:ch21:flow): the other positive zeros, all their mirror images, and the mirror of
# z_k itself. The first N = 600 zeros (heights up to 939.02) are integrated by an adaptive Runge Kutta method of order 8 (DOP853,
# rtol = atol = 1e-10) backward from t = 0 until the first gap z_{k+1} - z_k falls below 0.02, and forward to t = 0.3. The zeros
# 601 to 1700 of data/zeros_1700.txt stay frozen at their t = 0 positions and act as an explicit tail, the density
# log(z/4 pi)/(4 pi) continues beyond the last of them as a continuum, summed in closed form by a power series in (z/c)^2 (the
# chapter's script keeps only the leading term -4 z sum_{j>N} z_j^{-2}). Recorded at sample times: the minimum normalised gap
# s_min = min_k (z_{k+1} - z_k) log(z_{k+1}/4 pi)/(4 pi), and min_k Delta_k^2 g_k over the first 400 pairs (k = 1..400), with
# g_k = sum over the zeros of both signs other than z_k, z_{k+1} of (z_k - w)^-2 + (z_{k+1} - w)^-2 (the quantity of Csordas, Smith
# and Varga, Section sec:ch21:bounds), once with the tail (zeros 601 to 1700 and the continuum) and once over the 600 integrated
# zeros only (the chapter's convention); both are written. The stop pair, its ordinates and indices (from data/zeros_1700.txt),
# D_0 and the isolated pair time -D_0^2/8 give part (i); the crossing of 4/5 by min Delta^2 g is found by root finding on the
# dense output; the pair (71, 72) is tracked for the right panel of the figure and for the extrapolation to its collision.
# Error control. (a) The velocities at t = 0 of five zeros are compared with the exact H_zz/H_z of H_0 (200 to 360 digit
# trapezoid sum, as in dbn_collisions.py), which tests the formula, the mirror terms and the tail treatment. (b) The ODE
# positions of nine zeros at three times are compared with Newton's method on H_t itself started at the ODE position (the
# truncation error of the 600 zero system, largest near the top of the system where the zeros above are frozen). (c) The stop
# pair's double zero H_t = dH_t/dz = 0 is located by Newton's method in (z, t) on H_t (H_t = -H_zz holds exactly for the
# trapezoid sum), which gives the collision time that the stop only approximates (D^2 = 8 (t - t_*) gives t_* = t_stop - 0.02^2/8).
# (d) The backward leg is repeated at rtol = atol = 1e-12. Nothing here is a proof: the ODE is floating point and the monotonicity
# of (ii) is a statement about the sample times only.
# Output: data/phase_gas.json (values as strings; times, s_min, csv_min, D^2 of the pair (71, 72) at the sample times, the checks).
# Usage: python3 phase_gas.py [--quick] [--out FILE] [--workers W] [--zeros FILE]; --quick integrates 100 zeros (about 25 s) and
# writes to /tmp unless --out is given; the full case (600 zeros) is meant for the Hopper cluster (misc/reviews/nick-ch14/compute-C4-run.sh).
import argparse
import json
import math
import os
import sys
import tempfile
import time
from multiprocessing import get_context
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

T0 = time.time()
GSTOP = 0.02      # the backward run stops at the first gap below this (in z)
TFWD = 0.30
PI4 = 4 * math.pi


def fs(x, n=12):
    return '%.*g' % (n, x)


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def find_zeros(path):
    here = Path(__file__).resolve().parent
    for c in ([Path(path)] if path else []) + [Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), here / 'zeros_1700.txt',
                                                here.parent / 'data' / 'zeros_1700.txt']:
        if c.exists():
            return c
    sys.exit('zeros_1700.txt not found: pass --zeros')


# ---------------------------------------------------------------------------------------------------------------------
# The gas: N integrated zeros, an explicit frozen tail (zeros N+1..J) and a continuum beyond z_J
# ---------------------------------------------------------------------------------------------------------------------
class Cont:
    """The zeros above the last explicit zero z_J as a continuum with the density rho(w) = log(w/4 pi)/(4 pi), starting at the
    point c where the smooth counting function (c/4 pi) log(c/4 pi e) + 7/8 equals J + 1/2; closed forms as power series in (x/c)^2."""

    def __init__(self, J):
        c = 2 * math.pi * J / max(math.log(J), 1.0) * 2.0
        for _ in range(100):
            f = c / PI4 * (math.log(c / PI4) - 1) + 7.0 / 8 - (J + 0.5)
            c -= f / (math.log(c / PI4) / PI4)
        self.c = c
        L = math.log(c / PI4)
        m = np.arange(80)
        # a_m = (L/(2m+1) + 1/(2m+1)^2)/(4 pi c): x^{2m} int_c^oo log(w/4pi)/(4pi) w^{-2m-2} dw = a_m (x/c)^{2m}
        self.am = (L / (2 * m + 1) + 1.0 / (2 * m + 1) ** 2) / (PI4 * c)
        self.m = m

    def drift(self, x):      # int_c^oo rho(w) [1/(x-w) + 1/(x+w)] dw
        u = (x / self.c) ** 2
        return -2.0 * x * ((u[:, None] ** self.m[None, :]) @ self.am)

    def S(self, x):          # int_c^oo rho(w) [(x-w)^-2 + (x+w)^-2] dw
        u = (x / self.c) ** 2
        return 2.0 * ((u[:, None] ** self.m[None, :]) @ ((2 * self.m + 1) * self.am))


class Gas:
    def __init__(self, gam, N, J):
        z = 2.0 * gam
        self.N, self.J = N, J
        self.z0 = z[:N].copy()
        self.w = z[N:J].copy()
        self.w2 = self.w ** 2
        self.cont = Cont(J)
        self.c = self.cont.c

    def cont_drift(self, x):
        return self.cont.drift(x)

    def cont_S(self, x):
        return self.cont.S(x)

    def rhs(self, t, z):
        d = z[:, None] - z[None, :]
        np.fill_diagonal(d, np.inf)
        s = (1.0 / d).sum(1) + (1.0 / (z[:, None] + z[None, :])).sum(1)        # partners z_j, -z_j, and -z_k itself
        s += (2.0 * z[:, None] / (z[:, None] ** 2 - self.w2[None, :])).sum(1)  # frozen explicit tail
        s += self.cont_drift(z)
        return 2.0 * s

    def S_int(self, z):           # sum over the integrated zeros of both signs, w != z_k, of (z_k - w)^-2
        d = z[:, None] - z[None, :]
        np.fill_diagonal(d, np.inf)
        return (1.0 / d ** 2).sum(1) + (1.0 / (z[:, None] + z[None, :]) ** 2).sum(1)

    def S_tail(self, z):          # the same over the frozen tail and the continuum
        return (1.0 / (z[:, None] - self.w[None, :]) ** 2 + 1.0 / (z[:, None] + self.w[None, :]) ** 2).sum(1) + self.cont_S(z)


def diag(gas, z, npairs, nsafe):
    gaps = np.diff(z)
    s = gaps * np.log(z[1:] / PI4) / PI4
    S1 = gas.S_int(z)
    S2 = gas.S_tail(z)
    g_src = S1[:-1] + S1[1:] - 2.0 / gaps ** 2
    g_full = g_src + S2[:-1] + S2[1:]
    q_src = gaps ** 2 * g_src
    q_full = gaps ** 2 * g_full
    ks = int(np.argmin(s))
    kn = int(np.argmin(s[:nsafe]))
    kc = int(np.argmin(q_full[:npairs]))
    kd = int(np.argmin(q_src[:npairs]))
    return dict(s_min=float(s[ks]), s_arg=ks + 1, s_min_safe=float(s[kn]), csv=float(q_full[kc]), csv_arg=kc + 1,
                csv_src=float(q_src[kd]), csv_src_arg=kd + 1, gaps=gaps, q_full=q_full, g_full=g_full, s=s)


# ---------------------------------------------------------------------------------------------------------------------
# The exact layer: H_t(z) = (1/4) int_0^oo e^{t v^2/4} Phi(v) cos(z v/2) dv as a trapezoid sum (step 1/hinv on [0, 3.6]), mpmath
# ---------------------------------------------------------------------------------------------------------------------
_NODES = {}
_WC = {}


def phi(v):
    pi = mp.pi
    e2, a9, a5 = mp.exp(2 * v), mp.exp(9 * v / 2), mp.exp(5 * v / 2)
    eps = mp.mpf(10) ** (-(mp.mp.dps + 5))
    tot, n = mp.mpf(0), 1
    while True:
        term = (2 * pi * pi * n ** 4 * a9 - 3 * pi * n * n * a5) * mp.exp(-pi * n * n * e2)
        tot += term
        if n > 2 and abs(term) < abs(tot) * eps:
            break
        n += 1
    return 2 * tot


def get_nodes(hinv, dps):
    key = (hinv, dps)
    if key not in _NODES:
        with mp.workdps(dps + 20):
            h = mp.mpf(1) / hinv
            V = [j * h for j in range(int(round(3.6 * hinv)) + 1)]
            _NODES[key] = (h, V, [phi(v) for v in V], [v / 2 for v in V])
    return _NODES[key]


def get_w(t, hinv, dps):
    key = (hinv, dps)
    if _WC.get('key') != key or _WC.get('t') != t:
        h, V, P, U = get_nodes(hinv, dps)
        with mp.workdps(dps):
            tt = mp.mpf(t)
            W = [h * P[j] * mp.exp(tt * V[j] * V[j] / 4) / 4 for j in range(len(V))]
            W[0] = W[0] / 2
        _WC.update(key=key, t=t, W=W)
    return _WC['W']


def parts(z, t, hinv, dps):
    """H, H_z, H_zz, H_zzz of the trapezoid sum at real z and time t."""
    h, V, P, U = get_nodes(hinv, dps)
    W = get_w(t, hinv, dps)
    with mp.workdps(dps):
        x = mp.mpf(z) * h / 2
        c1, s1 = mp.cos(x), mp.sin(x)
        c, s = mp.mpf(1), mp.mpf(0)
        H = Hz = Hzz = Hzzz = mp.mpf(0)
        for j in range(len(W)):
            u = U[j]
            wc, ws = W[j] * c, W[j] * s
            H += wc
            Hz -= u * ws
            Hzz -= u * u * wc
            Hzzz += u * u * u * ws
            c, s = c * c1 - s * s1, s * c1 + c * s1
    return H, Hz, Hzz, Hzzz


def cfg(z):
    """(1/step, digits) so that the trapezoid error and the cancellation e^{-pi z/8} are covered; the digits are rounded up to
    a multiple of 40 so that the node tables are shared between tasks."""
    return (200 if z <= 800 else 400), 40 * int(math.ceil((0.1706 * z + 45) / 40))


def newton_zero(z, t, hinv, dps):
    with mp.workdps(dps):
        z = mp.mpf(z)
        for _ in range(40):
            H, Hz, Hzz, _ = parts(z, t, hinv, dps)
            dz = H / Hz
            z -= dz
            if abs(dz) < mp.mpf(10) ** -35:
                return z, Hz, Hzz
    return None, None, None


def double_zero(z, t, hinv, dps):
    with mp.workdps(dps):
        z, t = mp.mpf(z), mp.mpf(t)
        for _ in range(60):
            H, Hz, Hzz, Hzzz = parts(z, t, hinv, dps)
            det = -Hz * Hzzz + Hzz * Hzz
            dz = (-Hzzz * H + Hzz * Hz) / det
            dt = (-Hzz * H + Hz * Hz) / det
            z -= dz
            t -= dt
            if abs(dz) < mp.mpf(10) ** -30 and abs(dt) < mp.mpf(10) ** -30:
                return z, t
    return None, None


def task(a):
    kind = a[0]
    try:
        if kind == 'vel':          # exact velocity H_zz/H_z at the zero 2 gamma_k of H_0
            k, zk = a[1], a[2]
            hinv, dps = cfg(float(zk))
            with mp.workdps(dps):
                z, Hz, Hzz = newton_zero(mp.mpf(zk), mp.mpf(0), hinv, dps)
                return dict(kind=kind, k=k, v=mp.nstr(Hzz / Hz, 15), dz0=mp.nstr(z - mp.mpf(zk), 3))
        if kind == 'track':        # Newton on H_t from the ODE position
            k, t, zode = a[1], a[2], a[3]
            hinv, dps = cfg(zode)
            with mp.workdps(dps):
                z, Hz, Hzz = newton_zero(mp.mpf(zode), mp.mpf(t), hinv, dps)
                return dict(kind=kind, k=k, t=t, zode=zode, zex=mp.nstr(z, 18), diff=float(z - zode), dps=dps, hinv=hinv)
        if kind == 'coll':         # double zero near (z, t)
            z, t = a[1], a[2]
            hinv, dps = cfg(z)
            zc, tc = double_zero(z, t, hinv, dps)
            return dict(kind=kind, z=mp.nstr(zc, 20), t=mp.nstr(tc, 25), dps=dps, hinv=hinv)
    except Exception as e:      # pragma: no cover
        return dict(kind=kind, error=repr(e))
    return dict(kind=kind, error='unknown')


def run_tasks(tasks, workers):
    if workers > 1:
        try:
            with get_context('fork').Pool(workers) as pool:
                return pool.map(task, tasks, chunksize=1)
        except Exception as e:      # fall back to one process rather than lose the run
            log('the worker pool failed (%r): running the exact tasks in one process' % (e,))
    return [task(a) for a in tasks]


# ---------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__ if False else 'zero gas along the de Bruijn Newman flow (cp:prop:gas)')
    ap.add_argument('--quick', action='store_true', help='100 zeros, coarse sampling, writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/phase_gas.json; quick: a /tmp file)')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    ap.add_argument('--n', type=int, default=None, help='number of integrated zeros (default 600; quick 100)')
    a = ap.parse_args()
    quick = a.quick
    N = a.n or (100 if quick else 600)
    J = 400 if quick else 1700
    out_path = a.out or (os.path.join(tempfile.gettempdir(), 'phase_gas_quick.json') if quick else 'data/phase_gas.json')
    npairs = min(400, N - 2) if not quick else 80
    nsafe = int(0.85 * (N - 1))
    zf = find_zeros(a.zeros)
    lines = [ln.strip() for ln in zf.read_text().split() if ln.strip()]
    gam = np.array([float(x) for x in lines[:J]])
    gas = Gas(gam, N, J)
    z0 = gas.z0
    R = dict(script='phase_gas.py', mode='quick' if quick else 'full')
    R['parameters'] = dict(N=N, J=J, rtol='1e-10', atol='1e-10', method='DOP853', gap_stop_z=fs(GSTOP), t_forward=fs(TFWD),
                           csv_pairs=npairs, first_height_theta=lines[N - 1][:12], continuum_start_z=fs(gas.c, 10),
                           zeros_file=str(zf), safe_gaps_for_s_min_safe=nsafe)
    log('zeros', N, J, 'top theta', lines[N - 1][:10], 'continuum from z =', fs(gas.c, 8))

    # (0) the right hand side at t = 0 is tiny for the true zeros (the zeros are nearly in equilibrium)
    v0 = gas.rhs(0.0, z0)
    R['velocity_scale_t0'] = dict(max_abs_velocity_all=fs(np.abs(v0).max(), 5), max_abs_velocity_first_85_percent=fs(np.abs(v0[:nsafe]).max(), 5),
                                  note='the closest pairs move at about 2/D in z, so this is not small')
    log('max |dz/dt| at t=0 over the first 85 percent of the zeros', fs(np.abs(v0[:nsafe]).max(), 5))

    # (1) the integration
    def gap_event(t, z):
        return float(np.min(np.diff(z))) - GSTOP
    gap_event.terminal = True
    gap_event.direction = -1
    nfev = {}
    solb = solve_ivp(gas.rhs, (0.0, -0.4), z0, method='DOP853', rtol=1e-10, atol=1e-10, max_step=0.002, events=gap_event,
                     dense_output=True)
    nfev['backward'] = int(solb.nfev)
    if solb.status != 1:
        sys.exit('the backward run did not stop at a gap of %g (status %s)' % (GSTOP, solb.status))
    tstop = float(solb.t_events[0][0])
    zstop = solb.y_events[0][0]
    solf = solve_ivp(gas.rhs, (0.0, TFWD), z0, method='DOP853', rtol=1e-10, atol=1e-10, max_step=0.01, dense_output=True)
    nfev['forward'] = int(solf.nfev)
    log('backward stop at t =', fs(tstop, 10), 'nfev', nfev, 'forward status', solf.status)
    R['parameters']['nfev'] = nfev

    # (2) part (i): the stop pair
    gstop = np.diff(zstop)
    ks = int(np.argmin(gstop))
    order = np.argsort(gstop)
    k1 = ks + 1                       # 1-based index of the lower zero of the pair
    mp.mp.dps = 40
    g1, g2 = mp.mpf(lines[k1 - 1]), mp.mpf(lines[k1])
    D0 = 2 * (g2 - g1)
    iso = -D0 ** 2 / 8
    tcoll_est = tstop - float(gstop[ks]) ** 2 / 8
    R['part_i'] = dict(t_stop=fs(tstop, 10), stop_gap_z=fs(float(gstop[ks]), 8), pair_indices=[k1, k1 + 1],
                       ordinates_theta=[lines[k1 - 1][:14], lines[k1][:14]], D0_z=mp.nstr(D0, 12), isolated_estimate=mp.nstr(iso, 10),
                       collision_by_pair_law_from_stop=fs(tcoll_est, 10),
                       second_smallest_gap_at_stop=dict(pair=[int(order[1]) + 1, int(order[1]) + 2], gap=fs(float(gstop[order[1]]), 6)),
                       ratio_stop_over_isolated=fs(tstop / float(iso), 8))
    log('part (i): pair', k1, k1 + 1, 'theta', lines[k1 - 1][:10], lines[k1][:10], 'D0', mp.nstr(D0, 8), 'iso', mp.nstr(iso, 8))

    # (3) sampling: nb points on [t_stop, 0), nf+1 points on [0, 0.3]
    nb, nf = (8, 12) if quick else (50, 120)
    tb = np.linspace(tstop, 0.0, nb + 1)[:-1]
    tf = np.linspace(0.0, TFWD, nf + 1)
    ts = np.concatenate([tb, tf])
    states = [solb.sol(min(max(t, tstop), 0.0)) for t in tb] + [solf.sol(min(max(t, 0.0), TFWD)) for t in tf]
    rows = [diag(gas, zz, npairs, nsafe) for zz in states]
    i71 = 70
    D2 = [float((zz[i71 + 1] - zz[i71]) ** 2) for zz in states]
    smin = np.array([r['s_min'] for r in rows])
    ssafe = np.array([r['s_min_safe'] for r in rows])
    csv = np.array([r['csv'] for r in rows])
    csvs = np.array([r['csv_src'] for r in rows])
    mono = lambda x: bool(np.all(np.diff(x) > 0))
    R['samples'] = dict(t=[fs(x, 10) for x in ts], s_min=[fs(x, 10) for x in smin], s_min_arg=[r['s_arg'] for r in rows],
                        s_min_safe=[fs(x, 10) for x in ssafe], csv_min=[fs(x, 10) for x in csv], csv_min_arg=[r['csv_arg'] for r in rows],
                        csv_min_source_convention=[fs(x, 10) for x in csvs], D2_pair_71_72=[fs(x, 12) for x in D2])
    R['part_ii'] = dict(interval=[fs(tstop, 8), fs(TFWD)], n_samples=len(ts), s_min_first=fs(smin[0], 6), s_min_last=fs(smin[-1], 6),
                        csv_min_first=fs(csv[0], 6), csv_min_last=fs(csv[-1], 6), csv_min_first_source_convention=fs(csvs[0], 6),
                        csv_min_last_source_convention=fs(csvs[-1], 6), s_min_strictly_increasing=mono(smin),
                        s_min_safe_strictly_increasing=mono(ssafe), csv_min_strictly_increasing=mono(csv),
                        csv_min_source_convention_strictly_increasing=mono(csvs), min_increment_s_min=fs(np.diff(smin).min(), 4),
                        min_increment_csv_min=fs(np.diff(csv).min(), 4))
    r0 = rows[nb]                        # t = 0
    kk = r0['csv_arg'] - 1
    R['t0'] = dict(s_min=fs(r0['s_min'], 8), s_min_pair=[r0['s_arg'], r0['s_arg'] + 1],
                   csv_minimiser_pair=[r0['csv_arg'], r0['csv_arg'] + 1],
                   csv_minimiser_theta=[fs(z0[kk] / 2, 9), fs(z0[kk + 1] / 2, 9)], csv_min=fs(r0['csv'], 8),
                   csv_minimiser_pair_source_convention=[r0['csv_src_arg'], r0['csv_src_arg'] + 1],
                   csv_min_source_convention=fs(r0['csv_src'], 8))
    # the pair (71, 72) from all 1700 zeros of the file, independent of N: the book's Delta^2 g = 0.708, g = 0.33734, Lambda_71 = -0.3048
    za = 2.0 * np.array([float(x) for x in lines[:1700]])
    a_, b_ = za[70], za[71]
    ww = np.concatenate([-za, za])
    keep = np.ones(len(ww), bool)
    keep[[1700 + 70, 1700 + 71]] = False
    ww = ww[keep]
    g71_1700 = float(np.sum((a_ - ww) ** -2.0 + (b_ - ww) ** -2.0))
    c1700 = Cont(1700)
    g71 = g71_1700 + float(c1700.S(np.array([a_, b_])).sum())
    D71 = float(b_ - a_)
    lam = lambda g: ((1 - 1.25 * D71 ** 2 * g) ** 0.8 - 1) / (8 * g)
    R['book_check_pair_71_72'] = dict(Delta=fs(D71, 8), g_first_1700_zeros=fs(g71_1700, 8), g_with_continuum=fs(g71, 8),
                                      Delta2_g=fs(D71 ** 2 * g71_1700, 6), Lambda_71=fs(lam(g71_1700), 6),
                                      Lambda_71_with_continuum=fs(lam(g71), 6),
                                      book_values=dict(Delta='1.44863', g='0.33734', Delta2_g='0.708', Lambda_71='-0.3048'))
    # the crossing of 4/5 by min Delta^2 g
    f = lambda t: diag(gas, solf.sol(t), npairs, nsafe)['csv'] - 0.8
    fsrc = lambda t: diag(gas, solf.sol(t), npairs, nsafe)['csv_src'] - 0.8
    cross = {}
    for name, fn in (('with_tail', f), ('source_convention', fsrc)):
        cross[name] = fs(brentq(fn, 0.0, TFWD, xtol=1e-9), 8) if fn(0.0) * fn(TFWD) < 0 else None
    R['crossing_4_5'] = cross
    # the extrapolation of the chapter: linear fit of D^2 on the first 8 of 30 equally spaced points on [t_stop, 0]
    tt = np.linspace(tstop, 0.0, 30)[:8]
    dd = np.array([(solb.sol(t)[i71 + 1] - solb.sol(t)[i71]) ** 2 for t in tt])
    p = np.polyfit(tt, dd, 1)
    R['pair_71_72'] = dict(ordinates_theta_t0=[lines[70][:14], lines[71][:14]], D0_z=fs(D71, 10), isolated_estimate=fs(-D71 ** 2 / 8, 8),
                           extrapolated_collision_first8of30=fs(-p[1] / p[0], 8),
                           note='D^2 along the 600 zero ODE only on [t_stop, 0.3]; the exact pair is tracked in phase_collide.py')
    log('part (ii): s_min', fs(smin[0], 5), '->', fs(smin[-1], 5), 'csv', fs(csv[0], 5), '->', fs(csv[-1], 5), 'monotone', mono(smin), mono(csv),
        'cross', cross)

    # (4) the rtol check of the backward leg
    sol12 = solve_ivp(gas.rhs, (0.0, -0.4), z0, method='DOP853', rtol=1e-12, atol=1e-12, max_step=0.002, events=gap_event, dense_output=False)
    R['rtol_check'] = dict(t_stop_rtol_1e_12=fs(float(sol12.t_events[0][0]), 10), difference=fs(float(sol12.t_events[0][0]) - tstop, 3),
                           max_abs_diff_of_stop_state=fs(float(np.abs(sol12.y_events[0][0] - zstop).max()), 3))
    log('rtol check', R['rtol_check'])

    # (5) the exact layer: velocities, positions, the stop pair's double zero
    ks_vel = [1, 71, N // 2, int(0.9 * N)] if not quick else [1, 71, 90]
    kt = [10, 71, 72, 150, 300, k1, k1 + 1, 560, 590] if not quick else [10, 40, 71, 72, 95]
    kt = sorted(set(k for k in kt if k <= N))
    tchk = [0.8 * tstop, 0.1, 0.3]
    tasks = [('vel', k, mp.nstr(2 * mp.mpf(lines[k - 1]), 60)) for k in ks_vel]
    tasks += [('track', k, t, float((solb if t < 0 else solf).sol(t)[k - 1])) for t in tchk for k in kt]
    tasks += [('coll', 0.5 * float(zstop[ks] + zstop[ks + 1]), tcoll_est)]
    log('exact layer: %d tasks, %d workers' % (len(tasks), a.workers))
    res_all = run_tasks(tasks, max(1, a.workers))
    R['exact_layer_errors'] = [r for r in res_all if 'error' in r]
    res = [r for r in res_all if 'error' not in r]
    vel = []
    for r in res:
        if r['kind'] == 'vel':
            k = r['k']
            vel.append(dict(k=k, v_exact=r['v'], v_ode=fs(v0[k - 1], 12), abs_diff=fs(abs(float(r['v']) - v0[k - 1]), 3)))
    R['velocity_check_t0'] = vel
    trk = [dict(k=r['k'], t=fs(r['t'], 6), z_ode=fs(r['zode'], 12), z_exact=r['zex'], abs_diff=fs(abs(r['diff']), 3)) for r in res if r['kind'] == 'track']
    R['truncation_check'] = dict(rows=trk, max_abs_diff_first_85_percent=fs(max([abs(float(x['abs_diff'])) for x in trk if x['k'] <= nsafe] or [0]), 3),
                                 max_abs_diff_all=fs(max([abs(float(x['abs_diff'])) for x in trk] or [0]), 3),
                                 note='Newton on H_t from the ODE position; a difference far below the gap (about 3 in z) means the same zero')
    gaperr = []
    for r1 in [r for r in res if r['kind'] == 'track']:
        for r2 in [r for r in res if r['kind'] == 'track']:
            if r2['k'] == r1['k'] + 1 and r2['t'] == r1['t']:
                gaperr.append(dict(pair=[r1['k'], r2['k']], t=fs(r1['t'], 6), gap_ode=fs(r2['zode'] - r1['zode'], 10),
                                   gap_exact=fs(float(mp.mpf(r2['zex']) - mp.mpf(r1['zex'])), 10),
                                   abs_diff=fs(abs((r2['zode'] - r1['zode']) - float(mp.mpf(r2['zex']) - mp.mpf(r1['zex']))), 3)))
    R['truncation_check']['gap_errors_of_adjacent_pairs'] = gaperr
    cl = [r for r in res if r['kind'] == 'coll']
    if not cl:
        R['part_i']['exact_collision'] = dict(error='the double zero Newton step failed', details=R['exact_layer_errors'])
        R['seconds'] = round(time.time() - T0, 1)
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        json.dump(R, open(out_path, 'w'), indent=1)
        sys.exit('exact collision failed; partial JSON written')
    cres = cl[0]
    tex = float(mp.mpf(cres['t']))
    R['part_i']['exact_collision'] = dict(z=cres['z'], t=cres['t'], dps=cres['dps'], step_inverse=cres['hinv'],
                                          t_minus_stop=fs(tex - tstop, 4), t_minus_pair_law_from_stop=fs(tex - tcoll_est, 4),
                                          ratio_to_isolated=fs(tex / float(iso), 8),
                                          note='double zero H_t = H_t_z = 0 by Newton in (z, t) on the trapezoid sum; not certified')
    log('exact collision of the stop pair: t =', cres['t'][:14], ' pair law from stop', fs(tcoll_est, 10))
    R['seconds'] = round(time.time() - T0, 1)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    json.dump(R, open(out_path, 'w'), indent=1)
    log('wrote', out_path)


if __name__ == '__main__':
    main()
