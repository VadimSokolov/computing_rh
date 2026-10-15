# phase_collide.py: the collision of the pair (gamma_71, gamma_72) of zeros of H_t under the backward de Bruijn Newman flow, part (iii)
# of Proposition cp:prop:gas in Section cp:sec:gas of Chapter 14 (The Thorin phase diagram), and the drawing of Figure cp:fig:flow
# (fig/cp_flow.pdf). Added for the book in October 2026, written from the mathematics; the book's own value is t_* = -0.2880 with a
# grid error 6.1e-5 (Section sec:ch21:flow, and data/dbn_collisions.json: -0.287979837009 from a double zero Newton step at 200
# digits), which this script checks, certifies at one end and tightens.
# H_t(z) = (1/4) int_0^oo e^{t v^2/4} Phi(v) cos(z v/2) dv, z = 2 theta, evaluated by the trapezoid sum with step 1/200 on [0, 3.6]
# in Arb (python-flint) at 640 bits. Phi is summed in ball arithmetic with a rigorous tail bound; the sum S_t(z) is an entire
# function that satisfies the backward heat equation exactly (S_t = -S_zz), and H_t(z) - S_t(z) is enclosed by a rigorous bound:
# the trapezoid theorem |I_h - I| <= 2M/(e^{2 pi a/h} - 1) (Trefethen and Weideman, SIAM Review 2014, Theorem 5.1) on the strip
# |Im v| < a = 1/2 with M bounded in closed form (the cos factor by e^{z a/2}, Phi by its theta series on the strip), plus the
# cut-off at v = 3.6. At theta = 185 this bound is below 1e-200 against H_t of size 1e-63.
# (a) Certified end. Signs of H_t(2 theta) on the grid theta = 184.6 + 0.01 k, k = 0..130 (131 points), each decided by its ball;
# bisection on t of the number of sign changes (0 at t = -0.30, 2 at t = -0.27) gives t_lo < t_hi with 0 and 2 certified changes:
# two distinct real zeros at t_hi, so the pair is still real there and t_* < t_hi rigorously. At t_lo no change on a grid of step
# 0.01 does not exclude two real zeros less than 0.01 apart, so t_lo is NOT a certified lower end.
# (b) The other end. The pair law D^2 = 8 (t - t_*) (from H ~ c [(z - z_*)^2 - 2 (t - t_*)] near a double zero and H_t = -H_zz) says
# the changes vanish at most h^2/2 = 5e-5 above t_*; that is an asymptotic argument. Here it is replaced by a double zero of S_t:
# Newton's method in (z, t) on (S, S_z) = 0, then a Krawczyk test in Arb on the box of half-width r = 2^-400 (7.8e-121) around the result
# (the Jacobian [[S_z, -S_zz], [S_zz, -S_zzz]] has determinant S_zz^2 > 0 there), which proves that S has a real double zero in the box:
# t_* is located to +-8e-121 for the sum S, which differs from H by the enclosed bound. What is not proved: that this double zero is
# the collision of the zeros 71 and 72 (tracked here by Newton continuation from t = 0, not enclosed), and the transfer from S to H
# of the Krawczyk statement (the difference H - S is enclosed for values only, not for derivatives).
# (c) The exact pair is tracked from t_* to 0.3 (Newton continuation on S, 640 bits) for the right panel of the figure: D^2 against
# the isolated pair law D_0^2 + 8 t and the pair law 8 (t - t_*); the extrapolation of the integrated D^2 of phase_gas.py (first 8 of
# 30 points, as the chapter does) is compared with t_*.
# Output: data/phase_collide.json (values as strings), fig/cp_flow.pdf (left panel from data/phase_gas.json, copied into the JSON;
# grayscale safe: line styles and markers). Usage: python3 phase_collide.py [--quick] [--out FILE] [--fig FILE] [--gas FILE]
# [--book FILE] [--zeros FILE] [--workers W] [--figure]; --figure redraws from the JSON without computing; --quick uses 10 bisection steps and
# few track points (about 20 s) and writes to /tmp unless --out and --fig are given.
import argparse
import json
import math
import os
import sys
import tempfile
import time
from multiprocessing import get_context
from pathlib import Path

try:
    from flint import arb, ctx, fmpq
except ImportError:      # --figure needs only matplotlib
    arb = ctx = fmpq = None

T0 = time.time()
PREC = 640
HINV = 200
NJ = 720                      # nodes v_j = j/200, j = 0..720 (v up to 3.6)
if ctx is not None:
    ctx.prec = PREC
_T = {}


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def fl(x):
    """float of an arb (midpoint)."""
    return float(x.mid().str(20, radius=False))


def find_file(path, name):
    here = Path(__file__).resolve().parent
    for c in ([Path(path)] if path else []) + [Path('data') / name, Path(name), here / name, here.parent / 'data' / name]:
        if c.exists():
            return c
    return None


# ---------------------------------------------------------------------------------------------------------------------
# Arb layer
# ---------------------------------------------------------------------------------------------------------------------
def phi_ball(v):
    """Phi(v) = 2 sum_n (2 pi^2 n^4 e^{9v/2} - 3 pi n^2 e^{5v/2}) exp(-pi n^2 e^{2v}) for v >= 0, with a rigorous tail bound:
    |term_m| <= 4 pi^2 m^4 e^{9v/2} exp(-pi m^2 e^{2v}) =: B_m, and B_{m+1}/B_m <= 16 e^{-3 pi} < 1/2, so the sum over m >= n is <= 2 B_n."""
    pi = arb.pi()
    e2, a9, a5 = (2 * v).exp(), (9 * v / 2).exp(), (5 * v / 2).exp()
    thr = arb(2) ** (-(PREC + 60))
    tot, n = arb(0), 1
    while True:
        tot += (2 * pi * pi * n ** 4 * a9 - 3 * pi * n * n * a5) * (-pi * n * n * e2).exp()
        n += 1
        B = 4 * pi * pi * n ** 4 * a9 * (-pi * n * n * e2).exp()
        if B < thr:
            break
    return 2 * (tot + arb(0, (2 * B).upper()))


def tables():
    if not _T:
        ctx.prec = PREC
        v = [arb(j) / HINV for j in range(NJ + 1)]
        _T.update(h=arb(1) / HINV, v=v, u=[x / 2 for x in v], v2=[x * x / 4 for x in v], P=[phi_ball(x) for x in v])
    return _T


def weights(t):
    T = tables()
    W = [T['h'] * T['P'][j] * (t * T['v2'][j]).exp() / 4 for j in range(NJ + 1)]
    W[0] = W[0] / 2
    return W


def sums(z, W, order):
    """S, S_z, S_zz, S_zzz of the trapezoid sum (balls) at z for the weights W; order = highest derivative."""
    u = tables()['u']
    H = Hz = Hzz = Hzzz = arb(0)
    if order == 0:
        for j in range(NJ + 1):
            H += W[j] * (z * u[j]).cos()
        return H, None, None, None
    for j in range(NJ + 1):
        s, c = (z * u[j]).sin_cos()
        wc, ws = W[j] * c, W[j] * s
        H += wc
        Hz -= u[j] * ws
        if order >= 2:
            Hzz -= u[j] * u[j] * wc
        if order >= 3:
            Hzzz += u[j] * u[j] * u[j] * ws
    return H, Hz, Hzz, Hzzz


def err_bound(z, t):
    """Enclosure of |H_t(z) - S_t(z)| for real z and t <= 1 (see the header): E_quad + E_cutoff."""
    T = tables()
    a = arb(1) / 2
    pi = arb.pi()
    cc = (2 * a).cos()
    tu = t.upper()
    tp = tu if tu > 0 else arb(0)
    tab = abs(t).upper()
    Ssum = arb(0)
    for n in range(1, 15):
        an = pi * n * n * cc - tp / 8
        Ssum += 4 * pi * pi * n ** 4 * (-an).exp() / (an - arb(5) / 4)
    n = 15
    an = pi * n * n * cc - tp / 8
    Ssum += 2 * 4 * pi * pi * n ** 4 * (-an).exp() / (an - arb(5) / 4)
    M = (z.upper() * a / 2 + tab * a * a / 4).exp() * Ssum / 2
    Equad = M / ((2 * pi * a * HINV).exp() - 1)
    V = arb(36) / 10
    gV = (tp * V * V / 4).exp() * 16 * pi * pi * (9 * V / 2).exp() * (-pi * (2 * V).exp()).exp() / 4
    return Equad + T['h'] * gV


def enclosure(z, W, t):
    S = sums(z, W, 0)[0]
    return S + arb(0, err_bound(z, t).upper())


def sign_of(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


def grid_task(args):
    tp, tq, thetas = args
    t = arb(tp) / tq
    W = weights(t)
    out = []
    for (p, q) in thetas:
        val = enclosure(2 * arb(p) / q, W, t)
        out.append((sign_of(val), val.str(6)))
    return out


class Grid:
    def __init__(self, thetas, workers):
        self.thetas = thetas
        self.workers = workers
        tables()
        try:
            self.pool = get_context('fork').Pool(workers) if workers > 1 else None
        except Exception as e:      # fall back to one process rather than lose the run
            log('the worker pool failed (%r): one process' % (e,))
            self.pool = None

    def signs(self, tq):
        tp, tqq = int(tq.p), int(tq.q)
        if self.pool is None:
            res = grid_task((tp, tqq, self.thetas))
        else:
            n = len(self.thetas)
            chunks = [self.thetas[i::self.workers] for i in range(self.workers)]
            try:
                parts = self.pool.map(grid_task, [(tp, tqq, c) for c in chunks])
            except Exception as e:
                log('the worker pool failed (%r): one process' % (e,))
                self.pool = None
                return self.signs(tq)
            res = [None] * n
            for i in range(self.workers):
                for k, r in enumerate(parts[i]):
                    res[i + k * self.workers] = r
        return [r[0] for r in res], [r[1] for r in res]

    def count(self, tq):
        s, vals = self.signs(tq)
        if 0 in s:
            return None, s, vals
        return sum(1 for x, y in zip(s, s[1:]) if x != y), s, vals


def newton_double(z, t, iters=14):
    for _ in range(iters):
        H, Hz, Hzz, Hzzz = sums(z, weights(t), 3)
        det = -Hz * Hzzz + Hzz * Hzz
        dz = (-Hzzz * H + Hzz * Hz) / det
        dt = (-Hzz * H + Hz * Hz) / det
        z, t = (z - dz).mid(), (t - dt).mid()
    return z, t, dz, dt


def krawczyk(z0, t0, r):
    """Krawczyk test for F = (S, S_z) on the box x0 + [-r, r]^2; returns (ok, contraction bound, det enclosure, Szz enclosure)."""
    H, Hz, Hzz, Hzzz = sums(z0, weights(t0), 3)
    a, b, c, d = Hz, -Hzz, Hzz, -Hzzz
    det = a * d - b * c
    C = [[(d / det).mid(), (-b / det).mid()], [(-c / det).mid(), (a / det).mid()]]
    CF = [C[0][0] * H + C[0][1] * Hz, C[1][0] * H + C[1][1] * Hz]
    rb = arb(0, r)
    Zb, Tb = z0 + rb, t0 + rb
    HB, HzB, HzzB, HzzzB = sums(Zb, weights(Tb), 3)
    J = [[HzB, -HzzB], [HzzB, -HzzzB]]
    M = [[(1 if i == j else 0) - (C[i][0] * J[0][j] + C[i][1] * J[1][j]) for j in range(2)] for i in range(2)]
    disp = [-CF[i] + M[i][0] * rb + M[i][1] * rb for i in range(2)]
    bound = [abs(x).upper() for x in disp]
    ok = all(bd < arb(r) for bd in bound)
    return ok, [bd / arb(r) for bd in bound], det, HzzB, M


def newton_zero(z, t, W=None, iters=40):
    W = W if W is not None else weights(t)
    for _ in range(iters):
        H, Hz, Hzz, _ = sums(z, W, 2)
        dz = H / Hz
        z = (z - dz).mid()
        if abs(dz) < arb(10) ** -60:
            return z, Hz, Hzz
    return None, None, None


def track(zc, tc, ts, direction):
    """Follow one root of S_t from the double zero (zc, tc) for the increasing times ts (> tc); direction = +1 or -1."""
    out, zprev, tprev, vprev = [], None, tc, None
    for t in ts:
        s = t - tc
        if zprev is None:
            z = zc + direction * (2 * s).sqrt()
        else:
            z = zprev + (t - tprev) * vprev
        W = weights(t)
        zz, Hz, Hzz = newton_zero(z.mid(), t, W)
        if zz is None:
            sub = 2
            while zz is None and sub <= 256:       # retry with substeps from the previous solution
                zz = zprev
                tt = tprev
                for i in range(1, sub + 1):
                    ti = tprev + (t - tprev) * i / sub
                    Hh = sums(zz, weights(ti), 2)
                    zg = zz + (ti - tt) * Hh[2] / Hh[1]
                    zz, Hz, Hzz = newton_zero(zg.mid(), ti)
                    tt = ti
                    if zz is None:
                        break
                sub *= 2
            if zz is None:
                raise RuntimeError('tracking failed at t = %s' % t.str(10))
        out.append(zz)
        zprev, tprev, vprev = zz, t, (Hzz / Hz).mid()
    return out


# ---------------------------------------------------------------------------------------------------------------------
def draw(R, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    G = R.get('figure_data', {}).get('gas')
    if G:
        t = np.array([float(x) for x in G['t']])
        sm = np.array([float(x) for x in G['s_min']])
        cs = np.array([float(x) for x in G['csv_min']])
        k = max(1, len(t) // 14)
        a1.plot(t, sm, color='black', ls='-', lw=1.4, marker='o', ms=3.2, markevery=k, label=r'$s_{\min}(t)$')
        a1.plot(t, cs, color='0.35', ls='--', lw=1.4, marker='s', ms=3.2, markevery=k, mfc='white', label=r'$\min_k\Delta_k^2 g_k$')
        a1.axhline(0.8, color='0.5', ls=':', lw=1.2, label='Lehmer threshold 4/5')
        if G.get('crossing'):
            a1.axvline(float(G['crossing']), color='0.6', ls='-.', lw=0.7)
        a1.set_xlabel('$t$')
        a1.legend(frameon=False, fontsize=7, loc='upper left')
        a1.set_title('zero gas, 600 zeros' if G.get('N') == 600 else 'zero gas, %s zeros' % G.get('N'), fontsize=8)
    Tk = R['track']
    tt = np.array([float(x) for x in Tk['t']])
    d2 = np.array([float(x) for x in Tk['D2']])
    ts = float(R['double_zero']['t_star'])
    D0sq = float(R['pair']['D0_squared'])
    a2.plot(tt, d2, color='black', ls='-', lw=1.5, label=r'$D^2$, exact pair')
    tgrid = np.linspace(-D0sq / 8, tt.max(), 100)
    a2.plot(tgrid, D0sq + 8 * tgrid, color='0.4', ls='--', lw=1.2, label=r'isolated pair, $D_0^2+8t$')
    tp = np.linspace(ts, ts + 0.15, 60)
    a2.plot(tp, 8 * (tp - ts), color='0.55', ls='-.', lw=1.2, label=r'pair law, $8(t-t_*)$')
    og = R.get('figure_data', {}).get('ode_D2')
    if og:
        a2.plot([float(x) for x in og['t']], [float(x) for x in og['D2']], 'o', mfc='white', mec='0.2', ms=3, ls='none', label='ODE, %s zeros' % R.get('figure_data', {}).get('gas', {}).get('N', 600))
    a2.axvline(ts, color='0.6', ls=':', lw=0.8)
    a2.set_xlabel('$t$')
    a2.legend(frameon=False, fontsize=7, loc='upper left')
    a2.set_title('closest pair below 200 (zeros 71, 72)', fontsize=8)
    a2.set_ylim(0, None)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)


def main():
    ap = argparse.ArgumentParser(description='collision of the pair near 185 (cp:prop:gas (iii)) and Figure cp:fig:flow')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--fig', default=None)
    ap.add_argument('--gas', default=None, help='phase_gas.json (default data/phase_gas.json)')
    ap.add_argument('--book', default=None, help='dbn_collisions.json of the book, for the comparison (optional)')
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--figure', action='store_true', help='redraw from --out without computing')
    ap.add_argument('--steps', type=int, default=None, help='bisection steps (default 24; quick 10)')
    a = ap.parse_args()
    tmp = tempfile.gettempdir()
    out = a.out or (os.path.join(tmp, 'phase_collide_quick.json') if a.quick else 'data/phase_collide.json')
    fig = a.fig or (os.path.join(tmp, 'cp_flow_quick.pdf') if a.quick else 'fig/cp_flow.pdf')
    if a.figure:
        draw(json.load(open(out)), fig)
        print('drew', fig)
        return
    if arb is None:
        sys.exit('python-flint is needed for the computation (use the Hopper venv); --figure needs only matplotlib')
    steps = a.steps or (10 if a.quick else 24)
    zf = find_file(a.zeros, 'zeros_1700.txt')
    if zf is None:
        sys.exit('zeros_1700.txt not found: pass --zeros')
    lines = zf.read_text().split()
    R = dict(script='phase_collide.py', mode='quick' if a.quick else 'full')
    R['parameters'] = dict(precision_bits=PREC, step_inverse=HINV, nodes=NJ + 1, v_max='3.6', strip_half_width_a='1/2', bisection_steps=steps,
                           grid='theta = 184.6 + 0.01 k, k = 0..130', t_start=['-0.30', '-0.27'], krawczyk_half_width='2^-400')
    tables()
    log('tables built at', PREC, 'bits; python-flint nodes', NJ + 1)

    # (0) the evaluator at t = 0: the zeros 71, 72 of S_0 against the 70 digit values of the file
    W0 = weights(arb(0))
    chk = []
    for k in (71, 72):
        zk = arb(lines[k - 1]).mid() * 2
        zz, Hz, Hzz = newton_zero(zk, arb(0), W0)
        chk.append((k, zz, (zz - zk).mid()))
    R['evaluator_check_t0'] = dict(zeros=[dict(k=k, z=zz.str(40, radius=False), diff_to_file=abs(df).str(3, radius=False)) for k, zz, df in chk],
                                   note='Newton on S_0 from 2 gamma_k of data/zeros_1700.txt (70 digits): the difference is the rounding of the file')
    z71, z72 = chk[0][1], chk[1][1]
    D0 = z72 - z71
    R['pair'] = dict(indices=[71, 72], z_t0=[z71.str(30, radius=False), z72.str(30, radius=False)], theta_t0=[(z71 / 2).str(14, radius=False), (z72 / 2).str(14, radius=False)],
                     D0=D0.str(15, radius=False), D0_squared=(D0 * D0).str(15, radius=False), isolated_estimate=(-D0 * D0 / 8).str(12, radius=False))
    log('pair at t=0:', R['pair']['theta_t0'], 'D0^2/8 =', (D0 * D0 / 8).str(8, radius=False))

    # (a) the certified sign changes and the bisection
    thetas = [(18460 + k, 100) for k in range(131)]
    G = Grid(thetas, max(1, a.workers))
    lo, hi = fmpq(-30, 100), fmpq(-27, 100)
    c_lo, s_lo, v_lo = G.count(lo)
    c_hi, s_hi, v_hi = G.count(hi)
    log('sign changes at -0.30:', c_lo, ' at -0.27:', c_hi)
    if c_lo != 0 or c_hi != 2:
        sys.exit('unexpected sign change counts at the ends of the starting bracket')
    hist = []
    br14 = None
    for it in range(steps):
        mid = (lo + hi) / 2
        c, s, v = G.count(mid)
        if c is None:
            log('undecided ball at', mid)
            break
        hist.append((float(mid.p) / float(mid.q), c))
        if c >= 2:
            hi = mid
        else:
            lo = mid
        if it == 13:
            br14 = (lo, hi)       # the resolution of the chapter's bisection (14 steps from a bracket of width 0.03)
    c_lo, s_lo, v_lo = G.count(lo)
    c_hi, s_hi, v_hi = G.count(hi)
    tf = lambda q: (arb(int(q.p)) / int(q.q))
    changes = [k for k in range(130) if s_hi[k] != s_hi[k + 1]]
    R['certified'] = dict(
        t_lo=tf(lo).str(12, radius=False), t_hi=tf(hi).str(12, radius=False), width=(tf(hi) - tf(lo)).str(3, radius=False),
        changes_at_t_lo=c_lo, changes_at_t_hi=c_hi, all_131_balls_decided_at_both=bool(0 not in s_lo and 0 not in s_hi),
        sign_change_theta_intervals_at_t_hi=[[str(Fr) for Fr in (fmpq(18460 + k, 100), fmpq(18460 + k + 1, 100))] for k in changes],
        values_at_t_hi_around_the_pair=[dict(theta=str(fmpq(18460 + k, 100)), sign=s_hi[k], enclosure=v_hi[k]) for k in range(max(0, changes[0] - 1), min(131, changes[-1] + 3))] if changes else [],
        statement='at t_hi the enclosures of H_t(2 theta) show %d sign changes on the grid, so two distinct real zeros: t_* < t_hi. At t_lo '
                  'no change on a grid of step 0.01 is NOT a certificate of absence (two zeros closer than 0.01 would escape).' % c_hi,
        bracket_after_14_steps=[(arb(int(q.p)) / int(q.q)).str(12, radius=False) for q in br14] if br14 else None,
        bisection_history=[[('%.9f' % x), c] for x, c in hist])
    log('bracket [t_lo, t_hi] =', R['certified']['t_lo'], R['certified']['t_hi'], 'changes', c_lo, c_hi)

    # (b) the double zero and its Krawczyk certificate
    kc = changes[0]
    zstart = 2 * (arb(18460 + kc) / 100 + arb(1) / 200)
    zs, ts, dz, dt = newton_double(zstart, tf(hi))
    r = arb(2) ** -400
    ok, contr, det, Hzz, M = krawczyk(zs, ts, r)
    R['double_zero'] = dict(z_star=zs.str(60, radius=False), theta_star=(zs / 2).str(40, radius=False), t_star=ts.str(60, radius=False),
                            last_newton_step=[abs(dz).str(3, radius=False), abs(dt).str(3, radius=False)],
                            krawczyk=dict(half_width='2^-400 = 7.8e-121', contained_in_interior=bool(ok), contraction_bounds=[c_.str(3, radius=False) for c_ in contr],
                                          jacobian_determinant_enclosure=det.str(5), S_zz_enclosure=Hzz.str(8),
                                          statement='the trapezoid sum S has a real double zero (S = S_z = 0) in the box z_star +- r, t_star +- r, unique there'
                                          if ok else 'Krawczyk test FAILED'))
    ts_f = fl(ts)
    R['certified']['t_lo_below_t_star'] = bool(tf(lo) < ts)
    R['certified']['t_lo_below_t_star_note'] = ('False means the pair is real at t_lo but no grid point falls between its zeros: a deep grid bisection '
                                                 'converges to the visibility threshold t_* + 2 delta^2, not to t_*, so [t_lo, t_hi] is then not a bracket of the collision')
    # a local certificate just above the double zero: at the rational t_c = t_* + (1 to 2) e-12 the enclosures of H at z_* - e, z_*, z_* + e
    n12 = int(math.floor(ts_f * 1e12)) + 2
    tc = fmpq(n12, 10 ** 12)
    tca = arb(n12) / 10 ** 12
    eps = arb(1) / 250000
    Wc = weights(tca)
    vc = [enclosure(zs + k * eps, Wc, tca) for k in (-1, 0, 1)]
    sc = [sign_of(v) for v in vc]
    R['local_certificate'] = dict(t_c=str(tc), t_c_decimal=tca.str(15, radius=False), t_c_minus_t_star=(tca - ts).str(4, radius=False),
                                  offsets_z=['-4e-6', '0', '4e-6'], signs=sc, enclosures=[v.str(6) for v in vc],
                                  statement='two distinct real zeros of H_{t_c} lie in (z_* - 4e-6, z_*) and (z_*, z_* + 4e-6) (signs +,-,+ or -,+,-), '
                                            'so the pair is still real 1e-12 above the double zero' if sc[0] * sc[1] < 0 and sc[1] * sc[2] < 0 else 'NOT certified')
    log('local certificate at t_c = t_* +', (tca - ts).str(3, radius=False), 'signs', sc)
    R['double_zero']['t_star'] = ts.str(30, radius=False)
    R['double_zero']['t_star_minus_t_hi'] = (ts - tf(hi)).str(4, radius=False)
    R['double_zero']['t_lo_minus_t_star'] = (tf(lo) - ts).str(4, radius=False)
    log('double zero: z* =', zs.str(14, radius=False), 't* =', ts.str(20, radius=False), 'Krawczyk ok =', ok, 'contraction', [c_.str(3, radius=False) for c_ in contr])
    # the grid argument: a sign change is seen only when a grid point falls between the two zeros, D_theta/2 = sqrt(2 (t - t*))/2 > delta
    th_star = fl(zs / 2)
    dgrid = min(abs(th_star - (18460 + k) / 100.0) for k in range(131))
    R['grid_argument'] = dict(theta_star=repr(th_star), nearest_grid_distance=repr(dgrid), predicted_t_hi_minus_t_star_2delta2=repr(2 * dgrid ** 2),
                              worst_case_h2_over_2='5e-5', observed_t_hi_minus_t_star=(tf(hi) - ts).str(4, radius=False))

    # (c) the exact pair, for the right panel and for the comparisons
    gas = None
    gf = find_file(a.gas, 'phase_gas.json')
    if gf is not None:
        gas = json.load(open(gf))
    tlist = set()
    for sv in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.4):
        if ts_f + sv < 0.3:
            tlist.add(ts_f + sv)
    nuni = 12 if a.quick else 30
    for k in range(nuni + 1):
        tlist.add(ts_f + 0.005 + (0.3 - ts_f - 0.005) * k / nuni)
    tlist.add(0.0)
    extrap_t = []
    if gas:
        tst = float(gas['part_i']['t_stop'])
        extrap_t = [float(x) for x in __import__('numpy').linspace(tst, 0.0, 30)[:8]]
        tlist.update(extrap_t)
        gts = [float(x) for x in gas['samples']['t']]
        if not a.quick:
            tlist.update(x for x in gts if x > ts_f)
        else:
            tlist.update(x for x in gts[::3] if x > ts_f)
    tl = sorted(x for x in tlist if x > ts_f)
    ta = [arb(x) for x in tl]
    log('tracking the exact pair at', len(ta), 'times')
    zlo = track(zs, ts, ta, -1)
    zhi = track(zs, ts, ta, +1)
    D2 = [((h_ - l_) ** 2) for l_, h_ in zip(zlo, zhi)]
    R['track'] = dict(t=[x.str(15, radius=False) for x in ta], z_lo=[x.str(20, radius=False) for x in zlo], z_hi=[x.str(20, radius=False) for x in zhi],
                      D2=[x.str(14, radius=False) for x in D2],
                      D2_over_8_t_minus_tstar=[(d / (8 * (t_ - ts))).str(10, radius=False) for d, t_ in zip(D2, ta)])
    # t = 0 value of the tracked pair against the Newton zeros of S_0
    i0 = tl.index(0.0)
    R['pair']['tracked_D0_squared_at_t0'] = D2[i0].str(15, radius=False)
    R['pair']['tracked_minus_direct'] = abs(D2[i0] - D0 * D0).str(3, radius=False)
    iso = -D0 * D0 / 8
    rel = lambda x: ((x - ts) / abs(ts)).str(6, radius=False)
    R['comparisons'] = dict(isolated_estimate=iso.str(10, radius=False), isolated_minus_tstar_relative=rel(iso),
                            note='positive: the isolated pair estimate is later (closer to 0) than the collision; the chapter says 9 percent late')
    # the pair law check: D^2 / (8 (t - t*)) tends to 1
    R['pair_law'] = dict(ratio_D2_over_8s=[[ta[i].str(8, radius=False), R['track']['D2_over_8_t_minus_tstar'][i]] for i in range(0, len(ta), max(1, len(ta) // 12))])
    if gas:
        ex = gas['pair_71_72']['extrapolated_collision_first8of30']
        R['comparisons']['extrapolation_of_integrated_D2_ode'] = ex
        R['comparisons']['extrapolation_minus_tstar_relative'] = rel(arb(ex))
        import numpy as np
        idx = [tl.index(x) for x in extrap_t]
        p = np.polyfit(np.array(extrap_t), np.array([fl(D2[i]) for i in idx]), 1)
        R['comparisons']['extrapolation_of_exact_D2_same_8_times'] = repr(float(-p[1] / p[0]))
        R['comparisons']['extrapolation_exact_minus_tstar_relative'] = rel(arb(float(-p[1] / p[0])))
        # ODE against exact D^2 at the sample times
        od = []
        for sv, dv, tt in zip(gas['samples']['t'], gas['samples']['D2_pair_71_72'], gts):
            if tt > ts_f and tt in tl:
                od.append((tt, float(dv), fl(D2[tl.index(tt)])))
        if od:
            R['comparisons']['ode_vs_exact_D2'] = dict(n=len(od), max_abs_diff=repr(max(abs(x[1] - x[2]) for x in od)),
                                                       max_rel_diff=repr(max(abs(x[1] - x[2]) / x[2] for x in od)))
        R['figure_data'] = dict(gas=dict(t=gas['samples']['t'], s_min=gas['samples']['s_min'], csv_min=gas['samples']['csv_min'],
                                         crossing=(gas.get('crossing_4_5') or {}).get('with_tail'), N=gas['parameters']['N']),
                                ode_D2=dict(t=[repr(x[0]) for x in od], D2=[repr(x[1]) for x in od]))
    else:
        R['figure_data'] = {}
        log('phase_gas.json not found: the left panel of the figure is omitted')
    bf = find_file(a.book, 'dbn_collisions.json')
    if bf is not None:
        B = json.load(open(bf))
        rows = [x for x in B.get('collisions among the first 200 gaps on [-2.05, 0], latest first', []) if x.get('pair') == [71, 72]]
        if rows:
            bt = rows[0]['t_coll']
            R['comparisons']['book_dbn_collisions_t_coll'] = bt
            R['comparisons']['this_minus_book'] = (ts - arb(bt)).str(4, radius=False)
            R['comparisons']['book_z'] = rows[0]['z']
    R['seconds'] = round(time.time() - T0, 1)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(R, open(out, 'w'), indent=1)
    log('wrote', out)
    try:
        draw(R, fig)
        log('wrote', fig)
    except ImportError:
        log('matplotlib is missing: rerun with --figure where it is installed')


if __name__ == '__main__':
    main()
