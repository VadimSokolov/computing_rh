# ev_rouche5.py: the Rouche circle test for the five term folded kernel (E:425-450, tab:rouche5, M7 of review A2): for each T of the table, N(T), the number
# of sign changes of H^F_5(2t) on (0, T), and min over theta in [0, pi/2] of |H(2T e^(i theta))| / |E_5(2T e^(i theta))|; and the constants c = T/(N+1)^2.
# Claim supported: "every zero of H^F_5(2t) with |t| < 148.7 is real, c >= 4.13" and the certified constants 3.92, 3.62, 3.89, 4.13.
# Conventions. H(z) = int_0^oo Phi(u) cos(z u) du = Xi(z/2)/8 (ev_evt_lib); H^F_N(z) = (1/2)(F_N(z) + F_N(-z)), F_N(z) = sum_{n<=N} [ (1/2) a_n^(2-w9) Gamma(w9, a_n)
# - (3/4) a_n^(1-w5) Gamma(w5, a_n) ], a_n = pi n^2, w_j = (j + i z)/4 (the analytic form of the chapter's Re on the real axis); E_5 = sum_{n=6}^{14} of
# the same terms (the chapter truncates at n = 14; the neglected n = 15 term is checked to be below 10^-40 relative at theta = 0 here by a second run with n <= 18).
# Method. mpmath at dps digits, no interval arithmetic: the result is COMPUTED on the sampled circle (theta_j = (pi/2) j/J), not certified. The sign changes are counted
# on the grid t = h j, h = 0.03, up to 160. Check of the transform: H^F_N with N = 14 against Xi/8 at three real points.
# Error control. The minimum ratio is also recomputed at dps + 20 on the minimising theta and the relative difference is logged.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 -I ev_rouche5.py [--quick] [--rows I,J,..] [--thetas J] [--dps D] [--step h] [--no-signs] [--out FILE]
#   quick: dps 30, J = 24, sign grid step 0.12 (about 20 s); full: dps 40, J = 720 (721 points), step 0.03, on the cluster.
import argparse, json, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, mpc, pi, exp, gammainc, nstr, sin, cos, fabs
from ev_evt_lib import *

ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--rows', default=None)
ap.add_argument('--thetas', type=int, default=None); ap.add_argument('--dps', type=int, default=None); ap.add_argument('--step', type=float, default=None)
ap.add_argument('--no-signs', action='store_true'); ap.add_argument('--tmin', type=float, default=0); ap.add_argument('--tmax', type=float, default=160); ap.add_argument('--sdps', type=int, default=70); ap.add_argument('--out', default=None)
args = ap.parse_args()
dps = args.dps or (30 if args.quick else 40)
J = args.thetas or (24 if args.quick else 720)
step = args.step or (0.12 if args.quick else 0.03)
mp.dps = dps
ROWS = [('128.548', 42, '6.0e6'), ('130.333', 43, '7.7e5'), ('132.293', 44, '3.6e5'), ('134.127', 45, '3.6e4'), ('136.436', 46, '3.5e4'), ('138.926', 47, '1190'),
        ('140.430', 48, '247'), ('142.118', 49, '157'), ('144.556', 50, '57.1'), ('146.712', 51, '3.11'), ('148.738', 52, '1.36'), ('150.489', 53, '0.052')]
rows = ROWS if not args.rows else [ROWS[int(i)] for i in args.rows.split(',')]
logf = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.append(s)

def Fn(z, n):
    a = pi * n * n
    w9 = (9 + 1j * z) / 4; w5 = (5 + 1j * z) / 4
    return a ** (2 - w9) * gammainc(w9, a) / 2 - mpf(3) / 4 * a ** (1 - w5) * gammainc(w5, a)

def Hf(z, n0, n1):
    return sum((Fn(z, n) + Fn(-z, n)) / 2 for n in range(n0, n1 + 1))

def Hexact(z):
    return xi(mpf(1) / 2 + 1j * z / 2) / 8

t0 = time.time()
for tt in ('10', '25.3', '60.7'):
    z = 2 * mpf(tt)
    P('transform check t =', tt, 'rel diff of H^F_14 and Xi/8:', nstr(abs(Hf(z, 1, 14) - Hexact(z).real) / abs(Hexact(z).real), 3))
g = [mpf(x) for x in zero_strings(60)]
res = {}
for (Ts, N, chm) in (rows if not (args.no_signs is False and args.rows == '11' and J == 0) else []):
    T = mpf(Ts)
    cnt = sum(1 for x in g if x < T)
    between = (g[N - 1] < T < g[N]) if N < len(g) else None
    best = None; bestth = None
    for j in range(J + 1):
        th = pi / 2 * j / J
        z = 2 * T * exp(1j * th)
        H = Hexact(z)
        E = Hf(z, 6, 14)
        r = abs(H) / abs(E)
        if best is None or r < best: best, bestth = r, j
    # recompute at higher precision at the minimising theta
    mp.dps = dps + 20
    th = pi / 2 * bestth / J; z = 2 * T * exp(1j * th)
    r2 = abs(Hexact(z)) / abs(Hf(z, 6, 14))
    r3 = abs(Hexact(z)) / abs(Hf(z, 6, 18))
    mp.dps = dps
    P('T', Ts, 'N(T) from zeros file', cnt, 'chapter', N, 'T between gamma_N and gamma_N+1:', between, '| min |H|/|E5| =', nstr(best, 4), 'at theta index', bestth, 'of', J,
      '| chapter', chm, '| rel diff at dps+20', nstr(abs(best - r2) / best, 2), '| with n<=18', nstr(abs(best - r3) / best, 2), '| c = T/36 =', nstr(T / 36, 4))
    res[Ts] = {'N': cnt, 'min': nstr(best, 6), 'theta_index': bestth, 'chapter': chm}
cc = {'2': (('35.3'), 9), '3': ('57.9', 16), '4': ('97.3', 25), '5': ('148.738', 36)}
P('c = T/(N+1)^2 for the last passing rows:', {k: nstr(mpf(v[0]) / v[1], 4) for k, v in cc.items()})
if not args.no_signs:
    h = mpf(str(step)); prev = None; ch = []; t = mpf(args.tmin) + h / 10
    mp.dps = args.sdps
    xs = []
    while t < args.tmax:
        v = Hf(2 * t, 1, 5).real
        if prev is not None and v * prev < 0: xs.append(t)
        prev = v; t += h
    P('sign changes of H^F_5(2t) on (%g,%g) at %d digits, grid step' % (args.tmin, args.tmax, args.sdps), step, ':', len(xs), '(chapter 53 below 160); zeros of Xi in the range:', sum(1 for x in g if args.tmin < x < args.tmax))
    near = [min(abs(x - y) for y in g) for x in xs] or [mpf(0)]
    P('max distance of a sign change to an ordinate of Xi:', nstr(max(near), 3), '(chapter: within grid resolution 0.03)')
    res['signs'] = {'count': len(xs), 'maxdist': nstr(max(near), 3)}
P('time', round(time.time() - t0, 1), 's')
of = args.out or ('/tmp/ev_rouche5.json' if args.quick else os.path.join(HERE, '..', 'data', 'ev_rouche5.json'))
json.dump(res, open(of, 'w'), indent=1)
open(of.replace('.json', '.log'), 'w').write('\n'.join(logf) + '\n')
