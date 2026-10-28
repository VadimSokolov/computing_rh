# ev_dbn.py: kernel moments, the mean mu_1(lambda), Hankel determinants of the de Bruijn and Newman family and of the two term folded kernel, the pair at
# lambda = -4 and the masking table D_6 (E:380-400, 480-484, 487-512, 637-641 of Chapter 25 (ev:)).
# Claims supported: mu_1 = .0232010, .0231050, .0226358, .0190909 at lambda = .2, 0, -1, -10 and dmu_1/dlambda = 4.78e-4 at 0 (E:392); D_n, D_n' > 0 for n <= 5 on
# lambda in [-10, .2] (E:392); the Hankel signs at lambda = 0, -4, -8 to order 10 (E:484: D_10' < 0 at -4, D_8 to D_10 and primes < 0 at -8, all positive at 0);
# the folded N = 2 kernel positive to order 11 at 120 and 240 digits (E:400); H_{-4}(2t) < 0 on [28, 38] and the zero pair t = 33.1150940 +- .1505390 i (E:639-641);
# first rows of tab:mask and the first two zeros of H^F_2, t = 43.13891 - 3.28097 i, 47.52276 - 6.25251 i (E:398, 502-506).
# Conventions. Phi as in ev_evt_lib; H_lambda(z) = int_0^oo e^(lambda u^2) Phi(u) cos(z u) du; M_{2m} its even moments; with H_lambda(z) = H_lambda(0) prod (1 - z^2/(4 gamma_k^2)),
# H(2 i sqrt s)/H(0) = prod (1 + s/gamma_k^2), so the power sums p_j = sum_k gamma_k^(-2j) follow from the logarithm of sum_m M_{2m} 4^m s^m/((2m)! M_0);
# mu_1 = p_1 = 2 M_2/M_0 = 2 E U^2; D_n = det[p_{i+j+1}], D_n' = det[p_{i+j+2}] (i, j < n). Atoms x = 1/t^2 with t a zero of H(2t) (z = 2t).
# Method. Gauss-Legendre on [0, 1.5] (Phi(1.5) < e^(-1200)), nodes from mpmath, 3 * 2^(deg-1) nodes; the exact kernel summed over n <= 40 (terms below e^(-pi n^2) vanish);
# zeros of H^F_2(2t) by secant (mp.findroot) from the ordinates; D_6 of an atom set by the Hankel determinant of its power sums (equal to the Andreief sum).
# Error control. Each family is computed at two node counts and two precisions; the digits agreeing in D_n are logged (min over n and families).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Further parts: argp (argument principle for H^F_2(2t) on rectangles of the t plane, E:396-400, 489: counts of zeros in [36,150] x [-8,-0.05], [36,150] x [-40,-8],
# [0,36] x [-40,-0.05], by the accumulated phase of H^F_2 along the boundary, adaptive step, dps 30: computed, not certified) and counts (sign changes of H_lambda(2t) on a grid of t below tmax,
# E:482: 29, 29, 27, 19, 12, 8 real zeros below 100 at lambda = -0.25, -0.5, -1, -2, -4, -8).
# Usage: python3 -I ev_dbn.py [--quick] [--part moments|hankel|folded|pair|mask|argp|counts|all] [--rect X0,X1,Y0,Y1] [--lams L1,L2] [--tmax T] [--tstep h] [--dps D] [--deg G] [--out FILE]
import argparse, json, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, mpc, pi, exp, nstr, findroot, cos, sign, det, matrix, log10
from ev_evt_lib import *

ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--part', default='all'); ap.add_argument('--dps', type=int, default=None)
ap.add_argument('--deg', type=int, default=None); ap.add_argument('--rect', default='36,150,-8,-0.05'); ap.add_argument('--lams', default='-0.25,-0.5,-1,-2,-4,-8'); ap.add_argument('--tmax', type=float, default=100); ap.add_argument('--tstep', type=float, default=0.05); ap.add_argument('--out', default=None)
args = ap.parse_args()
dps = args.dps or (60 if args.quick else 130); deg = args.deg or (7 if args.quick else 9)
logf = []; out = {}
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.append(s)
def part(n): return args.part == n or (args.part == 'all' and n not in ('argp', 'counts'))
t0 = time.time()
mp.dps = dps
nodes = gl_nodes(deg, mpf('1.5'))
P('nodes', len(nodes), 'dps', dps, 'time', round(time.time() - t0, 1))
mmax = 21

def agree(a, b):
    d = abs(a - b)
    return 999 if d == 0 else float(-log10(d / abs(a)))

if part('moments') or part('hankel'):
    res = {}
    for lam in ('0.2', '0', '-1', '-10'):
        M = kernel_moments(nodes, mpf(lam), 40, 2)
        mu1 = 2 * M[1] / M[0]; dmu = 2 * (M[2] / M[0] - (M[1] / M[0]) ** 2)
        P('lambda', lam, 'mu_1 =', nstr(mu1, 10), ' dmu_1/dlambda =', nstr(dmu, 5))
        res[lam] = (nstr(mu1, 10), nstr(dmu, 5))
    out['mu1'] = res
if part('hankel'):
    # positivity n <= 5 on a lambda grid
    allpos = True; grid = [mpf(x) for x in ['0.2', '0', '-0.5', '-1', '-2', '-3', '-4', '-5', '-6', '-7', '-8', '-9', '-10']]
    for lam in grid:
        M = kernel_moments(nodes, lam, 40, 11); p = power_sums(M, 11)
        D, Dp = hankels(p, 5)
        ok = all(D[n] > 0 and Dp[n] > 0 for n in range(1, 6))
        allpos &= ok
        if not ok: P('  non positive at lambda', lam)
    P('D_n, D_n\' > 0 for n <= 5 on the lambda grid:', allpos)
    out['pos5'] = allpos
    for lam in ('0', '-4', '-8'):
        M = kernel_moments(nodes, mpf(lam), 40, mmax); p = power_sums(M, mmax)
        D, Dp = hankels(p, 10)
        s = ''.join('+' if D[n] > 0 else '-' for n in range(1, 11)); sp = ''.join('+' if Dp[n] > 0 else '-' for n in range(1, 11))
        P('lambda', lam, 'signs D_1..D_10', s, " D'_1..D'_10", sp)
        out['signs_' + lam] = (s, sp)
        # second node count
        mp2 = gl_nodes(deg - 1, mpf('1.5'))
        M2 = kernel_moments(mp2, mpf(lam), 40, mmax); p2 = power_sums(M2, mmax); D2, Dp2 = hankels(p2, 10)
        P('   digits of agreement D_n between %d and %d nodes (n = 1..10):' % (len(nodes), len(mp2)), [round(min(agree(D[n], D2[n]), 99), 0) for n in range(1, 11)])
if part('folded'):
    # folded N = 2 kernel: positivity to order 11 (needs p_j to 23)
    for N_dps, N_deg in ((dps, deg),) if args.quick else ((120, 9), (240, 10)):
        mp.dps = N_dps; nd = gl_nodes(N_deg, mpf('1.5'))
        M = kernel_moments(nd, mpf(0), 2, 24); p = power_sums(M, 24); D, Dp = hankels(p, 11)
        P('folded N = 2, dps', N_dps, 'nodes', len(nd), 'min D_n sign positive for n <= 11:', all(D[n] > 0 for n in range(1, 12)), "D' positive:", all(Dp[n] > 0 for n in range(1, 12)))
        P('   D_1..D_11 exponents', [nstr(D[n], 3) for n in range(1, 12)])
        out['folded_%d' % N_dps] = [nstr(D[n], 6) for n in range(1, 12)]
        D6_kernel = D[6]
    mp.dps = dps
    out['D6_kernel'] = nstr(D6_kernel, 6)
    P('D_6 from kernel moments (all zeros), chapter 8.10e-100:', nstr(D6_kernel, 4))
if part('pair'):
    mp.dps = 30
    # H_{-4}(2t) on [28, 38]
    nd30 = gl_nodes(8, mpf('1.5'))
    def Hl(z, lam=-4):
        return sum(w * exp(lam * u * u) * phi_terms(u, 12) * cos(z * u) for u, w in nd30)
    tt = [mpf(28) + mpf(j) / 4 for j in range(41)]
    vals = [Hl(2 * t) for t in tt]
    P('H_{-4}(2t) on [28, 38] step 0.25: max value', nstr(max(vals), 4), ' (negative: %s), min |value| near t = %s' % (all(v < 0 for v in vals), nstr(tt[min(range(41), key=lambda i: abs(vals[i]))], 4)))
    z = findroot(lambda t: Hl(2 * t), mpc('33.1', '0.15'))
    P('root of H_{-4}(2t) near 33.1 + 0.15 i:', nstr(z, 12), 'gamma^2 =', nstr(z * z, 10), ' chapter 33.1150940 + .1505390 i, 1096.587 + 9.970 i')
    out['root_m4'] = nstr(z, 12)
if part('mask'):
    mp.dps = 60
    # real zeros of H^F_2(2t) (7 below 70) and the first two off line zeros
    mp.dps = 40
    f2 = lambda t: HF(2 * t, 1, 2).real
    grid = [mpf(10) + mpf(j) / 20 for j in range(1201)]
    fv = [f2(t) for t in grid]
    real = []
    for j in range(1200):
        if fv[j] * fv[j + 1] < 0:
            real.append(findroot(f2, (grid[j], grid[j + 1]), solver='anderson', tol=1e-30))
    mp.dps = 60
    P('real zeros of H^F_2(2t) (sign changes on [10, 70], step 0.05):', len(real), [nstr(r, 8) for r in real])
    z1 = findroot(lambda t: HF(2 * t, 1, 2), mpc('43.14', '-3.28'), tol=1e-25)
    z2 = findroot(lambda t: HF(2 * t, 1, 2), mpc('47.52', '-6.25'), tol=1e-25)
    P('off line zeros of H^F_2(2t):', nstr(z1, 8), nstr(z2, 8), 'chapter 43.13891 - 3.28097 i, 47.52276 - 6.25251 i')
    mp.dps = 250
    def D6(atoms):
        p = [None] + [sum(a ** j for a in atoms) for j in range(1, 12)]
        return det(matrix([[p[i + j + 1] for j in range(6)] for i in range(6)]))
    ra = [1 / (r * r) for r in real][:7]
    P('7 real atoms: count', len(ra))
    d_a = D6(ra)
    ra1 = ra + [1 / (z1 * z1), 1 / (z1.conjugate() ** 2)]
    d_b = D6(ra1)
    ra2 = ra1 + [1 / (z2 * z2), 1 / (z2.conjugate() ** 2)]
    d_c = D6(ra2)
    P('D_6: 7 real atoms', nstr(d_a, 4), '(chapter 1.29e-105); plus first pair', nstr(d_b, 4), '(-1.81e-104); plus first two pairs', nstr(d_c, 4), '(-2.27e-103)')
    out['mask'] = [nstr(d_a, 6), nstr(d_b, 6), nstr(d_c, 6)]
if part('argp'):
    from mpmath import arg
    mp.dps = 30
    x0, x1, y0, y1 = [mpf(v) for v in args.rect.split(',')]
    f = lambda t: HF(2 * t, 1, 2)
    nev = [0]
    def seg(a, b, fa, fb, depth):
        nev[0] += 1
        da = arg(fb / fa)
        if abs(da) > 0.7 and depth < 14:
            m = (a + b) / 2; fm = f(m)
            return seg(a, m, fa, fm, depth + 1) + seg(m, b, fm, fb, depth + 1)
        return da
    def edge(a, b, n):
        tot = mpf(0); pts = [a + (b - a) * j / n for j in range(n + 1)]; fv = [f(t) for t in pts]
        for j in range(n): tot += seg(pts[j], pts[j + 1], fv[j], fv[j + 1], 0)
        return tot
    c = [mpc(x0, y0), mpc(x1, y0), mpc(x1, y1), mpc(x0, y1), mpc(x0, y0)]
    total = sum(edge(c[i], c[i + 1], max(2, int(abs(c[i + 1] - c[i]) / 0.2))) for i in range(4))
    P('argument principle, rectangle', args.rect, ': winding number', nstr(total / (2 * pi), 8), 'evaluations', nev[0] + 0)
    out['argp_' + args.rect] = nstr(total / (2 * pi), 8)
if part('counts'):
    mp.dps = dps
    nd = gl_nodes(deg, mpf('1.5'))
    for lam in [mpf(x) for x in args.lams.split(',')]:
        wu = [(u, w * exp(lam * u * u) * phi_terms(u, 40)) for u, w in nd]
        prev = None; n = 0; locs = []; t = mpf(args.tstep) / 3
        while t < args.tmax:
            v = sum(w * cos(2 * t * u) for u, w in wu)
            if prev is not None and v * prev < 0: n += 1; locs.append(t)
            prev = v; t += mpf(args.tstep)
        P('lambda', lam, 'sign changes of H_lambda(2t) for t <', args.tmax, 'step', args.tstep, ':', n, ' (first gaps:', [nstr(locs[i + 1] - locs[i], 3) for i in range(min(len(locs) - 1, 4))], ')')
        out['counts_' + str(lam)] = n
P('time', round(time.time() - t0, 1))
of = args.out or ('/tmp/ev_dbn.json' if args.quick else os.path.join(HERE, '..', 'data', 'ev_dbn_%s.json' % args.part))
json.dump(out, open(of, 'w'), indent=1)
open(of.replace('.json', '.log'), 'w').write('\n'.join(logf) + '\n')
