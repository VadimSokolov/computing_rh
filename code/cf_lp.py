# cf_lp.py: the numerical claims of Table cf:tab:lp of Chapter 34 (Cauchy deconvolution from U_1 to U_{1/2}, prefix cf): the ceiling of the
# scheme of Corollary cf:cor:ystar over all smoothings phi >= 0, phi(0) = 1, with positive definite even extension, computed by a discretised
# linear program (the text of the chapter, L:1421-1453 of the source): how the digits depend on the frequency grid, on Lambda and on the
# number of lambda points. Added for the book in October 2026 (group CF-LP), written from the mathematics (the collaborator's scripts were
# not read, run or copied).
# The program. For a fixed cosine polynomial p = sum_{k<=K} a_k cos(k phi) >= 0 (a_k >= 0, A = a_1 + ... + a_K) and a fixed y,
#   G(y) = max_phi { a_1 int_0^Lambda exp(-y l) phi - a_0 int_0^Lambda phi }
#   subject to 0 <= phi <= 1, phi(0) = 1, int_0^Lambda cos(u l) phi(l) dl >= 0 for u in a finite frequency set,
# and y is certified when G(y) >= A/2; G decreases in y, so the ceiling c is the root of G(y) = A/2, found by bisection (tolerance 1e-6, and
# 1e-8 in the error check). With p free (degree <= K) the two programs (phi for fixed p, p for fixed phi, with a_0 = 1, p >= 0 on 400 angles
# and 0 <= a_k <= 4) alternate inside the bisection from several starting polynomials (best single level of degree K, 3+4cos+cos2, a flat one).
# Two discretisations, because the text does not say how the integrals are discretised: lambda_i = i Lambda/(N-1), i < N.
#   trap  all integrals (the objective and the constraints) are trapezoid sums over the lambda grid;
#   pl    phi is the piecewise linear interpolant of its values: the objective by trapezoid sums (exact for the interpolant up to the
#         exponential, error below 1e-7 here) and the cosine transform exactly (closed form by parts), which is what a continuous
#         phi would give. For trap the constraints at u beyond pi/h (about 78 for Lambda = 12, N = 300) are aliased copies of lower ones.
# Frequency sets. unif: 460 equally spaced points on [0, 400] (the chapter's stated grid). c25: unif plus 239 equally spaced points on
# (0, 25] (699 in all); c40: unif plus 939 on (0, 40] (1399 in all). These counts reproduce the three grids described in the reconcile
# note; their exact shape in the source is unknown, so this is the grid of this script.
# Rows. P1: p = 3 + 4 cos + cos 2 (a = (3, 4, 1)); P4: the degree 4 polynomial of Table cf:tab:arith (the best single level polynomial,
# computed here by Nelder Mead on y* = 2 (sqrt a1 - sqrt a0)^2/A over Fejer Riesz b >= 0); both fixed. F4 and F8: polynomial of degree <= 4
# and <= 8 optimised together with phi. The cap of Proposition cf:prop:barrier2 is (2/A)(a1 - a0 - a0 log(a1/a0)).
# Error control. The stated row is recomputed with bisection tolerance 1e-8 and with the interior point method of HiGHS instead of the dual
# simplex, and the differences are logged; the sweep over grids, Lambda and N is itself the discretisation test.
# Output: data/cf_lp.json and the log data/cf_lp.log (or --out FILE and FILE with .log); --quick takes seconds.
import sys, json, math, time, argparse, os
from multiprocessing import Pool
import numpy as np
from scipy.optimize import linprog, minimize

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default='data/cf_lp.json')
ap.add_argument('--part', default='all', help='all, fixed or opt')
ap.add_argument('--procs', type=int, default=1)
args = ap.parse_args()
QUICK = args.quick
LOGF = open(args.out.rsplit('.', 1)[0] + '.log', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOGF.write(s + '\n'); LOGF.flush()
T0 = time.time()
RES = {'quick': QUICK, 'fixed': [], 'opt': []}

def coeffs(b):
    K = len(b) - 1
    a = np.zeros(K + 1); a[0] = np.dot(b, b)
    for k in range(1, K + 1): a[k] = 2 * np.dot(b[:-k], b[k:])
    return a
def single_closed(a):
    a0, a1, A = a[0], a[1], a[1:].sum()
    return 2 * (math.sqrt(a1) - math.sqrt(a0)) ** 2 / A if a1 > a0 else -(a0 - a1) / A
def best_single(K, nstart=8):
    rng = np.random.default_rng(1); best = (-1, None)
    for _ in range(nstart):
        v0 = rng.random(K + 1) + 0.05
        f = lambda v: -single_closed(coeffs(v ** 2))
        r = minimize(f, v0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 6000, 'maxfev': 6000})
        r = minimize(f, r.x, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-16, 'maxiter': 6000, 'maxfev': 6000})
        if -r.fun > best[0]: best = (-r.fun, coeffs(r.x ** 2))
    return best[1] / best[1][0]
def cap_new(a):
    return 2 / a[1:].sum() * (a[1] - a[0] - a[0] * math.log(a[1] / a[0]))

def freq_grid(name):
    q = QUICK   # quick: about a quarter of the points, only to test the code path
    unif = np.linspace(0, 400, 120 if q else 460)
    if name == 'unif': return unif
    if name == 'c25': return np.unique(np.concatenate([unif, np.linspace(0, 25, 40 if q else 240)]))
    if name == 'c40': return np.unique(np.concatenate([unif, np.linspace(0, 40, 80 if q else 940)]))
    raise ValueError(name)

class Disc:
    def __init__(self, Lam, N, fname, variant, nang=400):
        self.Lam, self.N, self.fname, self.variant = Lam, N, fname, variant
        lam = np.linspace(0, Lam, N); h = Lam / (N - 1); self.lam, self.h = lam, h
        w = np.full(N, h); w[0] = w[-1] = h / 2; self.w = w
        u = freq_grid(fname); self.u = u
        if variant == 'trap':
            C = np.cos(np.outer(u, lam)) * w[None, :]
        else:
            C = np.zeros((len(u), N)); z = u == 0
            C[z] = w[None, :]
            uu = u[~z]
            cl = np.cos(np.outer(uu, lam))                       # cos(u lam_k)
            E = cl[:, :-1] - cl[:, 1:]                           # E_k, k = 0..N-2
            Cc = np.zeros((len(uu), N))
            Cc[:, 1:-1] += (E[:, :-1] - E[:, 1:]) / h             # coefficient of phi_i, 1 <= i <= N-2
            Cc[:, 0] += -E[:, 0] / h
            Cc[:, -1] += E[:, -1] / h
            Cc *= (-1.0 / uu ** 2)[:, None]
            Cc[:, -1] += np.sin(uu * Lam) / uu
            C[~z] = Cc
        self.C = C
        self.ph = np.linspace(0, math.pi, nang)
        self.cosk = None
    def lp_phi(self, a, y, method='highs-ds'):
        e = np.exp(-y * self.lam) * self.w
        cost = -(a[1] * e - a[0] * self.w)                       # minimise -(a1 int e phi - a0 int phi)
        const = (a[1] * e[0] - a[0] * self.w[0])                 # phi_0 = 1
        S = 1 / np.abs(self.C).max(axis=1)                       # row scaling (the same feasible set)
        r = linprog(cost[1:], A_ub=-self.C[:, 1:] * S[:, None], b_ub=self.C[:, 0] * S, bounds=(0, 1), method=method)
        if r.status != 0: return None, -1e9
        phi = np.concatenate([[1.0], r.x])
        return phi, float(-r.fun + const)
    def lp_poly(self, phi, y, K, method='highs-ds'):
        E1 = float(np.sum(np.exp(-y * self.lam) * self.w * phi)); E0 = float(np.sum(self.w * phi))
        ks = np.arange(K + 1)
        Cos = np.cos(np.outer(self.ph, ks))
        cost = np.zeros(K + 1); cost[1:] = 0.5 * 1.0 - E1 * (np.arange(K + 1)[1:] == 1)   # minimise -(a1 E1 - a0 E0 - A/2) = -a1 E1 + a0 E0 + A/2
        cost[1] = 0.5 - E1; cost[0] = E0
        r = linprog(cost, A_ub=-Cos, b_ub=np.zeros(len(self.ph)), A_eq=np.eye(1, K + 1), b_eq=[1.0], bounds=[(None, None)] + [(0, 4)] * K, method=method)
        if r.status != 0: return None, -1e9
        a = r.x
        return a, float(a[1] * E1 - a[0] * E0 - a[1:].sum() / 2)

def newton_ceiling(D, a, y0=0.0, method='highs-ds', tol=1e-9, maxit=40):
    """Root of G(y) = A/2 by Newton from the left: G is convex and decreasing (a maximum of functions convex in y), G'(y) = -a1 int l e^{-yl} phi*
    (envelope theorem), so the iterates increase to the root and every iterate is certified (G >= A/2)."""
    A = a[1:].sum(); y = y0; phi = None
    for it in range(maxit):
        ph, g = D.lp_phi(a, y, method)
        if ph is None: return y, phi
        f = g - A / 2
        if f < -1e-9:
            if phi is None: return 0.0, None
            return y, phi
        phi = ph
        dG = -a[1] * float(np.sum(D.lam * np.exp(-y * D.lam) * D.w * phi))
        if f < tol or dG >= -1e-14: return y, phi
        y = y + f / (-dG)
    return y, phi

def ceiling_fixed(D, a, method='highs-ds', tol=1e-9):
    return newton_ceiling(D, a, 0.0, method, tol)

def ceiling_opt(D, K, starts, method='highs-ds', outer=14):
    best = (0.0, None, None)
    for a0 in starts:
        a = np.array(a0, float).copy()
        if len(a) < K + 1: a = np.concatenate([a, np.zeros(K + 1 - len(a))])
        a = a / a[0]
        y, phi = newton_ceiling(D, a, 0.0, method)
        for it in range(outer):
            if phi is None: break
            a2, g2 = D.lp_poly(phi, y, K, method)
            if a2 is None or g2 < -1e-9: break
            y2, phi2 = newton_ceiling(D, a2, y, method)
            if phi2 is None: break
            gain = y2 - y
            if gain >= 0: a, y, phi = a2, y2, phi2
            if gain < 1e-7: break
        if y > best[0]: best = (y, a / a[0], phi)
    return best

def run_task(t):
    t0 = time.time()
    D = Disc(t['Lambda'], t['N'], t['grid'], t['variant'])
    if t['kind'] == 'fixed':
        a = np.array(t['a'])
        c, _ = ceiling_fixed(D, a)
        out = dict(t, c=c, nfreq=len(D.u))
    elif t['kind'] == 'errcheck':
        a = np.array(t['a'])
        c1, _ = ceiling_fixed(D, a, tol=1e-6); c2, _ = ceiling_fixed(D, a, tol=1e-11); c3, _ = ceiling_fixed(D, a, 'highs-ipm')
        out = dict(t, c_1e6=c1, c_1e11=c2, c_ipm=c3, nfreq=len(D.u))
    else:
        c, a, _ = ceiling_opt(D, t['K'], t['starts'])
        out = dict(t, c=c, a=None if a is None else a.tolist(), nfreq=len(D.u))
    out['sec'] = time.time() - t0
    return out

def fmt(o):
    if o['kind'] == 'errcheck':
        return 'error check P1 stated grid, pl: Newton tol 1e-6 %.8f | tol 1e-11 %.8f (diff %.2e) | interior point %.8f (diff %.2e)' % (o['c_1e6'], o['c_1e11'], abs(o['c_1e6'] - o['c_1e11']), o['c_ipm'], abs(o['c_ipm'] - o['c_1e11']))
    c = o['c']
    if o['kind'] == 'fixed':
        lab = '%s fixed poly' % o['poly']
    else:
        lab = 'F%d polynomial and phi optimised' % o['K']
    s = '%s | grid %-4s (%4d freq) | %-4s | Lambda %2d N %3d | c = %.5f  1/c = %.3f' % (lab, o['grid'], o['nfreq'], o['variant'], o['Lambda'], o['N'], c, 1 / c if c else float('inf'))
    if o['kind'] == 'opt' and o['a'] is not None: s += '  a = %s' % np.round(o['a'], 3).tolist()
    return s + '  (%.0fs)' % o['sec']

if __name__ == '__main__':
    if QUICK:
        LAMN = [(12, 80)]; GRIDS = ['unif', 'c25']; VARS = ['trap', 'pl']; OPT_CFG = [(12, 80, 'c25', 'pl')]
    else:
        LAMN = [(8, 300), (12, 300), (16, 300), (8, 600), (12, 600), (16, 600), (20, 600)]
        GRIDS = ['unif', 'c25', 'c40']; VARS = ['trap', 'pl']
        OPT_CFG = [(Lm, N, g, v) for (Lm, N) in [(12, 300), (16, 600)] for g in GRIDS for v in VARS]
    a1p = np.array([3.0, 4.0, 1.0]); a4 = best_single(4)
    log('cf_lp: quick' if QUICK else 'cf_lp: full', ' procs', args.procs)
    log('P4 (degree 4 single level optimum): a =', np.round(a4, 5).tolist(), ' y* single = %.6f, cap_new P1 %.5f, P4 %.5f' % (single_closed(a4), cap_new(a1p), cap_new(a4)))
    RES['P4'] = a4.tolist(); RES['cap_new'] = {'P1': cap_new(a1p), 'P4': cap_new(a4)}
    tasks = []
    if args.part in ('all', 'fixed'):
        for name, a in [('P1', a1p), ('P4', a4)]:
            for fname in GRIDS:
                for variant in VARS:
                    for (Lm, N) in LAMN:
                        tasks.append(dict(kind='fixed', poly=name, a=a.tolist(), grid=fname, variant=variant, Lambda=Lm, N=N))
        tasks.append(dict(kind='errcheck', a=a1p.tolist(), grid='unif', variant='pl', Lambda=12, N=LAMN[1][1] if len(LAMN) > 1 else LAMN[0][1]))
    if args.part in ('all', 'opt'):
        flat = [1.0, 1.5, 0.8, 0.3, 0.05, 0.02, 0.01, 0.005, 0.002]
        for K in ([4] if QUICK else [4, 8]):
            starts = [best_single(K, 4).tolist(), a1p.tolist(), flat[:K + 1]]
            for (Lm, N, fname, variant) in OPT_CFG:
                tasks.append(dict(kind='opt', K=K, starts=starts, grid=fname, variant=variant, Lambda=Lm, N=N))
    tasks.sort(key=lambda t: -(t['N'] * (1 + 2 * (t['kind'] == 'opt'))))   # longest first
    if args.procs > 1:
        with Pool(args.procs) as pool:
            outs = []
            for o in pool.imap_unordered(run_task, tasks):
                log(fmt(o)); outs.append(o)
    else:
        outs = []
        for t in tasks:
            o = run_task(t); log(fmt(o)); outs.append(o)
    RES['tasks'] = outs
    RES['seconds'] = time.time() - T0
    json.dump(RES, open(args.out, 'w'), indent=1, default=float)
    log('done in %.1f s' % (time.time() - T0))
