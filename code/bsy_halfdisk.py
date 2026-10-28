# bsy_halfdisk.py: the numerical claims of Section bsy:sec:halfdisk of Chapter 35 (Theorem bsy:thm:halfdisk, the Balazard, Saias and Yor identity on a half disc, and
# Lemma bsy:lem:trap on the sorted trapezoid rule), group BSY-E. It fills Table bsy:tab:halfdisk and the numbers of the text after it (the Green term, the agreement of the
# Gauss and Legendre reference with the right side, the slopes, the ratios at N = 32768), and supplies the data of fig/bsy_halfdisk.pdf (code/bsy_figs.py).
# Setting. c = theta + i t0, R = 1, t0 = 85.7, H = {|s - c| < R, Re s > theta}, z = c + R/2, w = (s - c)/R. Three cases: zeta at theta = 3/4 (h = zeta, nothing to subtract);
# zeta at theta = 1/2 (h = zeta, the zero on the diameter, refined by Newton's method, subtracted); the Davenport and Heilbronn function f at theta = 1/2 (h = f, with the zeros of f on
# the diameter, if any, subtracted and the zero 0.8085... + 85.6993... i of data/cf_dhzeros.json in H, so that the right side contains G_H(z, rho)).
# Identity (Theorem bsy:thm:halfdisk, all multiplicities 1 here): with u = log|h| - sum_j log|s - rho_j| over the zeros rho_j on the diameter,
#   int_{dH} u d omega_z = log|h(z)| - sum_j log|z - rho_j| + sum_{rho in H} G_H(z, rho)          (the right side RS, the target).
# Densities (w = e^{i phi} on the arc, |phi| <= pi/2; w = i y on the diameter, |y| <= 1; G_D(a, w) = log|(1 - conj(a) w)/(w - a)|, w_z = (z - c)/R, w* = -conj(w_z)):
#   arc:       d omega/d phi = (1 - |w_z|^2)/(2 pi |e^{i phi} - w_z|^2) - (1 - |w*|^2)/(2 pi |e^{i phi} - w*|^2)
#   diameter:  d omega/d y   = (1/2 pi) d_x [G_D(w_z, x + i y) - G_D(w*, x + i y)] at x = 0, with d_x G_D(a, w) = Re[-conj(a)/(1 - conj(a) w) - 1/(w - a)].
# The checks of the configuration (all logged and stored under 'config'): (i) total harmonic measure 1 by Gauss and Legendre; (ii) the identity itself, with the code of the densities and
# of G_H, for the test function (s - c - 0.2 i)(s - c - 0.3 - 0.1 i) (one zero on the diameter, one in H) at a node count and at twice it; (iii) the argument principle: the winding number of
# q = h / prod (s - rho_j) along the closed boundary dH (counter clockwise, 40000 points per piece, every principal increment of the argument below 1 in absolute value) is the number of
# zeros of h in H off the diameter (0, 0, 1) and min |q| on the boundary is positive; (iv) the zeros of h on the diameter: for zeta the zero 1/2 + 84.7354929805... i by Newton's method
# (compared with data/zeros_1700.txt), for the Davenport and Heilbronn function a scan of |f| on the diameter at step 1e-3 with Newton's method from every local minimum below 0.05, where a zero
# counts if Re rho = 1/2 within 1e-9 and |f| < 1e-12 (none exists in 84.7 < t < 86.7, which the winding number of (iii) confirms); the DH zero in H is taken from data/cf_dhzeros.json (tag z85).
# Reference. Gauss and Legendre on each piece with n = 64, 128, ..., 4096 nodes, doubling until two successive values of the integral differ by less than 1e-14 (the first such n is logged);
# the reference is compared with RS (error 'ref_minus_RS'). Evaluator precision: h is evaluated by the Euler and Maclaurin code of bsy_zdh.py (default nf = 0.3, M = 28); the reference and
# the RS are recomputed with the second setting nf = 0.45, M = 36, and with mpmath at 30 digits at 50 points of the boundary, and the differences are logged ('prec').
# Estimators (N = 2^7, ..., 2^15, 200 replications, seeds SeedSequence([20261025, case, rep]) stored): plain Monte Carlo (N independent points from omega_z by inverse distribution
# function on a grid of 2^17 points per piece, the piece chosen with its mass), the sorted trapezoid rule with N/2 uniform nodes on each piece (variable phi on the arc, y on the diameter) and
# the end points added, and the same after the change of variable y = Phi(v) = v + sin(pi v)/pi on [-1, 1] (the arc through phi = (pi/2) y), the integrand g(Phi(v)) Phi'(v). Per N: the root mean square error
# (RMSE), the mean error against RS, and the slope of log RMSE against log N for N >= 1024 (least squares, with its standard error); the ratios plain over weighted at N = 32768.
# The mean error is also compared with the leading bias term of Lemma bsy:lem:trap, (b-a)^2 (g'(b) - g'(a)) / (2 (N/2+1)(N/2+2)), summed over the two pieces, as a check of the lemma (key 'bias_lemma').
# Written from the mathematics (the collaborator's scripts were not read, run or copied). Options: --quick (N up to 2^10, 6 replications, seconds), --out FILE (json; the log goes beside it
# with the extension .log), --procs P (processes for the replications).
# Run: python3 -I code/bsy_halfdisk.py [--quick] [--out data/bsy_halfdisk.json] [--procs 16]   (the full run is on Hopper: bash misc/tools/hopper_run.sh -c 16 -t 60 -n bsy_halfdisk -g bsy_halfdisk.json -g bsy_halfdisk.log code/bsy_halfdisk.py code/bsy_zdh.py data/cf_dhzeros.json data/zeros_1700.txt -- --procs 16)
import argparse, json, math, os, sys, time
import numpy as np
from numpy.polynomial.legendre import leggauss
import multiprocessing as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bsy_zdh as Z

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default=None)
ap.add_argument('--procs', type=int, default=1)
ap.add_argument('--data', default=None, help='folder of cf_dhzeros.json and zeros_1700.txt (default: ../data, then the script folder)')
args = ap.parse_args()
OUT = args.out or ('misc/reviews/nick-b4/quick/BSY-E/bsy_halfdisk.json' if args.quick else 'data/bsy_halfdisk.json')
LOGF = os.path.splitext(OUT)[0] + '.log'
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
LOG = open(LOGF, 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.write(s + '\n'); LOG.flush()


def finddata(name):
    for d in (args.data, os.path.join(HERE, '..', 'data'), HERE, '.'):
        if d and os.path.exists(os.path.join(d, name)):
            return os.path.join(d, name)
    raise FileNotFoundError(name)


R = 1.0
T0 = 85.7
KIND = {0: ('zeta', 0.75), 1: ('zeta', 0.5), 2: ('dh', 0.5)}
NAME = {0: 'zeta|0.75', 1: 'zeta|0.5', 2: 'dh|0.5'}
SEED0 = 20261025
NS = [2 ** k for k in range(7, 11 if args.quick else 16)]
REPS = 6 if args.quick else 200
NFIT = 1024
GRID = 2 ** 14 if args.quick else 2 ** 17


def hval(kind, s, nf=0.3, M=28):
    return Z.zeta(s, nf=nf, M=M) if kind == 'zeta' else Z.dh(s, nf=nf, M=M)


def hderiv(kind, s):
    return Z.zeta(s, deriv=True) if kind == 'zeta' else Z.dh(s, deriv=True)


# ----------------------------------------------------------------------------- geometry
def dens_arc(phi, wz):
    ws = -np.conj(wz)
    e = np.exp(1j * phi)
    return (1 - abs(wz) ** 2) / (2 * np.pi * np.abs(e - wz) ** 2) - (1 - abs(ws) ** 2) / (2 * np.pi * np.abs(e - ws) ** 2)


def dxG(a, w):
    return np.real(-np.conj(a) / (1 - np.conj(a) * w) - 1 / (w - a))


def dens_dia(y, wz):
    ws = -np.conj(wz)
    w = 1j * y
    return (dxG(wz, w) - dxG(ws, w)) / (2 * np.pi)


def GD(a, w):
    return np.log(np.abs((1 - np.conj(a) * w) / (w - a)))


def GH(wz, wr):
    return GD(wz, wr) - GD(-np.conj(wz), wr)


class Geo:
    def __init__(self, theta):
        self.theta = theta
        self.c = theta + 1j * T0
        self.z = self.c + 0.5 * R
        self.wz = (self.z - self.c) / R

    # piece 0: arc, variable phi in [-pi/2, pi/2]; piece 1: diameter, variable y in [-1, 1]
    def pts(self, piece, x):
        if piece == 0:
            return self.c + R * np.exp(1j * x)
        return self.c + 1j * R * x

    def dens(self, piece, x):
        return dens_arc(x, self.wz) if piece == 0 else dens_dia(x, self.wz)

    LO = (-np.pi / 2, -1.0)
    HI = (np.pi / 2, 1.0)


def gl_integral(f, lo, hi, n):
    x, w = leggauss(n)
    xm = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    return 0.5 * (hi - lo) * np.sum(w * f(xm))


# ----------------------------------------------------------------------------- zeros
def newton_zero(kind, s0, it=40):
    s = complex(s0)
    for _ in range(it):
        v, d = hderiv(kind, np.array([s]))
        ds = v[0] / d[0]
        s -= ds
        if abs(ds) < 1e-15:
            break
    return s


def find_diameter_zeros(kind, theta):
    t = np.arange(T0 - R, T0 + R + 1e-9, 1e-3)
    v = np.abs(hval(kind, theta + 1j * t))
    idx = [i for i in range(1, len(t) - 1) if v[i] <= v[i - 1] and v[i] <= v[i + 1] and v[i] < 0.05]
    zs = []
    for i in idx:
        s = newton_zero(kind, theta + 1j * t[i])
        fv = abs(hval(kind, np.array([s]))[0])
        if abs(s.real - theta) < 1e-9 and fv < 1e-12 and T0 - R < s.imag < T0 + R and all(abs(s - q) > 1e-8 for q in zs):
            zs.append(s)
    return zs, [float(t[i]) for i in idx], float(v.min())


def make_case(case):
    kind, theta = KIND[case]
    g = Geo(theta)
    zs, cand, vmin = find_diameter_zeros(kind, theta)
    if case == 1:
        # the known zeta zero by Newton from the first seed, and the cross-check with the table of zeros
        zs = [newton_zero('zeta', 0.5 + 84.7354929805j)]
    ins = []
    if kind == 'dh':
        with open(finddata('cf_dhzeros.json')) as f:
            d = json.load(f)
        zz = [r for r in d['T2']['zeros'] if r['tag'] == 'z85'][0]
        ins = [complex(float(zz['center_re']), float(zz['center_im']))]
    return g, kind, zs, ins


def u_of(kind, g, zs, piece, x, nf=0.3, M=28):
    s = g.pts(piece, x)
    hv = hval(kind, s, nf, M)
    u = np.log(np.abs(hv))
    for r in zs:
        u = u - np.log(np.abs(s - r))
    return u


def rhs(kind, g, zs, ins, nf=0.3, M=28):
    hz = hval(kind, np.array([g.z]), nf, M)[0]
    v = math.log(abs(hz)) - sum(math.log(abs(g.z - r)) for r in zs)
    gh = [float(GH(g.wz, (r - g.c) / R)) for r in ins]
    return v + sum(gh), gh


def reference(kind, g, zs, nf=0.3, M=28):
    prev = None
    n = 64
    while True:
        tot = 0.0
        parts = []
        for p in (0, 1):
            val = gl_integral(lambda x: u_of(kind, g, zs, p, x, nf, M) * g.dens(p, x), g.LO[p], g.HI[p], n)
            parts.append(val); tot += val
        if prev is not None and abs(tot - prev) < 1e-14:
            return tot, n, abs(tot - prev), parts
        if n >= 4096:
            return tot, n, abs(tot - prev) if prev is not None else float('nan'), parts
        prev = tot
        n *= 2


# ----------------------------------------------------------------------------- configuration checks
def winding(kind, g, zs, m=40000):
    # counter clockwise: the arc from phi = -pi/2 to pi/2, then the diameter from y = 1 to y = -1
    out = []
    x = np.linspace(-np.pi / 2, np.pi / 2, m + 1)
    out.append(g.pts(0, x))
    y = np.linspace(1, -1, m + 1)
    out.append(g.pts(1, y))
    s = np.concatenate(out)
    q = hval(kind, s)
    for r in zs:
        q = q / (s - r)
    ang = np.angle(q[1:] / q[:-1])  # the junction at the corner is a regular step; the path is closed (the last point is the first)
    tot = float(ang.sum())
    return tot / (2 * np.pi), float(np.abs(ang).max()), float(np.abs(q).min())


def identity_selftest(n):
    g = Geo(0.5)
    r1 = g.c + 0.2j
    r2 = g.c + 0.3 + 0.1j
    h = lambda s: (s - r1) * (s - r2)
    tot = 0.0
    for p in (0, 1):
        def f(x, p=p):
            s = g.pts(p, x)
            return (np.log(np.abs(h(s))) - np.log(np.abs(s - r1))) * g.dens(p, x)
        tot += gl_integral(f, g.LO[p], g.HI[p], n)
    rs = math.log(abs(h(g.z))) - math.log(abs(g.z - r1)) + float(GH(g.wz, (r2 - g.c) / R))
    return tot, rs


# ----------------------------------------------------------------------------- estimators
def prep_sampler(g):
    # per piece: grid, cumulative mass (Simpson-free: trapezoid on a fine grid of GRID points, a smooth integrand), and the mass
    S = []
    for p in (0, 1):
        x = np.linspace(g.LO[p], g.HI[p], GRID + 1)
        d = g.dens(p, x)
        c = np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(x))])
        S.append((x, c, c[-1]))
    return S


def sample_omega(g, S, rng, n):
    m0, m1 = S[0][2], S[1][2]
    piece = rng.random(n) < m0 / (m0 + m1)
    out = np.empty(n, dtype=complex)
    for p, mask in ((0, piece), (1, ~piece)):
        k = int(mask.sum())
        if k:
            x, c, m = S[p]
            q = rng.random(k) * m
            xs = np.interp(q, c, x)
            out[mask] = g.pts(p, xs)
    return out


def trap_est(kind, g, zs, rng, n_per, use_phi):
    tot = 0.0
    for p in (0, 1):
        lo, hi = g.LO[p], g.HI[p]
        # nodes on [-1, 1] (v), sorted, end points added; x = (hi) * Phi(v) for the change of variable, x = lo + (hi - lo) v' for the plain rule
        if use_phi:
            v = np.sort(rng.random(n_per) * 2 - 1)
            v = np.concatenate([[-1.0], v, [1.0]])
            Phi = v + np.sin(np.pi * v) / np.pi
            dPhi = 1 + np.cos(np.pi * v)
            x = hi * Phi
            gx = u_of(kind, g, zs, p, x) * g.dens(p, x) * hi * dPhi
            gx[0] = gx[-1] = 0.0
            tot += float(np.sum(0.5 * (gx[1:] + gx[:-1]) * np.diff(v)))
        else:
            x = np.sort(lo + (hi - lo) * rng.random(n_per))
            x = np.concatenate([[lo], x, [hi]])
            gx = u_of(kind, g, zs, p, x) * g.dens(p, x)
            gx[0] = gx[-1] = 0.0
            tot += float(np.sum(0.5 * (gx[1:] + gx[:-1]) * np.diff(x)))
    return tot


def one_rep(task):
    case, rep, ns = task
    g, kind, zs, ins = CTX[case]
    rng = np.random.default_rng(np.random.SeedSequence([SEED0, case, rep]))
    S = SAMP[case]
    res = []
    for N in ns:
        s = sample_omega(g, S, rng, N)
        hv = hval(kind, s)
        u = np.log(np.abs(hv))
        for r in zs:
            u = u - np.log(np.abs(s - r))
        plain = float(u.mean())
        tr = trap_est(kind, g, zs, rng, N // 2, False)
        tp = trap_est(kind, g, zs, rng, N // 2, True)
        res.append((plain, tr, tp))
    return case, rep, res


CTX = {}
SAMP = {}


def init():
    for case in KIND:
        CTX[case] = make_case(case)
        SAMP[case] = prep_sampler(CTX[case][0])


def fit(Ns, rm):
    Ns = np.array(Ns, float); rm = np.array(rm, float)
    sel = Ns >= NFIT
    X = np.log(Ns[sel]); Y = np.log(rm[sel])
    A = np.vstack([X, np.ones_like(X)]).T
    coef, res, *_ = np.linalg.lstsq(A, Y, rcond=None)
    n = len(X)
    r = Y - A @ coef
    s2 = float(r @ r) / max(n - 2, 1)
    se = math.sqrt(s2 / float(((X - X.mean()) ** 2).sum())) if n > 2 else float('nan')
    return float(coef[0]), se


def main():
    t_start = time.time()
    log('bsy_halfdisk.py', 'quick' if args.quick else 'full', 'Ns', NS, 'reps', REPS, 'procs', args.procs)
    init()
    R_ = {'quick': args.quick, 'Ns': NS, 'reps': REPS, 'seed': [SEED0, 'case', 'rep'], 'R': R, 't0': T0, 'nfit': NFIT, 'cases': {}}
    # ---- configuration
    st = {}
    for case in KIND:
        g, kind, zs, ins = CTX[case]
        nm = NAME[case]
        log('---- case', nm)
        c = {'theta': g.theta, 'z': [g.z.real, g.z.imag], 'diameter_zeros': [[r.real, r.imag] for r in zs], 'zeros_in_H': [[r.real, r.imag] for r in ins]}
        mass = [gl_integral(lambda x, p=p: g.dens(p, x), g.LO[p], g.HI[p], 400) for p in (0, 1)]
        c['mass_arc'], c['mass_diameter'] = mass
        c['mass_total_minus_1'] = sum(mass) - 1
        log('harmonic measure: arc %.15f diameter %.15f total-1 %.2e' % (mass[0], mass[1], sum(mass) - 1))
        w, mxinc, qmin = winding(kind, g, zs)
        c['winding'] = w; c['max_arg_increment'] = mxinc; c['min_abs_q_on_boundary'] = qmin
        log('argument principle: winding of h/prod(s-rho_j) on dH = %.6f (expected %d), max |arg increment| %.3f, min|q| on dH %.3e' % (w, len(ins), mxinc, qmin))
        if case == 1:
            tab = [float(l.split()[0]) for l in open(finddata('zeros_1700.txt')) if l.strip()]
            near = min(tab, key=lambda t: abs(t - zs[0].imag))
            c['zero_vs_table'] = abs(near - zs[0].imag)
            log('zeta zero on the diameter: %.15f + %.15f i, table value %.15f, difference %.2e' % (zs[0].real, zs[0].imag, near, abs(near - zs[0].imag)))
        elif kind == 'dh':
            log('DH diameter scan: %d candidate minima, min |f| on the diameter %.3e, zeros on the diameter found: %d' % (len(find_diameter_zeros(kind, g.theta)[1]), find_diameter_zeros(kind, g.theta)[2], len(zs)))
            log('DH zero in H: %.16f + %.16f i, inside: %s' % (ins[0].real, ins[0].imag, abs(ins[0] - g.c) < R and ins[0].real > g.theta))
        else:
            log('zeta at theta = 3/4: no zero on the diameter (min |zeta| on it %.3e)' % find_diameter_zeros(kind, g.theta)[2])
        # target and reference at two evaluator settings
        rs, gh = rhs(kind, g, zs, ins)
        rs2, _ = rhs(kind, g, zs, ins, 0.45, 36)
        ref, nn, dd, parts = reference(kind, g, zs)
        ref2, nn2, dd2, _ = reference(kind, g, zs, 0.45, 36)
        c['target'] = rs; c['green_terms'] = gh; c['reference'] = ref; c['ref_nodes'] = nn; c['ref_last_diff'] = dd
        c['ref_minus_RS'] = ref - rs
        c['prec'] = {'RS_second_evaluator': rs2 - rs, 'ref_second_evaluator': ref2 - ref, 'ref_second_nodes': nn2}
        log('target RS = %.15f (Green terms %s); reference GL = %.15f with %d nodes per piece (last doubling changed it by %.2e); reference - RS = %.3e' % (rs, ['%.10f' % x for x in gh], ref, nn, dd, ref - rs))
        log('precision: RS (nf .45, M 36) - RS = %.2e; reference (second evaluator) - reference = %.2e' % (rs2 - rs, ref2 - ref))
        # mpmath check at 50 points
        rng = np.random.default_rng(1)
        pp = [(rng.integers(0, 2), rng.uniform(-0.99, 0.99)) for _ in range(50)]
        import mpmath as mpm
        mpm.mp.dps = 30
        mx = 0.0
        for p, x in pp:
            s = complex(g.pts(p, np.array([x * g.HI[p]]))[0])
            a = hval(kind, np.array([s]))[0]
            b = complex(mpm.zeta(mpm.mpmathify(s))) if kind == 'zeta' else Z.mp_h('dh', g.theta, s, 30)
            mx = max(mx, abs(a - b) / abs(b))
        c['prec']['evaluator_vs_mpmath_max_rel'] = mx
        log('evaluator against mpmath (30 digits) at 50 boundary points: max relative difference %.2e' % mx)
        R_['cases'][nm] = c
        st[case] = c
    s1, s2 = identity_selftest(200)
    s1b, _ = identity_selftest(400)
    R_['selftest'] = {'lhs': s1, 'rhs': s2, 'diff': s1 - s2, 'lhs_n400_minus_n200': s1b - s1}
    log('identity self test (polynomial with one zero on the diameter and one in H): LHS %.15f RHS %.15f difference %.2e; 400 nodes minus 200 nodes %.2e' % (s1, s2, s1 - s2, s1b - s1))
    # ---- replications
    tasks = [(case, rep, NS) for case in KIND for rep in range(REPS)]
    t1 = time.time()
    if args.procs > 1:
        with mp.get_context('fork').Pool(args.procs) as pool:
            out = pool.map(one_rep, tasks, chunksize=1)
    else:
        out = [one_rep(t) for t in tasks]
    log('replications: %.1f s' % (time.time() - t1))
    for case in KIND:
        c = st[case]
        arr = np.array([r[2] for r in sorted([o for o in out if o[0] == case], key=lambda o: o[1])])  # reps x N x 3
        err = arr - c['target']
        rm = np.sqrt((err ** 2).mean(0))   # N x 3
        me = err.mean(0)
        sl = [fit(NS, rm[:, k]) for k in range(3)]
        c['rmse'] = {nm: rm[:, k].tolist() for k, nm in enumerate(('plain', 'trap', 'trap_phi'))}
        c['mean_error'] = {nm: me[:, k].tolist() for k, nm in enumerate(('plain', 'trap', 'trap_phi'))}
        c['slope'] = {nm: sl[k][0] for k, nm in enumerate(('plain', 'trap', 'trap_phi'))}
        c['slope_se'] = {nm: sl[k][1] for k, nm in enumerate(('plain', 'trap', 'trap_phi'))}
        i = NS.index(32768) if 32768 in NS else len(NS) - 1
        c['N_ratio'] = NS[i]
        c['ratio_plain_over_trap'] = rm[i, 0] / rm[i, 1]
        c['ratio_plain_over_trap_phi'] = rm[i, 0] / rm[i, 2]
        log('==== %s: RMSE at N = %d  plain %.3e (slope %.2f +- %.2f)  trap %.3e (%.2f +- %.2f)  trap+Phi %.3e (%.2f +- %.2f);  ratios plain/trap %.3g, plain/Phi %.3g' % (
            NAME[case], NS[i], rm[i, 0], sl[0][0], sl[0][1], rm[i, 1], sl[1][0], sl[1][1], rm[i, 2], sl[2][0], sl[2][1], c['ratio_plain_over_trap'], c['ratio_plain_over_trap_phi']))
        # the Lemma bias: leading term summed over the two pieces, g' at the ends by a 4 point one sided formula (step 1e-3)
        g_, kind, zs, ins = CTX[case]
        sl_end = []
        for p in (0, 1):
            f = lambda x, p=p: u_of(kind, g_, zs, p, x) * g_.dens(p, x)
            lo, hi = g_.LO[p], g_.HI[p]
            e = 1e-3 * (hi - lo) / 2
            xs = hi - e * np.arange(0, 4)
            vr = f(xs)
            dhi = (-11 * vr[0] + 18 * vr[1] - 9 * vr[2] + 2 * vr[3]) / (6 * e)  # derivative at hi from the left
            xs = lo + e * np.arange(0, 4)
            vl = f(xs)
            dlo = (-11 * vl[0] + 18 * vl[1] - 9 * vl[2] + 2 * vl[3]) / (6 * e) * 1.0  # forward difference at lo (derivative)
            sl_end.append((lo, hi, float(dlo), float(-dhi * -1.0) if False else float((11 * vr[0] - 18 * vr[1] + 9 * vr[2] - 2 * vr[3]) / (6 * e))))
        bl = []
        for N in NS:
            n = N // 2
            bl.append(float(sum((hi - lo) ** 2 * (dh - dl) / (2 * (n + 1) * (n + 2)) for lo, hi, dl, dh in sl_end)))
        c['bias_lemma'] = bl
        c['end_slopes'] = sl_end
        c['bias_ratio_trap'] = [float(me[k, 1] / bl[k]) if bl[k] != 0 else None for k in range(len(NS))]
        log('   Lemma bias (leading term) / measured mean error of the sorted trapezoid:', ['%.3f' % (bl[k] / me[k, 1]) for k in range(len(NS))])
        for k, N in enumerate(NS):
            log('   N %6d  mean error plain %+.3e trap %+.3e trap+Phi %+.3e   RMSE plain %.3e trap %.3e trap+Phi %.3e' % (N, me[k, 0], me[k, 1], me[k, 2], rm[k, 0], rm[k, 1], rm[k, 2]))
    R_['seconds'] = time.time() - t_start
    with open(OUT, 'w') as f:
        json.dump(R_, f)
    log('wrote', OUT, 'in %.1f s' % R_['seconds'])


if __name__ == '__main__':
    main()
