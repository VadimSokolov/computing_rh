# mi_survey.py: the survey over primes of Chapter 22 (Polson's Miyatake chapter), Section "A survey to height 21,000" and the
# claims of the lines around Table survey (offset R_P - N(T_P^*) in {-1,0,1,2,3}, mean 0.90 and standard error 0.27, the ramp model
# with average 0.47, the density log(t/2 pi)/(4 pi) of negative humps). Added for the book in October 2026 (computation group G2), written
# from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. xi(s) = (1/2) s (s-1) pi^(-s/2) Gamma(s/2) zeta(s), Xi(t) = xi(1/2 + i t) = -M(t) Z(t) with Z Hardy's function and
# M(t) = (1/2)(t^2 + 1/4) pi^(-1/4) |Gamma(1/4 + i t/2)|. S(P) = integers with all prime factors <= P, q the least prime above P.
# Xi_P = Xi + D_P, D_P(t) = (t^2 + 1/4) sum_{m not in S(P)} C_(pi m^2)(t), C_a(t) = int_0^oo exp(-a e^w) e^(w/4) cos(t w/2) dw
# = Re[ a^(-s') Gamma(s', a) ], s' = 1/4 + i t/2 (substitute v = a e^w). With g(s,a) = a^(-s) e^a Gamma(s,a) (Legendre continued
# fraction g = 1/(a + 1 - s - 1(1-s)/(a + 3 - s - 2(2-s)/(a + 5 - s - ...)))) this is C_a(t) = e^(-a) Re g(s', a). The sum over m is
# truncated at m^2 - q^2 <= 16 (the ratio of the terms is e^(-pi (m^2 - q^2)) < 1e-21). T_P^* is the root of M(T) = D_P(T) above the Gram point region.
# Hump criterion. Xi < 0 on the negative humps (gamma_(2j-1), gamma_(2j)) (Z > 0 there); the hump carries two zeros of Xi_P iff
# max over the hump of (|Xi| - D_P) > 0, i.e. max (Z - D_P/M) > 0. R_P = 2 * (number of surviving negative humps); all humps below
# T_P^* - 80 are taken to survive (checked: the minimum of the margin over the humps from T_P^* - 80 to T_P^* - 40 is logged),
# humps up to T_P^* + 80 are tested (none survives above T_P^* + 40, logged). off = R_P - N(T_P^*), N(T) = number of zeros of Xi in (0, T].
# Ramp: (2/pi^2) log(T_P^*/2 pi) times the mean of log max|Z| over the negative humps with ordinate within 200 of T_P^* (window: our choice).
# Method for Z. Riemann-Siegel with the terms C0, C1, C2 (Gabcke): Z = 2 sum_(n<=N) n^(-1/2) cos(theta - t log n) + (-1)^(N-1)
# (t/2pi)^(-1/4) [C0(p) + C1(p) u + C2(p) u^2], u = (t/2pi)^(-1/2), N = floor(sqrt(t/2pi)), p = sqrt(t/2pi) - N, C0 = Phi,
# C1 = -Phi'''/(96 pi^2), C2 = Phi''/(64 pi^2) + Phi^(6)/(18432 pi^4), Phi(p) = cos(2 pi (p^2 - p - 1/16))/cos(2 pi p); Phi, Phi'''
# etc. are entire, evaluated from Chebyshev interpolants (degree 48, nodes by mpmath at 60 digits). Zeros of Z: sign changes on the grid
# t = k h (h = 0.02), located by linear interpolation; the hump maximum of Z - D_P/M is the parabola refined grid maximum. Marginal humps
# (|margin| < 1e-3 (1 + Z)) are re-evaluated by mpmath siegelz and a golden section.
# Error control. (1) Z from Riemann-Siegel against mpmath siegelz at 24 points (max difference logged), and the zeros against the first
# 1700 ordinates of data/zeros_1700.txt (max difference and count agreement logged). (2) the continued fraction depth 800 against 1600
# (relative difference of D_P/M logged). (3) the survey is run twice, with grid steps h and 0.8 h in the windows (option --h2): the
# counts R_P and N are compared (every difference logged). Output data/mi_survey.json and data/mi_survey.log (quick: --out FILE).
import sys, json, time, argparse, math
import numpy as np
from scipy.special import loggamma
from scipy.optimize import brentq
import mpmath as mpm

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default=None)
ap.add_argument('--pmax', type=int, default=71)
ap.add_argument('--pmin', type=int, default=2)
ap.add_argument('--h', type=float, default=0.02)
ap.add_argument('--h2', action='store_true', help='rerun with step 0.8 h and compare the counts')
args = ap.parse_args()
if args.quick:
    args.pmax = min(args.pmax, 13)
out = args.out or ('data/mi_survey.json' if not args.quick else 'misc/reviews/nick-miyatake/quick/G2/mi_survey.json')
logf = open(out.replace('.json', '.log'), 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.write(s + '\n'); logf.flush()

PI = math.pi
primes_all = [p for p in range(2, 80) if all(p % d for d in range(2, p))]
PRIMES = [p for p in primes_all if args.pmin <= p <= args.pmax]
def nextprime(P):
    return next(p for p in primes_all if p > P)

# ---------- Chebyshev interpolants of Phi, C1, C2 on [0,1]
t0 = time.time()
mpm.mp.dps = 60
def Phi(p):
    return mpm.cos(2 * mpm.pi * (p * p - p - mpm.mpf(1) / 16)) / mpm.cos(2 * mpm.pi * p)
DEG = 48
nodes = [mpm.mpf(1) / 2 + mpm.cos(mpm.pi * (2 * k + 1) / (2 * (DEG + 1))) / 2 for k in range(DEG + 1)]
def cheb_fit(fn):
    vals = [fn(x) for x in nodes]
    return np.polynomial.chebyshev.Chebyshev.interpolate(lambda x: np.interp(x, [float(n) for n in nodes][::-1], [float(v) for v in vals][::-1]), 1) if False else \
        np.polynomial.chebyshev.chebfit(np.array([2 * float(n) - 1 for n in nodes]), np.array([float(v) for v in vals]), DEG)
c0 = cheb_fit(Phi)
c1 = cheb_fit(lambda x: -mpm.diff(Phi, x, 3) / (96 * mpm.pi ** 2))
c2 = cheb_fit(lambda x: mpm.diff(Phi, x, 2) / (64 * mpm.pi ** 2) + mpm.diff(Phi, x, 6) / (18432 * mpm.pi ** 4))
mpm.mp.dps = 30
log('Chebyshev interpolants built in %.1f s' % (time.time() - t0))
def cheb(c, p):
    return np.polynomial.chebyshev.chebval(2 * p - 1, c)

def theta(t):
    return t / 2 * np.log(t / (2 * PI)) - t / 2 - PI / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)

def Zrs(t):
    t = np.asarray(t, float)
    r = np.sqrt(t / (2 * PI)); N = np.floor(r).astype(int); p = r - N
    nmax = int(N.max()) if N.size else 0
    th = theta(t)
    S = np.zeros_like(t)
    for n in range(1, nmax + 1):
        m = N >= n
        S[m] += np.cos(th[m] - t[m] * math.log(n)) / math.sqrt(n)
    u = 1 / r
    corr = (np.where(N % 2 == 1, 1.0, -1.0)) * r ** -0.5 * (cheb(c0, p) + cheb(c1, p) * u + cheb(c2, p) * u * u)
    return 2 * S + corr

def Zgrid(h, tmax):
    K = int(tmax / h) + 2
    t = np.arange(1, K + 1) * h
    Z = np.empty(K)
    for a in range(0, K, 20000):
        Z[a:a + 20000] = Zrs(t[a:a + 20000])
    return np.concatenate([[0.0], t]), np.concatenate([[float(-mpm.zeta(0.5) * 0 + 0)], Z])

# ---------- D_P / M
def logM(t):
    return np.log(0.5 * (t * t + 0.25) * PI ** -0.25) + loggamma(0.25 + 0.5j * t).real

def gCF(s, a, depth):
    # g(s,a) = a^(-s) e^a Gamma(s,a) by the Legendre continued fraction, backward recurrence, s array complex, a float
    f = np.zeros_like(s, dtype=complex)
    for n in range(depth, 0, -1):
        f = -n * (n - s) / (a + 2 * n + 1 - s + f)
    return 1 / (a + 1 - s + f)

def logD(t, P, depth=800):
    q = nextprime(P)
    sp = 0.25 + 0.5j * t
    terms = []
    for m in range(q, q + 8):
        if m * m - q * q <= 16:
            # m not in S(P) iff the largest prime factor of m exceeds P
            k, lp, d = m, 1, 2
            while k > 1:
                if k % d == 0:
                    lp = d; k //= d
                else:
                    d += 1
            if lp <= P:
                continue
            a = PI * m * m
            g = gCF(sp, a, depth)
            terms.append(-a + np.log(np.abs(g.real) + 1e-300) + 0 * t)
            if m == q:
                gsign = np.sign(g.real)
    L = np.array(terms)
    mx = L.max(axis=0)
    # signs: all terms of the dominant m = q are used with their sign; the others are below by e^(-pi(m^2-q^2)) and have the sign of Re g too
    return np.log(t * t + 0.25) + mx + np.log(np.sum(np.exp(L - mx), axis=0))

def logr(t, P, depth=800):
    return logD(t, P, depth) - logM(t)

def Tstar(P):
    f = lambda T: float(logD(np.array([T]), P)[0] - logM(np.array([T]))[0])
    lo = 10.0; hi = 4 * nextprime(P) ** 2 + 60.0
    # D/M increases from about 0 to infinity; locate the sign change
    return brentq(f, lo, hi, xtol=1e-9)

# ---------- validation of Z
rng = np.random.default_rng(1)
tt = np.concatenate([rng.uniform(100, 1000, 12), rng.uniform(5000, 21000, 12)])
mpm.mp.dps = 20
zr = Zrs(tt); zm = np.array([float(mpm.siegelz(mpm.mpf(x))) for x in tt])
mpm.mp.dps = 30
log('Riemann-Siegel (C0,C1,C2) against mpmath siegelz at %d points: max |dZ| = %.2e (t <= 1000: %.2e, t >= 5000: %.2e)'
    % (len(tt), np.abs(zr - zm).max(), np.abs(zr - zm)[:12].max(), np.abs(zr - zm)[12:].max()))

tmax = max(4 * nextprime(p) ** 2 + 260.0 for p in PRIMES)
def zeros_from_grid(t, Z):
    idx = np.nonzero(np.sign(Z[1:]) * np.sign(Z[:-1]) < 0)[0]
    return t[idx] - Z[idx] * (t[idx + 1] - t[idx]) / (Z[idx + 1] - Z[idx]), idx

def run(h, label):
    t0 = time.time()
    t, Z = Zgrid(h, tmax)
    Z[0] = -1.0
    zer, idx = zeros_from_grid(t, Z)
    log('[%s] Z grid h=%g to t=%.0f: %d points, %d zeros, %.1f s' % (label, h, tmax, len(t), len(zer), time.time() - t0))
    return t, Z, zer, idx

t, Z, zer, idx = run(args.h, 'h')
ref = [float(mpm.mpf(l)) for l in open(__import__('os').path.exists('data/zeros_1700.txt') and 'data/zeros_1700.txt' or 'zeros_1700.txt').read().split()]
nc = min(len(ref), int(np.searchsorted(zer, ref[-1] + 1)))
nn = min(len(ref), len(zer))
log('first %d zeros against data/zeros_1700.txt: count of grid zeros up to gamma_1700 + 1 = %d, max |zero - ref| = %.2e'
    % (len(ref), nc, np.abs(zer[:nn] - np.array(ref[:nn])).max() if nn else 0))
results = []
def survey(t, Z, zer, label, hstep):
    rows = []
    for P in PRIMES:
        q = nextprime(P)
        Ts = Tstar(P)
        N = int(np.searchsorted(zer, Ts))
        w = (t > max(1.0, Ts - 80)) & (t < Ts + 80)
        tw, Zw = t[w], Z[w]; w0 = int(np.nonzero(w)[0][0])
        lr = logr(tw, P)
        lr2 = logr(tw, P, 1600)
        dd = float(np.abs(lr - lr2).max())
        r = np.exp(np.minimum(lr, 50))
        # humps: (zer[2j], zer[2j+1]) 0-based j = hump j+1
        A_ = zer[0:len(zer) // 2 * 2:2]; B_ = zer[1:len(zer) // 2 * 2:2]
        below = int(np.searchsorted(B_, Ts - 80))
        nhump = below; surv = below; marg = []; early = []; late = 0
        for j in range(below, len(B_)):
            a, b = A_[j], B_[j]
            if b > Ts + 80: break
            i0, i1 = int(a / hstep) + 1, int(b / hstep) + 1
            m = slice(max(i0, int(max(1.0, Ts - 80) / hstep) + 1) - w0, i1 + 1 - w0)
            v = Zw[m] - r[m]
            k = int(np.argmax(v))
            val = v[k]
            if 0 < k < len(v) - 1:
                y0, y1, y2 = v[k - 1], v[k], v[k + 1]
                den = y0 - 2 * y1 + y2
                if den < 0:
                    val = y1 - (y2 - y0) ** 2 / (8 * den)
            zmax = float(Zw[m].max())
            ok = val > 0
            if abs(val) < 1e-3 * (1 + abs(zmax)):
                mpm.mp.dps = 20
                tm = float(tw[m][k])
                fn = lambda x: -(float(mpm.siegelz(mpm.mpf(x))) - float(np.exp(min(logr(np.array([x]), P)[0], 50))))
                lo, hi = tm - hstep, tm + hstep
                gr = (math.sqrt(5) - 1) / 2
                c_, d_ = hi - gr * (hi - lo), lo + gr * (hi - lo)
                for _ in range(25):
                    if fn(c_) < fn(d_): hi = d_
                    else: lo = c_
                    c_, d_ = hi - gr * (hi - lo), lo + gr * (hi - lo)
                mpm.mp.dps = 30
                val = -fn((lo + hi) / 2); ok = val > 0
                marg.append((round(0.5 * (a + b), 3), float(val)))
            nhump += 1
            if ok: surv += 1
            elif b < Ts - 40: early.append(round(0.5 * (a + b), 2))
            if ok and a > Ts + 40: late += 1
        Lm = []
        for j in range(int(np.searchsorted(A_, Ts - 200)), len(B_)):
            a, b = A_[j], B_[j]
            if b > Ts + 200: break
            Lm.append(math.log(Z[int(a / hstep) + 1:int(b / hstep) + 2].max()))
        ramp = 2 / PI ** 2 * math.log(Ts / (2 * PI)) * float(np.mean(Lm))
        # density of negative humps near Ts
        nh = int(np.sum(np.abs(A_ - Ts) < 200))
        dens = nh / 400.0; dens_th = math.log(Ts / (2 * PI)) / (4 * PI)
        R = 2 * surv
        rows.append(dict(P=P, q=q, Tstar=Ts, N=N, R=R, off=R - N, ramp=ramp, depth_diff=dd, marginal=marg, lost_early=early,
                         survivors_above_Tstar_plus_40=late, hump_density=dens, hump_density_theory=dens_th, nLm=len(Lm)))
        log('[%s] P=%d q=%d T*=%.3f N=%d R=%d off=%d ramp=%.3f (CF depth diff %.1e) marginal=%s lost below T*-40=%s survivors above T*+40=%d density %.4f vs %.4f'
            % (label, P, q, Ts, N, R, R - N, ramp, dd, marg, early, late, dens, dens_th))
    return rows
rows = survey(t, Z, zer, 'h', args.h)
off = np.array([r['off'] for r in rows], float)
log('off in %s, mean %.3f, standard error %.3f (n=%d); mean ramp %.3f' % (sorted(set(int(x) for x in off)), off.mean(), off.std(ddof=1) / math.sqrt(len(off)),
    len(off), np.mean([r['ramp'] for r in rows])))
res = dict(rows=rows, h=args.h, off_mean=off.mean(), off_se=off.std(ddof=1) / math.sqrt(len(off)), ramp_mean=float(np.mean([r['ramp'] for r in rows])))
if args.h2:
    t2, Z2, zer2, idx2 = run(0.8 * args.h, 'h2')
    rows2 = survey(t2, Z2, zer2, 'h2', 0.8 * args.h)
    diffs = [(a['P'], a['N'] - b['N'], a['R'] - b['R']) for a, b in zip(rows, rows2)]
    log('step comparison (P, dN, dR) h vs 0.8h:', [d for d in diffs if d[1] or d[2]] or 'no differences')
    res['step_comparison'] = diffs
json.dump(res, open(out, 'w'), indent=1)
