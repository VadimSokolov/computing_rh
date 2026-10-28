# mi_depart.py: the departed zeros of Xi_P, Chapter 22 (Polson's Miyatake chapter), Section "Where the lost zeros go": the departure law
# |y_P(x)| = 2 log(D_P(x)/M(x))/log(x/2 pi) with D_P and with the limit c_P, Table dep (the six departed zeros of Xi_2 and Xi_3), the
# argument principle counts N_P(X) - N(X) against y_P(X)/4 (P = 2 at X = 51.5, 60, 68, 80; P = 3 at X = 120, 140, 150), the traced chain of
# Xi_2 from 41.19 + 1.46 i to 601.90 + 191.23 i with 189 zeros of real part below 600 against 191.4 predicted, and the density
# 1/(4 pi) log(x/2 pi) + pi/(16 log(x/2 pi)). Added for the book in October 2026 (computation group G2), written from the mathematics
# (the collaborator's scripts were not read, run or copied).
# Conventions. xi(s) = (1/2) s (s-1) pi^(-s/2) Gamma(s/2) zeta(s), Xi(t) = xi(1/2 + i t); Xi_P = Xi + D_P, D_P(t) = (t^2 + 1/4) sum over
# m not in S(P) of C_(pi m^2)(t), C_a(t) = int_0^oo exp(-a e^w) e^(w/4) cos(t w/2) dw = (1/2)[a^(-s') Gamma(s', a) + a^(-s'') Gamma(s'', a)],
# s' = 1/4 + i t/2, s'' = 1/4 - i t/2 (valid for complex t by analytic continuation, v = a e^w). M(t) = (1/2)(t^2 + 1/4) pi^(-1/4) |Gamma(1/4 + i t/2)|.
# c_P = lim D_P(t) = 4 sum over m not in S(P) of e^(-pi m^2)(pi m^2 - 1/4) (integration by parts of C_a at infinity). The sum over m is
# cut at m^2 - q^2 <= 30 (the ratio of the terms is e^(-pi (m^2 - q^2)) < 1e-40); the cut 120 is used for the check.
# Check of the evaluator: Xi_P from Proposition incgamma (sum over n in S(P)) at 60 digits against Xi + D_P at several complex points (logged).
# Zeros: complex Newton (mpmath findroot, secant) from the departure law. Argument principle: N_P(X) = number of zeros of Xi_P in
# 0 < Re t < X (both half planes), by the rectangle [0, X] x [-Y, Y]; since Xi_P(conj t) = conj Xi_P(t), the change of argument is twice that
# along (X, 0) -> (X, Y) -> (0, Y) -> (0, 0); the sides are first cut into steps of length 0.4 (the argument turns by less than 1.5 rad per unit on every side), then points are added by bisection until every step changes the argument by less than 0.5 rad, the
# integer is accepted if the total is within 0.02 of it, and the count is repeated with a second Y and a second precision (the differences are
# logged). On the top side |D_P/Xi| is logged. N(X) from data/zeros_1700.txt. The chain of Xi_2 is traced by continuation in the real part
# (prediction by the previous spacing, Newton, acceptance only if the real part increases and the step is within 0.5 to 2 times the
# predicted local spacing 4 pi/log(x/2 pi), else smaller predictions); completeness is checked against the argument principle count
# of the region (option --chain-check; the run to 600 is large and is left for Hopper when the laptop cannot do it in 30 s).
# Output data/mi_depart.json and data/mi_depart.log (quick: --out FILE). Parts: --part table|counts|chain|all (default all), --cases for counts.
import sys, json, time, argparse, math
import mpmath as mp

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default=None)
ap.add_argument('--part', default='all')
ap.add_argument('--cases', default='2:51.5', help='comma list P:X for the counts')
ap.add_argument('--chain-to', type=float, default=100.0)
ap.add_argument('--single', action='store_true', help='counts: only the first (Y, precision)')
ap.add_argument('--skip-check', action='store_true')
ap.add_argument('--seg', default=None, help='chain segment lo:hi traced from a seed at lo (a few zeros back to lo - 12)')
ap.add_argument('--seed', default=None, help='x:y of the seed of a segment')
ap.add_argument('--tmax', type=float, default=22.0, help='time limit of one trace (s)')
ap.add_argument('--chain-check', action='store_true')
ap.add_argument('--dps', type=int, default=20)
args = ap.parse_args()
out = args.out or ('data/mi_depart.json' if not args.quick else 'misc/reviews/nick-miyatake/quick/G2/mi_depart.json')
logf = open(out.replace('.json', '.log'), 'a')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.write(s + '\n'); logf.flush()
res = {}
mp.mp.dps = args.dps
PI = mp.pi
PR = [p for p in range(2, 200) if all(p % d for d in range(2, p))]
def lpf(m):
    k, d, lp = m, 2, 1
    while k > 1:
        if k % d == 0: lp = d; k //= d
        else: d += 1
    return lp
def nextprime(P): return next(p for p in PR if p > P)
def notS(P, cut):
    q = nextprime(P)
    return [m for m in range(q, 40) if m * m - q * q <= cut and lpf(m) > P]

def Xi(t):
    t = mp.mpmathify(t)
    s = mp.mpf(1) / 2 + (-1j if mp.im(t) > 0 else 1j) * t
    return s * (s - 1) / 2 * PI ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
def DP(t, P, cut=30):
    t = mp.mpmathify(t)
    sp = mp.mpf(1) / 4 + 0.5j * t; spp = mp.mpf(1) / 4 - 0.5j * t
    tot = 0
    for m in notS(P, cut):
        a = PI * m * m
        tot += 0.5 * (a ** (-sp) * mp.gammainc(sp, a) + a ** (-spp) * mp.gammainc(spp, a))
    return (t * t + mp.mpf(1) / 4) * tot
def XiP(t, P, cut=30):
    return Xi(t) + DP(t, P, cut)
def M(t):
    t = mp.mpmathify(t)
    return 0.5 * (t * t + 0.25) * PI ** (-0.25) * abs(mp.gamma(mp.mpf(1) / 4 + 0.5j * t))
def cP(P):
    return 4 * sum(mp.exp(-PI * m * m) * (PI * m * m - mp.mpf(1) / 4) for m in notS(P, 200))
def XiP_direct(t, P):
    t = mp.mpmathify(t); s = mp.mpf(1) / 2 + 1j * t
    tot = 0
    for n in range(1, 41):
        if lpf(n) <= P:
            a = PI * n * n
            tot += a ** (-s / 2) * mp.gammainc(s / 2, a) + a ** (-(1 - s) / 2) * mp.gammainc((1 - s) / 2, a)
    return mp.mpf(1) / 2 + s * (s - 1) / 2 * tot
def law(x, P, use_c=False):
    x = mp.mpmathify(x)
    D = cP(P) if use_c else mp.re(DP(x, P))
    return 2 * mp.log(D / M(x)) / mp.log(x / (2 * PI))

# ---- evaluator check
mp.mp.dps = 60
pts = [] if args.skip_check else [(45.418 + 4.421j, 2), (30.0 + 2j, 2), (70.0 + 19.0j, 2), (111.38 + 3.66j, 3), (148.3 + 21.5j, 3), (50.0, 2)]
dm = 0
for t, P in pts:
    a, b = XiP(t, P), XiP_direct(t, P)
    c = XiP(t, P, 120)
    log('evaluator check P=%d t=%s: Xi+D_P = %s, Prop incgamma = %s, relative diff %.1e, cut 60 vs 120 relative diff %.1e'
        % (P, t, mp.nstr(a, 8), mp.nstr(b, 8), float(abs(a - b) / abs(a)), float(abs(a - c) / abs(a))))
mp.mp.dps = args.dps

def newton(P, t0, cut=30):
    return mp.findroot(lambda z: XiP(z, P, cut), mp.mpmathify(t0), solver='secant', tol=mp.mpf(10) ** (-(args.dps - 8)), maxsteps=60, verify=False)

if args.part in ('table', 'all'):
    rows = []
    for P, seeds in ((2, [45, 58, 70]), (3, [111, 128, 148])):
        for xs in seeds:
            y0 = float(law(xs, P))
            z = newton(P, mp.mpc(xs, y0))
            z2 = newton(P, mp.mpc(xs, y0), 120)
            x, y = float(mp.re(z)), abs(float(mp.im(z)))
            lc, ld = float(law(x, P, True)), float(law(x, P))
            rows.append(dict(P=P, re=x, im=y, law_c=lc, law=ld, resid=float(abs(XiP(z, P))), cutdiff=float(abs(z - z2))))
            log('Table dep P=%d: z = %.4f + %.4fi (|Xi_P| = %.1e, cut difference %.1e), law with c_P %.2f, law with D_P %.2f (seed %d + %.2fi)' % (P, x, y, rows[-1]['resid'], rows[-1]['cutdiff'], lc, ld, xs, y0))
    res['table_dep'] = rows
    res['c_P'] = {P: float(cP(P)) for P in (2, 3)}
    log('c_2 = %.6e, c_3 = %.6e' % (res['c_P'][2], res['c_P'][3]))
    # D_P / c_P near T_P^*
    for P, T in ((2, 39.7), (3, 104.7)):
        log('D_P(T_P^*)/c_P at P=%d, T=%.1f: %.3f' % (P, T, float(mp.re(DP(T, P)) / cP(P))))

def arg_change(f, path, tol=0.5, maxdepth=14):
    # path: list of complex points (mp); returns the total change of arg f along the polyline, by bisection
    tot = 0; nev = 0
    for z0, z1 in zip(path[:-1], path[1:]):
        stack = [(z0, f(z0), z1, f(z1), 0)]; nev += 2
        while stack:
            a, fa, b, fb, d = stack.pop()
            dth = mp.arg(fb / fa)
            if abs(dth) < tol or d >= maxdepth:
                if d >= maxdepth: log('WARNING: maximal depth at', a)
                tot += dth
            else:
                m = (a + b) / 2; fm = f(m); nev += 1
                stack.append((m, fm, b, fb, d + 1)); stack.append((a, fa, m, fm, d + 1))
    return tot, nev

ref = [mp.mpf(l) for l in open(__import__('os').path.exists('data/zeros_1700.txt') and 'data/zeros_1700.txt' or 'zeros_1700.txt').read().split()]
def Ncount(X): return sum(1 for g in ref if g < X)

def count(P, X, Yfac, dps):
    mp.mp.dps = dps
    yP = float(law(X, P))
    Y = Yfac * max(yP, 2.0) + 6
    f = lambda z: XiP(z, P)
    path = [mp.mpc(X, 0), mp.mpc(X, Y), mp.mpc(0, Y), mp.mpc(0, 0)]
    # subdivide the long edges to keep the bisection stack small
    pts = [mp.mpc(X, 0)]
    for (za, zb) in ((mp.mpc(X, 0), mp.mpc(X, Y)), (mp.mpc(X, Y), mp.mpc(0, Y)), (mp.mpc(0, Y), mp.mpc(0, 0))):
        n = int(abs(zb - za) / 0.4) + 1
        pts += [za + (zb - za) * k / n for k in range(1, n + 1)]
    t0 = time.time()
    ch, nev = arg_change(f, pts)
    top = max(float(abs(DP(mp.mpc(X * k / 6, Y), P) / Xi(mp.mpc(X * k / 6, Y)))) for k in range(0, 7))
    val = 2 * ch / (2 * PI)
    mp.mp.dps = args.dps
    return float(val), Y, top, nev, time.time() - t0, yP

if args.part in ('counts', 'all'):
    cres = []
    for case in args.cases.split(','):
        P, X = case.split(':'); P = int(P); X = mp.mpf(X)
        v1, Y1, top1, nev1, tm1, yP = count(P, X, 2.5, args.dps)
        v2, Y2, top2, nev2, tm2, _ = (v1, 0, 0, 0, 0, 0) if args.single else count(P, X, 3.5, args.dps + 8)
        NX = Ncount(X)
        n1 = round(v1)
        log('count P=%d X=%s: N_P = %.4f (Y=%.1f, dps %d, %d evals, %.1fs, max|D_P/Xi| top %.1e) and %.4f (Y=%.1f, dps %d); integer %d, fractional parts %.1e %.1e; N(X) = %d; N_P - N = %d; y_P/4 = %.2f'
            % (P, X, v1, Y1, args.dps, nev1, tm1, top1, v2, Y2, args.dps + 8, n1, abs(v1 - n1), abs(v2 - round(v2)), NX, n1 - NX, yP / 4))
        cres.append(dict(P=P, X=float(X), NP=v1, NP2=v2, N=NX, diff=n1 - NX, yP4=yP / 4, Y=Y1, top=top1))
    res['counts'] = cres

def dens(x):
    L = math.log(x / (2 * math.pi))
    return L / (4 * math.pi) + math.pi / (16 * L)
def trace(P, z0, sign, xstop, tmax):
    # continuation of the chain from z0 in the direction sign (+1 up, -1 down) until Re passes xstop; returns the list without z0
    chain = [z0]; t0 = time.time(); outl = []
    while (mp.re(chain[-1]) < xstop if sign > 0 else mp.re(chain[-1]) > xstop):
        zc = chain[-1]; x = float(mp.re(zc))
        sp = 1.0 / dens(x)
        if len(chain) >= 2:
            d = chain[-1] - chain[-2]
        else:
            slope = float(law(x + 1, P) - law(x, P))
            d = sign * sp * mp.mpc(1, slope)
        found = None
        for sc in (1.0, 0.7, 1.3, 0.5, 1.6):
            try:
                zn = newton(P, zc + d * sc)
            except Exception:
                continue
            dx = sign * float(mp.re(zn - zc))
            if 0.4 * sp < dx < 1.8 * sp and abs(XiP(zn, P)) < 1e-6 * abs(XiP(zn + 0.05, P)):
                found = zn; break
        if found is None:
            log('chain STOP at', mp.nstr(zc, 8)); break
        chain.append(found); outl.append(found)
        if time.time() - t0 > tmax:
            log('chain time limit at', mp.nstr(found, 8)); break
    return outl

if args.part in ('chain', 'all'):
    P = 2
    t0 = time.time(); Xend = args.chain_to
    if args.seg:
        lo, hi = [float(v) for v in args.seg.split(':')]
        if args.seed:
            sx, sy = [float(v) for v in args.seed.split(':')]; zs = newton(P, mp.mpc(sx, sy))
        else:
            y0 = float(law(lo, P)); zs = newton(P, mp.mpc(lo, y0))
        if not abs(XiP(zs, P)) < 1e-6 * abs(XiP(zs + 0.05, P)):
            log('seed NOT converged', mp.nstr(zs, 8)); raise SystemExit(1)
        up = trace(P, zs, 1, hi, 24); dn = trace(P, zs, -1, lo - 12, 6)
        chain = sorted([z for z in [zs] + up + dn if float(mp.im(z)) > 1e-3], key=lambda z: float(mp.re(z)))
        log('segment %s: seed %s, %d zeros from %s to %s, %.1fs' % (args.seg, mp.nstr(zs, 7), len(chain), mp.nstr(chain[0], 7), mp.nstr(chain[-1], 7), time.time() - t0))
    else:
        z = newton(P, mp.mpc(41.19, 1.46))
        log('chain start', mp.nstr(z, 8))
        chain = [z] + trace(P, z, 1, Xend, args.tmax)
    xs = [float(mp.re(c)) for c in chain]
    log('chain of %d zeros from %s to %s; zeros with real part below %.0f: %d; time %.1fs'
        % (len(chain), mp.nstr(chain[0], 7), mp.nstr(chain[-1], 7), Xend, sum(1 for v in xs if v < Xend), time.time() - t0))
    import scipy.integrate as si
    pred = si.quad(dens, 41.19, Xend)[0]
    log('integral of the density from 41.19 to %.0f: %.2f' % (Xend, pred))
    res['chain'] = [(mp.nstr(mp.re(c), 20), mp.nstr(mp.im(c), 20)) for c in chain]
    res['chain_pred'] = pred
    gaps = [b - a for a, b in zip(xs[:-1], xs[1:])]
    log('chain spacings: min %.3f max %.3f' % (min(gaps) if gaps else 0, max(gaps) if gaps else 0))
json.dump(res, open(out, 'w'), indent=1, default=str)
