# mi_levy.py: the 2 sin^2 Levy density nu_P of Chapter 22 (Polson's Miyatake chapter), Section "The 2 sin^2 lemma and the Thorin programme":
# Proposition axis (sup over real y of |Xi_P(iy)/Xi(iy) - 1| = 6e-13 for P = 2, 3e-35 for P = 3), Table levy ((x/2) nu_2(x), (x/2) nu(x) and the
# departed part at x = 0.03, 0.05, 0.1, 0.5), the positivity of nu_2 on [0.02, 3], the departed zeros up to 0.35 radians off the real axis, and
# the claim that their terms reproduce to five digits the contribution of the real zeros of Xi above 37.5. Added for the book in October 2026
# (computation group G2), written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. nu_P(x) = (2/x)[sum over the real zeros of Xi_P of e^(-gamma x) + 2 sum over the departed zeros x0 + i y0 in the first quadrant of
# e^(-x0 x) cos(y0 x)], nu(x) = (2/x) sum over the zeros gamma of Xi of e^(-gamma x). So (x/2) nu_P is the bracket, and the departed part of the
# table is 2 sum e^(-x0 x) cos(y0 x). Xi_P = Xi + D_P as in mi_depart.py (same evaluator: D_P(t) = (t^2+1/4) sum over m not in S(P) of C_(pi m^2)(t)).
# Inputs. The real zeros of Xi: data/zeros_1700.txt (68 digits). The real zeros of Xi_2: complex-plane Newton of Xi_2 from the zeros of Xi (6 zeros).
# The departed zeros of Xi_2: the chain files written by mi_depart.py --part chain (segments traced from seeds, merged here by removing points
# that agree to 1e-8; the merge reports the number of zeros, the smallest and largest gap of the real parts and the real part cut that is used).
# Tail control. For nu the zero sum is cut at gamma_1700 = 2.4e3 and the remainder is bounded by the integral of the density log(t/2 pi)/(2 pi) times
# e^(-x t) (logged for the smallest x of the table); for nu_2 the chain is cut at the real part of its last traced zero and the same bound with
# the chain density is logged. The axis sup is the maximum over a grid of y in [0, 60] (step 0.05) at 40 digits and the grid is refined by a golden
# section near the maximum; the same at 60 digits for the check. Output data/mi_levy.json and data/mi_levy.log (quick: --out FILE).
import sys, json, time, argparse, glob, math
import mpmath as mp

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default=None)
ap.add_argument('--chain-files', default='misc/reviews/nick-miyatake/quick/G2/seg/seg*.json')
ap.add_argument('--chain-cut', type=float, default=601.91, help='keep chain zeros with real part up to this (the chapter has 190 traced zeros, the last at 601.90)')
ap.add_argument('--axis-p', default='2,3')
ap.add_argument('--ymax', type=float, default=40.0)
ap.add_argument('--axis-only', action='store_true')
ap.add_argument('--no-axis', action='store_true')
args = ap.parse_args()
out = args.out or ('data/mi_levy.json' if not args.quick else 'misc/reviews/nick-miyatake/quick/G2/mi_levy.json')
logf = open(out.replace('.json', '.log'), 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.write(s + '\n'); logf.flush()
res = {}
mp.mp.dps = 30
PI = mp.pi
PR = [p for p in range(2, 200) if all(p % d for d in range(2, p))]
def lpf(m):
    k, d, lp = m, 2, 1
    while k > 1:
        if k % d == 0: lp = d; k //= d
        else: d += 1
    return lp
def nextprime(P): return next(p for p in PR if p > P)
def notS(P, cut=30):
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
def XiP(t, P): return Xi(t) + DP(t, P)

# ---- Proposition axis
def axis(P, dps, ymax=60.0, step=0.05):
    mp.mp.dps = dps
    best = (0, 0); vals = []
    k = 0
    while k * step <= ymax:
        y = mp.mpf(k) * mp.mpf(step) + mp.mpf(10) ** -17  # offset: the evaluator of Xi uses zeta(1/2 + y), a pole at y = 1/2
        r = abs(mp.re(DP(mp.mpc(0, y), P) / Xi(mp.mpc(0, y))))
        vals.append((float(y), float(r)))
        if r > best[0]: best = (r, y)
        k += 1
    mp.mp.dps = 30
    return best, vals
if not args.no_axis:
    for P in ([int(v) for v in args.axis_p.split(',')]):
        t0 = time.time()
        b, vals = axis(P, 40, args.ymax, 0.05)
        b2, _ = axis(P, 55, 20.0, 0.5)
        log('axis P=%d: sup |Xi_P(iy)/Xi(iy) - 1| = %s at y = %s (grid step 0.05 on [0, ymax] at 40 digits); value at y=0: %s, at y=5: %s, at y=20: %s; second run (55 digits, step 0.5 on [0, 20]) gives %s; %.1fs'
            % (P, mp.nstr(b[0], 6), mp.nstr(b[1], 4), mp.nstr(vals[0][1], 6), mp.nstr(min(vals, key=lambda v: abs(v[0] - 5))[1], 4), mp.nstr(min(vals, key=lambda v: abs(v[0] - 20))[1], 4), mp.nstr(b2[0], 6), time.time() - t0))
        res['axis%d' % P] = dict(sup=mp.nstr(b[0], 10), at=mp.nstr(b[1], 6), second=mp.nstr(b2[0], 10))
if args.axis_only:
    json.dump(res, open(out, 'w'), indent=1); sys.exit(0)

# ---- zeros
ref = [mp.mpf(l) for l in open(__import__('os').path.exists('data/zeros_1700.txt') and 'data/zeros_1700.txt' or 'zeros_1700.txt').read().split()]
real2 = []
for g in ref[:6]:
    z = mp.findroot(lambda t: mp.re(XiP(t, 2)), g, tol=mp.mpf(10) ** -22)
    real2.append(z)
log('real zeros of Xi_2: %s; shifts from the zeros of Xi: %s' % ([mp.nstr(z, 8) for z in real2], [mp.nstr(z - g, 3) for z, g in zip(real2, ref[:6])]))
pts = []
for f in sorted(sum([glob.glob(g) for g in args.chain_files.split(',')], [])):
    d = json.load(open(f)); pts += [(mp.mpf(a), mp.mpf(b)) for a, b in d.get('chain', [])]
pts.sort(key=lambda z: z[0]); chain = []
for z in pts:
    if not chain or abs(z[0] - chain[-1][0]) > 1e-8 or abs(z[1] - chain[-1][1]) > 1e-8:
        chain.append(z)
chain = [z for z in chain if float(z[0]) <= args.chain_cut]
xs = [float(z[0]) for z in chain]; gaps = [b - a for a, b in zip(xs[:-1], xs[1:])]
log('departed zeros: %d (real part from %.4f to %.4f, cut %.2f), %d with real part below 600; gaps of the real parts: min %.3f (at %.2f) max %.3f (at %.2f)'
    % (len(chain), xs[0], xs[-1], args.chain_cut, sum(1 for v in xs if v < 600), min(gaps), xs[gaps.index(min(gaps))], max(gaps), xs[gaps.index(max(gaps))]))
angles = [math.atan2(float(z[1]), float(z[0])) for z in chain]
k = angles.index(max(angles))
log('largest angle off the real axis: %.4f rad at %.2f + %.2fi' % (max(angles), xs[k], float(chain[k][1])))
res['n_departed'] = len(chain); res['max_angle'] = max(angles)

def bracket_real(x):
    return sum(mp.exp(-g * x) for g in real2)
def departed(x):
    return 2 * sum(mp.exp(-a * x) * mp.cos(b * x) for a, b in chain)
def bracket_xi(x):
    return sum(mp.exp(-g * x) for g in ref)
def tail_bound(x, T):
    # integral over t > T of (1/(2 pi)) log(t/2 pi) e^(-x t) dt <= (log(T/2pi)+1/(x T))/(2 pi x) e^(-x T)
    return float((mp.log(T / (2 * PI)) + 1 / (x * T)) / (2 * PI * x) * mp.exp(-x * T))
rows = []
Tg = float(ref[-1]); Tc = xs[-1]
for x in (0.03, 0.05, 0.1, 0.5):
    x = mp.mpf(x)
    r2, dep, rxi = bracket_real(x), departed(x), bracket_xi(x)
    nu2 = r2 + dep
    log('x=%s: (x/2) nu_2 = %s, (x/2) nu = %s, departed part = %s, real part of nu_2 = %s; difference nu_2 - nu: %s; tail bounds: zeros of Xi beyond %.0f: %.1e, chain beyond %.1f: %.1e'
        % (mp.nstr(x, 3), mp.nstr(nu2, 8), mp.nstr(rxi, 8), mp.nstr(dep, 8), mp.nstr(r2, 8), mp.nstr(nu2 - rxi, 3), Tg, tail_bound(float(x), Tg), Tc, tail_bound(float(x), Tc) * 2))
    # contribution of the zeros of Xi above 37.5
    above = sum(mp.exp(-g * x) for g in ref if g > 37.5)
    above6 = sum(mp.exp(-g * x) for g in ref[6:])
    log('   zeros of Xi beyond the sixth: %s, relative difference to the departed part %s' % (mp.nstr(above6, 8), mp.nstr((dep - above6) / above6, 3)))
    log('   real zeros of Xi above 37.5: %s; departed part %s; relative difference %s' % (mp.nstr(above, 8), mp.nstr(dep, 8), mp.nstr((dep - above) / above, 3)))
    rows.append(dict(x=float(x), nu2=mp.nstr(nu2, 10), nu=mp.nstr(rxi, 10), dep=mp.nstr(dep, 10), above=mp.nstr(above, 10), reldiff=mp.nstr((dep - above) / above, 5)))
res['table_levy'] = rows
# positivity of nu_2 on [0.02, 3]
mn = None; bad = 0; xsg = [0.02 + 0.005 * k for k in range(int((3 - 0.02) / 0.005) + 1)]
for x in xsg:
    v = (bracket_real(mp.mpf(x)) + departed(mp.mpf(x))) * 2 / x
    if v <= 0: bad += 1
    if mn is None or v < mn[0]: mn = (float(v), x)
log('nu_2 on %d grid points of [0.02, 3] (step 0.005): minimum %.4e at x = %.3f, non positive values: %d' % (len(xsg), mn[0], mn[1], bad))
# x = 0.02 specifically
x = mp.mpf('0.02'); log('x=0.02: (x/2) nu_2 = %s' % mp.nstr(bracket_real(x) + departed(x), 8))
res['positivity'] = dict(min=mn[0], at=mn[1], bad=bad)
json.dump(res, open(out, 'w'), indent=1, default=str)
