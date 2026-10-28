# ht_lcert.py: ball arithmetic (python-flint, Arb) check of the log concavity certificate of Nick Polson's residue note (R:360 to 420, 540 to 560):
# Proposition explicit (constants M_j), the pair identity ratios, the Taylor model windows of tab:cert, and g = f f'' - f'^2 at y = 0.005.
# Claims supported: R:363 (c, gamma_120, gamma_121), R:370 (log M_0 = 206.26, M_(j+1)/M_j = 1.0005 c), R:385 and R:553 (pair ratios 2.83, 1.39,
# 0.622, 0.034, 2.2e-5; the enclosure [0.62185566984, 0.62185566985]), R:404 to 407 and R:553 (Taylor steps 25856, 73, 74; 18487, 4827, 2091, 360,
# 94), R:418 (g ~ -6.6e-40 at y = 0.005).
# Conventions. F(s) = xi(1/2 + sqrt s)/xi(1/2), xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s). lambda_n = gamma_n^2, n <= 120, w_n = 1/F'(-lambda_n) =
# 2 i gamma_n xi(1/2)/xi'(1/2 + i gamma_n), f_120(y) = sum w_n exp(-lambda_n y), c = (lambda_120 + lambda_121)/2, and
# M_j = (1/2 pi) int |-c + it|^j /|F(-c + it)| dt = (1/pi) int_0^inf (c^2 + t^2)^(j/2)/|F(-c + it)| dt by the symmetry t -> -t.
# Zeros: acb.zeta_zeros (certified balls). xi'(rho) = (rho (rho-1)/2) pi^(-rho/2) Gamma(rho/2) zeta'(rho), zeta'(rho) by the trapezoid rule
# on the circle of radius 1/2 with N = 64 nodes (the aliasing error is below (1/540)^64 times the maximum of |zeta| on the circle); the
# difference to N = 96 is logged and the ball radius of every w_n is widened by 100 times that difference.
# Enclosure of M_j. |1/F(-c + it)|^2 = 1/(F(-c+it) F(-c-it)) is an analytic function of t, so the integral over [0, T] is enclosed by Arb's rigorous
# Gauss Legendre integration (acb.integral, with the root taken as nan where a ball meets the branch cut, so that the integrator subdivides),
# on the chunks [0, 64], [64, 256], ..., [64 * 4^8, 64 * 4^9]. The tail beyond T = 64 * 4^9 is NOT enclosed: the log reports the integrand at T
# relative to its peak and the share of the last chunk; an analytic tail bound (Stirling) is still needed, and the status of M_j is: enclosed on
# [0, T], tail unproved. The line Re s = -c is taken at the exact midpoint (a number between lambda_120 and lambda_121), which any line between the poles serves.
# Taylor model. As in R:395 to 410 with K = 6: at each step y_i the derivatives g^(m)(y_i), m < 6, are formed by Leibniz from balls for
# f^(j) = sum + [-e_j, e_j], e_j = M_j^hi exp(-c y_i); |g^(6)| is bounded with A_j(y_i) = sum |w_n| lambda_n^j exp(-lambda_n y_i) + e_j; the
# step h is the largest of a doubling ladder (2 x previous step, halved until the test passes). The step counts are therefore not the
# chapter's (its rule for h is not given); what is tested is that the windows are covered.
# Options: --quick (fewer chunks, 200 bits, one window), --out FILE, --mj J (stage 1: only M_J, written to --out; run once for each J = 0..8, 2 minutes each), --mfiles PATTERN (stage 2: read M_J from PATTERN % J, then pair ratios, Taylor windows, g(0.005); 80 s). Output json and log.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, time, argparse
from flint import arb, acb, ctx, fmpq
ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/ht_lcert.json')
ap.add_argument('--part', default='all'); ap.add_argument('--mj', type=int, default=-1); ap.add_argument('--mfiles', default='')
a = ap.parse_args()
LOG = open(a.out.replace('.json', '.log'), 'w')
def lg(*x):
    s = ' '.join(str(t) for t in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
ctx.prec = 200 if a.quick else 400
R = {}
t0 = time.time()
pi = arb.pi(); half = arb(1) / 2
def xi(s):
    return s * (s - 1) / 2 * (-(s / 2) * pi.log()).exp() * (s / 2).gamma() * s.zeta()
xi_half = xi(acb(half))
xi_half = xi_half.real
lg('xi(1/2) =', xi_half.str(20))
NZ = 121
Z = acb.zeta_zeros(1, NZ)
gam = [z.imag for z in Z]
lam = [g * g for g in gam]
lg('gamma_120 =', gam[119].str(15), ' gamma_121 =', gam[120].str(15))
c = (lam[119] + lam[120]) / 2
c = arb(c.mid())   # an exact number between lambda_120 and lambda_121 (any line between the poles serves)
lg('c =', c.str(15))
R['c'] = c.str(15); R['gamma120'] = gam[119].str(15); R['gamma121'] = gam[120].str(15)

def zeta_deriv(rho, N):
    r = arb(1) / 2; tot = acb(0)
    for k in range(N):
        th = 2 * pi * k / N
        e = acb(th.cos(), th.sin())
        tot += (rho + r * e).zeta() * e.conjugate() if False else (rho + r * e).zeta() / e
    return tot / (N * r)
nz = 30 if a.quick else 120
w = []; dmax = arb(0)
for n in range(nz):
    rho = acb(half, gam[n])
    d1 = zeta_deriv(rho, 64 if not a.quick else 32); d2 = zeta_deriv(rho, 96 if not a.quick else 48)
    diff = abs(d1 - d2) / abs(d1)
    dmax = max(dmax, diff, key=lambda x: float(x.mid())) if False else (diff if float(diff.mid()) > float(dmax.mid()) else dmax)
    r = rho
    xp = r * (r - 1) / 2 * (-(r / 2) * pi.log()).exp() * (r / 2).gamma() * d1
    wn = 2 * acb(0, 1) * gam[n] * xi_half / xp
    wr = wn.real
    # widen by 100 x the discretisation difference
    err = abs(wr) * diff.upper() * 100
    wr = wr + (err.union(-err) if True else 0)
    w.append(wr)
lg('w_1 =', w[0].str(12), ' w_2 =', w[1].str(12), ' w_5 =', w[4].str(12))
lg('largest relative difference of zeta\' between N=64 and N=96 circles:', dmax.str(5))
R['w1'] = w[0].str(12); R['w2'] = w[1].str(12); R['w5'] = w[4].str(12)
if nz >= 120:
    lg('w_120 =', w[119].str(10), ' log|w_120| =', abs(w[119]).log().str(8))
    R['w120'] = w[119].str(10); R['logabsw120'] = abs(w[119]).log().str(8)
lg('weights done, time %.0f s' % (time.time() - t0))

# ------------------------------------------------------------ M_j (rigorous analytic integration of the real analytic integrand)
JMAX = 8
sqc = c.sqrt()
I_ = acb(0, 1)
def sq(z, an):
    # principal square root; as an analytic function only where the ball stays in the right half plane (no branch cut), else nan so that the integrator subdivides
    if an and not (z.real > 0):
        return acb(float("nan"), float("nan"))
    return z.sqrt()
def Fabs2(z, an):
    # |F(-c + i t)|^2 = F(-c + i t) F(-c - i t) as an analytic function of t = z (F even in sqrt s, so either branch of the root serves)
    Wp = I_ * sqc * sq(1 - I_ * z / c, an)
    Wm = I_ * sqc * sq(1 + I_ * z / c, an)
    return xi(acb(half) + Wp) * xi(acb(half) + Wm) / (xi_half * xi_half)
def make_f(j):
    def f(z, an):
        return sq(c * c + z * z, an) ** j / sq(Fabs2(z, an), an)
    return f
edges = [0, 64, 256] if a.quick else [0] + [64 * 4 ** k for k in range(8)]
RT = arb(10) ** (-8 if a.quick else -12)
peak0 = make_f(0)(acb(0), False).real
lg('integrand at t = 0: %s' % peak0.str(8, radius=False))
M = [None] * (JMAX + 1)
if a.mfiles:
    for j in range(JMAX + 1):
        M[j] = arb(json.load(open(a.mfiles % j))['M'])
    lg('M_j read from', a.mfiles, [m.str(8, radius=False) for m in M])
else:
    for j in ([a.mj] if a.mj >= 0 else range(JMAX + 1)):
        f = make_f(j)
        peak = peak0 * c ** j
        tot = arb(0); parts = []
        for e0, e1 in zip(edges[:-1], edges[1:]):
            val = acb.integral(f, e0, e1, rel_tol=RT, abs_tol=peak * arb(10) ** -40, eval_limit=10 ** 7, depth_limit=200)
            parts.append(val.real); tot += val.real
            lg('  j=%d chunk [%d, %d]: %s ; time %.0f s' % (j, e0, e1, (val.real / pi).str(8, radius=False), time.time() - t0))
        Mj = tot / pi
        M[j] = Mj
        ref = f(acb(edges[-1]), False).real
        lg('M_%d = %s (log M_%d = %s); last chunk share %s; integrand at t=%d relative to its peak %s; time %.0f s' % (j, Mj.str(12), j, Mj.log().str(10, radius=False), (parts[-1] / tot).str(3, radius=False), edges[-1], (ref / peak).str(3, radius=False), time.time() - t0))
        if a.mj >= 0:
            json.dump({'j': j, 'M': Mj.str(30, radius=True), 'T': edges[-1]}, open(a.out, 'w'))
            lg('written', a.out); sys.exit(0)
lg('integration range [0, %d], tail beyond NOT enclosed (chunk shares above show how much sits near the end)' % edges[-1])
for j in range(3):
    lg('M_%d/M_%d / c = %s' % (j + 1, j, (M[j + 1] / M[j] / c).str(8, radius=False)))
R['M'] = [m.str(14) for m in M]
R['T'] = edges[-1]
Mhi = [m.upper() for m in M]
# ------------------------------------------------------------ pair identity
def pair_ratio(y, with_trunc):
    wl = [(w[i], lam[i]) for i in range(nz)]
    E = [(-lam[i] * y).exp() for i in range(nz)]
    lead = abs(w[0] * w[1]) * (lam[1] - lam[0]) ** 2 * (-(lam[0] + lam[1]) * y).exp()
    tot = arb(0)
    for i in range(nz):
        for j in range(i + 1, nz):
            if i == 0 and j == 1: continue
            tot += abs(w[i] * w[j]) * (lam[i] - lam[j]) ** 2 * E[i] * E[j]
    ratio = tot / lead
    if with_trunc:
        A = [sum((abs(w[i]) * lam[i] ** j * E[i] for i in range(nz)), arb(0)) for j in range(3)]
        e = [Mhi[j] * (-c * y).exp() for j in range(3)]
        d = e[0] * A[2] + A[0] * e[2] + 2 * A[1] * e[1] + e[0] * e[2] + e[1] ** 2
        ratio = (tot + d) / lead
    return ratio
if a.part in ('all', 'pair'):
    for ys in ['0.02', '0.022', '0.025', '0.04', '0.08']:
        y = arb(fmpq(int(round(float(ys) * 1000)), 1000))
        r0 = pair_ratio(y, False); r1 = pair_ratio(y, True)
        lg('pair ratio at y = %s: %s ; with the truncation bound: %s' % (ys, r0.str(12, radius=True), r1.str(12, radius=True)))
        R['pair_' + ys] = [r0.str(12), r1.str(12)]
    # monotone on a grid
    prev = None; mono = True
    for k_ in range(25, 301, 5):
        y = arb(fmpq(k_, 1000)); r = pair_ratio(y, False)
        if prev is not None and not (r < prev): mono = False
        prev = r
    lg('pair ratio strictly decreasing and below 1 on the grid y = 0.025:0.005:0.3 (balls):', mono, ' last value', prev.str(5))
    R['pair_monotone_grid'] = bool(mono)

# ------------------------------------------------------------ Taylor model
from math import comb
def step_data(y):
    E = [(-lam[i] * y).exp() for i in range(nz)]
    fj = []; Aj = []
    for j in range(JMAX + 1):
        s = arb(0); sa = arb(0)
        for i in range(nz):
            t = w[i] * lam[i] ** j * E[i]
            s += t if j % 2 == 0 else -t
            sa += abs(w[i]) * lam[i] ** j * E[i]
        e = Mhi[j] * (-c * y).exp()
        fj.append(s + e.union(-e)); Aj.append(sa + e)
    return fj, Aj
def taylor_window(ya, yb, h0=arb(1) / 10 ** 6):
    y = ya; steps = 0; h = h0; hmin = None
    while y < yb:
        fj, Aj = step_data(y)
        G = []
        for m in range(6):
            s = arb(0)
            for kk in range(m + 1):
                s += comb(m, kk) * (fj[kk] * fj[m - kk + 2] - fj[kk + 1] * fj[m - kk + 1])
            G.append(s)
        Bd = sum((comb(6, kk) * (Aj[kk] * Aj[8 - kk] + Aj[kk + 1] * Aj[7 - kk]) for kk in range(7)), arb(0))
        if not (G[0] < 0):
            return steps, False, y
        hh = min(2 * h, yb - y)
        fact = [1, 1, 2, 6, 24, 120, 720]
        while True:
            P = G[0] + sum((abs(G[m]) * hh ** m / fact[m] for m in range(1, 6)), arb(0)) + Bd * hh ** 6 / 720
            if P < 0: break
            hh = hh / 2
            if hh < arb(10) ** -30: return steps, False, y
        y = y + hh; h = hh; steps += 1
        hmin = hh if hmin is None or float(hh.mid()) < float(hmin.mid()) else hmin
    return steps, True, hmin
wins = [('0.020', '0.080'), ('0.012', '0.020'), ('0.011', '0.012'), ('0.010', '0.011'), ('0.009', '0.010'), ('0.0085', '0.009'), ('0.008', '0.0085')]
if a.quick: wins = [('0.020', '0.080')]
if a.part in ('all', 'taylor'):
    for wa, wb in wins:
        ya = arb(fmpq(int(round(float(wa) * 10000)), 10000)); yb = arb(fmpq(int(round(float(wb) * 10000)), 10000))
        t1 = time.time()
        steps, ok, info = taylor_window(ya, yb)
        lg('Taylor window [%s, %s]: covered %s in %d steps (smallest step %s), %.0f s' % (wa, wb, ok, steps, info.str(3) if ok else ('failed at y=' + info.str(8)), time.time() - t1))
        R['taylor_%s_%s' % (wa, wb)] = [bool(ok), steps]
    # g at y = 0.005
    y = arb(fmpq(5, 1000))
    fj, Aj = step_data(y)
    g = fj[0] * fj[2] - fj[1] ** 2
    lg('g(0.005) from the 120 terms with the truncation balls = %s ; largest term of f at 0.005 over f: see A_0/|f_0| = %s' % (g.str(8), (Aj[0] / abs(fj[0])).str(5, radius=False)))
    R['g_0005'] = g.str(8)
lg('time %.0f s' % (time.time() - t0))
json.dump(R, open(a.out, 'w'), indent=1)
