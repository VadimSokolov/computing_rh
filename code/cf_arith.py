# cf_arith.py: the numerical claims of Chapter 34 (Cauchy deconvolution from U_1 to U_{1/2}, prefix cf) about the zero free constants of
# admissible pairs (p, w): Table cf:tab:arith (every row from Corollary cf:cor:ystar, y* = 0 when the equation has no positive root), the
# cap of Proposition cf:prop:barrier2 per row (the sharper cap (2/A)(a1 - a0 - a0 log(a1/a0)) and the source's 2(a1 - a0)/A), the
# Stechkin row of Lemma cf:lem:stechkin and Proposition cf:prop:stechkin (0.05195 for 3+4cos+cos2; 0.06809, 0.09800, 0.10367 for the
# best polynomials of degree 2, 3, 4; the 1/9.646 of degree 4). Added for the book in October 2026 (group CF-LP), written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. p(phi) = sum_{k=0}^K a_k cos(k phi) >= 0 with a_k >= 0, A = a_1 + ... + a_K; the level kernel w has Laplace transform
# phi_w(lambda) = int exp(-x lambda) w(dx); I = int_0^inf phi_w, J(y) = int_0^inf exp(-y lambda) phi_w. y*(p, w) is the root of
# a_1 J(y) = a_0 I + (A/2) phi_w(0) when (a_1 - a_0) I > (A/2) phi_w(0), and 0 otherwise (Corollary cf:cor:ystar). Rows:
#   single level   phi = exp(-x l):              J = 1/(x+y), I = 1/x. Closed form of the optimum over x: y* = 2 (sqrt a1 - sqrt a0)^2 / A at
#                  x = 2 (sqrt(a0 a1) - a0)/A (derived in the log by calculus and checked against a numerical maximisation over x).
#   two levels     phi = exp(-x l)(1 + x l):     J = 1/(x+y) + x/(x+y)^2, I = 2/x (the limit kappa = x1/x2 -> 1 of delta_x1 - kappa delta_x2).
#   derivatives    phi = exp(-x l) q(l), q(0) = 1, q >= 0 on [0, inf), sum_j q_j j! Re(x+iu)^(-j-1) >= 0 for all u (Lemma cf:lem:admissible
#                  (ii)); J = sum_j q_j j!/(x+y)^(j+1), I = sum_j q_j j!/x^(j+1). The constraints are imposed on grids (s = x l in [0, 60],
#                  t = u/x on a log grid up to 1e4), so this row is a computed estimate, as the chapter's wording says.
# Method. Polynomials: Fejer Riesz, p = |sum_j b_j e^{ij phi}|^2 with b real; mode b >= 0 (b = v^2, as the source says) with Nelder Mead
# from many starts then a polish, and mode general (b real, a_k >= 0 imposed) with SLSQP; the best of the starts is reported and the spread
# of the starts is logged. Derivative rows: alternating linear programs (q for fixed p, p for fixed q) inside a bisection on y, and a
# bounded scalar search over x; started from the single level optimum so that the result is at least the single level one.
# Error control. Every closed form row is recomputed with mpmath at 30 digits (root of the defining equation at the optimised parameters)
# and the difference is logged; the optimisation is repeated with two seeds and the difference of the optima logged; the derivative rows
# are repeated on a finer constraint grid. Lemma cf:lem:stechkin is tested on a random grid (minimum of F_theta - kappa F_theta1 over 1e6
# points). Output: data/cf_arith.json and the log data/cf_arith.log (or --out FILE and FILE with .log); --quick takes seconds.
import sys, json, math, time, argparse
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq, linprog
import mpmath as mp

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default='data/cf_arith.json')
ap.add_argument('--skip-deriv', action='store_true')
args = ap.parse_args()
QUICK = args.quick
LOGF = open(args.out.rsplit('.', 1)[0] + '.log', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOGF.write(s + '\n'); LOGF.flush()
T0 = time.time()
RES = {'quick': QUICK, 'rows': {}}
EXTRA_STARTS = []

def coeffs(b):
    K = len(b) - 1
    a = np.zeros(K + 1)
    a[0] = np.dot(b, b)
    for k in range(1, K + 1):
        a[k] = 2 * np.dot(b[:-k], b[k:])
    return a

def cap_new(a):
    a0, a1, A = a[0], a[1], a[1:].sum()
    return 2 / A * (a1 - a0 - a0 * math.log(a1 / a0)) if a1 > a0 else 0.0
def cap_old(a):
    return 2 * (a[1] - a[0]) / a[1:].sum()

# ---- y* by Corollary cf:cor:ystar for a given kernel
def ystar_gen(a, I, J, phi0):
    """J callable; root of a1 J(y) = a0 I + A/2 phi0, zero when no positive root."""
    a0, a1, A = a[0], a[1], a[1:].sum()
    if not (a1 * J(0.0) > a0 * I + A / 2 * phi0):
        return 0.0
    f = lambda y: a1 * J(y) - a0 * I - A / 2 * phi0
    hi = 1.0
    while f(hi) > 0: hi *= 2
    return brentq(f, 0.0, hi, xtol=1e-15, rtol=1e-14)

def ystar_single(a, x): return ystar_gen(a, 1 / x, lambda y: 1 / (x + y), 1.0)
def ystar_two(a, x): return ystar_gen(a, 2 / x, lambda y: 1 / (x + y) + x / (x + y) ** 2, 1.0)
def ystar_mp(a, x, kind, dps=30):
    mp.mp.dps = dps
    a0, a1, A = mp.mpf(a[0]), mp.mpf(a[1]), mp.mpf(sum(a[1:])); x = mp.mpf(x)
    if kind == 'single':
        I = 1 / x; J = lambda y: 1 / (x + y)
    else:
        I = 2 / x; J = lambda y: 1 / (x + y) + x / (x + y) ** 2
    f = lambda y: a1 * J(y) - a0 * I - A / 2
    if f(mp.mpf(0)) <= 0: return mp.mpf(0)
    return mp.findroot(f, mp.mpf(0.05) , tol=mp.mpf(10) ** (-dps + 5))

# ---- optimisation of the polynomial
def single_closed(a):
    a0, a1, A = a[0], a[1], a[1:].sum()
    return 2 * (math.sqrt(a1) - math.sqrt(a0)) ** 2 / A if a1 > a0 else -(a0 - a1) / A  # negative branch only guides the search

def optimise_poly(obj, K, nstart, seed, general=False):
    """maximise obj(a) over p = |sum b_j e^{ij phi}|^2, b >= 0 (b = v^2), or b real with a_k >= 0 (general)."""
    rng = np.random.default_rng(seed)
    best = (-1, None); vals = []
    for s in range(nstart):
        v0 = rng.random(K + 1) + 0.05
        if general:
            f = lambda b: -obj(coeffs(b))
            cons = [{'type': 'ineq', 'fun': (lambda b, k=k: coeffs(b)[k])} for k in range(1, K + 1)]
            b0 = v0 * np.where(rng.random(K + 1) < 0.2, -1, 1)
            r = minimize(f, b0 / np.linalg.norm(b0), method='SLSQP', constraints=cons, options={'ftol': 1e-14, 'maxiter': 400})
            b = r.x
            if min(coeffs(b)[1:]) < -1e-12: continue
        else:
            f = lambda v: -obj(coeffs(v ** 2))
            r = minimize(f, v0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 6000, 'maxfev': 6000})
            r = minimize(f, r.x, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-16, 'maxiter': 6000, 'maxfev': 6000})
            b = r.x ** 2
        val = obj(coeffs(b)); vals.append(val)
        if val > best[0]: best = (val, coeffs(b))
    a = best[1] / best[1][0]
    return best[0], a, (min(vals), max(vals))

nst = 6 if QUICK else 30
degs_single = [2, 3, 4] if QUICK else [2, 3, 4, 5, 6, 7, 8]
degs_two = [3, 4]

log('cf_arith: quick' if QUICK else 'cf_arith: full')
# ================= R18: y* = 0 when there is no root
a345 = np.array([3.0, 4.0, 1.0])
log('--- R18, Corollary cf:cor:ystar: no root cases')
r18 = {}
r18['3+4cos+cos2, delta_1: (a1-a0) I = 1 vs (A/2) phi(0) = 2.5: y*'] = ystar_single(a345, 1.0)
r18['3+4cos+cos2, delta_x, x = 0.186: y*'] = ystar_single(a345, 0.186)
r18['p = 1 (a1 = 0), delta_1: y*'] = ystar_gen(np.array([1.0, 0.0]), 1.0, lambda y: 1 / (1 + y), 1.0)
r18['3+4cos+cos2, curvature kernel x = 1: y*'] = ystar_two(a345, 1.0)
for k, v in r18.items(): log('  ', k, '=', v)
RES['r18'] = r18

# ================= single level rows
log('--- single level: y* = 2 (sqrt a1 - sqrt a0)^2/A, max over x of a1/(a0/x + A/2) - x')
rows = {}
a = a345
x_opt = 2 * (math.sqrt(a[0] * a[1]) - a[0]) / a[1:].sum()
xg = np.exp(np.linspace(-6, 1, 400)); xb = xg[np.argmax([ystar_single(a, x) for x in xg])]
sc = minimize_scalar(lambda x: -ystar_single(a, x), bounds=(xb / 1.1, xb * 1.1), method='bounded', options={'xatol': 1e-12})
ymp = ystar_mp(a, x_opt, 'single')
log('3+4cos+cos2: closed form', single_closed(a), 'at x =', x_opt, '; numerical max over x', -sc.fun, 'at', sc.x, '; mpmath 30 digits', mp.nstr(ymp, 20),
    '; difference', abs(float(ymp) - single_closed(a)))
rows['single, 3+4cos+cos2'] = {'a': list(a), 'y': single_closed(a), 'x': x_opt, 'cap_new': cap_new(a), 'cap_old': cap_old(a), 'mp_diff': abs(float(ymp) - single_closed(a))}
single_best = {}
for K in degs_single:
    v1, a1v, sp = optimise_poly(single_closed, K, nst, 1)
    v2, a2v, _ = optimise_poly(single_closed, K, max(3, nst // 3), 2)
    vg, ag, spg = optimise_poly(single_closed, K, max(3, nst // 3), 3, general=True)
    xo = 2 * (math.sqrt(a1v[0] * a1v[1]) - a1v[0]) / a1v[1:].sum()
    ymp = ystar_mp(a1v, xo, 'single')
    log('single level, best of degree %d: %.6f  (1/c = %.3f) a = %s | seed 2: %.6f | general b, a_k>=0: %.6f | start spread %.6f..%.6f | cap_new %.5f cap_old %.4f | mp diff %.2e' %
        (K, v1, 1 / v1, np.round(a1v, 4).tolist(), v2, vg, sp[0], sp[1], cap_new(a1v), cap_old(a1v), abs(float(ymp) - v1)))
    rows['single, best of degree %d' % K] = {'a': a1v.tolist(), 'y': v1, 'y_seed2': v2, 'y_general': vg, 'x': xo, 'cap_new': cap_new(a1v), 'cap_old': cap_old(a1v), 'spread': sp}
    single_best[K] = (v1, a1v)
RES['rows'].update(rows)

# ================= Stechkin rows
log('--- Stechkin (Proposition cf:prop:stechkin): constant = single level optimum / (1 - 1/sqrt5)')
ks = 1 / math.sqrt(5)
st = {}
st['3+4cos+cos2'] = single_closed(a345) / (1 - ks)
for K in [2, 3, 4]:
    st['best of degree %d' % K] = single_best[K][0] / (1 - ks)
for k, v in st.items(): log('  Stechkin', k, '= %.6f   1/c = %.5f' % (v, 1 / v))
RES['stechkin'] = st
for K in [4] + ([5, 6, 8] if not QUICK else []):
    v = single_best[K][0] / (1 - ks)
    log('  Stechkin degree %d: 1/c = %.6f (Rosser and Schoenfeld R = 9.645908801, recalled, not in the repo)' % (K, 1 / v))
    RES['stechkin']['inverse_deg%d' % K] = 1 / v
a4 = single_best[4][1]
xs = 2 * (math.sqrt(a4[0] * a4[1]) - a4[0]) / ((1 - ks) * a4[1:].sum())
log('  degree 4: a =', np.round(a4, 4).tolist(), ' x (Stechkin) = %.5f' % xs, ' caps: new %.5f old %.5f' % (cap_new(a4), cap_old(a4)))
RES['stechkin']['a4'] = a4.tolist(); RES['stechkin']['x4'] = xs
# Lemma cf:lem:stechkin: F_theta(b, tau) >= kappa F_theta1(b, tau)
rng = np.random.default_rng(7)
nn = 200000 if QUICK else 2000000
worst = 1e9
for theta in [1.0000001, 1.001, 1.01, 1.1, 1.5, 3.0, 10.0, 100.0]:
    th1 = 0.5 * (1 + math.sqrt(1 + 4 * theta ** 2)); kap = (2 * theta - 1) / (2 * th1 - 1)
    b = rng.random(nn // 8); tau = np.exp(rng.uniform(-8, 8, nn // 8)) * rng.choice([-1, 1], nn // 8)
    F = lambda T, bb, t: (1 / ((T - bb) + 1j * t) + 1 / ((T - 1 + bb) + 1j * t)).real
    d = (F(theta, b, tau) - kap * F(th1, b, tau)) / F(theta, b, tau)
    worst = min(worst, d.min())
log('  Lemma cf:lem:stechkin: min over %d random points (theta in 1..100, b in [0,1], |tau| in e^-8..e^8) of (F - kappa F1)/F = %.4g (must be >= 0)' % (nn, worst))
RES['lemma_stechkin_min_rel'] = worst

# ================= two levels, curvature limit
log('--- two levels, curvature limit: phi = exp(-x l)(1 + x l)')
def two_obj_factory():
    def best_x(a):
        lg = np.linspace(-6, 3, 60); vals = [ystar_two(a, math.exp(l)) for l in lg]; i = int(np.argmax(vals))
        if vals[i] <= 0: return 0.0, math.exp(lg[i])
        r = minimize_scalar(lambda lx: -ystar_two(a, math.exp(lx)), bounds=(lg[max(i - 1, 0)], lg[min(i + 1, len(lg) - 1)]), method='bounded', options={'xatol': 1e-10})
        return -r.fun, math.exp(r.x)
    return best_x
bx = two_obj_factory()
def two_total(a):
    return bx(a)[0]
def two_rows(label, K, startpoly=None):
    rng = np.random.default_rng(11)
    best = (-1, None)
    for s in range(max(4, nst // 3)):
        v0 = (rng.random(K + 1) + 0.05) if startpoly is None or s > 0 else np.sqrt(np.sqrt(np.abs(np.array(startpoly))))
        f = lambda v: -two_total(coeffs(v ** 2))
        r = minimize(f, v0, method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 400 if QUICK else 1500, 'maxfev': 400 if QUICK else 1500})
        if -r.fun > best[0]: best = (-r.fun, coeffs(r.x ** 2))
    a = best[1] / best[1][0]; y, x = bx(a)
    return y, a, x
y, x = bx(a345)
ymp = ystar_mp(a345, x, 'two')
log('3+4cos+cos2: y* = %.6f at x = %.5f (1/c = %.2f); cap_new %.4f cap_old %.4f; mp diff %.2e' % (y, x, 1 / y, cap_new(a345), cap_old(a345), abs(float(ymp) - y)))
RES['rows']['two levels, 3+4cos+cos2'] = {'y': y, 'x': x, 'a': list(a345), 'cap_new': cap_new(a345), 'cap_old': cap_old(a345)}
for K in degs_two:
    y, a, x = two_rows('two', K, single_best[K][1])
    ymp = ystar_mp(a, x, 'two')
    log('two levels, best of degree %d: y* = %.6f (1/c = %.2f) at x = %.5f, a = %s; cap_new %.5f cap_old %.4f; mp diff %.2e' %
        (K, y, 1 / y, x, np.round(a, 4).tolist(), cap_new(a), cap_old(a), abs(float(ymp) - y)))
    RES['rows']['two levels, best of degree %d' % K] = {'y': y, 'x': x, 'a': a.tolist(), 'cap_new': cap_new(a), 'cap_old': cap_old(a)}
    EXTRA_STARTS.append(a.tolist())
EXTRA_STARTS.append([1.0, 1.5, 0.8, 0.3, 0.05])

# ================= level with derivatives (alternating LPs)
def make_grids(fine):
    f = int(fine)   # 0 coarse, 1 fine, 2 extra fine
    ns = (150, 300, 1500)[f]
    s = np.linspace(0, 80 if f == 2 else 60, ns)
    t = np.concatenate([[0], np.logspace(-3, 5 if f == 2 else 4, (200, 400, 3000)[f])])
    ph = np.linspace(0, math.pi, (400, 800, 800)[f])
    return s, t, ph

def audit_c(c, npts=2000000):
    """dense check of a computed q: min over s in [0, 400] of q(s) e^{-s}-free polynomial (relative to q(0) = 1) and min over t of the
    positive definiteness function R(t) = sum_j c_j j! Re(1+it)^(-j-1) relative to R(0) = sum c_j j!."""
    d = len(c) - 1
    s = np.concatenate([np.linspace(0, 400, npts), np.logspace(-6, 6, npts)]); q = np.polyval(c[::-1], s)
    t = np.concatenate([[0.0], np.logspace(-6, 6, npts)]); z = 1 + 1j * t
    R = sum(c[j] * math.factorial(j) * (z ** (-(j + 1))).real for j in range(d + 1))
    R0 = sum(c[j] * math.factorial(j) for j in range(d + 1))
    return float(q.min()), float((R * (1 + t * t) ** 0.5).min() / R0)

def deriv_row(d, a_start, fine=False, xs=None, nx=22, fixed_x=None):
    s, t, ph = make_grids(fine)
    fact = [math.factorial(j) for j in range(d + 1)]
    th = np.arctan(t); r = np.sqrt(1 + t * t)
    PDrow = np.array([[fact[j] * math.cos((j + 1) * th[i]) / r[i] ** (j + 1) * r[i] for j in range(d + 1)] for i in range(len(t))])
    Pos = np.array([[s[i] ** j for j in range(d + 1)] for i in range(len(s))])
    Cos = np.array([[math.cos(k * p) for k in range(len(a_start))] for p in ph])
    K = len(a_start) - 1
    def S(eta):
        return np.array([fact[j] / (1 + eta) ** (j + 1) for j in range(d + 1)]), np.array([float(fact[j]) for j in range(d + 1)])
    def lp_c(a, x, y):
        eta = y / x; S1, S0 = S(eta)
        cost = -(a[1] * S1 - a[0] * S0)
        Aub = np.vstack([-Pos, -PDrow]); bub = np.zeros(Aub.shape[0])
        bounds = [(1, 1)] + [(-1e5, 1e5)] * (d - 1) + ([(0, 1e5)] if d >= 1 else [])
        if d == 0: bounds = [(1, 1)]
        r_ = linprog(cost, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs')
        if r_.status != 0: return None, -1e9
        c = r_.x
        return c, float(a[1] * S1 @ c - a[0] * S0 @ c - a[1:].sum() * x / 2 * c[0])
    def lp_a(c, x, y):
        eta = y / x; S1, S0 = S(eta)
        s1, s0 = float(S1 @ c), float(S0 @ c)
        cost = np.zeros(K + 1); cost[1] = -s1; cost[0] = s0; cost[1:] += x / 2 * c[0]  # minimise -(a1 s1 - a0 s0 - A x c0/2)
        r_ = linprog(cost, A_ub=-Cos, b_ub=np.zeros(len(ph)), A_eq=np.eye(1, K + 1), b_eq=[1.0], bounds=[(None, None)] + [(0, 4)] * K, method='highs')
        if r_.status != 0: return None, -1e9
        a = r_.x
        return a, float(a[1] * s1 - a[0] * s0 - a[1:].sum() * x / 2 * c[0])
    def certify(x, y, a0poly, iters=10):
        a = a0poly.copy(); c = None; best = -1e9
        for it in range(iters):
            c, g = lp_c(a, x, y)
            if c is None: return False, a, c
            if g >= -1e-10: return True, a, c
            a2, g2 = lp_a(c, x, y)
            if a2 is None: return False, a, c
            if g2 >= -1e-10: return True, a2, c
            if g2 <= best + 1e-12: break
            best = g2; a = a2
        return False, a, c
    def ystar_x(x, a0poly, tol=1e-6):
        cands = [a0poly] + [np.asarray(z, float) for z in EXTRA_STARTS if len(z) == len(a0poly)]
        outs = [ystar_x1(x, z, tol) for z in cands]
        return max(outs, key=lambda o: o[0])
    def ystar_x1(x, a0poly, tol=1e-6):
        lo, hi = 0.0, 0.6
        ok, a, c = certify(x, hi, a0poly)
        if ok: return hi, a, c
        bestac = (a0poly, None)
        for _ in range(16 if QUICK else 22):
            mid = 0.5 * (lo + hi)
            ok, a, c = certify(x, mid, a0poly)
            if ok: lo = mid; bestac = (a, c); a0poly = a
            else: hi = mid
            if hi - lo < tol: break
        return lo, bestac[0], bestac[1]
    if fixed_x is not None:
        return (fixed_x,) + ystar_x(fixed_x, a_start)
    xs_ = np.linspace(0.55, 1.2, nx) if xs is None else xs
    best = (-1, None, None, None)
    for x in xs_:
        y, a, c = ystar_x(x, a_start)
        if y > best[0]: best = (y, x, a, c)
    # refine x
    y0, x0 = best[0], best[1]; step = (xs_[1] - xs_[0]) if len(xs_) > 1 else 0.05
    r = minimize_scalar(lambda x: -ystar_x(x, best[2])[0], bounds=(max(0.2, x0 - step), x0 + step), method='bounded', options={'xatol': 2e-3, 'maxiter': 12})
    y1, a1, c1 = ystar_x(r.x, best[2])
    if y1 > best[0]: best = (y1, r.x, a1, c1)
    return best[1], best[0], best[2], best[3]

if not args.skip_deriv:
    log('--- level with derivatives (alternating linear programs, constraints on grids)')
    a4s = single_best[4][1]
    dr = {}
    dlist = [2, 3] if QUICK else [2, 3, 4, 5, 6]
    for d in dlist:
        if QUICK:
            x, y, a, c = deriv_row(d, a4s, fixed_x=0.88)[:4]
        else:
            x, y, a, c = deriv_row(d, a4s)
        log('level with derivatives, deg q = %d, best of degree 4: y* = %.5f (1/c = %.2f) at x = %.4f; a = %s; q coefficients (c_j = q_j/x^j) = %s; cap_new %.5f; t=%.0fs' %
            (d, y, 1 / y, x, np.round(a / a[0], 4).tolist(), np.round(c, 4).tolist(), cap_new(a), time.time() - T0))
        dr[d] = {'y': y, 'x': float(x), 'a': a.tolist(), 'c': c.tolist(), 'cap_new': cap_new(a), 'cap_old': cap_old(a)}
        if not QUICK:
            mq, mR = audit_c(np.asarray(c))
            log('   dense audit of the computed q (4e6 points): min q(s) = %.3e (q(0) = 1), min R(t) sqrt(1+t^2)/R(0) = %.3e (both must be >= 0)' % (mq, mR))
            dr[d]['audit'] = [mq, mR]
            xf, yf, af, cf_ = deriv_row(d, a4s, fine=2, fixed_x=x)[:4]
            mq2, mR2 = audit_c(np.asarray(cf_))
            log('   extra fine constraint grids (1500 s points to 80, 3000 t points to 1e5) at the same x: y* = %.5f (difference %.2e); audit %.2e %.2e' % (yf, abs(yf - y), mq2, mR2))
            dr[d]['y_fine'] = yf; dr[d]['audit_fine'] = [mq2, mR2]
            if d >= 4:
                xs2 = np.linspace(0.5, 1.3, 9)
                x2, y2, a2_, c2_ = deriv_row(d, a4s, fine=2, xs=xs2)[:4]
                mq3, mR3 = audit_c(np.asarray(c2_))
                log('   extra fine grids, scan over x in 0.5..1.3 then refined: y* = %.5f at x = %.4f; audit %.2e %.2e' % (y2, x2, mq3, mR3))
                dr[d]['y_fine_scan'] = y2; dr[d]['x_fine_scan'] = float(x2)
    RES['deriv'] = dr

# ================= caps and y* <= cap check
bad = [k for k, v in RES['rows'].items() if v['y'] > v['cap_new'] + 1e-12]
log('--- check y* <= cap_new in every row: violations', bad)
RES['violations'] = bad
RES['seconds'] = time.time() - T0
json.dump(RES, open(args.out, 'w'), indent=1, default=float)
log('done in %.1f s' % (time.time() - T0))
