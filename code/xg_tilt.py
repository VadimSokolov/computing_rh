# code/xg_tilt.py (October 2026, from misc/reviews/nick-xiggc/: review-A2.md, table (b) rows prop:gauss, thm:thorintilt, tab:tilt, thm:weak, thm:ggcdebruijn and the lessons, items (d)2 to (d)5, (d)7, errors E6 to E9, E14, E15): the numbers of Section sec:xg:tilt ("The Thorin atoms of the tilted clock", the end of Section 8.4), written from the mathematics; the scripts of the authors' chapter were neither read nor run.
# Definitions. Kernel and flow in the normalisation of Section sec:ch21:flow: Phi_dBN(u) = sum_{n>=1} (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u)) exp(-pi n^2 e^(4u)) (all terms positive for u >= 0), H_t(z) = int_0^oo e^(t u^2) Phi_dBN(u) cos(z u) du, H_0(z) = Xi(z/2)/8, R_t(s) = H_t(0)/H_t(2 i sqrt s), a_j(t) = z_j(t)^2/4 for the positive zeros z_j(t) of H_t (the Thorin atoms of R_t for t >= Lambda; a_j(0) = gamma_j^2), M_k(t) = int_R v^(2k) e^(t v^2/4) Phi(v) dv with Phi_dBN(u) = Phi(2u)/2 (Section sec:ch21:hankel), so that M_k/M_0 = 4^k int u^(2k) e^(t u^2) Phi_dBN du / int e^(t u^2) Phi_dBN du. T_t is the variable with Laplace transform R_t. By the Newton identities of Section sec:ch21:hankel, E T_t = sum_j 1/a_j = M_1/(2 M_0) and Var T_t = sum_j a_j^(-2) = (M_1/M_0)^2/4 - (M_2/M_0)/12 (the two sums are the power sums m_0(t), m_1(t) of that section).
# (1) The table of Section sec:xg:tilt. The four smallest zeros z_j of H_t for t = 0 and 0.2 (Newton's method on the integral, started at 2 gamma_j; at t = 0 compared with 2 gamma_j from zeros_8000.txt), the atoms a_j = z_j^2/4, E T_t and Var T_t; at t = 0 compared with the constants p_1 = 0.0231049931154189707889 and p_2 = 0.0000371725992852696862 of Section sec:ch18:moments and with the sums over the 8000 zeros completed by the tail (log(T/2 pi) + 1)/(2 pi T), which is accurate to about 1/T^2. The rates da_j/dt at t = 0 from the dynamics da_j/dt = 1/2 + 2 sum_{k != j} a_j/(a_j - a_k) (the sum over the 8000 zeros a_k = gamma_k^2, plus the tail -2 a_j (log(T/2 pi) + 1)/(2 pi T), T = gamma_8000, which is the first term of sum_{k>8000} 2 a_j/(a_j - a_k) = -2 a_j sum_{k>8000} a_k^(-1) - ...), against the central difference of the zeros at t = +-0.01 and +-0.005 (Richardson) and against the difference quotient (a_j(0.2) - a_j(0))/0.2.
# (2) The mean identity dE T_t/dt = (E T_t)^2 - (3/2) Var T_t at t = 0 and 0.2. Both sides from the integrals: the left side is (M_2/M_0 - (M_1/M_0)^2)/8, because dM_k/dt = M_(k+1)/4, the right side from E T_t and Var T_t above; and the left side from the central difference of E T_t at t = +-0.01, +-0.005 (Richardson), which does not use the formula. The difference quotient of E T over [0, 0.2] against the mean of the two end values.
# (3) The shadow of the clock (the paragraph after the mean identity): R_0.2(s)/R_0(s) - 1 at s = 1, 10, 100 and E T_0.2 - E T_0. For s = 1, 10, 100 also 8 H_0(2 i sqrt s) against xi(1/2 + sqrt s), which checks the quadrature at imaginary arguments.
# (4) The Gaussian representation (Proposition xg:prop:gauss): 8 H_t(-2 i sigma) = E xi(1/2 + sigma + sqrt(t/2) Z), Z standard normal, at t = 0.3 and sigma = 2 (s = sigma^2 = 4), and the same at sigma = 0, with xi(s) = s (s - 1) pi^(-s/2) Gamma(s/2) zeta(s)/2 from mpmath (a different route from the kernel) and a Gauss Hermite rule of 80 and of 120 nodes; the number of agreeing digits is reported.
# (5) The large tilt (Proposition xg:prop:large): for t = 300, 1000, 3000 the mode u_t of e^(t u^2) Phi_dBN(u) (the zero of 2 t u + Phi_dBN'(u)/Phi_dBN(u), with Phi_dBN' from the series), the standard deviation of the density proportional to e^(t u^2) Phi_dBN(u) (integrated over a window around the mode where the integrand exceeds e^(-105) times its maximum), the zeros z_1, z_2, z_3 of H_t and u_t z_j against pi (j - 1/2), the ratio e^(4 u_t) 2 pi/(t u_t) (which tends to 1) and R_t(s/u_t^2) cosh(2 sqrt s) at s = 1. Also that 2 t u + Phi_dBN'/Phi_dBN stays positive on a grid below the mode.
# Error control. Composite Gauss Legendre rules: for t of order 1 the rule on [0, 1.5] (Phi_dBN(1.5) < 1e-540) with 6 pieces of 48 nodes, and a second one with 5 pieces of 60 nodes, whose difference from the first is reported for every quantity of (1) to (3); for (5) a window of pieces of 32 nodes and a second rule of 40 nodes. Nodes and weights are polished by Newton's method at the working precision (40 digits), as in tf_heat.py. Nothing here is a proof.
# Output: xg_tilt.json in the working directory (values as strings, with the parameters and the precision), copied to data/xg_tilt.json; the log is data/xg_tilt.log. Run on Hopper from the book root: bash misc/tools/hopper_run.sh -t 30 -n xg_tilt -g xg_tilt.json code/xg_tilt.py data/zeros_8000.txt (one core, about a minute); python3 xg_tilt.py [zeros file].
import json
import os
import sys
import time

import mpmath as mp
import numpy as np

mp.mp.dps = 40
T0 = time.time()
PREC = 40


def say(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def ns(x, n=15):
    return mp.nstr(x, n)


def agree(a, b):
    """Number of agreeing decimal digits of a and b (relative to b), as a float."""
    d = abs(mp.mpf(a) - mp.mpf(b))
    if d == 0:
        return float(PREC)
    return float(-mp.log10(d / abs(mp.mpf(b))))


def phi_dphi(u):
    """Phi_dBN(u) and Phi_dBN'(u) from the series, to the working precision."""
    u = mp.mpf(u)
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    stop = mp.mpf(10) ** (-(mp.mp.dps + 6))
    tot = mp.mpf(0)
    dtot = mp.mpf(0)
    n = 1
    while True:
        a = 2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5
        da = 18 * mp.pi ** 2 * n ** 4 * e9 - 15 * mp.pi * n ** 2 * e5
        ex = mp.exp(-mp.pi * n * n * e4)
        term = a * ex
        tot += term
        dtot += (da - 4 * mp.pi * n * n * e4 * a) * ex
        if term < stop * tot:
            return tot, dtot
        n += 1


def xi_func(s):
    s = mp.mpf(s)
    if s < mp.mpf(1) / 2:
        s = 1 - s
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


_GL = {}


def gl_nodes(n):
    """All nodes and weights of the n point Gauss Legendre rule on [-1, 1], polished by Newton's method at the working precision."""
    if n in _GL:
        return _GL[n]
    x0 = np.polynomial.legendre.leggauss(n)[0]
    tol = mp.mpf(10) ** (-mp.mp.dps - 4)
    xs, ws = [], []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(12):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < tol:
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    _GL[n] = (xs, ws)
    return _GL[n]


def gh_nodes(n):
    """All nodes and weights of the n point Gauss Hermite rule for the weight exp(-x^2), polished at the working precision."""
    x0 = np.polynomial.hermite.hermgauss(n)[0]
    tol = mp.mpf(10) ** (-mp.mp.dps - 4)
    xs, ws = [], []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(12):
            h0, h1 = mp.mpf(1), 2 * x
            for k in range(1, n):
                h0, h1 = h1, 2 * x * h1 - 2 * k * h0
            dx = h1 / (2 * n * h0)
            x -= dx
            if abs(dx) < tol:
                break
        h0, h1 = mp.mpf(1), 2 * x
        for k in range(1, n):
            h0, h1 = h1, 2 * x * h1 - 2 * k * h0
        xs.append(x)
        ws.append(mp.mpf(2) ** (n - 1) * mp.factorial(n) * mp.sqrt(mp.pi) / (n * n * h0 * h0))
    return xs, ws


def composite(a, b, pieces, n):
    xs, ws = gl_nodes(n)
    a, b = mp.mpf(a), mp.mpf(b)
    h = (b - a) / pieces
    u, W = [], []
    for p in range(pieces):
        lo = a + p * h
        for x, w in zip(xs, ws):
            u.append(lo + h * (1 + x) / 2)
            W.append(h * w / 2)
    return u, W


class Rule:
    """Quadrature of int_0^1.5 f(u) Phi_dBN(u) du with the values of Phi_dBN stored; H_t is computed with the tilt e^(t u^2)."""

    def __init__(self, pieces, n, b=mp.mpf('1.5')):
        self.u, self.W = composite(0, b, pieces, n)
        self.P = [phi_dphi(x)[0] for x in self.u]
        self.cache = {}

    def w(self, t):
        t = mp.mpf(t)
        key = mp.nstr(t, 20)
        if key not in self.cache:
            self.cache[key] = [W * mp.exp(t * x * x) * p for x, W, p in zip(self.u, self.W, self.P)]
        return self.cache[key]

    def H(self, t, z):
        return mp.fsum(wi * mp.cos(z * x) for wi, x in zip(self.w(t), self.u))

    def H_imag(self, t, y):
        """H_t(i y) = int e^(t u^2) Phi_dBN(u) cosh(y u) du."""
        return mp.fsum(wi * mp.cosh(y * x) for wi, x in zip(self.w(t), self.u))

    def ratio(self, t, k):
        """M_k(t)/M_0(t) = 4^k int u^(2k) e^(t u^2) Phi_dBN du / int e^(t u^2) Phi_dBN du."""
        w = self.w(t)
        return 4 ** k * mp.fsum(wi * x ** (2 * k) for wi, x in zip(w, self.u)) / mp.fsum(w)

    def zero(self, t, z0):
        w = self.w(t)
        z = mp.mpf(z0)
        tol = mp.mpf(10) ** (-mp.mp.dps + 5)
        for _ in range(60):
            f = mp.fsum(wi * mp.cos(z * x) for wi, x in zip(w, self.u))
            df = -mp.fsum(wi * x * mp.sin(z * x) for wi, x in zip(w, self.u))
            dz = f / df
            z -= dz
            if abs(dz) < tol:
                break
        return z


def mean_var(rule, t):
    r1, r2 = rule.ratio(t, 1), rule.ratio(t, 2)
    return r1 / 2, r1 ** 2 / 4 - r2 / 12, r1, r2


def read_zeros(path):
    return [mp.mpf(line.strip()) for line in open(path) if line.strip()]


zpath = None
for cand in (sys.argv[1] if len(sys.argv) > 1 else None, 'zeros_8000.txt', 'data/zeros_8000.txt', '../data/zeros_8000.txt'):
    if cand and os.path.exists(cand):
        zpath = cand
        break
if zpath is None:
    sys.exit('zeros_8000.txt not found')
gam = read_zeros(zpath)
NZ = len(gam)
say('zeros read from %s: %d ordinates, gamma_1 = %s, gamma_%d = %s' % (zpath, NZ, ns(gam[0], 20), NZ, ns(gam[-1], 12)))

out = dict(script='xg_tilt.py', working_precision_digits=PREC, zeros_file=os.path.basename(zpath), number_of_zeros=NZ)

# ---------------------------------------------------------------- (1) the table
ruleA = Rule(6, 48)
ruleB = Rule(5, 60)
say('rules built: A has %d nodes, B has %d nodes' % (len(ruleA.u), len(ruleB.u)))

P1 = mp.mpf('0.0231049931154189707889')
P2 = mp.mpf('0.0000371725992852696862')


def tail_M(T):
    return (mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T)


rows = {}
for t in ('0', '0.2'):
    zs, zsB = [], []
    for j in range(4):
        z0 = 2 * gam[j]
        zs.append(ruleA.zero(t, z0))
        zsB.append(ruleB.zero(t, z0))
    atoms = [z * z / 4 for z in zs]
    ET, VT, r1, r2 = mean_var(ruleA, t)
    ETB, VTB, _, _ = mean_var(ruleB, t)
    rows[t] = dict(z=zs, a=atoms, ET=ET, Var=VT)
    say('t = %s: a_1..a_4 = %s' % (t, ', '.join(ns(a, 12) for a in atoms)))
    say('t = %s: E T = %s, Var T = %s   (M_1/M_0 = %s, M_2/M_0 = %s)' % (t, ns(ET, 20), ns(VT, 15), ns(r1, 15), ns(r2, 15)))
    say('t = %s: rule B differs from rule A by %.2e in a_j, %.2e in E T (relative), %.2e in Var T (relative)' % (
        t, float(max(abs(z * z / 4 - a) for z, a in zip(zsB, atoms))), float(abs(ETB - ET) / ET), float(abs(VTB - VT) / VT)))

say('check at t = 0: z_j against 2 gamma_j from the file (digits of agreement): ' + ', '.join('%.1f' % agree(rows['0']['z'][j], 2 * gam[j]) for j in range(4)))
say('check at t = 0: E T against p_1 of Section sec:ch18:moments: %.1f digits; Var T against p_2: %.1f digits' % (agree(rows['0']['ET'], P1), agree(rows['0']['Var'], P2)))
S1 = mp.fsum(g ** -2 for g in gam) + tail_M(gam[-1])
S2 = mp.fsum(g ** -4 for g in gam)
say('check at t = 0: sum_{k<=%d} gamma_k^-2 + tail = %s against E T = %s (absolute difference %.2e, the size of 1/T^2 = %.1e, T = gamma_%d, which is the error of the tail); sum_{k<=%d} gamma_k^-4 = %s (the tail is about 4e-14) against Var T = %s' % (
    NZ, ns(S1, 15), ns(rows['0']['ET'], 15), float(abs(S1 - rows['0']['ET'])), float(1 / gam[-1] ** 2), NZ, NZ, ns(S2, 15), ns(rows['0']['Var'], 15)))

# rates da_j/dt at t = 0 from the dynamics
a_all = [g * g for g in gam]
TM = tail_M(gam[-1])
pred = []
for j in range(4):
    aj = a_all[j]
    s = mp.fsum(aj / (aj - a_all[k]) for k in range(NZ) if k != j)
    pred.append(mp.mpf(1) / 2 + 2 * s - 2 * aj * TM)
say('predicted rates da_j/dt at t = 0 (sum over %d zeros plus tail): %s' % (NZ, ', '.join(ns(p, 8) for p in pred)))
say('   the tail term alone: %s' % ', '.join(ns(-2 * a_all[j] * TM, 6) for j in range(4)))

# central differences of the zeros
cd = []
for j in range(4):
    D = {}
    for h in ('0.01', '0.005'):
        hp, hm = mp.mpf(h), -mp.mpf(h)
        zp = ruleA.zero(hp, 2 * gam[j])
        zm = ruleA.zero(hm, 2 * gam[j])
        D[h] = (zp * zp / 4 - zm * zm / 4) / (2 * hp)
    cd.append((4 * D['0.005'] - D['0.01']) / 3)
say('central difference (Richardson, h = 0.01 and 0.005) of a_j at t = 0: %s' % ', '.join(ns(c, 8) for c in cd))
say('   difference from the predicted rate: %s' % ', '.join('%.1e' % float(c - p) for c, p in zip(cd, pred)))
obs = [(rows['0.2']['a'][j] - rows['0']['a'][j]) / mp.mpf('0.2') for j in range(4)]
say('difference quotient (a_j(0.2) - a_j(0))/0.2: %s' % ', '.join(ns(o, 8) for o in obs))

out['table'] = {t: dict(zeros_z=[ns(z, 25) for z in r['z']], atoms=[ns(a, 20) for a in r['a']], E_T=ns(r['ET'], 25), Var_T=ns(r['Var'], 20)) for t, r in rows.items()}
out['table_check_t0'] = dict(E_T_digits_vs_p1=agree(rows['0']['ET'], P1), Var_T_digits_vs_p2=agree(rows['0']['Var'], P2),
                             zero_sum_plus_tail=ns(S1, 15), zero_sum_gamma_minus4=ns(S2, 15))
out['rates'] = dict(predicted_at_0=[ns(p, 12) for p in pred], tail_term=[ns(-2 * a_all[j] * TM, 8) for j in range(4)],
                    central_difference_at_0=[ns(c, 12) for c in cd], difference_quotient_0_to_0p2=[ns(o, 12) for o in obs])

# ---------------------------------------------------------------- (2) the mean identity
mid = {}
for t in ('0', '0.2'):
    ET, VT, r1, r2 = mean_var(ruleA, t)
    lhs = (r2 - r1 ** 2) / 8
    rhs = ET ** 2 - mp.mpf(3) / 2 * VT
    mid[t] = (lhs, rhs)
    say('t = %s: dE T/dt from the integrals (M_2/M_0 - (M_1/M_0)^2)/8 = %s ; (E T)^2 - 1.5 Var T = %s ; relative difference %.1e' % (t, ns(lhs, 12), ns(rhs, 12), float(abs(lhs - rhs) / rhs)))
ETs = {}
for h in ('0.01', '0.005'):
    ETs[h] = (mean_var(ruleA, mp.mpf(h))[0] - mean_var(ruleA, -mp.mpf(h))[0]) / (2 * mp.mpf(h))
rich = (4 * ETs['0.005'] - ETs['0.01']) / 3
say('central difference of E T at t = 0 (Richardson): %s against %s (relative difference %.1e)' % (ns(rich, 12), ns(mid['0'][1], 12), float(abs(rich - mid['0'][1]) / mid['0'][1])))
quot = (rows['0.2']['ET'] - rows['0']['ET']) / mp.mpf('0.2')
say('difference quotient of E T over [0, 0.2]: %s ; mean of the two end values: %s' % (ns(quot, 10), ns((mid['0'][1] + mid['0.2'][1]) / 2, 10)))
out['mean_identity'] = dict(lhs_t0=ns(mid['0'][0], 15), rhs_t0=ns(mid['0'][1], 15), lhs_t0p2=ns(mid['0.2'][0], 15), rhs_t0p2=ns(mid['0.2'][1], 15),
                            central_difference_t0=ns(rich, 15), difference_quotient=ns(quot, 12))

# ---------------------------------------------------------------- (3) the shadow
sh = []
for s in (1, 10, 100):
    y = 2 * mp.sqrt(s)
    R0 = ruleA.H(0, 0) / ruleA.H_imag(0, y)
    R2 = ruleA.H(mp.mpf('0.2'), 0) / ruleA.H_imag(mp.mpf('0.2'), y)
    R0B = ruleB.H(0, 0) / ruleB.H_imag(0, y)
    R2B = ruleB.H(mp.mpf('0.2'), 0) / ruleB.H_imag(mp.mpf('0.2'), y)
    rel = R2 / R0 - 1
    relB = R2B / R0B - 1
    xi_route = 8 * ruleA.H_imag(0, y) / xi_func(mp.mpf(1) / 2 + mp.sqrt(s))
    say('s = %d: R_0.2(s)/R_0(s) - 1 = %s (rule B %s); 8 H_0(2 i sqrt s)/xi(1/2 + sqrt s) = %s' % (s, ns(rel, 10), ns(relB, 10), ns(xi_route, 20)))
    sh.append((s, rel, xi_route))
dmean = rows['0.2']['ET'] - rows['0']['ET']
say('E T_0.2 - E T_0 = %s (relative to E T_0: %s)' % (ns(dmean, 10), ns(dmean / rows['0']['ET'], 6)))
out['shadow'] = dict(rel_diff=[dict(s=s, value=ns(r, 12), xi_route_ratio=ns(x, 25)) for s, r, x in sh], mean_difference=ns(dmean, 12))

# ---------------------------------------------------------------- (4) the Gaussian representation
tg = mp.mpf('0.3')
gauss = {}
for nh in (80, 120):
    xs, ws = gh_nodes(nh)
    for sigma in (2, 0):
        avg = mp.fsum(w * xi_func(mp.mpf(1) / 2 + sigma + mp.sqrt(tg) * x) for x, w in zip(xs, ws)) / mp.sqrt(mp.pi)
        gauss[(nh, sigma)] = avg
lhs2 = 8 * ruleA.H_imag(tg, 4)
lhs2B = 8 * ruleB.H_imag(tg, 4)
lhs0 = 8 * ruleA.H(tg, 0)
say('t = 0.3, sigma = 2 (s = 4): 8 H_t(-4 i) = %s (rule B %s)' % (ns(lhs2, 25), ns(lhs2B, 25)))
for nh in (80, 120):
    say('   Gauss Hermite %3d nodes: E xi(1/2 + 2 + sqrt(0.15) Z) = %s, agreement %.1f digits' % (nh, ns(gauss[(nh, 2)], 25), agree(gauss[(nh, 2)], lhs2)))
say('t = 0.3, sigma = 0: 8 H_t(0) = %s' % ns(lhs0, 25))
for nh in (80, 120):
    say('   Gauss Hermite %3d nodes: E xi(1/2 + sqrt(0.15) Z) = %s, agreement %.1f digits' % (nh, ns(gauss[(nh, 0)], 25), agree(gauss[(nh, 0)], lhs0)))
ratio_lhs = ruleA.H_imag(tg, 4) / ruleA.H(tg, 0)
ratio_rhs = gauss[(120, 2)] / gauss[(120, 0)]
say('1/R_0.3(4) = H(2 i 2)/H(0) = %s ; ratio of the two averages = %s ; agreement %.1f digits' % (ns(ratio_lhs, 20), ns(ratio_rhs, 20), agree(ratio_rhs, ratio_lhs)))
out['gaussian'] = dict(t='0.3', sigma=2, lhs_8H=ns(lhs2, 30), average_80=ns(gauss[(80, 2)], 30), average_120=ns(gauss[(120, 2)], 30),
                       digits_80=agree(gauss[(80, 2)], lhs2), digits_120=agree(gauss[(120, 2)], lhs2),
                       sigma0_lhs=ns(lhs0, 30), sigma0_average_120=ns(gauss[(120, 0)], 30), sigma0_digits_120=agree(gauss[(120, 0)], lhs0),
                       ratio_digits=agree(ratio_rhs, ratio_lhs))

# ---------------------------------------------------------------- (5) the large tilt
large = []


def large_tilt(tt):
    t = mp.mpf(tt)
    u = mp.mpf(1)
    for _ in range(80):
        u = mp.log((2 * t * u + 9) / (4 * mp.pi)) / 4

    def gp(x):
        p, dp = phi_dphi(x)
        return 2 * t * x + dp / p

    um = mp.findroot(gp, (u, u + mp.mpf('1e-3')), solver='secant', tol=mp.mpf(10) ** (-30), maxsteps=100)

    def g(x):
        return t * x * x + mp.log(phi_dphi(x)[0])

    gt = g(um)
    g2 = mp.diff(gp, um)
    sig = 1 / mp.sqrt(-g2)
    L = mp.mpf(105)
    step = sig / 2
    a = um
    while a > 0 and g(a) > gt - L:
        a -= step
    a = max(a, mp.mpf(0))
    b = um
    while g(b) > gt - L:
        b += step
    res = {}
    for npts in (32, 40):
        pieces = int(mp.ceil((b - a) / (2 * sig))) + 1
        u_, W_ = composite(a, b, pieces, npts)
        w = [Wi * mp.exp(g(x) - gt) for x, Wi in zip(u_, W_)]
        Z0 = mp.fsum(w)
        mean = mp.fsum(wi * x for wi, x in zip(w, u_)) / Z0
        var = mp.fsum(wi * (x - mean) ** 2 for wi, x in zip(w, u_)) / Z0
        zz = []
        for j in range(1, 4):
            z = (j - mp.mpf(1) / 2) * mp.pi / um
            for _ in range(60):
                f = mp.fsum(wi * mp.cos(z * x) for wi, x in zip(w, u_))
                df = -mp.fsum(wi * x * mp.sin(z * x) for wi, x in zip(w, u_))
                dz = f / df
                z -= dz
                if abs(dz) < mp.mpf(10) ** (-30):
                    break
            zz.append(z)
        Rs = Z0 / mp.fsum(wi * mp.cosh(2 * x / um) for wi, x in zip(w, u_))
        res[npts] = dict(mean=mean, sd=mp.sqrt(var), z=zz, R1=Rs, pieces=pieces)
    grid = [um * mp.mpf(k) / 400 for k in range(1, 400) if um * mp.mpf(k) / 400 < um - 3 * sig]
    gmin = min(gp(x) / t for x in grid)
    return dict(t=t, um=um, sig_laplace=sig, a=a, b=b, res=res, gmin=gmin, ratio4=mp.exp(4 * um) * 2 * mp.pi / (t * um))


for tt in (300, 1000, 3000):
    r = large_tilt(tt)
    r1, r2 = r['res'][32], r['res'][40]
    say('t = %d: mode u_t = %s, e^(4u) 2pi/(t u) = %s, Laplace s.d. %s, s.d. %s (rule 40: %s), mean - u_t = %s' % (
        tt, ns(r['um'], 12), ns(r['ratio4'], 10), ns(r['sig_laplace'], 6), ns(r1['sd'], 10), ns(r2['sd'], 10), ns(r1['mean'] - r['um'], 6)))
    say('   window [%s, %s]; z_1, z_2, z_3 of H_t: %s ; u_t z_j: %s ; against pi (j - 1/2): %s' % (
        ns(r['a'], 8), ns(r['b'], 8), ', '.join(ns(z, 10) for z in r1['z']), ', '.join(ns(r['um'] * z, 10) for z in r1['z']),
        ', '.join(ns(r['um'] * z / (mp.pi * (j + mp.mpf(1) / 2)) - 1, 4) for j, z in enumerate(r1['z']))))
    say('   rule 40 gives u_t z_1 = %s ; R_t(1/u_t^2) cosh(2) = %s ; min of g\'/t on the grid below the mode = %s' % (
        ns(r['um'] * r2['z'][0], 10), ns(r1['R1'] * mp.cosh(2), 10), ns(r['gmin'], 6)))
    large.append(dict(t=tt, mode=ns(r['um'], 15), ratio_e4u=ns(r['ratio4'], 12), sd=ns(r1['sd'], 12), sd_rule40=ns(r2['sd'], 12),
                      mean_minus_mode=ns(r1['mean'] - r['um'], 8), z=[ns(z, 15) for z in r1['z']], u_z=[ns(r['um'] * z, 15) for z in r1['z']],
                      u_z1_rule40=ns(r['um'] * r2['z'][0], 15), R_cosh_ratio=ns(r1['R1'] * mp.cosh(2), 12), min_gprime_over_t=ns(r['gmin'], 8)))
out['large_tilt'] = large
say('pi/2 = %s' % ns(mp.pi / 2, 12))

out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open('xg_tilt.json', 'w'), indent=1)
say('wrote xg_tilt.json in %.1f s' % (time.time() - T0))
