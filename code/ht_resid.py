# ht_resid.py: the numerical claims of Nick Polson's residue note (R, Chapter 24 material) about the residues w_n = 1/F'(-gamma_n^2) of
# 1/F, F(s) = xi(1/2 + sqrt s)/xi(1/2), the density f of the reciprocal Xi law as the 120 term residue sum, the truncation tables, the
# Taylor coefficients c_n of F and their Toeplitz minors, the explicit truncation constants c, gamma_120, gamma_121, and the zero z_* of the
# continuation of the 120 term sum. R:56, R:99 (not here), R:182 (tab:f), R:190 (mass, mean), R:191 and R:203 to 209 (tab:trunc), R:192
# (log f)'' range, R:251 (w_1, w_2), R:318 to 327 (z_*, the zeros in the box), R:363 (c, gamma_120, gamma_121), R:426 (c_1, c_2, c_3), R:446
# (F(1), Toeplitz minima, 400 random minors), R:473 (k <= floor(pi H) - 1), R:737 (gamma_1^2).
# Conventions. xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s); zeros 1/2 + i gamma_n from data/zeros_1700.txt (68 digits), polished by Newton
# on zeta at the working precision; s_n = -gamma_n^2 = -lambda_n; w_n = 1/F'(s_n) = 2 i gamma_n xi(1/2)/xi'(1/2 + i gamma_n) (real).
# f(y) = sum_{n<=120} w_n exp(-lambda_n y) (the truncated residue sum; the explicit remainder constants M_j are in ht_lcert.py).
# c_n = int Phi(u) u^(2n) du /((2n)! int Phi du), Phi(u) = 2 sum_m (2 pi^2 m^4 e^(9u/2) - 3 pi m^2 e^(5u/2)) exp(-pi m^2 e^(2u)) (Titchmarsh 10.1)
# by the trapezoid rule of step 1/100 on [-3, 3] (the strip of analyticity has half width pi/4, so the error is below exp(-4.9 * 100)).
# Method. mpmath at dps D; every family is also computed at dps D + 30 (and the Polya sums with step 1/80) and the maximum relative difference
# is logged and written: that difference is the error of the printed numbers. Minors are normalised two ways because the chapter does not
# say how: by c_n^k (diag) and by the product over rows of the largest entry of the row (rowmax).
# Options: --quick (dps 40/60, 30 zeros, small boxes, seconds), --out FILE (json; the log goes next to it with the extension .log).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, time, argparse, random
from mpmath import mp, mpf, mpc, zeta, gamma, pi, exp, log, sqrt, findroot, matrix, det, arg, fabs, nstr, binomial, atan, mpmathify

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default='data/ht_resid.json')
a = ap.parse_args()
LOG = open(a.out.replace('.json', '.log'), 'w')
def lg(*x):
    s = ' '.join(str(t) for t in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
R = {}
D = 40 if a.quick else 110
NZ = 30 if a.quick else 121
t0 = time.time()
import os
zfile = [l.strip() for l in open('data/zeros_1700.txt' if os.path.exists('data/zeros_1700.txt') else 'zeros_1700.txt') if l.strip()]

def xi(s): return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)
def dxi(r): return r * (r - 1) / 2 * pi ** (-r / 2) * gamma(r / 2) * zeta(r, 1, 1)   # xi'(rho) at a zero of zeta

def weights(dps, nz):
    mp.dps = dps
    gam, w = [], []
    h = mpf(1) / 2
    x0 = xi(h)
    for n in range(nz):
        r = mpc(h, mpmathify(zfile[n]))
        for _ in range(3):
            r = r - zeta(r) / zeta(r, 1, 1)
        g = r.imag
        wn = 2j * g * x0 / dxi(r)
        gam.append(g); w.append(wn.real)
        assert abs(wn.imag) < mpf(10) ** (-dps // 2) * abs(wn.real), (n, wn)
    return gam, w, x0

gam, w, x0 = weights(D, NZ)
gam2, w2, x02 = weights(D + 30, NZ)
mp.dps = D
dw = max(abs(w[i] / w2[i] - 1) for i in range(NZ)); dg = max(abs(gam[i] - gam2[i]) for i in range(NZ))
lg('weights: %d zeros, dps %d vs %d, max rel diff of w %s, max abs diff of gamma %s, time %.1f' % (NZ, D, D + 30, nstr(dw, 3), nstr(dg, 3), time.time() - t0))
R['weights_reldiff'] = nstr(dw, 3)
lam = [g * g for g in gam]
S = lambda n: n - 1
R['gamma_1_sq'] = nstr(lam[0], 8)
lg('N31 gamma_1^2 =', nstr(lam[0], 8))
lg('N9 w_1 =', nstr(w[0], 8), ' w_2 =', nstr(w[1], 8), ' w_5 =', nstr(w[4], 8))
R['w1'] = nstr(w[0], 12); R['w2'] = nstr(w[1], 12); R['w5'] = nstr(w[4], 12)
if NZ >= 120:
    lg('N1 w_120 =', nstr(w[119], 8), ' log|w_120| =', nstr(log(abs(w[119])), 8))
    R['w120'] = nstr(w[119], 10); R['logabsw120'] = nstr(log(abs(w[119])), 8)
    c = (lam[119] + lam[120]) / 2
    lg('N14 gamma_120 =', nstr(gam[119], 8), ' gamma_121 =', nstr(gam[120], 8), ' c =', nstr(c, 10))
    R['gamma120'] = nstr(gam[119], 10); R['gamma121'] = nstr(gam[120], 10); R['c'] = nstr(c, 12)
lg('N21 floor(pi H) - 1 for H = 3e12:', int(pi * 3 * 10 ** 12) - 1, ' pi H =', nstr(pi * 3 * 10 ** 12, 10))
R['rungs'] = str(int(pi * 3 * 10 ** 12) - 1)

# ---- Taylor coefficients c_n by the Polya integral
def polya_cn(dps, hinv, nmax):
    mp.dps = dps
    h = mpf(1) / hinv
    J = 3 * hinv
    Phi = []
    for j in range(J + 1):
        u = j * h; e2 = exp(2 * u); s = mpf(0); m = 1
        while True:
            t = (2 * pi ** 2 * m ** 4 * exp(9 * u / 2) - 3 * pi * m ** 2 * exp(5 * u / 2)) * exp(-pi * m * m * e2)
            s += t
            if m > 2 and abs(t) < mpf(10) ** (-dps - 30): break
            m += 1
        Phi.append(2 * s)
    def mom(k):
        tot = Phi[0] * (1 if k == 0 else 0) if True else 0
        tot = (Phi[0] if k == 0 else mpf(0)) + 2 * sum(Phi[j] * (j * h) ** k for j in range(1, J + 1))
        return tot * h
    m0 = mom(0)
    cs = [mom(2 * n) / (m0 * mp.factorial(2 * n)) for n in range(nmax + 1)]
    return m0, cs
NC = 36
m0, cs = polya_cn(70, 100, NC)
m0b, csb = polya_cn(90, 80, NC)
mp.dps = 70
dc = max(abs(cs[n] / csb[n] - 1) for n in range(NC + 1))
lg('Polya: xi(1/2) =', nstr(m0, 15), ' closed form', nstr(xi(mpf(1) / 2), 15), ' c_n max rel diff between (70 digits, h 1/100) and (90, 1/80):', nstr(dc, 3))
R['polya_reldiff'] = nstr(dc, 3)
lg('N19 c_1 =', nstr(cs[1], 8), ' c_2 =', nstr(cs[2], 8), ' c_3 =', nstr(cs[3], 8))
R['c1'] = nstr(cs[1], 12); R['c2'] = nstr(cs[2], 12); R['c3'] = nstr(cs[3], 12)
S120 = sum(1 / l for l in lam[:120]) if NZ >= 120 else sum(1 / l for l in lam)
lg('N4 sum_{n<=%d} gamma^-2 = %s ; c_1 - that = %s' % (min(NZ, 120), nstr(S120, 8), nstr(cs[1] - S120, 8)))
R['S120'] = nstr(S120, 10)
F1 = sum(cs); lg('N20 F(1) = sum c_n (n<=36) =', nstr(F1, 12), ' xi(3/2)/xi(1/2) =', nstr(xi(mpf(3) / 2) / xi(mpf(1) / 2), 12))
R['F1'] = nstr(F1, 12)

def minors(cs, N, ks, off):
    out = {}
    for k in ks:
        best_d, best_r, nd, nr = None, None, None, None
        for n in range(k - 1, N - k + 1 - off):
            M = matrix(k, k)
            for i in range(k):
                for j in range(k):
                    M[i, j] = cs[n + i - j]
            d = det(M)
            dd = d / cs[n] ** k
            rm = mpf(1)
            for i in range(k): rm *= max(M[i, j] for j in range(k))
            rr = d / rm
            if best_d is None or dd < best_d: best_d, nd = dd, n
            if best_r is None or rr < best_r: best_r, nr = rr, n
        out[k] = (best_d, nd, best_r, nr)
    return out
for off in (0, 1):
    mn = minors(cs, NC, range(2, 8), off)
    for k, (bd, nd, br, nr) in mn.items():
        lg("N20 consecutive Toeplitz minors order %d (n from k-1 to %d): min diag-normalised %s at n=%d, min rowmax-normalised %s at n=%d" % (k, 36 - k - off, nstr(bd, 3), nd, nstr(br, 3), nr))
        R["toep_%d_%d" % (k, off)] = [nstr(bd, 3), nd, nstr(br, 3), nr]
random.seed(1)
negs = 0; worst = mpf(1); tot = 0
for k in (3, 4, 5):
    for _ in range(0 if False else 133):
        rows = sorted(random.sample(range(0, 40), k)); cols = sorted(random.sample(range(0, 40), k))
        M = matrix(k, k)
        for i in range(k):
            for j in range(k):
                idx = rows[i] - cols[j]
                M[i, j] = cs[idx] if 0 <= idx <= NC else mpf(0)
        d = det(M); rm = mpf(1)
        for i in range(k):
            rm *= max(abs(M[i, j]) for j in range(k))
        v = d / rm if rm != 0 else mpf(0)
        tot += 1; worst = min(worst, v)
        if v < -mpf(10) ** (-40): negs += 1
lg('N20 random minors (orders 3,4,5, %d of them, rows/cols from 0..39): negative beyond -1e-40: %d, smallest normalised %s' % (tot, negs, nstr(worst, 3)))
R['random_minors'] = [tot, negs, nstr(worst, 3)]

# ---- the density f from the residue sum
mp.dps = D
def fj(y, j, W=w, L=lam, n=None):
    n = n or len(W)
    return sum(W[i] * (-L[i]) ** j * exp(-L[i] * y) for i in range(n))
n120 = min(NZ, 120)
if NZ >= 120:
    ys = ['0.010', '0.015', '0.020', '0.025', '0.030', '0.040', '0.050', '0.060', '0.100']
    tab = []
    for ystr in ys:
        y = mpf(ystr); v = fj(y, 0, n=120)
        mp.dps = D + 30; v2 = sum(w2[i] * exp(-lam[i] * y) for i in range(120)); mp.dps = D
        lg('N3 f(%s) = %s  (second precision rel diff %s)' % (ystr, nstr(v, 8), nstr(abs(v / v2 - 1), 3)))
        tab.append((ystr, nstr(v, 8)))
    R['tab_f'] = tab
    # mass and mean of f restricted to [a, inf), a = 0.008
    aa = mpf('0.008')
    mass = sum(w[i] / lam[i] * exp(-lam[i] * aa) for i in range(120))
    mean = sum(w[i] * exp(-lam[i] * aa) * (aa / lam[i] + 1 / lam[i] ** 2) for i in range(120))
    lg('N4 int_{0.008}^inf f = %s (1 - mass below 0.008 = %s); int_{0.008}^inf y f = %s (c_1 = %s, deficit %s)' % (nstr(mass, 12), nstr(1 - mass, 4), nstr(mean, 12), nstr(cs[1], 12), nstr(cs[1] - mean, 4)))
    R['mass_tail'] = nstr(mass, 14); R['mean_tail'] = nstr(mean, 14)
    # (log f)'' on [0.008, 0.2] step 0.001
    lo, hi, ylo, yhi = None, None, None, None
    for k in range(8, 201):
        y = mpf(k) / 1000
        f0, f1, f2 = fj(y, 0, n=120), fj(y, 1, n=120), fj(y, 2, n=120)
        v = (f0 * f2 - f1 ** 2) / f0 ** 2
        if lo is None or v < lo: lo, ylo = v, y
        if hi is None or v > hi: hi, yhi = v, y
    lg('N6 (log f)\'\' on grid 0.008:0.001:0.2: min %s at y=%s, max %s at y=%s' % (nstr(lo, 5), nstr(ylo, 4), nstr(hi, 5), nstr(yhi, 4)))
    R['logf2'] = [nstr(lo, 5), nstr(ylo, 4), nstr(hi, 5), nstr(yhi, 4)]
    for ystr in ('0.1', '0.15', '0.2'):
        y = mpf(ystr); f0, f1, f2 = fj(y, 0, n=120), fj(y, 1, n=120), fj(y, 2, n=120)
        lg('N6 (log f)\'\'(%s) = %s' % (ystr, nstr((f0 * f2 - f1 ** 2) / f0 ** 2, 5)))

# ---- tab:trunc
if NZ >= 120:
    rows = {}
    for n0 in (1, 2, 5):
        lamn = lam[n0 - 1]
        for N in (10, 30, 60, 120):
            prod = mpf(1)
            for m in range(N):
                if m != n0 - 1: prod *= (1 - lamn / lam[m])
            raw = lamn / prod
            tail = cs[1] - sum(1 / l for l in lam[:N])
            corr = raw * exp(lamn * tail)
            # refined: exact factors for zeros N+1..NZ available (here up to 121) are not enough; use the file's 1700 zeros
            rows[(n0, N)] = (raw, corr)
        lg('N5 zero %d: residue %s; raw/corrected at N=10,30,60,120: %s' % (n0, nstr(w[n0 - 1], 8), '; '.join('%s / %s' % (nstr(rows[(n0, N)][0], 5), nstr(rows[(n0, N)][1], 5)) for N in (10, 30, 60, 120))))
    R['trunc'] = {'%d_%d' % k: [nstr(v[0], 6), nstr(v[1], 6)] for k, v in rows.items()}
    n0, N = 1, 120
    lamn = lam[0]; tail = cs[1] - sum(1 / l for l in lam[:N])
    lg('N5/M1 first zero N=120: tail sum %s (chapter 0.0027993), exp(+lambda tail) = %s (chapter 1.7494), exp(-lambda tail) = %s, raw x exp(+) = %s, raw x exp(-) = %s' % (nstr(tail, 6), nstr(exp(lamn * tail), 6), nstr(exp(-lamn * tail), 6), nstr(rows[(1, 120)][1], 6), nstr(rows[(1, 120)][0] * exp(-lamn * tail), 6)))
    # refinement with the exact factors up to the 1700 zeros of the file (68 digits) and the second order term beyond
    mp.dps = 30
    lam1700 = [mpmathify(z) ** 2 for z in zfile[:1700]]
    for n0 in (1, 2, 5):
        lamn = lam[n0 - 1]
        for N in (120,):
            corr_exact = mpf(1)
            for m in range(N, 1700): corr_exact /= (1 - lamn / lam1700[m])
            s_rest = cs[1] - sum(1 / l for l in lam1700)
            corr_exact *= exp(lamn * s_rest)
            lg('N5 refinement zero %d N=120: correction with exact factors to 1700 zeros and exp tail beyond: %s versus exp-only %s versus w/raw %s' % (n0, nstr(corr_exact, 8), nstr(exp(lamn * (cs[1] - sum(1 / l for l in lam[:N]))), 8), nstr(w[n0 - 1] / rows[(n0, 120)][0], 8)))
    mp.dps = D

# ---- z_*: continuation of the 120 term sum, Newton, box counts
def fz(z, n=None):
    n = n or n120
    return sum(w[i] * exp(-lam[i] * z) for i in range(n))
def fzn(z):
    return fz(z) * exp(lam[0] * z) / w[0]
mp.dps = 50 if not a.quick else 30
zs0 = mpc('0.013831874754506', '0.020001801712034')
try:
    zs = findroot(fzn, zs0, tol=mpf(10) ** (-(mp.dps - 8)), maxsteps=40)
    lg('N12 Newton from the printed z_*: z_* = %s' % nstr(zs, 16))
    R['zstar'] = nstr(zs, 16)
    terms = [abs(w[i] * exp(-lam[i] * zs) * exp(lam[0] * zs) / w[0]) for i in range(n120)]
    lg('N12 |normalised sum at printed z_*| = %s ; largest term modulus %s' % (nstr(abs(fzn(zs0)), 3), nstr(max(terms), 3)))
    R['zstar_terms'] = [nstr(abs(fzn(zs0)), 3), nstr(max(terms), 3)]
except Exception as e:
    lg('N12 Newton failed:', e)
    zs = zs0

def winding(z0, z1, z2, z3, eps=mpf('0.5')):
    # count zeros of fzn inside the rectangle with corners z0 (lower left), z1 (lower right), z2 (upper right), z3 (upper left)
    pts = [z0, z1, z2, z3, z0]
    total = mpf(0); calls = 0
    for a_, b_ in zip(pts[:-1], pts[1:]):
        t = mpf(0); prev = fzn(a_); dt = mpf('1e-3')
        while t < 1:
            t2 = min(mpf(1), t + dt)
            cur = fzn(a_ + (b_ - a_) * t2); calls += 1
            d = arg(cur / prev)
            if abs(d) > eps and dt > mpf('1e-9'):
                dt /= 2; continue
            total += d; prev = cur; t = t2; dt = min(dt * 2, mpf('0.02'))
    return total / (2 * pi), calls
if True:
    r = mpf('1e-3')
    wn, calls = winding(zs - r - 1j * r, zs + r - 1j * r, zs + r + 1j * r, zs - r + 1j * r)
    lg('N12 winding number of the 120 term sum on a box of half width 1e-3 around z_*: %s (calls %d)' % (nstr(wn, 8), calls)); R['wind_small'] = nstr(wn, 8)
    if not a.quick:
        wn, calls = winding(mpc('0.012', '1e-12'), mpc('0.15', '1e-12'), mpc('0.15', '0.05'), mpc('0.012', '0.05'))
        lg('N12 winding number on the box 0.012<=Re z<=0.15, 0<Im z<=0.05: %s (calls %d)' % (nstr(wn, 8), calls)); R['wind_big'] = nstr(wn, 8)
        # the same box for the truncation of the sum to 100 zeros, to see how the count depends on the number of terms
        wn, calls = winding(mpc('0.012', '1e-12'), mpc('0.15', '1e-12'), mpc('0.15', '0.05'), mpc('0.012', '0.05'), eps=mpf('0.5'))
lg('time %.1f s' % (time.time() - t0))
json.dump(R, open(a.out, 'w'), indent=1)
