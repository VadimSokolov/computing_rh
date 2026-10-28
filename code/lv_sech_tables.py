# lv_sech_tables.py: numerical claims of Chapter 29 (sech completeness and the Mobius problem), part 1: Section "Numerical check" of the sech
# chapter (the secant table ||S_N - e^-x||_c, the zeta kernel table ||T_N - e^-x||_c, the Mobius sum M_chi), the first zero of beta (barrier
# claim), the Plancherel identity behind Proposition (segment) of the open problems chapter (the factor 2 in 2 D_N - 1/beta), and the
# Thorin check log xi(1+sqrt s)/xi(1) against the sum over the first 1000 zeros plus a smooth tail. Added for the book in October 2026,
# written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. H_c = L^2((0,inf), x^(2c-1) dx), ||f||_c^2 = int |f|^2 x^(2c-1) dx. sech x = 2 sum_m chi_4(m) e^(-m x), so
# sum_(n<=N) mu(n) chi_4(n) sech(n x) -> 2 e^-x and S_N = (1/2) sum_(n<=N) mu(n) chi_4(n) sech(n x). T_N = sum_(n<=N) mu(n) phi(n x),
# phi(y) = 1/(e^y - 1) - 1/y (phi = -1/2 + y/12 - y^3/720 + ... for y < 0.1). beta = L(s, chi_4).
# Method. The norms are trapezoid sums in u = log x (integrand |f|^2 x^(2c), analytic in a strip, so geometric convergence) over
# u in [-30, log 60]; below e^-30 the error is (value at 0)^2 x^(2c)/(2c) < 1e-26 and is dropped. Above x = 60: S_N - e^-x is below 1e-25;
# for T_N the tail is exact: T_N(x) = -m_N/x + O(e^-x), m_N = sum_(n<=N) mu(n)/n, so the tail of the squared norm is
# m_N^2 60^(2c-2)/(2-2c). Error control: every norm is computed on two grids (steps h and h/2) and the maximum relative difference is logged.
# Segment check: ||e^-x - sum a_n sech(n x)||_c^2 against (1/2pi) int |Gamma(c+it)|^2 |1 - 2 beta(c+it) D_N(c+it)|^2 dt (the Mellin Plancherel
# formula), a_n = mu(n) chi_4(n)/2, N = 8, c = 0.75; beta from mpmath. Compared also with the printed version |D_N - 1/beta|^2 |K_c|^2.
# Thorin check: zero sum = sum_(k<=1000) log(((1/2 + sqrt s)^2 + g_k^2)/(1/4 + g_k^2)) (data/zeros_1700.txt), tail = int_(T0)^inf
# (1/2pi) log(t/2pi) log(((1/2+sqrt s)^2 + t^2)/(1/4 + t^2)) dt with T0 the midpoint of g_1000 and g_1001 (a smooth count, not a bound).
# Output: data/lv_sech_tables.json and data/lv_sech_tables.log. Usage: python3 lv_sech_tables.py [--quick] [--out FILE] [--zeros FILE]
import argparse, json, math, time, sys
import numpy as np
import mpmath as mp

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/lv_sech_tables.json')
ap.add_argument('--zeros', default='data/zeros_1700.txt')
a = ap.parse_args()
t0 = time.time()
logf = open(a.out.replace('.json', '.log'), 'w')
def log(*x):
    s = ' '.join(str(i) for i in x); print(s); logf.write(s + '\n'); logf.flush()
res = {}
NS = [10, 100, 1000, 10000] if not a.quick else [10, 100]
NMAX = max(NS)
# mu and chi4
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
isp = np.ones(NMAX + 1, dtype=bool)
for p in range(2, NMAX + 1):
    if isp[p]:
        isp[2 * p::p] = False
        mu[p::p] *= -1
        if p * p <= NMAX: mu[p * p::p * p] = 0
chi = np.zeros(NMAX + 1, dtype=np.int64); chi[1::4] = 1; chi[3::4] = -1
mc = mu * chi
def sech_arr(y):
    e = np.exp(-np.minimum(y, 700.0)); return 2 * e / (1 + e * e)
def phi(y):
    out = np.empty_like(y); s = y < 0.1
    ys = y[s]; out[s] = -0.5 + ys / 12 - ys**3 / 720 + ys**5 / 30240 - ys**7 / 1209600
    yl = y[~s]; out[~s] = 1 / np.expm1(yl) - 1 / yl
    return out
def grid(h):
    ulo, uhi = -30.0, math.log(60.0); m = int(round((uhi - ulo) / h))
    u = np.linspace(ulo, uhi, m + 1); w = np.full(m + 1, (uhi - ulo) / m); w[0] *= .5; w[-1] *= .5
    return u, w
def errS(N, c, h):
    u, w = grid(h); x = np.exp(u); n = np.nonzero(mc[:N + 1])[0]
    S = np.zeros_like(x)
    for i in range(0, len(x), 500):
        xs = x[i:i + 500]
        S[i:i + 500] = 0.5 * (sech_arr(np.outer(xs, n)) @ mc[n].astype(float))
    return math.sqrt(np.sum(w * (S - np.exp(-x))**2 * x**(2 * c)))
def errT(N, c, h):
    u, w = grid(h); x = np.exp(u); n = np.nonzero(mu[:N + 1])[0]
    T = np.zeros_like(x)
    for i in range(0, len(x), 500):
        xs = x[i:i + 500]
        T[i:i + 500] = phi(np.outer(xs, n).ravel()).reshape(len(xs), len(n)) @ mu[n].astype(float)
    I = np.sum(w * (T - np.exp(-x))**2 * x**(2 * c))
    mN = float(np.sum(mu[1:N + 1] / np.arange(1, N + 1)))
    tail = mN**2 * 60.0**(2 * c - 2) / (2 - 2 * c)
    return math.sqrt(I + tail), math.sqrt(I)
h1 = 0.01 if not a.quick else 0.04
log('# secant table ||S_N - e^-x||_c, grid steps', h1, h1 / 2)
res['secant'] = {}; mxd = 0
for N in NS:
    row = []
    for c in [0.60, 0.75, 1.00, 1.25]:
        e1 = errS(N, c, h1); e2 = errS(N, c, h1 / 2); mxd = max(mxd, abs(e1 - e2) / e1)
        row.append(e2); res['secant'][f'{N},{c}'] = e2
    log('N=%d' % N, ' '.join('%.3e' % v for v in row))
res['secant_maxreldiff'] = mxd; log('secant max relative difference between grids', mxd)
log('# zeta kernel table ||T_N - e^-x||_c (with exact tail beyond 60; second value without tail)')
res['zeta'] = {}; mxd = 0
for N in NS:
    row = []
    for c in [0.75, 0.875, 0.95]:
        e1, _ = errT(N, c, h1); e2, e2n = errT(N, c, h1 / 2); mxd = max(mxd, abs(e1 - e2) / e1)
        row.append((e2, e2n)); res['zeta'][f'{N},{c}'] = e2
    log('N=%d' % N, ' '.join('%.3e (no tail %.3e)' % v for v in row))
res['zeta_maxreldiff'] = mxd; log('zeta max relative difference between grids', mxd)
# Mobius errors of tab:dN
res['mobius_dN'] = {}
for c in [0.60, 0.75, 1.00]:
    row = []
    for N in ([4, 8, 16, 32] if not a.quick else [4, 8]):
        row.append(errS(N, c, h1)); res['mobius_dN'][f'{N},{c}'] = row[-1]
    log('Mobius error ||S_N - e^-x||_c at c=%.2f, N=4,8,16,32:' % c, ' '.join('%.3e' % v for v in row))
# M_chi
M = np.cumsum(mc); r = np.abs(M[1:]) / np.sqrt(np.arange(1, NMAX + 1)); k = int(np.argmax(r)) + 1
log('M_chi(N) at', NS, [int(M[n]) for n in NS]); log('max |M_chi(n)|/sqrt n over n <=', NMAX, '=', r[k - 1], 'at n =', k)
n2 = min(NMAX, 20000)
res['Mchi'] = {'max_ratio': float(r[k - 1]), 'at_n': k, 'M_at_NS': {str(n): int(M[n]) for n in NS}}
# first zero of beta
mp.mp.dps = 25
beta = lambda s: mp.dirichlet(s, [0, 1, 0, -1])
Z = lambda t: mp.re(mp.exp(1j * (mp.arg(mp.gamma(0.25 + 0.5j * t)) * 0 )) * 0)  # unused placeholder
z = mp.findroot(lambda t: mp.re(mp.exp(1j * mp.im(mp.loggamma(mp.mpf(3) / 4 + 0.5j * t)) - 0.5j * t * mp.log(mp.pi)) * beta(mp.mpf(1) / 2 + 1j * t)), 6.0)
log('first zero of beta: 1/2 +', mp.nstr(z, 15), 'i; |beta| there', mp.nstr(abs(beta(0.5 + 1j * z)), 3))
res['beta_zero'] = mp.nstr(z, 15)
# segment check (Plancherel)
N, c = 8, 0.75
mp.mp.dps = 20
u, w = grid(0.01 if not a.quick else 0.04); x = np.exp(u); n = np.nonzero(mc[:N + 1])[0]
S = 0.5 * (sech_arr(np.outer(x, n)) @ mc[n].astype(float))
lhs = float(np.sum(w * (S - np.exp(-x))**2 * x**(2 * c)))
hh = 0.1; ts = np.arange(0, 60 + hh / 2, hh); A = B = 0.0
for i, t in enumerate(ts):
    s = mp.mpf(c) + 1j * t; D = sum(mp.mpf(int(mc[m])) / 2 * mp.mpf(int(m))**(-s) for m in n)
    bt = beta(s); g2 = abs(mp.gamma(s))**2
    wt = hh * (0.5 if i == 0 else 1.0) * 2
    A += float(wt * g2 * abs(1 - 2 * bt * D)**2) / (2 * math.pi)
    B += float(wt * g2 * abs(bt)**2 * abs(D - 1 / bt)**2) / (2 * math.pi)
log('Plancherel check N=8 c=0.75: x-space %.10e; (1/2pi) int |K|^2 |2D-1/beta|^2 %.10e; with printed |D-1/beta|^2 %.10e' % (lhs, A, B))
res['segment'] = {'xspace': lhs, 'plancherel_2D': A, 'printed_D': B}
# Thorin
zs = []
for line in open(a.zeros):
    line = line.strip()
    if line: zs.append(mp.mpf(line))
g = zs[:1000]; T0 = (g[999] + zs[1000]) / 2
res['thorin'] = {}
for s in [1, 4, 16]:
    r_ = mp.sqrt(s); xi = lambda z: z * (z - 1) / 2 * mp.pi**(-z / 2) * mp.gamma(z / 2) * mp.zeta(z)
    ex = mp.log(xi(1 + r_) / mp.mpf(0.5)); zsum = sum(mp.log(((0.5 + r_)**2 + t * t) / (0.25 + t * t)) for t in g)
    tail = mp.quad(lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) * mp.log(((0.5 + r_)**2 + t * t) / (0.25 + t * t)), [T0, 10 * T0, mp.inf])
    log('s=%d: log xi ratio %s, zero sum %s, smooth tail %s, total %s, remainder exact - zero sum %s' % (s, mp.nstr(ex, 8), mp.nstr(zsum, 6), mp.nstr(tail, 4), mp.nstr(zsum + tail, 6), mp.nstr(ex - zsum, 4)))
    res['thorin'][str(s)] = [mp.nstr(v, 10) for v in (ex, zsum, tail, zsum + tail)]
log('elapsed', time.time() - t0)
json.dump(res, open(a.out, 'w'), indent=1)
