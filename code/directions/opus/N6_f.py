"""Mean field and hierarchical checks; Gaussian damping thresholds. Run: timeout 500 python3 -I N6_f.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
import mpmath as mp
from scipy.special import gammaln
print('F1 Curie Weiss with Ising spins, law of S_N/N^(3/4)')
def cw_moments(N, beta):
    k = np.arange(0, N + 1)
    S = 2.0 * k - N
    lw = gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1) + beta * S * S / (2 * N)
    lw -= lw.max()
    w = np.exp(lw)
    w /= w.sum()
    x = S / N ** 0.75
    m2 = (w * x * x).sum()
    m4 = (w * x ** 4).sum()
    return m2, m4 / m2 ** 2 - 3
g14 = math.gamma(0.25); g34 = math.gamma(0.75)
ex = 0.25 * (g14 / g34) ** 2 - 3
print('  excess of exp(-x^4): %.4f; of the convolution of two: %.4f' % (ex, ex / 2))
for N in (10 ** 4, 10 ** 6, 10 ** 7):
    print('  N=%d critical beta=1: excess %.4f; two independent blocks: %.4f' % (N, cw_moments(N, 1.0)[1], cw_moments(N, 1.0)[1] / 2))
t = 2.0
from scipy.integrate import quad
Mom = [quad(lambda x: x ** p * math.exp(-x ** 4 / 12 - t * x * x / 2), 0, 12)[0] for p in (0, 2, 4)]
exw = (Mom[2] / Mom[0]) / (Mom[1] / Mom[0]) ** 2 - 3
print('  critical window beta=1-%g/sqrt(N): N=10^7 excess %.4f; density exp(-x^4/12-%gx^2/2) has excess %.4f' % (t, cw_moments(10 ** 7, 1 - t / math.sqrt(10 ** 7))[1], t, exw))
print('F2 exact convolution at c=2: (rho*rho)(2x)=(2/lam)K_0(2a cosh lam x)')
a, lam, x = 1.0, 1.0, 0.7
num = quad(lambda v: math.exp(-a * math.cosh(lam * (x + v)) - a * math.cosh(lam * (x - v))), -15, 15)[0]
print('  a=1, lam=1, x=0.7: quadrature %.10f, formula %.10f' % (num, 2 / lam * float(mp.besselk(0, 2 * a * mp.cosh(lam * x)))))
print('F3 Dyson map from an Ising spin, x_n=S_n/c^n, rho_{n+1}(x)=C exp(beta x^2/2)(rho_n*rho_n)(c x)')
def iterate(c, beta, n):
    p = np.array([0.5, 0.5])
    out = []
    for lev in range(1, n + 1):
        p = np.convolve(p, p)
        S = 2.0 * np.arange(len(p)) - (len(p) - 1)
        x = S / c ** lev
        lw = beta * x * x / 2
        p = p * np.exp(lw - lw.max())
        p = p / p.sum()
        out.append((x * x * p).sum())
    S = 2.0 * np.arange(len(p)) - (len(p) - 1)
    x = S / c ** n
    return p, x, out
print('  c=2 (support stays in [-1,1] by construction)')
for beta in (1.0, 2.0, 4.0, 8.0, 16.0):
    p, x, var = iterate(2.0, beta, 13)
    m2 = (x * x * p).sum(); m4 = (x ** 4 * p).sum()
    print('   beta=%-5g Var x_n at n=5,9,13: %.4g %.4g %.4g; ratio Var_13/Var_12 %.4f; excess at n=13: %.3f; max |x| with p>1e-300: %.3f' % (beta, var[4], var[8], var[12], var[12] / var[11], m4 / m2 ** 2 - 3, np.abs(x[p > 1e-300]).max()))
def tuned(c, n):
    lo, hi = 0.2, 6.0
    for _ in range(46):
        mid = 0.5 * (lo + hi)
        _, _, var = iterate(c, mid, n)
        if var[-1] > var[n // 2 + 1] * 1.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)
for d, n in ((3, 13), (3.5, 13)):
    c = 2 ** ((d + 2) / (2 * d))
    b = tuned(c, n)
    p, x, var = iterate(c, b, n)
    keep = (p > 1e-290) & (x > 0)
    i0 = len(p) // 2
    u = -np.log(p[keep] / p[i0])
    xx = x[keep]
    sel = [np.argmin(np.abs(u - target)) for target in (10, 20, 40, 80, 160, 320, 600)]
    slopes = []
    for s in sel:
        if 1 <= s < len(u) - 1 and u[s] > 0 and u[s + 1] > 0:
            slopes.append(round(float((np.log(u[s + 1]) - np.log(u[s - 1])) / (np.log(xx[s + 1]) - np.log(xx[s - 1]))), 2))
    m2 = (x * x * p).sum(); m4 = (x ** 4 * p).sum()
    print('  d=%g c=%.4f tuned beta=%.8f: Var at n=5,9,13: %.4f %.4f %.4f; excess %.3f; x_max=%.3f; local tail exponent dlog U/dlog x at U=10,20,40,80,160,320,600: %s (formula 2d/(d-2)=%.2f)' % (d, c, b, var[4], var[8], var[12], m4 / m2 ** 2 - 3, x[-1], slopes, 2 * d / (d - 2)))
print('F4 independent sums: X_alpha=sum sigma_j j^(-alpha): local slope of log Lambda(z) against log z, should tend to 1/alpha')
for alpha in (0.75, 5.0 / 6):
    J = 2000000
    jj = np.arange(1, J + 1, dtype=float) ** (-alpha)
    def Lam(z):
        v = z * jj
        return (np.where(v < 30, np.log(np.cosh(np.minimum(v, 30))), v - math.log(2))).sum()
    s = []
    for z in (30.0, 300.0, 3000.0):
        s.append(round((math.log(Lam(z * 1.1)) - math.log(Lam(z))) / math.log(1.1), 4))
    print('  alpha=%.4f slopes at z=30,300,3000: %s; 1/alpha=%.4f; tail exponent p=1/(1-alpha)=%.2f' % (alpha, s, 1 / alpha, 1 / (1 - alpha)))
print('F5 Gaussian damping of the three spins: kernel rho(u)exp(t u^2), t<0. Lee Yang forces 3 k2 k6 >= 10 k4^2')
mp.mp.dps = 60
kernels = {
    'base': lambda u: mp.exp(-2 * mp.pi * mp.cosh(4 * u)),
    'Polya': lambda u: mp.cosh(9 * u) * mp.exp(-2 * mp.pi * mp.cosh(4 * u)),
    'two term': lambda u: (2 * mp.pi ** 2 * mp.cosh(9 * u) - 3 * mp.pi * mp.cosh(5 * u)) * mp.exp(-2 * mp.pi * mp.cosh(4 * u)),
}
h = mp.mpf('0.005')
nodes = [j * h for j in range(0, 301)]
def cumul(vals, t):
    w = [v * mp.exp(t * u * u) * (h if j else h / 2) for j, (u, v) in enumerate(zip(nodes, vals))]
    Z = sum(w)
    m = {p: sum(wi * nodes[j] ** p for j, wi in enumerate(w)) / Z for p in (2, 4, 6)}
    k2 = m[2]; k4 = m[4] - 3 * m[2] ** 2; k6 = m[6] - 15 * m[4] * m[2] + 30 * m[2] ** 3
    return k2, k4, k6
for name, fn in kernels.items():
    vals = [fn(u) for u in nodes]
    def test(t):
        k2, k4, k6 = cumul(vals, t)
        return 3 * k2 * k6 - 10 * k4 ** 2
    def t6(t):
        return cumul(vals, t)[2]
    res = []
    for f in (test, t6):
        hi = mp.mpf(0); lo = mp.mpf(-50)
        while f(lo) > 0 and lo > -20000:
            hi = lo; lo = lo * 2
        for _ in range(60):
            mid = (hi + lo) / 2
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        res.append(float(hi))
    print('  %-9s Lee Yang test fails for t<%.3f (b=|t|/16=%.3f); kappa_6 changes sign at t=%.3f' % (name, res[0], -res[0] / 16, res[1]))
    h2 = mp.mpf('0.01')
    nd = [j * h2 for j in range(0, 141)]
    vv = [fn(u) for u in nd]
    prev = None
    zgrid = [mp.mpf(j) / 4 for j in range(1, 1001)]
    cosm = [[mp.cos(z * u) for u in nd] for z in zgrid[::1]] if False else None
    rows = []
    for t in (-10, -20, -40, -60, -80, -100, -150, -200, -300, -400, -600, -800):
        wts = [v * mp.exp(t * u * u) * (h2 if j else h2 / 2) for j, (u, v) in enumerate(zip(nd, vv))]
        sig = []
        for z in zgrid:
            sig.append(mp.sign(sum(wi * mp.cos(z * u) for wi, u in zip(wts, nd))))
        cnt = sum(1 for a_, b_ in zip(sig[:-1], sig[1:]) if a_ * b_ < 0)
        fz = [float(zgrid[i]) for i in range(len(sig) - 1) if sig[i] * sig[i + 1] < 0]
        rows.append((t, cnt, fz[:2]))
    print('    real zeros with z<250 at t: ' + ' '.join('%d:%d' % (r[0], r[1]) for r in rows))
    print('    first two zeros (coarse) at t=-10,-200,-800: %s %s %s' % (rows[0][2], rows[7][2], rows[11][2]))
