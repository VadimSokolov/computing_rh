"""Family rho_{*,N}, second pass in ball arithmetic: lowest nonreal zero of N=3, zeros below 100, ratios, larger N, polygon test, DH control family. Run: timeout 560 python3 -I N6_d.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 256
pi = arb.pi()
def Fa(t, N):
    """F_N(t)/((1/2)gamma_inf(1/2+it)) as a complex float, t python complex."""
    tb = acb(t.real, t.imag)
    v = sum((T_ball(tb, n) for n in range(2, N + 1)), T_ball(tb, 1)) / (gamma_inf_ball(acb(arb(1) / 2) + I_ * tb) / 2)
    return v
def Psi_ball(t, N):
    tb = acb(t.real, t.imag)
    return psiN_ball(tb, N) * 4 / (gamma_inf_ball(acb(arb(1) / 2) + I_ * tb) / 2)
def c_(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))
def gr(t, N):
    """real function with the sign of psi_N(2t) up to a fixed sign, scaled, real t."""
    return float(psi_N_real(t, N))
def psi_N_real(t, N):
    tb = acb(t)
    s = T_ball(tb, 1)
    for n in range(2, N + 1):
        s = s + T_ball(tb, n)
    return (s.real * arb(PI * t / 4).exp()).mid()
print('D1 ratio |T_2|/|T_1| at single heights in ball arithmetic')
print('  ' + ' '.join('t=%g:%.3g' % (t, float(abs(T_ball(acb(t), 2)) / abs(T_ball(acb(t), 1)))) for t in (20, 30, 40, 42, 44, 46, 50, 55, 60, 70, 80, 100, 120, 200, 400, 2000)))
print('D2 zeros below 100, ball arithmetic values')
tt = np.arange(0.2, 100.0, 0.01)
zz = {}
for N in (1, 2, 3):
    v = np.array([float(psi_N_real(float(t), N)) for t in tt])
    z, _ = sign_zeros(tt, v)
    zz[N] = z
    print('  N=%d: %d zeros; first 5 %s; last 3 %s' % (N, len(z), np.round(z[:5], 3), np.round(z[-3:], 3)))
for N in (2, 3):
    if len(zz[N]) == len(zz[1]):
        d = np.abs(zz[N] - zz[1])
        k = [np.abs(d[(zz[1] > a) & (zz[1] <= b)]).max() for a, b in ((0, 30), (30, 40), (40, 50), (50, 70), (70, 100))]
        print('  N=%d against N=1: max|diff| on (0,30] (30,40] (40,50] (50,70] (70,100]: %s' % (N, ['%.2g' % x for x in k]))
print('D3 small blocks below 300 for N=3: real zeros against rectangle count')
out = []
hits = []
for lo in range(0, 300, 10):
    hi = lo + 10
    if lo == 0:
        lo_ = 0.5
    else:
        lo_ = lo
    t1 = np.arange(lo_, hi, 0.01)
    v = np.array([float(psi_N_real(float(t), 3)) for t in t1])
    nreal = len(sign_zeros(t1, v)[0])
    pts = rect_path(lo_, hi, -0.75, 0.75, 0.04)
    vals = [c_(Psi_ball(complex(p), 3)) for p in pts]
    va = np.array(vals)
    w = np.angle(va[1:] / va[:-1]).sum() / (2 * PI)
    if abs(w - nreal) > 0.5:
        hits.append((lo, nreal, round(w, 2)))
print('  blocks with a deficit (lo, real, rectangle): %s' % hits)
print('D4 nonreal zeros of rho_{*,3}: Newton from a grid')
sols = []
for lo, nr, w in hits[:6]:
    for x0 in np.arange(lo + 0.5, lo + 10, 1.0):
        for y0 in (-0.05, -0.15, -0.3, -0.5):
            t = complex(x0, y0)
            ok = False
            for _ in range(30):
                h = 1e-6
                f0 = c_(Psi_ball(t, 3))
                f1 = c_(Psi_ball(t + h, 3))
                dt = f0 / ((f1 - f0) / h)
                if abs(dt) > 2:
                    break
                t = t - dt
                if abs(dt) < 1e-11:
                    ok = True
                    break
            if ok and t.imag < -1e-5 and not any(abs(t - q) < 1e-6 for q in sols):
                sols.append(t)
sols = sorted(sols, key=lambda z: z.real)
print('  lowest: ' + ', '.join('%.6f%+.6fi' % (z.real, z.imag) for z in sols[:8]))
if sols:
    tm = sols[0]
    r = 0.25 * min(abs(tm.imag), 0.1)
    w, ok, ns = cert_winding(lambda t: psiN_ball(t, 3), tm.real - r, tm.real + r, tm.imag - r, tm.imag + r, r / 4, 22)
    print('  certified box centre %.6f%+.6fi half width %.4f: winding %s, ok %s, %d pieces' % (tm.real, tm.imag, r, w.str(6) if w is not None else 'NA', ok, ns))
print('D5 larger N on 200<t<700: real sign changes against rectangle (ball values)')
for N in (4, 5, 6, 10):
    t1 = np.arange(200.0, 700.0, 0.02)
    v = np.array([float(psi_N_real(float(t), N)) for t in t1])
    nreal = len(sign_zeros(t1, v)[0])
    pts = rect_path(200, 700, -1, 1, 0.1)
    va = np.array([c_(Psi_ball(complex(p), N)) for p in pts])
    d = np.angle(va[1:] / va[:-1])
    print('  N=%d real %d, rectangle %.3f (max step %.2f)' % (N, nreal, d.sum() / (2 * PI), np.abs(d).max()))
print('D6 polygon test: is sum_{n<=N} chi(n) n^(-sigma)=0 solvable at sigma=0.52, chi completely multiplicative unimodular?')
def primes(n):
    return [p for p in range(2, n + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
rng = np.random.default_rng(1)
fails = []
for N in range(3, 201):
    sigma = 0.52
    P = primes(N)
    p = P[-1]
    Q = P[:-1]
    ns = [n for n in range(1, N + 1) if n != p]
    E = np.zeros((len(ns), len(Q)))
    for i, n in enumerate(ns):
        m = n
        for j, q in enumerate(Q):
            while m % q == 0:
                m //= q
                E[i, j] += 1
    co = np.array(ns, dtype=float) ** (-sigma)
    target = p ** (-sigma)
    best = 1e9
    for trial in range(30):
        ph = np.pi * np.ones(len(Q)) if trial == 0 else rng.uniform(0, 2 * np.pi, len(Q))
        val = abs((co * np.exp(1j * E @ ph)).sum()) if len(Q) else abs(co.sum())
        for sweep in range(25):
            for j in range(len(Q)):
                grid = np.linspace(0, 2 * np.pi, 49)[:-1]
                base = E @ ph - E[:, j] * ph[j]
                vals = [abs((co * np.exp(1j * (base + E[:, j] * g))).sum()) for g in grid]
                k = int(np.argmin(vals))
                if vals[k] < val - 1e-12:
                    val = vals[k]
                    ph[j] = grid[k]
            if val <= target:
                break
        best = min(best, val)
        if best <= target:
            break
    if best > target:
        fails.append((N, round(best, 3), round(target, 3)))
print('  N in 3..200 where the search did not close the polygon: %s' % fails)
print('D7 Davenport Heilbronn control family, sum a_n n cosh(3u) exp(-(2 pi n^2/5) cosh 4u)')
K_ = KAPPA_DH
a5 = [0.0, 1.0, K_, -K_, -1.0]
def TD(t, n):
    a = a5[n % 5]
    if a == 0:
        return acb(0)
    return arb(a) * n * K_ball(2 * pi * n * n / 5, acb(arb(3) / 4) + I_ * t / 2)
def psiD_ball(t, N):
    t = acb(t)
    s = acb(0)
    for n in range(1, N + 1):
        s = s + TD(t, n) + TD(-t, n)
    return s
def norm_D(t):
    s = acb(arb(1) / 2) + I_ * t
    return (acb(5) / pi).pow((s + 1) / 2) * ((s + 1) / 2).gamma()
print('  sum of |a_n| n^(-1/2), 2<=n<=N: ' + ' '.join('N=%d:%.3f' % (N, sum(abs(a5[n % 5]) * n ** -0.5 for n in range(2, N + 1))) for n_ in [0] for N in (4, 5, 6, 7, 9, 11)))
for N in (4, 6, 7, 9, 11):
    t1 = np.arange(100.0, 500.0, 0.02)
    v = np.array([float((sum((TD(acb(float(t)), n) for n in range(1, N + 1)), acb(0)).real * arb(PI * t / 4).exp()).mid()) for t in t1])
    nreal = len(sign_zeros(t1, v)[0])
    pts = rect_path(100, 500, -1, 1, 0.1)
    va = []
    for p in pts:
        p = complex(p)
        tb = acb(p.real, p.imag)
        va.append(c_(psiD_ball(tb, N) / norm_D(tb)))
    va = np.array(va)
    d = np.angle(va[1:] / va[:-1])
    print('  N=%d real %d, rectangle %.3f (max step %.2f)' % (N, nreal, d.sum() / (2 * PI), np.abs(d).max()))
print('D8 N=2 further up (double precision) and minimal phase derivative for N=1,2,3')
for lo, hi in ((2100, 3100), (3100, 4100), (4100, 6100)):
    t1 = np.arange(lo, hi, 0.005)
    v = np.real(np.exp(1j * theta(t1)) * A_N(t1 + 0j, 2))
    z, _ = sign_zeros(t1, v)
    w, j = winding(lambda t: Psi_N(t, 2), rect_path(lo, hi, -1, 1, 0.02))
    print('  N=2 [%d,%d]: real %d, rectangle %.3f' % (lo, hi, len(z), w))
t1 = np.arange(60.0, 6000.0, 0.01)
for N in (1, 2, 3):
    ph = theta(t1) + np.unwrap(np.angle(A_N(t1 + 0j, N)))
    dp = np.diff(ph) / 0.01
    print('  N=%d min phase derivative on [60,6000]: %.4f at t=%.2f; fraction of t with negative derivative %.4f' % (N, dp.min(), t1[int(np.argmin(dp))], np.mean(dp < 0)))
