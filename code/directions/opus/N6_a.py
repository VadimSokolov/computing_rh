"""Checks of Sections 1 to 3 and of the Gram theorems of cosh_spin.tex (v7). Run: timeout 300 python3 -I N6_a.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
import mpmath as mp
mp.mp.dps = 40
ctx.prec = 256
print('A1 Proposition bessel: direct quadrature against the Bessel form')
def kern(u):
    return (2 * mp.pi ** 2 * mp.cosh(9 * u) - 3 * mp.pi * mp.cosh(5 * u)) * mp.exp(-2 * mp.pi * mp.cosh(4 * u))
def psi_mp(x):
    return 0.5 * mp.re(2 * mp.pi ** 2 * mp.besselk(mp.mpf(9) / 4 + 1j * x / 4, 2 * mp.pi) - 3 * mp.pi * mp.besselk(mp.mpf(5) / 4 + 1j * x / 4, 2 * mp.pi))
for x in (0, 10, 37.5):
    d = 2 * mp.quad(lambda u: kern(u) * mp.cos(x * u), mp.linspace(0, 2.5, 11))
    print(' x=%g direct %s form %s reldiff %s' % (x, mp.nstr(d, 25), mp.nstr(psi_mp(x), 25), mp.nstr((d - psi_mp(x)) / d, 3)))
a = 2 * mp.pi
for x in (0, 10):
    d = 2 * mp.quad(lambda u: mp.exp(-a * mp.cosh(4 * u)) * mp.cos(x * u), mp.linspace(0, 2.5, 11))
    f = 0.5 * mp.re(mp.besselk(1j * mp.mpf(x) / 4, a))
    print(' base x=%g direct %s form %s' % (x, mp.nstr(d, 25), mp.nstr(f, 25)))

print('A2 Proposition interval at a=2pi: P(A<B)=a K_1(a), E(B-A)_+=2K_0(a)')
pAB = mp.quad(lambda x: mp.exp(-x - a * a / (4 * x)), [0, 1, 5, 20, 100, mp.inf])
print(' P(A<B) quad %s, a K1(a) %s; 2K0(a) %s' % (mp.nstr(pAB, 15), mp.nstr(a * mp.besselk(1, a), 15), mp.nstr(2 * mp.besselk(0, a), 15)))

print('A3 U_*\'\'\' on (0,2]')
def Uppp(u):
    u = mp.mpf(u)
    p0 = 2 * mp.pi ** 2 * mp.cosh(9 * u) - 3 * mp.pi * mp.cosh(5 * u)
    p1 = 18 * mp.pi ** 2 * mp.sinh(9 * u) - 15 * mp.pi * mp.sinh(5 * u)
    p2 = 162 * mp.pi ** 2 * mp.cosh(9 * u) - 75 * mp.pi * mp.cosh(5 * u)
    p3 = 1458 * mp.pi ** 2 * mp.sinh(9 * u) - 375 * mp.pi * mp.sinh(5 * u)
    return 128 * mp.pi * mp.sinh(4 * u) - (p3 / p0 - 3 * p2 * p1 / p0 ** 2 + 2 * p1 ** 3 / p0 ** 3)
for u in (1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1, 2):
    print('  u=%-7g U3=%s' % (u, mp.nstr(Uppp(u), 8)))
print('  U4(0)=', mp.nstr(mp.diff(Uppp, mp.mpf(0)), 8))
gr = [i / 2000 for i in range(1, 4001)]
vals = [float(Uppp(u)) for u in gr]
imin = int(np.argmin(vals[100:])) + 100
print('  grid 0.0005 on (0,2]: first value %.4f, global grid min %.6f at u=%.4f; local min beyond u=0.05: %.4f at u=%.4f' % (vals[0], min(vals), gr[int(np.argmin(vals))], vals[imin], gr[imin]))
lm = [(gr[i], vals[i]) for i in range(1, 3999) if vals[i] < vals[i - 1] and vals[i] < vals[i + 1]]
lx = [(gr[i], vals[i]) for i in range(1, 3999) if vals[i] > vals[i - 1] and vals[i] > vals[i + 1]]
print('  local minima', lm, 'local maxima', lx)

print('A4 cumulants of rho_*')
mom = {}
Z = mp.quad(kern, mp.linspace(0, 2.5, 11))
for j in range(1, 6):
    mom[2 * j] = mp.quad(lambda u: u ** (2 * j) * kern(u), mp.linspace(0, 2.5, 11)) / Z
m = mom
k2 = m[2]; k4 = m[4] - 3 * m[2] ** 2; k6 = m[6] - 15 * m[4] * m[2] + 30 * m[2] ** 3
k8 = m[8] - 28 * m[6] * m[2] - 35 * m[4] ** 2 + 420 * m[4] * m[2] ** 2 - 630 * m[2] ** 4
k10 = m[10] - 45 * m[8] * m[2] - 210 * m[6] * m[4] + 1260 * m[6] * m[2] ** 2 + 3150 * m[4] ** 2 * m[2] - 18900 * m[4] * m[2] ** 3 + 22680 * m[2] ** 5
print('  kappa2..10', [mp.nstr(x, 5) for x in (k2, k4, k6, k8, k10)], 'k4/k2^2', mp.nstr(k4 / k2 ** 2, 5))
def G(z):
    z = mp.mpf(z)
    A = mp.mpf(9) / 4; B = mp.mpf(5) / 4
    return (2 * mp.pi ** 2 * (mp.besselk(A + z / 4, a) + mp.besselk(A - z / 4, a)) - 3 * mp.pi * (mp.besselk(B + z / 4, a) + mp.besselk(B - z / 4, a))) / 4
def K3(z):
    h = mp.mpf('1e-3')
    lg = [mp.log(G(z + j * h)) for j in (-2, -1, 1, 2)]
    return (lg[3] - 2 * lg[2] + 2 * lg[1] - lg[0]) / (2 * h ** 3)
vs = [(z, K3(z)) for z in [0.5 * i for i in range(1, 241)] + [200]]
print('  K3<0 on grid of 240 points in (0,120] and z=200:', all(v < 0 for _, v in vs), 'max over (0,120]', mp.nstr(max(v for z, v in vs if z <= 120), 4), 'at z=', [z for z, v in vs if v == max(w for y, w in vs if y <= 120)], 'K3(0.5)', mp.nstr(vs[0][1], 4), 'K3(120)', mp.nstr(vs[239][1], 4), 'K3(200)', mp.nstr(vs[240][1], 4))

print('A5 zeros of psi(2t), N=1')
ts = np.arange(0.2, 705.0, 0.01)
gv = g_real(ts, 1)
z0, idx = sign_zeros(ts, gv)
zz = np.array([refine_zero(lambda x: float(g_real(np.array([x]), 1)[0]), ts[i], ts[i + 1]) for i in idx])
print('  first 10 zeros', np.round(zz[:10], 3))
th0 = np.array([float(arb(0))])
# check against ball arithmetic at sample points
for t in (7.6, 18.0, 100.3, 699.7):
    v = psi_ball(acb(t), 1)
    A = acb(0)
    d = float(g_real(np.array([t]), 1)[0])
    s = acb(0.5, t)
    nor = gamma_inf_ball(s)
    ratio = (F_ball(acb(t), 1)) / (nor / 2)
    print('  check A_1 at t=%g: double %s ball %s' % (t, A_N(np.array([t + 0j]), 1)[0], ratio.str(12)))
for T in (100, 260):
    print('  N(%d)=%d main %.3f' % (T, np.sum(zz < T), main_term(T)))
mid = np.sum((zz > 200) & (zz < 700))
print('  N(200,700)=%d main diff %.3f' % (mid, main_term(700) - main_term(200)))
w1, jump = winding(lambda t: Psi_N(t, 1), rect_path(200, 700, -1, 1, 0.01))
print('  argument principle, t in [200,700], |Im t|<=1: %.6f (max step %.3f)' % (w1, jump))
print('Gram')
th = theta(ts)
kmax = int(th[-1] / PI)
mask = ts > 6.5
gk = {}
for k in range(-1, kmax + 1):
    x = np.interp(k * PI, th[mask], ts[mask])
    for _ in range(4):
        x = x - (theta(np.array([x]))[0] - k * PI) / (0.5 * (np.real(__import__('scipy.special', fromlist=['digamma']).digamma(0.25 + 0.5j * x))) - 0.5 * math.log(PI))
    gk[k] = x
print('  g_-1=%.5f g_0=%.5f g_1=%.5f g_2=%.5f g_3=%.5f g_4=%.5f g_5=%.5f' % tuple(gk[k] for k in (-1, 0, 1, 2, 3, 4, 5)))
print('  theta min on (0,10): %.4f at t=%.3f' % (th[ts < 12].min(), ts[ts < 12][np.argmin(th[ts < 12])]))
print('  zeros below g_-1: %d; in [g_-1,g_0): %d' % (np.sum(zz < gk[-1]), np.sum((zz >= gk[-1]) & (zz < gk[0]))))
ks = [k for k in range(0, kmax) if gk[k + 1] <= 700]
cnt = [int(np.sum((zz >= gk[k]) & (zz < gk[k + 1]))) for k in ks]
print('  Gram intervals [g_k,g_k+1) with g_k+1<=700: %d, counts min %d max %d' % (len(ks), min(cnt), max(cnt)))
f = np.mod(theta(zz) / PI, 1.0)
for t0 in (36.6, 258):
    i = np.argmin(abs(zz - t0))
    print('  zero %.3f frac(theta/pi)=%.4f (1/2-frac)*t=%.3f' % (zz[i], f[i], (0.5 - f[i]) * zz[i]))
sel = zz > 200
print('  mean (1/2-frac)*t for t>200: %.4f; 2pi=%.4f' % (np.mean((0.5 - f[sel]) * zz[sel]), 2 * PI))
EE = A_N(ts + 0j, 1) - 1
absE = np.abs(EE)
print('  largest t with |E|>=1: %.3f' % ts[absE >= 1].max())
for t0 in (17.846, 18.5, 20, 25, 30, 31, 35, 40, 700):
    i = np.argmin(abs(ts - t0)); y = ts[i] / 2
    eps = math.exp(PI ** 2 / y) - 1
    pref = (abs(1.25 + 1j * y) + 1.5) / abs(-0.25 + 1j * y)
    print('  t=%-7g |E|=%.5f t|E|=%.3f note bound e^(2pi^2/t)-1=%.4f corrected bound %.4f delta=%.4f  delta/(2pi^2/t)=%.4f' % (t0, absE[i], t0 * absE[i], eps, pref * eps, np.angle(1 + EE[i]), np.angle(1 + EE[i]) * t0 / (2 * PI ** 2)))
print('  phase at g_0: cos(theta+delta) sign, |1+E| at g_0: %s' % abs(A_N(np.array([gk[0] + 0j]), 1)[0]))
print('  (-1)^k cos(theta+delta) at g_k, k=0..6:', [round(float((-1) ** k * np.real(np.exp(1j * (k * PI)) * A_N(np.array([gk[k] + 0j]), 1)[0]) / abs(A_N(np.array([gk[k] + 0j]), 1)[0])), 4) for k in range(0, 7)])

print('A6 Polya spin and base spin on Gram intervals')
def scan(F, tmax=705, h=0.02, tmin=0.2):
    tt = np.arange(tmin, tmax, h)
    out = np.empty(len(tt))
    for i, t in enumerate(tt):
        out[i] = float(F(float(t)).real.mid()) * math.exp(PI * t / 4)
    return tt, out
pi_b = arb.pi()
Fp = lambda t: 2 * pi_b * pi_b * K_ball(2 * pi_b, acb(arb(9) / 4, t / 2))
Fb = lambda t: K_ball(2 * pi_b, acb(0, t / 2))
for name, F in (('Polya', Fp), ('base', Fb)):
    tt, vv = scan(F)
    zs, _ = sign_zeros(tt, vv)
    fr = np.mod(theta(zs) / PI, 1.0)
    cnt = [int(np.sum((zs >= gk[k]) & (zs < gk[k + 1]))) for k in ks]
    print('  %s: first zeros %s; below g_-1 %d; in [g_-1,g_0) %d; Gram counts min %d max %d over %d intervals; N(100)=%d N(260)=%d N(200,700)=%d; mean frac(theta/pi) t>500: %.4f' % (name, np.round(zs[:8], 3), np.sum(zs < gk[-1]), np.sum((zs >= gk[-1]) & (zs < gk[0])), min(cnt), max(cnt), len(ks), np.sum(zs < 100), np.sum(zs < 260), np.sum((zs > 200) & (zs < 700)), np.mean(fr[zs > 500])))

print('A7 Liouville wall at a=2pi: zeros nu_k of K_{i nu}(2pi)')
nus = np.arange(0.5, 425.0, 0.01)
vv = np.array([float(K_ball(2 * pi_b, acb(0, float(n))).real.mid()) * math.exp(PI * n / 2) for n in nus])
_, ii = sign_zeros(nus, vv)
nz = np.array([refine_zero(lambda x: float(K_ball(2 * pi_b, acb(0, x)).real.mid()) * math.exp(PI * x / 2), nus[i], nus[i + 1], 1e-10) for i in ii])
print('  count up to 425: %d; nu_1,2,3=%s; nu_490=%.4f' % (len(nz), np.round(nz[:3], 4), nz[489]))
phi0 = nz * np.log(nz / PI) - nz + PI / 4
off = phi0 / PI - np.arange(1, len(nz) + 1)
print('  WKB rule offsets at k=1,2,10,100,490: %s; book prediction -(pi^2-1/12)/(pi nu): %s' % (np.round(off[[0, 1, 9, 99, 489]], 4), np.round(-(PI ** 2 - 1 / 12) / (PI * nz[[0, 1, 9, 99, 489]]), 4)))
for u in (1, 10, 50):
    exact = float((K_ball(2 * pi_b, acb(0)) / K_ball(2 * pi_b, acb(math.sqrt(u)))).real.mid())
    pr = -np.sum(np.log1p(u / nz[:490] ** 2))
    V = nz[489]
    x = np.linspace(V, V * 4000, 2000001)
    tail = np.trapezoid(np.log1p(u / x ** 2) * np.log(x / PI) / PI, x) + u / (x[-1] * PI) * (np.log(x[-1] / PI) + 1) - 0.5 * math.log1p(u / V ** 2)
    print('  lambda=%d exact %.7f; product with 490 zeros alone %.7f; with Weyl tail %.7f' % (u, exact, math.exp(pr), math.exp(pr - tail)))

print('A8 DH control: cosh spin of the Davenport Heilbronn function, Re K_{3/4+it/2}(2pi/5)')
aDH = 2 * pi_b / 5
tt = np.arange(0.2, 505.0, 0.02)
vv = np.array([float(K_ball(aDH, acb(arb(3) / 4, float(t) / 2)).real.mid()) * math.exp(PI * t / 4) for t in tt])
zd, _ = sign_zeros(tt, vv)
thd = np.imag(loggamma(0.75 + 0.5j * tt)) + 0.5 * tt * math.log(5 / PI)
ktop = int(thd[-1] / PI)
md = tt > 3
gd = {k: np.interp(k * PI, thd[md], tt[md]) for k in range(int(np.ceil(thd[md][0] / PI)), ktop)}
keys = sorted(gd)
cnts = [int(np.sum((zd >= gd[k]) & (zd < gd[k + 1]))) for k in keys[:-1]]
frd = np.mod(np.interp(zd, tt, thd) / PI, 1.0)
print('  first zeros %s; DH Gram points from k=%d (t=%.3f); counts per interval min %d max %d over %d intervals; first exceptional k: %s; mean frac(theta_DH/pi) t>300: %.4f' % (np.round(zd[:8], 3), keys[0], gd[keys[0]], min(cnts), max(cnts), len(cnts), [keys[i] for i, c in enumerate(cnts) if c != 1][:5], np.mean(frd[zd > 300])))
