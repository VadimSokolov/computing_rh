"""Rouche certificate: the limit kernel sum over all n has a nonreal zero; then N=2,3,4 below 200 and N=6. Run: timeout 560 python3 -I N6_h2.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
from N6_common import _seg_ball
ctx.prec = 512
pi = arb.pi()
def c_(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))
def Psi_b(t, N):
    tb = acb(t.real, t.imag)
    return psiN_ball(tb, N) * 4 / (gamma_inf_ball(acb(arb(1) / 2) + I_ * tb) / 2)
t0 = complex(188.43447686, -0.19249917)
r = 0.02
a, b, c, d = t0.real - r, t0.real + r, t0.imag - r, t0.imag + r
f4 = lambda z: psiN_ball(z, 4)
w, ok, ns = cert_winding(f4, a, b, c, d, r / 4, 20)
print('R1 N=4 box centre %.8f%+.8fi half width %.3f: winding %s ok %s, %d pieces' % (t0.real, t0.imag, r, w.str(6) if w is not None else 'NA', ok, ns), flush=True)
# rigorous lower bound of |psi_4| on the boundary
low = None
cs = [complex(a, c), complex(b, c), complex(b, d), complex(a, d), complex(a, c)]
for u, v in zip(cs[:-1], cs[1:]):
    n = 400
    for j in range(n):
        B = f4(_seg_ball(u + (v - u) * j / n, u + (v - u) * (j + 1) / n))
        val = abs(B).lower()
        low = val if low is None else (val if val < low else low)
print('   lower bound of |psi_4| on the boundary: %s' % arb(low).str(5))
# majorant |K_nu(z)|<=exp(-y theta)K_b(z cos theta), nu=b+iy or b-iy, b=Re nu, 0<=theta<pi/2
ymin = (t0.real - r) / 2
bmax1 = 9 / 4 + (abs(t0.imag) + r) / 2
bmax2 = 5 / 4 + (abs(t0.imag) + r) / 2
tot = arb(0)
terms = []
for n in range(5, 61):
    z = 2 * math.pi * n * n
    th = math.asin(min(1.0, ymin / z))
    zz = arb(2) * pi * n * n * arb(th).cos()
    e = (-arb(ymin) * arb(th)).exp()
    bd = 2 * pi * pi * n ** 4 * e * acb(zz).bessel_k(bmax1).real + 3 * pi * n * n * e * acb(zz).bessel_k(bmax2).real
    # psi=(F(t)+F(-t))/4 so each T_n enters as a quarter of two terms of this size
    tot = tot + bd / 2
    terms.append(bd / 2)
print('   majorant of |sum_{5<=n<=60} T_n contribution to psi|: %s; first term %s; ratio to lower bound of |psi_4|: %s' % (tot.str(5), terms[0].str(5), (tot / arb(low)).str(5)), flush=True)
print('H4 N=2,3,4 below 200 (N=4 to 190): rectangle (ball values) against real zeros')
for N, top in ((2, 200.0), (3, 180.0), (4, 188.0)):
    tt = np.arange(0.5, top, 0.01)
    v = np.array([float((sum((T_ball(acb(float(x)), n) for n in range(2, N + 1)), T_ball(acb(float(x)), 1)).real * arb(PI * float(x) / 4).exp()).mid()) for x in tt])
    nreal = len(sign_zeros(tt, v)[0])
    pts = rect_path(0.5, top, -0.75, 0.75, 0.02)
    va = np.array([c_(Psi_b(complex(p), N)) for p in pts])
    dd = np.angle(va[1:] / va[:-1])
    print('  N=%d to %g: real %d rectangle %.3f (max step %.2f)' % (N, top, nreal, dd.sum() / (2 * PI), np.abs(dd).max()), flush=True)
print('H5 N=6 on 200<t<700')
N = 6
t1 = np.arange(200.0, 700.0, 0.02)
v = np.array([float((sum((T_ball(acb(float(x)), n) for n in range(2, N + 1)), T_ball(acb(float(x)), 1)).real * arb(PI * float(x) / 4).exp()).mid()) for x in t1])
nreal = len(sign_zeros(t1, v)[0])
pts = rect_path(200, 700, -1, 1, 0.1)
va = np.array([c_(Psi_b(complex(p), N)) for p in pts])
dd = np.angle(va[1:] / va[:-1])
print('  N=6 real %d, rectangle %.3f (max step %.2f)' % (nreal, dd.sum() / (2 * PI), np.abs(dd).max()), flush=True)
