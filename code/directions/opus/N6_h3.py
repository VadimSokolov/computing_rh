"""Rouche certificate for the limit kernel of the family, with a Lipschitz walk: exact point values in Arb, Lipschitz constant from the contour majorant
|K_nu(z)|<=exp(-|y|th)K_b(z cos th), |d_nu K_nu(z)|<=exp(-|y|th)[th K_b(x)+(K_{b+1}(x)+K_{b-1}(x))/2], x=z cos th, nu=b+iy. Run: timeout 300 python3 -I N6_h3.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 512
pi = arb.pi()
t0 = complex(188.43447686, -0.19249917)
r = 0.02
a, b, c, d = t0.real - r, t0.real + r, t0.imag - r, t0.imag + r
ymin = (t0.real - r) / 2
db = (abs(t0.imag) + r) / 2
def Kr(bb, x):
    return acb(x).bessel_k(arb(bb)).real
def maj_val(n, bb):
    best = None
    z = 2 * math.pi * n * n
    for th in np.linspace(0.02, 1.5707, 400):
        x = arb(z) * arb(float(th)).cos()
        v = (-arb(ymin) * arb(float(th))).exp() * Kr(bb, x)
        if best is None or v < best:
            best = v
    return best
def maj_der(n, bb):
    best = None
    z = 2 * math.pi * n * n
    for th in np.linspace(0.02, 1.5707, 400):
        x = arb(z) * arb(float(th)).cos()
        v = (-arb(ymin) * arb(float(th))).exp() * (arb(float(th)) * Kr(bb, x) + (Kr(bb + 1, x) + Kr(bb - 1, x)) / 2)
        if best is None or v < best:
            best = v
    return best
L = arb(0)
for n in range(1, 5):
    L = L + (2 * pi * pi * n ** 4 * maj_der(n, 9 / 4 + db) + 3 * pi * n * n * maj_der(n, 5 / 4 + db)) / 4
L = arb(float(L.upper()) * 1.001)
print('Lipschitz constant of psi_4 on the box: L<=%s' % L.str(5))
f4 = lambda z: psiN_ball(z, 4)
print('value at the centre: %s; check at N=3 style point' % f4(acb(t0.real, t0.imag)).str(5))
total = arb(0); steps = 0; low = None
cs = [complex(a, c), complex(b, c), complex(b, d), complex(a, d), complex(a, c)]
for u, v_ in zip(cs[:-1], cs[1:]):
    length = abs(v_ - u); e = (v_ - u) / length; pos = 0.0
    v = f4(acb(u.real, u.imag))
    while pos < length:
        step = float((abs(v) / (2 * L)).lower()) * 0.99
        if not step > 0:
            raise RuntimeError('boundary too close to a zero')
        low = abs(v) / 2 if low is None or abs(v) / 2 < low else low
        pos2 = min(length, pos + step)
        z = u + e * pos2
        v2 = f4(acb(z.real, z.imag))
        total = total + (v2 / v).arg()
        v = v2; pos = pos2; steps += 1
print('winding of psi_4 on the boundary: %s in %d steps; |psi_4|>=%s on the boundary' % ((total / (2 * pi)).str(8), steps, low.str(5)))
tot = arb(0)
for n in range(5, 61):
    z = 2 * math.pi * n * n
    th = math.asin(min(1.0, ymin / z))
    x = arb(2) * pi * n * n * arb(th).cos()
    e_ = (-arb(ymin) * arb(th)).exp()
    tot = tot + (2 * pi * pi * n ** 4 * e_ * Kr(9 / 4 + db, x) + 3 * pi * n * n * e_ * Kr(5 / 4 + db, x)) / 2
print('majorant of the terms 5<=n<=60 on the box: %s; ratio to the lower bound: %s' % (tot.str(5), (tot / low).str(5)))
