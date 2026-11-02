"""Last numbers: derivative of the reflected term at t=34, DH folded with exact kappa, persistence of the nonreal zero as N grows, N=2 below 200, N=6. Run: timeout 560 python3 -I N6_h.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 512
pi = arb.pi(); PI2 = pi * pi
print('H1 reflected term on the circle |t-34|=1/4')
def refl_ball(t):
    t = acb(t); y = I_ * t / 2
    n1 = acb(arb(9) / 4) + y; n2 = acb(arb(5) / 4) + y
    z = acb(2 * pi)
    r = 2 * PI2 * (-pi / (2 * (pi * n1).sin())) * z.bessel_i(n1) - 3 * pi * (-pi / (2 * (pi * n2).sin())) * z.bessel_i(n2)
    return r / (gamma_inf_ball(acb(arb(1) / 2) + I_ * t) / 2)
sup = 0.0
M = 160
for j in range(M):
    ang = 2 * math.pi * (j + 0.5) / M
    c = complex(34 + 0.25 * math.cos(ang), 0.25 * math.sin(ang))
    b = acb(arb(c.real, 0.25 * math.pi / M * 1.01), arb(c.imag, 0.25 * math.pi / M * 1.01))
    sup = max(sup, float(abs(refl_ball(b)).upper()))
print('  sup|R|<=%.3e so |R\'(34)|<=4 sup=%.3e; with |E_m\'|<=0.0352, |E|<=0.858: |delta\'(34)|<=%.4f; theta\'(34)=0.8442' % (sup, 4 * sup, (0.0352 + 4 * sup) / (1 - 0.858)))
print('H2 DH folded truncations with exact kappa')
s5 = arb(5).sqrt()
kap = ((10 - 2 * s5).sqrt() - 2) / (s5 - 1)
print('  kappa=%s' % kap.str(20))
A5 = [arb(0), arb(1), kap, -kap, arb(-1)]
def XD(t, N):
    s = acb(0); t = acb(t); w = I_ * t / 2
    for n in range(1, N + 1):
        a = A5[n % 5]
        if n % 5 == 0:
            continue
        z = pi * n * n / 5
        s = s + a * n * acb(z).pow(-(acb(arb(3) / 4) + w)) * acb(z).gamma_upper(acb(arb(3) / 4) + w)
    return s.real * arb(PI * float(t.real) / 4).exp()
print('  sum of all first derivatives: %s' % sum((A5[n % 5] * n * (-pi * n * n / 5).exp() * (3 - 4 * pi * n * n / 5) for n in range(1, 60)), arb(0)).str(5))
for N in (8, 9, 11):
    M_ = N + 1 if (N + 1) % 5 else N + 2
    pred = 4 * M_ * M_ / 5
    kD = sum((A5[n % 5] * n * (-pi * n * n / 5).exp() * (3 - 4 * pi * n * n / 5) for n in range(1, N + 1)), arb(0))
    tt = np.arange(0.25, pred + 40, 0.02)
    v = np.array([float(XD(float(t), N).mid()) for t in tt])
    z, _ = sign_zeros(tt, v)
    print('  N=%d real zeros %d, largest %.3f, prediction 4M^2/5=%.1f (M=%d), kink %.4g' % (N, len(z), z[-1], pred, M_, float(kD)))
print('H3 the nonreal zero of N=3 as N grows')
def Psi_b(t, N):
    tb = acb(t.real, t.imag)
    return psiN_ball(tb, N) * 4 / (gamma_inf_ball(acb(arb(1) / 2) + I_ * tb) / 2)
def c_(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))
t = complex(188.433681, -0.274981)
for N in (3, 4, 5, 6, 8):
    for _ in range(25):
        h = 1e-7
        f0 = c_(Psi_b(t, N)); f1 = c_(Psi_b(t + h, N))
        dt = f0 / ((f1 - f0) / h)
        t = t - dt
        if abs(dt) < 1e-11:
            break
    print('  N=%d zero at %.8f%+.8fi' % (N, t.real, t.imag), flush=True)
ctx.prec = 768
r = 0.02
w, ok, ns = cert_winding(lambda z: psiN_ball(z, 8), t.real - r, t.real + r, t.imag - r, t.imag + r, r / 4, 20)
mn = min(float(abs(psiN_ball(acb(t.real + r * math.cos(a), t.imag + r * math.sin(a)), 8)).lower()) * math.exp(PI * 188 / 4) for a in np.linspace(0, 6.283, 41))
print('  N=8 certified box half width %.3f: winding %s ok %s; on a circle inside |psi_8|e^(pi 188/4)>=%.3e' % (r, w.str(6), ok, mn))
tail = sum(2 * math.pi ** 2 * n ** 4 * math.exp(-2 * math.pi * n * n + 2.75 * math.log(2)) for n in range(9, 12))
print('  tail n>=9: bound by 2 pi^2 n^4 K_{Re nu}(2 pi n^2) with K_b(z)<=e^{-z}e^{b^2/z}: log10 about %.0f, against psi at 10^(%.0f)' % (math.log10(tail), -PI * 188 / 4 / math.log(10)))
ctx.prec = 512
print('H4 N=2 and N=3 below 200: rectangle count (ball values) against real zeros')
for N in (2, 3):
    tt = np.arange(0.5, 200.0, 0.01)
    v = np.array([float((sum((T_ball(acb(float(x)), n) for n in range(2, N + 1)), T_ball(acb(float(x)), 1)).real * arb(PI * float(x) / 4).exp()).mid()) for x in tt])
    nreal = len(sign_zeros(tt, v)[0])
    pts = rect_path(0.5, 200, -0.75, 0.75, 0.02)
    va = np.array([c_(Psi_b(complex(p), N)) for p in pts])
    d = np.angle(va[1:] / va[:-1])
    print('  N=%d real %d rectangle %.3f (max step %.2f)' % (N, nreal, d.sum() / (2 * PI), np.abs(d).max()), flush=True)
print('H5 N=6 on 200<t<700')
N = 6
t1 = np.arange(200.0, 700.0, 0.02)
v = np.array([float((sum((T_ball(acb(float(x)), n) for n in range(2, N + 1)), T_ball(acb(float(x)), 1)).real * arb(PI * float(x) / 4).exp()).mid()) for x in t1])
nreal = len(sign_zeros(t1, v)[0])
pts = rect_path(200, 700, -1, 1, 0.1)
va = np.array([c_(Psi_b(complex(p), N)) for p in pts])
d = np.angle(va[1:] / va[:-1])
print('  N=6 real %d, rectangle %.3f (max step %.2f)' % (nreal, d.sum() / (2 * PI), np.abs(d).max()), flush=True)
