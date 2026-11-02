"""Rerun of D5 and D7 at 1024 bits with radius control. Run: timeout 560 python3 -I N6_d2.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 1024
pi = arb.pi()
def c_(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))
def relrad(z):
    return float((z.real.rad() + z.imag.rad()) / abs(z.mid() if hasattr(z, 'mid') else acb(z.real.mid(), z.imag.mid())))
print('T_5 at t=200+i: 256 bits against 1024 bits')
for prec in (256, 512, 1024):
    ctx.prec = prec
    v = T_ball(acb(200, 1), 5)
    print('  prec %d: %s' % (prec, v.str(8)))
ctx.prec = 1024
def psi_real(t, N):
    tb = acb(t)
    s = acb(0)
    for n in range(1, N + 1):
        s = s + T_ball(tb, n)
    return float((s.real * arb(PI * t / 4).exp()).mid())
def Psi_b(t, N):
    tb = acb(t.real, t.imag)
    return psiN_ball(tb, N) * 4 / (gamma_inf_ball(acb(arb(1) / 2) + I_ * tb) / 2)
for N in (5, 6, 10):
    t1 = np.arange(200.0, 700.0, 0.02)
    v = np.array([psi_real(float(t), N) for t in t1])
    nreal = len(sign_zeros(t1, v)[0])
    pts = rect_path(200, 700, -1, 1, 0.05)
    bad = 0
    vals = []
    for p in pts:
        z = Psi_b(complex(p), N)
        vals.append(c_(z))
    va = np.array(vals)
    d = np.angle(va[1:] / va[:-1])
    print('N=%d real %d, rectangle %.3f (max step %.2f)' % (N, nreal, d.sum() / (2 * PI), np.abs(d).max()), flush=True)
K_ = KAPPA_DH
a5 = [0.0, 1.0, K_, -K_, -1.0]
def TD(t, n):
    a = a5[n % 5]
    if a == 0:
        return acb(0)
    return arb(a) * n * K_ball(2 * pi * n * n / 5, acb(arb(3) / 4) + I_ * t / 2)
def norm_D(t):
    s = acb(arb(1) / 2) + I_ * t
    return (acb(5) / pi).pow((s + 1) / 2) * ((s + 1) / 2).gamma()
for N in (4, 9, 11):
    t1 = np.arange(100.0, 500.0, 0.02)
    v = []
    for t in t1:
        tb = acb(float(t))
        v.append(float((sum((TD(tb, n) for n in range(1, N + 1)), acb(0)).real * arb(PI * float(t) / 4).exp()).mid()))
    nreal = len(sign_zeros(t1, np.array(v))[0])
    pts = rect_path(100, 500, -1, 1, 0.05)
    vals = []
    for p in pts:
        p = complex(p)
        tb = acb(p.real, p.imag)
        s = acb(0)
        for n in range(1, N + 1):
            s = s + TD(tb, n) + TD(-tb, n)
        vals.append(c_(s / norm_D(tb)))
    va = np.array(vals)
    d = np.angle(va[1:] / va[:-1])
    print('DH N=%d real %d, rectangle %.3f (max step %.2f)' % (N, nreal, d.sum() / (2 * PI), np.abs(d).max()), flush=True)
