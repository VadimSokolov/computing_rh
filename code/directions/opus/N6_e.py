"""Folded truncations (Theorem cusp), DH folded control, Robin family of the Liouville wall. Run: timeout 500 python3 -I N6_e.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 768
pi = arb.pi()
print('E0 calling convention, Gamma(2,1)=2/e=0.7357589:', acb(1).gamma_upper(2).str(8), acb(2).gamma_upper(1).str(8))
def Gu(s, z):
    return acb(z).gamma_upper(s)
def X_term(t, n):
    t = acb(t)
    z = pi * n * n
    w = I_ * t / 2
    return 2 * pi * pi * n ** 4 * acb(z).pow(-(acb(arb(9) / 4) + w)) * Gu(acb(arb(9) / 4) + w, z) - 3 * pi * n * n * acb(z).pow(-(acb(arb(5) / 4) + w)) * Gu(acb(arb(5) / 4) + w, z)
def Xi_real(t, N):
    s = acb(0)
    for n in range(1, N + 1):
        s = s + X_term(t, n)
    return s.real * arb(PI * t / 4).exp()
def kink(N):
    s = arb(0)
    for n in range(1, N + 1):
        s = s + pi * n * n * (-pi * n * n).exp() * (30 * pi * n * n - 15 - 8 * pi * pi * n ** 4)
    return s
print('E1 phi_1\'(0)=%.5f phi_2\'(0)=%.5f' % (float(kink(1)), float(kink(2) - kink(1))))
print('E2 test: N=6 first zeros (should be the zeta zeros)')
tt = np.arange(0.5, 45.0, 0.02)
v = np.array([float(Xi_real(float(t), 6).mid()) for t in tt])
print('   ', np.round(sign_zeros(tt, v)[0][:9], 4))
for N in (1, 2, 3, 4):
    Tm = {1: 40, 2: 70, 3: 95, 4: 135}[N]
    tt = np.arange(0.25, Tm, 0.01)
    v = np.array([float(Xi_real(float(t), N).mid()) for t in tt])
    z, _ = sign_zeros(tt, v)
    lz = z[-1] if len(z) else float('nan')
    far = [(t, float(Xi_real(float(t), N).mid()) * math.exp(-PI * t / 4) * t * t / float(kink(N))) for t in (lz + 2, lz + 20, 300, 1000, 4000)]
    grid = np.arange(lz + 0.5, 1200, 0.5)
    sg = [np.sign(float(Xi_real(float(t), N).mid())) for t in grid]
    print('E3 N=%d real zeros %d, largest %.3f, 4(N+1)^2=%d; zeros in the last 15 units %s; constant sign on (largest,1200) on a 0.5 grid: %s (sign %d); kink %.4g; t^2 Xi/kink at t=%s: %s' % (N, len(z), lz, 4 * (N + 1) ** 2, np.round(z[-4:], 3), len(set(sg)) == 1, sg[0], float(kink(N)), [round(a, 1) for a, _ in far], [round(b, 4) for _, b in far]))
print('E4 Davenport Heilbronn folded truncations')
K_ = KAPPA_DH
a5 = [0.0, 1.0, K_, -K_, -1.0]
def XD(t, N):
    s = acb(0)
    t = acb(t)
    w = I_ * t / 2
    for n in range(1, N + 1):
        a = a5[n % 5]
        if a == 0:
            continue
        z = pi * n * n / 5
        s = s + arb(a) * n * acb(z).pow(-(acb(arb(3) / 4) + w)) * Gu(acb(arb(3) / 4) + w, z)
    return s.real * arb(PI * float(t.real) / 4).exp()
print('   sum over all n of the first derivatives (should vanish): %s' % (sum((arb(a5[n % 5]) * n * (-pi * n * n / 5).exp() * (3 - 4 * pi * n * n / 5) for n in range(1, 60)), arb(0))).str(5))
for N in (1, 2, 3, 4, 6, 7, 8, 9):
    M = N + 1 if a5[(N + 1) % 5] != 0 else N + 2
    pred = 4 * M * M / 5
    kD = sum((arb(a5[n % 5]) * n * (-pi * n * n / 5).exp() * (3 - 4 * pi * n * n / 5) for n in range(1, N + 1)), arb(0))
    tt = np.arange(0.25, 2.5 * pred + 25, 0.01)
    v = np.array([float(XD(float(t), N).mid()) for t in tt])
    z, _ = sign_zeros(tt, v)
    print('   N=%d real zeros %d, largest %.3f, prediction 4M^2/5=%.1f (M=%d), kink %.4g' % (N, len(z), z[-1] if len(z) else float('nan'), pred, M, float(kD)))
print('E5 Robin family of the wall, a=2pi: g(nu)=-a K\'_{i nu}(a)+sigma K_{i nu}(a), transform of (a cosh t+sigma)exp(-a cosh t)')
a_ = 2 * pi
def Kp(nu):
    """K'_nu(a)=-(K_{nu-1}(a)+K_{nu+1}(a))/2"""
    return -(acb(a_).bessel_k(nu - 1) + acb(a_).bessel_k(nu + 1)) / 2
def g(nu, sigma):
    return -a_ * Kp(nu) + arb(sigma) * acb(a_).bessel_k(nu)
sc = -float(a_ * acb(a_).bessel_k(1).real / acb(a_).bessel_k(0).real)
print('   sigma_c=-aK_1(a)/K_0(a)=%.6f, a=%.6f' % (sc, float(a_)))
for sigma in (0.0, 3.0, -6.5, -6.7, sc, -7.0, -9.0):
    nus = np.arange(0.02, 60.0, 0.01)
    vv = np.array([float((g(acb(0, float(n)), sigma).real * arb(PI * n / 2).exp()).mid()) for n in nus])
    z, _ = sign_zeros(nus, vv)
    mu = np.arange(0.0, 8.0, 0.002)
    vm = np.array([float(g(acb(float(m), 0), sigma).real.mid()) for m in mu])
    zm, _ = sign_zeros(mu, vm)
    # phase at higher levels
    nh = np.arange(300.0, 340.0, 0.01)
    vh = np.array([float((g(acb(0, float(n)), sigma).real * arb(PI * n / 2).exp()).mid()) for n in nh])
    zh, _ = sign_zeros(nh, vh)
    fr = np.mod(theta(2 * zh) / PI, 1.0)
    print('   sigma=%-9.5g nu_1..3=%s; imaginary zeros mu in (0,8): %s; negative kernel at t=0: %s; mean frac(theta/pi) at nu near 320: %.4f' % (sigma, np.round(z[:3], 3), np.round(zm, 4), a_.__float__() + sigma < 0, fr.mean()))
print('E6 count of Robin levels with a signed kernel (sigma=-6.5): rectangle in nu against real zeros')
sigma = -6.5
hh = 0.05
pts = rect_path(0.3, 100.3, -2, 2, hh)
va = np.array([complex(float(x.real.mid()), float(x.imag.mid())) for x in [(g(I_ * acb(complex(p).real, complex(p).imag), sigma) * arb(PI * complex(p).real / 2).exp()) for p in pts]])
d = np.angle(va[1:] / va[:-1])
nus = np.arange(0.3, 100.3, 0.01)
vv = np.array([float((g(acb(0, float(n)), sigma).real * arb(PI * n / 2).exp()).mid()) for n in nus])
print('   rectangle %.3f (max step %.2f), real %d' % (d.sum() / (2 * PI), np.abs(d).max(), len(sign_zeros(nus, vv)[0])))
