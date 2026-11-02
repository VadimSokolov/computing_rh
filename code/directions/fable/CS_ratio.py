# Track CS. The ratio Phi/rho* (is Phi = rho* times a coupling multiplier?), the normalising
# constants, and the total variation between the Polya law and the cosh spin. mpmath, "computed".
from mpmath import mp, mpf, pi, exp, cosh, quad, besselk, diff
mp.dps = 30

def Phi(u):
    u = mpf(u); s = mpf(0)
    for n in range(1, 40):
        s += (2 * pi**2 * n**4 * exp(9 * u) - 3 * pi * n**2 * exp(5 * u)) * exp(-pi * n**2 * exp(4 * u))
    return s
def rho_star(u):
    u = mpf(u)
    return (2 * pi**2 * cosh(9 * u) - 3 * pi * cosh(5 * u)) * exp(-2 * pi * cosh(4 * u))

ZPhi = 2 * quad(Phi, [0, 1, 3])
Zstar = 2 * quad(rho_star, [0, 1, 3])
# closed form for Zstar: 2 pi^2 * (1/4) * 2 K_{9/4}(2 pi) - 3 pi * (1/4) * 2 K_{5/4}(2 pi)
Zstar_cf = pi**2 * besselk(mpf(9) / 4, 2 * pi) - (3 * pi / 2) * besselk(mpf(5) / 4, 2 * pi)
print("Z_Phi = %s (H_0(0) = 0.1242802), Z_star = %s, closed form %s" % (mp.nstr(ZPhi, 10), mp.nstr(Zstar, 10), mp.nstr(Zstar_cf, 10)))
print("u, Phi/rho*, p_V/p* (normalised), d/du log(p_V/p*)")
prev = None; mono = True
for k in range(0, 31):
    u = mpf(k) / 10
    r = Phi(u) / rho_star(u)
    rn = r * Zstar / ZPhi
    d = diff(lambda v: mp.log(Phi(v) / rho_star(v)), u)
    if prev is not None and r > prev: mono = False
    prev = r
    print("  %.1f  %s  %s  %s" % (float(u), mp.nstr(r, 8), mp.nstr(rn, 8), mp.nstr(d, 6)))
print("ratio Phi/rho* decreasing on the grid [0, 3]:", mono, "; limit at infinity: 2 (n = 1 term over pi^2 e^{9u - pi e^{4u}})")
tv = quad(lambda u: abs(Phi(u) / ZPhi - rho_star(u) / Zstar), [0, 0.2, 0.5, 1, 3])
print("total variation between the Polya law and rho*: %s" % mp.nstr(tv, 6))
# variances
vV = 2 * quad(lambda u: u**2 * Phi(u), [0, 1, 3]) / ZPhi
vS = 2 * quad(lambda u: u**2 * rho_star(u), [0, 1, 3]) / Zstar
print("Var V = %s, Var rho* = %s" % (mp.nstr(vV, 8), mp.nstr(vS, 8)))
print("U*''(0) = %s" % mp.nstr(diff(lambda v: -mp.log(rho_star(v)), 0, 2), 8))
