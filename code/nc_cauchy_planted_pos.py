# Is the planted function F the characteristic function of a positive law? f(x) = (1/pi) int_0^inf F(t) cos(xt) dt,
# compared with the Polya kernel version p(x) from Xi. Floating point, mpmath.
import mpmath as mp
mp.mp.dps = 20
g1, g2 = mp.im(mp.zetazero(1)), mp.im(mp.zetazero(2)); m = (g1+g2)/2; eta = mp.mpf('0.3')
def Xi(z):
    s = mp.mpf(1)/2 + 1j*z
    return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
c = m + 1j*eta
def R(t): return mp.re((t*t-c*c)*(t*t-mp.conj(c)**2))/((t*t-g1*g1)*(t*t-g2*g2))
X0 = Xi(0)
def dens(x, planted):
    f = lambda t: (Xi(t)*(R(t) if planted else 1))*mp.cos(x*t)
    pts = [0, 10, g1, 17.6, g2, 30, 45, 60, 80, 110]
    return mp.quad(f, pts)/(mp.pi*X0)
for x in [0, 0.25, 0.5, 0.75, 1.0, 1.25]:
    a, b = dens(x, False), dens(x, True)
    print(f"x={x}: Polya {mp.nstr(a,8)}  planted {mp.nstr(b,8)}  ratio {mp.nstr(b/a,6)}")
