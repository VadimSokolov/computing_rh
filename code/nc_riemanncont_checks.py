"""Checks for ch/nc_riemanncont.tex (Riemann's continuation and the GGC question).
All computations in mpmath floating point (computed, not certified)."""
import mpmath as mp
mp.mp.dps = 40
pi = mp.pi

def xi(s):
    if s == 0 or s == 1:
        return mp.mpf(1)/2
    return s*(s-1)/2*pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)

print("== 1. Riemann's kernel as Gamma(2) minus Gamma(1) mixtures ==")
def mellin_series(s, N=4000):
    # sum_n (pi n^2)^{-s/2} (4 Gamma(2+s/2) - 6 Gamma(1+s/2)), Re s > 2
    return (4*mp.gamma(2+s/2)-6*mp.gamma(1+s/2))*pi**(-s/2)*mp.zeta(s)
for s in [mp.mpf(3), mp.mpc(4,2), mp.mpf(6)]:
    print("s=", s, " series:", mp.nstr(mellin_series(s),20), " 2xi(s):", mp.nstr(2*xi(s),20))
def fV(x, N=120):
    # density of V = pi S_2/2
    return mp.fsum(4*(pi*n*n)**2*x*mp.exp(-pi*n*n*x) - 6*pi*n*n*mp.exp(-pi*n*n*x) for n in range(1, N+1))
mp.mp.dps = 60
for s in [mp.mpf(0), mp.mpf(1), mp.mpc(0.5, 3)]:
    I = mp.quad(lambda x: x**(s/2)*fV(x), [mp.mpf('0.02'), 0.1, 0.5, 2, 10, 40])
    print("Mellin of f_V at s=", s, ":", mp.nstr(I, 18), " 2xi(s):", mp.nstr(2*xi(s), 18))
print("min of f_V on grid [0.03,5]:", mp.nstr(min(fV(mp.mpf(k)/100) for k in range(3, 501)), 6))
mp.mp.dps = 40

print("== 2. Logistic, Gumbel and the sinh clock ==")
for t in [mp.mpf('0.7'), mp.mpf(2), mp.mpf(5)]:
    cf_log = pi*t/mp.sinh(pi*t)
    gg = mp.gamma(1-1j*t)*mp.gamma(1+1j*t)           # difference of two Gumbels
    prod = mp.nprod(lambda k: 1/(1+t**2/k**2), [1, mp.inf])
    print("t=", t, mp.nstr(cf_log, 18), mp.nstr(mp.re(gg), 18), mp.nstr(prod, 18))
for s in [mp.mpf('0.5'), mp.mpf(3)]:
    lhs = mp.gamma(1+s)*mp.exp(mp.euler*s)
    rhs = mp.nprod(lambda k: mp.exp(s/k)/(1+s/k), [1, mp.inf])
    print("Gumbel compensated product s=", s, mp.nstr(lhs, 18), mp.nstr(rhs, 18))

print("== 3. The sine squared image of the positive half ==")
NZ = 400
gam = [mp.im(mp.zetazero(k)) for k in range(1, NZ+1)]
def K(x):
    return mp.fsum(mp.exp(-g*x) for g in gam)
def m_series(th):
    return mp.fsum(th**2/(g*(g*g+th**2)) for g in gam)/pi
for th in [mp.mpf(5), mp.mpf('14.134725'), mp.mpf(30)]:
    q = mp.quad(lambda x: (1-mp.cos(x*th))*K(x), [0, mp.mpf(1)/th, 0.2, 1, 3, 8])/pi
    print("theta=", th, " quad:", mp.nstr(q, 15), " series:", mp.nstr(m_series(th), 15))
# contribution of a hypothetical off line pair c = gamma - i delta (and conjugate)
for (g, d) in [(mp.mpf(1000), mp.mpf('0.4')), (mp.mpf(100), mp.mpf('0.45'))]:
    c = mp.mpc(g, -d)
    f = lambda th: 2*mp.re(th**2/(c*(c*c+th**2)))/pi
    vals = [f(mp.mpf(k)/10) for k in range(1, 30000)]
    print("off line pair gamma=", g, "delta=", d, " min of contribution on (0,3000):", mp.nstr(min(vals), 8),
          " value at theta=gamma:", mp.nstr(f(g), 8))
    a = c*c
    print("   pole of m(theta)/theta^2 in w=theta^2 at w=-a=", mp.nstr(-a, 10))

print("== 4. Tilts of the kernel ==")
def Phi(u):
    u = abs(u)
    return 2*mp.fsum((2*pi**2*n**4*mp.exp(4.5*u)-3*pi*n*n*mp.exp(2.5*u))*mp.exp(-pi*n*n*mp.exp(2*u)) for n in range(1, 12))
Phi0 = Phi(0)
print("Phi(0) =", mp.nstr(Phi0, 20))
for u in [mp.mpf('0.3'), mp.mpf(2)]:
    lhs = mp.exp(-u)
    rhs = mp.quad(lambda t: (pi*t)**mp.mpf(-0.5)*mp.exp(-t)*mp.exp(-u*u/(4*t)), [0, 1, mp.inf])
    print("e^{-|u|} as normal scale mixture, u=", u, mp.nstr(lhs, 18), mp.nstr(rhs, 18))
mp.mp.dps = 30
def F(x):
    return 2*mp.quad(lambda u: mp.exp(-u)*Phi(u)*mp.cos(x*u), [0, 0.5, 1, 1.5, 2.5])
xs = [mp.mpf(k)/4 for k in range(0, 4*80+1)]
Fv = [F(x) for x in xs]
ch = [mp.nstr(xs[i], 6) for i in range(len(xs)-1) if Fv[i]*Fv[i+1] < 0]
print("sign changes of F on [0,80]:", ch)
for x in [mp.mpf(40), mp.mpf(80)]:
    print("x=", x, " F(x)=", mp.nstr(F(x), 10), " 2Phi(0)/x^2=", mp.nstr(2*Phi0/x**2, 10))
def Xi(z):
    return xi(mp.mpf(1)/2+1j*z)
lam = mp.mpf('0.05')
z = mp.mpf(10)
lhs = 2*mp.quad(lambda u: Phi(u)*mp.exp(lam*u*u)*mp.cos(z*u), [0, 0.5, 1, 1.5, 2.5])
rhs = mp.quad(lambda y: Xi(z+1j*mp.sqrt(2*lam)*y)*mp.exp(-y*y/2)/mp.sqrt(2*pi), [-12, 0, 12])
print("de Bruijn tilt lambda=0.05 at z=10: kernel side", mp.nstr(lhs, 15), " vertical Gaussian average", mp.nstr(mp.re(rhs), 15))
