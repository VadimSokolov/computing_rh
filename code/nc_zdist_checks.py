"""Checks for Chapter nc_zdist: z distributions, Meixner duality and the reciprocal Xi function.

Run: python3 code/nc_zdist_checks.py
All computations in mpmath at 30 significant digits (floating point; nothing here is certified).
"""
from mpmath import mp, mpf, gamma, beta, quad, inf, exp, log, pi, zeta, digamma, cos, cosh, sqrt, re, zetazero, nprod, fabs

mp.dps = 30


def xi(s):
    return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def Xi(t):
    return re(xi(mpf(1) / 2 + 1j * t))


def show(name, a, b):
    a, b = mpf(re(a)), mpf(re(b))
    print("%-48s %s %s  rel.diff %.2e" % (name, mp.nstr(a, 15), mp.nstr(b, 15), float(fabs(a / b - 1))))


# 1. z law: E e^{sX} = B(a+s,b-s)/B(a,b)
a, b, s = mpf("0.7"), mpf("1.9"), mpf("0.4")
f = lambda x: exp(a * x) / (beta(a, b) * (1 + exp(x)) ** (a + b))
show("1 z mgf, integral vs Beta ratio", quad(lambda x: f(x) * exp(s * x), [-inf, 0, inf]), beta(a + s, b - s) / beta(a, b))

# 2. variance mean mixture: product over (a+k)(b+k)
prod = nprod(lambda k: (a + k) * (b + k) / ((a + k + s) * (b + k - s)), [0, inf])
show("2 mixing product (a+k)(b+k) vs Beta ratio", prod, beta(a + s, b - s) / beta(a, b))
mu = (a - b) / 2
prod2 = nprod(lambda k: (a + k) * (b + k) / ((a + k) * (b + k) - 2 * mu * s - s ** 2), [0, inf])
show("2b E exp(W(mu s+s^2/2)), W atoms (a+k)(b+k)", prod2, beta(a + s, b - s) / beta(a, b))

# 3. Meixner duality: char fn of z(1/4,1/4) vs |Gamma(1/4+iv)|^2/Gamma(1/4)^2 vs rz2 Meixner density
q = mpf(1) / 4
fz = lambda x: exp(q * x) / (beta(q, q) * (1 + exp(x)) ** (2 * q))
v = mpf("0.8")
cf = quad(lambda x: fz(x) * cos(v * x), [-inf, 0, inf])
show("3a z(1/4,1/4) char fn vs |Gamma|^2 ratio", cf, abs(gamma(q + 1j * v)) ** 2 / gamma(q) ** 2)
p = lambda tau: abs(gamma(q + 1j * tau / 2)) ** 2 / (2 ** 1.5 * pi ** 1.5)
show("3b same vs 2^{3/2}pi^{3/2} p(2v)/Gamma(1/4)^2", cf, 2 ** 1.5 * pi ** 1.5 * p(2 * v) / gamma(q) ** 2)
print("   B(1/4,1/4) =", mp.nstr(beta(q, q), 15))

# 4. X/2 = sum L_k/c_k, c_k = 2k+1/2
s4 = mpf(3)
show("4 prod (1+s^2/c_k^2)^-1 vs Gamma ratio at s=3", nprod(lambda k: 1 / (1 + s4 ** 2 / (2 * k + q * 2) ** 2), [0, inf]),
     abs(gamma(q + 1j * s4 / 2)) ** 2 / gamma(q) ** 2)
print("   Var(X/2) = psi'(1/4)/2 =", mp.nstr(mp.psi(1, q) / 2, 12))

# 5. archimedean identity on the critical line
s5 = mpf(5)
lhs = abs(zeta(mpf(1) / 2 + 1j * s5) / zeta(mpf(1) / 2)) ** 2
phi = abs(gamma(q + 1j * s5 / 2)) ** 2 / gamma(q) ** 2
rhs = (Xi(s5) / Xi(0)) ** 2 * (q / (s5 ** 2 + q)) ** 2 / phi
show("5 |zeta ratio|^2 vs Xi^2, rational, 1/phi_X", lhs, rhs)

# zeros
N = 300
gam = [zetazero(n).imag for n in range(1, N + 1)]
G = gam[-1]
ET_exact = mp.diff(lambda x: xi(x), mpf(1) / 2, 2) / (2 * xi(mpf(1) / 2))   # E T = sum_k 1/a_k
print("   E T = xi''(1/2)/(2 xi(1/2)) =", mp.nstr(ET_exact, 12))
tail2 = ET_exact - sum(1 / g ** 2 for g in gam)   # sum_{gamma>G} gamma^-2, using RH below gamma_300 (verified)
tail4 = log(G / (2 * pi)) / (2 * pi * 3 * G ** 3)  # Weyl approximation to sum_{gamma>G} gamma^-4

# 6. Xi(0)/Xi(s) = E e^{s^2 T} = prod (1-s^2/gamma^2)^-1 (Hadamard), s=5
s6 = mpf(5)
pr = mpf(1)
for g in gam:
    pr *= 1 / (1 - s6 ** 2 / g ** 2)
pr *= exp(s6 ** 2 * tail2 + s6 ** 4 * tail4 / 2)
show("6 Xi(0)/Xi(5) vs product over 300 zeros+tail", Xi(0) / Xi(s6), pr)

# 7. tilted mixing law: Xi(beta)/Xi(sqrt(beta^2-u)) = prod (1+u/(gamma^2-beta^2))^-1
bt, u = mpf(5), mpf(3)
pr = mpf(1)
for g in gam:
    pr *= 1 / (1 + u / (g ** 2 - bt ** 2))
pr *= exp(-u * (tail2 + bt ** 2 * tail4))
show("7 Xi(5)/Xi(sqrt(22)) vs tilted atoms", Xi(bt) / Xi(sqrt(bt ** 2 - u)), pr)

# 8. 2 sin^2 lemma for the gamma function
a8, u8 = q, mpf(2)
m = lambda th: re(digamma(a8 + 1j * th)) - digamma(a8)
I = quad(lambda th: log(1 + u8 ** 2 / th ** 2) * m(th), [0, 1, 10, 100, inf])
show("8 e^{psi(a)u}Gamma(a)/Gamma(a+u) vs Thorin form", exp(digamma(a8) * u8) * gamma(a8) / gamma(a8 + u8), exp(-I / pi))
a8b = mpf(3)
m = lambda th: re(digamma(a8b + 1j * th)) - digamma(a8b)
I = quad(lambda th: log(1 + u8 ** 2 / th ** 2) * m(th), [0, 1, 10, 100, inf])
show("8b same at a=3", exp(digamma(a8b) * u8) * gamma(a8b) / gamma(a8b + u8), exp(-I / pi))

# 9. variance of Y and of the partial sums over zeros
ET = sum(1 / g ** 2 for g in gam) + tail2
print("9 2*sum gamma^-2 (300 zeros + Weyl tail) =", mp.nstr(2 * ET, 8), " against 2 E T = 0.046210")


# 10. heat flow: R_t(1) / R_0(1) - 1 at t = 0.2, from the Gaussian smoothing of xi on the real axis
def Ecs(sig, t):
    c = sqrt(t / 2)
    return quad(lambda z: exp(-z ** 2 / 2) / sqrt(2 * pi) * re(xi(mpf(1) / 2 + sig + c * z)), [-inf, -3, 0, 3, inf])


t = mpf("0.2")
R02 = Ecs(0, t) / Ecs(1, t)
R00 = xi(mpf(1) / 2) / xi(mpf(3) / 2)
print("10 R_0.2(1)/R_0(1) - 1 =", mp.nstr(R02 / R00 - 1, 6))

# 11. Cauchy's beta integral: first member of the z van Dantzig pair, c = 0.8, v = 0.7
c, v = mpf("0.8"), mpf("0.7")
lhs = quad(lambda th: cos(th) ** (2 * c - 2) * cos(2 * v * th), [-pi / 2, 0, pi / 2])
lhs0 = quad(lambda th: cos(th) ** (2 * c - 2), [-pi / 2, 0, pi / 2])
show("11 E e^{2iv Theta} vs Gamma(c)^2/(Gamma(c+v)Gamma(c-v))", lhs / lhs0, gamma(c) ** 2 / (gamma(c + v) * gamma(c - v)))
for c in [mpf("0.25")]:
    vals = [gamma(c) ** 2 / (gamma(c + v) * gamma(c - v)) for v in [mpf(10), mpf("10.5"), mpf(100), mpf("100.5")]]
    print("   c=1/4: Gamma(c)^2/(Gamma(c+v)Gamma(c-v)) at v=10,10.5,100,100.5:", [mp.nstr(x, 6) for x in vals])

# 12. real axis form of the archimedean factor, s = 0.3
s = mpf("0.3")
lhs = xi(mpf(1) / 2 + s) ** 2 / xi(mpf(1) / 2) ** 2
q = mpf(1) / 4
EzX = gamma(q + s / 2) * gamma(q - s / 2) / gamma(q) ** 2
rhs = (1 - 4 * s ** 2) ** 2 * EzX * zeta(mpf(1) / 2 + s) * zeta(mpf(1) / 2 - s) / zeta(mpf(1) / 2) ** 2
show("12 xi(1/2+s)^2/xi(1/2)^2 vs (1-4s^2)^2 E e^{sX/2} zeta zeta", lhs, rhs)

# 13. Meixner density 2^{2c}|Gamma(c+ix)|^2/(2 pi Gamma(2c)) has char fn cosh(u/2)^{-2c}
c, uu = mpf("0.7"), mpf("1.1")
pc = lambda x: 2 ** (2 * c) * abs(gamma(c + 1j * x)) ** 2 / (2 * pi * gamma(2 * c))
show("13a Meixner(c) mass", quad(pc, [-inf, 0, inf]), 1)
show("13b Meixner(c) char fn vs cosh(u/2)^{-2c}", quad(lambda x: pc(x) * cos(uu * x), [-inf, 0, inf]), cosh(uu / 2) ** (-2 * c))
show("13c Meixner char fn vs product over (2k+1)pi", nprod(lambda k: (1 + uu ** 2 / ((2 * k + 1) * pi) ** 2) ** (-2 * c), [0, inf]), cosh(uu / 2) ** (-2 * c))

# 14. the constant of the line family bound
T0 = mpf("3e12")
print("14 2(1+log T0)/(pi T0) =", mp.nstr(2 * (1 + log(T0)) / (pi * T0), 6))

# 15. mean of D_0 = T + T' - H_inf two ways
lz2 = mp.diff(lambda x: log(zeta(x)), mpf(1) / 2, 2)
print("15 2ET - psi'(1/4)/4 =", mp.nstr(2 * ET_exact - mp.psi(1, mpf(1) / 4) / 4, 12),
      "  -8 + (log zeta)''(1/2) =", mp.nstr(-8 + lz2, 12), "  E H_inf =", mp.nstr(mp.psi(1, mpf(1) / 4) / 4, 10))

# 16. Laplace transform of D_0 at w = 0.1: 1/G(w)^2 * E e^{w H_inf} vs zeta expression
w = mpf("0.1")
lhs = (xi(mpf(1) / 2) / xi(mpf(1) / 2 + sqrt(w))) ** 2 * gamma(mpf(1) / 4 + sqrt(w) / 2) * gamma(mpf(1) / 4 - sqrt(w) / 2) / gamma(mpf(1) / 4) ** 2
rhs = zeta(mpf(1) / 2) ** 2 / ((1 - 4 * w) ** 2 * zeta(mpf(1) / 2 + sqrt(w)) * zeta(mpf(1) / 2 - sqrt(w)))
show("16 E e^{-w D_0} two ways at w=0.1", lhs, rhs)

# 17. the line family: error constant at T' = T0 and m(T') ; partial check at T' = gamma_300
def eps(Tp):
    return 2 / pi * (log(Tp) / 3 + mpf(1) / 9) / Tp ** 3
print("17 eps(T0) =", mp.nstr(eps(T0), 6), " eps(gamma_300) =", mp.nstr(eps(G), 6), " gamma_300 =", mp.nstr(G, 10))
mT = ET_exact - sum(1 / g ** 2 for g in gam)
for u in [mpf(5), mpf(20)]:
    approx = exp(-u ** 2 * mT)
    for g in gam:
        approx /= (1 + u ** 2 / g ** 2)
    exact = xi(mpf(1) / 2) / xi(mpf(1) / 2 + u)
    print("   u=%s: |log ratio| = %.3e, bound u^4 eps = %.3e" % (u, float(abs(log(exact / approx))), float(u ** 4 * eps(G))))
print("   m(gamma_300) =", mp.nstr(mT, 8))

# 18. Malmsten's L\'evy form of log Gamma used in the 2 sin^2 section, c = 0.3, s = 1.7 (50 digits, expm1 form)
mp.dps = 50
from mpmath import expm1, loggamma
c, s = mpf("0.3"), mpf("1.7")
fm = lambda x: (expm1(-s * x) + s * x) * exp(-c * x) / (x * (-expm1(-x)))
show("18 Malmsten: logGamma ratio vs Levy integral", loggamma(c + s) - loggamma(c) - s * digamma(c),
     quad(fm, [0, 0.5, 1, 2, 5, 10, 20, 40, 80, 160, inf]))
