"""Track E, the Davenport and Heilbronn control (mpmath, computed, no enclosure).
a = (1, kappa, -kappa, -1, 0) mod 5, kappa = (sqrt(10-2 sqrt 5) - 2)/(sqrt 5 - 1), f(s) = sum a_n n^-s,
Lambda(s) = (5/pi)^{s/2} Gamma((s+1)/2) f(s) = Lambda(1-s); density p_DH(t) proportional to e^{3t/2} theta(e^{2t}),
theta(y) = sum a_n n exp(-pi n^2 y/5). Tests: functional equation, the off line zero, evenness and positivity,
GHS curve (log Lambda)'''(1/2+h), V' convexity, cumulant signs and Hankel determinants of p_m for DH and Polya,
by Cauchy integrals of log Lambda(1/2+z) and log xi(1/2+z) on |z| = 3."""
import mpmath as mp, time
mp.mp.dps = 40
pi = mp.pi
kappa = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
coef = {1: mp.mpf(1), 2: kappa, 3: -kappa, 4: mp.mpf(-1)}
def f(s):
    return mp.mpf(5)**(-s)*mp.fsum(c*mp.zeta(s, mp.mpf(r)/5) for r, c in coef.items())
def Lam(s):
    return (5/pi)**(s/2)*mp.gamma((s+1)/2)*f(s)
def xi(s):
    return mp.mpf(1)/2*s*(s-1)*pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
print("kappa = %s" % mp.nstr(kappa, 12))
s = mp.mpc('0.3', '2.0')
print("functional equation: Lambda(s) = %s, Lambda(1-s) = %s" % (mp.nstr(Lam(s), 12), mp.nstr(Lam(1-s), 12)))
rho = mp.mpc('0.808517', '85.699348')
print("off line zero: |Lambda(rho)| = %s against |Lambda(rho+0.01)| = %s" % (mp.nstr(abs(Lam(rho)), 5), mp.nstr(abs(Lam(rho + mp.mpf('0.01'))), 5)))
Lh = Lam(mp.mpf(1)/2)
print("Lambda(1/2) = %s" % mp.nstr(Lh, 15))
def theta(y, N=80):
    y = mp.mpf(y)
    return mp.fsum(coef[n % 5]*n*mp.exp(-pi*n*n*y/5) for n in range(1, N+1) if n % 5)
norm = 2/(mp.sqrt(5/pi)*Lh)
def pDH(t):
    t = mp.mpf(t)
    return norm*mp.exp(3*t/2)*theta(mp.exp(2*t))
for t in ['0.2', '0.5', '0.8']:
    print("  evenness: p_DH(%s) = %s  p_DH(-%s) = %s" % (t, mp.nstr(pDH(t), 14), t, mp.nstr(pDH('-' + t), 14)))
grid = [mp.mpf(k)/40 for k in range(0, 121)]
vals = [pDH(t) for t in grid]
print("  positivity on [0,3] step 1/40: min p_DH = %s at t = %s; p_DH(0) = %s" % (mp.nstr(min(vals), 6), mp.nstr(grid[vals.index(min(vals))], 4), mp.nstr(vals[0], 8)))
print("  int p_DH = %s" % mp.nstr(2*mp.quad(pDH, [0, 0.5, 1, 2, 3]), 12))

print("GHS: (log Lambda)'''(1/2+h) for DH, (log xi)'''(1/2+h) for Polya")
for h in ['0.05', '0.1', '0.2', '0.45', '0.55', '1', '1.5', '2', '3', '5', '8', '12', '16', '20', '24', '28', '32', '40', '50', '60', '80', '100']:
    hh = mp.mpf(h)
    gD = mp.diff(lambda u: mp.log(Lam(mp.mpf(1)/2 + u)), hh, 3)
    gP = mp.diff(lambda u: mp.log(xi(mp.mpf(1)/2 + u)), hh, 3)
    print("  h=%5s  DH: %s   Polya: %s   Lambda(1/2+h)=%s" % (h, mp.nstr(gD, 8), mp.nstr(gP, 8), mp.nstr(Lam(mp.mpf(1)/2 + hh), 6)))

print("V_DH = -log p_DH: V'(t)/t and V''' on a grid")
V = lambda t: -mp.log(pDH(t))
prev = None
for t in ['0.05', '0.1', '0.2', '0.4', '0.6', '0.8', '1.0', '1.3', '1.6', '2.0']:
    t = mp.mpf(t)
    r = mp.diff(V, t)/t; v3 = mp.diff(V, t, 3)
    print("  t=%s  V'(t)/t=%s %s  V'''=%s" % (mp.nstr(t, 3), mp.nstr(r, 10), "" if prev is None or r > prev else "NOT INCREASING", mp.nstr(v3, 8)))
    prev = r

print("Cumulants by Cauchy integrals on |z| = 3, N = 256 points (64 computed, symmetry), dps 200")
mp.mp.dps = 200
R = mp.mpf(3); N = 256
def coefs(F, kmax):
    t0 = time.time()
    vals = [None]*N
    for j in range(N//4 + 1):
        z = R*mp.expjpi(2*mp.mpf(j)/N)
        vals[j] = F(z)
    for j in range(N//4 + 1):          # F(conj z) = conj F(z), F(-z) = F(z)
        vals[(N//2 - j) % N] = mp.conj(vals[j])      # angle pi - phi: conj(F(-z)) = conj F(z)
        vals[(N - j) % N] = mp.conj(vals[j]); vals[(N//2 + j) % N] = vals[j]   # angle -phi and pi + phi
    c = []
    for k in range(kmax + 1):
        ck = mp.fsum(vals[j]*mp.expjpi(-2*mp.mpf(j)*k/N) for j in range(N))/(N*R**k)
        c.append(ck.real)
    print("  done in %.0f s" % (time.time() - t0))
    return c
half = mp.mpf(1)/2
FD = lambda z: mp.log(Lam(half + z))
FP = lambda z: mp.log(xi(half + z))
cD = coefs(FD, 48); cP = coefs(FP, 48)
print("  check: c_0 (log of the value at the centre) DH %s Polya %s; odd coefficients should vanish: max |c_odd| DH %s Polya %s" % (mp.nstr(cD[0], 10), mp.nstr(cP[0], 10), mp.nstr(max(abs(cD[k]) for k in range(1, 49, 2)), 3), mp.nstr(max(abs(cP[k]) for k in range(1, 49, 2)), 3)))
def pm(c):
    return {m: (-1)**(m+1)*m*c[2*m] for m in range(1, 25)}
pD, pP = pm(cD), pm(cP)
print("  p_m = (-1)^{m+1} m kappa_{2m}/(2m)!  (m, DH, Polya); all should be positive if all zeros are on the line")
for m in range(1, 25):
    print("   m=%2d  DH %s   Polya %s" % (m, mp.nstr(pD[m], 8), mp.nstr(pP[m], 8)))
print("  Hankel determinants H_n = det(p_{2+i+j})_{i,j<n}, raw values and signs")
for n in range(1, 12):
    HD = mp.det(mp.matrix([[pD[2+i+j] for j in range(n)] for i in range(n)]))
    HP = mp.det(mp.matrix([[pP[2+i+j] for j in range(n)] for i in range(n)]))
    print("   n=%2d  DH %s (%s)   Polya %s (%s)" % (n, mp.nstr(HD, 6), "+" if HD > 0 else "-", mp.nstr(HP, 6), "+" if HP > 0 else "-"))
print("  zero counts to height 85.7 (main terms): zeta (T/2pi)log(T/2pi e) = %s, DH conductor 5 (T/2pi)log(5T/2pi e) = %s" % (mp.nstr(mp.mpf('85.7')/(2*pi)*mp.log(mp.mpf('85.7')/(2*pi*mp.e)), 4), mp.nstr(mp.mpf('85.7')/(2*pi)*mp.log(5*mp.mpf('85.7')/(2*pi*mp.e)), 4)))
