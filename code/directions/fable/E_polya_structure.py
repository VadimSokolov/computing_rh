"""Track E, structure of the Polya law (mpmath, computed, no enclosure).
(a) Kennedy's law of the excursion maximum M: E[M^s] = 2 (pi/2)^{s/2} xi(s).
(b) Zero free Euler form: p(t) = sum_n n^{-1/2} p1(|t| + log n), p1 >= 0 on [0,inf); kink identities for
    P-smooth truncations; TV of the smooth truncations.
(c) GHS for the Polya law: (log xi)'''(1/2+h) on a grid; cross check with the Hadamard sum over 200 zeros
    (the cross check uses zeros, nothing else does); Stirling asymptotics.
(d) V = -log Phi: V'(t)/t increasing, V' convex, on a grid.
(e) The law in L that violates GHS: L(z) = e^{z^2/2}(1+z^2/c^2)."""
import mpmath as mp, time
mp.mp.dps = 40
pi = mp.pi
def xi(s):
    s = mp.mpf(s) if not isinstance(s, mp.mpc) else s
    return mp.mpf(1)/2*s*(s-1)*pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
xi_half = xi(mp.mpf(1)/2)
print("xi(1/2) = %s" % mp.nstr(xi_half, 15))

print("(a) Kennedy: E[M^s] by quadrature against 2 (pi/2)^{s/2} xi(s)")
def fM(x, N=200):
    return 8*mp.fsum((4*n**4*x**3 - 3*n**2*x)*mp.exp(-2*n**2*x**2) for n in range(1, N+1))
for s in [mp.mpf(1)/2, mp.mpf('2.3'), mp.mpf(4)]:
    lhs = mp.quad(lambda x: x**s*fM(x), [0.02, 0.3, 1, 3, 10])
    rhs = 2*(pi/2)**(s/2)*xi(s)
    print("  s=%s  E[M^s]=%s  2(pi/2)^{s/2}xi(s)=%s" % (mp.nstr(s, 3), mp.nstr(lhs, 14), mp.nstr(rhs, 14)))

print("(b) zero free Euler form")
def p1(v):
    v = mp.mpf(v)
    return (4*pi**2*mp.exp(9*v/2) - 6*pi*mp.exp(5*v/2))*mp.exp(-pi*mp.exp(2*v))/xi_half
def p(t, N=60):
    t = abs(mp.mpf(t))
    return mp.fsum(mp.sqrt(n)**-1*p1(t + mp.log(n)) for n in range(1, N+1))
def Phi_over_xi(t, N=60):
    t = mp.mpf(t)
    return 2*mp.fsum((2*pi**2*n**4*mp.exp(9*t/2) - 3*pi*n**2*mp.exp(5*t/2))*mp.exp(-pi*n**2*mp.exp(2*t)) for n in range(1, N+1))/xi_half
print("  sign change of p1 at v* = (1/2)log(3/(2pi)) = %s; p1(0) = %s" % (mp.nstr(mp.log(3/(2*pi))/2, 8), mp.nstr(p1(0), 10)))
print("  mode of p1: v_mode = %s; p1'(0.01) = %s (negative: decreasing from 0.01 on)" % (mp.nstr(mp.findroot(lambda v: mp.diff(p1, v), 0.0), 8), mp.nstr(mp.diff(p1, mp.mpf('0.01')), 6)))
for t in ['0.3', '-0.3', '1.0']:
    print("  t=%s: sum_n n^-1/2 p1(|t|+log n) = %s   Phi(t)/xi(1/2) = %s" % (t, mp.nstr(p(t), 15), mp.nstr(Phi_over_xi(t), 15)))
print("  int_0^inf p = %s (should be 1/2)" % mp.nstr(mp.quad(lambda t: p(t), [0, 0.5, 1, 2, 3]), 12))
dp1 = lambda v: mp.diff(p1, v)
kink_all = mp.fsum(mp.sqrt(n)**-1*dp1(mp.log(n)) for n in range(1, 61))
print("  kink of the full sum: sum_n n^-1/2 p1'(log n) = %s (Jacobi: 0);  n=1 term alone p1'(0) = %s" % (mp.nstr(kink_all, 6), mp.nstr(dp1(0), 8)))
def smooth(P, N=2000):
    out = []
    for n in range(1, N+1):
        m = n
        for q in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
            if q > P: break
            while m % q == 0: m //= q
        if m == 1: out.append(n)
    return out
tail = lambda n: mp.quad(lambda v: p1(v), [mp.log(n), mp.log(n)+1, mp.log(n)+3])
print("  P-smooth truncations p^P (normalised): P, number of n <= 2000, kink (p^P)'(0+), TV(p^P, p)")
for P in [1, 2, 3, 5, 7, 11, 13]:
    S = smooth(P)
    mass = 2*mp.fsum(mp.sqrt(n)**-1*tail(n) for n in S if n <= 40)
    c = 1/mass
    kink = c*mp.fsum(mp.sqrt(n)**-1*dp1(mp.log(n)) for n in S if n <= 40)
    pP = lambda t: c*mp.fsum(mp.sqrt(n)**-1*p1(abs(t) + mp.log(n)) for n in S if n <= 40)
    tv = mp.quad(lambda t: abs(pP(t) - p(t)), [0, 0.3, 0.6, 1, 1.5, 2.5])
    print("  P=%2d  #=%4d  kink=%s  TV=%s" % (P, len(S), mp.nstr(kink, 6), mp.nstr(tv, 4)))

print("(c) GHS: g3(h) = (log xi)'''(1/2+h); Hadamard sum over 200 zeros + tail; -1/(2h^2)")
t0 = time.time()
with mp.workdps(20):
    zeros = [mp.mpf(mp.zetazero(k).imag) for k in range(1, 201)]
print("  200 zeros in %.0f s, gamma_200 = %s" % (time.time()-t0, mp.nstr(zeros[-1], 8)))
def g3(h):
    h = mp.mpf(h)
    return mp.diff(lambda u: mp.log(xi(mp.mpf(1)/2 + u)), h, 3)
def had(h):
    h = mp.mpf(h)
    s = mp.fsum(-4*h*(3*g*g - h*h)/(g*g + h*h)**3 for g in zeros)
    T = zeros[-1]
    tail = mp.quad(lambda g: -4*h*(3*g*g - h*h)/(g*g + h*h)**3*mp.log(g/(2*pi))/(2*pi), [T, 2*T, 10*T, mp.inf])
    return s + tail
for h in ['0.05', '0.1', '0.2', '0.45', '0.55', '1', '1.5', '2', '3', '5', '8', '12', '16', '20', '24', '28', '32', '40', '50', '60']:
    print("  h=%5s  g3=%s  hadamard=%s  -1/(2h^2)=%s" % (h, mp.nstr(g3(h), 8), mp.nstr(had(h), 8), mp.nstr(-1/(2*mp.mpf(h)**2), 5)))

print("(d) V = -log Phi: V'(t)/t and V''' on a grid")
V = lambda t: -mp.log(Phi_over_xi(t))
prev = None
for t in ['0.05', '0.1', '0.2', '0.4', '0.6', '0.8', '1.0', '1.3', '1.6', '2.0', '2.4']:
    t = mp.mpf(t)
    r = mp.diff(V, t)/t; v3 = mp.diff(V, t, 3)
    print("  t=%s  V'(t)/t=%s %s  V'''=%s" % (mp.nstr(t, 3), mp.nstr(r, 10), "" if prev is None or r > prev else "NOT INCREASING", mp.nstr(v3, 8)))
    prev = r

print("(e) L(z)=e^{z^2/2}(1+z^2/c^2), c=1: (log L)''' at h = 1, 1.5, 2, 3")
for h in [1, 1.5, 2, 3]:
    h = mp.mpf(h)
    print("  h=%s  (log L)'''=%s" % (h, mp.nstr(4*h*(h*h - 3)/(1 + h*h)**3, 6)))
