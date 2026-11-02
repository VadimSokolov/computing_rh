#!/usr/bin/env python3
# Track E. Exact identities: Euler factors as splines, Renyi type identity for the range, Euler exponents of
# the theta functions, Laplace transforms for odd periodic coefficients (DH and real characters),
# PF minors in balls, Lebowitz moment tests.
# Run: timeout 300 python3 -I E_spline.py
import sys, math
sys.path.insert(0, '/tmp/e_lib')
from flint import arb, acb, arb_mat, ctx
from fractions import Fraction
ctx.prec = 256
pi = arb.pi()

print('--- (a) hyperbolic spline: (1/4+th^2) int g_p(x) e^{ith x} dx = 1+1/p-2 p^{-1/2} cos(th log p)')
def gp_transform_quad(p, th, N=4000):
    # composite Gauss free check: exact antiderivative is used in E.md; here a Simpson check in balls of the midpoint only
    l = arb(p).log()
    h = l / N
    s = arb(0)
    for k in range(N + 1):
        x = h * k
        w = 1 if k in (0, N) else (4 if k % 2 else 2)
        s += w * ((l - x) / 2).sinh() * (arb(th) * x).cos()
    return 2 * s * h / 3 * 2 * arb(p).sqrt().__rtruediv__(1)
for p, th in ((2, 3.0), (3, 14.0), (7, 21.0)):
    lhs = (arb(1) / 4 + arb(th)**2) * gp_transform_quad(p, th)
    rhs = 1 + arb(1) / p - 2 / arb(p).sqrt() * (arb(th) * arb(p).log()).cos()
    print(' p=%d th=%g  lhs=%.12f rhs=%.12f (Simpson, no remainder bound)' % (p, th, float(lhs), float(rhs)))

print('--- (b) Renyi type identity: E exp(-u Y^2) = (sqrt(pi u)/sinh sqrt(pi u))^2')
def F_Y_laplace(u, N=60, T=40.0, M=20000):
    # u int_0^inf e^{-ux} F(x) dx with F(x)=sum_{n in Z}(1-2 pi n^2 x) e^{-pi n^2 x}: closed partial fractions, absolutely convergent
    s = arb(1)
    for n in range(1, N + 1):
        a = pi * n**2
        s += 2 * arb(u) * (1 / (arb(u) + a) - 2 * a / (arb(u) + a)**2)
    # tail: terms are u(u-a)/(u+a)^2*2 ~ -2u/a; bound the tail by the integral of 2u/(pi n^2)
    tail = 2 * arb(u) / (pi * N)
    return s + arb(0, 1) * tail * 2
for u in (0.5, 1.0, 5.0):
    r = (pi * u).sqrt()
    print(' u=%g  partial fractions %s  closed form %s' % (u, F_Y_laplace(u, 2000).str(6), ((r / r.sinh())**2).str(12)))

print('--- (c) Euler exponents: f(q)=prod_n (1-q^n)^{a_n}')
def euler_exponents(coeffs, K):
    # coeffs: f=1+sum c_n q^n, returns a_1..a_K with f=prod(1-q^n)^{a_n}, exact Fractions
    f = [Fraction(0)] * (K + 1)
    f[0] = Fraction(1)
    for n, c in coeffs.items():
        if n <= K:
            f[n] = Fraction(c)
    # log derivative: q f'/f = -sum_n n a_n q^n/(1-q^n) = -sum_m q^m sum_{d|m} d a_d
    b = [Fraction(0)] * (K + 1)  # b_m = coefficient of q^m in q f'/f
    for m in range(1, K + 1):
        s = m * f[m]
        for j in range(1, m):
            s -= b[j] * f[m - j]
        b[m] = s
    a = [Fraction(0)] * (K + 1)
    for m in range(1, K + 1):
        s = -b[m]
        for d in range(1, m):
            if m % d == 0:
                s -= d * a[d]
        a[m] = s / m
    return [int(x) if x.denominator == 1 else x for x in a[1:]]
K = 24
print(' theta_3 = sum q^{n^2}          :', euler_exponents({n*n: 2 for n in range(1, 10)}, K))
print(' theta_4 = sum (-1)^n q^{n^2}   :', euler_exponents({n*n: 2*(-1)**n for n in range(1, 10)}, K))

print('--- (d) odd periodic coefficients: Laplace transform of theta_a(x)=sum a_n n exp(-pi n^2 x/q)')
def analyse(name, q, a):
    # a: dict j->a_j for 1<=j<q, a_{q-j}=-a_j. Numerator of R(t)=sum_j a_j*2 sin(2 pi j/q)/(t-cos(2 pi j/q)), j<q/2.
    import numpy as np
    js = [j for j in range(1, (q + 1) // 2) if abs(a.get(j, 0)) > 1e-14]
    s = {j: 2 * a[j] * math.sin(2 * math.pi * j / q) for j in js}
    c = {j: math.cos(2 * math.pi * j / q) for j in js}
    num = np.poly1d([0.0])
    for j in js:
        prod = np.poly1d([1.0])
        for k in js:
            if k != j:
                prod = prod * np.poly1d([1.0, -c[k]])
        num = num + s[j] * prod
    roots = num.r
    ys = []
    for r in roots:
        if abs(complex(r).imag) < 1e-12 and abs(complex(r).real) < 1:
            y = q / (2 * math.pi) * math.acos(complex(r).real)
            for k in range(0, 6):
                for z in (y + k * q, -y + k * q):
                    if z > 0:
                        ys.append(z)
    ys = sorted(ys)
    print(' %-10s q=%d: numerator degree %d, roots in t %s; first zeros y=%s' % (
        name, q, len(roots), [round(complex(r).real, 6) if abs(complex(r).imag) < 1e-12 else complex(round(r.real, 6), round(r.imag, 6)) for r in roots],
        [round(y, 4) for y in ys[:5]]))
kap = float(((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1))
analyse('chi_-3', 3, {1: 1, 2: -1})
analyse('chi_-4', 4, {1: 1, 3: -1})
analyse('DH', 5, {1: 1, 2: kap, 3: -kap, 4: -1})
analyse('chi_-7', 7, {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1})
analyse('chi_-8', 8, {1: 1, 3: 1, 5: -1, 7: -1})
analyse('chi_-11', 11, {1: 1, 2: -1, 3: 1, 4: 1, 5: 1, 6: -1, 7: -1, 8: -1, 9: 1, 10: -1})
print('--- (d2) check of T_a(u)=sum_{j<q/2} a_j sin(2 pi j/q)/(cosh(2 sqrt(pi u/q))-cos(2 pi j/q)) against the block summed series of a_n n/(u+pi n^2/q)')
def check_T(name, q, a, u, blocks=300000):
    import numpy as np
    n = np.arange(1, q * blocks + 1, dtype=float)
    coef = np.tile(np.array([float(a.get(k % q, 0.0)) for k in range(1, q + 1)]), blocks)
    lhs = float(np.sum(coef * n / (u + np.pi * n**2 / q)))
    rhs = sum(a[j] * math.sin(2 * math.pi * j / q) / (math.cosh(2 * math.sqrt(math.pi * u / q)) - math.cos(2 * math.pi * j / q)) for j in range(1, (q + 1) // 2) if abs(a.get(j, 0)) > 1e-14)
    print(' %-8s u=%g  block sum %.10f  closed form %.10f' % (name, u, lhs, rhs))
check_T('chi_-3', 3, {1: 1, 2: -1}, 1.0)
check_T('chi_-4', 4, {1: 1, 3: -1}, 1.0)
check_T('DH', 5, {1: 1, 2: kap, 3: -kap, 4: -1}, 1.0)
check_T('chi_-7', 7, {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}, 0.3)

print('--- (e) PF minors: det[p((i-j)h)], in balls')
def Phi(t):
    t = abs(t)
    s = arb(0)
    for n in range(1, 8):
        y = pi * n**2 * (2 * t).exp()
        s += (2 * pi**2 * n**4 * (arb(9) / 2 * t).exp() - 3 * pi * n**2 * (arb(5) / 2 * t).exp()) * (-y).exp()
    return 2 * s + arb(0, 1e-80)  # tail n>=8 is below 1e-82 for every t, see E.md
sq5 = arb(5).sqrt()
kappa = ((10 - 2 * sq5).sqrt() - 2) / (sq5 - 1)
def DH(t):
    t = abs(t)
    s = arb(0)
    for n in range(1, 40):
        r = n % 5
        if r == 0:
            continue
        a = [arb(0), arb(1), kappa, -kappa, arb(-1)][r]
        s += a * n * (-pi * n**2 * (2 * t).exp() / 5).exp()
    return (arb(3) / 2 * t).exp() * s + arb(0, 1e-300)  # tail n>=40 is below 1e-400
def minors(f, name, hs, ks):
    for h in hs:
        row = []
        f0 = f(arb(0))
        for k in ks:
            M = arb_mat(k, k, [f(arb(i - j) * arb(h)) / f0 for i in range(k) for j in range(k)])
            row.append('k=%d %s' % (k, M.det().str(3)))
        print(' %s h=%g: ' % (name, h) + '; '.join(row))
minors(Phi, 'Polya', [0.05, 0.1], [2, 3, 4, 5, 6])
minors(DH, 'DH', [0.1, 0.2, 0.3], [2, 3, 4, 5, 6])

print('--- (f) Lebowitz moment test: E[T^{2a}Q^{2b}] <= E[T^{2a}] E[Q^{2b}], T=(X+X\')/sqrt2, Q=(X-X\')/sqrt2')
def moments_from_density(f, T, N, M):
    # trapezoid on [-T,T] with step T/N; double exponentially accurate; the ball carries no quadrature remainder
    h = arb(T) / N
    vals = [f(h * k) for k in range(N + 1)]
    out = []
    for m in range(0, M + 1):
        s = vals[0] if m == 0 else arb(0)
        for k in range(1, N + 1):
            s += 2 * vals[k] * (h * k)**(2 * m)
        out.append(s)
    return [x / out[0] for x in out]
def binom(n, k):
    return math.comb(n, k)
def mixed(mu, a, b):
    # E[(X+Y)^{2a}(X-Y)^{2b}]/2^{a+b}, mu[m]=E X^{2m}
    n = 2 * a + 2 * b
    tot = arb(0)
    for i in range(0, 2 * a + 1):
        for j in range(0, 2 * b + 1):
            ex = i + j          # power of X
            ey = n - ex          # power of Y
            if ex % 2 or ey % 2:
                continue
            sign = (-1)**(2 * b - j)
            tot += sign * binom(2 * a, i) * binom(2 * b, j) * mu[ex // 2] * mu[ey // 2]
    return tot / arb(2)**(a + b)
def lebowitz(mu, name, D=6):
    worst = None
    fails = []
    for a in range(1, D):
        for b in range(1, D - a + 1):
            lhs = mixed(mu, a, b)
            rhs = mixed(mu, a, 0) * mixed(mu, 0, b)
            r = (lhs / rhs)
            if not (r.upper() < 1):
                fails.append((a, b, float(r)))
            v = float(r)
            if worst is None or v > worst[0]:
                worst = (v, a, b)
    print(' %-22s largest ratio lhs/rhs %.6f at (a,b)=(%d,%d); (a,b) with ratio not certified below 1: %s' % (name, worst[0], worst[1], worst[2], fails[:6]))
mu_P = moments_from_density(Phi, 3, 3000, 14)
lebowitz(mu_P, 'Polya')
mu_D = moments_from_density(DH, 4, 4000, 14)
lebowitz(mu_D, 'DH')
# Maxwell law, density x^2 exp(-x^2/2)/sqrt(2 pi): E X^{2m} = (2m+1)!!
dbl = [arb(1)]
for m in range(1, 15):
    dbl.append(dbl[-1] * (2 * m + 1))
lebowitz(dbl, 'Maxwell (LY, not GHS)')
# f_1: L=(1+s^2/g^2) exp(beta s^2), g=gamma_1, beta=Var/2-1/g^2: moments by series
g = 14.134725141734693790
var = float(mu_P[1])
beta = var / 2 - 1 / g**2
c = [Fraction(0)] * 30
def f1_moments(M):
    # coefficients of L in s^2: e_j=beta^j/j! + beta^{j-1}/((j-1)! g^2); E X^{2j}=(2j)! e_j
    out = []
    for j in range(0, M + 1):
        e = arb(beta)**j / math.factorial(j)
        if j >= 1:
            e += arb(beta)**(j - 1) / math.factorial(j - 1) / arb(g)**2
        out.append(e * math.factorial(2 * j))
    return out
lebowitz(f1_moments(14), 'zero completed f_1')
