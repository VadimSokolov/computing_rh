#!/usr/bin/env python3
# Track E. Hankel test of the Lee and Yang property from moments alone, no zero as input. Proposal script for Hopper; small sizes run locally.
# mu_{2m}=E X^{2m} by the trapezoid rule on the theta series, p_m=sum_gamma gamma^{-2m} by Newton, H_n=[p_{i+j+1}]_{0<=i,j<n}.
# If every zero is on the line, H_n is positive definite for every n, and a negative pivot refutes the line. Pivots are the ratios of leading minors.
# The balls carry rounding only: the quadrature remainder (about exp(-pi^2/(2h))) and the series tail are not in the balls, so the output is "computed";
# stability under a change of digits is the check. Run: timeout 170 python3 -I E_hankel.py P 12 300   (law P or DH, size n, decimal digits)
import sys, math
sys.path.insert(0, '/tmp/e_lib')
from flint import arb, ctx
law, n, D = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
ctx.prec = int(3.33 * D) + 64
pi = arb.pi()
M = 2 * n + 1
sq5 = arb(5).sqrt()
kappa = ((10 - 2 * sq5).sqrt() - 2) / (sq5 - 1)
A = [arb(0), arb(1), kappa, -kappa, arb(-1)]
c = pi if law == 'P' else pi / 5
L = 1.2 * D * math.log(10)
tmax = 0.5 * math.log(L / float(c)) + 0.1
nmax = int(math.sqrt(L / float(c))) + 2
h = arb(2.2) / D


def dens(t):
    w = (2 * t).exp()
    base = (-c * w).exp()
    qn, step, sq = base, base * base * base, base * base
    s = arb(0)
    e1 = (arb(9) / 2 * t).exp()
    e2 = (arb(5) / 2 * t).exp()
    for k in range(1, nmax + 1):
        if law == 'P':
            s += (2 * pi**2 * k**4 * e1 - 3 * pi * k**2 * e2) * qn
        elif k % 5:
            s += A[k % 5] * k * qn
        qn *= step
        step *= sq
    return s if law == 'P' else s * (arb(3) / 2 * t).exp()


K = int(tmax * D / 2.2) + 1
mom = [arb(0)] * (M + 1)
for k in range(0, K + 1):
    t = h * k
    f = dens(t) * (1 if k == 0 else 2)
    tp = arb(1)
    t2 = t * t
    for m in range(M + 1):
        mom[m] += f * tp
        tp *= t2
mom = [x / mom[0] for x in mom]
a = [arb(1)]
fact = arb(1)
for m in range(1, M + 1):
    fact *= (2 * m - 1) * (2 * m)
    a.append((-1)**m * mom[m] / fact)
p = [arb(0)] * (M + 1)
for m in range(1, M + 1):
    s = -m * a[m]
    for k in range(1, m):
        s -= p[k] * a[m - k]
    p[m] = s
print(law, 'n=%d digits=%d points=%d terms=%d' % (n, D, K, nmax))
print(' p_1=%s  p_2=%s' % (p[1].str(15), p[2].str(15)))
H = [[p[i + j + 1] for j in range(n)] for i in range(n)]
piv = []
for i in range(n):
    d = H[i][i]
    piv.append(d)
    for r in range(i + 1, n):
        q = H[r][i] / d
        for s2 in range(i, n):
            H[r][s2] = H[r][s2] - q * H[i][s2]
print(' pivots, log10 of the midpoint, sign certified by the ball:')
print(' ' + ' '.join('%.2f%s' % (float(abs(d).log() / arb(10).log()), '+' if d.lower() > 0 else ('-' if d.upper() < 0 else '?')) for d in piv))
print(' last pivot', piv[-1].str(10))
