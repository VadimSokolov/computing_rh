"""Last refinements: largest real zeros of the folded truncations, N=3 against N=2 below 100, DH friable Dirichlet series at sigma=0.55. Run: timeout 200 python3 -I N6_i.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 768
pi = arb.pi()
def X_term(t, n):
    t = acb(t); z = pi * n * n; w = I_ * t / 2
    return 2 * pi * pi * n ** 4 * acb(z).pow(-(acb(arb(9) / 4) + w)) * acb(z).gamma_upper(acb(arb(9) / 4) + w) - 3 * pi * n * n * acb(z).pow(-(acb(arb(5) / 4) + w)) * acb(z).gamma_upper(acb(arb(5) / 4) + w)
def Xi_real(t, N):
    s = acb(0)
    for n in range(1, N + 1):
        s = s + X_term(t, n)
    return float((s.real * arb(PI * t / 4).exp()).mid())
print('I1 largest real zeros of the folded truncations')
for N, g in ((1, 14.045), (2, 39.532), (3, 65.032), (4, 103.368)):
    a, b = g - 0.05, g + 0.05
    print('  N=%d: %.5f' % (N, refine_zero(lambda x: Xi_real(x, N), a, b, 1e-7)))
ctx.prec = 512
print('I2 zeros of N=2 and N=3 below 100')
def pr(t, N):
    s = acb(0)
    for n in range(1, N + 1):
        s = s + T_ball(acb(t), n)
    return float((s.real * arb(PI * t / 4).exp()).mid())
tt = np.arange(30.0, 100.0, 0.01)
z2 = sign_zeros(tt, np.array([pr(float(t), 2) for t in tt]))[0]
z3 = sign_zeros(tt, np.array([pr(float(t), 3) for t in tt]))[0]
d = np.abs(z2 - z3)
print('  counts %d %d; max|z2-z3| on (30,70]: %.5f; (70,90]: %.5f; (90,100): %.5f' % (len(z2), len(z3), d[z2 <= 70].max(), d[(z2 > 70) & (z2 <= 90)].max(), d[z2 > 90].max()))
print('I3 DH friable Dirichlet series f=sum over 3-friable n of a_n n^-s, phases x=2^-it, y=3^-it free')
K_ = KAPPA_DH
a5 = [0.0, 1.0, K_, -K_, -1.0]
fr = []
for i in range(0, 15):
    for j in range(0, 10):
        n = 2 ** i * 3 ** j
        if n <= 10 ** 5:
            fr.append((n, i, j, a5[n % 5]))
fr = [q for q in fr if q[3] != 0]
ns = np.array([q[0] for q in fr], dtype=float); ii = np.array([q[1] for q in fr]); jj = np.array([q[2] for q in fr]); aa = np.array([q[3] for q in fr])
def f(sig, p2, p3):
    return (aa * ns ** (-sig) * np.exp(1j * (ii * p2 + jj * p3))).sum()
for sig in (0.55, 0.65, 0.75):
    best = 9
    g = np.linspace(0, 2 * np.pi, 121)[:-1]
    arg = None
    for p2 in g:
        for p3 in g:
            v = abs(f(sig, p2, p3))
            if v < best:
                best, arg = v, (p2, p3)
    # Newton on the torus
    p = np.array(arg, dtype=float)
    for _ in range(50):
        F0 = f(sig, p[0], p[1])
        h = 1e-6
        F1 = (f(sig, p[0] + h, p[1]) - F0) / h
        F2 = (f(sig, p[0], p[1] + h) - F0) / h
        M = np.array([[F1.real, F2.real], [F1.imag, F2.imag]])
        try:
            dp = np.linalg.solve(M, [F0.real, F0.imag])
        except Exception:
            break
        p = p - dp
    print('  sigma=%.2f grid min |f|=%.4f, after Newton |f|=%.2e at phases (%.4f,%.4f)' % (sig, best, abs(f(sig, p[0], p[1])), p[0] % (2 * np.pi), p[1] % (2 * np.pi)))
print('  zeta_3 has no zero: |zeta_3(s)|>=1/((1+2^-sigma)(1+3^-sigma))')
print('I4 phase criterion thresholds 2 pi exp(2 sum_{p<=P} log p/(sqrt p-1)): ' + ' '.join('P=%d:%.3g' % (P, 2 * PI * math.exp(2 * sum(math.log(p) / (math.sqrt(p) - 1) for p in (2, 3, 5, 7) if p <= P))) for P in (2, 3, 5, 7)))
