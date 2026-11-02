"""The family rho_{*,N}: note's numbers, zeros of f_3, first nonreal zero of rho_{*,3}. Run: timeout 500 python3 -I N6_b.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 256
print('B0 double against ball, N=3')
for t in (20.0, 60.0, 100.3, 150.7 - 0.2j, 500.1 + 0.5j):
    tb = acb(complex(t).real, complex(t).imag)
    v = psiN_ball(tb, 3) * 4 / (gamma_inf_ball(acb(0.5) + I_ * tb) / 2)
    print('  t=%s double %s ball %s' % (t, np.round(Psi_N(np.array([complex(t)]), 3)[0], 9), v.str(10, radius=False)))
print('B1 term size ratio n=2 against n=1')
for lo, hi in ((20, 40), (60, 80), (100, 120), (400, 420)):
    ts = np.arange(lo, hi + 1e-9, 0.05)
    t1 = np.abs(A_term(ts + 0j, 1)) * np.abs(np.exp(np.real(loggamma(0.25 + 0.5j * ts))) ) * (0.25 + ts ** 2)
    t2 = np.abs(A_term(ts + 0j, 2)) * np.abs(np.exp(np.real(loggamma(0.25 + 0.5j * ts))) ) * (0.25 + ts ** 2)
    print('  [%d,%d]: max|T2|/max|T1|=%.3g, max of pointwise ratio=%.3g, mean %.3g' % (lo, hi, t2.max() / t1.max(), (t2 / t1).max(), (t2 / t1).mean()))
print('B2 zeros below 100: N=1,2,3')
tt = np.arange(0.2, 100.0, 0.005)
zs = {}
for N in (1, 2, 3):
    v = np.real(Psi_N(tt + 0j, N))
    v = np.real(np.exp(1j * theta(tt)) * A_N(tt + 0j, N))
    z, _ = sign_zeros(tt, v)
    zs[N] = z
    print('  N=%d count %d, last %s' % (N, len(z), np.round(z[-3:], 3)))
if len(zs[1]) == len(zs[2]):
    d = zs[2] - zs[1]
    print('  max|z2-z1| overall %.4f; below t=50: %.2e; 50 to 70: %.3f; 70 to 100: %.3f' % (np.abs(d).max(), np.abs(d[zs[1] < 50]).max(), np.abs(d[(zs[1] > 50) & (zs[1] < 70)]).max(), np.abs(d[zs[1] > 70]).max()))
if len(zs[1]) == len(zs[3]):
    d = zs[3] - zs[2]
    print('  max|z3-z2| below 100: %.4f; below 70 %.2e' % (np.abs(d).max(), np.abs(d[zs[1] < 70]).max()))
print('B3 counts on 200<t<700 and rectangle')
for N in (1, 2, 3):
    tt = np.arange(200, 700, 0.005)
    v = np.real(np.exp(1j * theta(tt)) * A_N(tt + 0j, N))
    z, _ = sign_zeros(tt, v)
    w, j = winding(lambda t: Psi_N(t, N), rect_path(200, 700, -1, 1, 0.01))
    print('  N=%d real sign changes %d, argument principle %.5f (max step %.3f)' % (N, len(z), w, j))
print('B4 zeros of f_3=1+2^-s+3^-s in sigma>1/2')
fz = lambda s: 1 + 2.0 ** (-s) + 3.0 ** (-s)
dfz = lambda s: -math.log(2) * 2.0 ** (-s) - math.log(3) * 3.0 ** (-s)
tt = np.arange(0, 3000, 0.01)
vals = np.abs(fz(0.62 + 1j * tt))
idx = np.nonzero((vals[1:-1] < vals[:-2]) & (vals[1:-1] < vals[2:]) & (vals[1:-1] < 0.5))[0] + 1
found = []
for i in idx:
    s = 0.62 + 1j * tt[i]
    for _ in range(60):
        s = s - fz(s) / dfz(s)
    if abs(fz(s)) < 1e-12 and s.real > 0.5 and abs(s.imag - tt[i]) < 3:
        if not any(abs(s - f) < 1e-6 for f in found):
            found.append(s)
found = sorted(found, key=lambda s: s.imag)
print('  %d zeros below 3000 with sigma>1/2; first ones (sigma,t): %s' % (len(found), [(round(s.real, 4), round(s.imag, 3)) for s in found[:14]]))
print('  largest sigma found %.4f (limit 0.7878)' % max(s.real for s in found))
print('B5 blocks: real sign changes against rectangle count for N=2,3')
for N in (2, 3):
    line = []
    first = None
    for lo in range(100, 2100, 200):
        hi = lo + 200
        tt = np.arange(lo, hi, 0.004)
        v = np.real(np.exp(1j * theta(tt)) * A_N(tt + 0j, N))
        z, _ = sign_zeros(tt, v)
        w, j = winding(lambda t: Psi_N(t, N), rect_path(lo, hi, -1, 1, 0.01))
        line.append('%d:%d/%.2f' % (lo, len(z), w))
        if first is None and abs(w - len(z)) > 0.5:
            first = lo
    print('  N=%d lo:real/rect %s; first deficit block %s' % (N, ' '.join(line), first))
print('B6 nonreal zeros of rho_{*,3} near zeros of f_3')
def newton(f, t0):
    t = t0
    for _ in range(40):
        h = 1e-5
        fp = (f(np.array([t + h]))[0] - f(np.array([t - h]))[0]) / (2 * h)
        dt = f(np.array([t]))[0] / fp
        t = t - dt
        if abs(dt) < 1e-11:
            return t, True
    return t, False
sols = []
for s0 in found[:40]:
    t0 = s0.imag - 1j * (s0.real - 0.5)
    t, ok = newton(lambda x: Psi_N(x, 3), t0)
    if ok and abs(t.imag) > 1e-6 and abs(t - t0) < 1.5:
        if not any(abs(t - q) < 1e-6 for q, _ in sols):
            sols.append((t, s0))
print('  from the first 40 zeros of f_3: %d converge to nonreal zeros of Psi_3' % len(sols))
for t, s0 in sols[:8]:
    print('   t=%.6f%+.6fi  (sigma=%.4f) from f_3 zero sigma=%.4f, t=%.3f' % (t.real, t.imag, 0.5 - t.imag, s0.real, s0.imag))
if sols:
    tm = min(sols, key=lambda q: q[0].real)[0]
    a, b, c, d = tm.real - 0.02, tm.real + 0.02, tm.imag - 0.02, tm.imag + 0.02
    w, ok, ns = cert_winding(lambda t: psiN_ball(t, 3), a, b, c, d, 0.01)
    print('  certified box around the lowest: [%.5f,%.5f]x[%.5f,%.5f]: winding %s ok %s pieces %d' % (a, b, c, d, w.str(6) if w is not None else 'NA', ok, ns))
