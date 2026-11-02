#!/usr/bin/env python3
# Track E. Control for E_coupled.py: the same family of ferromagnetic coupled uniforms fitted to the Davenport and Heilbronn law,
# whose transform has an off line zero, and a finer scan for the Polya law. Double precision, "computed".
# Run: timeout 500 python3 -I E_coupled_dh.py
import numpy as np, time
from scipy.optimize import least_squares, brentq
from scipy.special import polygamma, digamma
sq5 = np.sqrt(5.0)
kappa = (np.sqrt(10 - 2 * sq5) - 2) / (sq5 - 1)
A = [0.0, 1.0, kappa, -kappa, -1.0]
gx, gw = np.polynomial.legendre.leggauss(80)
NN = 18


def setup(L, N, dens, b, J0=7, K=3000):
    dx = 2 * L / N
    x = (np.arange(N) - N // 2) * dx
    theta = 2 * np.pi * np.fft.rfftfreq(N, dx)
    p = dens(x)
    p = p / (p.sum() * dx)
    tail = np.zeros_like(theta)
    for k in range(J0 + 1, K + 1):
        tail += np.log(np.abs(np.sinc(theta / (2 * (k + b)) / np.pi)) + 1e-300)
    tail += -0.5 * float(polygamma(1, K + 1 + b)) / 12 * theta**2
    return dict(dx=dx, x=x, theta=theta, p=p, tail=tail, N=N)


def law(S, c, J):
    theta, x, dx, N = S['theta'], S['x'], S['dx'], S['N']
    logphi = S['tail'].copy()
    sign = np.ones_like(theta)
    for ck in c[1:]:
        s = np.sinc(theta * ck / np.pi)
        logphi += np.log(np.abs(s) + 1e-300)
        sign *= np.sign(s)
    g = np.fft.fftshift(np.fft.irfft(sign * np.exp(np.minimum(logphi, 0.0)), N)) / dx
    t = c[0] * gx
    tw = gw / 2
    E = np.exp(-1j * np.outer(theta, t))
    total = np.zeros_like(theta, dtype=complex)
    pref = 1.0
    gn = g.copy()
    tn = np.ones_like(t)
    for n in range(NN + 1):
        total += pref * (E @ (tn * tw)) * (np.fft.rfft(np.fft.ifftshift(gn)) * dx)
        pref *= J / (n + 1)
        gn = gn * x
        tn = tn * t
    m = np.fft.fftshift(np.fft.irfft(total, N)) / dx
    return m / (m.sum() * dx)


def tv(S, c, J):
    return 0.5 * np.abs(law(S, c, J) - S['p']).sum() * S['dx']


def fit(S, J, c_start, rounds=5, nfev=40, cut=None):
    x, p = S['x'], S['p']
    mask = np.abs(x) < cut
    c = c_start.copy()
    w = np.ones(mask.sum())
    best = (tv(S, c, J), c.copy())
    for it in range(rounds):
        r = least_squares(lambda z: (law(S, z, J) - p)[mask] * w, c, x_scale=np.abs(c_start) * 1e-3, max_nfev=nfev, xtol=1e-13, ftol=1e-13, gtol=1e-15)
        c = r.x
        val = tv(S, c, J)
        if val < best[0]:
            best = (val, c.copy())
        e = np.abs((law(S, c, J) - p)[mask])
        w = 1 / np.sqrt(np.maximum(e, 1e-3 * e.max()) / e.max())
    return best


def zeros(S, c, J, a, b, h=0.01):
    m = law(S, c, J)
    ths = np.arange(a, b, h)
    ph = np.array([(m * np.cos(th * S['x'])).sum() * S['dx'] for th in ths])
    return [round(float(ths[i] - ph[i] * h / (ph[i + 1] - ph[i])), 4) for i in range(len(ths) - 1) if ph[i] * ph[i + 1] < 0]


def dens_P(t):
    t = np.abs(t)
    s = np.zeros_like(t)
    for n in range(1, 8):
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * t) - 3 * np.pi * n**2 * np.exp(2.5 * t)) * np.exp(-np.pi * n**2 * np.exp(2 * t))
    return s


def dens_DH(t):
    t = np.abs(t)
    s = np.zeros_like(t)
    for n in range(1, 60):
        a = A[n % 5]
        if a != 0:
            s += a * n * np.exp(np.minimum(-np.pi * n**2 * np.exp(2 * t) / 5, 0.0))
    return np.exp(1.5 * t) * s


t0 = time.time()
print('--- Polya, finer scan near J=4')
S = setup(1.6, 2**13, dens_P, 1.0214628512344885)
cs = np.pi / np.array([16.2688, 44.9588, 35.7307, 21.1272, 25.5363, 35.6944, 45.8734])
for J in (3.0, 3.5, 4.0, 4.5, 5.0):
    val, c2 = fit(S, J, cs, cut=1.25)
    print('J=%-5g TV=%.4e pi/c=%s zeros in (13,34): %s' % (J, val, np.round(np.pi / c2, 3), zeros(S, c2, J, 13, 34)), flush=True)
    if J == 4.0:
        print('   parameters c=', [float(v) for v in c2])
print('elapsed %.0f' % (time.time() - t0))
print('--- Davenport and Heilbronn. Tail constant pi/5 needs psi(1+b)-gamma_E=-log 5')
bdh = brentq(lambda b: digamma(1 + b) - 0.5772156649 + np.log(5.0), -0.9, 3.0)
print(' b=%.5f' % bdh)
S2 = setup(3.2, 2**14, dens_DH, bdh)
print(' variance of the law on the grid %.6f (0.19614 expected)' % (S2['p'] * S2['x']**2).sum().__mul__(S2['dx']))
zs = [5.094159, 8.939914, 12.133545, 14.404003, 17.130239, 19.308, 21.8]
cs = np.pi / np.array(zs)
for J in (0.0, 0.3, 1.0, 2.0, 3.0, 4.0, 6.0):
    val, c2 = fit(S2, J, cs, cut=2.4, nfev=30)
    print('J=%-5g TV=%.4e pi/c=%s zeros in (3,20): %s' % (J, val, np.round(np.pi / c2, 3), zeros(S2, c2, J, 3, 20)), flush=True)
    cs = c2
    if time.time() - t0 > 400:
        break
