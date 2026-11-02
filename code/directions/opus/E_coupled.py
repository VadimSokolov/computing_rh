#!/usr/bin/env python3
# Track E. Does a ferromagnetic coupling escape the floor TV >= 3.04e-7 of independent sums of uniforms (opus/C Theorem 11)?
# Law: M = X_1 + R, X_1 uniform on [-c_1,c_1], R = sum_{k=2}^7 c_k U_k + sum_{k>7} U_k/(2(k+b)), joint density weighted by exp(J x_1 r), J>=0.
# Every such law is ferromagnetic for J>=0 (Newman, Lieb and Sokal: uniform spins have the Lee and Yang property), with no zero as input.
# Double precision on a grid, no enclosure: "computed". Run: timeout 500 python3 -I E_coupled.py
import numpy as np, math, time, sys
from scipy.optimize import least_squares
from scipy.special import polygamma
N = 2**13
L = 1.6
dx = 2 * L / N
x = (np.arange(N) - N // 2) * dx
theta = 2 * np.pi * np.fft.rfftfreq(N, dx)
xi12 = 0.4971207781883141


def Phi(t):
    t = np.abs(t)
    s = np.zeros_like(t)
    for n in range(1, 8):
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * t) - 3 * np.pi * n**2 * np.exp(2.5 * t)) * np.exp(-np.pi * n**2 * np.exp(2 * t))
    return 2 * s


p = Phi(x) / xi12
print('mass of p on the grid', p.sum() * dx)
b = 1.0214628512344885
K = 3000
tail_log = np.zeros_like(theta)
for k in range(8, K + 1):
    a = theta / (2 * (k + b))
    tail_log += np.log(np.abs(np.sinc(a / np.pi)) + 1e-300)
vrest = float(polygamma(1, K + 1 + b)) / 12
tail_log += -0.5 * vrest * theta**2
gx, gw = np.polynomial.legendre.leggauss(80)
NN = 16


def law(params):
    c = params[:7]
    J = params[7]
    logphi = tail_log.copy()
    for ck in c[1:]:
        logphi += np.log(np.abs(np.sinc(theta * ck / np.pi)) + 1e-300)
    sign = np.ones_like(theta)
    for ck in c[1:]:
        sign *= np.sign(np.sinc(theta * ck / np.pi))
    phi2 = sign * np.exp(np.minimum(logphi, 0.0))
    g = np.fft.fftshift(np.fft.irfft(phi2, N)) / dx
    t = c[0] * gx
    tw = c[0] * gw / (2 * c[0])
    total = np.zeros_like(theta, dtype=complex)
    pref = 1.0
    gn = g.copy()
    tn = np.ones_like(t)
    # the phase convention of rfft with the origin at index N//2 is removed by ifftshift
    E = np.exp(-1j * np.outer(theta, t))
    for n in range(NN + 1):
        F2 = np.fft.rfft(np.fft.ifftshift(gn)) * dx
        F1 = E @ (tn * tw)
        total += pref * F1 * F2
        pref *= J / (n + 1)
        gn = gn * x
        tn = tn * t
    m = np.fft.fftshift(np.fft.irfft(total, N)) / dx
    return m / (m.sum() * dx)


def tv(params):
    return 0.5 * np.abs(law(params) - p).sum() * dx


c0 = np.array([0.2222600795, 0.0866676631, 0.081748211, 0.1494331606, 0.1239904344, 0.0727049767, 0.0523685307])
t0 = time.time()
print('J=0, the independent fit of opus/C: TV=%.4e (reported there 5.13e-7), one evaluation %.3f s' % (tv(np.r_[c0, 0.0]), time.time() - t0))
mask = np.abs(x) < 1.25


def fit(J, c_start, rounds=4):
    c = c_start.copy()
    w = np.ones(mask.sum())
    best = (tv(np.r_[c, J]), c.copy())
    for it in range(rounds):
        def res(z):
            return (law(np.r_[z, J]) - p)[mask] * w
        r = least_squares(res, c, x_scale=np.abs(c0) * 1e-3, max_nfev=40, xtol=1e-12, ftol=1e-12, gtol=1e-14)
        c = r.x
        val = tv(np.r_[c, J])
        if val < best[0]:
            best = (val, c.copy())
        # iteratively reweighted least squares towards the L1 norm
        e = np.abs((law(np.r_[c, J]) - p)[mask])
        w = 1 / np.sqrt(np.maximum(e, 1e-3 * e.max()) / e.max())
    return best


cs = c0.copy()
print('J         TV (best of IRLS)   pi/c_1     pi/c_k')
for J in (0.0, 1e-3, 1e-2, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0):
    val, cs2 = fit(J, cs)
    mm = law(np.r_[cs2, J])
    ths = np.arange(13.0, 22.0, 0.01)
    ph = np.array([(mm * np.cos(th * x)).sum() * dx for th in ths])
    zz = [round(float(ths[i] - ph[i] * 0.01 / (ph[i + 1] - ph[i])), 5) for i in range(len(ths) - 1) if ph[i] * ph[i + 1] < 0]
    print('%-8g  %.4e   %s  zeros of the law in (13,22): %s  elapsed %.0f s' % (J, val, np.round(np.pi / cs2, 4), zz, time.time() - t0), flush=True)
    cs = cs2
    if time.time() - t0 > 240:
        break
