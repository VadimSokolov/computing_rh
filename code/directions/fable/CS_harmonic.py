# Track CS. The random harmonic series X = sum sigma_k / k: density by FFT of prod cos(t/k),
# cumulants, tails against exp(-c e^x), and the best cosh law e^{-a cosh(x/s)} with the same
# variance and fourth cumulant (total variation between them). Floating point, "computed".
import numpy as np
from math import pi, log, exp, gamma
from scipy import integrate, optimize, special

K = 4000                      # explicit factors; tail k > K as a Gaussian with variance sum 1/k^2
T = 6000.0; N = 2 ** 17
t = np.linspace(-T, T, N, endpoint=False)
dt = t[1] - t[0]
logphi = np.zeros(N)
for k0 in range(1, K + 1, 100):
    ks = np.arange(k0, min(k0 + 100, K + 1))
    logphi += np.log(np.abs(np.cos(np.outer(t, 1.0 / ks)))).sum(axis=1)
# sign of the product: count negative cos factors
sgn = np.ones(N)
for k0 in range(1, K + 1, 100):
    ks = np.arange(k0, min(k0 + 100, K + 1))
    sgn *= np.prod(np.sign(np.cos(np.outer(t, 1.0 / ks))), axis=1)
vtail = special.zeta(2) - np.sum(1.0 / np.arange(1, K + 1) ** 2)
phi = sgn * np.exp(logphi - 0.5 * vtail * t ** 2)
# density f(x) = (1/2pi) int phi(t) e^{-itx} dt on x grid
x = np.fft.fftshift(np.fft.fftfreq(N, d=dt)) * 2 * pi
f = np.real(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(phi)))) * dt / (2 * pi)
f = np.maximum(f, 0)
dx = x[1] - x[0]
mask = np.abs(x) < 12
x, f = x[mask], f[mask]
Z = np.trapezoid(f, x); f /= Z
m2 = np.trapezoid(x ** 2 * f, x); m4 = np.trapezoid(x ** 4 * f, x)
k4 = m4 - 3 * m2 ** 2
print("random harmonic series: f(0)=%.6f f(2)=%.9f  Var=%.6f (pi^2/6=%.6f)  k4=%.5f (-2 zeta(4)=%.5f)  excess kurtosis=%.4f"
      % (np.interp(0, x, f), np.interp(2, x, f), m2, pi ** 2 / 6, k4, -2 * special.zeta(4), k4 / m2 ** 2))
# tails
print("tails: x, P(X>x), -log P, e^x, (log 2) e^(x - gamma)")
for xv in [1, 2, 3, 3.5, 4]:
    p = np.trapezoid(f[x > xv], x[x > xv])
    print("  %.1f  %.3e  %.3f  %.3f  %.3f" % (xv, p, -log(max(p, 1e-300)), exp(xv), log(2) * exp(xv - 0.5772156649)))

# cosh law e^{-a cosh y}: moments by quadrature
def cosh_moments(a):
    L = float(np.arccosh(1 + 80.0 / a))          # e^{-a(cosh y - 1)} = e^{-80} at the ends
    w = lambda y: exp(-a * (np.cosh(y) - 1))
    Z = integrate.quad(w, -L, L, points=[0], limit=400)[0]
    m2 = integrate.quad(lambda y: y ** 2 * w(y), -L, L, points=[0], limit=400)[0] / Z
    m4 = integrate.quad(lambda y: y ** 4 * w(y), -L, L, points=[0], limit=400)[0] / Z
    return m2, m4
def excess(a):
    m2c, m4c = cosh_moments(a)
    return (m4c - 3 * m2c ** 2) / m2c ** 2
for a in [0.01, 0.1, 0.5, 1, 2, 5, 20]:
    print("cosh law a=%g: excess kurtosis %.4f" % (a, excess(a)))
target = k4 / m2 ** 2
a_star = optimize.brentq(lambda a: excess(a) - target, 1e-3, 50)
m2c, _ = cosh_moments(a_star)
s = np.sqrt(m2 / m2c)           # X ~ s * Y matches variances
print("matching a = %.6f, scale s = %.6f" % (a_star, s))
Zc = 2 * special.kv(0, a_star)  # exact normaliser of e^{-a cosh y}
g = np.exp(-a_star * np.cosh(x / s)) / (Zc * s)
tv = 0.5 * np.trapezoid(np.abs(f - g), x)
print("total variation between the harmonic law and the matched cosh law: %.5f" % tv)
print("sup |f - g| = %.5f, at x = %.3f" % (np.max(np.abs(f - g)), x[np.argmax(np.abs(f - g))]))
