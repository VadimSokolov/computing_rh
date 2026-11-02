"""Shared helpers for the N6 review of cosh_spin.tex. Nothing here reads misc/reviews/."""
import sys, os, math
sys.path.insert(0, '/tmp/cs_lib')
import numpy as np
from scipy.special import loggamma
from flint import arb, acb, ctx

PI = math.pi
ROOT = '.'
KAPPA_DH = 0.2840790438404122960282918

def zeta_zeros(n=8000):
    out = []
    with open(os.path.join(ROOT, 'data', 'zeros_8000.txt')) as f:
        for l in f:
            l = l.strip()
            if l:
                out.append(float(l))
            if len(out) >= n:
                break
    return np.array(out)

def S_series(nu, z2, kmax=400):
    """sum_k z2^k/(k!(1-nu)_k), the factor in K_nu = (1/2)Gamma(nu)(z/2)^(-nu) S(nu) + reflected; z2=z^2/4."""
    nu = np.asarray(nu, dtype=complex)
    term = np.ones_like(nu)
    s = np.ones_like(nu)
    for k in range(kmax):
        term = term * z2 / ((k + 1) * (1 - nu + k))
        s = s + term
        if k > 8 and np.max(np.abs(term)) < 1e-18 * max(1.0, np.max(np.abs(s))):
            break
    return s

def A_term(t, n=1, reflect=True):
    """Term n of F_N(t)/(gamma_inf(s)/2), F_N the bracket of Prop. bessel at x=2t; t complex array."""
    t = np.asarray(t, dtype=complex)
    w = 0.25 + 0.5j * t
    z2 = PI ** 2 * n ** 4
    s = 0.5 + 1j * t
    main = ((1 + w) * S_series(2 + w, z2) - 1.5 * S_series(1 + w, z2)) / (w - 0.5)
    if n > 1:
        main = main * np.exp(-s * math.log(n))
    if not reflect:
        return main
    L = math.log(PI * n ** 2)
    lg = loggamma(1 + w) + 0j
    e1 = np.exp(loggamma(-2 - w) + (2 + w) * L + w * math.log(PI) - lg)
    e2 = np.exp(loggamma(-1 - w) + (1 + w) * L + w * math.log(PI) - lg)
    refl = (PI ** 2 * n ** 4 * e1 * S_series(-2 - w, z2) - 1.5 * PI * n ** 2 * e2 * S_series(-1 - w, z2)) / (w - 0.5)
    return main + refl

def A_N(t, N=1):
    t = np.asarray(t, dtype=complex)
    out = np.zeros_like(t)
    for n in range(1, N + 1):
        out = out + A_term(t, n)
    return out

def theta(t):
    """Riemann Siegel theta, continuous, real t (or complex: analytic continuation of the same expression)."""
    t = np.asarray(t)
    return np.imag(loggamma(0.25 + 0.5j * t)) - 0.5 * t * math.log(PI)

def chi(t):
    t = np.asarray(t, dtype=complex)
    return np.exp(1j * t * math.log(PI) + loggamma(0.25 - 0.5j * t) - loggamma(0.25 + 0.5j * t))

def Psi_N(t, N=1):
    """4 psi_N(2t)/((1/2)gamma_inf(1/2+it)) = A_N(t)+chi(t)A_N(-t), for complex t. Same zeros as psi_N(2t) away from t=+-i/2."""
    t = np.asarray(t, dtype=complex)
    return A_N(t, N) + chi(t) * A_N(-t, N)

def g_real(t, N=1):
    t = np.asarray(t, dtype=float)
    return np.real(np.exp(1j * theta(t)) * A_N(t + 0j, N))

def main_term(T):
    return T / (2 * PI) * np.log(T / (2 * PI * math.e)) + 7 / 8

def winding(f, path_pts):
    """(1/2pi) total change of arg f along a closed polygon given by its vertices; f takes complex arrays."""
    z = np.asarray(path_pts, dtype=complex)
    v = f(z)
    d = np.angle(v[1:] / v[:-1])
    return d.sum() / (2 * PI), np.max(np.abs(d))

def rect_path(a, b, c, d, h):
    """counterclockwise rectangle a<=Re<=b, c<=Im<=d, step about h."""
    def seg(p, q):
        n = max(2, int(abs(q - p) / h) + 1)
        return p + (q - p) * np.linspace(0, 1, n, endpoint=False)
    return np.concatenate([seg(a + 1j * c, b + 1j * c), seg(b + 1j * c, b + 1j * d), seg(b + 1j * d, a + 1j * d), seg(a + 1j * d, a + 1j * c), np.array([a + 1j * c])])

def sign_zeros(ts, vals):
    idx = np.nonzero(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    return ts[idx] - vals[idx] * (ts[idx + 1] - ts[idx]) / (vals[idx + 1] - vals[idx]), idx

def refine_zero(f, a, b, tol=1e-11):
    fa, fb = f(a), f(b)
    for _ in range(100):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
        if b - a < tol:
            break
    return 0.5 * (a + b)

# ball arithmetic, Arb trusted
I_ = acb(0, 1)

def K_ball(z, nu):
    return acb(z).bessel_k(nu)

def F_ball(t, N=1, b1=None, b2=None):
    """bracket 2pi^2 n^4 K_{9/4+it/2}(2pi n^2) - 3 pi n^2 K_{5/4+it/2}(2 pi n^2), summed over n<=N; t acb."""
    pi = arb.pi()
    b1 = arb(9) / 4 if b1 is None else b1
    b2 = arb(5) / 4 if b2 is None else b2
    t = acb(t)
    out = acb(0)
    for n in range(1, N + 1):
        z = 2 * pi * n * n
        out = out + 2 * pi * pi * n ** 4 * K_ball(z, acb(b1) + I_ * t / 2) - 3 * pi * n * n * K_ball(z, acb(b2) + I_ * t / 2)
    return out

def psi_ball(t, N=1):
    """psi_N(2t) as an entire function of t."""
    t = acb(t)
    return (F_ball(t, N) + F_ball(-t, N)) / 4

def gamma_inf_ball(s):
    s = acb(s)
    pi = arb.pi()
    return s * (s - 1) / 2 * acb(pi).pow(-s / 2) * (s / 2).gamma()

# certified winding number by ball arithmetic: each boundary piece is enclosed in one ball whose image avoids a sector
def _seg_ball(p, q):
    mx, my = (p.real + q.real) / 2, (p.imag + q.imag) / 2
    rx, ry = abs(p.real - q.real) / 2, abs(p.imag - q.imag) / 2
    return acb(arb(mx, rx * (1 + 1e-9) + 1e-30), arb(my, ry * (1 + 1e-9) + 1e-30))

def cert_winding(fun, a, b, c, d, h=0.1, maxdepth=14):
    """Winding number of fun on the boundary of [a,b]x[c,d] (variable t), certified. Returns (ball, ok, nseg)."""
    total = arb(0)
    nseg = 0
    pi2 = 2 * arb.pi()
    def seg(p, q, depth):
        nonlocal total, nseg
        B = fun(_seg_ball(p, q))
        m = acb(B.real.mid(), B.imag.mid())
        rr = B.real.rad() + B.imag.rad()
        if (abs(m) * arb(7) / 10 - rr) > 0:
            vp = fun(acb(p.real, p.imag))
            vq = fun(acb(q.real, q.imag))
            total = total + (vq / m).arg() - (vp / m).arg()
            nseg += 1
            return True
        if depth >= maxdepth:
            return False
        r = (p + q) / 2
        return seg(p, r, depth + 1) and seg(r, q, depth + 1)
    cs = [complex(a, c), complex(b, c), complex(b, d), complex(a, d), complex(a, c)]
    for u, v in zip(cs[:-1], cs[1:]):
        n = max(1, int(abs(v - u) / h))
        for j in range(n):
            if not seg(u + (v - u) * j / n, u + (v - u) * (j + 1) / n, 0):
                return None, False, nseg
    return total / pi2, True, nseg

def T_ball(t, n):
    pi = arb.pi()
    z = 2 * pi * n * n
    return 2 * pi * pi * n ** 4 * K_ball(z, acb(arb(9) / 4) + I_ * t / 2) - 3 * pi * n * n * K_ball(z, acb(arb(5) / 4) + I_ * t / 2)

def psiN_ball(t, N):
    t = acb(t)
    return sum((T_ball(t, n) + T_ball(-t, n) for n in range(2, N + 1)), T_ball(t, 1) + T_ball(-t, 1)) / 4
