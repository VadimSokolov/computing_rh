# ev_evt_lib.py: shared helpers of the ev_*.py scripts that reproduce the numbers of Nick Polson's extreme value hitting times chapter
# (Chapter 25, prefix ev:). Zeros of zeta from data/zeros_1700.txt (68 digits), refinement to higher precision by Newton on xi, xi'(rho)
# at a zero, the kernel Phi of the Riemann Xi function and its folded truncations, Gauss-Legendre moments, power sums of inverse atoms.
# Written from the mathematics (the collaborator's scripts were not read, run or copied). Conventions: xi(s) = s(s-1) pi^(-s/2) Gamma(s/2)
# zeta(s)/2, Xi(t) = xi(1/2 + i t); H(z) = int_0^oo Phi(u) cos(z u) du = Xi(z/2)/8, Phi(u) = sum_n (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u))
# exp(-pi n^2 e^(4u)); folded Phi^F_N keeps n <= N.
import os
from mpmath import mp, mpf, mpc, zeta, gamma, pi, exp, log, loggamma, psi, sqrt, cos, cosh, factorial, matrix, det
from mpmath.calculus.quadrature import GaussLegendre

HERE = os.path.dirname(os.path.abspath(__file__))
ZEROS = os.path.join(HERE, '..', 'data', 'zeros_1700.txt')
if not os.path.exists(ZEROS):
    ZEROS = os.path.join(HERE, 'zeros_1700.txt')  # on the cluster the input file sits beside the script


def zero_strings(n=1700):
    return [l.split()[0] for l in open(ZEROS) if l.strip()][:n]


def xi(s):
    return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def xi_logder(s):
    # xi'/xi at s
    return 1 / s + 1 / (s - 1) - log(pi) / 2 + psi(0, s / 2) / 2 + zeta(s, 1, 1) / zeta(s)


def refine_zero(g, iters=3):
    # Newton on t -> Xi(t) real; Xi'(t) = i xi'(1/2+it)
    t = mpf(g)
    for _ in range(iters):
        s = mpc(mpf(1) / 2, t)
        x = xi(s)
        d = 1j * x * xi_logder(s)
        t = t - (x / d).real
    return t


def xi_prime_at_zero(t):
    rho = mpc(mpf(1) / 2, t)
    return rho * (rho - 1) / 2 * pi ** (-rho / 2) * gamma(rho / 2) * zeta(rho, 1, 1)


def ck_from_zero(t, xi_half):
    # c_k = 2 i xi(1/2) / (gamma_k xi'(1/2 + i gamma_k))
    return (2j * xi_half / (t * xi_prime_at_zero(t))).real


def phi_terms(u, nmax, nmin=1):
    s = mpf(0)
    for n in range(nmin, nmax + 1):
        a = pi * n * n
        s += (2 * a * a * exp(9 * u) - 3 * a * exp(5 * u)) * exp(-a * exp(4 * u))
    return s


def gl_nodes(deg, upper):
    # Gauss-Legendre nodes and weights on [0, upper] (degree deg -> 3*2^(deg-1) nodes), at the current precision
    gl = GaussLegendre(mp)
    nodes = gl.calc_nodes(deg, mp.prec)
    h = mpf(upper) / 2
    return [(h * (x + 1), h * w) for x, w in nodes]


def power_sums(M, kmax):
    # M[m] = int u^(2m) kernel (m = 0..kmax). H(2 i sqrt s)/H(0) = sum_m M_m 4^m s^m / ((2m)! M_0) = prod_k (1 + s/lam_k);
    # returns p_1..p_kmax with p_j = sum_k lam_k^(-j), from log of the series: log = sum_j (-1)^(j-1) p_j s^j / j
    a = [mpf(1)] + [M[m] * mpf(4) ** m / (factorial(2 * m) * M[0]) for m in range(1, kmax + 1)]
    l = [mpf(0)] * (kmax + 1)
    for n in range(1, kmax + 1):
        l[n] = a[n] - sum(mpf(k) / n * l[k] * a[n - k] for k in range(1, n))
    return [None] + [(-1) ** (j - 1) * j * l[j] for j in range(1, kmax + 1)]


def hankels(p, nmax):
    # D_n = det[p_{i+j+1}], D'_n = det[p_{i+j+2}], i,j < n
    D, Dp = [None], [None]
    for n in range(1, nmax + 1):
        D.append(det(matrix(n, n, ) if False else matrix([[p[i + j + 1] for j in range(n)] for i in range(n)])))
        Dp.append(det(matrix([[p[i + j + 2] for j in range(n)] for i in range(n)])))
    return D, Dp


def Fn_inc(z, n):
    # n-th term of the transform of phi_n (half line), analytic in z: (1/2) a^(2-w9) Gamma(w9, a) - (3/4) a^(1-w5) Gamma(w5, a), w_j = (j + i z)/4
    from mpmath import gammainc
    a = pi * n * n
    w9 = (9 + 1j * z) / 4; w5 = (5 + 1j * z) / 4
    return a ** (2 - w9) * gammainc(w9, a) / 2 - mpf(3) / 4 * a ** (1 - w5) * gammainc(w5, a)


def HF(z, n0, n1):
    # sum_{n0 <= n <= n1} (F_n(z) + F_n(-z))/2; n0 = 1, n1 = N is H^F_N(z)
    return sum((Fn_inc(z, n) + Fn_inc(-z, n)) / 2 for n in range(n0, n1 + 1))


def kernel_moments(nodes, lam, nmax, mmax):
    # M[m] = int_0^oo u^(2m) e^(lam u^2) Phi_nmax(u) du with Phi_nmax = sum_{n <= nmax} phi_n (nmax = 40 for the exact kernel at the precision used here)
    vals = [(u, w * exp(lam * u * u) * phi_terms(u, nmax)) for u, w in nodes]
    M = []
    for m in range(mmax + 1):
        M.append(sum(w * u ** (2 * m) for u, w in vals))
    return M
