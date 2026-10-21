# tf_xi.py: the functions shared by tf_deltac.py and tf_window.py (Chapter tf, Three Families of Generalised Gamma Convolutions and the Deconvolution Obstruction). It is imported, never run.
# Added for the book in October 2026, written from the mathematics (the collaborator's three_families.py was neither read nor run).
# With P(s) = (1/2) s (s - 1) pi^(-s/2) Gamma(s/2) one has xi(s) = P(s) zeta(s), and on the line s = 1/2 + i theta
#   Xi(theta) = xi(s) = -(1/2)(theta^2 + 1/4) pi^(-1/4) |Gamma(1/4 + i theta/2)| Z(theta)      (Z is Hardy's function).
# Put L1 = (log P)' = 1/s + 1/(s - 1) - (1/2) log pi + (1/2) psi(s/2) and L2 = (log P)'' = -1/s^2 - 1/(s - 1)^2 + (1/4) psi'(s/2). Then
#   xi'/xi(s) = L1(s) + zeta'/zeta(s),
#   S(theta) = -(log|Xi|)''(theta) = (xi'/xi)'(s) = L2(s) + zeta''/zeta(s) - (zeta'/zeta(s))^2      (real on the line),
# because Xi'(theta) = i (xi'/xi)(s) Xi(theta), so that (Xi'/Xi)' = -(xi'/xi)'. This is the review's S = -[(log(theta^2 + 1/4))'' - (1/4) Re psi'(1/4 + i theta/2) + (log|Z|)''],
# the first term being -2 Re s^(-2) in the sum -1/s^2 - 1/(s - 1)^2, the second (1/4) Re psi', the third Re (zeta'/zeta)'.
# Near a simple zero rho = 1/2 + i gamma, with u = s - rho, zeta = zeta'(rho) u (1 + a u + b u^2 + ...), a = zeta''/(2 zeta'), b = zeta'''/(6 zeta'), so
# (zeta'/zeta)' = -1/u^2 + 2b - a^2 + O(u), and -1/u^2 = (theta - gamma)^(-2) for u = i (theta - gamma). Hence
#   S_0(gamma) = lim (S(theta) - (theta - gamma)^(-2)) = L2(rho) + zeta'''(rho)/(3 zeta'(rho)) - (zeta''(rho)/zeta'(rho))^2/4,
# which equals (Xi''^2/4 - Xi' Xi'''/3)/Xi'^2 at gamma, the closed form of the review (no numerical averaging, no cancellation).
# All derivatives of zeta are mpmath's zeta(s, derivative=k) at the working precision of the caller.
import os

import mpmath as mp


def half():
    return mp.mpf(1) / 2


def find_zeros_file(explicit=None):
    """data/zeros_1700.txt, wherever it is: --zeros FILE, beside the script (hopper_run.sh copies it there), or in ../data."""
    here = os.path.dirname(os.path.abspath(__file__))
    cands = ([explicit] if explicit else []) + [os.path.join(here, 'zeros_1700.txt'), os.path.join(here, '..', 'data', 'zeros_1700.txt'),
                                                 'zeros_1700.txt', os.path.join('data', 'zeros_1700.txt')]
    for c in cands:
        if os.path.isfile(c):
            return c
    raise FileNotFoundError('zeros_1700.txt not found; pass --zeros FILE')


def read_zero_table(path):
    """The first 1700 ordinates, 68 digits each, as strings (converted at the working precision on use)."""
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def ordinate(n, table):
    """gamma_n, from the table for n <= len(table), else by mpmath's zetazero (working precision)."""
    if n <= len(table):
        return mp.mpf(table[n - 1])
    return mp.zetazero(n).imag


def L1(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.psi(0, s / 2) / 2


def L2(s):
    return -1 / s ** 2 - 1 / (s - 1) ** 2 + mp.psi(1, s / 2) / 4


def xi_logder(s):
    """xi'/xi(s) = L1(s) + zeta'(s)/zeta(s)."""
    return L1(s) + mp.zeta(s, derivative=1) / mp.zeta(s)


def S_line(theta):
    """S(theta) = -(log|Xi|)''(theta), as a complex number whose imaginary part must vanish (it is a check)."""
    s = mp.mpc(half(), theta)
    z0 = mp.zeta(s)
    z1 = mp.zeta(s, derivative=1)
    z2 = mp.zeta(s, derivative=2)
    return L2(s) + z2 / z0 - (z1 / z0) ** 2


def S0_zero(gamma):
    """S_0(gamma) = lim (S(theta) - (theta - gamma)^(-2)) at an on line zero 1/2 + i gamma, in closed form (see the header)."""
    s = mp.mpc(half(), gamma)
    z1 = mp.zeta(s, derivative=1)
    z2 = mp.zeta(s, derivative=2)
    z3 = mp.zeta(s, derivative=3)
    return L2(s) + z3 / (3 * z1) - (z2 / z1) ** 2 / 4


def xi_line(theta):
    """Xi(theta) = xi(1/2 + i theta), real, from the definition of xi (used only to test S by a numerical second derivative)."""
    s = mp.mpc(half(), theta)
    return (s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)).real


def median(xs):
    """The median of a list of numbers; for an even count the mean of the two middle values."""
    ys = sorted(xs)
    n = len(ys)
    return ys[n // 2] if n % 2 else (ys[n // 2 - 1] + ys[n // 2]) / 2


def upper_median(xs):
    ys = sorted(xs)
    return ys[len(ys) // 2]


def fmt(x, n=20):
    return mp.nstr(x, n)
