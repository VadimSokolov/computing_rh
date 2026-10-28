# cf_eval.py: shared Arb evaluation of the logarithmic derivatives used by the CF-XI scripts (cf_moments.py, cf_xi.py, cf_info.py, cf_ns.py)
# for the new Chapter 34 (prefix cf). Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, so xi'/xi(s) = 1/s + 1/(s-1) - (log pi)/2 + psi(s/2)/2 + zeta'/zeta(s).
# Davenport and Heilbronn: f(s) = 5^(-s) [zeta(s,1/5) + kappa zeta(s,2/5) - kappa zeta(s,3/5) - zeta(s,4/5)], kappa = (sqrt(10-2 sqrt5)-2)/(sqrt5-1),
# Lambda(s) = (5/pi)^(s/2) Gamma((s+1)/2) f(s), so Lambda'/Lambda = (log 5 - log pi)/2 + psi((s+1)/2)/2 - log 5 + F'/F with F = bracket.
# Density at level theta on the line: u_theta(v) = pi^(-1) Re (L'/L)(theta + i v); the measure of the chapter, U_theta on the line, has this density
# (no factor 1/2; the source's m_sigma is the same function). zeta and its s derivative come from one acb_series jet of length 2, whose
# enclosure is rigorous for ball arguments, so the same code serves point evaluation and acb.integral (analytic = True) calls.
import os
import sys

_p = os.environ.get('CF_FLINT_PATH')   # laptop tests only: a pip --target folder with python-flint (the isolated mode ignores PYTHONPATH)
if _p:
    sys.path.insert(0, _p)
from flint import acb, arb, acb_series, ctx  # noqa: E402,F401


def zjet(s, a=None):
    """(zeta(s), zeta'(s)) or the Hurwitz pair (zeta(s,a), d/ds zeta(s,a)) as acb."""
    ser = acb_series([s, 1], prec=2)
    z = ser.zeta() if a is None else ser.zeta(a)
    return z[0], z[1]


def xi_logder(s):
    z, dz = zjet(s)
    return 1 / s + 1 / (s - 1) - arb.pi().log() / 2 + (s / 2).digamma() / 2 + dz / z


def kappa():
    r5 = arb(5).sqrt()
    return ((10 - 2 * r5).sqrt() - 2) / (r5 - 1)


def dh_logder(s):
    k = kappa()
    F = acb(0)
    dF = acb(0)
    for num, c in ((1, 1), (2, k), (3, -k), (4, -1)):
        z, dz = zjet(s, acb(arb(num) / 5))
        F += c * z
        dF += c * dz
    return (arb(5).log() - arb.pi().log()) / 2 + ((s + 1) / 2).digamma() / 2 - arb(5).log() + dF / F


def density(logder, theta, v):
    """u_theta(v) = pi^(-1) Re logder(theta + i v) as arb; theta, v are arb/acb-compatible."""
    return (logder(acb(theta, v))).real / arb.pi()


def fl(x, nd=12):
    """decimal string of an arb/acb midpoint, nd significant digits"""
    if isinstance(x, acb):
        return fl(x.real, nd) + ('+' if x.imag.mid() >= 0 else '-') + fl(abs(x.imag), nd) + 'i'
    return x.mid().str(nd, radius=False)


def rad(x):
    return float(x.rad()) if isinstance(x, arb) else max(float(x.real.rad()), float(x.imag.rad()))
