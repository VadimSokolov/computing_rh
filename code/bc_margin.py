# bc_margin.py: the margin b_alpha - 2 N_alpha of the review M2 (misc/reviews/two-volumes/review-barrier-qrh.md; ch/basepoint.tex, formula for x Pi_alpha(dx) and v_alpha)
# in ball arithmetic (python-flint). v_alpha(theta) = b_alpha + int 2 sin^2(x theta / 2) x Pi_alpha(dx) >= b_alpha - 2 N_alpha, because the prime atoms and the positive
# part of the continuous density only add, and 2 sin^2 <= 2 on the negative part, whose mass is
#   N_alpha = int_{x0}^oo e^{-(alpha-1)x} (1 - g(x)) dx,  g(x) = e^{-3x} / (1 - e^{-2x}),  x0 the root of e^{-3x} + e^{-2x} = 1.
# Certified:
#  * x0 = -log y0, y0 the root of y^3 + y^2 = 1 in (0, 1), enclosed by exact rational bisection (signs of an increasing cubic are exact).
#  * N_alpha in closed form: N_alpha = e^{-(alpha-1) x0}/(alpha-1) - sum_{k>=0} e^{-(alpha+2+2k) x0}/(alpha+2+2k); the sum is truncated at K terms with the geometric
#    tail term_K / (1 - e^{-2 x0}) added (the terms decrease with ratio at most e^{-2 x0}).
#  * b_alpha = xi'/xi(alpha) = 1/alpha + 1/(alpha-1) - (1/2) log pi + (1/2) psi(alpha/2) + zeta'/zeta(alpha), psi by Arb (digamma), and
#    -zeta'/zeta(alpha) = sum_n Lambda(n) n^{-alpha} summed to n <= NN with the tail bound sum_{n > NN} ln n n^{-alpha} <= NN^{1-alpha} (ln NN/(alpha-1) + 1/(alpha-1)^2) (Lambda <= ln).
# Cross checks (floating point, mpmath at 40 digits): zeta'/zeta by mpmath, b_alpha by differentiating log xi, N_alpha by quadrature.
# Output: ../data/bc_margin.json.
import json
import os
from fractions import Fraction

from flint import arb, ctx
import mpmath as mp

ctx.prec = 200
NN = 200000


def root_y0():
    lo, hi = Fraction(3, 4), Fraction(4, 5)
    f = lambda y: y ** 3 + y ** 2 - 1
    assert f(lo) < 0 < f(hi)
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def sieve_lambda(n):
    lam = [0.0] * (n + 1)
    isp = bytearray([1]) * (n + 1)
    isp[0:2] = b"\x00\x00"
    for p in range(2, n + 1):
        if isp[p]:
            for m in range(p * p, n + 1, p):
                isp[m] = 0
    return isp


def mzl_sum(alpha, isp):
    """sum_{n <= NN} Lambda(n) n^{-alpha} as a ball (Lambda(p^k) = ln p)."""
    tot = arb(0)
    a = arb(alpha)
    p = 2
    while p <= NN:
        if isp[p]:
            lp = arb(p).log()
            q = p
            while q <= NN:
                tot += lp * (arb(q) ** (-a))
                q *= p
        p += 1
    return tot


def run():
    isp = sieve_lambda(NN)
    yl, yh = root_y0()
    x0 = arb(-1) * arb(yl.numerator) / arb(yl.denominator)
    x0 = (-(arb(yh.numerator) / arb(yh.denominator)).log()).union(-(arb(yl.numerator) / arb(yl.denominator)).log())
    out = {"x0": x0.str(40), "x0_check_root": (arb(-3) * x0).exp().str(10) + " + " + (arb(-2) * x0).exp().str(10)}
    res = {}
    for alpha in (4, 5, 10):
        a = arb(alpha)
        e2 = (-2 * x0).exp()
        K = 400
        s = arb(0)
        for k in range(K):
            c = a + 2 + 2 * k
            s += (-(c) * x0).exp() / c
        last = (-(a + 2 + 2 * K) * x0).exp() / (a + 2 + 2 * K)
        s += last / (1 - e2)
        Nal = (-(a - 1) * x0).exp() / (a - 1) - s
        s_lam = mzl_sum(alpha, isp)
        tail = arb(NN) ** (1 - a) * (arb(NN).log() / (a - 1) + 1 / (a - 1) ** 2)
        zeta_log = -(s_lam + arb(0, tail.upper()))      # zeta'/zeta(alpha), with the tail as radius
        b = 1 / a + 1 / (a - 1) - arb.pi().log() / 2 + (a / 2).digamma() / 2 + zeta_log
        margin = b - 2 * Nal
        # floating point cross checks
        mp.mp.dps = 40
        al = mp.mpf(alpha)
        xi = lambda s_: s_ * (s_ - 1) / 2 * mp.pi ** (-s_ / 2) * mp.gamma(s_ / 2) * mp.zeta(s_)
        b_mp = mp.diff(lambda s_: mp.log(xi(s_)), al)
        x0m = mp.findroot(lambda x: mp.exp(-3 * x) + mp.exp(-2 * x) - 1, 0.28)
        N_mp = mp.quad(lambda x: mp.exp(-(al - 1) * x) * (1 - mp.exp(-3 * x) / (1 - mp.exp(-2 * x))), [x0m, 5, 20, mp.inf])
        res[str(alpha)] = dict(b=b.str(25), N=Nal.str(25), margin=margin.str(25), margin_lower=margin.lower().str(15, radius=False), positive=bool(margin > 0),
                               zeta_log=zeta_log.str(25), b_mpmath=mp.nstr(b_mp, 25), N_mpmath=mp.nstr(N_mp, 25), margin_mpmath=mp.nstr(b_mp - 2 * N_mp, 25))
    out["alpha"] = res
    here = os.path.dirname(os.path.abspath(__file__))
    dd = os.path.join(here, "..", "data")
    json.dump(out, open(os.path.join(dd if os.path.isdir(dd) else here, "bc_margin.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    run()
