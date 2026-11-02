"""Arb Bromwich inversion with an explicit analytic tail bound.

The probability enclosure includes quadrature and the infinite tail.
Run: python3 ball_inversion.py. See report.md for the tail proof.
"""
from pathlib import Path
import json
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
from flint import arb, acb, ctx

ctx.dps = 65
BASE = Path(__file__).resolve().parent


def xi(z):
    return z * (z - 1) / 2 * (-z * arb.pi().log() / 2).exp() * (z / 2).gamma() * z.zeta()


def compute(estr, xstr):
    start = time.monotonic()
    eps = arb(estr)
    # The contour can use any positive x. These exact decimal values do not
    # require an enclosed saddle, since contour inversion is an identity.
    x = arb(xstr)
    s0 = x*x
    z0 = 1 + x
    xi0 = xi(acb(z0)).real
    I = xi0.log() + arb(2).log() - s0*eps
    # For t >= T=8*u0^2, Re(sqrt(s0+it))/2 >= u0.
    u0 = arb(300)
    T = 8*u0*u0
    vmax = T/s0
    c = (u0 / (arb.pi() * arb("1.5").exp())).log()
    C = 2 * arb("0.5").exp() / (2*arb.pi()).sqrt()
    tail_integral = xi0 * C * u0**arb("-1.5") * 2 * (-c*u0).exp() / (c*u0)
    # Normalized integrand: the final prefactor is exp(-I)/pi.
    def integrand(v, analytic):
        q = 1 + acb(0, 1)*v
        z = 1 + x*q.sqrt(analytic=analytic)
        return (acb(0, 1)*s0*eps*v).exp() * xi0 / (q*xi(z))

    total = acb(0)
    edges = [arb(0), arb("0.125"), arb("0.25"), arb("0.5"), arb(1), arb(2), arb(4), arb(8)]
    while edges[-1] < vmax:
        edges.append(edges[-1]*2)
    edges[-1] = vmax
    for lo, hi in zip(edges, edges[1:]):
        val = acb.integral(integrand, acb(lo), acb(hi),
                           abs_tol=arb("1e-40"), rel_tol=arb("1e-40"),
                           eval_limit=300000, depth_limit=40)
        if not val.is_finite():
            raise ArithmeticError("Quadrature returned a nonfinite ball")
        total += val
        print(estr, "through", hi.str(5), "integral", total.real.str(12), flush=True)
    enclosed = total.real + arb(0, tail_integral.upper())
    probability = (-I).exp() / arb.pi() * enclosed
    assert probability > 0
    return {"status": "computed", "arithmetic": "Arb ball quadrature with analytic infinite-tail bound",
            "epsilon": estr, "x_contour": xstr, "T": str(T),
            "normalized_tail_bound": str(tail_integral),
            "normalized_integral_ball": str(enclosed),
            "probability_ball": str(probability),
            "seconds": round(time.monotonic()-start, 3)}


if __name__ == "__main__":
    results=[]
    for eps, x in (("0.004", "220"), ("0.003", "350")):
        results.append(compute(eps, x))
        (BASE / "ball_inversion_results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2), flush=True)
