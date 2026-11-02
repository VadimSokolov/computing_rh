"""Independent audit of note.tex, Section 4. Run: python3 smallball_audit.py.

The saddle and Lambert values use Arb balls. Inverse Laplace values are
mpmath computations, with independent method and degree comparisons.
"""
from pathlib import Path
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import mpmath as mp
from flint import arb, arb_series, ctx

mp.mp.dps = 90
ctx.dps = 100
BASE = Path(__file__).resolve().parent


def logxi(z):
    if z == 1 or z == 0:
        return -mp.log(2)
    return (mp.log(z) + mp.log(z - 1) - mp.log(2)
            - z * mp.log(mp.pi) / 2 + mp.loggamma(z / 2)
            + mp.log(mp.zeta(z)))


def logxi_derivative(z):
    return (1 / z + 1 / (z - 1) - mp.log(mp.pi) / 2
            + mp.digamma(z / 2) / 2 + mp.diff(mp.zeta, z) / mp.zeta(z))


def blogxi_series(z):
    v = arb_series([z, 1], prec=3)
    return (v.log() + (v - 1).log() - arb(2).log()
            - v * arb.pi().log() / 2 + (v / 2).lgamma()
            + v.zeta().log())


def blogxi_at_alpha(alpha):
    if alpha == 1:
        return -arb(2).log()
    return blogxi_series(arb(alpha))[0]


def saddle(alpha, estr):
    eps = mp.mpf(estr)
    x0 = -mp.lambertw(-8 * mp.pi * eps, -1).real / (4 * eps)
    xhat = mp.findroot(lambda x: logxi_derivative(alpha + x) - 2 * eps * x,
                       (x0, x0 + 2), tol=mp.mpf("1e-75"))
    radius = mp.mpf("1e-55") * xhat
    lo = arb(mp.nstr(xhat - radius, 85))
    hi = arb(mp.nstr(xhat + radius, 85))
    E = arb(estr)
    # Concavity in s follows from the GGC proof in report.md.
    assert blogxi_series(lo + alpha)[1] - 2 * E * lo > 0
    assert blogxi_series(hi + alpha)[1] - 2 * E * hi < 0
    X = arb(mp.nstr(xhat, 85), mp.nstr(2 * radius, 85))
    v = blogxi_series(X + alpha)
    I = v[0] - blogxi_at_alpha(alpha) - E * X**2
    var = (v[1] - 2 * X * v[2]) / (4 * X**3)
    assert var > 0
    bx0 = -(-8 * arb.pi() * E).lambertw(-1) / (4 * E)
    lambert = E * bx0**2 - bx0 / 2 + arb(alpha + 3) * bx0.log() / 2
    c = arb.pi().log() / 2 - arb(alpha) * (2 * arb.pi()).log() / 2 - blogxi_at_alpha(alpha)
    leading = (1 / E).log()**2 / (16 * E)
    pref = (-I).exp() / (X**2 * (2 * arb.pi() * var).sqrt())
    return {
        "alpha": alpha, "epsilon": estr,
        "x_bracket_lower": str(lo), "x_bracket_upper": str(hi),
        "I_ball": str(I), "Lambert_ball": str(lambert),
        "I_minus_Lambert_ball": str(I - lambert),
        "c_alpha_ball": str(c), "leading_ball": str(leading),
        "I_over_leading_ball": str(I / leading),
        "saddle_probability_ball": str(pref),
        "variance_under_tilt_ball": str(var),
        "I_mid": str(I.mid()), "saddle_probability_mid": str(pref.mid()),
    }


def invert(alpha, estr, method, degree):
    eps = mp.mpf(estr)
    started = time.monotonic()
    value = mp.invertlaplace(
        lambda s: mp.exp(logxi(alpha) - logxi(alpha + mp.sqrt(s))) / s,
        eps, method=method, degree=degree)
    return {"alpha": alpha, "epsilon": estr, "method": method,
            "degree": degree, "probability": mp.nstr(value, 65),
            "seconds": round(time.monotonic() - started, 3)}


def main():
    rows = [saddle(a, e) for a in (1, 2)
            for e in ("0.004", "0.003", "0.001", "0.0001", "0.000001")]
    inversions = []
    data = {"status_saddles": "computed", "arithmetic_saddles": "Arb balls, 100 decimal digits",
            "status_inversions": "computed", "rows": rows, "inversions": inversions}
    (BASE / "smallball_results.json").write_text(json.dumps(data, indent=2) + "\n")
    for e in ("0.004", "0.003"):
        for method, degree in (("talbot", 160), ("talbot", 240), ("dehoog", 160)):
            rec = invert(1, e, method, degree)
            inversions.append(rec)
            print(json.dumps(rec), flush=True)
            (BASE / "smallball_results.json").write_text(json.dumps(data, indent=2) + "\n")
    print("Wrote smallball_results.json", flush=True)


if __name__ == "__main__":
    main()
