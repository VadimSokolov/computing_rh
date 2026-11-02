"""Floating point cross-check of the stated cosine-grid experiment.

This file makes no ball-arithmetic claim. It records trapezoidal, Simpson,
and independent oscillatory quadrature results with their method labels.
Run: python3 grid_recheck.py
"""

from pathlib import Path
import json
import math
import warnings

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad, simpson
from scipy.special import gammaln, zeta

ROOT = Path(__file__).resolve().parent


def logxi_real(s):
    if s == 1:
        return -math.log(2)
    return (math.log(s * (s - 1) / 2) - s / 2 * math.log(math.pi)
            + gammaln(s / 2) + math.log(zeta(s, 1)))


def main():
    u = np.arange(60001, dtype=float) / 1000
    out = {"status": "computed", "numpy": np.__version__, "scipy": scipy.__version__,
           "mpmath": mp.__version__, "step": 0.001, "upper_limit": 60, "rows": {}}
    for a in (1, 1.5, 2, 3):
        s = a + u
        logv = np.empty_like(s)
        good = s != 1
        sg = s[good]
        logv[good] = (np.log(sg * (sg - 1) / 2) - sg / 2 * np.log(np.pi)
                     + gammaln(sg / 2) + np.log(zeta(sg, 1)))
        logv[~good] = -np.log(2)
        phi = np.exp(logxi_real(a) - logv)
        row = {"phi_at_60": float(phi[-1])}
        for x in (50, 100, 200):
            values = phi * np.cos(x * u)
            trapezoidal = np.trapezoid(values, dx=0.001) * x * x
            simpson_value = simpson(values, dx=0.001) * x * x
            with warnings.catch_warnings(record=True) as messages:
                warnings.simplefilter("always")
                q, error = quad(lambda v: math.exp(logxi_real(a) - logxi_real(a + v)),
                                0, np.inf, weight="cos", wvar=x, epsabs=1e-14,
                                limlst=300, limit=300)
            row[str(x)] = {"pi_x2_density_trapezoidal": float(trapezoidal),
                           "pi_x2_density_simpson": float(simpson_value),
                           "pi_x2_density_weighted_quad": float(q * x * x),
                           "quad_internal_error_scaled": float(error * x * x),
                           "quad_warnings": [str(message.message) for message in messages]}
        out["rows"][str(a)] = row
    mp.mp.dps = 70
    identities = []
    for theta, gamma, s in ((1, 2, 3), (1, 0, 2), (mp.mpf("0.25"), 3, 4)):
        theta, gamma, s = map(mp.mpf, (theta, gamma, s))
        density = lambda w: (theta / (theta ** 2 + (w - gamma) ** 2)
                             + theta / (theta ** 2 + (w + gamma) ** 2)) / mp.pi
        integral = mp.quad(lambda w: mp.log(1 + s / w ** 2) * density(w),
                           [0, 1, max(mp.mpf(2), gamma), max(mp.mpf(4), gamma + 1), mp.inf])
        target = mp.log(((theta + mp.sqrt(s)) ** 2 + gamma ** 2) / (theta ** 2 + gamma ** 2))
        identities.append({"theta": str(theta), "gamma": str(gamma), "s": str(s),
                           "integral": mp.nstr(integral, 50),
                           "closed_form": mp.nstr(target, 50),
                           "difference": mp.nstr(integral - target, 10)})
    out["folded_Cauchy_integral_checks"] = identities
    (ROOT / "results" / "grid.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
