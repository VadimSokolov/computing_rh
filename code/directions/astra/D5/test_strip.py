"""Bounded D5 experiment, written independently from the book's scripts.

Run: timeout 30 python3 -I test_strip.py
Only floating point arithmetic is used. No numerical certification is claimed.
"""

import json
import platform
import time
from pathlib import Path

import mpmath as mp
import numpy as np


OUT = Path(__file__).resolve().parent


def spin(scale, delta=1):
    scale, delta = mp.mpf(scale), mp.mpf(delta)
    a = mp.cosh(scale * delta)
    return [(-scale, 1 / (2 * (1 + a))),
            (mp.mpf(0), a / (1 + a)),
            (scale, 1 / (2 * (1 + a)))]


def laplace(h, scale, delta=1):
    a = mp.cosh(scale * delta)
    return (a + mp.cosh(scale * h)) / (a + 1)


def partition(h, j, eps, delta=1):
    """Closed form, with both input laws normalized to probability one."""
    a, b = mp.cosh(delta), mp.cosh(eps * delta)
    return (a * b + a * mp.cosh(eps * h) + b * mp.cosh(h)
            + mp.cosh(j * eps) * mp.cosh(h) * mp.cosh(eps * h)
            + mp.sinh(j * eps) * mp.sinh(h) * mp.sinh(eps * h)
            ) / ((a + 1) * (b + 1))


def direct_partition(h, j, eps, delta=1):
    return sum(px * py * mp.exp(j * x * y + h * (x + y))
               for x, px in spin(1, delta)
               for y, py in spin(eps, delta))


def tracked_zero(j, eps, delta=1):
    r = mp.mpc(delta, mp.pi)
    seed = r - eps**2 * (j * r / 2 + j**2 * mp.coth(delta) / 4)
    f = lambda h: partition(h, j, eps, delta)
    return mp.findroot(f, (seed, seed + mp.mpf("0.001") * (1 + 1j)),
                       tol=mp.mpf(10)**(-mp.mp.dps + 15), maxsteps=80)


def identical_width(j, scale=1, delta=1):
    """All zeros follow from a quadratic in cosh(scale * h)."""
    a = mp.cosh(scale * delta)
    u = mp.exp(j * scale**2)
    disc = (u - 1) * ((u + 1) / 2 - a**2)
    cs = [(-a + sign * mp.sqrt(disc)) / u for sign in (-1, 1)]
    return max(abs(mp.re(mp.acosh(c))) / scale for c in cs)


def polynomial_check(n, j):
    """Every zero in one period: t=exp(h/n), degree 2*n+2.

    NumPy supplies initial guesses. Each is refined using the nine-term
    partition sum at 70 digits, independently of the closed-form evaluator.
    """
    eps = mp.mpf(1) / n
    coeff = np.zeros(2 * n + 3)
    sx, sy = spin(1), spin(eps)
    for ix, (x, px) in enumerate(sx):
        for iy, (y, py) in enumerate(sy):
            exponent = n + 1 + n * (ix - 1) + iy - 1
            coeff[exponent] += float(px * py * mp.exp(j * x * y))
    roots = np.roots(coeff[::-1])
    refined = []
    residuals = []
    for t in roots:
        guess = n * mp.log(mp.mpc(complex(t)))
        f = lambda h: direct_partition(h, j, eps)
        z = mp.findroot(f, (guess, guess + mp.mpf("1e-7")),
                        tol=mp.mpf("1e-55"), maxsteps=50)
        z = mp.mpc(mp.re(z), mp.im(z) % (2 * n * mp.pi))
        refined.append(z)
        residuals.append(abs(f(z)))
    # Count distinct roots on the periodic cylinder. These test parameters
    # avoid multiple roots. This is numerical bookkeeping, not a root count proof.
    period = 2 * n * mp.pi
    distinct = []
    for z in refined:
        if all(min(abs(z - w), abs(z - w + 1j * period),
                   abs(z - w - 1j * period)) > mp.mpf("1e-25")
               for w in distinct):
            distinct.append(z)
    assert len(distinct) == 2 * n + 2
    return {"n": n, "degree": 2 * n + 2,
            "distinct_refined_roots": len(distinct),
            "global_width_computed": mp.nstr(max(abs(mp.re(z)) for z in refined), 25),
            "max_direct_residual": mp.nstr(max(residuals), 8)}


def main():
    start = time.monotonic()
    rows = []
    differences = []
    identity_errors = []
    for n in (2, 5, 10, 20, 50, 100, 1000):
        zs = []
        for dps in (60, 100):
            with mp.workdps(dps):
                eps, j = mp.mpf(1) / n, mp.mpf(1)
                z = tracked_zero(j, eps)
                zs.append(z)
                if dps == 100:
                    gap = 1 - mp.re(z)**2
                    rows.append({"J": "1", "delta": "1", "n": n,
                                 "epsilon": mp.nstr(eps, 15),
                                 "real_zero": mp.nstr(mp.re(z), 30),
                                 "imag_zero": mp.nstr(mp.im(z), 30),
                                 "squared_width_gap_at_zero": mp.nstr(gap, 25),
                                 "gap_divided_by_epsilon_squared": mp.nstr(gap / eps**2, 25),
                                 "direct_residual": mp.nstr(abs(direct_partition(z, j, eps)), 8)})
                    probe = mp.mpc("0.73", "2.29")
                    identity_errors.append(abs(partition(probe, j, eps)
                                               - direct_partition(probe, j, eps)))
        with mp.workdps(100):
            differences.append(abs(zs[0] - zs[1]))

    with mp.workdps(70):
        global_checks = [polynomial_check(n, mp.mpf(1)) for n in (2, 5, 10, 20)]
        iid = [{"J": str(j), "width": mp.nstr(identical_width(mp.mpf(str(j))), 25)}
               for j in (0, 0.01, 0.1, 0.5, 1, 1.3, 1.4, 2)]
        scaled_iid = [{"scale": str(s),
                       "closing_J": mp.nstr(mp.log(mp.cosh(2 * mp.mpf(str(s))))
                                             / mp.mpf(str(s))**2, 25)}
                      for s in (1, 0.5, 0.1, 0.01)]
        j_checks = []
        for j in (mp.mpf("0.25"), mp.mpf(1), mp.mpf(4)):
            eps = mp.mpf("0.001")
            z = tracked_zero(j, eps)
            j_checks.append({"J": str(j), "epsilon": "0.001",
                             "gap_over_epsilon_squared": mp.nstr((1-mp.re(z)**2)/eps**2, 25),
                             "asymptotic_coefficient": mp.nstr(j+j**2*mp.coth(1)/2, 25)})
        initial_rates = []
        for j in (mp.mpf("1e-2"), mp.mpf("1e-4"), mp.mpf("1e-6")):
            initial_rates.append({"J": str(j),
                                  "iid_gap_over_J": mp.nstr((1-identical_width(j)**2)/j, 25)})
        rho = mp.mpf("0.1")
        m = mp.sinh(mp.mpf(7)/8) / (2*(1+mp.cosh(1)))
        hmax = mp.sqrt(1+mp.pi**2) + mp.mpf(1)/4
        k = hmax + 1
        eps_explicit = min(mp.mpf("0.5"),
                           mp.sqrt(2*m*rho)/(k*mp.exp((k+mp.mpf(5)/4)/2)))
        summary = {
            "status": "computed, not ball-arithmetic certified",
            "python": platform.python_version(), "mpmath": mp.__version__,
            "numpy": np.__version__, "precision_dps": [60, 100],
            "asymptotic_coefficient_J1": mp.nstr(1+mp.coth(1)/2, 35),
            "max_precision_repeat_difference": mp.nstr(max(differences), 8),
            "max_nine_term_identity_error": mp.nstr(max(identity_errors), 8),
            "counterexample_rows": rows,
            "global_polynomial_checks": global_checks,
            "identical_spin_widths": iid,
            "identical_spin_initial_rates": initial_rates,
            "scaled_identical_spin_closing_thresholds": scaled_iid,
            "coupling_checks": j_checks,
            "explicit_rouche_J1_rho_point1_epsilon": mp.nstr(eps_explicit, 25),
            "wall_seconds": round(time.monotonic()-start, 4),
        }
    (OUT / "numerical_results.json").write_text(json.dumps(summary, indent=2)+"\n")
    lines = ["D5 numerical experiment: computed, not certified",
             "Run: timeout 30 python3 -I test_strip.py",
             "J=delta=1, both input spins have exact defect 1",
             "n epsilon Re(root) Im(root) 1-Re(root)^2 gap/epsilon^2"]
    for row in rows:
        lines.append(" ".join(str(row[k]) for k in
                              ("n", "epsilon", "real_zero", "imag_zero",
                               "squared_width_gap_at_zero", "gap_divided_by_epsilon_squared")))
    lines.extend(["coefficient_J1="+summary["asymptotic_coefficient_J1"],
                  "precision_difference="+summary["max_precision_repeat_difference"],
                  "identity_error="+summary["max_nine_term_identity_error"],
                  "All remaining values: numerical_results.json"])
    (OUT / "numbers.txt").write_text("\n".join(lines)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
