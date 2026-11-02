"""Bounded synthetic diagnostics for track B. Standard library only.

Run: timeout 30 python3 -I diagnostic.py
These are floating point calculations, not ball arithmetic certificates.
"""

from decimal import Decimal, localcontext
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PI_DEC = Decimal(
    "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
)


def pair_difference(a, delta, x):
    """P_(a-delta)(x) + P_(a+delta)(x) - 2 P_a(x), without cancellation."""
    return (
        2 * a * delta**2 * (a*a - 3*x*x - delta*delta)
        / (math.pi * (a*a + x*x)
           * (x*x + (a-delta)**2) * (x*x + (a+delta)**2))
    )


def p_second(a, x):
    return 2*a*(a*a-3*x*x)/(math.pi*(a*a+x*x)**3)


def quartet_difference(a, delta, height, t):
    return pair_difference(a, delta, t-height) + pair_difference(a, delta, t+height)


def decimal_difference(height, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        a, delta, x = Decimal("0.5"), Decimal("0.1"), Decimal(str(height))
        return +(4*a*delta**2*(a*a-3*x*x-delta*delta)
                 / (PI_DEC*(a*a+x*x)*(x*x+(a-delta)**2)*(x*x+(a+delta)**2)))


def posterior_row(tau, steps):
    a, gamma, delta_max = 0.5, 20.0, 0.25
    kappa = p_second(a, 0.0) + p_second(a, 2*gamma)
    scale = math.sqrt(tau/kappa)
    upper = min(delta_max/scale, 6.0)
    du = upper/steps
    weights = []
    for j in range(steps+1):
        delta = j*du*scale
        shift = pair_difference(a, delta, 0.0) + pair_difference(a, delta, 2*gamma)
        weights.append(math.exp(-0.5*(shift/tau)**2))
    cells = [(weights[j]+weights[j+1])*du/2 for j in range(steps)]
    total = math.fsum(cells)
    target, cumulative, quantile = 0.95*total, 0.0, None
    for j, cell in enumerate(cells):
        if cumulative+cell >= target:
            quantile = (j+(target-cumulative)/cell)*du*scale
            break
        cumulative += cell
    marginal_slab = scale*total/delta_max
    return {
        "tau": tau,
        "q95_slab": quantile,
        "q95_over_sqrt_tau": quantile/math.sqrt(tau),
        "posterior_null_equal_prior_odds": 1/(1+marginal_slab),
        "kappa": kappa,
        "quadrature_cells": steps,
    }


def main():
    a, delta, width, count, tau = 0.5, 0.1, 10.0, 201, 1e-8
    grid = [-width+2*width*j/(count-1) for j in range(count)]
    rows = []
    for height in [20, 40, 80, 160, 1000, 10**13]:
        shifts = [quartet_difference(a, delta, height, t) for t in grid]
        norm = math.sqrt(math.fsum(x*x for x in shifts))
        maximum = max(abs(x) for x in shifts)
        bound = 12*delta*delta/(math.pi*(height-width)**4)
        assert maximum <= bound
        rows.append({
            "height": height,
            "max_abs_shift": maximum,
            "l2_shift": norm,
            "kl_at_tau_1e_minus_8": norm*norm/(2*tau*tau),
            "optimal_sum_of_testing_errors": math.erfc(norm/(2*math.sqrt(2)*tau)),
            "analytic_upper_bound": bound,
            "absolute_fractional_bits_for_one_point_signal": -math.log2(maximum),
        })
    posteriors = []
    max_relative_quad_difference = 0.0
    for error in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]:
        coarse, fine = posterior_row(error, 20000), posterior_row(error, 40000)
        relative = abs(coarse["q95_slab"]-fine["q95_slab"])/fine["q95_slab"]
        max_relative_quad_difference = max(max_relative_quad_difference, relative)
        fine["relative_q95_mesh_difference"] = relative
        posteriors.append(fine)
    decimal60, decimal90 = decimal_difference(10**13, 60), decimal_difference(10**13, 90)
    with localcontext() as ctx:
        ctx.prec = 90
        decimal_relative = abs(decimal60-decimal90)/abs(decimal90)
    global_rows = [
        {"verified_shells": n, "prior_null": p,
         "posterior_null": p/(p+(1-p)/(n+1))}
        for n, p in [(0, 0.5), (10, 0.5), (100, 0.5), (1000, 0.5),
                     (1000000, 0.5), (100, 0.01), (100, 0.99)]
    ]
    out = {
        "status": "synthetic, non-certified, deterministic calculations",
        "settings": {"a": a, "delta": delta, "window_half_width": width,
                     "grid_size": count, "noise_sd": tau},
        "remote_quartet": rows,
        "local_posterior_zero_residual": posteriors,
        "max_relative_q95_mesh_difference": max_relative_quad_difference,
        "decimal_remote_shift_90_digits": str(decimal90),
        "decimal_60_vs_90_relative_difference": str(decimal_relative),
        "first_defect_prior": global_rows,
    }
    (ROOT/"diagnostic_results.json").write_text(json.dumps(out, indent=2)+"\n")
    (ROOT/"results").mkdir(exist_ok=True)
    lines = ["Track B synthetic diagnostics", "Command: timeout 30 python3 -I diagnostic.py",
             "Remote quartet: H, max shift, KL, optimal error sum, signal bits"]
    for row in rows:
        lines.append("{height}: {max_abs_shift:.10e}, {kl_at_tau_1e_minus_8:.10e}, "
                     "{optimal_sum_of_testing_errors:.10e}, "
                     "{absolute_fractional_bits_for_one_point_signal:.6f}".format(**row))
    lines.append("Local slab posterior: tau, q95, q95/sqrt(tau), spike probability")
    for row in posteriors:
        lines.append("{tau:.1e}: {q95_slab:.10e}, {q95_over_sqrt_tau:.10e}, "
                     "{posterior_null_equal_prior_odds:.10e}".format(**row))
    lines.extend([f"Maximum relative q95 mesh change: {max_relative_quad_difference:.10e}",
                  f"Decimal t=0 remote shift: {decimal90}",
                  f"Decimal relative check: {decimal_relative}"])
    for row in global_rows:
        lines.append("First defect prior n={verified_shells}, p={prior_null}: "
                     "P(null|data)={posterior_null:.12g}".format(**row))
    (ROOT/"results"/"numbers.txt").write_text("\n".join(lines)+"\n")
    (ROOT/"empirical.md").write_text(
        "Track B bounded numerical notebook\n\n"
        "Run: `timeout 30 python3 -I diagnostic.py`. All inputs are fixed in the script.\n\n"
        "The remote quartet experiment evaluates 201 rational kernel differences per height. "
        "The local posterior uses one scalar datum with zero residual at a double line zero, "
        "a uniform displacement prior, and trapezoid quadrature on a scaled coordinate. "
        "It compares 20000 and 40000 cells; the numerical upper cutoff is six in that coordinate. "
        "The Gaussian testing errors are calculated exactly from their formula using erfc. "
        "No zeta values, zero lists, network training, or Monte Carlo samples are used. "
        "Binary64 and decimal arithmetic are not ball arithmetic.\n\n"
        "Every reported synthetic value is in `results/numbers.txt` and "
        "`diagnostic_results.json`; formulas and limitations are in `report.md`.\n"
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
