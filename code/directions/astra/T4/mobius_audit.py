#!/usr/bin/env python3
"""Independent reconstruction of note.tex, lines 134-201.

Arb balls enclose the scalar constants and zeros. Table quadratures and
optimizations use double precision, with independent resolution runs; these
table entries are computed, not certified. No author code or zero data used.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "pylib"))
from flint import acb, arb, arb_series, ctx


def mobius_sieve(limit):
    mu = [1] * (limit + 1)
    for p in range(2, limit + 1):
        if all(p % d for d in range(2, int(p**0.5) + 1)):
            for j in range(p, limit + 1, p):
                mu[j] = -mu[j]
            for j in range(p * p, limit + 1, p * p):
                mu[j] = 0
    mu[0] = 0
    return mu


def xi_ball(s):
    s = acb(s)
    return s * (s - 1) * (-s * arb.pi().log() / 2).exp() * (s / 2).gamma() * s.zeta() / 2


def xi_mp(s):
    return s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def scalar_constants():
    ctx.dps = 60
    ctx.cap = 5
    s = arb_series([arb("0.5"), 1])
    xi = s * (s - 1) * (-s * arb.pi().log() / 2).exp() * (s / 2).gamma() * s.zeta() / 2
    logxi = xi.log()
    variance = 2 * logxi[2]
    total_s4 = -2 * logxi[4]
    zeros_ball = [acb.zeta_zero(k).imag for k in range(1, 30)]
    s2 = sum(g**-2 for g in zeros_ball)
    s4 = sum(g**-4 for g in zeros_ball)
    sigma = variance - 2 * s2
    algebraic_tail = 2 * (total_s4 - s4)
    absolute_tail_bound = algebraic_tail * arb(1600080001) / 1599760001
    g1, g2 = zeros_ball[:2]
    counter_sigma = variance - 2/g1**2
    counter_psi = (-counter_sigma*g2**2/2).exp()
    for k, power in [(1, 1), (2, -1), (3, -1)]:
        y = arb.pi()*g2/(k*g1)
        counter_psi *= (y.sin()/y)**power
    g30 = acb.zeta_zero(30).imag
    k1_omitted_zero_psi = (-sigma*g30**2/2).exp()
    for g in zeros_ball:
        y = arb.pi()*g30/g
        k1_omitted_zero_psi *= y.sin()/y
    f0_partial = sum((4 * arb.pi()**2 * n**4 - 6 * arb.pi() * n**2) * (-arb.pi() * n*n).exp() for n in range(1, 10)) / xi[0]
    # For n>=10: 4*pi^2*n^4*exp(-pi*n^2), the ratio of successive terms is
    # less than exp(-60). Add a symmetric 1e-130 ball to cover this tail.
    f0 = f0_partial + arb("0 +/- 1e-130")
    p2_at_pi = 8*arb.pi()**3 - 56*arb.pi()**2 + arb("82.5")*arb.pi() - arb("18.75")
    density_n1_second_derivative_upper = 2*arb.pi()*arb("-3.2").exp()*p2_at_pi
    density_other_second_derivative_upper = arb(".005").exp()*sum(
        (16*(arb.pi()*n*n)**4 + 165*(arb.pi()*n*n)**2)*(-arb.pi()*n*n).exp() for n in range(2, 10)
    ) + arb("0 +/- 1e-120")
    # A global elementary envelope for the omitted Fourier tail uses
    # |zeta(1/2+it)| <= 2t+2 and |Gamma(1/4+it/2)| <= sqrt(2)*Gamma(1/4)*exp(-pi*t/4).
    a, cut = arb.pi()/4, arb(160)
    def exponential_moment(degree):
        return (-a*cut).exp()*sum(arb(math.factorial(degree))/math.factorial(degree-j)*cut**(degree-j)/a**(j+1) for j in range(degree+1))
    fourier_far_tail = (arb(2).sqrt()*arb(".25").gamma()*arb.pi()**arb("-.25")/(arb.pi()*xi[0])) * (exponential_moment(3)+exponential_moment(2)+(exponential_moment(1)+exponential_moment(0))/4)
    mu = mobius_sieve(10000)
    tails = {k: 6 / arb.pi()**2 - sum(arb(mu[j]) / j**2 for j in range(1, k + 1)) for k in [1, 2, 3, 5, 10, 30]}
    mp.mp.dps = 65
    t = mp.mpf(100)
    integral_bounds = {}
    for name, a, b, c in [("arxiv", ".111", ".275", "2.450"), ("published", ".112", ".278", "2.510")]:
        a, b, c = map(mp.mpf, (a, b, c))
        one_no_boundary = (2 / (3 * mp.pi * t**3)) * (mp.log(t/(2*mp.pi)) - mp.mpf(2)/3)
        one_no_boundary += (mp.mpf(7)/8 + c + mp.mpf(".2")/t + a*(mp.log(t)+mp.mpf(1)/4) + b*mp.log(mp.log(t))) / t**4
        one_no_boundary += b * mp.e1(4 * mp.log(t))
        integral_bounds[name] = {
            "one_sided_no_boundary": str(one_no_boundary),
            "two_sided_no_boundary": str(2*one_no_boundary),
            "one_sided_with_N100_boundary": str(one_no_boundary - 29/t**4),
            "two_sided_with_N100_boundary": str(2*one_no_boundary - 58/t**4),
        }
    constants = {
        "arithmetic": "certified scalar enclosures from python-flint Arb at 60 decimal digits; counting integral evaluations use mpmath at 65 digits",
        "xi_half": str(xi[0]), "variance": str(variance), "S_T": str(s2),
        "S4_T": str(s4), "sigma_T_squared": str(sigma),
        "algebraic_one_sided_total_S4": str(total_s4),
        "algebraic_two_sided_tail_S4": str(algebraic_tail),
        "algebraic_tail_note": "This equals the absolute tail only conditional on RH above T. The algebraic identity itself is unconditional.",
        "unconditional_absolute_tail_upper_ball": str(absolute_tail_bound),
        "absolute_tail_factor": "1600080001/1599760001, from the strip and T=100",
        "lemma_counterexample_T_gamma1_K3_t_gamma2_psi": str(counter_psi),
        "K1_T100_gaussian_counterexample_at_gamma30": str(k1_omitted_zero_psi),
        "N100": str(arb(100).zeta_nzeros()), "f_X_zero": str(f0),
        "density_monotonicity_n1_second_derivative_upper": str(density_n1_second_derivative_upper),
        "density_monotonicity_remaining_second_derivative_upper": str(density_other_second_derivative_upper),
        "fourier_integral_tail_beyond_160_upper_ball": str(fourier_far_tail),
        "gaussian_smoothing_variance_K1": str(arb.pi()**2 / 3 * s2 * abs(tails[1])),
        "mobius_tails": {str(k): str(v) for k, v in tails.items()},
        "zero_count_integrals": integral_bounds,
        "zeros_balls": [str(g) for g in zeros_ball],
    }
    return constants, np.array([float(g.mid()) for g in zeros_ball]), float(xi[0].mid()), float(variance.mid()), float(s2.mid()), float(s4.mid()), float(sigma.mid()), {k: float(v.mid()) for k, v in tails.items()}, mu


def density_kernel(x, xi_half):
    x = np.abs(np.asarray(x))
    answer = np.zeros_like(x, dtype=float)
    q = np.pi * np.exp(2*x)
    for n in range(1, 10):
        answer += (4*np.pi**2*n**4*np.exp(4.5*x) - 6*np.pi*n*n*np.exp(2.5*x)) * np.exp(-q*n*n)
    return answer / xi_half


def quadrature(a, b, split_zeros, degree):
    breaks = [a] + [x for x in split_zeros if a < x < b] + [b]
    q, w = leggauss(degree)
    nodes, weights = [], []
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        # Split at most 3 units apart, even when there are no zeros nearby.
        fine = np.linspace(lo, hi, int(np.ceil((hi-lo)/3)) + 1)
        for left, right in zip(fine[:-1], fine[1:]):
            nodes.append((left+right)/2 + (right-left)*q/2)
            weights.append((right-left)*w/2)
    return np.concatenate(nodes), np.concatenate(weights)


def phi_values(nodes, xi_half):
    ctx.dps = 30
    return np.array([float((xi_ball(acb(0.5, float(t))) / xi_half).real.mid()) for t in nodes])


def psi_values(nodes, zeros, sigma, k, mu):
    logabs = -sigma*nodes*nodes/2
    sign = np.ones_like(nodes)
    for g in zeros:
        for j in range(1, k+1):
            if mu[j]:
                factor = np.sinc(nodes/(j*g))
                sign *= np.sign(factor)
                logabs += mu[j]*np.log(np.abs(factor))
    return sign*np.exp(logabs)


def table_run(degree, data):
    constants, zeros, xi_half, variance, s2, s4, sigma, tails, mu = data
    allzeros = [float(acb.zeta_zero(k).imag.mid()) for k in range(1, 59)]
    rows = []
    for k in [1, 2, 3, 5, 10, 30]:
        u = min((k+1)*zeros[0]/2, 50.0)
        t, w = quadrature(0.0, u, allzeros, degree)
        phi = phi_values(t, xi_half)
        psi = psi_values(t, zeros, sigma, k, mu)
        tfar, wfar = quadrature(u, 160.0, allzeros, degree)
        phifar = phi_values(tfar, xi_half)
        fourier_tail = float(np.dot(wfar, abs(phifar))/np.pi)
        a4 = (2*np.pi**4/90)/(9*k**3)*s4
        a2 = np.pi**2/6 * abs(tails[k]) * s2
        crude_a2 = np.pi**2/6 / k * s2
        b4 = (2.0/3)*4.87e-7
        refined = float(np.dot(w, abs(phi)*np.expm1(a2*t*t+(a4+b4)*t**4))/np.pi + fourier_tail)
        crude = float(np.dot(w, abs(phi)*np.expm1(crude_a2*t*t+(a4+b4)*t**4))/np.pi + fourier_tail)
        alternative_crude = float(np.dot(w, abs(phi)*np.expm1((4/3)*crude_a2*t*t+b4*t**4))/np.pi + fourier_tail)
        corrected_tail = float(constants["zero_count_integrals"]["published"]["two_sided_with_N100_boundary"])
        corrected_refined = float(np.dot(w, abs(phi)*np.expm1(a2*t*t+(a4+corrected_tail/3)*t**4))/np.pi + fourier_tail)
        strip_tail = float(arb(constants["unconditional_absolute_tail_upper_ball"]).upper())
        strip_refined = float(np.dot(w, abs(phi)*np.expm1(a2*t*t+(a4+strip_tail/3)*t**4))/np.pi + fourier_tail)
        # The tighter coefficient is 1/3 for a two-sided tail, since the product
        # itself uses one representative from each plus/minus pair.
        weighted_psi = w*psi/np.pi
        def signed_error(x):
            return float(np.dot(weighted_psi, np.cos(t*x)) - density_kernel(x, xi_half))
        grid = np.linspace(0, 3, 3001)
        errors = np.array([signed_error(x) for x in grid])
        locs = [0, len(grid)-1]
        locs += [i for i in range(1, len(grid)-1) if abs(errors[i]) >= abs(errors[i-1]) and abs(errors[i]) >= abs(errors[i+1])]
        candidates = [(abs(errors[i]), grid[i]) for i in [0, len(grid)-1]]
        for i in locs[2:]:
            opt = minimize_scalar(lambda x: -abs(signed_error(x)), bounds=(grid[i-1], grid[i+1]), method="bounded", options={"xatol":1e-13})
            candidates.append((-float(opt.fun), float(opt.x)))
        actual_error, argmax = max(candidates)
        pole = next((n for n in range(2, 10001) if sum(mu[d] for d in range(1,k+1) if n%d == 0) < 0), None)
        rows.append({"K": k, "U": u, "m_K": tails[k], "actual_error_on_abs_x_le_3": actual_error, "argmax_abs_x": argmax,
                     "signed_error_at_zero": signed_error(0), "refined_bound_as_printed": refined, "crude_bound_as_printed": crude,
                     "crude_from_full_log_bound_with_4_over_3": alternative_crude,
                     "refined_bound_corrected_tighter_tail": corrected_refined,
                     "refined_bound_strip_derivative_tail": strip_refined,
                     "fourier_tail_U_to_160": fourier_tail,
                     "first_pole_multiple_up_to_10000": pole,
                     "refined_over_actual": refined/actual_error,
                     "percent_of_f0": 100*actual_error/float(density_kernel(0, xi_half)),
                     "coefficient_t2": a2, "coefficient_t4": a4+b4})
        print("computed", degree, k, actual_error, refined, crude, flush=True)
    return rows


def fixed_T_limit(data):
    constants, zeros, xi_half, variance, s2, s4, sigma, tails, mu = data
    t,w = quadrature(0,50,zeros,64)
    psi = np.exp(-sigma*t*t/2)*np.prod(1-t[:,None]**2/zeros[None,:]**2, axis=1)
    return {"U": 50, "T": 100, "density_error_at_zero_as_K_tends_to_infinity": float(np.dot(w,psi)/np.pi-density_kernel(0,xi_half))}


def identity_checks(mu):
    mp.mp.dps = 70
    records = []
    for z in [mp.mpf(".3"), mp.mpc(".3", ".2"), mp.mpc("1.3", ".2")]:
        kmax = 10000
        logprod = sum(mu[k]*mp.log(mp.sin(mp.pi*z/k)/(mp.pi*z/k)) for k in range(1,kmax+1) if mu[k])
        prod = mp.exp(logprod)
        records.append({"z": str(z), "K": kmax, "product": str(prod), "target": str(1-z*z), "absolute_error": str(abs(prod-(1-z*z)))})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--degrees", type=int, nargs="+", default=[24, 48])
    args = parser.parse_args()
    started = time.time()
    data = scalar_constants()
    result = {"status": "computed", "script": str(Path(__file__).name), "scalar_constants": data[0], "table_runs": {}, "identity_checks": identity_checks(data[-1]), "fixed_T_limit": fixed_T_limit(data)}
    for degree in args.degrees:
        result["table_runs"][str(degree)] = table_run(degree, data)
    if len(args.degrees) > 1:
        fields = ["actual_error_on_abs_x_le_3", "refined_bound_as_printed", "crude_bound_as_printed", "crude_from_full_log_bound_with_4_over_3"]
        result["max_change_between_first_last_resolution"] = {
            field: max(abs(a[field]-b[field]) for a,b in zip(result["table_runs"][str(args.degrees[0])], result["table_runs"][str(args.degrees[-1])])) for field in fields
        }
    result["elapsed_seconds"] = time.time() - started
    (ROOT / "mobius_results.json").write_text(json.dumps(result, indent=2) + "\n")
    out = ["computed: independent reconstruction of Table 2", "Scalar constants below are certified Arb balls where so stated.", json.dumps(data[0], indent=2),
           "K U m_K actual_error refined_as_printed crude_as_printed alternative_4_over_3_crude corrected_tighter_refined argmax refined/actual percent_peak"]
    for row in result["table_runs"][str(args.degrees[-1])]:
        out.append(" ".join(str(row[key]) for key in ["K", "U", "m_K", "actual_error_on_abs_x_le_3", "refined_bound_as_printed", "crude_bound_as_printed", "crude_from_full_log_bound_with_4_over_3", "refined_bound_corrected_tighter_tail", "argmax_abs_x", "refined_over_actual", "percent_of_f0"]))
    out.append("Table values are computed quadrature/optimization, not certified enclosures.")
    (ROOT / "mobius_results.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
