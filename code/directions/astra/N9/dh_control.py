"""Independent, bounded Davenport and Heilbronn increment covariance control.

Run with timeout 30 python3 -I dh_control.py (or a larger explicit timeout).
The double precision eigensolver provides evidence, not a sign certificate.
"""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh, toeplitz
from scipy.optimize import minimize_scalar
from scipy.special import digamma, polygamma


KAPPA = (math.sqrt(10.0 - 2.0 * math.sqrt(5.0)) - 2.0) / (math.sqrt(5.0) - 1.0)
ROOT_BETA = 0.8085171824566374
ROOT_T = 85.6993484853775922


def log_derivative_coefficients(limit):
    """Solve (sum c_n/n^s)(sum a_n/n^s) = sum a_n log(n)/n^s."""
    n = np.arange(limit + 1, dtype=np.int64)
    pattern = np.array([0.0, 1.0, KAPPA, -KAPPA, -1.0], dtype=np.longdouble)
    a = pattern[n % 5]
    logs = np.zeros(limit + 1, dtype=np.longdouble)
    logs[1:] = np.log(n[1:].astype(np.longdouble))
    c = a * logs
    for d in range(2, limit // 2 + 1):
        if c[d] != 0:
            last = limit // d
            c[2 * d : (last + 1) * d : d] -= c[d] * a[2 : last + 1]
    roots = np.sqrt(np.maximum(1, n).astype(np.longdouble))
    weights = c / roots
    return a, c, logs, np.cumsum(weights), np.cumsum(weights * logs)


def archimedean(t):
    if t == 0.0:
        return 0.0
    # The omitted positive exponential tail is below
    # exp(-2*t*(K+3/4))/(4*(K+3/4)^2*(1-exp(-2*t))).
    cutoff = max(1, int(math.ceil(40.0 / (2.0 * t))))
    k = np.arange(cutoff, dtype=np.float64) + 0.75
    exponential = np.sum(np.exp(-2.0 * k * t) / (k * k))
    return 0.5 * t * (math.log(5.0 / math.pi) + digamma(0.75)) + 0.25 * (polygamma(1, 0.75) - exponential)


def psi_grid(delta, max_n, data):
    _, _, logs, cw, cwl = data
    out = np.zeros(max_n + 1)
    for j in range(1, max_n + 1):
        t = np.longdouble(str(delta)) * j
        # searchsorted uses high precision logarithms, avoiding exp cutoff errors.
        upper = int(np.searchsorted(logs[1:], t, side="right"))
        arithmetic = t * cw[upper] - cwl[upper]
        out[j] = archimedean(float(t)) - float(arithmetic)
    return out


def increment_covariances(psi):
    out = np.empty(len(psi) - 1)
    out[0] = psi[1]
    out[1:] = 0.5 * (psi[2:] + psi[:-2] - 2.0 * psi[1:-1])
    return out


def lowest(cov, n, vector=False):
    answer = eigh(toeplitz(cov[:n]), subset_by_index=[0, 0], check_finite=False, eigvals_only=not vector, driver="evr")
    if vector:
        return float(answer[0][0]), answer[1][:, 0]
    return float(answer[0])


def quartet_covariance(delta, n, beta=ROOT_BETA):
    gamma = complex(ROOT_T, beta - 0.5)
    m = np.arange(n)
    return 4.0 * np.real((1.0 - np.cos(gamma * delta)) / (gamma * gamma) * np.cos(gamma * delta * m))


def vector_diagnostics(cov, delta, n, eigenvalue, v):
    matrix = toeplitz(cov[:n])
    residual = np.linalg.norm(matrix @ v - eigenvalue * v)
    times = delta * np.arange(n)
    freq = np.linspace(0.0, math.pi / delta, 12001)
    power = np.empty(len(freq))
    for j in range(0, len(freq), 1000):
        power[j : j + 1000] = np.abs(np.exp(1j * freq[j : j + 1000, None] * times) @ v) ** 2
    peak = int(np.argmax(power))
    left, right = freq[max(0, peak - 1)], freq[min(len(freq) - 1, peak + 1)]
    optimize = minimize_scalar(lambda w: -abs(np.dot(np.exp(1j * w * times), v)) ** 2, bounds=(left, right), method="bounded")
    q4 = float(v @ toeplitz(quartet_covariance(delta, n)) @ v)
    q4line = float(v @ toeplitz(quartet_covariance(delta, n, beta=0.5)) @ v)
    return {
        "eigenpair_residual": float(residual),
        "matrix_frobenius_norm": float(np.linalg.norm(matrix)),
        "vector_l1_norm": float(np.linalg.norm(v, 1)),
        "discrete_fourier_power_peak": float(optimize.x),
        "known_quartet_quadratic_form": q4,
        "quartet_moved_to_line_quadratic_form": q4line,
        "same_witness_quadratic_form_after_quartet_replacement": eigenvalue - q4 + q4line,
    }


def find_transition(delta, data, max_a=12.0):
    max_n = int(round(max_a / delta))
    psi = psi_grid(delta, max_n, data)
    cov = increment_covariances(psi)
    samples = {}

    def evaluate(n):
        if n not in samples:
            samples[n] = lowest(cov, n)
        return samples[n]

    lower, upper = 1, min(max_n, max(2, int(round(2.0 / delta))))
    evaluate(lower)
    while evaluate(upper) >= 0 and upper < max_n:
        lower, upper = upper, min(max_n, upper * 2)
    if evaluate(upper) >= 0:
        return {"delta": delta, "first_negative": None, "samples": samples}
    while upper - lower > 1:
        middle = (lower + upper) // 2
        if evaluate(middle) < 0:
            upper = middle
        else:
            lower = middle
    value, vector = lowest(cov, upper, vector=True)
    results = {
        "delta": delta,
        "last_nonnegative_n": lower,
        "last_nonnegative_a": lower * delta,
        "last_nonnegative_lambda_min": evaluate(lower),
        "first_negative_n": upper,
        "first_negative_a": upper * delta,
        "first_negative_lambda_min": value,
        "samples": samples,
        "diagnostics": vector_diagnostics(cov, delta, upper, value, vector),
    }
    np.savez(Path(__file__).with_name("dh_witness_" + str(delta).replace(".", "p") + ".npz"), delta=delta, n=upper, vector=vector, covariance=cov[:upper], psi=psi[:upper + 1])
    return results


def independent_zero_check():
    try:
        import mpmath as mp
    except ImportError:
        return {"available": False}
    mp.mp.dps = 45
    kappa = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)

    def dh(s):
        return mp.power(5, -s) * (mp.zeta(s, mp.mpf(1) / 5) + kappa * mp.zeta(s, mp.mpf(2) / 5) - kappa * mp.zeta(s, mp.mpf(3) / 5) - mp.zeta(s, mp.mpf(4) / 5))

    def completed(s):
        return mp.power(5 / mp.pi, s / 2) * mp.gamma((s + 1) / 2) * dh(s)

    s = mp.mpc("0.8085171824566374", "85.6993484853775922")
    root = mp.findroot(dh, (s, s + mp.mpf("0.00001")), tol=mp.mpf("1e-40"))
    test = mp.mpc("0.8", "4.5")
    return {"available": True, "mpmath_version": mp.__version__, "working_dps": mp.mp.dps, "input_zero_absolute_residual": str(abs(dh(s))), "refined_root_real": str(root.real), "refined_root_imaginary": str(root.imag), "refined_root_absolute_residual": str(abs(dh(root))), "functional_equation_relative_defect": str(abs(completed(test) - completed(1 - test)) / abs(completed(test))), "certified": False}


def main():
    started = time.monotonic()
    data = log_derivative_coefficients(int(math.ceil(math.exp(12.0))))
    output = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "kappa": KAPPA,
        "coefficient_limit": len(data[1]) - 1,
        "coefficients_2_to_12": {str(n): float(data[1][n]) for n in range(2, 13)},
        "coefficient_identity_max_residual_n_le_100": max(abs(float(sum(data[1][d] * data[0][n // d] for d in range(1, n + 1) if n % d == 0) - data[0][n] * data[2][n])) for n in range(2, 101)),
        "signs_certified": False,
        "grids": [],
    }
    for delta in (0.025, 0.0125, 0.005):
        result = find_transition(delta, data)
        output["grids"].append(result)
        print(json.dumps(result), flush=True)
        Path(__file__).with_name("dh_results.json").write_text(json.dumps(output, indent=2) + "\n")
    output["zero_check"] = independent_zero_check()
    output["elapsed_seconds"] = time.monotonic() - started
    Path(__file__).with_name("dh_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"elapsed_seconds": output["elapsed_seconds"], "zero_check": output["zero_check"]}), flush=True)


if __name__ == "__main__":
    main()
