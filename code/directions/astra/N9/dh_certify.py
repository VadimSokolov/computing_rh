"""Arb certificates using a floating approximate congruence and residual bounds.

certify_matrix(c) takes Arb balls for a real symmetric Toeplitz first row.
All floating computations only select approximations. Bounds use Arb and
the standard IEEE 754 dot-product forward error bound.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
import sys
sys.path.insert(0, '/private/tmp/d1lib')
from flint import arb
import numpy as np
from scipy.linalg import eigh, toeplitz


def float_arb(x):
    num, den = float(x).as_integer_ratio()
    return arb(num) / den


def frobenius_upper(a):
    # An exact sum of squared dyadic upper bounds for absolute entries.
    # The small integer rounding is independent of hardware accumulation.
    bits = 256
    total = 0
    for x in np.asarray(a).flat:
        numerator, denominator = abs(float(x)).as_integer_ratio()
        scaled = ((numerator << bits) + denominator - 1) // denominator
        total += scaled * scaled
    return (arb(total) / (arb(2) ** (2 * bits))).sqrt()


def certify_matrix(c):
    n = len(c)
    row = np.array([float(x.mid()) for x in c])
    # These dyadic size checks exclude product overflow and underflow.
    # First products have denominator at most 2^490, second at most 2^530.
    assert n <= 2**20
    assert all(np.isfinite(x) and abs(x) <= 2.0**400 and x.as_integer_ratio()[1] <= 2**450 for x in row)
    matrix = toeplitz(row)
    diagonal, basis = eigh(matrix, check_finite=False)
    scale = 1 << 40
    integers = np.rint(basis * scale).astype(np.int64)
    assert all(abs(int(x)) <= 2 * scale for x in integers.flat)
    basis = integers.astype(float) / scale
    # The dyadic denominator of either matrix product is at most the row
    # denominator times 2^80. These checks exclude product underflow and
    # overflow; the final residual subtraction is bounded separately.
    assert max(float(x).as_integer_ratio()[1].bit_length() - 1 for x in row) + 80 < 1000
    assert np.max(np.abs(row)) * n * n < 1e200
    assert np.max(np.abs(basis)) <= 2
    qnorm2 = arb(sum(int(x) * int(x) for x in integers.flat)) / (arb(scale) * scale)
    dnorm = frobenius_upper(matrix)
    product = basis.T @ (matrix @ basis)
    gram = basis.T @ basis
    residual = product.copy()
    residual[np.diag_indices(n)] -= diagonal
    gram_residual = gram - np.eye(n)
    unit = arb(2) ** -53
    gamma = n * unit / (1 - n * unit)
    # The last subtraction can round, including gradual underflow.
    underflow = n * (arb(2) ** -1074)
    residual_norm = (frobenius_upper(residual) + underflow) / (1 - unit)
    gram_norm = (frobenius_upper(gram_residual) + underflow) / (1 - unit)
    multiply_error = gamma * (2 + gamma) * qnorm2 * dnorm
    input_norm_squared = arb(0)
    for k in range(n):
        error = arb(abs(c[k] - float_arb(row[k])).upper())
        weight = n if k == 0 else 2 * (n - k)
        input_norm_squared += weight * error * error
    input_error = input_norm_squared.sqrt() * qnorm2
    error = residual_norm + multiply_error + input_error
    rho = gram_norm + gamma * qnorm2
    if not rho < 1:
        return {'n': n, 'certified': False, 'reason': 'Gram residual does not prove nonsingularity', 'rho': str(rho)}
    mu = float_arb(diagonal[0])
    lower_bound = mu - (error + abs(mu) * rho) / (1 - rho)
    # A separate exact dyadic Rayleigh vector gives an upper bound.
    v = [int(x) for x in integers[:, 0]]
    correlation = [sum(v[j] * v[j + k] for j in range(n - k)) for k in range(n)]
    rayleigh = (c[0] * correlation[0] + 2 * sum((c[k] * correlation[k] for k in range(1, n)), arb(0))) / correlation[0]
    return {
        'n': n,
        'certified': bool(lower_bound > 0 or rayleigh < 0),
        'positive_definite': bool(lower_bound > 0),
        'negative_rayleigh': bool(rayleigh < 0),
        'eigenvalue_lower_bound_ball': str(lower_bound),
        'eigenvalue_upper_bound_ball': str(rayleigh),
        'lower_bound_display': float(lower_bound.lower()),
        'upper_bound_display': float(rayleigh.upper()),
        'approximate_eigenvalue': float(diagonal[0]),
        'congruence_residual_bound': str(error),
        'gram_residual_bound': str(rho),
        'computed_congruence_residual_frobenius_bound': str(residual_norm),
        'floating_matrix_multiplication_bound': str(multiply_error),
        'analytic_input_ball_bound': str(input_error),
        'qnorm_squared': str(qnorm2),
        'assumptions': 'IEEE binary64 arithmetic, round to nearest, BLAS dot accumulation with at most n products and n additions; Arb elementary enclosures',
    }
