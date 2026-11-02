"""Independent numerical controls for referee_validity.md.

These high-precision calculations are controls, not interval certificates.
The mathematical counterexamples are proved analytically in the report.
"""

import json
from pathlib import Path

import mpmath as mp


mp.mp.dps = 60


def xi(s):
    if s == 0 or s == 1:
        return mp.mpf("0.5")
    return s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def phi(u):
    return 2 * mp.fsum(
        (2 * mp.pi**2 * n**4 * mp.exp(9 * u / 2) - 3 * mp.pi * n**2 * mp.exp(5 * u / 2))
        * mp.exp(-mp.pi * n**2 * mp.exp(2 * u))
        for n in range(1, 41)
    )


def symmetric_mgf(z, h):
    return mp.cosh(z) * mp.cosh(z - 2 * h) * mp.cosh(z + 2 * h) / mp.cosh(2 * h) ** 2


def boundary_derivative(w, h):
    z = h + 1j * w
    return mp.re(mp.tanh(z) + mp.tanh(z - 2 * h) + mp.tanh(z + 2 * h))


results = {"precision_digits": mp.mp.dps, "classification": "computed"}
results["phi_normalization"] = []
for level in [mp.mpf("0.5"), mp.mpf("1"), mp.mpf("2")]:
    integral = 2 * mp.quad(
        lambda u: phi(u) * mp.cosh((level - mp.mpf("0.5")) * u),
        [0, mp.mpf("0.25"), mp.mpf("0.5"), 1, 2, 3],
    )
    target = xi(level)
    results["phi_normalization"].append(
        {"level": str(level), "integral": str(integral), "xi": str(target), "error": str(integral - target)}
    )

h = mp.mpf("0.1")
grid = [mp.pi * j / 1000 for j in range(1001)]
errors = [abs(boundary_derivative(w, h) - mp.re(mp.tanh(3 * h + 1j * w))) for w in grid]
results["positive_boundary_with_right_zeros"] = {
    "h": str(h),
    "exact_minimum": str(mp.tanh(3 * h)),
    "grid_minimum": str(min(boundary_derivative(w, h) for w in grid)),
    "maximum_cancellation_error": str(max(errors)),
    "right_zero_real_part": str(2 * h),
    "right_zero_mgf_absolute_value": str(abs(symmetric_mgf(2 * h + 1j * mp.pi / 2, h))),
}

gamma = mp.mpf("2")
theta = mp.mpf("0.3")
v = mp.mpf("0.7")
f = lambda z: mp.log(1 + z**2 / gamma**2)
lhs = f(theta + v) + f(theta - v) - 2 * f(theta)
rhs_undoubled_kernel = 2 * mp.quad(
    lambda x: (1 - mp.cos(v * x)) * mp.cos(theta * x) * mp.exp(-gamma * x) / x,
    [0, 1, mp.inf],
)
results["symmetric_kernel_factor"] = {
    "left_side": str(lhs),
    "integral_with_kernel_exp_minus_gamma_x": str(rhs_undoubled_kernel),
    "left_minus_twice_integral": str(lhs - 2 * rhs_undoubled_kernel),
}

output = Path(__file__).with_name("referee_validity_checks.json")
output.write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
