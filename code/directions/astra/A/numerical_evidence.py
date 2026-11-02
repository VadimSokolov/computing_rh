"""Bounded checks for report.md. All numerical values are floating point.

Run: timeout 30 python3 -I numerical_evidence.py
No book scripts or other research reports are imported.
"""

import json
import math
from pathlib import Path
import mpmath as mp


def q(width, x):
    return width / (math.pi * (width * width + x * x))


def pair(width, delta, x):
    return q(width - delta, x) + q(width + delta, x)


def pair_difference(width, delta, x):
    # Algebraic form avoids subtracting nearly equal Cauchy tails.
    z = x * x
    d2 = delta * delta
    a2 = width * width
    return (2 * width * d2 * (a2 - 3 * z - d2)
            / (math.pi * (a2 + z)
               * ((width - delta) ** 2 + z)
               * ((width + delta) ** 2 + z)))


out = {}
a = 0.5
out["pair_quadratic"] = [
    {"delta": d, "data_sup": 2 * d * d / (math.pi * a * (a * a - d * d)),
     "ratio_to_delta_squared": 2 / (math.pi * a * (a * a - d * d))}
    for d in (0.1, 0.01, 0.001)
]
out["pair_noise_threshold"] = [
    {"eps": eps, "delta_star": math.sqrt(math.pi * a ** 3 * eps / (1 + math.pi * a * eps))}
    for eps in (1e-4, 1e-8, 1e-12)
]
T, delta, y = 10.0, 0.2, 0.1
grid = [-T + 2 * T * i / 4000 for i in range(4001)]
out["far_quartet"] = []
for gamma in (30.0, 60.0, 120.0, 240.0):
    differences = [abs(pair_difference(a, delta, t - gamma)
                       + pair_difference(a, delta, t + gamma)) for t in grid]
    out["far_quartet"].append({
        "T": T, "gamma": gamma, "delta": delta, "y": y,
        "grid_data_sup": max(differences),
        "proved_data_bound": 12 * a * delta ** 2 / (math.pi * (gamma - T) ** 4),
        "gamma_fourth_times_center_difference": gamma ** 4 * abs(2 * pair_difference(a, delta, gamma)),
        "target_at_gamma": pair(y, delta, 0) + pair(y, delta, 2 * gamma),
    })
out["far_quartet_asymptotic_constant"] = 12 * a * delta ** 2 / math.pi

out["pole_crossing"] = []
gamma, y = 5.0, 0.2
for h in (1e-2, 1e-4, 1e-6):
    dm, dp = y - h, y + h
    diffs = [abs(pair(a, dp, t-gamma) + pair(a, dp, t+gamma)
                 - pair(a, dm, t-gamma) - pair(a, dm, t+gamma)) for t in grid]
    out["pole_crossing"].append({
        "h": h, "grid_data_sup": max(diffs),
        "target_left": pair(y, dm, 0) + pair(y, dm, 2*gamma),
        "target_right": pair(y, dp, 0) + pair(y, dp, 2*gamma),
    })

out["holder_and_bandwidth"] = [
    {"y": b, "exponent": b / a, "amplitude_at_noise_1e_12": (1e-12) ** (b / a),
     "omega_to_make_exp_minus_y_omega_1e_6": math.log(1e6) / b,
     "log10_X": math.log10(1e6) / b,
     "decimal_digits_lost": (a-b) * math.log10(1e6) / b}
    for b in (0.25, 0.1, 0.01)
]

def dh(s):
    k = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    return 5 ** (-s) * (mp.zeta(s, mp.mpf(1)/5)
                        + k * mp.zeta(s, mp.mpf(2)/5)
                        - k * mp.zeta(s, mp.mpf(3)/5)
                        - mp.zeta(s, mp.mpf(4)/5))


def log_derivative_density(s, is_dh):
    if is_dh:
        z = mp.log(5/mp.pi)/2 + mp.digamma((s+1)/2)/2 + mp.diff(dh, s)/dh(s)
    else:
        z = 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.diff(mp.zeta, s)/mp.zeta(s)
    return mp.re(z)/mp.pi


out["dh_comparison"] = []
for digits in (40, 60):
    with mp.workdps(digits):
        root = mp.findroot(dh, (mp.mpc("0.8085", "85.6993"), mp.mpc("0.8086", "85.6994")))
        t = mp.im(root)
        out["dh_comparison"].append({
            "digits": digits, "root": mp.nstr(root, 30),
            "residual": mp.nstr(abs(dh(root)), 5),
            "values": [{"theta": theta,
                        "dh_u": mp.nstr(log_derivative_density(mp.mpc(theta,t), True), 20),
                        "zeta_u": mp.nstr(log_derivative_density(mp.mpc(theta,t), False), 20)}
                       for theta in ("1.0", "0.9", "0.8", "0.7", "0.6")]
        })

# Algebraic and asymptotic checks, independent of grid extrema.
assert out["dh_comparison"][0]["root"] == out["dh_comparison"][1]["root"]
assert out["dh_comparison"][0]["values"] == out["dh_comparison"][1]["values"]
for d in (0.01, 0.1, 0.2):
    for x in (0, 0.5, 2, 10):
        assert math.isclose(pair_difference(a,d,x), pair(a,d,x)-2*q(a,x), rel_tol=1e-8, abs_tol=1e-15)
for row in out["far_quartet"]:
    assert row["grid_data_sup"] <= row["proved_data_bound"]
    assert row["target_at_gamma"] < 0
for row in out["pole_crossing"]:
    assert row["target_left"] > 0 > row["target_right"]

folder = Path(__file__).resolve().parent / "results"
folder.mkdir(exist_ok=True)
payload = json.dumps(out, indent=2)
(folder / "numbers.txt").write_text(payload + "\n")
print(payload)
