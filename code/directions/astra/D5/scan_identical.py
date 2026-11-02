"""Small floating point search for an identical-spin obstruction.

Run: timeout 30 python3 -I scan_identical.py
The candidate c(J)=J/2 is tested, not assumed or proved.
"""

import json
from pathlib import Path

import numpy as np


def width(coeff):
    roots = np.roots(np.asarray(coeff)[::-1])
    return float(np.max(np.abs(np.log(np.abs(roots)))))


def main():
    rng = np.random.default_rng(20261101)
    closest = None
    violations = []
    active = 0
    cases = 3000
    for k in range(cases):
        m = int(rng.integers(1, 5))
        weights = np.exp(rng.uniform(-6, 0, m + 1))
        weights /= weights.sum()
        probs = np.r_[weights[:0:-1] / 2, weights[0], weights[1:] / 2]
        xs = np.arange(-m, m + 1)
        j = float(np.exp(rng.uniform(np.log(1e-4), np.log(3))))
        delta = width(probs)
        pair = np.zeros(4 * m + 1)
        for i, x in enumerate(xs):
            for l, y in enumerate(xs):
                pair[i + l] += probs[i] * probs[l] * np.exp(j * x * y)
        coupled = width(pair)
        allowed2 = max(delta**2 - j / 2, 0)
        violation = coupled**2 - allowed2
        item = {"case": k, "m": m, "J": j, "probabilities": probs.tolist(),
                "single_width": delta, "coupled_width": coupled,
                "candidate_bound_violation_squared": violation,
                "squared_decrease_over_J": (delta**2-coupled**2)/j}
        if coupled > 1e-6 and delta > 1e-6:
            active += 1
            if closest is None or item["squared_decrease_over_J"] < closest["squared_decrease_over_J"]:
                closest = item
        if violation > 1e-7:
            violations.append(item)
    summary = {"status": "exploratory double precision, not certified",
               "seed": 20261101, "cases": cases, "cases_off_axis_after_coupling": active,
               "violation_tolerance_squared": 1e-7,
               "candidate": "W^2 <= max(delta^2-J/2, 0)",
               "closest_off_axis_case": closest, "violations": violations}
    Path(__file__).with_name("identical_scan.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
