"""Outward dyadic ball verification of the independent DH control.

The ball arithmetic is the shared, newly written local ball_review.py.
No source-review script is imported. Invoke with timeout and Python -I.
"""

import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ball_review as br
import mpmath as mp
import numpy as np

B = br.B


def record(ball):
    return {"midpoint": ball.mid(), "radius": ball.radius(), "bounds_display_only": ball.bounds(), "exact_midpoint_integer": str(ball.m), "exact_radius_integer": str(ball.r), "scale_bits": br.BITS}


def psi_balls(d_string, size, bits):
    br.set_precision(bits)
    d = B(d_string)
    max_t = float(d_string) * size
    limit = int(math.ceil(math.exp(max_t))) + 1
    logs = [B(0)] + [br.blog(B(n)) for n in range(1, limit + 1)]
    assert (logs[-1] - size * d).positive()
    kappa = (br.bsqrt(10 - 2 * br.bsqrt(5)) - 2) / (br.bsqrt(5) - 1)
    pattern = [B(0), B(1), kappa, -kappa, B(-1)]
    a = [pattern[n % 5] for n in range(limit + 1)]
    coeff = [a[n] * logs[n] for n in range(limit + 1)]
    for n in range(2, limit // 2 + 1):
        for multiple in range(2 * n, limit + 1, n):
            coeff[multiple] = coeff[multiple] - coeff[n] * a[multiple // n]
    weights = [B(0), B(0)] + [coeff[n] / br.bsqrt(n) for n in range(2, limit + 1)]
    pi = B.from_iv(mp.iv.pi)
    euler = B.from_iv(mp.iv.euler)
    catalan = B.from_iv(mp.iv.catalan)
    linear = (br.blog(5 / pi) - euler + pi / 2 - 3 * br.blog(2)) / 2
    trigamma = pi * pi - 8 * catalan
    out = [B(0)]
    for j in range(1, size + 1):
        t = d * j
        # This formula chooses K. The subsequent ball tail formula proves
        # the actual tail enclosure, independently of this float choice.
        cutoff = max(1, math.ceil((bits * math.log(2.0) + 40.0) / (2 * float(d_string) * j)))
        ratio = br.bexp(-2 * t)
        exponential = br.bexp(-B("1.5") * t)
        power = exponential
        exponential_sum = B(0)
        for k in range(cutoff):
            r = B(k) + B("0.75")
            exponential_sum = exponential_sum + power / (r * r)
            power = power * ratio
        tail = power / (4 * (B(cutoff) + B("0.75")) * (B(cutoff) + B("0.75")) * (1 - ratio))
        value = linear * t + (trigamma - exponential_sum) / 4
        value = value + B.raw(0, max(0, tail.m + tail.r))
        for n in range(2, limit + 1):
            difference = t - logs[n]
            if difference.positive():
                value = value - weights[n] * difference
            elif not difference.negative():
                raise ArithmeticError("Unresolved arithmetic cutoff")
        out.append(value)
    covariance = [out[1]] + [(out[j + 1] + out[j - 1] - 2 * out[j]) / 2 for j in range(1, size)]
    return covariance


def rayleigh(covariance, vector):
    # The witness is an exact dyadic rational vector, obtained by rounding
    # the exploratory eigenvector. Its sign is independent of normalization.
    scale = 1 << 40
    integers = [int(round(float(x) * scale)) for x in vector]
    autocorrelation = [sum(integers[j] * integers[j + lag] for j in range(len(integers) - lag)) for lag in range(len(integers))]
    total = covariance[0] * autocorrelation[0]
    for lag in range(1, len(covariance)):
        total = total + 2 * covariance[lag] * autocorrelation[lag]
    return total / autocorrelation[0]


def main():
    started = time.monotonic()
    output = {"arithmetic": "Custom exact-integer dyadic midpoint-radius balls, elementary functions and constants enclosed by mpmath.iv; not Arb", "grids": []}
    for d, n, bits, do_levinson in [("0.025", 218, 2048, True), ("0.0125", 317, 256, False), ("0.005", 743, 256, False)]:
        covariance = psi_balls(d, n, bits)
        witness_path = Path(__file__).with_name("dh_witness_" + d.replace(".", "p") + ".npz")
        witness = np.load(witness_path)["vector"]
        q = rayleigh(covariance, witness)
        item = {"delta": d, "n": n, "bits": bits, "negative_witness_certified": q.negative(), "rayleigh": record(q)}
        if do_levinson:
            pivots, reflections, fail = br.ball_levinson(covariance)
            item["levinson"] = {"processed": len(pivots), "first_nonpositive_pivot_index": fail, "all_preceding_pivots_positive": all(v.positive() for v in pivots[:-1]), "last_pivot_negative": pivots[-1].negative(), "last_pivot": record(pivots[-1]), "preceding_pivot": record(pivots[-2]) if len(pivots) > 1 else None}
        output["grids"].append(item)
        Path(__file__).with_name("dh_ball_results.json").write_text(json.dumps(output, indent=2) + "\n")
        print(json.dumps({"delta": d, "n": n, "bits": bits, "rayleigh_midpoint": q.mid(), "rayleigh_radius": q.radius(), "negative_witness_certified": q.negative(), "levinson": item.get("levinson", {}).get("first_nonpositive_pivot_index")}), flush=True)
    output["elapsed_seconds"] = time.monotonic() - started
    Path(__file__).with_name("dh_ball_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"elapsed_seconds": output["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
