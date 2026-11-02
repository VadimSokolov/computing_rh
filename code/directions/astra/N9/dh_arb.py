"""Arb construction and finite-matrix certificates for the DH control."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
import sys
import math
import json
import time
from pathlib import Path
sys.path.insert(0, '/private/tmp/d1lib')
sys.path.insert(0, str(Path(__file__).resolve().parent))
from flint import arb, ctx
from dh_certify import certify_matrix


def covariance(d_string, n, bits=192):
    ctx.prec = bits
    d = arb(d_string)
    limit = math.ceil(math.exp(float(d_string) * n)) + 1
    logs = [arb(0)] + [arb(k).log() for k in range(1, limit + 1)]
    assert logs[-1] > d * n
    kappa = ((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
    pattern = [arb(0), arb(1), kappa, -kappa, arb(-1)]
    a = [pattern[k % 5] for k in range(limit + 1)]
    coefficients = [a[k] * logs[k] for k in range(limit + 1)]
    for k in range(2, limit // 2 + 1):
        for multiple in range(2 * k, limit + 1, k):
            coefficients[multiple] -= coefficients[k] * a[multiple // k]
    weights = [arb(0), arb(0)] + [coefficients[k] / arb(k).sqrt() for k in range(2, limit + 1)]
    pi = arb.pi()
    linear = ((5 / pi).log() - arb.const_euler() + pi / 2 - 3 * arb(2).log()) / 2
    trigamma = pi * pi - 8 * arb.const_catalan()
    values = [arb(0)]
    for j in range(1, n + 1):
        t = d * j
        cutoff = max(1, math.ceil((bits * math.log(2) + 40) / (2 * float(d_string) * j)))
        ratio = (-2 * t).exp()
        power = (-arb('1.5') * t).exp()
        exponential_sum = arb(0)
        for k in range(cutoff):
            r = arb(k) + arb('0.75')
            exponential_sum += power / (r * r)
            power *= ratio
        endpoint = arb(cutoff) + arb('0.75')
        tail = power / (4 * endpoint * endpoint * (1 - ratio))
        value = linear * t + (trigamma - exponential_sum) / 4 + arb(0, tail.upper())
        for k in range(2, limit + 1):
            difference = t - logs[k]
            if difference > 0:
                value -= weights[k] * difference
            elif not difference < 0:
                raise ArithmeticError('Unresolved arithmetic cutoff')
        values.append(value)
    return [values[1]] + [(values[j + 1] + values[j - 1] - 2 * values[j]) / 2 for j in range(1, n)]


def main():
    started = time.monotonic()
    results = {'arithmetic': 'python-flint 0.9.0 Arb, 192 bits, with explicit IEEE dot-product residual bounds', 'grids': []}
    for d, n in [('0.025', 218), ('0.0125', 317), ('0.005', 743)]:
        c = covariance(d, n)
        entry = {'delta': d, 'last_positive': certify_matrix(c[:-1]), 'first_negative': certify_matrix(c)}
        results['grids'].append(entry)
        Path(__file__).with_name('dh_arb_results.json').write_text(json.dumps(results, indent=2) + '\n')
        print(json.dumps({'delta': d, 'positive_n': n - 1, 'positive_certified': entry['last_positive']['positive_definite'], 'negative_n': n, 'negative_certified': entry['first_negative']['negative_rayleigh'], 'positive_eigenvalue_lower_bound': entry['last_positive']['lower_bound_display'], 'negative_rayleigh_upper_bound': entry['first_negative']['upper_bound_display']}), flush=True)
    results['elapsed_seconds'] = time.monotonic() - started
    Path(__file__).with_name('dh_arb_results.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps({'elapsed_seconds': results['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
