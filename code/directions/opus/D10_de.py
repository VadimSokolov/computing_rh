#!/usr/bin/env python3
"""D10: global search (differential evolution) for F = (z^2+1) * prod lower pairs * prod reals * exp(beta z).
Top pair fixed at (0,1); other pairs have height <= 1.  Objective: w*/1, the largest w<1 where the tilt has a
nonreal zero (coarse grid during the search, checked on a fine grid at the end).  Floating point, exploratory.
Usage: timeout 580 python3 -I D10_de.py secs_per_config
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution


def nonreal(z, beta, w, tol=1e-6):
    co = (np.poly(z - 1j * w) * np.exp(1j * beta * w)).real
    co = co / np.max(np.abs(co))
    k = 0
    while k < len(co) - 1 and abs(co[k]) < 1e-13:
        k += 1
    co = co[k:]
    if len(co) <= 1:
        return False
    r = np.roots(co)
    return np.any(np.abs(r.imag) > tol * np.maximum(1.0, np.abs(r)))


def wstar(z, beta, n):
    grid = np.linspace(1.0, 0.0, n + 1)[1:-1]
    last = 1.0
    for w in grid:
        if nonreal(z, beta, w):
            lo, hi = w, last
            for _ in range(22):
                mid = 0.5 * (lo + hi)
                if nonreal(z, beta, mid):
                    lo = mid
                else:
                    hi = mid
            return hi
        last = w
    return 0.0


def build(x, np_, nr):
    beta = x[0]
    z = [1j, -1j]
    k = 1
    for _ in range(np_):
        m, h = x[k], x[k + 1]
        k += 2
        z += [m + 1j * h, m - 1j * h]
    for _ in range(nr):
        z.append(x[k])
        k += 1
    return np.array(z, dtype=complex), beta


def main():
    secs = float(sys.argv[1]) if len(sys.argv) > 1 else 60
    configs = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 3)]
    for np_, nr in configs:
        bounds = [(-6, 6)] + [(-5, 5), (0.05, 1.0)] * np_ + [(-5, 5)] * nr
        t0 = time.time()

        class Stop(Exception):
            pass

        best = [9.0, None]

        def f(x):
            z, beta = build(x, np_, nr)
            v = wstar(z, beta, 80)
            if v < best[0] and v > 0:
                best[0], best[1] = v, x.copy()
            if time.time() - t0 > secs:
                raise Stop
            return v

        try:
            differential_evolution(f, bounds, popsize=25, maxiter=100000, tol=0, seed=np_ * 10 + nr,
                                   polish=False, init="sobol", updating="immediate", mutation=(0.4, 1.0), recombination=0.8)
        except Stop:
            pass
        z, beta = build(best[1], np_, nr)
        fine = wstar(z, beta, 4000)
        print(f"pairs+1={np_+1} reals={nr}: coarse {best[0]:.4f} fine {fine:.4f} beta {beta:.3f} x {np.array2string(best[1][1:],precision=3)}", flush=True)


if __name__ == "__main__":
    main()
