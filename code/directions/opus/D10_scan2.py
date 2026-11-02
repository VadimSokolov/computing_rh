#!/usr/bin/env python3
"""D10: grid scans and local refinement for F = (z^2+1) ((z-m)^2+h^2) exp(beta z), then more pairs.

Top pair at height 1 (abscissa 0).  Reports the least ratio w*/1 found.  Floating point, exploratory.
Usage: timeout 300 python3 -I D10_scan2.py
"""
import time
import numpy as np


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


def wstar(z, beta, n=300):
    delta = np.max(np.abs(z.imag))
    grid = np.linspace(delta, 0.0, n + 1)[1:-1]
    last = delta
    for w in grid:
        if nonreal(z, beta, w):
            lo, hi = w, last
            for _ in range(25):
                mid = 0.5 * (lo + hi)
                if nonreal(z, beta, mid):
                    lo = mid
                else:
                    hi = mid
            return hi
        last = w
    return 0.0


def cfg(pairs, reals=()):
    z = list(reals)
    for m, h in pairs:
        z += [m + 1j * h, m - 1j * h]
    return np.array(z, dtype=complex)


def scan_two():
    best = (9, None)
    for m in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.6, 3, 4, 5, 7, 10]:
        for sign in (1, -1):
            for h in np.linspace(0.3, 1.0, 15):
                for beta in np.linspace(0.5, 5.0, 46):
                    r = wstar(cfg([(0, 1), (sign * m, h)]), beta, 120)
                    if r < best[0]:
                        best = (r, (sign * m, h, beta))
    return best


def refine(x0, build, secs, seed=0):
    rng = np.random.default_rng(seed)
    x = np.array(x0, float)
    f = build(x)
    t0 = time.time()
    step = 0.1
    stall = 0
    while time.time() - t0 < secs and step > 1e-4:
        y = x + step * rng.normal(size=len(x))
        g = build(y)
        if g < f:
            x, f = y, g
            stall = 0
        else:
            stall += 1
            if stall > 40:
                step *= 0.7
                stall = 0
    return f, x


def main():
    t0 = time.time()
    b = scan_two()
    print("two pairs grid best", b, "t", round(time.time() - t0, 1), flush=True)

    def f2(x):
        m, h, beta = x
        h = min(abs(h), 1.0)
        return wstar(cfg([(0, 1), (m, h)]), beta, 400)

    f, x = refine(b[1], f2, 40)
    print("two pairs refined", f, x, flush=True)

    # three pairs: the top pair and two lower ones, random restarts
    rng = np.random.default_rng(5)

    def f3(x):
        m1, h1, m2, h2, beta = x
        return wstar(cfg([(0, 1), (m1, min(abs(h1), 1.0)), (m2, min(abs(h2), 1.0))]), beta, 150)

    best = (9, None)
    t1 = time.time()
    while time.time() - t1 < 60:
        x0 = [rng.uniform(-3, 3), rng.uniform(0.2, 1), rng.uniform(-3, 3), rng.uniform(0.2, 1), rng.uniform(0, 5)]
        f, x = refine(x0, f3, 4, seed=int(rng.integers(1 << 30)))
        if f < best[0]:
            best = (f, x)
            print("three pairs", f, x, flush=True)

    # pair plus real zeros near the pair, with field: n real zeros and top pair
    def f4(x):
        beta = x[0]
        return wstar(cfg([(0, 1)], reals=tuple(x[1:])), beta, 200)

    for nr in (1, 2, 3, 4):
        best = (9, None)
        t1 = time.time()
        while time.time() - t1 < 25:
            x0 = np.concatenate(([rng.uniform(0, 5)], rng.uniform(-3, 3, nr)))
            f, x = refine(x0, f4, 3, seed=int(rng.integers(1 << 30)))
            if f < best[0]:
                best = (f, x)
        print("pair + reals", nr, best, flush=True)


if __name__ == "__main__":
    main()
