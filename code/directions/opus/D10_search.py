#!/usr/bin/env python3
"""D10: random search for the local model F(z) = P(z) exp(phi z).

For real P with zeros c_j and real phi, T_w F(x) = 2 exp(phi x) Re[exp(i phi w) P(x+iw)].
w*(F) = inf{w0 : T_w F has only real zeros for every w > w0}, delta(F) = max |Im c_j|.
We minimise the ratio w*/delta over the configuration.  Floating point, exploratory.

Usage: timeout 300 python3 -I D10_search.py npairs nreal seconds [seed]
"""
import sys, time
import numpy as np


def zeros_from(params, npairs, nreal):
    """params: pair abscissae, log heights, real zeros, log phi."""
    m = params[:npairs]
    h = np.exp(params[npairs:2 * npairs])
    a = params[2 * npairs:2 * npairs + nreal]
    phi = np.exp(params[2 * npairs + nreal])
    z = list(a)
    for mj, hj in zip(m, h):
        z += [mj + 1j * hj, mj - 1j * hj]
    return np.array(z, dtype=complex), phi


def nonreal(z, phi, w, tol=1e-6):
    """True if Re[exp(i phi w) P(x+iw)] has a nonreal root."""
    co = np.poly(z - 1j * w) * np.exp(1j * phi * w)
    co = co.real
    # scale
    scale = np.max(np.abs(co))
    co = co / scale
    # drop leading coefficients that vanish to rounding
    k = 0
    while k < len(co) - 1 and abs(co[k]) < 1e-13:
        k += 1
    co = co[k:]
    if len(co) <= 1:
        return False
    r = np.roots(co)
    size = np.maximum(1.0, np.abs(r))
    return np.any(np.abs(r.imag) > tol * size)


def wstar(z, phi, n=400):
    """Largest w below delta at which the tilt has a nonreal zero (grid and bisection)."""
    delta = np.max(np.abs(z.imag))
    grid = np.linspace(delta, 0.0, n + 1)[1:-1]
    last = delta
    for w in grid:
        if nonreal(z, phi, w):
            lo, hi = w, last
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                if nonreal(z, phi, mid):
                    lo = mid
                else:
                    hi = mid
            return hi
        last = w
    return 0.0


def ratio(params, npairs, nreal, n=400):
    z, phi = zeros_from(params, npairs, nreal)
    delta = np.max(np.abs(z.imag))
    if delta <= 0:
        return 1.0
    # tol of the nonreal test is relative to max(1,|r|): normalise the top height to 1
    return wstar(z / delta, phi * delta, n)


def main():
    npairs = int(sys.argv[1])
    nreal = int(sys.argv[2])
    secs = float(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    rng = np.random.default_rng(seed)
    dim = 2 * npairs + nreal + 1
    t0 = time.time()
    best, bestx = 2.0, None
    restarts = 0
    while time.time() - t0 < secs:
        restarts += 1
        x = rng.normal(size=dim) * 1.5
        x[npairs:2 * npairs] = rng.normal(size=npairs) * 0.5
        f = ratio(x, npairs, nreal, 150)
        step = 0.5
        stall = 0
        while stall < 60 and time.time() - t0 < secs:
            y = x + step * rng.normal(size=dim)
            g = ratio(y, npairs, nreal, 150)
            if g < f:
                x, f = y, g
                stall = 0
            else:
                stall += 1
                if stall % 15 == 0:
                    step *= 0.6
        f = ratio(x, npairs, nreal, 1200)
        if 0 < f < best:
            best, bestx = f, x.copy()
            z, phi = zeros_from(x, npairs, nreal)
            d = np.max(np.abs(z.imag))
            print(f"restart {restarts} ratio {best:.5f} phi*delta {phi*d:.4f} zeros/delta:",
                  np.array2string(np.sort_complex((z - z.real.mean()) / d)[z.imag >= 0] if False else (z - 0) / d, precision=4), flush=True)
    print("best", best, "restarts", restarts)


if __name__ == "__main__":
    main()
