#!/usr/bin/env python3
"""D10: laws on the integers, L(s) = sum c_k cosh(k s), c_k >= 0, characteristic function
F(z) = sum c_k cos(k z) = polynomial in y = cos z.  The cosine tilt multiplies c_k by cosh(k w).
Ferromagnetic at w iff all n roots of sum c_k cosh(k w) T_k(y) lie in [-1,1].
(1) Random search on the coefficients for the least ratio w*/delta at each degree n.
(2) The explicit family (cos z + B)^n ((cos z + A)^2 + eps^2).
Floating point, exploratory.
Usage: timeout 400 python3 -I D10_law_search.py secs_per_degree
"""
import sys, time
import numpy as np
from numpy.polynomial import chebyshev as Ch


def is_ferro(c, w, tol=1e-7):
    k = np.arange(len(c))
    cw = c * np.cosh(k * w)
    cw = cw / np.max(np.abs(cw))
    r = Ch.chebroots(cw)
    if np.any(np.abs(r.imag) > tol):
        return False
    return np.all(np.abs(r.real) <= 1 + tol)


def delta_of(c):
    r = Ch.chebroots(c / np.max(np.abs(c)))
    z = np.arccos(r.astype(complex))
    return np.max(np.abs(z.imag))


def wstar(c, ngrid=240):
    d = delta_of(c)
    if d < 1e-5:
        return 0.0, d
    grid = np.linspace(d, 0, ngrid + 1)[1:-1]
    last = d
    for w in grid:
        if not is_ferro(c, w):
            lo, hi = w, last
            for _ in range(22):
                mid = 0.5 * (lo + hi)
                if is_ferro(c, mid):
                    hi = mid
                else:
                    lo = mid
            return hi, d
        last = w
    return 0.0, d


def obj(x, ngrid=120):
    c = np.exp(x)
    ws, d = wstar(c, ngrid)
    if d < 1e-3:
        return 1.0
    return ws / d


def search(n, secs, rng):
    t0 = time.time()
    best, bx = 2.0, None
    while time.time() - t0 < secs:
        # start from a product of three point factors, with random log perturbation
        x = np.log(1e-3 + rng.uniform(0, 1, n + 1)) + rng.normal(size=n + 1)
        f = obj(x)
        step = 0.7
        stall = 0
        while stall < 80 and time.time() - t0 < secs:
            y = x + step * rng.normal(size=n + 1)
            g = obj(y)
            if g < f and g > 0:
                x, f = y, g
                stall = 0
            else:
                stall += 1
                if stall % 20 == 0:
                    step *= 0.6
        f = obj(x, 480)
        if 0 < f < best:
            best, bx = f, x.copy()
    return best, bx


def family(n, A, B, eps):
    """Chebyshev coefficients of (y+B)^n ((y+A)^2+eps^2) with y = cos z, as c_k in cos(kz)."""
    p = Ch.cheb2poly(np.array([0.0, 1.0]))  # just to make sure the module is loaded
    q = np.array([A * A + eps * eps + 0.5, 2 * A, 0.5])
    f = Ch.chebpow(np.array([B, 1.0]), n) if n > 0 else np.array([1.0])
    return Ch.chebmul(f, q)


def main():
    secs = float(sys.argv[1]) if len(sys.argv) > 1 else 15
    rng = np.random.default_rng(1)
    print("(1) random search over laws sum c_k cos(kz)")
    for n in (2, 3, 4, 5, 6, 8, 10):
        b, x = search(n, secs, rng)
        c = np.exp(x)
        print(f"n={n} best ratio {b:.4f} delta {delta_of(c):.4f} c/c_max {np.array2string(c/np.max(c),precision=3)}", flush=True)
    print("(2) family (y+B)^n ((y+A)^2+eps^2): best ratio over B, eps at A=0.5")
    for n in (1, 2, 4, 8, 16, 32):
        best = (9, None)
        for B in (0.0, 0.1, 0.2, 0.3, 0.4):
            for eps in np.geomspace(0.002, 0.6, 40):
                c = family(n, 0.5, B, eps)
                ws, d = wstar(c, 240)
                if d > 1e-4 and ws / d < best[0] and ws > 0:
                    best = (ws / d, (B, eps, d, ws))
        print(f"n={n} best {best}", flush=True)


if __name__ == "__main__":
    main()
