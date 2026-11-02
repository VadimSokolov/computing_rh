#!/usr/bin/env python3
"""D10 numerical evidence (floating point; nothing here is certified).

E0  constants of Theorem 1.
E1  the laws X_{n,eps} with mgf cosh(s)^n ((cosh s + A)^2 + eps^2), A = 1/2: ratio w*/delta against n,
    by the up crossing test of Theorem 2, cross checked by Chebyshev roots for small n.
E3  the fixed point of Remark cm:prop:gap against the true gap, for one pair among real zeros in a field.
E5  the exact example of Proposition 3 (the set of couplings is not an interval).
Usage: timeout 580 python3 -I D10_evidence.py
"""
from fractions import Fraction as Fr
import numpy as np
from numpy.polynomial import chebyshev as Ch

T_C = 2.798386045783887
BETA_C = np.sqrt(1 + T_C ** 2)
KAPPA1 = 0.3927607377918873


def r_star_formula(beta):
    """r*(beta) = inverse of arccos(r)/r, valid for beta <= beta_c."""
    lo, hi = 1e-6, 1.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if np.arccos(m) / m > beta:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


# ---------- E1 ----------

def delta_of(A, eps):
    """Height of the zeros cos z = -A -+ i eps: solve A^2/cosh^2 y + eps^2/sinh^2 y = 1."""
    lo, hi = 1e-12, 5.0
    for _ in range(200):
        y = 0.5 * (lo + hi)
        if A ** 2 / np.cosh(y) ** 2 + eps ** 2 / np.sinh(y) ** 2 > 1:
            lo = y
        else:
            hi = y
    return 0.5 * (lo + hi)


def eps_of(A, delta):
    return np.sinh(delta) * np.sqrt(1 - A ** 2 / np.cosh(delta) ** 2)


def upcross(n, A, eps, w, delta, npts=4001):
    """True if Phi_w rises through a level near the cluster."""
    x0 = np.arccos(-A / np.cosh(delta))
    u = np.linspace(-3 * delta, 3 * delta, npts)
    x = x0 + u
    x = x[(x > 0) & (x < np.pi)]
    z = x + 1j * w
    y = np.cos(z)
    theta = np.angle(-y)            # arg(-cos(x+iw)), the distance of Theta_w above -pi
    psi = np.angle((y + A) ** 2 + eps ** 2)  # principal, by the proof of Theorem 2
    phi = n * theta + psi           # Phi_w + n pi
    lev = np.floor((phi - np.pi / 2) / np.pi)
    runmin = np.minimum.accumulate(lev)
    return np.any(lev > runmin)


def ratio_law(n, A, delta, grid=np.linspace(0.999, 0.05, 1900)):
    eps = eps_of(A, delta)
    last = 1.0
    for r in grid:
        if not upcross(n, A, eps, r * delta, delta):
            return last
        last = r
    return 0.0


def cheb_coeffs(n, A, eps):
    # exact product by repeated Chebyshev multiplication
    c = np.array([1.0])
    for _ in range(n):
        c = Ch.chebmul(c, np.array([0.0, 1.0]))
    return Ch.chebmul(c, np.array([A * A + eps * eps + 0.5, 2 * A, 0.5]))


def ferro_cheb(c, w, tol=1e-6):
    k = np.arange(len(c))
    cw = c * np.cosh(k * w)
    r = Ch.chebroots(cw / np.max(np.abs(cw)))
    return np.all(np.abs(r.imag) < tol) and np.all(np.abs(r.real) <= 1 + tol)


def e1():
    print("E1: X_{n,eps}, A=1/2, tan(omega0)=sqrt 3; beta = n sqrt(3) delta")
    print("   n   best_beta  delta      ratio(w*/delta)   limit r*(beta)")
    for n in (5, 10, 20, 50, 100, 300, 1000, 10000, 100000):
        best = (9, None)
        for beta in np.linspace(1.5, 3.05, 32):
            d = beta / (n * np.sqrt(3))
            if d > 0.9:
                continue
            r = ratio_law(n, 0.5, d, np.linspace(0.999, 0.05, 479))
            if 0 < r < best[0]:
                best = (r, beta, d)
        # refine beta
        beta0 = best[1]
        for beta in np.linspace(beta0 - 0.05, beta0 + 0.05, 21):
            d = beta / (n * np.sqrt(3))
            r = ratio_law(n, 0.5, d)
            if 0 < r < best[0]:
                best = (r, beta, d)
        print(f"{n:7d}  {best[1]:.4f}  {best[2]:.4e}  {best[0]:.4f}   {r_star_formula(min(best[1], BETA_C)):.4f}", flush=True)
    print("cross check by Chebyshev roots (n small): largest w/delta at which a root leaves [-1,1]")
    for n, beta in ((5, 2.2), (10, 2.5), (14, 2.6)):
        d = beta / (n * np.sqrt(3))
        eps = eps_of(0.5, d)
        c = cheb_coeffs(n, 0.5, eps)
        last = 1.0
        for r in np.linspace(0.999, 0.05, 1900):
            if not ferro_cheb(c, r * d):
                break
            last = r
        print(f"  n={n} beta={beta} up-crossing ratio {ratio_law(n,0.5,d):.4f}  roots ratio {last:.4f}", flush=True)


# ---------- E3 ----------

def true_gap(reals, phi):
    """w* for (z^2+1) prod(z-a) exp(phi z), height 1, by the roots (grid then bisection)."""
    z = np.concatenate(([1j, -1j], np.array(reals, dtype=complex)))

    def nonreal(w):
        co = (np.poly(z - 1j * w) * np.exp(1j * phi * w)).real
        co = co / np.max(np.abs(co))
        k = 0
        while k < len(co) - 1 and abs(co[k]) < 1e-13:
            k += 1
        co = co[k:]
        if len(co) <= 1:
            return False
        r = np.roots(co)
        return np.any(np.abs(r.imag) > 1e-7 * np.maximum(1.0, np.abs(r)))

    last = 1.0
    for w in np.linspace(0.9995, 0.0005, 2000):
        if nonreal(w):
            lo, hi = w, last
            for _ in range(40):
                m = 0.5 * (lo + hi)
                if nonreal(m):
                    lo = m
                else:
                    hi = m
            return 1 - hi
        last = w
    return 1.0


def fixed_point(reals, phi):
    """Solve eps = dist(arg G(iw), pi Z)^2 / (4 (b + 1/(2-eps))), w = 1-eps, b = -Im G'/G(iw), from eps=0 upward."""
    reals = np.array(reals, dtype=float)

    def rhs(eps):
        w = 1 - eps
        z = 1j * w
        arg = np.sum(np.angle(z - reals)) + phi * w
        dist = abs((arg + np.pi / 2) % np.pi - np.pi / 2)
        b = -np.sum(1 / (z - reals)).imag
        return dist ** 2 / (4 * (b + 1 / (2 - eps)))

    grid = np.linspace(0, 0.999, 4000)
    v = np.array([rhs(e) - e for e in grid])
    idx = np.where(v[:-1] * v[1:] <= 0)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    lo, hi = grid[i], grid[i + 1]
    for _ in range(50):
        m = 0.5 * (lo + hi)
        if (rhs(lo) - lo) * (rhs(m) - m) <= 0:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def e3():
    print("E3: book fixed point against the true gap, eta=1")
    print("pure field: phi, exact gap 1-w (w=cos(phi w)), fixed point")
    for phi in (0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 2.97):
        g = true_gap([], phi)
        f = fixed_point([], phi)
        print(f"  phi={phi:4.2f} gap {g:.5f} fp {f:.5f} fp/gap {f/g:.4f}")
    print("central real zero z(z^2+1): gap", round(true_gap([0.0], 0.0), 5), "fp", round(fixed_point([0.0], 0.0), 5), "exact 1-1/sqrt3 =", round(1 - 1 / np.sqrt(3), 5))
    rng = np.random.default_rng(3)
    ratios = []
    rec = []
    for trial in range(1500):
        n = rng.integers(0, 7)
        reals = rng.uniform(-6, 6, n)
        reals = reals[np.abs(reals) > 0.05]
        phi = rng.uniform(-3, 3)
        g = true_gap(reals, phi)
        f = fixed_point(reals, phi)
        if np.isfinite(f) and g > 1e-5 and g < 0.999:
            ratios.append(f / g)
            rec.append((f / g, g, f, phi, reals))
    ratios = np.array(ratios)
    print(f"random configurations {len(ratios)}: fp/gap min {ratios.min():.4f} median {np.median(ratios):.4f} max {ratios.max():.4f}; "
          f"fraction with fp>gap*(1+1e-3) {np.mean(ratios>1.001):.4f}")
    rec.sort(key=lambda t: -t[0])
    for t in rec[:3]:
        print("   largest fp/gap:", np.round(t[0], 4), "gap", np.round(t[1], 5), "fp", np.round(t[2], 5), "phi", np.round(t[3], 3), "reals", np.round(t[4], 3))


# ---------- E5 ----------

def e5():
    print("E5: X = Y1+Y2, A1=13/12 (eta1=log 3/2), A2=5/4 (eta2=log 2), exact arithmetic")
    A1, A2 = Fr(13, 12), Fr(5, 4)
    S, P = A1 + A2, A1 * A2

    def Q(C):
        return (2 * C * C - 1, S * C, P + 1 - C * C)

    for C in (Fr(21, 20), Fr(13, 12), Fr(11, 10), Fr(6, 5), Fr(5, 4), Fr(3, 2)):
        a, b, c = Q(C)
        disc = b * b - 4 * a * c
        val = a - b + c
        roots = None
        if disc >= 0:
            r1 = (-float(b) - float(disc) ** 0.5) / (2 * float(a))
            r2 = (-float(b) + float(disc) ** 0.5) / (2 * float(a))
            roots = (round(r1, 5), round(r2, 5))
        ferro = disc >= 0 and val >= 0 and (a + b + c) >= 0 and b <= 2 * a
        print(f"  C=cosh w={C}: Q=({a},{b},{c}) Q(-1)={val} disc={float(disc):.5f} roots {roots} ferromagnetic {ferro}")
    # lower end of the first window
    lo, hi = 1.0, 13 / 12
    for _ in range(100):
        m = 0.5 * (lo + hi)
        a, b, c = 2 * m * m - 1, float(S) * m, float(P) + 1 - m * m
        if b * b - 4 * a * c >= 0:
            hi = m
        else:
            lo = m
    print("  first window: cosh w from", hi, "to 13/12; w from", np.arccosh(hi), "to", np.log(1.5))


if __name__ == "__main__":
    print("beta_c", BETA_C, "kappa_1", KAPPA1, "r_c", np.cos(np.pi - T_C))
    e5()
    e3()
    e1()
