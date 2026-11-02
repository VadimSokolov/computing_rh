#!/usr/bin/env python3
"""D2: the equilibrium identities and the explicit test function of Theorem 5, against the Weil grid.

Part 1 checks numerically: int_b^inf q = N - N0(b), U(t) = 0 for t > b, U(0) = E/2.
Part 2 builds the zero set (ordinates up to b, band zeros beta_j from the smooth count) and evaluates
its Rayleigh quotient on the zero side (under RH): first 8000 ordinates, then the density with mean sin^2 = 1/2.
Run: timeout 300 python3 -I D2_testfunction.py
Status: floating point in log space; the rational form is exact, the tail is a model above height 8148.
"""
import json, math
import numpy as np
from scipy import integrate
ROOT = '.'


def params(N, L):
    T = math.pi * N / L
    lo, hi = 1e-12, math.pi / 2
    g = lambda p: (4 * L / math.pi) * p - math.log(T / (4 * math.pi * math.sin(p)))
    if g(hi) <= 0:
        return math.pi / 2, 4 * math.pi * math.exp(2 * L), T
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if g(m) > 0 else (m, hi)
    p = 0.5 * (lo + hi)
    return p, T / math.sin(p), T


def q_band(x, phi, b, L):
    """q(x) = (1/pi)[ (1/2) log(1+w) + (2L/pi)(phi - psi) ], w = exp(-2 eta), psi = arctan(tan phi tanh eta)."""
    x = np.asarray(x, dtype=float)
    r = x / b
    eta_exp = r + np.sqrt(np.maximum(r * r - 1, 0.0))   # e^eta
    w = eta_exp ** -2
    tanh = (1 - w) / (1 + w)
    tp = math.tan(min(phi, math.pi / 2 - 1e-14))
    # phi - psi by the subtraction formula for arctan, stable for large x
    phi_minus_psi = np.arctan(tp * (2 * w / (1 + w)) / (1 + tp * tp * tanh))
    return (0.5 * np.log1p(w) + (2 * L / math.pi) * phi_minus_psi) / math.pi


def N0(s):
    return s / (2 * math.pi) * math.log(s / (2 * math.pi * math.e))


def check_identities(N, L):
    phi, b, T = params(N, L)
    mass = integrate.quad(lambda u: float(q_band(b * math.exp(u), phi, b, L)) * b * math.exp(u), 0, 60, limit=400)[0]
    rho0 = lambda s: math.log(s / (2 * math.pi)) / (2 * math.pi)

    def U(t):
        pts_lo = [x for x in (t, T, 2 * math.pi) if 0 < x < b]
        a = integrate.quad(lambda s: math.log(abs(t * t - s * s) + 1e-300) * rho0(s), 0, b, points=pts_lo, limit=800)[0]
        c = integrate.quad(lambda s: math.log(abs(t * t - s * s) + 1e-300), 0, T, points=[t] if t < T else None, limit=400)[0] * L / math.pi
        f = lambda u: math.log(abs(t * t - (b * math.exp(u)) ** 2) + 1e-300) * float(q_band(b * math.exp(u), phi, b, L)) * b * math.exp(u)
        pt = [math.log(t / b)] if t > b else None
        d = integrate.quad(f, 0, 60, points=pt, limit=800)[0]
        return a + d - c
    E = b + 4 * N * math.log(1 / math.tan(phi / 2))
    return dict(N=N, L=L, b=b, phi=phi, mass=mass, N_minus_N0=N - N0(b), U0=U(1e-9), E_half=E / 2, U_15b=U(1.5 * b), U_3b=U(3 * b),
                U_half_b=U(0.5 * b), U_half_b_closed=closed_U(0.5 * b, phi, b, L))


def closed_U(t, phi, b, L):
    def u(s):
        th = math.asin(min(s / b, 1.0))
        return -0.5 * math.acos(min(s / b, 1.0)) + (L / math.pi) * math.log(abs(math.sin(phi - th) / math.sin(phi + th)))
    return -integrate.quad(u, t, b, points=[b * math.sin(phi)] if t < b * math.sin(phi) else None, limit=400)[0]


def smooth_count_table(N, L, phi, b):
    """A(s) = N - D(s) on s >= b with D(s) = int_s^inf q."""
    u = np.linspace(0, 40, 400001)
    s = b * np.exp(u)
    qs = q_band(s, phi, b, L) * s
    du = u[1] - u[0]
    tail = np.concatenate([np.cumsum(0.5 * (qs[1:] + qs[:-1])[::-1] * du)[::-1], [0.0]])
    return s, N - tail - qs[-1] * 1e-40


def zero_set(N, L, gam):
    phi, b, T = params(N, L)
    K = int(np.searchsorted(gam, b, side='right'))
    if K > N:
        K = N
    s, A = smooth_count_table(N, L, phi, b)
    J = N - K
    targets = K + np.arange(1, J + 1) - 0.5
    betas = np.interp(targets, A, s)
    low = targets <= A[0]
    betas[low] = b + np.arange(1, low.sum() + 1) * (2 * math.pi / math.log(b / (2 * math.pi)))
    return np.concatenate([gam[:K], np.sort(betas)]), K, phi, b, T


def rayleigh(N, L, gam):
    z, K, phi, b, T = zero_set(N, L, gam)
    Tn = math.pi * np.arange(1, N + 1) / L

    def logG(t):
        t = np.atleast_1d(t)
        out = np.zeros_like(t)
        for i in range(0, len(t), 500):
            tt = t[i:i + 500][:, None]
            out[i:i + 500] = (np.log(np.abs(1 - tt ** 2 / z ** 2)).sum(1) - np.log(np.abs(1 - tt ** 2 / Tn ** 2)).sum(1))
        return out
    g = gam[K:]
    lg = logG(g)
    terms = 2 * lg + 2 * np.log(np.abs(np.sin(L * g)) + 1e-300) - 2 * np.log(g)
    # tail above the last ordinate: density with mean sin^2 = 1/2
    tgrid = g[-1] * np.exp(np.linspace(0, 16, 1601))
    lgt = logG(tgrid)
    dens = np.log(tgrid / (2 * math.pi)) / (2 * math.pi)
    integrand = 0.5 * dens * np.exp(2 * lgt) / tgrid ** 2 * tgrid
    du = 16 / 1600
    tail = np.sum(0.5 * (integrand[1:] + integrand[:-1])) * du + math.exp(2 * lgt[-1]) * 0.5 * (math.log(tgrid[-1] / (2 * math.pi)) + 1) / (2 * math.pi * tgrid[-1])
    m = terms.max()
    lnsum = m + math.log(np.exp(terms - m).sum() + tail * math.exp(-m)) if tail > 0 else m + math.log(np.exp(terms - m).sum())
    # exact norm through the values at the lattice points
    logF = np.zeros(N)
    for j in range(N):
        y = Tn[j] ** 2
        others = np.delete(Tn, j)
        logF[j] = -math.log(2) + np.log(np.abs(1 - y / z ** 2)).sum() - np.log(np.abs(1 - y / others ** 2)).sum()
    nm = np.logaddexp.reduce(np.concatenate([[math.log(1.0)], math.log(2) + 2 * logF]))
    return -(math.log(4 / L) + lnsum - nm), -(math.log(4 / L) + lnsum), K, phi, b, z


def main():
    gam = np.array([float(x) for x in open(f'{ROOT}/data/zeros_8000.txt') if x.strip()])
    print('identities')
    for N, X in [(40, 13), (100, 1000), (600, 1000), (60, 50), (20, 1e6)]:
        L = 0.5 * math.log(X)
        r = check_identities(N, L)
        print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
    grid = {(r['X'], r['N']): r for r in json.load(open(f'{ROOT}/misc/research/directions/opus/D2_law_vs_grid.json'))}
    print('    X    N   K     b  -ln eps(grid)    E      -ln RQ(exact norm)  -ln RQ(F(0) norm)')
    out = []
    for X, N in [(13, 8), (13, 20), (13, 40), (13, 60), (50, 40), (50, 100), (50, 200), (200, 100), (200, 300), (1000, 30), (1000, 100), (1000, 200), (1000, 600), (100000, 100), (100000, 600)]:
        L = 0.5 * math.log(X)
        a, c, K, phi, b, z = rayleigh(N, L, gam)
        e = grid[(X, N)]
        out.append(dict(X=X, N=N, K=K, b=b, grid=e['minus_ln_eps'], E=e['E'], rq_exact=a, rq_F0=c))
        print(f"{X:7d} {N:4d} {K:4d} {b:7.1f} {e['minus_ln_eps']:10.2f} {e['E']:9.2f} {a:10.2f} {c:10.2f}")
    json.dump(out, open(f'{ROOT}/misc/research/directions/opus/D2_testfunction.json', 'w'), indent=0)


if __name__ == '__main__':
    main()
