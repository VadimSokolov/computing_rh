#!/usr/bin/env python3
"""D2: the pencil form (Theorem 1) against the grid, and the fixed N law (Theorem 3) as L grows.

Zero side under RH: first 8000 ordinates of data/zeros_8000.txt (50 digits), above them the density with mean sin^2 = 1/2.
Run: timeout 600 python3 -I D2_fixedN.py
Status: mpmath floating point at 150 digits; the tail is a model, so the agreement with the grid is computed, not certified.
"""
import json, math
import mpmath as mp
mp.mp.dps = 150
ROOT = '.'
OPUS = ROOT + '/misc/research/directions/opus'
GAM = [mp.mpf(x.strip()) for x in open(ROOT + '/data/zeros_8000.txt') if x.strip()]
TOP = GAM[-1]


def tail_moment(j):
    """(1/2) int_TOP^inf rho0(t) t^(-2-2j) dt, rho0 = log(t/2pi)/(2pi)."""
    p = mp.mpf(2 + 2 * j)
    return mp.mpf(1) / 2 * (1 / (2 * mp.pi)) * TOP ** (1 - p) * (mp.log(TOP / (2 * mp.pi)) / (p - 1) + 1 / (p - 1) ** 2)


def moments(L, jmax, N=0):
    """m_j = sum sin^2(L gamma) r^2 gamma^(-2-2j), r = prod_{n=1}^N (1 - T_n^2/t^2)^(-1) (N=0: the pure moments of Theorem 3)."""
    Tn2 = [(mp.pi * n / L) ** 2 for n in range(1, N + 1)]
    m = [mp.mpf(0)] * (jmax + 1)
    for g in GAM:
        w = 1 / g ** 2
        s2 = mp.sin(L * g) ** 2
        r = mp.mpf(1)
        for y in Tn2:
            r /= (1 - y * w)
        val = s2 * r * r * w
        for j in range(jmax + 1):
            m[j] += val
            val *= w
    return [m[j] + tail_moment(j) for j in range(jmax + 1)]


def pencil_eps(N, L):
    """Least eigenvalue on V_N from the pencil of Theorem 1, basis p(w), a_n = (-1)^n P(y_n)/prod (y_n - y_j)."""
    m = moments(L, 2 * N, N)
    A = mp.matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            A[i, j] = 8 * m[i + j]
    y = [(mp.pi * n / L) ** 2 for n in range(N + 1)]
    Z = mp.matrix(N + 1, N + 1)
    for n in range(N + 1):
        den = mp.mpf(1)
        for j in range(N + 1):
            if j != n:
                den *= (y[n] - y[j])
        for i in range(N + 1):
            yp = mp.mpf(1) if N - i == 0 else y[n] ** (N - i)
            Z[n, i] = yp / den * (-1) ** n
    D = mp.matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            D[i, j] = L * sum((2 if n == 0 else 1) * Z[n, i] * Z[n, j] for n in range(N + 1))
    C = mp.cholesky(D)
    Ci = mp.inverse(C)
    S = Ci * A * Ci.T
    S = (S + S.T) / 2
    E = mp.eigsy(S, eigvals_only=True)
    return min(E)


def c_N(N):
    return 4 * mp.pi ** (4 * N) * mp.factorial(2 * N) ** 4 / mp.factorial(4 * N)


def hankel_ratio(m, N):
    def det(n):
        M = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                M[i, j] = m[i + j]
        return mp.det(M)
    return det(N + 1) / (det(N) if N > 0 else 1)


def main():
    grid = {(r['X'], r['N']): r for r in json.load(open(OPUS + '/D2_law_vs_grid.json'))}
    print('Theorem 1 against the grid: X, N, -ln eps(pencil), -ln eps(grid)')
    for X, N in [(5, 8), (13, 8), (13, 12), (50, 10), (1000, 8), (10000, 8), (1000000, 8), (13, 16), (1000, 16)]:
        if (X, N) not in grid:
            continue
        e = pencil_eps(N, mp.log(X) / 2)
        print(' ', X, N, mp.nstr(-mp.log(e), 8), round(grid[(X, N)]['minus_ln_eps'], 4))
    print('Theorem 3: ratio eps_N / (c_N L^(-4N-1) h_N(L)), and L^2 (ratio - 1)')
    for N in (0, 1, 2, 3, 4, 6):
        for L in (3, 5, 10, 30, 100):
            if N > 0 and L < 2 * mp.pi * N / GAM[0]:
                continue
            L = mp.mpf(L)
            e = pencil_eps(N, L)
            m0 = moments(L, 2 * N, 0)
            h = hankel_ratio(m0, N)
            pred = c_N(N) * L ** (-4 * N - 1) * h
            ra = e / pred
            print(f"  N={N} L={int(L):3d} eps={mp.nstr(e, 5):>12} ratio={mp.nstr(ra, 8)} L^2(ratio-1)={mp.nstr(L * L * (ra - 1), 5)}")
    print('annealed comparison: h_N(L) against (1/2) H_N/H_{N-1}, H_N the Hankel determinant of the power sums')
    Mj = [sum(g ** (-2 - 2 * j) for g in GAM) + 2 * tail_moment(j) for j in range(13)]
    for N in (1, 3, 6):
        ann = hankel_ratio(Mj, N) / 2
        vals = []
        for L in [mp.mpf(40) + mp.mpf(k) * mp.mpf('1.37') for k in range(16)]:
            vals.append(hankel_ratio(moments(L, 2 * N, 0), N) / ann)
        vals = [float(v) for v in vals]
        print(f"  N={N} annealed={mp.nstr(ann, 5)} ratio h_N/annealed over 16 values of L in [40,61]: min {min(vals):.3f} median {sorted(vals)[8]:.3f} max {max(vals):.3f}")


if __name__ == '__main__':
    main()
