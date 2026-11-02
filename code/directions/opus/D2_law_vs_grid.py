#!/usr/bin/env python3
"""D2: compare the closed form margin law and the proved upper bound with the book's Weil grid.

Reads (never writes) data/wl_grid*.json, data/wl_float*.json, data/zeros_8000.txt of the book.
Run: timeout 120 python3 -I D2_law_vs_grid.py
Status of the output: floating point, heuristic law against computed eigenvalues.
"""
import json, math, os, sys
ROOT = '.'
OUT = os.path.join(ROOT, 'misc/research/directions/opus/D2_law_vs_grid.json')


def load_grid():
    cells = {}
    for fn in ['wl_float', 'wl_float2', 'wl_grid', 'wl_grid2', 'wl_grid3']:
        for e in json.load(open(f'{ROOT}/data/{fn}.json')):
            if 'eps_up_log10' in e:
                v = float(e['eps_up_log10'])
            elif 'eps' in e:
                m = e['eps']
                v = float(m.split('e')[0]) if 'e' not in m else math.log10(float(m.split('e')[0])) + float(m.split('e')[1])
            else:
                continue
            cells[(e['X'], e['N'])] = v
    return cells


def load_zeros(n=8000):
    z = []
    with open(f'{ROOT}/data/zeros_8000.txt') as fh:
        for line in fh:
            line = line.strip()
            if line:
                z.append(float(line))
            if len(z) >= n:
                break
    return z


def law(N, L):
    """Closed form E(N,L) = -log eps_N to leading order, with the mode count N+1/2 optional."""
    T = math.pi * N / L
    b_inf = 4 * math.pi * math.exp(2 * L)
    if T >= b_inf:
        return b_inf, math.pi / 2, b_inf, 0.0
    lo, hi = 1e-12, math.pi / 2
    g = lambda phi: (4 * L / math.pi) * phi - math.log(T / (4 * math.pi * math.sin(phi)))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0:
            hi = mid
        else:
            lo = mid
    phi = 0.5 * (lo + hi)
    b = T / math.sin(phi)
    slope = 4 * math.log(1 / math.tan(phi / 2))
    return b + N * slope, phi, b, slope


def bound_thmB(N, L, z):
    """log of (4 S_N / L) prod_{k<=N} (T_k/gamma_k)^4, or None if gamma_k < T_k for some k<=N."""
    s = 0.0
    for k in range(1, N + 1):
        Tk = math.pi * k / L
        if z[k - 1] < Tk:
            return None
        s += math.log(Tk / z[k - 1])
    # S_N: sum over k>N of gamma_k^-2, first 8000 ordinates plus an integral tail
    SN = sum(1 / g ** 2 for g in z[N:])
    T = z[-1]
    SN += (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T)
    return math.log(4 * SN / L) + 4 * s


def main():
    cells = load_grid()
    z = load_zeros()
    rows = []
    for (X, N), lg in sorted(cells.items()):
        L = 0.5 * math.log(X)
        E, phi, b, slope = law(N, L)
        B = bound_thmB(N, L, z) if N <= 8000 else None
        rows.append(dict(X=X, N=N, L=L, minus_ln_eps=-lg * math.log(10), E=E, phi=phi, b=b, slope=slope,
                         offset=-lg * math.log(10) - E, bound_ln=B))
    json.dump(rows, open(OUT, 'w'), indent=0)
    print('cells', len(rows))
    print('    X    N    -ln eps      E_law    offset   ratio   b      phi   slope   ThmB(-ln)')
    want = {2, 5, 13, 50, 200, 1000, 10000, 100000, 1000000}
    for r in rows:
        if r['X'] in want and r['N'] in (8, 12, 16, 20, 30, 40, 60, 80, 100, 200, 300, 400, 600):
            B = r['bound_ln']
            print(f"{r['X']:7d} {r['N']:4d} {r['minus_ln_eps']:9.2f} {r['E']:9.2f} {r['offset']:8.2f} {r['minus_ln_eps']/r['E']:6.3f} "
                  f"{r['b']:8.1f} {r['phi']:5.3f} {r['slope']:6.3f} {(-B if B is not None else float('nan')):9.2f}")


if __name__ == '__main__':
    main()
