#!/usr/bin/env python3
"""Linearised invisibility test: moving zeros in a window against the primes below X = e^x.

Real matrix B_{jk} = x_j sin(gamma_k x_j), x_j on (0, x_max]: the first order change of the real
zero sum 2 Re D(x) under real displacements a_k of the zeros gamma_k. The smallest singular value
measures the best move of the zeros that the band |x| <= x_max cannot see. Double precision, not certified.
Usage: python3 -I nyquist_svd.py   (reads data/zeros_8000.txt)
"""
import os
import numpy as np
from math import pi, log, sqrt, exp

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", "..", "..", "..", ".."))
Z = np.array([float(l) for l in open(os.path.join(root, "data", "zeros_8000.txt"))])
for H, ell in ((300.0, 120.0), (1000.0, 300.0), (3000.0, 300.0), (7000.0, 300.0)):
    g = Z[(Z >= H) & (Z <= H + ell)]
    c = 0.5 * (g[0] + g[-1]); n = g.size
    rho = n / (g[-1] - g[0]); rho_th = log((H + ell / 2) / (2 * pi)) / (2 * pi)
    print(f"H={H:.0f} ell={ell:.0f} n={n} empirical density {rho:.4f}, formula {rho_th:.4f}, pi rho={pi*rho:.3f}, e^(pi rho)={exp(pi*rho):.2f}, sqrt((H+ell/2)/2 pi)={sqrt((H+ell/2)/(2*pi)):.2f}")
    row = []
    for f in (0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5):
        xm = f * pi * rho
        J = int(2.5 * xm * g[-1] / pi) + 2 * n    # samples resolve the oscillation at the true ordinates
        x = (np.arange(J) + 0.5) * xm / J
        A = np.sin(np.outer(x, g)) * (x[:, None] / xm)
        s = np.linalg.svd(A, compute_uv=False)
        rank = int(np.sum(s > 1e-6 * s[0]))
        row.append((f, float(np.log10(s[-1] / s[0])), rank))
    print("   x_max/(pi rho), log10 sigma_min/sigma_max, numerical rank (1e-6) of", n, "columns:")
    print("   " + "  ".join(f"{a:.2f}:{b:.1f}:{r}" for a, b, r in row))
