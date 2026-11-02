#!/usr/bin/env python3
"""Arithmetic for the tables of B.md (exact formulas evaluated in floating point). python3 -I tables.py"""
from math import acos, sqrt, pi, log, e, exp, erf

def g(nu):          # exponent function of Theorem 1(f)
    return sqrt(1 - nu * nu) - nu * acos(nu)
def N_count(T):     # main term of Riemann and von Mangoldt
    return T / (2 * pi) * (log(T / (2 * pi)) - 1) + 7 / 8
def rho(H):
    return log(H / (2 * pi)) / (2 * pi)
Phi = lambda z: 0.5 * (1 + erf(z / sqrt(2)))

print("g(nu):", {nu: round(g(nu), 4) for nu in (0, .25, .5, .75, .9, .99)})
T = 3e12
for nu in (0.25, 0.5, 0.75):
    N = sqrt(T / (2 * pi * nu))
    print(f"T=3e12, nu={nu}: N={N:.4g}, primes below N about {N/log(N):.3g} (x/log x), pi N^2 g(nu)={pi*N*N*g(nu):.4g}, log10 size={-pi*N*N*g(nu)/log(10):.4g}")
# Theorem 2: cluster length and count, bits
for H, Tw in ((1e3, 500.), (1e6, 1.), (1e6, 5e5), (3e12, 1.), (3e12, 1.5e12)):
    ell = (H - Tw) / e
    n = ell * rho(H + ell / 2)
    print(f"H={H:.3g} T={Tw:.3g}: ell={ell:.4g}, n about {n:.4g}, bits n/ln2 = {n/log(2):.4g}, RS length sqrt(H/2pi)={sqrt(H/(2*pi)):.4g}, X proved H/2pi^0.088={(H/(2*pi))**0.088:.3g}")
# constant c_star: 1 + c/2 = log(2/c)
lo, hi = 0.1, 1.0
for _ in range(80):
    c = 0.5 * (lo + hi)
    if 1 + c / 2 - log(2 / c) > 0: hi = c
    else: lo = c
print(f"c_star={c:.5f}, exponent a_star=c_star/(2 pi)={c/(2*pi):.5f}")
# prior tail table
NT = N_count(T)
print(f"N(3e12) about {NT:.4g}")
for name, S in (("exchangeable zeros, Laplace", 1 / (NT + 1)), ("1/(1+log T)", 1 / (1 + log(T))), ("1/(1+log log T)", 1 / (1 + log(log(T)))), ("1/(1+log log log T)", 1 / (1 + log(log(log(T)))))):
    print(f"{name:32s} S={S:.3e}  P(RH | verified, w=1/2)={1/(1+S):.14f}")
# two point error: Neyman Pearson total error 2 Phi(-kappa/2)
for kappa in (1e-6, 1e-2, 1, 2, 4, 8):
    print(f"kappa={kappa:g}: minimal sum of errors {2*Phi(-kappa/2):.6f}")
print("||d_eta psi||_L2 at a=1/2:", sqrt(3 / (4 * pi * 0.5 ** 5)))
