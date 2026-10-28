"""lo_common.py: shared helpers of the computation group LO (Chapter 30, the Thorin measure at level one and the backward Cauchy flow).

What it provides: the sampled level one density u_1 (import of the book's cx_em.u1_half, python-flint balls at 128 bits, checked at 192 bits), the
analytic value u_1(0) = (2 + gamma_E - log 4 pi)/(2 pi) that replaces the pole sample at v = 0, the zeros of data/zeros_1700.txt, Cauchy kernels
and a peak matcher.
Conventions: xi(s) = s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)/2, u_1(v) = pi^{-1} Re (xi'/xi)(1 + i v) = sum over zeros rho of
P_{1-sigma}(v - gamma), P_a(v) = (a/pi)/(a^2 + v^2), the sum over both signs of gamma (H's measure, no factor 1/2).
Written from the mathematics (the collaborator's scripts were not read, run or copied).
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import cx_em

GAMMA_E = 0.5772156649015328606
U10 = (2 + GAMMA_E - math.log(4 * math.pi)) / (2 * math.pi)

def zeros(n=1700):
    return cx_em.zeros(n)

def u1_grid(vmax, h, procs=1, log=print):
    """u_1 on 0, h, .., vmax with the exact value at v = 0 (cx_em.u1_half puts 1e-4 there)."""
    u = cx_em.u1_half(vmax, h, procs, log)
    u[0] = U10
    return u

def P(a, v):
    return (a / np.pi) / (a * a + v * v)

def match(peaks, zs):
    """for each zero the distance to the nearest peak; returns (max distance, number of zeros whose nearest peak is used by no other zero)."""
    peaks = np.asarray(peaks); zs = np.asarray(zs)
    if len(peaks) == 0: return float('nan'), 0
    d = np.abs(zs[:, None] - peaks[None, :])
    j = d.argmin(1)
    return float(d.min(1).max()), int(len(np.unique(j)))
