#!/usr/bin/env python3
"""D1 exploration, part 2 (double precision, NOT a certificate).
First order and second order effect of the prime 2 term, thresholds of the
archimedean form, weight interval, reduced form in the Legendre basis, and the
Davenport Heilbronn rung. Imports the closed forms from D1_explore.py.
Usage: timeout 500 python3 -I -u D1_explore2.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.linalg as sla
from scipy.special import spherical_jn, digamma, sici
from numpy.polynomial.legendre import leggauss
import D1_explore as E

Ltop = 0.5 * np.log(3.0)
Lbot = 0.5 * np.log(2.0)
K2 = E.K2


def eig(M):
    return sla.eigh(0.5 * (M + M.T))


print("=== table N=300: sector L mu1 mu2 (A+P) | -T00 second lam1 lam2 (W) overlap")
for L in [Lbot + 1e-9, 0.36, 0.38, 0.40, 0.41, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, Ltop]:
    d = E.matrices(L, 300)
    for sec in "eo":
        A = d["A" + sec] + d["P" + sec]
        T = d["T" + sec]
        mu, V = eig(A)
        lam, U = eig(A - T)
        t00 = V[:, 0] @ T @ V[:, 0]
        ov = (V[:, 0] @ U[:, 0]) ** 2
        print("%s %.5f mu %.4e %.4e | -T00 %+.4e second %+.4e lam %.4e %.4e ov %.3f | mu2-|T| %+.3f" % (sec, L, mu[0], mu[1], -t00, lam[0] - mu[0] + t00, lam[0], lam[1], ov, mu[1] - K2 / 2))

