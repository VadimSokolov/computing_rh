#!/usr/bin/env python3
"""D1: test of the plan's first lemma as stated (double precision, NOT a certificate).

M = A + P is the archimedean form with the pole, T the prime 2 form, W = M - T.
For each L and sector this prints
  lam1, lam2 of M - T and of M + T.
  By Cauchy interlacing, min over any hyperplane v^perp of a form is at most its lam2.
  So if lam2(M + T) < 0 then |T| <= M fails on v^perp for EVERY v, however chosen.
  min of W and of M + T on the complement of three candidate ground states:
    gs:  ground state of M itself
    pz:  prolate state, top eigenvector of the sinc kernel on [-L,L] with band Omega = 2 pi exp(2L)
    p1:  prolate state with band Omega = gamma_1 = 14.1347
  min of W and of M + T on the Connes and Consani subspace F(0) = F(i/2) = F(-i/2) = 0.

Usage: timeout 300 python3 -I D1_explore4.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.linalg as sla
from scipy.special import roots_legendre
from D1_explore import matrices

N = 300
NQ = 1100
xq, wq = roots_legendre(NQ)


def restricted_min(M, C):
    """Least eigenvalue of M on the orthogonal complement of the columns of C."""
    Q, _ = np.linalg.qr(C, mode="complete")
    B = Q[:, C.shape[1]:]
    R = B.T @ M @ B
    return sla.eigvalsh(0.5 * (R + R.T), subset_by_index=[0, 0])[0]


def prolate(L, Omega):
    x = L * xq
    w = L * wq
    d = x[:, None] - x[None, :]
    K = np.where(np.abs(d) < 1e-14, Omega / np.pi, np.sin(Omega * d) / (np.pi * np.where(np.abs(d) < 1e-14, 1, d)))
    s = np.sqrt(w)
    S = s[:, None] * K * s[None, :]
    vals, vecs = sla.eigh(0.5 * (S + S.T), subset_by_index=[NQ - 2, NQ - 1])
    out = {}
    for k in range(2):
        psi = vecs[:, k] / s
        if abs(psi @ (w * psi[::-1]) - 1) < 1e-6:
            out["e"] = (psi, vals[k])
        else:
            out["o"] = (psi, vals[k])
    om = np.pi * np.arange(N + 1) / L
    pe, le = out["e"]
    ce = np.empty(N + 1)
    ce[0] = (w * pe).sum() / np.sqrt(2 * L)
    ce[1:] = (np.cos(np.outer(om[1:], x)) * (w * pe)).sum(axis=1) / np.sqrt(L)
    po, lo = out["o"]
    co = (np.sin(np.outer(om[1:], x)) * (w * po)).sum(axis=1) / np.sqrt(L)
    return ce / np.linalg.norm(ce), co / np.linalg.norm(co), le, lo


print("sector L | lam1(W) lam2(W) | lam1(M+T) lam2(M+T) | on gs^perp: W, M+T | on pz^perp: W, M+T, top overlap with W ground | on p1^perp: W, M+T, overlap | on CC subspace: M, W, M+T")
for L in [0.5 * np.log(2), 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.5 * np.log(3)]:
    d = matrices(L, N)
    pze, pzo, _, _ = prolate(L, 2 * np.pi * np.exp(2 * L))
    p1e, p1o, _, _ = prolate(L, 14.134725)
    for sec in "eo":
        A = d["Ae"] if sec == "e" else d["Ao"]
        P = d["Pe"] if sec == "e" else d["Po"]
        T = d["Te"] if sec == "e" else d["To"]
        M = A + P
        W = M - T
        Wp = M + T
        lw = sla.eigvalsh(0.5 * (W + W.T), subset_by_index=[0, 1])
        vW = sla.eigh(0.5 * (W + W.T), subset_by_index=[0, 0])[1][:, 0]
        lp = sla.eigvalsh(0.5 * (Wp + Wp.T), subset_by_index=[0, 1])
        gs = sla.eigh(0.5 * (M + M.T), subset_by_index=[0, 0])[1]
        pz = (pze if sec == "e" else pzo)[:, None]
        p1 = (p1e if sec == "e" else p1o)[:, None]
        if sec == "e":
            e0 = np.zeros((N + 1, 1))
            e0[0, 0] = 1.0
            CC = np.hstack([e0, d["ap"][:, None]])
        else:
            CC = d["am"][:, None]
        print("%s %.5f | %+.4e %+.4e | %+.4e %+.4e | %+.3e %+.3e | %+.3e %+.3e %.4f | %+.3e %+.3e %.4f | %+.3e %+.3e %+.3e" % (
            sec, L, lw[0], lw[1], lp[0], lp[1],
            restricted_min(W, gs), restricted_min(Wp, gs),
            restricted_min(W, pz), restricted_min(Wp, pz), abs(pz[:, 0] @ vW),
            restricted_min(W, p1), restricted_min(Wp, p1), abs(p1[:, 0] @ vW),
            restricted_min(M, CC), restricted_min(W, CC), restricted_min(Wp, CC)), flush=True)
