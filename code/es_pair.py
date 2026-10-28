# es_pair.py: the pair model of "The functional equation makes the problem singular" (explicit sparse chapter): ||t_perp||, the profile likelihood in delta, and the
# Bayes factors of Theorem "Bayes factors" for the zero at 201.264752 and for the pair added at 60.089, on the level 1 window of the chapter's section "The computation"
# (c = 200 or 60, taper ell = 10, omega_j = 5 + 0.1 j, j < 400, a0 = 1/2; the chapter does not define the window of the pair experiment, this is the assumption).
# Conventions: f(a)_j = exp(-q^2/(2 ell^2) + i omega_j q), q = gamma - c + i a (Lemma "evenness"); data vector in R^(2J); noise N(0, tau^2) per real coordinate
# (Theorem "rate", Fisher information 1/tau^2) or tau^2/2 per real coordinate (the ReLU experiment's convention); zeros = data/zeros_1700.txt (68 digits, read to double).
# Nuisance: (gamma, a) of every zero with |gamma - c| <= W other than the central one, the central ordinate, and optionally (kappa=1) the complex amplitude of the
# central pair (2 real parameters; see note on the simple zero versus the double zero in the checks file).
# Method: linearization at the truth (nuisance shifts are of order tau, so the model is linear in them); the profile chi^2(delta) = ||P_perp(delta) (y - Phi(delta))||^2/tau^2
# with P_perp the projector orthogonal to the nuisance tangent columns at delta (truncated SVD, relative cutoff --rcond, swept for error control).
# BF01 = exp(-chi2(0)/2) / int pi(delta) exp(-chi2(delta)/2) d delta (nuisance Jacobians cancel to first order), single zero alternative: width a0+eps, eps uniform on (0, 1/2).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import numpy as np, sys
import os
A0 = float(os.environ.get('ES_A0', 0.5)); ELL = float(os.environ.get('ES_ELL', 10.0))
OM = float(os.environ.get('ES_OM0', 5)) + float(os.environ.get('ES_DOM', 0.1)) * np.arange(int(os.environ.get('ES_J', 400)))  # env vars select the source design (om0 2, Delta 0.2, J 100, ell 3, a0 0.1)
def zeros():
    return np.array([float(l.split()[0]) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_1700.txt')) if l.strip()])
def f(g, a): q = g + 1j * a; return np.exp(-q ** 2 / (2 * ELL ** 2) + 1j * OM * q)
def da(g, a): return -(OM + 1j * (g + 1j * a) / ELL ** 2)
def dg(g, a): return -(g + 1j * a) / ELL ** 2 + 1j * OM
def r2(v): return np.concatenate([v.real, v.imag])
def proj_perp(B, rcond):
    U, S, _ = np.linalg.svd(B, full_matrices=False); U = U[:, S > rcond * S[0]]
    return lambda v: v - U @ (U.T @ v)
def setup(c, extra=None, W=15.0, central=None):
    z = zeros(); g = z[abs(z - c) <= W] - c
    if extra is not None: g = np.sort(np.append(g, extra - c))
    k0 = int(np.argmin(abs(g - (central - c)))) if central is not None else int(np.argmin(abs(g)))
    return g, k0
def pair(g0, d): return f(g0, A0 - d) + f(g0, A0 + d)
def tv(g0): ff = f(g0, A0); return ff * (da(g0, A0) ** 2 + 1 / ELL ** 2)
def nuis(g, k0, kappa, pairmode=True, dd=0.0):
    g0 = g[k0]; cols = []
    pc = pair(g0, dd)
    cols.append(r2(f(g0, A0 - dd) * dg(g0, A0 - dd) + f(g0, A0 + dd) * dg(g0, A0 + dd)))
    if kappa: cols += [r2(pc), r2(1j * pc)]
    for k, x in enumerate(g):
        if k == k0: continue
        ff = f(x, A0); cols += [r2(ff * dg(x, A0)), r2(ff * da(x, A0))]
    return np.array(cols).T
def tperp(g, k0, kappa, rcond):
    t = r2(tv(g[k0])); P = proj_perp(nuis(g, k0, kappa), rcond); return np.linalg.norm(P(t)), np.linalg.norm(t)
def single_tangent(g, k0, rcond):   # the single zero alternative: derivative in a of one zero (amplitude of the simple zero fixed)
    g0 = g[k0]; t = r2(f(g0, A0) * da(g0, A0)); cols = [r2(f(g0, A0) * dg(g0, A0))]
    for k, x in enumerate(g):
        if k == k0: continue
        ff = f(x, A0); cols += [r2(ff * dg(x, A0)), r2(ff * da(x, A0))]
    P = proj_perp(np.array(cols).T, rcond); return np.linalg.norm(P(t))
