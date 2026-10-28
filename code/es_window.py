# es_window.py: the window model of Theorem "tapered transform" of the explicit sparse chapter at level theta = 0.6, taper ell, J samples, and the
# quantities of Theorem "the rate for the real part": the tangent t of the pair at 1/2 +- delta (Lemma "evenness"), its component t_perp orthogonal to the
# tangent vectors of the other parameters, the bound delta* = (z_alpha tau/||t_perp||)^(1/2). Supports the claims of the subsection "A ReLU network"
# (median ||t_perp|| = 318, bound 0.060 tau^(1/2)). Conventions: f(a)_j = exp(-q^2/(2 ell^2) + i omega_j q), q = (gamma - c) + i a, omega_j = omega0 + j Delta
# (omega0, Delta are NOT defined in the chapter's description of the experiment: they are parameters here), vectors in R^(2J) (real and imaginary parts),
# ordinates = eigenangles of a Haar unitary matrix (CUE), scaled to mean spacing 1/d with d = log(200/(2 pi))/(2 pi), centre at 0.
# Nuisance parameters: (gamma, a) of every other zero within +-W of the centre, and gamma of the central pair (the central pair's mean width is fixed by the
# functional equation). Written from the mathematics (the collaborator's scripts were not read, run or copied).
import numpy as np
THETA = 0.6; A0 = THETA - 0.5
D200 = np.log(200 / (2 * np.pi)) / (2 * np.pi)

def cue_ordinates(rng, N=80):
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z); Q = Q * (np.diag(R) / abs(np.diag(R)))
    th = np.sort(np.angle(np.linalg.eigvals(Q)))
    x = th * N / (2 * np.pi) / D200           # mean spacing 1/d
    x = x - rng.uniform(-1, 1) * 0.0
    return x

def omegas(J, om0, dom): return om0 + dom * np.arange(J)

def f(om, g, a, ell):          # g = gamma - c
    q = g + 1j * a
    return np.exp(-q ** 2 / (2 * ell ** 2) + 1j * om * q)

def dlogf_a(om, g, a, ell): return -(om + 1j * (g + 1j * a) / ell ** 2)
def dlogf_g(om, g, a, ell): return -(g + 1j * a) / ell ** 2 + 1j * om

def real2(v): return np.concatenate([v.real, v.imag], axis=0)

def tvec(om, g0, ell):
    ff = f(om, g0, A0, ell)
    return ff * (dlogf_a(om, g0, A0, ell) ** 2 + 1 / ell ** 2)

def tperp(om, xs, ell, W):
    k0 = np.argmin(abs(xs)); g0 = xs[k0]
    t = real2(tvec(om, g0, ell))
    cols = [real2(f(om, g0, A0, ell) * 2 * dlogf_g(om, g0, A0, ell))]
    for k, x in enumerate(xs):
        if k == k0 or abs(x) > W: continue
        ff = f(om, x, A0, ell)
        cols += [real2(ff * dlogf_g(om, x, A0, ell)), real2(ff * dlogf_a(om, x, A0, ell))]
    B = np.array(cols).T
    Q, _ = np.linalg.qr(B)
    tp = t - Q @ (Q.T @ t)
    return np.linalg.norm(tp), np.linalg.norm(t), B.shape[1]
