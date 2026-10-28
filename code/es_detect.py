# es_detect.py: the detection experiment of the subsection "A ReLU network" (explicit sparse chapter), pencil and plug-in matched filter rows (the two ReLU rows are not
# reproduced: they need training). Window model of es_window.py: level 0.6 (a0 = 0.1), taper ell = 3, J = 100 samples omega_j = omega0 + j Delta (omega0, Delta are not given in the
# chapter: defaults omega0 = 5, Delta = 0.404, i.e. Omega = 45), CUE ordinates at the density at height 200 (N = 80), zeros with |x| <= W = 10 in the model, the zero nearest the
# centre either a double zero on the line (2 f(a0)) or the pair f(a0 - delta) + f(a0 + delta). Noise: complex Gaussian with tau/sqrt2 per real coordinate (the chapter's convention).
# Detectors, each calibrated to 5 per cent false alarms on --nnull null windows: pencil statistic T = max |a_hat - a0| over the pencil components within 0.5 of the central ordinate
# (T = 1 if none), model order = number of singular values of H0 above 10 times the median of the second half of the spectrum (Hankel (J - L) x L, L = J/2); plug-in matched filter
# <t_perp, F>/||t_perp|| with t the second derivative of the central zero's f at the pencil's ordinate and a0, t_perp orthogonal to the tangent columns (amplitude and pole, complex) of
# every pencil component (this is my reading of "built from the pencil's components", the chapter does not define it further). delta* = smallest delta with power 1/2, log-log interpolation.
# Error control: two seeds. Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import es_window as W0
ELL = 3.0; J = 100; L = 50; WIN = 10.0; A0 = 0.1
def run_pencil(F, dom):
    Hk = np.array([F[i:i + L] for i in range(J - L + 1)]); H0 = Hk[:-1]; H1 = Hk[1:]
    U, S, Vh = np.linalg.svd(H0); fl = np.median(S[L // 2:]); M = max(1, int(np.sum(S > 10 * fl)))
    Um = U[:, :M]; Vm = Vh[:M].conj().T
    z = np.linalg.eigvals(Um.conj().T @ H1 @ Vm / S[:M]); return z, M
def window(rng, om):
    xs = W0.cue_ordinates(rng); xs = xs[abs(xs) <= WIN]; k0 = int(np.argmin(abs(xs))); return xs, k0
def data(om, xs, k0, delta, ell=ELL):
    F = np.zeros(len(om), complex)
    for k, x in enumerate(xs):
        F += W0.f(om, x, A0, ell) if k != k0 else (W0.f(om, x, A0 - delta, ell) + W0.f(om, x, A0 + delta, ell) if delta > 0 else 2 * W0.f(om, x, A0, ell))
    return F
def stats(F, om, dom, g0):
    z, M = run_pencil(F, dom); gam = np.angle(z) / dom; a = -np.log(abs(z)) / dom
    near = abs(gam - g0) <= 0.5; T = float(np.max(abs(a[near] - A0))) if near.any() else 1.0
    V = np.exp(np.outer(np.arange(J), np.log(z.astype(complex)))); nrm = np.linalg.norm(V, axis=0)
    cols = []
    for m in range(len(z)):
        v = V[:, m] / nrm[m]; w = np.arange(J) * v; cols += [W0.real2(v), W0.real2(1j * v), W0.real2(w), W0.real2(1j * w)]
    ic = np.argmin(abs(gam - g0)) if len(z) else 0
    t = W0.real2(W0.tvec(om, gam[ic], ELL))
    B = np.array(cols).T; U, S, _ = np.linalg.svd(B, full_matrices=False); U = U[:, S > 1e-10 * S[0]]
    tp = t - U @ (U.T @ t); n = np.linalg.norm(tp)
    P = float(tp @ W0.real2(F) / n) if n > 1e-300 else 0.0
    return T, P, M
def main():
    import argparse; ap = argparse.ArgumentParser(); ap.add_argument('--tau', type=float); ap.add_argument('--nnull', type=int, default=1200); ap.add_argument('--npow', type=int, default=150)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--om0', type=float, default=5.0); ap.add_argument('--Om', type=float, default=45.0); ap.add_argument('--out'); a = ap.parse_args()
    dom = (a.Om - a.om0) / (J - 1); om = a.om0 + dom * np.arange(J); rng = np.random.default_rng(a.seed); tau = a.tau; sd = tau / np.sqrt(2)
    def noisy(F): return F + sd * (rng.standard_normal(J) + 1j * rng.standard_normal(J))
    nullT = []; nullP = []; Ms = []
    for _ in range(a.nnull):
        xs, k0 = window(rng, om); T, P, M = stats(noisy(data(om, xs, k0, 0.0)), om, dom, xs[k0]); nullT.append(T); nullP.append(P); Ms.append(M)
    thT = np.quantile(nullT, 0.95); thP = np.quantile(nullP, 0.95)
    deltas = np.logspace(-7.5, np.log10(0.09), 16); powT = []; powP = []
    for d in deltas:
        cT = cP = 0
        for _ in range(a.npow):
            xs, k0 = window(rng, om); T, P, M = stats(noisy(data(om, xs, k0, d)), om, dom, xs[k0]); cT += T > thT; cP += P > thP
        powT.append(cT / a.npow); powP.append(cP / a.npow)
    def dstar(p):
        p = np.array(p)
        for i in range(len(p) - 1):
            if p[i] < 0.5 <= p[i + 1]:
                lo, hi = np.log(deltas[i]), np.log(deltas[i + 1]); return float(np.exp(lo + (0.5 - p[i]) / (p[i + 1] - p[i]) * (hi - lo)))
        return None
    r = dict(tau=tau, om0=a.om0, Om=a.Om, seed=a.seed, median_order_null=float(np.median(Ms)), thr_pencil=float(thT), dstar_pencil=dstar(powT), dstar_plugin=dstar(powP), power_pencil=powT, power_plugin=powP, deltas=deltas.tolist())
    print(json.dumps(r)); 
    if a.out: json.dump(r, open(a.out, 'w'))
main()
