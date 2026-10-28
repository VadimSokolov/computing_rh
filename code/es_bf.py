# es_bf.py: profile likelihood of the pair model near 60 and Bayes factors for the real zero 201.264752 (explicit sparse chapter, "Pairs, not zeros" and "Bayes factors").
# Uses es_pair.py (level 1 window, ell = 10, 400 samples; model, nuisance and linearization described there). Claims: profile flat below delta ~ 1e-5, detected from 3e-5 at
# tau = 1.7e-11, t_perp = 1.41, delta* = 4.4e-6; log10 BF01 = 9.1 (single zero alternative) and 3.9 (pair alternative, pi uniform on (0,1/2)) for 201.264752, 6.9 to 9.3 for the ten
# other central zeros; the formula of Theorem "Bayes factors" with s = 1.7e11 gives 5.3 at Z = 0 and 4.5 at Z = 1.7. Output JSON and a log to stdout.
# Error control: ||t_perp|| for relative SVD cutoffs 1e-8..1e-14 and nuisance windows W = 10, 15, 25 (es_pair.py output); BF grids of 300 and 600 points are compared.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import es_pair as E
from scipy.integrate import quad
from scipy.stats import norm
TAU = 1.7e-11
def G(Z): return quad(lambda u: u ** -0.5 * np.exp(Z * u - u * u / 2), 0, 1, limit=200)[0] + quad(lambda u: u ** -0.5 * np.exp(Z * u - u * u / 2), 1, 60, limit=200)[0]
def formula_log10(s, Z, pi0=2.0): return np.log10(2 * np.sqrt(s) / (pi0 * G(Z)))
def Ps(g, k0, kappa, grid, rcond=1e-10):
    out = []
    for d in grid:
        B = E.nuis(g, k0, kappa, dd=d); U, S, _ = np.linalg.svd(B, full_matrices=False); out.append(U[:, S > rcond * S[0]])
    return out
def bf_pair(r0, noises, g, k0, kappa, grid, scale, tau):
    """r(delta) = r0(0) - scale*(pair(delta)) ... returns log10 BF01 per noise draw (pi uniform on (0,1/2), density 2) and the profile chi2 of the first draw"""
    Us = Ps(g, k0, kappa, grid); g0 = g[k0]; pc0 = E.r2(E.pair(g0, 0.0))
    res = []; prof = None
    for n in noises:
        chi = np.empty(len(grid))
        for i, d in enumerate(grid):
            r = r0 - scale * (E.r2(E.pair(g0, d)) - pc0) + n; U = Us[i]; r = r - U @ (U.T @ r); chi[i] = r @ r / tau ** 2
        if prof is None: prof = chi.copy()
        L = np.exp(-(chi - chi[0]) / 2)   # L(delta)/L(0)... chi[0] is the grid's first point, ~ delta = 0
        integ = 2 * (np.sum(0.5 * (L[1:] + L[:-1]) * np.diff(grid)))
        res.append(-np.log10(integ)); 
    return np.array(res), prof
def main():
    out = {}
    grid = np.concatenate([[0.0], np.logspace(-9, np.log10(0.5), int(sys.argv[2]) if len(sys.argv) > 2 else 300)])
    rng = np.random.default_rng(7)
    # window about 60: pair added at 60.089 (synthetic, fixed amplitudes), data = model + noise
    g, k0 = E.setup(60.0, 60.089, 15.0, 60.089); print("window 60: zeros in window", len(g), "central", g[k0] + 60)
    for kappa in (0, 1):
        tp = E.tperp(g, k0, kappa, 1e-10)[0]; print("kappa", kappa, "||t_perp|| = %.4g" % tp, "delta* (z=1.645, tau=1.7e-11, tau^2 I) = %.3g" % np.sqrt(1.645 * TAU / tp),
              " with tau/sqrt2: %.3g" % np.sqrt(1.645 * TAU / np.sqrt(2) / tp))
    rows = []
    for kappa in (0, 1):
        for dt in (0.0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
            nz = rng.standard_normal(2 * len(E.OM)) * TAU
            r0 = E.r2(E.pair(g[k0], dt) - E.pair(g[k0], 0.0))      # data minus model at delta = 0 (central part only), truth delta = dt
            Us = Ps(g, k0, kappa, grid[:1]); U = Us[0]
            chi_grid = []
            for d in grid[1:]:
                U1 = Ps(g, k0, kappa, [d])[0]; r = r0 + nz - (E.r2(E.pair(g[k0], d)) - E.r2(E.pair(g[k0], 0.0))); r = r - U1 @ (U1.T @ r); chi_grid.append(r @ r / TAU ** 2)
            chi_grid = np.array(chi_grid); r = r0 + nz; r = r - U @ (U.T @ r); chi0 = r @ r / TAU ** 2
            i = int(np.argmin(chi_grid)); lg = -(chi_grid - chi0) / 2; m = lg.max(); L = np.exp(lg - m)
            integ = 2 * (np.sum(0.5 * (L[1:] + L[:-1]) * np.diff(grid[1:])) + L[0] * grid[1])
            lb = -(np.log10(integ) + m / np.log(10))
            rows.append(dict(kappa=kappa, delta_true=dt, chi2_at_0=chi0, chi2_min=float(chi_grid[i]), delta_hat=float(grid[1:][i]), dchi2=chi0 - float(chi_grid[i]), log10BF01=lb))
            print("profile kappa", kappa, "delta_true %.0e" % dt, "chi2(0) %.3g min %.3g at %.3g  dchi2 %.3g log10BF01 %.2f" % (chi0, chi_grid[i], grid[1:][i], chi0 - chi_grid[i], lb))
    out["profile60"] = rows
    # real zero 201.264752: simple zero (data = f(a0) at the central zero), pair alternative with the model 1/2 pair(delta) and kappa free (2 extra amplitude parameters)
    g, k0 = E.setup(200.0, None, 15.0, 201.264752); print("window 200: zeros", len(g), "central", g[k0] + 200)
    tp1 = E.tperp(g, k0, 1, 1e-10)[0]; tp0 = E.tperp(g, k0, 0, 1e-10)[0]; s1 = E.single_tangent(g, k0, 1e-10)
    print("||t_perp|| kappa free %.4g, kappa fixed %.4g, single zero ||P d_a f|| %.4g" % (tp1, tp0, s1))
    noises = [rng.standard_normal(2 * len(E.OM)) * TAU for _ in range(40)]
    zero = np.zeros(2 * len(E.OM))
    r0 = zero
    for scale, lab in ((0.5, "pair model 1/2 pair(delta), kappa free"), (1.0, "pair model pair(delta) unhalved (chapter's t), kappa free")):
        lg, prof = bf_pair(r0, [zero] + noises, g, k0, 1, grid, scale, TAU)
        print(lab, "log10BF01 Asimov (Z=0) %.2f, noise draws median %.2f range %.2f..%.2f" % (lg[0], np.median(lg[1:]), lg[1:].min(), lg[1:].max()))
        out["pair201_" + str(scale)] = dict(asimov=float(lg[0]), median=float(np.median(lg[1:])), lo=float(lg[1:].min()), hi=float(lg[1:].max()))
    # single zero alternative (width a0 + eps, eps uniform (0, 1/2)): linear Gaussian closed form, Z1 standard normal under the null
    for tau, lab in ((TAU, "tau=1.7e-11 per real coordinate"),):
        zs = np.array([0.0, 1.0, 1.7, -1.0]); bf = s1 / (2 * tau * np.sqrt(2 * np.pi) * np.exp(zs ** 2 / 2) * norm.cdf(zs))
        print("single zero alternative", lab, "Z1", zs, "log10BF01", np.round(np.log10(bf), 2))
        out["single201"] = dict(Z=zs.tolist(), log10BF01=np.log10(bf).tolist(), s1=s1)
    # ten other central zeros: the single zero alternative at Z1 = 0 (Asimov) and the median over Z1 ~ N(0,1) (under the null)
    z = E.zeros(); cen = z[np.argsort(abs(z - 200))[:11]]; vals = []
    for c0 in sorted(cen):
        g, k0 = E.setup(200.0, None, 15.0, c0); s1c = E.single_tangent(g, k0, 1e-10); zz = rng.standard_normal(4001)
        bfd = s1c / (2 * TAU * np.sqrt(2 * np.pi) * np.exp(zz ** 2 / 2) * norm.cdf(zz)); vals.append((c0, s1c, np.log10(s1c / (TAU * 2.5066282746310002)), float(np.median(np.log10(bfd)))))
        print("central zero %.6f  ||P d_a f|| %.4g  log10BF01(Z1=0) %.2f  median over Z1 %.2f" % vals[-1])
    out["central11"] = vals
    # the formula of Theorem bf as printed
    for s in (1.7e11, 2 * 1.7e11, 8.3e10, tp1 / TAU, 2.816 / TAU):
        print("formula s=%.3g: log10BF01 Z=0 %.2f, Z=1.7 %.2f" % (s, formula_log10(s, 0.0), formula_log10(s, 1.7)))
    print("G(0)=%.5f exact 2^(-3/4) Gamma(1/4)=%.5f, G(1.7)=%.4f" % (G(0.0), 2 ** -0.75 * 3.625609908221908, G(1.7)))
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
main()
