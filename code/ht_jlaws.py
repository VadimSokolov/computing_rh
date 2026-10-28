# ht_jlaws.py: the distance between the Jensen laws mu_d and mu, Section "Rate of convergence" of Chapter 24 (Nick Polson's hitting times chapter), claim H:176 ("Inverting 1/(s J_d(s)) by Talbot's
# method and comparing with the residue expansion of mu on a grid of 18 points in [0.007, 0.08], the largest difference of distribution functions is 0.155, 0.0548, 0.0154, 0.00397, 0.00100 for
# d = 10, 40, 160, 640, 2560, so d sup|F_d - F| tends to about 2.57; the maximum sits at y ~ 0.016"). The grid of 18 points is not defined in the chapter: two grids are used, 18 equispaced points
# y_i = 0.007 + i 0.073/17 and the fine grid y = 0.007 + 0.001 i, i = 0..73 (the largest difference on the fine grid is the better estimate of a supremum).
# Conventions. F_d(y) = P(mu_d <= y), 1/J_d the Laplace transform of mu_d (J_d(s) = sum_{k<=d} c_k prod_{j<k}(1 - j/d) s^k, c_k = m_{2k}/((2k)! m_0) from the Polya kernel by the trapezoid rule, as in ht_jensen.py);
# F(y) = P(mu <= y) = 1 - sum_{gamma>0} w_gamma gamma^-2 e^(-gamma^2 y) with w_gamma = 2 i gamma xi(1/2)/xi'(1/2 + i gamma) (the residue series of Section 1; 200 ordinates of data/zeros_1700.txt, the terms for
# gamma > 200 are below 1e-100 at y >= 0.007).
# Method. F_d by the fixed Talbot rule of Abate and Valko with M nodes (r = 2M/(5y)) applied to 1/(s J_d(s)); J_d is evaluated by Horner on the exact polynomial (all d + 1 coefficients, so no truncation)
# at 400 digits; the 2 coefficient sets (trapezoid steps 1/500 and 1/600) are compared (the largest difference of F_d over the grids is logged), and M = 40 against M = 50 (the largest difference of F_d on the
# 18 point grid is logged); the series for F is checked against 150 instead of 200 ordinates. Output: for each d the largest |F_d - F| on both grids, its argmax, d times it, and F_d, F on the 18 point grid.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 ht_jlaws.py [--quick] [--out FILE] [--zeros FILE] [--cpus N]; quick: d = 10, 40, 160, step 1/200 and 1/250, 120 digits, M = 30 and 40, fine grid thinned (smoke test).
import argparse, json, os, time
from multiprocessing import Pool
import mpmath as mp

def phi(v, dps):
    e2 = mp.exp(2 * v); a = mp.exp(9 * v / 2); b = mp.exp(5 * v / 2); s = mp.mpf(0); n = 1; eps = mp.mpf(10) ** (-dps - 20)
    while True:
        t = (2 * mp.pi ** 2 * n ** 4 * a - 3 * mp.pi * n ** 2 * b) * mp.exp(-mp.pi * n * n * e2); s += t
        if n > 2 and abs(t) < eps * (abs(s) + eps): break
        n += 1
    return 2 * s

def coeffs(K, hinv, dps, V=5):
    mp.mp.dps = dps; h = mp.mpf(1) / hinv; J = V * hinv
    vs = [j * h for j in range(J + 1)]; w = [h * phi(v, dps) for v in vs]; w[0] /= 2; m0 = 2 * sum(w)
    c = [mp.mpf(1)]; pw = [mp.mpf(1)] * (J + 1); v2 = [v * v for v in vs]; fact = mp.mpf(1)
    for k in range(1, K):
        pw = [p * q for p, q in zip(pw, v2)]; fact *= (2 * k - 1) * (2 * k)
        c.append(2 * sum(x * y for x, y in zip(w, pw)) / (fact * m0))
    return c

G = {}

def talbot(d, y, M, which):
    c = G["c%d" % which]; mp.mp.dps = G["dps"]
    co = []; pr = mp.mpf(1)
    for k in range(d + 1): co.append(c[k] * pr); pr *= (1 - mp.mpf(k) / d)
    co.reverse()
    def Ld(s):
        f = mp.mpf(0)
        for a in co: f = f * s + a
        return 1 / (s * f)
    y = mp.mpf(y); r = 2 * M / (5 * y); tot = Ld(mp.mpf(r)) * mp.exp(r * y) / 2
    for k in range(1, M):
        th = k * mp.pi / M; cot = mp.cot(th); s = r * th * mp.mpc(cot, 1); sg = th + (th * cot - 1) * cot
        tot += mp.re(mp.exp(y * s) * mp.mpc(1, sg) * Ld(s))
    return r / M * tot

def job(args):
    return args, talbot(*args)

def F_zeros(ys, gam, nz, dps=90):
    mp.mp.dps = dps
    xi0 = -mp.pi ** (-mp.mpf(1) / 4) * mp.gamma(mp.mpf(1) / 4) * mp.zeta(mp.mpf(1) / 2) / 8
    W = []
    for g in gam[:nz]:
        rho = mp.mpc(0.5, g); xip = rho * (rho - 1) * mp.pi ** (-rho / 2) * mp.gamma(rho / 2) * mp.zeta(rho, derivative=1) / 2
        W.append((2j * g * xi0 / xip).real)
    return [1 - sum(w / g ** 2 * mp.exp(-g ** 2 * y) for w, g in zip(W, gam[:nz])) for y in ys]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--zeros", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "zeros_1700.txt"))
    a = ap.parse_args(); t0 = time.time(); L = []; R = {}
    def log(s): L.append(s); print(s, flush=True)
    if a.quick: ds, dps, hs, Ms, thin = (10, 40, 160), 120, (200, 250), (30, 40), 6
    else: ds, dps, hs, Ms, thin = (10, 40, 160, 640, 2560), 400, (500, 600), (40, 50), 1
    K = max(ds) + 1
    mp.mp.dps = dps + 20; gam = [mp.mpf(l) for l in open(a.zeros)][:250]; mp.mp.dps = dps
    G["dps"] = dps
    G["c0"] = coeffs(K, hs[0], dps); log("coefficients 1/%d done %.1f s" % (hs[0], time.time() - t0))
    G["c1"] = coeffs(K, hs[1], dps); log("coefficients 1/%d done %.1f s" % (hs[1], time.time() - t0))
    mp.mp.dps = dps
    R["coeff_max_rel_diff_steps"] = mp.nstr(max(abs(x - y) / abs(y) for x, y in zip(G["c0"], G["c1"])), 3)
    log("max relative difference of the coefficients between the two steps: " + R["coeff_max_rel_diff_steps"])
    g18 = [mp.mpf("0.007") + i * mp.mpf("0.073") / 17 for i in range(18)]
    gfine = [mp.mpf("0.007") + mp.mpf("0.001") * i for i in range(0, 74, thin)]
    ys = g18 + gfine
    Fz = F_zeros(ys, gam, 200); Fz150 = F_zeros(ys, gam, 150)
    mp.mp.dps = 40
    R["F_zero_series_diff_200_vs_150"] = mp.nstr(max(abs(x - y) for x, y in zip(Fz, Fz150)), 3)
    log("F(y) from 200 vs 150 ordinates: max diff %s; F at 18 grid: %s" % (R["F_zero_series_diff_200_vs_150"], ", ".join(mp.nstr(x, 4) for x in Fz[:18])))
    tasks = [(d, y, Ms[0], 0) for d in ds for y in ys]
    tasks2 = [(d, y, Ms[1], 0) for d in ds[:3] for y in g18] + [(d, y, Ms[0], 1) for d in ds[:3] for y in g18]
    with Pool(a.cpus) as P:
        res = dict(P.map(job, tasks, chunksize=1)); log("main tasks done %.1f s" % (time.time() - t0))
        res2 = dict(P.map(job, tasks2, chunksize=1))
    mp.mp.dps = 40
    R["talbot_M_diff_(d<=160,18pt)"] = mp.nstr(max(abs(res[(d, y, Ms[0], 0)] - res2[(d, y, Ms[1], 0)]) for d in ds[:3] for y in g18), 3)
    R["step_diff_(d<=160,18pt)"] = mp.nstr(max(abs(res[(d, y, Ms[0], 0)] - res2[(d, y, Ms[0], 1)]) for d in ds[:3] for y in g18), 3)
    log("Talbot M=%d vs %d max diff %s; coefficient step 1/%d vs 1/%d max diff %s (d <= 160, 18 point grid)" % (Ms[0], Ms[1], R["talbot_M_diff_(d<=160,18pt)"], hs[0], hs[1], R["step_diff_(d<=160,18pt)"]))
    for d in ds:
        e18 = [abs(res[(d, y, Ms[0], 0)] - Fz[i]) for i, y in enumerate(g18)]
        ef = [abs(res[(d, y, Ms[0], 0)] - Fz[18 + i]) for i, y in enumerate(gfine)]
        i18 = max(range(18), key=lambda i: e18[i]); jf = max(range(len(gfine)), key=lambda i: ef[i])
        R["d%d" % d] = {"sup18": mp.nstr(e18[i18], 5), "argmax18": mp.nstr(g18[i18], 4), "supfine": mp.nstr(ef[jf], 5), "argmaxfine": mp.nstr(gfine[jf], 4), "d_sup18": mp.nstr(d * e18[i18], 5), "d_supfine": mp.nstr(d * ef[jf], 5)}
        log("d=%d: sup(18 pt) = %s at y=%s (d*sup = %s); sup(fine) = %s at y=%s (d*sup = %s)" % (d, R["d%d" % d]["sup18"], R["d%d" % d]["argmax18"], R["d%d" % d]["d_sup18"], R["d%d" % d]["supfine"], R["d%d" % d]["argmaxfine"], R["d%d" % d]["d_supfine"]))
    R["F_d_on_18pt_grid"] = {str(d): [mp.nstr(res[(d, y, Ms[0], 0)], 8) for y in g18] for d in ds}; R["F_on_18pt_grid"] = [mp.nstr(x, 8) for x in Fz[:18]]
    log("elapsed %.1f s" % (time.time() - t0))
    if a.out:
        json.dump(R, open(a.out + ".json", "w"), indent=1); open(a.out + ".log", "w").write("\n".join(L) + "\n")
if __name__ == "__main__":
    import multiprocessing; multiprocessing.set_start_method("fork"); main()
