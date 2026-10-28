# ht_ev.py: the extreme value and prime side numbers of Chapter 24 (Nick Polson's hitting times chapter): the weights b_k and the Gram midpoints (H:263, H:266), the smoothed zero laws of tab:smooth
# (H:272 to H:290), the Toeplitz minors of the Polya density (H:293), the table of log E exp(-s Z_P) and rho_4, rho_5 (H:318 to H:337). Claims: b_1..b_5 = 0.0050, 0.0045, 0.0048, 0.0043, 0.0046; t^_1,t^_2,t^_3 = 14.52,
# 20.65, 25.49; tab:smooth (N = 0, 1, 2, 5, 10, 30, 100: sup|F_N - F| = 1.2e-2, 3.2e-3, 2.6e-3, 8.5e-4, 2.5e-5, 2.7e-6, 6.8e-7 and mean errors -2.2e-4, 4.6e-5, -3.6e-5, 1.1e-5, -3.2e-7, 3.4e-8, -8.4e-9);
# "the smoothed law alone is within 1.2% of mu ... a hundred times better than the Jensen law of degree 2560" (1.00e-3 from H:176); the Toeplitz determinants det[p((i-j)h)]_{i,j<k} of p = Phi(|u|)/Phi(0): nonnegative for k <= 4,
# -4.0e-11 and -1.4e-7 at k = 5, h = 0.03 and 0.05 (90 digits), positive from h = 0.1; the table of log E e^(-s Z_P) for P = 1e3..1e7, s = 0.5, 1, 2 (H:321), the limits 0.1667, 0.4714, 1.3333, rho_4 = 0.54, 0.37, 0.26,
# 0.19, 0.13, the expansion within 2e-2 at P = 1e3 and 2e-4 at 1e7 and 0.44986 against 0.44983 at P = 1e7, s = 1 (H:337), and the asymptotic 2^(-1/3) P^(-1/6) (log P)^(1/3) of the review (0.48, 0.36, 0.26, 0.19, 0.14).
# Conventions. Gram midpoints: theta(t^_k) = (k - 3/2) pi, theta the Riemann Siegel theta (mpmath siegeltheta for k <= 200, the Stirling series of theta with Newton for larger k); b_k = k/gamma_k^2.
# mu^[N]: sum E_k/lambda_k with lambda_k = gamma_k^2 for k <= N and t^_k^2 for k > N, k up to 4e5 (the first 4e5 smoothed points; quick: 2e4) plus the integral int_{t_L}^inf log(1 + s/t^2) theta'(t)/pi dt, theta(t_L) = (4e5 - 1) pi (the
# midpoint rule in k), theta'/pi = (log(t/2 pi))/(2 pi) + O(t^-2); F_N by the fixed Talbot rule of Abate and Valko (double precision, M = 18, 20, 22 (20 reported); the difference of the sup is logged) applied to exp(-sum log(1 + s/lambda_k))/s; F by the
# residue series with 200 ordinates (mpmath, 90 digits, as in ht_residues.py); two grids as in ht_jlaws.py (18 equispaced points on [0.007, 0.08] and the fine grid step 0.001); the mean error is E mu^[N] - c_1 with c_1 the exact 0.02310499311541...
# computed from the sum 1/gamma^2 plus tail inside ht_jensen.py's value (here from the (log xi)'' second derivative). Toeplitz: Phi(u) = sum_n (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u)) exp(-pi n^2 e^(4u)), dets by LU in mpmath at 90 and 130 digits
# (difference logged). Euler: primes to 1e7 (numpy sieve), q = p^(-1/2), K_P(w) = sum_p [log((1 - q)/(1 - q p^(-w))) + w log p q/(1 - q)], log E e^(-sZ) = s kappa_2/c - K_P(sqrt(2 s/c)), c = kappa_3^(2/3), kappa_j = sum_p (log p)^j
# sum_m m^(j-1) q^m (Eulerian closed forms for j = 2..5), in long double; rho_j = kappa_j/kappa_3^(j/3); expansion (2s)^(3/2)/6 - rho_4 (2s)^2/24 + rho_5 (2s)^(5/2)/120.
# Written from the mathematics (the collaborator's scripts were not read, run or copied). Usage: python3 ht_ev.py [--quick] [--out FILE] [--zeros FILE] [--part smooth,pf,euler]
import argparse, json, os, time
import numpy as np
import mpmath as mp
from scipy import integrate

def theta_np(t): return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)
def dtheta_np(t): return 0.5 * np.log(t / (2 * np.pi)) - 1 / (48 * t ** 2) - 7 / (1920 * t ** 4)

def gram_mid(Mn, nexact=200):
    mp.mp.dps = 30
    out = np.zeros(Mn)
    for k in range(1, min(nexact, Mn) + 1):
        out[k - 1] = float(mp.findroot(lambda t: mp.siegeltheta(t) - (k - mp.mpf(3) / 2) * mp.pi, 2 * mp.pi * k / mp.log(k + 3) + 10))
    if Mn > nexact:
        ks = np.arange(nexact + 1, Mn + 1, dtype=float); t = np.maximum(out[nexact - 1] + (ks - nexact) * 2 * np.pi / np.log(out[nexact - 1] / (2 * np.pi)), 30.0)
        for _ in range(60): t = t - (theta_np(t) - (ks - 1.5) * np.pi) / dtheta_np(t)
        out[nexact:] = t
    return out

def talbot_nodes(y, M):
    r = 2 * M / (5 * y); th = np.arange(1, M) * np.pi / M; cot = 1 / np.tan(th)
    s = r * th * (cot + 1j); sg = th + (th * cot - 1) * cot
    return r, s, sg

def part_smooth(a, log, R):
    quick = a.quick; Mn = 20000 if quick else 400000
    mp.mp.dps = 90
    gam = [mp.mpf(l) for l in open(a.zeros)][:200]
    g2 = np.array([float(g) ** 2 for g in gam[:100]])
    # b_k and midpoints
    tt = gram_mid(Mn); R["b_1_5"] = [mp.nstr(k / gam[k - 1] ** 2, 4) for k in range(1, 6)]
    R["tgram_1_3"] = [float(x) for x in tt[:3]]
    log("b_1..5 = %s; Gram midpoints t^_1..3 = %s; t^_%d = %.3f" % (", ".join(R["b_1_5"]), ", ".join("%.4f" % x for x in tt[:3]), Mn, tt[-1]))
    # F by residue series
    xi0 = -mp.pi ** (-mp.mpf(1) / 4) * mp.gamma(mp.mpf(1) / 4) * mp.zeta(mp.mpf(1) / 2) / 8
    W = []
    for g in gam:
        rho = mp.mpc(0.5, g); xip = rho * (rho - 1) * mp.pi ** (-rho / 2) * mp.gamma(rho / 2) * mp.zeta(rho, derivative=1) / 2
        W.append((2j * g * xi0 / xip).real)
    g18 = [0.007 + i * 0.073 / 17 for i in range(18)]; gfine = [0.007 + 0.001 * i for i in range(0, 74, 6 if quick else 1)]
    ys = g18 + gfine
    Fz = np.array([float(1 - sum(w / g ** 2 * mp.exp(-g ** 2 * mp.mpf(y)) for w, g in zip(W, gam))) for y in ys])
    xi = lambda s: s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2
    mp.mp.dps = 40; c1 = float(mp.diff(lambda s: mp.log(xi(s)), mp.mpf(1) / 2, 2) / 2)
    gall = np.array([float(l) for l in open(a.zeros)]); gam2all = gall ** 2
    NS = [0, 1, 2, 5, 10, 30, 100]
    t2 = tt ** 2; tL = float(tt[-1]) + 0.0
    # t_L: theta(t_L) = (Mn - 1) pi
    tl = tt[-1]
    for _ in range(100): tl = tl - (theta_np(tl) - (Mn - 1) * np.pi) / dtheta_np(tl)
    # mean tail: sum_{k>Mn} t^_k^-2 ~ int_{tl}^inf t^-2 dtheta/pi
    tailmean = (np.log(tl / (2 * np.pi)) + 1) / (2 * np.pi * tl)
    means = {N: float(np.sum(1 / gam2all[:N]) if N else 0) + float(np.sum(1 / t2[N:])) + tailmean for N in NS}
    R["mean_error"] = {N: means[N] - c1 for N in NS}
    log("mean errors (E mu^[N] - c_1): " + ", ".join("N=%d: %.3e" % (N, R["mean_error"][N]) for N in NS))
    out = {}
    for M in (18, 20, 22):
        FN = {N: np.zeros(len(ys)) for N in NS}
        for iy, y in enumerate(ys):
            r, s, sg = talbot_nodes(y, M)
            nodes = np.concatenate([[r + 0j], s])
            # all-smoothed log transform
            LA = np.zeros(len(nodes), dtype=complex); B = {N: np.zeros(len(nodes), dtype=complex) for N in NS}
            for i, sv in enumerate(nodes):
                lg = np.log1p(sv / t2)
                if abs(sv) < 0.05 * tl ** 2:
                    lt = np.log(tl / (2 * np.pi)); tail = 0j
                    for n_ in range(1, 6):
                        m_ = 2 * n_; tail += (-1) ** (n_ + 1) * sv ** n_ / n_ * tl ** (1 - m_) * (lt / (m_ - 1) + 1 / (m_ - 1) ** 2) / (2 * np.pi)
                else:
                    tail = integrate.quad(lambda t: (np.log1p(sv / t ** 2)).real * dtheta_np(t) / np.pi, tl, np.inf, limit=300)[0] + 1j * integrate.quad(lambda t: (np.log1p(sv / t ** 2)).imag * dtheta_np(t) / np.pi, tl, np.inf, limit=300)[0]
                LA[i] = lg.sum() + tail
                cs_s = np.cumsum(lg[:100]); cs_g = np.cumsum(np.log1p(sv / g2))
                for N in NS:
                    B[N][i] = (cs_s[N - 1] - cs_g[N - 1]) if N else 0
            for N in NS:
                Gs = np.exp(-(LA - B[N])) / nodes
                tot = 0.5 * Gs[0].real * np.exp(r * y) + np.sum((np.exp(y * s) * (1 + 1j * sg) * Gs[1:]).real)
                FN[N][iy] = r / M * tot
        out[M] = FN
        log("Talbot M=%d done" % M)
    for N in NS:
        e18 = np.abs(out[20][N][:18] - Fz[:18]); ef = np.abs(out[20][N][18:] - Fz[18:])
        ya = np.array(ys); mk = ya >= 0.011   # Talbot in double is unreliable where F < 1e-3 (y < 0.011): masked sup is the result
        d_all = np.max(np.abs(out[22][N] - out[20][N])); d_mk = np.max(np.abs(out[22][N] - out[20][N])[mk])
        e18m = e18[mk[:18]]; efm = ef[mk[18:]]
        R["sup_N%d" % N] = {"sup18": float(e18.max()), "supfine": float(ef.max()), "sup18_y>=0.011": float(e18m.max()), "supfine_y>=0.011": float(efm.max()), "argmaxfine_y>=0.011": float(ya[18:][mk[18:]][int(efm.argmax())]), "M22_vs_M20_all": float(d_all), "M22_vs_M20_y>=0.011": float(d_mk)}
        log("N=%d: sup|F_N - F| 18pt %.3e, fine %.3e; for y >= 0.011: 18pt %.3e, fine %.3e (at y=%.3f); mean error %.3e; Talbot M22 vs M20 %.1e (all), %.1e (y >= 0.011)" % (N, e18.max(), ef.max(), e18m.max(), efm.max(), ya[18:][mk[18:]][int(efm.argmax())], R["mean_error"][N], d_all, d_mk))

def Phi_dbn(u):
    e4 = mp.exp(4 * u); a = mp.exp(9 * u); b = mp.exp(5 * u); s = mp.mpf(0); n = 1
    while True:
        t = (2 * mp.pi ** 2 * n ** 4 * a - 3 * mp.pi * n ** 2 * b) * mp.exp(-mp.pi * n * n * e4); s += t
        if n > 2 and abs(t) < mp.mpf(10) ** (-mp.mp.dps - 30) * abs(s): break
        n += 1
    return s

def part_pf(a, log, R):
    hs = ["0.01", "0.02", "0.03", "0.05", "0.07", "0.1", "0.15", "0.2", "0.3", "0.5", "1"]
    res = {}
    for dps in (90, 130):
        mp.mp.dps = dps
        for h in hs:
            hh = mp.mpf(h); p0 = Phi_dbn(mp.mpf(0)); vals = {m: Phi_dbn(m * hh) / p0 for m in range(0, 7)}
            for k in (3, 4, 5, 6):
                Mx = mp.matrix(k, k)
                for i in range(k):
                    for j in range(k): Mx[i, j] = vals[abs(i - j)]
                res[(dps, h, k)] = mp.det(Mx)
    R["toeplitz"] = {}
    for h in hs:
        row = {}
        for k in (3, 4, 5, 6):
            d1 = res[(90, h, k)]; d2 = res[(130, h, k)]
            row[k] = (mp.nstr(d2, 6), mp.nstr(abs(d1 - d2) / abs(d2), 2))
        R["toeplitz"][h] = row
        log("h=%s: " % h + "; ".join("k=%d: %s (rel diff dps 90/130 %s)" % (k, row[k][0], row[k][1]) for k in (3, 4, 5, 6)))

def primes_upto(n):
    sv = np.ones(n + 1, dtype=bool); sv[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]: sv[i * i::i] = False
    return np.nonzero(sv)[0]

def part_euler(a, log, R):
    Pm = 100000 if a.quick else 10 ** 7
    p = primes_upto(Pm).astype(np.longdouble); lp = np.log(p); q = p ** (-np.longdouble(0.5))
    e = {2: q / (1 - q) ** 2, 3: q * (1 + q) / (1 - q) ** 3, 4: q * (1 + 4 * q + q * q) / (1 - q) ** 4, 5: q * (1 + 11 * q + 11 * q ** 2 + q ** 3) / (1 - q) ** 5}
    # check the closed forms against direct sums for p=2
    qq = np.longdouble(2) ** (-np.longdouble(0.5)); m = np.arange(1, 400, dtype=np.longdouble)
    R["eulerian_check"] = [float(abs(np.sum(m ** (j - 1) * qq ** m) - e[j][0]) / e[j][0]) for j in (2, 3, 4, 5)]
    ss = [0.5, 1.0, 2.0]; table = {}
    Ps = [10 ** j for j in range(3, 8) if 10 ** j <= Pm]
    for P in Ps:
        idx = int(np.searchsorted(p, P, side="right")); pl = p[:idx]; lpl = lp[:idx]; ql = q[:idx]
        kap = {j: float(np.sum(lpl ** j * e[j][:idx])) for j in (2, 3, 4, 5)}
        c = kap[3] ** (2 / 3); rho4 = kap[4] / kap[3] ** (4 / 3); rho5 = kap[5] / kap[3] ** (5 / 3)
        row = {}
        for s in ss:
            w = np.longdouble(np.sqrt(2 * s / c))
            K = np.sum(np.log((1 - ql) / (1 - ql * pl ** (-w))) + w * lpl * ql / (1 - ql))
            val = float(s * kap[2] / c - K)
            lim = (2 * s) ** 1.5 / 6; ex = lim - rho4 * (2 * s) ** 2 / 24 + rho5 * (2 * s) ** 2.5 / 120
            row[s] = (val, lim, ex)
        table[P] = (row, rho4, rho5, float(2 ** (-1 / 3) * P ** (-1 / 6) * np.log(P) ** (1 / 3)), kap)
        log("P=1e%d: log E e^{-sZ} = %s; limit %s; expansion %s; rho_4 = %.4f (asymptotic %.4f), rho_5 = %.4f; kappa_2 %.4e kappa_3 %.4e" % (int(round(np.log10(P))), ", ".join("%.5f" % row[s][0] for s in ss), ", ".join("%.4f" % row[s][1] for s in ss), ", ".join("%.5f" % row[s][2] for s in ss), rho4, table[P][3], rho5, kap[2], kap[3]))
    R["euler"] = {str(P): {"values": [table[P][0][s][0] for s in ss], "expansion": [table[P][0][s][2] for s in ss], "rho4": table[P][1], "rho5": table[P][2], "rho4_asympt": table[P][3]} for P in Ps}
    R["euler_max_expansion_error"] = {str(P): max(abs(table[P][0][s][0] - table[P][0][s][2]) for s in ss) for P in Ps}
    log("max |exact - expansion| over s: " + ", ".join("P=%g: %.2e" % (P, v) for P, v in ((float(k), v) for k, v in R["euler_max_expansion_error"].items())))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--part", default="smooth,pf,euler")
    ap.add_argument("--zeros", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "zeros_1700.txt"))
    a = ap.parse_args(); t0 = time.time(); L = []; R = {}
    def log(s): L.append(s); print(s, flush=True)
    for part in a.part.split(","):
        {"smooth": part_smooth, "pf": part_pf, "euler": part_euler}[part](a, log, R); log("part %s done %.1f s" % (part, time.time() - t0))
    if a.out:
        json.dump(R, open(a.out + ".json", "w"), indent=1, default=str); open(a.out + ".log", "w").write("\n".join(L) + "\n")
if __name__ == "__main__":
    main()
