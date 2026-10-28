# ht_residues.py: the residue weights, the total mass and the mean of the reciprocal Xi law mu, Section 1 of Chapter 24 (Nick Polson's hitting times chapter),
# and the cancellation factor of the density series. Claims supported (H = chapter_hitting_times.tex): H:51 w_{gamma_1} = 1.016e4, w_{gamma_2} = -1.18e6, |w| ~ 4e89 at the
# 120th zero, "with the first 120 zeros the series reproduces the total mass 1.0000000 and the mean sum gamma^-2 = 0.0231049931"; H:83 "near y = 0.008 the largest single
# term exceeds f(y) by a factor of about e^28".
# Conventions. F(s) = xi(1/2 + sqrt s)/xi(1/2), xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, w_gamma = 1/F'(-gamma^2) = 2 i gamma xi(1/2)/xi'(1/2 + i gamma),
# xi'(rho) = rho(rho-1) pi^(-rho/2) Gamma(rho/2) zeta'(rho)/2 at a zero rho. f(y) = sum_{gamma>0} w_gamma e^(-gamma^2 y) (RH, simple zeros, y > 0); the survival function is
# P(mu > y) = sum w_gamma gamma^-2 e^(-gamma^2 y) and int_y^inf t f(t) dt = sum w_gamma e^(-gamma^2 y)(y/gamma^2 + gamma^-4).
# Method. Zeros: data/zeros_1700.txt (68 digits). Mass and mean use y0 where the K-th term is below 1e-45 (y0 = (pi gamma_K/4 + 45)/gamma_K^2 ~ 0.0038 for K = 120), plus the
# mass below y0, which is below exp(-(log(1/y0))^2/(16 y0)) < 1e-100 (the small time law of H:270, an upper bound here only by order of magnitude); the sum is done for
# K = 100, 120, 150, the differences between K are the truncation error. The exact mean sum gamma^-2 = c_1 = (log xi)''(1/2)/2 by numerical differentiation at 60 digits, and
# (second check) the 1700 ordinates plus the smooth tail int_{gamma_N}^inf (log(t/2pi)/2pi) dt/t^2. The weight w_1, w_2 are also checked against the Hadamard product
# gamma^2 prod (1 - gamma^2/gamma'^2)^(-1) over the 1700 zeros with the tail exp(gamma^2 * smooth tail). Working precision dps and dps + 20, difference logged.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 ht_residues.py [--quick] [--out FILE] [--zeros FILE]; --quick uses K = 20, 30, 40 and 40 digits (about 2 s).
import argparse, json, time, os
import mpmath as mp

def load(path, n):
    return [mp.mpf(l.strip()) for l in open(path) if l.strip()][:n]

def weights(gam, dps):
    mp.mp.dps = dps
    xi0 = mp.mpf(1) / 2 * mp.mpf(-1) / 4 * mp.pi ** (-mp.mpf(1) / 4) * mp.gamma(mp.mpf(1) / 4) * mp.zeta(mp.mpf(1) / 2)
    ws = []
    for g in gam:
        rho = mp.mpc(0.5, g)
        xip = rho * (rho - 1) * mp.pi ** (-rho / 2) * mp.gamma(rho / 2) * mp.zeta(rho, derivative=1) / 2
        w = 2j * g * xi0 / xip
        ws.append(w)
    return ws

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
    ap.add_argument("--zeros", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "zeros_1700.txt"))
    a = ap.parse_args(); t0 = time.time(); R = {}; L = []
    Ks = (20, 30, 40) if a.quick else (100, 120, 150); dps = 40 if a.quick else 90
    mp.mp.dps = dps + 20
    gam = load(a.zeros, 1700)
    mp.mp.dps = dps
    K = max(Ks)
    w1 = weights(gam[:K], dps); w2 = weights(gam[:K], dps + 20)
    mp.mp.dps = dps
    L.append("max |Im w| / |w| = %s" % mp.nstr(max(abs(w.imag) / abs(w) for w in w1), 3))
    L.append("max relative difference of w between dps %d and %d: %s" % (dps, dps + 20, mp.nstr(max(abs(x - y) / abs(y) for x, y in zip(w1, w2)), 3)))
    w = [x.real for x in w1]
    R["w1"] = mp.nstr(w[0], 10); R["w2"] = mp.nstr(w[1], 10); R["w3"] = mp.nstr(w[2], 10)
    for n in (K,):
        R["w_%d" % n] = mp.nstr(w[n - 1], 6)
    if not a.quick:
        R["w_120"] = mp.nstr(w[119], 6); R["gamma_120"] = mp.nstr(gam[119], 12)
    # Hadamard product check for w_1, w_2
    N = len(gam); tailsum = lambda g: (mp.log(gam[-1] / (2 * mp.pi)) + 1) / (2 * mp.pi * gam[-1])
    chk = []
    for n in (0, 1):
        g2 = gam[n] ** 2; pr = g2
        for m in range(N):
            if m != n: pr /= (1 - g2 / gam[m] ** 2)
        pr *= mp.exp(g2 * tailsum(0))
        chk.append(pr)
    R["w1_product_with_tail"] = mp.nstr(chk[0], 8); R["w2_product_with_tail"] = mp.nstr(chk[1], 8)
    # exact c1
    xi = lambda s: s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2
    mp.mp.dps = 60
    c1 = mp.diff(lambda s: mp.log(xi(s)), mp.mpf(1) / 2, 2) / 2
    mp.mp.dps = dps
    S1 = sum(1 / g ** 2 for g in gam) + tailsum(0)
    R["c1_exact_logxi_second_derivative"] = mp.nstr(c1, 15); R["sum_gamma_minus2_1700_plus_tail"] = mp.nstr(S1, 10)
    # mass and mean for K in Ks
    for k in Ks:
        gk = gam[k - 1]; y0 = (mp.pi * gk / 4 + 45) / gk ** 2
        mass = sum(w[j] / gam[j] ** 2 * mp.exp(-gam[j] ** 2 * y0) for j in range(k))
        mean = sum(w[j] * mp.exp(-gam[j] ** 2 * y0) * (y0 / gam[j] ** 2 + 1 / gam[j] ** 4) for j in range(k))
        R["K%d_y0" % k] = mp.nstr(y0, 6); R["K%d_mass" % k] = mp.nstr(mass, 25); R["K%d_mean" % k] = mp.nstr(mean, 15)
        # same at a different y0 (1.3 y0) for the same K
        y1 = 1.3 * y0
        mass1 = sum(w[j] / gam[j] ** 2 * mp.exp(-gam[j] ** 2 * y1) for j in range(k))
        R["K%d_mass_at_1.3y0" % k] = mp.nstr(mass1, 25)
    # cancellation at y = 0.008
    y = mp.mpf("0.008"); kk = min(K, 120)
    terms = [w[j] * mp.exp(-gam[j] ** 2 * y) for j in range(kk)]
    fy = sum(terms); jm = max(range(kk), key=lambda j: abs(terms[j]))
    R["f(0.008)"] = mp.nstr(fy, 10); R["largest_term_0.008"] = mp.nstr(terms[jm], 10); R["largest_term_index_gamma"] = mp.nstr(gam[jm], 6)
    R["ratio_largest_term_over_f"] = mp.nstr(abs(terms[jm]) / fy, 6); R["log_ratio"] = mp.nstr(mp.log(abs(terms[jm]) / fy), 6)
    R["f(0.008)_with_K_minus_20"] = mp.nstr(sum(terms[:kk - 20]), 10)
    L.append("elapsed %.1f s" % (time.time() - t0))
    L += ["%s = %s" % (k, v) for k, v in R.items()]
    log = "\n".join(L); print(log)
    if a.out:
        json.dump(R, open(a.out + ".json", "w"), indent=1); open(a.out + ".log", "w").write(log + "\n")
main()
