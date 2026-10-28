# ht_jensen.py: the Taylor coefficients c_k of F(s) = xi(1/2 + sqrt s)/xi(1/2), the first order shifts A_n, the zeros of the Jensen polynomials J_d (tab:jensen), the Hermite
# identity and its asymptotics, and the height arithmetic of the Jensen section of Chapter 24 (Nick Polson's hitting times chapter, Section "The family" to "Assessment").
# Claims supported (H = chapter_hitting_times.tex): H:159 (600 coefficients c_k to 220 digits; A_1..A_6 = -8.69, -16.37, -15.15, -29.57, -15.83, -32.02), H:161 to H:174 (tab:jensen: relative
# shifts l_{d,n}/gamma_n^2 - 1 at d = 80, 640, 10240 against A_n/(2d)), H:140 (kappa_2 = c_1^2 - 2 c_2 (1 - 1/d)), H:196 to H:204 (Hermite identity J_d(-z^2) = E* C_d(zU) with
# C_d(x) = 1F1(-d; 1/2; x^2/(4d)) = (-1)^d d!/(2d)! H_{2d}(x/(2 sqrt d)) to twelve digits at five pairs (d, z), and the asymptotic C_d(x) ~ e^(x^2/(8d)) cos(x sqrt(1 + 1/(4d))) at x = 7, d = 10, 100, 1000
# with differences 5e-2, 3e-4, 3e-6), H:192 and H:215 to H:219 (arithmetic: 4.5e8, 1e25, D(T) = T^2/(2 Lambda_0) = 3e17 and 2e25).
# Conventions. xi(1/2 + s) = int_R Phi(v) e^(s v) dv with Phi(v) = 2 sum_{n>=1} (2 pi^2 n^4 e^(9v/2) - 3 pi n^2 e^(5v/2)) exp(-pi n^2 e^(2v)) (even, the Polya kernel in the variable v of
# the book, v = U); hence F(s) = sum c_k s^k with c_k = m_{2k}/((2k)! m_0), m_{2k} = int_R v^(2k) Phi(v) dv, and E* U^(2k) = m_{2k}/m_0. J_d(s) = sum_{k<=d} c_k prod_{j<k}(1 - j/d) s^k.
# A_n = s F''(s)/F'(s) at s = -gamma_n^2. Method. Moments by the trapezoid rule of step h on [0, 4] (Phi is analytic in |Im v| < pi/4; the integrand v^(2k) Phi has width ~ 0.03 at k = 600, so the
# error is about exp(-2 pi^2 (0.03/h)^2), below 1e-300 at h = 1/200); the coefficients are computed at two steps (full: h = 1/200 and 1/250, 260 digits; quick: 1/100 and 1/125, 80 digits) and the
# largest relative difference over k is logged (it is the error of every c_k). c_1 against (log xi)''(1/2)/2 (mpmath differentiation) and c_2 against (c_1^2 - sum gamma^-4)/2 with the 1700 ordinates of
# data/zeros_1700.txt plus the smooth tail. A_n also from -2 gamma_n^2 sum_{m != n} (gamma_m^2 - gamma_n^2)^-1 with the smooth tail int (log(t/2 pi)/2 pi) dt/(t^2 - gamma_n^2). The zeros of J_d by Newton
# (start from the first order shift) at two precisions and two truncations (K and K - 100 terms), and for d = 80 all 80 roots by polyroots (hyperbolic, the six smallest are the Newton roots).
# The Hermite identity: E* C_d(zU) by Gauss Legendre of mp.quad on [0, 4] (independent of the trapezoid moments), compared with J_d(-z^2) from c_k and with the Hermite formula of mp.hermite.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 ht_jensen.py [--quick] [--out FILE] [--zeros FILE] [--cout FILE]  (--cout writes the coefficients as strings to a json for ht_jlaws.py); --quick: K = 100 coefficients, 80 digits.
import argparse, json, os, time, sys
import mpmath as mp

def phi(v, dps):
    e2 = mp.exp(2 * v); a = mp.exp(9 * v / 2); b = mp.exp(5 * v / 2); s = mp.mpf(0); n = 1; eps = mp.mpf(10) ** (-dps - 20)
    while True:
        t = (2 * mp.pi ** 2 * n ** 4 * a - 3 * mp.pi * n ** 2 * b) * mp.exp(-mp.pi * n * n * e2)
        s += t
        if n > 2 and abs(t) < eps * (abs(s) + eps): break
        n += 1
    return 2 * s

def coeffs(K, hinv, dps, V=4):
    mp.mp.dps = dps
    h = mp.mpf(1) / hinv; J = V * hinv
    vs = [j * h for j in range(J + 1)]; ph = [phi(v, dps) for v in vs]
    w = [h * p for p in ph]; w[0] = w[0] / 2
    m0 = 2 * sum(w); c = [mp.mpf(1)]
    pw = [mp.mpf(1)] * (J + 1); v2 = [v * v for v in vs]; fact = mp.mpf(1)
    for k in range(1, K):
        pw = [p * q for p, q in zip(pw, v2)]; fact *= (2 * k - 1) * (2 * k)
        c.append(2 * sum(a * b for a, b in zip(w, pw)) / (fact * m0))
    return c, ph, m0

def jd_eval(c, d, s, K):
    # J_d(s) and derivative by Horner on the coefficients c_k prod_{j<k}(1-j/d), k <= min(d, K-1)
    top = min(d, K - 1); co = []; pr = mp.mpf(1)
    for k in range(top + 1):
        co.append(c[k] * pr); pr *= (1 - mp.mpf(k) / d)
    return co

def horner(co, s):
    f = mp.mpf(0) * s; g = mp.mpf(0) * s
    for a in reversed(co):
        g = g * s + f; f = f * s + a
    return f, g

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--cout", default=None)
    ap.add_argument("--zeros", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "zeros_1700.txt"))
    a = ap.parse_args(); t0 = time.time(); L = []; R = {}
    def log(s):
        L.append(s); print(s, flush=True)
    if a.quick: K, h1, h2, dps = 100, 100, 125, 80
    else: K, h1, h2, dps = 600, 200, 250, 260
    c, ph, m0 = coeffs(K, h1, dps); c2, _, _ = coeffs(K, h2, dps)
    mp.mp.dps = dps
    rel = [abs(x - y) / abs(y) for x, y in zip(c, c2)]
    R["coeff_max_rel_diff_between_steps"] = mp.nstr(max(rel), 3); R["coeff_rel_diff_at_k=K-1"] = mp.nstr(rel[-1], 3)
    log("coefficients: K=%d, steps 1/%d and 1/%d, %d digits, max rel diff %s (at last k %s), %.1f s" % (K, h1, h2, dps, R["coeff_max_rel_diff_between_steps"], R["coeff_rel_diff_at_k=K-1"], time.time() - t0))
    digits_agree = int(-mp.log10(max(rel))) if max(rel) > 0 else dps
    R["coeff_agreeing_digits"] = digits_agree
    # evenness check of phi
    R["phi_evenness"] = mp.nstr(abs(phi(mp.mpf("0.3"), 40) - phi(mp.mpf("-0.3"), 40)) / phi(mp.mpf("0.3"), 40), 3)
    R["c1"] = mp.nstr(c[1], 25); R["c2"] = mp.nstr(c[2], 25); R["c10"] = mp.nstr(c[10], 20); R["c_last"] = mp.nstr(c[K - 1], 20)
    mp.mp.dps = 60
    xi = lambda s: s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2
    c1ex = mp.diff(lambda s: mp.log(xi(s)), mp.mpf(1) / 2, 2) / 2
    R["c1_from_logxi"] = mp.nstr(c1ex, 20); R["c1_rel_diff"] = mp.nstr(abs(c1ex - c[1]) / c1ex, 3)
    mp.mp.dps = dps
    gam = [mp.mpf(l) for l in open(a.zeros)]
    mp.mp.dps = 40
    gam = [mp.mpf(g) for g in gam]
    N = len(gam); gN = gam[-1]
    tail2 = (lambda: mp.quad(lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) / t ** 4, [gN, mp.inf]))()
    p2 = sum(1 / g ** 4 for g in gam) + tail2
    c2z = (c[1] ** 2 - p2) / 2
    R["c2_from_zeros"] = mp.nstr(c2z, 12); R["c2_rel_diff"] = mp.nstr(abs(c2z - c[2]) / c[2], 3)
    log("c1 = %s (logxi %s), c2 = %s (zeros %s)" % (R["c1"][:16], R["c1_from_logxi"][:16], R["c2"][:14], R["c2_from_zeros"]))
    # A_n
    mp.mp.dps = dps
    co_full = c[:K]
    An = []
    for n in range(6):
        s = -mp.mpf(gam[n]) ** 2
        # derivative coefficients
        f1 = mp.mpf(0); f2 = mp.mpf(0)
        for k in range(K - 1, 0, -1):
            f1 = f1 * s + k * c[k]
        for k in range(K - 1, 1, -1):
            f2 = f2 * s + k * (k - 1) * c[k]
        An.append(s * f2 / f1)
    mp.mp.dps = 40
    An_z = []
    for n in range(6):
        g2 = gam[n] ** 2
        S = sum(1 / (gam[m] ** 2 - g2) for m in range(N) if m != n)
        S += mp.quad(lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) / (t * t - g2), [gN, 2 * gN, 10 * gN, mp.inf])
        An_z.append(-2 * g2 * S)
    R["A_series"] = [mp.nstr(x, 8) for x in An]; R["A_zeros"] = [mp.nstr(x, 8) for x in An_z]
    log("A_n (series) " + ", ".join(mp.nstr(x, 6) for x in An) + "; (zeros sum + tail) " + ", ".join(mp.nstr(x, 6) for x in An_z))
    # tab:jensen: zeros of J_d
    mp.mp.dps = dps
    tab = {}
    for d in (80, 640, 10240):
        row = []; row2 = []
        for n in range(6):
            g2 = mp.mpf(gam[n]) ** 2
            for (Kuse, pd) in ((K, dps), (K - K // 6, dps // 2)):
                mp.mp.dps = pd
                co = jd_eval(c, d, 0, Kuse); s = -g2 * (1 + An[n] / (2 * d))
                for it in range(60):
                    f, g = horner(co, s); ds = f / g; s -= ds
                    if abs(ds) < abs(s) * mp.mpf(10) ** (-(pd - 25)): break
                if Kuse == K: row.append(s / (-g2) - 1) if False else row.append(-s / g2 - 1)
                else: row2.append(-s / g2 - 1)
        mp.mp.dps = 40
        tab[d] = (row, row2)
        R["tab_d%d" % d] = [mp.nstr(x, 5) for x in row]
        R["tab_d%d_maxrelchange_truncK_and_precision" % d] = mp.nstr(max(abs(x - y) / abs(x) for x, y in zip(row, row2)), 3)
        log("d=%d: " % d + ", ".join(mp.nstr(x, 4) for x in row) + "  (max rel change with K-100 terms and half precision %s)" % R["tab_d%d_maxrelchange_truncK_and_precision" % d])
    R["A_over_2d_at_10240"] = [mp.nstr(x / (2 * 10240), 4) for x in An]
    log("A_n/(2d), d=10240: " + ", ".join(mp.nstr(x / 20480, 4) for x in An))
    # d = 80: hyperbolicity by sign changes of J_80 on a log grid of the negative axis (80 sign changes = 80 real simple roots)
    mp.mp.dps = 160
    co = jd_eval(c, 80, 0, K); co = [mp.mpf(x) for x in co]
    ng = 12000 if a.quick else 48000; prev = None; nchg = 0
    for i in range(ng + 1):
        sx = -mp.mpf(10) ** (mp.mpf(40) * i / ng)
        val = horner(co, sx)[0]
        if prev is not None and (val > 0) != (prev > 0): nchg += 1
        prev = val
    R["J80_sign_changes_on_(-1e40,-1)"] = nchg
    log("J_80: %d sign changes on the grid (-1e40, -1) of %d points (hyperbolic if 80)" % (nchg, ng + 1))
    # kappa_2 identity: variance of mu_d minus variance of mu = 2 c2/d
    mp.mp.dps = 40
    R["var_excess_2c2_over_d_at_d=100"] = mp.nstr(2 * c[2] / 100, 8)
    # Hermite identity
    mp.mp.dps = 40
    def Cd(d, x): return mp.hyp1f1(-d, mp.mpf(1) / 2, x * x / (4 * d))
    def CdH(d, x): return (-1) ** d * mp.factorial(d) / mp.factorial(2 * d) * mp.hermite(2 * d, x / (2 * mp.sqrt(d)))
    pairs = [(5, "0.7"), (10, "2.0"), (20, "3.5"), (40, "5.0"), (100, "7.0")]
    hm = []
    mp.mp.dps = 40
    m0q = 2 * mp.quad(lambda v: phi(v, 40), [0, 0.5, 1, 2, 3, 4])
    for d, z in pairs:
        z = mp.mpf(z)
        Jd = sum(c[k] * mp.fprod([(1 - mp.mpf(j) / d) for j in range(k)]) * (-z * z) ** k for k in range(min(d, K - 1) + 1))
        E = 2 * mp.quad(lambda v: phi(v, 40) * Cd(d, z * v), [0, 0.25, 0.5, 1, 1.5, 2, 3, 4]) / m0q
        hmv = max(abs(Cd(d, z * v) - CdH(d, z * v)) / (abs(Cd(d, z * v)) + 1e-30) for v in (mp.mpf("0.4"), mp.mpf("1.1")))
        hm.append((d, z, Jd, E, abs(Jd - E) / abs(Jd), hmv))
        log("Hermite (d=%d, z=%s): J_d(-z^2) = %s, E* C_d(zU) = %s, rel diff %s; 1F1 vs Hermite formula rel diff %s" % (d, z, mp.nstr(Jd, 14), mp.nstr(E, 14), mp.nstr(abs(Jd - E) / abs(Jd), 3), mp.nstr(hmv, 3)))
    R["hermite"] = [[d, str(z), mp.nstr(Jd, 14), mp.nstr(E, 14), mp.nstr(r, 3), mp.nstr(h, 3)] for d, z, Jd, E, r, h in hm]
    x = mp.mpf(7)
    for d in (10, 100, 1000):
        lhs = Cd(d, x); rhs = mp.exp(x * x / (8 * d)) * mp.cos(x * mp.sqrt(1 + mp.mpf(1) / (4 * d)))
        R["asym_d%d" % d] = {"C": mp.nstr(lhs, 10), "approx": mp.nstr(rhs, 10), "abs_diff": mp.nstr(abs(lhs - rhs), 4), "x3_over_96d2": mp.nstr(x ** 3 / (96 * d * d), 4), "x4_over_d2_times_rhs": mp.nstr(abs(rhs) * x ** 4 / (d * d), 4)}
        log("asymptotic x=7, d=%d: C_d = %s, e^(x^2/8d) cos(.) = %s, |diff| = %s (x^3/(96 d^2) = %s)" % (d, mp.nstr(lhs, 8), mp.nstr(rhs, 8), mp.nstr(abs(lhs - rhs), 3), mp.nstr(x ** 3 / (96 * d * d), 3)))
    # arithmetic
    R["sqrt_2e17"] = mp.nstr(mp.sqrt(2 * mp.mpf(10) ** 17), 6)
    R["D_T=5.45e8_L0=0.5"] = mp.nstr((mp.mpf("5.45e8")) ** 2 / (2 * mp.mpf("0.5")), 5)
    R["D_T=3e12_L0=0.2"] = mp.nstr((mp.mpf("3e12")) ** 2 / (2 * mp.mpf("0.2")), 5)
    R["gamma2_over_4delta2_DH"] = mp.nstr(mp.mpf("85.6993") ** 2 / (4 * mp.mpf("0.3085") ** 2), 6)
    R["s_rho_DH"] = [mp.nstr(mp.mpc("0.3085", "85.6993") ** 2, 8)]; R["abs_s_rho_sq"] = mp.nstr(abs(mp.mpc("0.3085", "85.6993") ** 2) ** 2, 4)
    R["d*_t_c=0.1765"] = mp.nstr(mp.mpf("85.699") ** 2 / (2 * mp.mpf("0.1765")), 6)
    log("arithmetic: sqrt(2e17) = %s; D(5.45e8, 1/2) = %s; D(3e12, 0.2) = %s; 85.6993^2/(4 0.3085^2) = %s; |s_rho|^2 = %s; d*(t_c = 0.1765) = %s" % (R["sqrt_2e17"], R["D_T=5.45e8_L0=0.5"], R["D_T=3e12_L0=0.2"], R["gamma2_over_4delta2_DH"], R["abs_s_rho_sq"], R["d*_t_c=0.1765"]))
    log("elapsed %.1f s" % (time.time() - t0))
    if a.cout:
        json.dump({"K": K, "dps": dps, "c": [mp.nstr(x, dps - 5) for x in c]}, open(a.cout, "w"))
    if a.out:
        json.dump(R, open(a.out + ".json", "w"), indent=1); open(a.out + ".log", "w").write("\n".join(L) + "\n")
main()
