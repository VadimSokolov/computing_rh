# ht_kent.py: the Brownian exit check of Kent's formula in Section "Kent's formula" of Chapter 24 (Nick Polson's hitting times chapter, H:67 to H:79).
# Claims supported (H = chapter_hitting_times.tex): H:67 (exit time of standard Brownian motion from (-1,1): Laplace transform 1/cosh sqrt(2s),
# eigenvalues (k+1/2)^2 pi^2/2, residues (-1)^k pi (k+1/2)), mass 1, mean 1, transform at s = 1 equal to 1/cosh(sqrt 2) = 0.459098131085425,
# and the table of the truncated Kent products with N = 50 and N = 400 terms against the residues (k = 0, 1, 5), H:69 to H:79.
# Conventions. lambda_m = (m + 1/2)^2 pi^2 / 2 for m = 0, 1, ...; F(s) = cosh(sqrt(2 s)) = prod_m (1 + s/lambda_m); Kent weight of eigenvalue n with N terms
# W_n^N = lambda_n prod_{m != n, m < N} lambda_m/(lambda_m - lambda_n) (m = 0..N-1, n < N), compared with 1/F'(-lambda_n) = (-1)^n pi (n + 1/2).
# Method. Products at 40 digits (second run at 60 digits, difference logged). The density f(t) = sum_k (-1)^k pi (k+1/2) exp(-(k+1/2)^2 pi^2 t/2) is also
# computed from its dual (small time) form f(t) = 2 (2 pi t^3)^(-1/2) sum_{m >= 0} (-1)^m (2 m + 1) exp(-(2m+1)^2/(2 t)) (the first passage density of Brownian motion from the
# centre of (-1,1)); the two series are compared at t = 0.05 .. 1, the mass, mean and transform are integrated with the theta form on t >= 0.5 and the dual form on t < 0.5.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 ht_kent.py [--quick] [--out FILE]; --quick = 25 digits, runs in about 1 s; the full run takes a few seconds.
import argparse, json, time
import mpmath as mp

def run(dps, N_list=(50, 400)):
    mp.mp.dps = dps
    lam = lambda m: (mp.mpf(m) + mp.mpf(1) / 2) ** 2 * mp.pi ** 2 / 2
    out = {}
    for N in N_list:
        L = [lam(m) for m in range(N)]
        for k in (0, 1, 5):
            w = L[k]
            for m in range(N):
                if m != k: w *= L[m] / (L[m] - L[k])
            out["W_N%d_k%d" % (N, k)] = w
    for k in (0, 1, 5):
        out["res_k%d" % k] = (-1) ** k * mp.pi * (k + mp.mpf(1) / 2)
    theta = lambda t: mp.nsum(lambda k: (-1) ** int(k) * mp.pi * (k + mp.mpf(1) / 2) * mp.exp(-(k + mp.mpf(1) / 2) ** 2 * mp.pi ** 2 * t / 2), [0, mp.inf]) if False else \
        sum((-1) ** k * mp.pi * (k + mp.mpf(1) / 2) * mp.exp(-(k + mp.mpf(1) / 2) ** 2 * mp.pi ** 2 * t / 2) for k in range(int(mp.sqrt(2 * dps * 2.31 / t) / mp.pi) + 5))
    dual = lambda t: 2 * (2 * mp.pi * t ** 3) ** -0.5 * sum((-1) ** m * (2 * m + 1) * mp.exp(-(2 * m + 1) ** 2 / (2 * t)) for m in range(int(mp.sqrt(2 * t * dps * 2.31)) + 6))
    f = lambda t: theta(t) if t >= 0.5 else dual(t)
    out["series_dual_maxdiff"] = max(abs(theta(mp.mpf(t)) - dual(mp.mpf(t))) for t in (0.1, 0.2, 0.3, 0.5, 0.8, 1.0)) if True else 0
    # integrals: on [0.5, T] with theta form, the tail beyond T is the first term only (exact), on (0, 0.5] with the dual form
    def integ(g):
        a = mp.quad(lambda t: g(t) * dual(t), [0, 0.05, 0.1, 0.25, 0.5])
        b = mp.quad(lambda t: g(t) * theta(t), [0.5, 1, 2, 4, 8, 16])
        # tail t > 16: theta series with terms k: int_16^inf g(t) e^{-l t}: g = 1 -> e^{-16 l}/l; g = t -> e^{-16 l}(16/l + 1/l^2); g = e^{-t} -> e^{-17 l'}/(l+1) form
        return a + b
    out["mass_noTail"] = integ(lambda t: 1)
    tail = lambda T, sh=0, p=0: sum((-1) ** k * mp.pi * (k + mp.mpf(1) / 2) * mp.exp(-(lam(k) + sh) * T) * (1 / (lam(k) + sh) if p == 0 else (T / (lam(k) + sh) + 1 / (lam(k) + sh) ** 2)) for k in range(6))
    out["mass"] = out["mass_noTail"] + tail(16)
    out["mean"] = integ(lambda t: t) + tail(16, 0, 1)
    out["laplace_s1"] = integ(lambda t: mp.exp(-t)) + tail(16, 1, 0)
    out["laplace_s1_exact"] = 1 / mp.cosh(mp.sqrt(2))
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
    a = ap.parse_args(); t0 = time.time()
    d1, d2 = (20, 30) if a.quick else (40, 60)
    r1, r2 = run(d1), run(d2)
    res = {}
    for k in r1:
        res[k] = mp.nstr(r2[k], 20)
    diff = max(abs(r1[k] - r2[k]) / max(abs(r2[k]), 1) for k in r1)
    res["max_rel_diff_between_precisions"] = mp.nstr(diff, 5)
    lines = ["ht_kent: precisions %d and %d digits, %.1f s" % (d1, d2, time.time() - t0)]
    for k, v in res.items(): lines.append("%s = %s" % (k, v))
    lines.append("W_N50: " + ", ".join(mp.nstr(r2["W_N50_k%d" % k], 5) for k in (0, 1, 5)) + "; W_N400: " + ", ".join(mp.nstr(r2["W_N400_k%d" % k], 5) for k in (0, 1, 5)) + "; residues: " + ", ".join(mp.nstr(r2["res_k%d" % k], 5) for k in (0, 1, 5)))
    log = "\n".join(lines); print(log)
    if a.out:
        json.dump(res, open(a.out + ".json", "w"), indent=1); open(a.out + ".log", "w").write(log + "\n")
