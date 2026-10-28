# bsy_negw.py: every number of Remark bsy:rem:negW (a negative trace) of Chapter 35 (ch/bsy1.tex, "a negative trace"), group BSY-C.
# Claims supported: for B = L(s, chi_4), Lambda_B(s) = (4/pi)^(s/2) Gamma((s+1)/2) B(s), odd primes p < q < 2p, a = log(pq)/2, b = log(p/q)/2, P(s) = 1 + 2 sqrt p p^-s + 2 sqrt q q^-s + sqrt(pq) (pq)^-s:
#   (i)   the identity (4pq/pi)^(s/2) Gamma((s+1)/2) P(s) B(s) = 2 (pq)^(1/4) Lambda_B(s) [cosh(a (s - 1/2)) + 2 cosh(b (s - 1/2))]  (checked at random complex s, and its consequences:
#         real and positive on R, invariant under s -> 1 - s);
#   (ii)  P B is a real linear combination of B, B_p, B_q, B_pq (L functions modulo 4, 4p, 4q, 4pq), B_m = B prod_(l | m) (1 - chi_4(l) l^-s)  (checked numerically, with the coefficients printed);
#   (iii) H_eps(z) = cosh z + 2 cosh(eps z), z = a (s - 1/2), eps = b/a: its zeros have |Re z| < 3, exactly two zeros with 0 <= Im z <= 2 pi (argument principle on the rectangle), near
#         +-arcosh 2 + i pi; for p = 101, q = 103 the pair 1/2 +- delta_0 + i gamma_0 with delta_0 = Re z/a, gamma_0 = Im z/a, which lies below the first zero of Lambda_B;
#   (iv)  the trace K(u) = sum over the zeros of the completed function with Im > 0 of e^(-gamma u) cos(delta u) (zeros of H_eps with m <= 5 and the first five zeros of L(s, chi_4)),
#         K(u) = 2 e^(-gamma_0 u) cos(delta_0 u) + ..., and K(pi/delta_0) < 0; printed values: delta_0 = 0.28475, gamma_0 = 0.67928, K(pi/delta_0) = -1.11e-3 (and W(11.03) = -1.1e-3).
# Conventions. The trace of the zeros rho = 1/2 + delta + i gamma of the completed function G(s) (zeros of Lambda_B and of H_eps) over those with gamma > 0, as in bsy_levy.py (a zero on the
# line contributes e^(-gamma u), a symmetric pair 1/2 +- delta + i gamma contributes 2 e^(-gamma u) cos(delta u)). "m <= 5" is read as the zeros of H_eps near +-arcosh 2 + i pi (2k + 1),
# k = 0 .. m (six pairs for m = 5), and the sum is also logged for k = 0..4 and for 60 pairs. The zeros of L(s, chi_4) are found from the sign changes of the real function Lambda_B(1/2 + it) on a grid of step
# 0.02 on (0, 19) and refined by findroot; their number up to height 18.3 and the position of the sixth zero are logged. All at mp.dps = 30, repeated at dps = 50 (differences logged).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Output: data/bsy_negw.json and the log data/bsy_negw.log. Usage: python3 -I bsy_negw.py [--quick] [--out FILE]; --quick: dps 20 and a coarser grid, to misc/reviews/nick-b4/quick/BSY-C/.
import argparse, json, math, random
from pathlib import Path
import mpmath as mp

ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
here = Path(__file__).resolve().parent
out = Path(args.out) if args.out else (Path("/tmp/bsy_negw_quick.json") if args.quick else here.parent / "data" / "bsy_negw.json")
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)

p, q = 101, 103
chi = lambda n: [0, 1, 0, -1][n % 4]


def run(dps, step, npairs):
    mp.mp.dps = dps
    res = {}
    B = lambda s: mp.dirichlet(s, [0, 1, 0, -1])
    LamB = lambda s: (4 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * B(s)
    a = mp.log(p * q) / 2; b = mp.log(mp.mpf(p) / q) / 2; eps = b / a
    Pf = lambda s: 1 + 2 * mp.sqrt(p) * mp.mpf(p) ** (-s) + 2 * mp.sqrt(q) * mp.mpf(q) ** (-s) + mp.sqrt(p * q) * mp.mpf(p * q) ** (-s)
    H = lambda z: mp.cosh(z) + 2 * mp.cosh(eps * z)
    # (i) identity
    random.seed(7); worst = 0
    for _ in range(4):
        s = mp.mpc(random.uniform(-1, 2), random.uniform(-6, 6))
        lhs = (4 * p * q / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * Pf(s) * B(s)
        rhs = 2 * (p * q) ** mp.mpf("0.25") * LamB(s) * (mp.cosh(a * (s - mp.mpf(1) / 2)) + 2 * mp.cosh(b * (s - mp.mpf(1) / 2)))
        worst = max(worst, abs(lhs - rhs) / abs(rhs))
    res["identity_max_rel"] = float(worst)
    # real positive on R, invariance
    xs = [mp.mpf(x) / 4 for x in range(-12, 17)]
    vals = [(4 * p * q / mp.pi) ** (x / 2) * mp.gamma((x + 1) / 2) * Pf(x) * B(x) for x in xs if x not in (-1, -3, -5)]
    res["real_axis_min_positive"] = bool(all(v.real > 0 and abs(v.imag) < mp.mpf(10) ** (-dps + 8) * abs(v) for v in vals))
    G = lambda s: (4 * p * q / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * Pf(s) * B(s)
    s0 = mp.mpc("0.3", "1.7")
    res["functional_eq_rel"] = float(abs(G(s0) - G(1 - s0)) / abs(G(s0)))
    # (ii) combination
    s1 = mp.mpc("0.7", "2.3")
    Bp = lambda s: B(s) * (1 - chi(p) * mp.mpf(p) ** (-s)); Bq = lambda s: B(s) * (1 - chi(q) * mp.mpf(q) ** (-s)); Bpq = lambda s: B(s) * (1 - chi(p) * mp.mpf(p) ** (-s)) * (1 - chi(q) * mp.mpf(q) ** (-s))
    # P B = B + 2 sqrt p chi(p)(B - Bp) + 2 sqrt q chi(q)(B - Bq) + sqrt(pq) chi(pq) (B - Bp - Bq + Bpq)
    cp, cq, cpq = chi(p), chi(q), chi(p) * chi(q)
    c0 = 1 + 2 * mp.sqrt(p) * cp + 2 * mp.sqrt(q) * cq + mp.sqrt(p * q) * cpq
    c1 = -(2 * mp.sqrt(p) * cp + mp.sqrt(p * q) * cpq); c2 = -(2 * mp.sqrt(q) * cq + mp.sqrt(p * q) * cpq); c3 = mp.sqrt(p * q) * cpq
    comb = c0 * B(s1) + c1 * Bp(s1) + c2 * Bq(s1) + c3 * Bpq(s1)
    res["combination_rel"] = float(abs(comb - Pf(s1) * B(s1)) / abs(comb))
    res["combination_coeffs"] = [float(c) for c in (c0, c1, c2, c3)]
    # (iii) zeros of H_eps
    zs = []
    for k in range(npairs):
        for sg in (1, -1):
            z0 = sg * mp.acosh(2) + 1j * mp.pi * (2 * k + 1)
            z = mp.findroot(H, z0, tol=mp.mpf(10) ** (-dps + 5))
            zs.append(z)
    res["H_zero_resid_max"] = float(max(abs(H(z)) for z in zs))
    res["H_max_abs_re"] = float(max(abs(z.real) for z in zs))
    # argument principle on [-3,3] x [0, 2 pi] and on [-3, 3] x [0, 12 pi]
    def count(top):
        pts = [mp.mpc(-3, 0), mp.mpc(3, 0), mp.mpc(3, top), mp.mpc(-3, top), mp.mpc(-3, 0)]
        tot = 0
        for u, v in zip(pts[:-1], pts[1:]):
            n = 400 if top > 7 else 120
            prev = H(u)
            for j in range(1, n + 1):
                cur = H(u + (v - u) * j / n)
                tot += mp.arg(cur / prev); prev = cur
        return float(tot / (2 * mp.pi))
    res["count_rect_2pi"] = count(2 * mp.pi)
    res["count_rect_12pi"] = count(12 * mp.pi)
    first = min((z for z in zs if z.real > 0), key=lambda z: z.imag)
    res["eps"] = float(eps); res["a"] = float(a); res["acosh2"] = float(mp.acosh(2))
    res["z0"] = [float(first.real), float(first.imag)]
    delta0 = first.real / a; gamma0 = first.imag / a
    res["delta0"] = float(delta0); res["gamma0"] = float(gamma0)
    # zeros of Lambda_B on the line
    Lreal = lambda t: LamB(mp.mpc(0.5, t))
    ts = [mp.mpf(j) * step for j in range(1, int(19 / step) + 1)]
    lz = []; prev = Lreal(ts[0]); res["Lambda_B_max_imag_rel"] = 0.0
    for t in ts[1:]:
        cur = Lreal(t)
        res["Lambda_B_max_imag_rel"] = max(res["Lambda_B_max_imag_rel"], float(abs(cur.imag) / abs(cur)))
        if cur.real * prev.real < 0:
            lz.append(mp.findroot(lambda x: Lreal(x).real, (t - step, t), solver="anderson", tol=mp.mpf(10) ** (-dps + 5)))
        prev = cur
    res["Lambda_B_zeros_below_19"] = [float(x) for x in lz]
    res["first_zero_Lambda_B_above_gamma0"] = bool(gamma0 < lz[0])
    # trace
    def K(u, npr):
        u = mp.mpf(u)
        sK = sum(mp.exp(-z.imag / a * u) * mp.cos(z.real / a * u) for z in zs if z.imag > 0 and z in zs[: 2 * npr])
        return sK + sum(mp.exp(-g * u) for g in lz[:5])
    u0 = mp.pi / delta0
    res["u0"] = float(u0)
    res["K_u0_m5"] = float(K(u0, 6)); res["K_u0_k4"] = float(K(u0, 5)); res["K_u0_all"] = float(K(u0, npairs)); res["K_u0_first_pair"] = float(2 * mp.exp(-gamma0 * u0) * mp.cos(delta0 * u0))
    res["K_11.03"] = float(K(mp.mpf("11.03"), 6))
    # scan of K on (0, 40]
    us = [mp.mpf(j) / 4 for j in range(1, 161)]
    Ks = [K(u, 6) for u in us]
    jm = min(range(len(us)), key=lambda j: Ks[j])
    res["K_scan_min"] = [float(us[jm]), float(Ks[jm])]
    res["K_scan_sign_changes"] = sum(1 for j in range(len(us) - 1) if Ks[j] * Ks[j + 1] < 0)
    res["K_at_odd_multiples"] = {str(2 * k + 1): float(K((2 * k + 1) * mp.pi / delta0, 6)) for k in range(0, 4)}
    return res


r1 = run(20 if args.quick else 30, mp.mpf("0.05") if args.quick else mp.mpf("0.02"), 8 if args.quick else 60)
for k, v in r1.items():
    P(k, v)
if not args.quick:
    r2 = run(50, mp.mpf("0.02"), 12)
    for k in ("delta0", "gamma0", "K_u0_m5", "K_11.03"):
        P("difference dps 30 vs 50:", k, r1[k] - r2[k])
    for k in ("count_rect_2pi", "count_rect_12pi"):
        P("dps 50", k, r2[k])
else:
    r2 = None
P("printed in the remark: delta_0 = 0.28475, gamma_0 = 0.67928, K(pi/delta_0) = -1.11e-3; computed:", f"delta_0 = {r1['delta0']:.5f}, gamma_0 = {r1['gamma0']:.5f}, K(pi/delta_0) = {r1['K_u0_m5']:.3e}")
out.write_text(json.dumps({"p": p, "q": q, "dps30": r1, "dps50": r2}, indent=1))
Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
