# cx_slope_X.py: Arb enclosures for the end of Proposition prop:slope of the chapter "complexity" (Nick Polson's notes, source2/lvjbook/
# chapter_complexity.tex), task C3 of misc/reviews/nick-b4/reconcile-I.md; written from the mathematics (the collaborator's scripts were not read,
# run or copied). Claims: X'(17.4978) = -1.04e-5 and X'(17.4980) = +7.25e-6 (floating, under RH), and inf X in [0.183897, 0.183907] with X at the
# minimiser 0.18389874 (floating), where, under the Riemann hypothesis (all zeros 1/2 + i gamma, gamma real),
#   X(v) = sum over all zeros rho of (1/(v-gamma)^2 - 1/gamma^2) = sum_{gamma>0} f_v(gamma), f_v(g) = 1/(g-v)^2 + 1/(g+v)^2 - 2/g^2 (the chapter's definition).
# Hence X'(v) = -2 sum_{gamma>0} [(v-gamma)^-3 + (v+gamma)^-3] and X''(v) = 6 sum [(v-gamma)^-4 + (v+gamma)^-4] > 0, so X is convex on every interval free of zeros.
# Method. X = Q + F, Q the sum over the first K = 1000 ordinates (certified Arb enclosures, acb.zeta_zeros at 320 bits, K+1 zeros, T0 = 1419.9 lies
# between gamma_1000 and gamma_1001, checked, and N(T0) = K) and F the sum over gamma > T0, 0 <= F, F' (all terms >= 0 for v < T0).
# Tail bound (rigorous): f_v(g) <= 6 v^2/(g-v)^4 and 2[(g-v)^-3 - (g+v)^-3] <= 12 v/(g-v)^4 for g > v; for h(t) = q (t-v)^-4, partial summation gives
# sum_{g>T0} h(g) = int_{T0}^inf 4 q (t-v)^-5 (N(t) - K) dt, and N(t) <= t/(2 pi) log(t/(2 pi e)) + E(t), E(t) = 0.1038 log t + 0.2573 log log t + 9.3675
# (Hasanalizade, Shen, Wong 2022, Cor. 1.2, t >= e; as quoted in the book), with log log t <= log t, t <= kappa w, log t <= log w + log kappa
# (w = t - v, kappa = T0/(T0-v)); the integrals of w^-4, w^-4 log w, w^-5, w^-5 log w are elementary and evaluated in Arb.
# Minimum: Q convex with Q'(a) < 0 < Q'(b) on [a, b] = [17.4978, 17.4980]; the tangent lines at a and b meet at a point whose value is a lower bound of
# min_[a,b] Q <= inf X (F >= 0); an upper bound is X(c) <= Q(c) + tail for any c, here c = 17.4979. X' changes sign in [a, b] (X' = Q' + F', F' in [0, tail']).
# Real axis facts used: none (everything is at real v on the critical line under RH).
# Error control: the whole computation is made at two precisions (--prec 160 and 2x) and the maximum difference of the enclosure midpoints is logged.
# Usage: python3 cx_slope_X.py [--quick] [--prec P] [--out FILE]; --quick: K = 100, T0 = 237 (gamma_100 = 236.52, gamma_101 = 237.77), ten seconds on Hopper.
import argparse, json, time
from flint import arb, acb, ctx


def tail(q, v, T0, K):
    """upper bound of sum_{gamma>T0} q (gamma - v)^-4 given N(T0) = K"""
    w0 = T0 - v
    kappa = T0 / w0
    lk = kappa.log()
    A = 1 / (2 * arb.pi())
    c1 = lk - 1 - (2 * arb.pi()).log()
    ev = arb(1038) / 10000 + arb(2573) / 10000
    e0 = arb(93675) / 10000 + ev * lk - K
    lw = w0.log()
    I4 = w0 ** -3 / 3
    I4l = w0 ** -3 * (lw / 3 + arb(1) / 9)
    I5 = w0 ** -4 / 4
    I5l = w0 ** -4 * (lw / 4 + arb(1) / 16)
    return 4 * q * (A * kappa * (I4l + c1 * I4) + e0 * I5 + ev * I5l)


def run(prec, K, T0f, log):
    ctx.prec = prec
    T0 = arb(T0f[0]) / T0f[1]
    t0 = time.time()
    zs = acb.zeta_zeros(1, K + 1)
    gam = [z.imag for z in zs]
    assert gam[K - 1] < T0 < gam[K], "T0 not between gamma_K and gamma_(K+1)"
    assert all(abs(z.real - arb(1) / 2).abs_upper() < 1e-20 for z in zs)
    log("  zeros at prec %d: radius of gamma_K %.2e, time %.1f s" % (prec, float(gam[K - 1].rad()), time.time() - t0))

    def Q(v):
        return sum((1 / (v - g) ** 2 + 1 / (v + g) ** 2 - 2 / g ** 2 for g in gam), arb(0))

    def Qp(v):
        return -2 * sum((1 / (v - g) ** 3 + 1 / (v + g) ** 3 for g in gam), arb(0))

    a = arb(174978) / 10000
    b = arb(17498) / 1000
    c = arb(174979) / 10000
    out = {}
    Fp = {}
    for nm, v in (("a", a), ("b", b)):
        qp = Qp(v)
        tp = tail(12 * v, v, T0, K)
        lo = qp
        hi = qp + tp
        out["Xp_" + nm] = (lo, hi, tp)
        log("  X'(%s) in [Q', Q'+tail'] with Q' = %s, tail' <= %.3e" % (nm, qp.str(12, radius=True), float(tp.upper())))
    Qa, Qb, Qc = Q(a), Q(b), Q(c)
    pa, pb = out["Xp_a"][0], out["Xp_b"][0]
    assert pa < 0 and pb > 0 and (out["Xp_a"][1] < 0) and (out["Xp_b"][0] - 0 > 0)
    w = (Qb - Qa + a * pa - b * pb) / (pa - pb)
    LB = Qa + pa * (w - a)
    tX = tail(6 * c * c, c, T0, K)
    UB = Qc + tX
    log("  Q(a)=%s Q(c)=%s Q(b)=%s" % (Qa.str(12, radius=True), Qc.str(12, radius=True), Qb.str(12, radius=True)))
    log("  tangent intersection w = %s ; lower bound of inf X (Q only) = %s" % (w.str(12, radius=True), LB.str(12, radius=True)))
    log("  tail of F at c: %.3e ; X(c) in [Q(c), Q(c)+tail] ; upper bound of inf X = %s" % (float(tX.upper()), UB.str(12, radius=True)))
    return {"Xp_a": (out["Xp_a"][0], out["Xp_a"][1]), "Xp_b": (out["Xp_b"][0], out["Xp_b"][1]), "Qc": Qc, "LB": LB, "UB": UB, "tailXp_a": out["Xp_a"][2], "tailXp_b": out["Xp_b"][2], "tailX": tX, "w": w}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--prec", type=int, default=160)
    ap.add_argument("--out", default=None)
    A = ap.parse_args()
    out = A.out or ("/tmp/cx_slope_X_quick.json" if A.quick else "cx_slope_X.json")
    logf = open(out.replace(".json", ".log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        logf.write(s + "\n")
        logf.flush()

    K, T0f = (100, (237, 1)) if A.quick else (1000, (14199, 10))
    log("cx_slope_X: K=%d T0=%s" % (K, T0f))
    r1 = run(A.prec, K, T0f, log)
    r2 = run(2 * A.prec, K, T0f, log)
    ctx.prec = 2 * A.prec
    diff = max(abs(r1[k][i].mid() - r2[k][i].mid()) if isinstance(r1[k], tuple) else abs(r1[k].mid() - r2[k].mid())
               for k in ("Xp_a", "Xp_b", "Qc", "LB", "UB", "tailX", "w") for i in ((0, 1) if isinstance(r1[k], tuple) else (0,)))
    log("max difference of midpoints between precision %d and %d: %.3e" % (A.prec, 2 * A.prec, float(diff)))
    ctx.prec = 2 * A.prec
    res = {"K": K, "T0": str(T0f), "prec": [A.prec, 2 * A.prec], "maxdiff_mid": float(diff)}
    for nm in ("a", "b"):
        lo, hi = r2["Xp_" + nm]
        enc = lo.union(hi)
        res["Xp_" + nm] = enc.str(8, radius=True)
        res["Xp_%s_lo_hi" % nm] = [float(lo.mid()), float(hi.mid())]
        log("X'(%s) enclosure (union of [Q', Q'+tail']):" % nm, enc.str(8, radius=True), "mid-range [%.4e, %.4e]" % (float(lo.mid()), float(hi.mid())))
    Xc = r2["Qc"].union(r2["Qc"] + r2["tailX"])
    res["X_c"] = Xc.str(10, radius=True)
    res["inf_X_lower"] = r2["LB"].str(10, radius=True)
    res["inf_X_upper"] = r2["UB"].str(10, radius=True)
    log("X(17.4979) in", Xc.str(10, radius=True))
    log("inf X over [17.4978,17.4980]: lower", r2["LB"].str(10, radius=True), " upper", r2["UB"].str(10, radius=True))
    log("contains floating 0.18389874:", bool(Xc.contains(arb(18389874) / 10 ** 8)))
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
