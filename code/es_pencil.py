# es_pencil.py: Hua and Sarkar's matrix pencil applied to the tapered Fourier transform of the level theta data u_theta, the computation of the
# subsection "Above the verified range" of the explicit sparse chapter (zeros of zeta recovered, ordinates and real parts, from the data u_theta alone),
# run at a height above the verified range of Platt and Trudgian (2021, RH verified for 0 < gamma <= 3 000 175 332 800): the window centred at
# c = 10^13 + 100. Added for the book in October 2026, written from the mathematics (the author's scripts were not read, run or copied).
# Conventions. u_theta(v) = pi^(-1) Re (xi'/xi)(theta + i v) = sum_rho a_rho/(pi (a_rho^2 + (v - gamma_rho)^2)) with a_rho = theta - sigma_rho, the sum over all
# zeros rho = sigma_rho + i gamma_rho (xi'/xi = sum 1/(s - rho)); theta = 0.6 here. xi'/xi(s) = 1/s + 1/(s - 1) - (1/2) log pi + (1/2) psi(s/2) + zeta'/zeta(s).
# Data. u_theta(v_k), v_k = c - 10 + 0.02 k, k = 0..1000 (1001 points, 5 ell on both sides of c for the taper ell = 2): zeta(s) and zeta'(s) from the
# first two Taylor coefficients of zeta at s = 0.6 + i v_k (python-flint acb_series, zeta, 160 bits), psi(s/2) from acb.digamma; every value is a ball
# and its radius is recorded. Arrays of tasks (--task I --ntasks N) each compute a slice of k and write a JSON file of the midpoints (40 digits) and radii.
# Method. F(omega) = int u_theta(v) w(v) exp(i omega (v - c)) dv with w(v) = exp(-(v - c)^2/(2 ell^2)); the trapezoid sum with step h = 0.02 equals F up
# to aliasing of order exp(-2 pi a/h) (a >= 0.1: below 1e-26) and to the truncation of the taper at +-5 ell (relative weight e^(-12.5) = 3.7e-6, the limit of
# the accuracy for the zeros at the edge). Theorem (tapered transform): for omega_j = omega_0 + j Delta, F(omega_j) = sum_rho A_rho z_rho^j + O(exp(-ell^2 omega_0^2/2)),
# z_rho = exp(Delta(-a_rho + i (gamma_rho - c))), A_rho = w(gamma_rho + i a_rho) exp(-a_rho omega_0) exp(i omega_0 (gamma_rho - c)), |A_rho| = w(gamma_rho)
# exp(a_rho^2/(2 ell^2) - a_rho omega_0). The pencil: omega in [5, 85] step Delta = 0.2 (401 samples), the Hankel matrix H[i, k] = F(omega_{i+k}), 0 <= i <= 200,
# 0 <= k < 200, H_0 = H without its last row, H_1 = H without its first row, the SVD H_0 = U S V*, model order M = the number of singular values above 10 times the plateau
# of the spectrum (the median of its second half; --order overrides; the validation sweeps M with --sweep), poles = eigenvalues of U_M* H_1 V_M S_M^(-1). gamma = c + arg(z)/Delta, a = -log|z|/Delta, sigma = theta - a (the real part is
# read from the modulus of the pole, nothing is assumed about it). Amplitudes by least squares of F against the Vandermonde matrix of the poles; a pole is a zero
# when |A| lies within a factor 3 of the profile |A_rho| above computed with the pole's OWN a (not theta - 1/2); poles with |gamma - c| > 7.5 are dropped.
# Error control. The data of every task are re-evaluated at 1.25 times the precision at the first point of the slice (--check) and the difference is logged;
# the pencil is rerun at the model orders M - 5, M and M + 5 and at two precisions of the transform (mpmath-free: double and the same double sums in
# longdouble) and the maximum change of the ordinates and of sigma - 1/2 of the zeros present in all runs is logged.
# Modes. gen (a slice of the data), pencil (aggregate data files, run the pencil, optionally compare with a reference), and --quick (all in one process:
# c = 200, 96 bits, validation against data/zeros_1700.txt, seconds).
# Reference for the comparison. --ref data/zeros_1700.txt (68 digits, one ordinate a line) or --ref a JSON of the sign scan (data/es_zeros1e13.json,
# zeros with "ordinate", "lo", "hi"); nearest neighbour matching inside the window |gamma - c| <= 6.
# Output: the JSON of the poles (ordinate, real part, amplitude) given by --out, and the comparison printed and stored in it.
import sys, json, argparse, time, math
import numpy as np
from decimal import Decimal

THETA = 0.6; ELL = 2.0; H = 0.02; NPTS = 1001; HALF = 10.0
SWEEP = False
OMEGA0 = 5.0; DOMEGA = 0.2; NOMEGA = 401; ORDER_H = 200; WIN = 6.0; WIN_KEEP = 7.5


def gen(c, k0, k1, prec, check=False):
    from flint import acb, arb, ctx, acb_series, fmpq
    out = []
    t0 = time.time()
    maxdiff = 0.0

    def u_at(k, p):
        ctx.prec = p
        v = arb(fmpq(c * 50 - 500 + k, 50))  # v = c - 10 + 0.02 k exactly
        s = acb(arb(fmpq(3, 5)), v)
        ser = acb_series([s, 1], prec=2).zeta()
        zt, dz = ser[0], ser[1]
        xi = 1 / s + 1 / (s - 1) - arb.pi().log() / 2 + acb.digamma(s / 2) / 2 + dz / zt
        return xi.real / arb.pi()

    for k in range(k0, k1):
        u = u_at(k, prec)
        rec = dict(k=k, mid=u.mid().str(40, radius=False), rad=float(u.rad()))
        if check and k == k0:
            u2 = u_at(k, int(prec * 1.25))
            maxdiff = abs(float(u.mid() - u2.mid()))
            rec["check_diff"] = maxdiff; rec["check_prec"] = int(prec * 1.25)
        out.append(rec)
        if (k - k0) % 10 == 0:
            print("k", k, "u", rec["mid"][:20], "rad", rec["rad"], "t=%.0fs" % (time.time() - t0), flush=True)
    return dict(c=c, theta=THETA, prec=prec, k0=k0, k1=k1, seconds=time.time() - t0, check_diff=maxdiff, data=out)


def load_data(files):
    d = {}; meta = None; mr = 0.0; cd = 0.0
    for fn in files:
        j = json.load(open(fn))
        for r in j["data"]:
            d[r["k"]] = float(r["mid"]); mr = max(mr, r["rad"])
        cd = max(cd, j.get("check_diff", 0.0)); meta = j
    assert sorted(d) == list(range(NPTS)), "data incomplete: %d points" % len(d)
    return np.array([d[k] for k in range(NPTS)]), meta, mr, cd


def transform(u, dtype=np.float64):
    x = (-HALF + H * np.arange(NPTS)).astype(dtype)  # v - c, exact to rounding
    w = np.exp(-x ** 2 / (2 * ELL ** 2))
    om = OMEGA0 + DOMEGA * np.arange(NOMEGA)
    cdt = np.complex128 if dtype == np.float64 else np.clongdouble
    E = np.exp(1j * np.outer(om.astype(dtype), x))
    return (H * (E @ (u.astype(dtype) * w))).astype(np.complex128), om


def pencil(F, c, order=None, kmin=10):
    J = len(F); L = ORDER_H
    Hk = np.array([[F[i + k] for k in range(L)] for i in range(J - L)])  # (J-L) x L = 201 x 200
    H0 = Hk[:-1]; H1 = Hk[1:]
    U, S, Vh = np.linalg.svd(H0)
    # model order: the singular values fall to a plateau (the noise floor, the median of the second half of the spectrum); M is the number of
    # singular values above 10 times that floor
    floor = float(np.median(S[L // 2:]))
    M = int(np.sum(S > 10 * floor)) if order is None else int(order)
    Um = U[:, :M]; Vm = Vh[:M].conj().T
    A = Um.conj().T @ H1 @ Vm / S[:M]
    z = np.linalg.eigvals(A)
    gam = c_float_offset(z) ; a = -np.log(np.abs(z)) / DOMEGA
    V = np.exp(np.outer(np.arange(J), np.log(z.astype(complex))))
    nrm = np.linalg.norm(V, axis=0)  # column scaling: poles off the unit circle make the raw Vandermonde matrix hopeless
    cf, *_ = np.linalg.lstsq(V / nrm, F, rcond=1e-13)
    amp = cf / nrm
    return dict(M=M, S=S, z=z, dgam=gam, a=a, amp=amp)


def c_float_offset(z):
    return np.angle(z) / DOMEGA  # gamma - c


def profile(dg, a):
    return np.exp(-dg ** 2 / (2 * ELL ** 2) + a ** 2 / (2 * ELL ** 2) - a * OMEGA0)


def select(res):
    dg, a, amp = res["dgam"], res["a"], res["amp"]
    ok = (np.abs(dg) <= WIN_KEEP) & (np.abs(a) < 2)
    pr = profile(dg, a)
    ratio = np.abs(amp) / pr
    ok &= (ratio > 1 / 3) & (ratio < 3)
    idx = np.where(ok)[0]
    idx = idx[np.argsort(dg[idx])]
    return idx, ratio


def read_ref(fn, c):
    from fractions import Fraction
    if fn.endswith(".json"):
        j = json.load(open(fn))
        return [(Fraction(z["ordinate"]) - c, Fraction(z["lo"]) - c, Fraction(z["hi"]) - c) for z in j["zeros"]]
    zs = []
    for ln in open(fn):
        ln = ln.split()
        if ln:
            g = Fraction(ln[0])
            if abs(g - c) <= WIN + 1: zs.append((g - c, g - c, g - c))
    return zs


def compare(dgs, ref):
    # dgs: list of ordinates (offsets from c) of the pencil's zeros within the window; ref: list of (mid, lo, hi) offsets
    rd = [float(r[0]) for r in ref if abs(r[0]) <= WIN]
    pd = [d for d in dgs if abs(d) <= WIN]
    rows = []; used = set()
    for d in pd:
        j = min(range(len(rd)), key=lambda i: abs(rd[i] - d)) if rd else None
        if j is None: continue
        rows.append((d, rd[j], d - rd[j])); used.add(j)
    unmatched_ref = [rd[i] for i in range(len(rd)) if i not in used]
    # the exact bracket test
    inbr = sum(1 for d in pd if any(float(r[1]) - 1e-7 <= d <= float(r[2]) + 1e-7 for r in ref if abs(r[0]) <= WIN))
    return dict(n_pencil=len(pd), n_ref=len(rd), max_abs_diff=max((abs(r[2]) for r in rows), default=None),
                median_abs_diff=float(np.median([abs(r[2]) for r in rows])) if rows else None,
                n_inside_bracket_pm1e7=inbr, unmatched_ref=unmatched_ref, rows=rows)


def run_pencil(u, c, order, ref, out, meta=None, extra=None):
    t0 = time.time()
    F, om = transform(u)
    res = pencil(F, c, order)
    idx, ratio = select(res)
    dg = res["dgam"]; a = res["a"]; sig = THETA - a
    zeros = [dict(ordinate="%.9f" % (Decimal(c) + Decimal(float(dg[i]))), offset=float(dg[i]), sigma=float(sig[i]),
                  sigma_minus_half=float(sig[i] - 0.5), amplitude=float(abs(res["amp"][i])), amp_over_profile=float(ratio[i])) for i in idx]
    inwin = [z for z in zeros if abs(z["offset"]) <= WIN]
    print("model order", res["M"], "poles kept", len(zeros), "in window", len(inwin),
          "largest |sigma-1/2| in window %.3g" % max((abs(z["sigma_minus_half"]) for z in inwin), default=float("nan")))
    # stability: other orders and long double transform
    stab = {}
    base = {round(z["offset"], 5) for z in inwin}
    runs = {}
    for name, kw in (("M-5", dict(order=res["M"] - 5)), ("M+5", dict(order=res["M"] + 5)), ("longdouble", dict(order=res["M"]))):
        FF = transform(u, np.longdouble)[0] if name == "longdouble" else F
        r2 = pencil(FF, c, kw["order"]); i2, rt2 = select(r2)
        runs[name] = {float(r2["dgam"][i]): float(THETA - r2["a"][i]) for i in i2 if abs(r2["dgam"][i]) <= WIN}
    for name, d in runs.items():
        pairs = []
        for z in inwin:
            if d:
                g2 = min(d, key=lambda g: abs(g - z["offset"]))
                if abs(g2 - z["offset"]) < 1e-3: pairs.append((abs(g2 - z["offset"]), abs((d[g2] - 0.5) - z["sigma_minus_half"])))
        stab[name] = dict(count=len(d), common=len(pairs), max_ordinate_change=max((p[0] for p in pairs), default=None),
                          max_sigma_change=max((p[1] for p in pairs), default=None))
    print("stability", stab)
    cmpd = None
    if ref:
        cmpd = compare([z["offset"] for z in zeros], read_ref(ref, c))
        print("compare: pencil in window", cmpd["n_pencil"], "reference", cmpd["n_ref"], "max |diff| %s" % cmpd["max_abs_diff"],
              "median %s" % cmpd["median_abs_diff"], "inside bracket +-1e-7:", cmpd["n_inside_bracket_pm1e7"], "unmatched ref", cmpd["unmatched_ref"])
    sweep = {}
    if ref and SWEEP:
        for M2 in range(10, 101, 5):
            r3 = pencil(F, c, M2); i3, _ = select(r3)
            cc = compare([float(r3["dgam"][i]) for i in i3], read_ref(ref, c))
            sweep[M2] = (cc["n_pencil"], cc["max_abs_diff"], max((abs(THETA - r3["a"][i] - 0.5) for i in i3 if abs(r3["dgam"][i]) <= WIN), default=None))
            print("sweep M", M2, "n", cc["n_pencil"], "ref", cc["n_ref"], "maxdiff", cc["max_abs_diff"], "max|sigma-1/2|", sweep[M2][2])
    sv = [float(x) for x in res["S"][:140]]
    d = dict(c=c, theta=THETA, ell=ELL, h=H, omega0=OMEGA0, domega=DOMEGA, nomega=NOMEGA, hankel=ORDER_H, model_order=res["M"],
             singular_values_first140=sv, kept=len(zeros), in_window=len(inwin), window=WIN,
             max_abs_sigma_minus_half_in_window=max((abs(z["sigma_minus_half"]) for z in inwin), default=None),
             data_radius_max=meta, stability=stab, order_sweep=sweep, comparison=cmpd, seconds=time.time() - t0, zeros=zeros)
    if extra: d.update(extra)
    if out: json.dump(d, open(out, "w"), indent=1)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="pencil", choices=["gen", "pencil"])
    ap.add_argument("--quick", action="store_true"); ap.add_argument("--c", type=int, default=None)
    ap.add_argument("--task", type=int, default=0); ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--prec", type=int, default=160); ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", default=None); ap.add_argument("--data", nargs="*"); ap.add_argument("--ref", default=None)
    ap.add_argument("--order", type=int, default=None); ap.add_argument("--sweep", action="store_true")
    a = ap.parse_args()
    global SWEEP
    SWEEP = a.sweep
    if a.quick:
        c = a.c or 200
        j = gen(c, 0, NPTS, 96, check=True)
        u = np.array([float(r["mid"]) for r in j["data"]])
        mr = max(r["rad"] for r in j["data"])
        print("data: %d points, max radius %.2g, precision check difference %.2g, %.1fs" % (NPTS, mr, j["check_diff"], j["seconds"]))
        run_pencil(u, c, a.order, a.ref or "data/zeros_1700.txt", a.out, meta=mr, extra=dict(prec=96, check_diff=j["check_diff"]))
        return
    if a.mode == "gen":
        c = a.c
        k0 = a.task * NPTS // a.ntasks; k1 = (a.task + 1) * NPTS // a.ntasks
        j = gen(c, k0, k1, a.prec, check=a.check)
        json.dump(j, open(a.out, "w"))
        print("task", a.task, "points", k1 - k0, "seconds", j["seconds"], "check_diff", j["check_diff"])
        return
    u, meta, mr, cd = load_data(a.data)
    run_pencil(u, a.c, a.order, a.ref, a.out, meta=mr, extra=dict(prec=meta["prec"], check_diff=cd))


main()
