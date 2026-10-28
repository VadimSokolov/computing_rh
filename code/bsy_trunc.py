# bsy_trunc.py: Table bsy:tab:trunc of Chapter 35, Section bsy:sec:trunc (truncation of the Cauchy average on the whole line), task N2 of reconcile-M3.md. Added for the book in
# October 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Claims supported: the truncated Cauchy average of log|zeta(theta + i tau)| over |tau - t_0| <= T on a randomly shifted grid of step h (N = 2T/h nodes, tau_j = t_0 - T + h (j + U))
# has, minus log|h_theta(1 + i t_0)|, a mean error equal to minus the Cauchy tail beyond T (exactly, up to the spread over shifts, because a grid of window length N h with one uniform
# shift is unbiased for the window integral); theta = 3/4 with nothing deflated, theta = 1/2 with every zero rho_k = 1/2 + i gamma_k, |gamma_k - t_0| <= T + 1, deflated by
# u_k(s) = log|(s - rho_k)(s - rho_k + 2)/(s - rho_k + 1)^2| (on the line (1/2) log(x^2 (4 + x^2)/(1 + x^2)^2), x = tau - gamma_k, harmonic on Re s > 1/2, O(|s|^-2), integral over R equal to 0).
# Conventions. t_0 = 85.7, Poisson weight P(d) = b/(pi (b^2 + d^2)) with b = 1 - theta (the point 1 + i t_0 is at distance 1 - theta from the line Re s = theta), d = t_0 - tau,
# h_theta(s) = zeta(s) (s - 1)/(s - 2 theta + 1). The shifts U_i = numpy default_rng(20261025 + i).random(), i = 0..n-1 (n = 10), the same U for every T and theta, so that the grids of
# the rows are nested (the grid of T is a sub-grid of the grid of the largest T of the same theta and the same U; each row is still a random shifted grid of window length exactly N h).
# Variants at theta = 1/2. Window: estimate = h sum_j [(log|zeta| - sum_k u_k) P](tau_j) + sum_k W_k, W_k = int_{|tau - t_0| <= T} u_k P d tau by mpmath quadrature with breakpoints
# at t_0, gamma_k and t_0 (resp. gamma_k) +- 2^m. Whole: the same with sum_k u_k(1 + i t_0) in place of sum_k W_k (the whole line mean, exact since u_k is harmonic with u_k = O(|s|^-2)),
# so that Error(whole) = Error(window) + sum_k tail_k, tail_k = int_{|tau - t_0| > T} u_k P d tau (the deflation tails, quadrature to +-infinity, independent of W_k).
# Per row: mean error over the shifts, spread (standard deviation, ddof 1), T^2 |mean error|, the ratio spread / |mean error|, mean / (spread/sqrt n), the zero count, the deflation tails.
# Zero range convention: the zeros with |gamma - t_0| <= T + 1, gamma over all ordinates of the zeros (negative ordinates by conjugation: the zeros are +-g, g in data/zeros_1700.txt, 68 digits);
# the count is checked against mpmath.nzeros at both ends of the range (N(b) - N(a), N(-x) = -N(x)) and the first ordinates against mpmath.zetazero.
# Error control. (1) log|zeta| on the grid: python-flint acb.zeta at 64 bits (the ball radius of log|zeta| is recorded for every node and its maximum is logged), without flint mpmath at
# 18 digits; a second evaluator (flint at 128 bits and mpmath at 30 digits) on a random sample of nodes: max difference logged. (2) Quadratures at 20 digits and again at 30 digits: the difference of
# every sum (sum W_k, sum tails) is logged, and for every zero W_k + tail_k (left and right) is compared with the exact u_k(1 + i t_0); the maximum is logged. (3) The reference
# log|h_theta(1 + i t_0)| at 30 and 50 digits (difference logged). Sums over the grid are in double precision (about 1e-12 for 1.3e5 nodes).
# Output: data/bsy_trunc.json and the log bsy_trunc.log (via --merge); a run covers one or two thetas. Usage: python3 bsy_trunc.py [--quick] [--out FILE] [--theta 0.75|0.5|both] [--Ts a,b,c] [--nshift n]
# [--zeros FILE] [--procs p]; python3 bsy_trunc.py --merge A.json B.json --out FILE writes the merged json and its log. --quick: T = 10, 20, 3 shifts, a few seconds, writes to /tmp unless --out.
# Cluster: bash misc/tools/hopper_run.sh -c 8 -t 120 -g bsy_trunc_t75.json -g bsy_trunc_t75.log code/bsy_trunc.py data/zeros_1700.txt -- --theta 0.75 --Ts 10,40,160,640,1280 --out bsy_trunc_t75.json
# and the same with --theta 0.5 --Ts 10,40,160,640 --out bsy_trunc_t50.json; then merge locally.
import argparse, json, math, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import mpmath as mp
try:
    from flint import acb, arb, ctx
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False

T0S = "85.7"
T0 = float(T0S)
here = Path(__file__).resolve().parent
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)

# ---------------------------------------------------------------- log|zeta| on the grid
def eval_chunk(args):
    sigma, taus, prec = args
    out = np.empty(len(taus)); mr = 0.0
    if HAVE_FLINT:
        ctx.prec = prec
        for i, t in enumerate(taus):
            z = acb.zeta(acb(sigma, abs(float(t))))
            l = (z.real * z.real + z.imag * z.imag).log() / 2
            out[i] = float(l.mid()); r = float(l.rad())
            if not (r == r): r = float("inf")
            mr = max(mr, r)
    else:
        mp.mp.dps = prec
        for i, t in enumerate(taus):
            z = mp.zeta(mp.mpc(sigma, abs(float(t))))
            out[i] = float(mp.log(abs(z)))
    return out, mr

def eval_many(pool, sigma, taus, prec, csize=1500):
    chunks = [(sigma, taus[i:i + csize], prec) for i in range(0, len(taus), csize)]
    res = pool.map(eval_chunk, chunks)
    return np.concatenate([r[0] for r in res]), max(r[1] for r in res)

# ---------------------------------------------------------------- the deflation functions and their quadratures
def u_line(x):
    return mp.log(abs(x)) + mp.log(4 + x * x) / 2 - mp.log(1 + x * x)

def u_exact(gam, t0):          # u_k(1 + i t_0), s - rho = 1/2 + i (t_0 - gamma)
    z = mp.mpc(mp.mpf(1) / 2, t0 - gam)
    return mp.log(abs(z * (z + 2) / (z + 1) ** 2))

def breaks(a, b, centers):
    s = {a, b}
    for c in centers:
        if a < c < b: s.add(c)
        k = 0
        while 2 ** k < (b - a):
            for p in (c + 2 ** k, c - 2 ** k):
                if a < p < b: s.add(p)
            k += 1
    return sorted(s)

def zero_task(args):
    gstr, T, dps, bstr = args
    mp.mp.dps = dps
    g = mp.mpf(gstr); t0 = mp.mpf(T0S); b = mp.mpf(bstr); T = mp.mpf(T)
    f = lambda tau: u_line(tau - g) * b / (mp.pi * (b * b + (t0 - tau) ** 2))
    lo, hi = t0 - T, t0 + T
    W = mp.quad(f, breaks(lo, hi, [t0, g]))
    big = mp.mpf(2) ** 40
    pr = breaks(hi, hi + big, [g]); pr.append(mp.inf)
    pl = breaks(lo - big, lo, [g]); pl.insert(0, -mp.inf)
    Tr = mp.quad(f, pr); Tl = mp.quad(f, pl)
    ex = u_exact(g, t0)
    return float(W), float(Tl + Tr), float(W + Tl + Tr - ex)

def href(theta, dps):
    mp.mp.dps = dps
    s = mp.mpc(1, mp.mpf(T0S))
    return mp.log(abs(mp.zeta(s) * (s - 1) / (s - 2 * mp.mpf(theta) + 1)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
    ap.add_argument("--theta", default="both"); ap.add_argument("--Ts", default=None)
    ap.add_argument("--nshift", type=int, default=None); ap.add_argument("--zeros", default=None)
    ap.add_argument("--procs", type=int, default=None); ap.add_argument("--merge", nargs="+", default=None)
    args = ap.parse_args()
    if args.merge:
        return merge(args)
    out = Path(args.out) if args.out else (Path("/tmp/bsy_trunc_quick.json") if args.quick else here.parent / "data" / "bsy_trunc.json")
    nsh = args.nshift or (3 if args.quick else 10)
    thetas = [0.75, 0.5] if args.theta == "both" else [float(args.theta)]
    procs = args.procs or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count())
    Tdef = {0.75: [10, 20] if args.quick else [10, 40, 160, 640, 1280], 0.5: [10, 20] if args.quick else [10, 40, 160, 640]}
    hh = {0.75: 0.02, 0.5: 0.01}
    prec1, prec2 = (64, 128) if HAVE_FLINT else (18, 30)
    zf = Path(args.zeros) if args.zeros else None
    if zf is None:
        for c in (here.parent / "data" / "zeros_1700.txt", here / "zeros_1700.txt", Path("zeros_1700.txt")):
            if c.exists(): zf = c; break
    zstr = [l.strip() for l in open(zf) if l.strip()]
    gpos = np.array([float(s) for s in zstr])
    P(f"bsy_trunc: t_0 = {T0S}, evaluator {'python-flint at %d bits' % prec1 if HAVE_FLINT else 'mpmath at %d digits' % prec1}, second evaluator at {prec2}, shifts {nsh}, processes {procs}, zeros {zf} ({len(zstr)} ordinates)")
    U = [float(np.random.default_rng(20261025 + i).random()) for i in range(nsh)]
    P("shifts U_i (seed 20261025 + i):", " ".join("%.6f" % u for u in U))
    res = {"t0": T0S, "nshift": nsh, "seeds": [20261025 + i for i in range(nsh)], "U": U, "rows": [], "checks": {}}
    chk = res["checks"]
    pool = Pool(procs)
    tstart = time.time()
    for theta in thetas:
        h = hh[theta]; b = 1.0 - theta; Ts = [int(x) for x in args.Ts.split(",")] if args.Ts else Tdef[theta]
        Tmax = max(Ts); Kmax = int(round(Tmax / h))
        r30 = href(theta, 30); r50 = href(theta, 50)
        ref = float(r30)
        P(f"theta {theta}: reference log|h_theta(1 + i t_0)| = {mp.nstr(r30, 17)} (30 digits minus 50 digits: {mp.nstr(r30 - r50, 3)}), h = {h}, b = {b}, T in {Ts}")
        chk[f"ref_{theta}"] = {"value": mp.nstr(r30, 20), "diff30_50": mp.nstr(r30 - r50, 3)}
        # grid values at the largest window, per shift
        ks = np.arange(-Kmax, Kmax)
        F = []; maxrad = 0.0; t1 = time.time()
        for i in range(nsh):
            d = h * (ks + U[i]); taus = T0 + d
            f1, mr = eval_many(pool, theta, taus, prec1)
            F.append((d, f1)); maxrad = max(maxrad, mr)
        P(f"theta {theta}: evaluated log|zeta| at {nsh} x {len(ks)} nodes in {time.time() - t1:.1f} s; max ball radius of log|zeta| {maxrad:.2e}")
        rng = np.random.default_rng(777)
        samp = []
        for i in range(nsh):
            idx = rng.choice(len(ks), size=min(150, len(ks)), replace=False); samp.append((i, idx))
        sd = 0.0
        for i, idx in samp:
            d, f1 = F[i]; taus = T0 + d[idx]
            f2, _ = eval_many(pool, theta, taus, prec2, csize=200)
            sd = max(sd, float(np.max(np.abs(f2 - f1[idx]))))
        if HAVE_FLINT:     # independent evaluator: mpmath at 30 digits on 100 nodes
            i, idx = samp[0]; d, f1 = F[i]
            mp.mp.dps = 30
            sd2 = max(abs(float(mp.log(abs(mp.zeta(mp.mpc(theta, abs(float(T0 + d[j]))))))) - f1[j]) for j in idx[:100])
        else:
            sd2 = sd
        P(f"theta {theta}: second evaluator on {sum(len(ix) for _, ix in samp)} sample nodes: max |difference| {sd:.2e}" + (f"; mpmath 30 digits on 100 nodes: {sd2:.2e}" if HAVE_FLINT else ""))
        chk[f"eval_{theta}"] = {"max_radius": maxrad, "second_eval_maxdiff": sd, "mpmath30_maxdiff": sd2, "nodes": len(ks) * nsh}
        for T in Ts:
            K = int(round(T / h)); sl = slice(Kmax - K, Kmax + K); N = 2 * K
            row = {"theta": theta, "h": h, "T": T, "N": N, "nshift": nsh}
            if theta == 0.75:
                errs = []
                for i in range(nsh):
                    d, f1 = F[i]; dd = d[sl]; w = b / (math.pi * (b * b + dd * dd))
                    errs.append(h * float(np.sum(f1[sl] * w)) - ref)
                row["err"] = errs
            else:
                lo, hi = T0 - T - 1, T0 + T + 1
                gam = np.concatenate([gpos, -gpos]); gam = np.sort(gam[(gam >= lo) & (gam <= hi)])
                sel = [s for s, g in zip(zstr, gpos) if lo <= g <= hi] + ["-" + s for s, g in zip(zstr, gpos) if lo <= -g <= hi]
                nz = len(sel); assert nz == len(gam)
                mp.mp.dps = 20
                nb = mp.nzeros(hi); na = mp.nzeros(lo) if lo > 0 else -mp.nzeros(-lo)
                row["zeros"] = nz; row["zeros_nzeros"] = int(nb - na)
                P(f"theta {theta} T {T}: zeros with |gamma - t_0| <= T + 1 (gamma in [{lo:.1f}, {hi:.1f}]): {nz}; mpmath.nzeros: N({hi:.1f}) - N({lo:.1f}) = {int(nb - na)}")
                assert nz == int(nb - na)
                tq = time.time()
                r20 = pool.map(zero_task, [(s, T, 20, "0.5") for s in sel], chunksize=4)
                r30_ = pool.map(zero_task, [(s, T, 30, "0.5") for s in sel], chunksize=4)
                Wk = np.array([r[0] for r in r20]); Tk = np.array([r[1] for r in r20]); cl = max(abs(r[2]) for r in r20)
                W30 = sum(r[0] for r in r30_); T30 = sum(r[1] for r in r30_); cl30 = max(abs(r[2]) for r in r30_)
                Wsum = float(math.fsum(Wk)); Tsum = float(math.fsum(Tk)); Usum = Wsum + Tsum
                row.update({"sum_W": Wsum, "defl_tail": Tsum, "sum_u_line": Usum, "quad_diff_W": float(W30 - Wsum), "quad_diff_tail": float(T30 - Tsum), "closure_max": max(cl, cl30)})
                P(f"theta {theta} T {T}: sum W_k {Wsum:.12e}, deflation tail sum {Tsum:.6e} (quad 20 vs 30 digits: {W30 - Wsum:.1e}, {T30 - Tsum:.1e}); closure W + tails - u_k(1+it_0) max {max(cl, cl30):.1e}; quads {time.time() - tq:.1f} s")
                errw = []; errl = []
                for i in range(nsh):
                    d, f1 = F[i]; dd = d[sl]; taus = T0 + dd; w = b / (math.pi * (b * b + dd * dd))
                    S = np.zeros(len(dd))
                    for j in range(0, len(dd), 2000):
                        x = taus[j:j + 2000, None] - gam[None, :]; x2 = x * x
                        S[j:j + 2000] = (0.5 * (np.log(x2) + np.log(4 + x2) - 2 * np.log1p(x2))).sum(axis=1)
                    base = h * float(np.sum((f1[sl] - S) * w))
                    errw.append(base + Wsum - ref); errl.append(base + Usum - ref)
                row["err"] = errw; row["err_whole"] = errl
            for key in ("err", "err_whole"):
                if key in row:
                    e = np.array(row[key]); m = float(e.mean()); sdv = float(e.std(ddof=1))
                    row[key + "_stats"] = {"mean": m, "spread": sdv, "T2abs": T * T * abs(m), "spread_over_abs": sdv / abs(m) if m else None, "t": m / (sdv / math.sqrt(nsh)) if sdv else None}
            res["rows"].append(row)
            st = row["err_stats"]
            P(f"theta {theta} T {T} N {N}: mean error {st['mean']:+.4e}, spread {st['spread']:.2e}, T^2|err| {st['T2abs']:.4f}, spread/|err| {st['spread_over_abs']:.3f}" + (f", zeros {row['zeros']}; whole line add-back: mean {row['err_whole_stats']['mean']:+.4e}, spread {row['err_whole_stats']['spread']:.2e}, T^2|err| {row['err_whole_stats']['T2abs']:.4f}" if "err_whole" in row else ""))
    pool.close(); pool.join()
    # first ordinates against mpmath.zetazero (full runs only)
    mp.mp.dps = 25
    nchk = [1, 10, 100] if args.quick else [1, 10, 100, 500, 742]
    dz = max(abs(mp.zetazero(n).imag - mp.mpf(zstr[n - 1])) for n in nchk)
    chk["zetazero_maxdiff"] = mp.nstr(dz, 3); P(f"ordinates n = {nchk} against mpmath.zetazero (25 digits): max |difference| {mp.nstr(dz, 3)}")
    res["seconds"] = time.time() - tstart
    json.dump(res, open(out, "w"), indent=1)
    open(out.with_suffix(".log"), "w").write("\n".join(log) + "\n")
    P(f"wrote {out}, total {time.time() - tstart:.1f} s")

def merge(args):
    out = Path(args.out); R = None
    for f in args.merge:
        d = json.load(open(f))
        if R is None: R = d
        else:
            R["rows"] += d["rows"]; R["checks"].update(d["checks"]); R["seconds"] += d["seconds"]
    R["rows"].sort(key=lambda r: (-r["theta"], r["T"]))
    json.dump(R, open(out, "w"), indent=1)
    P("merged table (t_0 = %s, %d shifts, seeds %s..)" % (R["t0"], R["nshift"], R["seeds"][0]))
    for r in R["rows"]:
        s = r["err_stats"]
        P("theta %.2f h %.2f T %5d N %7d: err %+.3e spread %.2e T2|err| %.4f spread/|err| %.3f t %.2f%s" % (r["theta"], r["h"], r["T"], r["N"], s["mean"], s["spread"], s["T2abs"], s["spread_over_abs"], s["t"] or 0,
          (" zeros %d tail %.3e | whole: err %+.3e spread %.2e T2|err| %.4f" % (r["zeros"], r["defl_tail"], r["err_whole_stats"]["mean"], r["err_whole_stats"]["spread"], r["err_whole_stats"]["T2abs"])) if "zeros" in r else ""))
    for k, v in sorted(R["checks"].items()): P("check", k, v)
    open(out.with_suffix(".log"), "w").write("\n".join(log) + "\n")

if __name__ == "__main__":
    main()
