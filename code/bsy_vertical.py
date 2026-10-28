# bsy_vertical.py: data of Figure bsy:fig:vertical of Chapter 35 (ch/bsy2.tex, Section bsy:sec:mc, "The vertical representation"), group BSY-C.
# Claims supported: at theta = 1/2 the layer cake curves P(f > lambda) and P(f < -lambda) of f(tau) = log|F(1/2 + i tau)| under the Cauchy law P_{1/2}(t0 - .) and the running balance
# A+(L) - A-(L) - log|h(1 + i t0)| = E[clip(f, -L, L)] - log|h(1 + i t0)| (A+(L) = int_0^L P(f > l) dl, A-(L) = int_0^L P(f < -l) dl), which settles at 0 for zeta (t0 = 0 and t0 = 85.7) once L passes
# about 4 and at about 1.43 (the Green potential log(0.8085/0.1915) = 1.44 of the zero off the line, plus about 2e-4 of the others) for the Davenport and Heilbronn function (DH) at t0 = 85.7.
# Conventions. P_a(x) = a/(pi (a^2 + x^2)), a = 1 - theta = 1/2. F = zeta with h(s) = zeta(s)(s - 1)/(s - 2 theta + 1) (h(1) = 1/(2 - 2 theta) = 1 at t0 = 0; log|h(1 + i t0)| = log|zeta(1 + i t0)| +
# log|t0/(1 + i t0)| for t0 != 0), or F = DH with h = F. The window is |tau - t0| <= W (W = 200); the integral is by cells of width hh (midpoint value, weight the exact Cauchy mass of the cell,
# 2/pi arctan differences); the neglected Cauchy mass 1 - (2/pi) arctan(2 W) = 0.0016 outside the window is logged and not added back (it is the same truncation as the quadrature
# columns of Table bsy:tab:defects). The curves are on lambda = 0, 0.05, ..., 12.
# Error control. The cell width is hh = 5e-4 and again 1e-3 (the balance at L = 4 and L = 8 and the tails at lambda = 2 and 6 are compared and the differences logged); the evaluator
# (code/bsy_zdh.py, Euler and Maclaurin) runs at (nf, M) = (0.3, 28) and at (0.45, 36) on 4000 random nodes (maximum difference of f logged).
# Output: data/bsy_vertical.json, the log data/bsy_vertical.log. Usage: python3 bsy_vertical.py [--quick] [--out FILE] [--procs P]; --quick: W = 20, hh = 2e-3 (seconds), to
# misc/reviews/nick-b4/quick/BSY-C/; full run: bash misc/tools/hopper_run.sh -c 10 -t 60 -g bsy_vertical.json -g bsy_vertical.log code/bsy_vertical.py code/bsy_zdh.py -- --out bsy_vertical.json
import argparse, json, math, os, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_zdh as Z

A = 0.5
CASES = [("zeta", 0.0), ("zeta", 85.7), ("dh", 85.7)]


def fval(args):
    kind, t, nf, M = args
    s = 0.5 + 1j * t
    F = Z.zeta(s, nf=nf, M=M) if kind == "zeta" else Z.dh(s, nf=nf, M=M)
    return np.log(np.abs(F))


def logh(kind, t0):
    if kind == "dh":
        return float(np.log(abs(Z.dh(np.array([1 + 1j * t0]))[0])))
    if t0 == 0:
        return 0.0
    s = np.array([1 + 1j * t0])
    return float(np.log(abs(Z.zeta(s)[0])) + math.log(abs(t0) / abs(1 + 1j * t0)))


def curves(f, w, t0, lam, lh):
    order = np.argsort(f); fs = f[order]; cw = np.cumsum(w[order]); tot = cw[-1]
    up = np.array([tot - (cw[np.searchsorted(fs, l, side="right") - 1] if np.searchsorted(fs, l, side="right") > 0 else 0.0) for l in lam])
    lo = np.array([(cw[np.searchsorted(fs, -l, side="left") - 1] if np.searchsorted(fs, -l, side="left") > 0 else 0.0) for l in lam])
    bal = np.array([float((w * np.clip(f, -l, l)).sum()) - lh for l in lam])
    return up, lo, bal, float(tot)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--procs", type=int, default=0)
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    out = Path(args.out) if args.out else (Path("/tmp/bsy_vertical_quick.json") if args.quick else here.parent / "data" / "bsy_vertical.json")
    log = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
    W = 20.0 if args.quick else 200.0
    hhs = [2e-3, 4e-3] if args.quick else [5e-4, 1e-3]
    lam = np.arange(0, 12.0001, 0.05)
    procs = args.procs or (1 if args.quick else (os.cpu_count() or 1))
    pool = None
    if procs > 1:
        import multiprocessing as mpc
        pool = mpc.Pool(procs)
    res = {}
    for kind, t0 in CASES:
        lh = logh(kind, t0)
        res[f"{kind}|{t0}"] = {"log_h": lh}
        for hh in hhs:
            t_edges = t0 - W + hh * np.arange(int(round(2 * W / hh)) + 1)
            mid = 0.5 * (t_edges[1:] + t_edges[:-1])
            w = (2 / math.pi) * 0 + (np.arctan((t_edges[1:] - t0) / A) - np.arctan((t_edges[:-1] - t0) / A)) / math.pi
            # the Cauchy mass of a cell is (atan((b - t0)/a) - atan((a_ - t0)/a))/pi
            chunks = np.array_split(mid, 200 if pool else 1)
            t1 = time.time()
            if pool:
                f = np.concatenate(pool.map(fval, [(kind, c, 0.3, 28) for c in chunks]))
            else:
                f = np.concatenate([fval((kind, c, 0.3, 28)) for c in chunks])
            up, lo, bal, tot = curves(f, w, t0, lam, lh)
            res[f"{kind}|{t0}"][f"hh={hh}"] = {"P_up": up.tolist(), "P_lo": lo.tolist(), "balance": bal.tolist(), "mass": tot, "n_cells": len(mid), "seconds": time.time() - t1, "min_f": float(f.min()), "max_f": float(f.max())}
            P(f"{kind} t0={t0} hh={hh}: cells {len(mid)}, window mass {tot:.6f} (neglected {1 - tot:.6f}), log|h| = {lh:.10f}, min f {f.min():.2f}, max f {f.max():.2f}, balance(4) {bal[80]:.5f}, balance(8) {bal[160]:.5f}, balance(12) {bal[240]:.5f}, {time.time() - t1:.0f} s")
        a_, b_ = res[f"{kind}|{t0}"][f"hh={hhs[0]}"], res[f"{kind}|{t0}"][f"hh={hhs[1]}"]
        P(f"  difference hh {hhs[0]} vs {hhs[1]}: balance(4) {a_['balance'][80] - b_['balance'][80]:.2e}, balance(8) {a_['balance'][160] - b_['balance'][160]:.2e}, P_lo(2) {a_['P_lo'][40] - b_['P_lo'][40]:.2e}, P_lo(6) {a_['P_lo'][120] - b_['P_lo'][120]:.2e}, P_up(2) {a_['P_up'][40] - b_['P_up'][40]:.2e}")
    rg = np.random.default_rng(3)
    for kind in ("zeta", "dh"):
        ts = 85.7 + rg.uniform(-200, 200, 4000)
        P(kind, "evaluator (0.3, 28) vs (0.45, 36), max |difference of f|:", float(np.abs(fval((kind, ts, 0.3, 28)) - fval((kind, ts, 0.45, 36))).max()))
    out.write_text(json.dumps({"W": W, "hh": hhs, "lambda": lam.tolist(), "cases": res}))
    Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
