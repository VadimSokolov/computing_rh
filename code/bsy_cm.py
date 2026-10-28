# bsy_cm.py: Table tab:cm of Chapter 35 (Balazard, Saias and Yor meet Thorin), subsection "The complex logarithm: a conformal martingale", group BSY-A, task C2 of
# misc/reviews/nick-b4/reconcile-M.md. Claims supported: Theorem bsy:thm:cmart (E X = x0, E Y = y0, Var X = Var Y, Cov(X, Y) = 0 under QRH(theta)) and the six rows of the table
# (zeta at theta = 1/2 for t0 = 0, 85.7, 1000; zeta at theta = 3/4, t0 = 85.7; the Davenport and Heilbronn function f at theta = 1/2 and 3/4, t0 = 85.7), M2:355 to M2:372.
# Normalisation: Chapters 30 to 35 use the measure without the factor 1/2 (misc/reviews/nick-b4/decisions.md); not used here, the table needs only the harmonic measure.
# Conventions. h_theta(s) = zeta(s) (s-1)/(s-2 theta+1) for zeta (h = f for the entire function f of Davenport and Heilbronn), a = 1 - theta, tau has the law
# P(tau) = (1/pi) a/(a^2 + (tau - t0)^2) (harmonic measure of {Re s > theta} seen from 1 + i t0), Z = X + iY = L(theta + i tau) with L the logarithm of h along the horizontal ray
# from Re s = +infinity (L -> 0 there): X = log|h|, Y = arg h continued along the ray (NOT the raw principal argument of zeta: that fails the identities, reconcile-M NEW-cmart),
# z0 = x0 + i y0 = L(1 + i t0) (for t0 = 0, theta = 1/2: h(1) = 1, so z0 = 0). Rows are lines of the table: E X - x0, E Y - y0, Var X, Var Y, Cov(X, Y), with the full
# variance corrected for the tail, Var X = I_2 - I_1^2 + r_2 - 2 I_1 r_1 - r_1^2 (I_j the integral over the window, r_j the contribution of the part beyond it).
# Method. Window |tau - t0| <= 13000 (quick: 400). On |tau - t0| <= 1000 (quick: 100) panels between breakpoints with the tanh-sinh rule (step 0.125 and 0.0625 for the two meshes;
# weights = Poisson density times rule weight): the breakpoints are the zeros of zeta (ordinates of data/zeros_1700.txt, both signs, theta = 1/2: log singularities of X, jumps pi of Y),
# the zeros of f on the line (sign changes of Lambda(1/2 + it), found as in bsy_levy.py) and the heights of its zeros off the line (jumps 2 pi of Y; below height 300 only), for theta = 3/4 a
# grid of step 1 and the heights of the zeros with beta > 3/4. Beyond, midpoint cells of step 0.5 (0.25 for the finer mesh) with the exact Poisson mass of every cell. Quality of the
# weights: |sum of weights + m_out - 1| is logged (mass defect). Y is computed by arg_ray (bsy_zdh.py): 24 levels from Re s = 2 to Re s = theta, bisected where the increment of the
# argument exceeds 0.9. Beyond the window the exact mass is m_out = (2/pi) arctan(a/T); r_j = m_out times the unweighted mean of X^j (Y^j, XY) over the outer half of the midpoint
# region (the Hoelder bound is replaced by this estimate, flagged as such); the error column is the larger of |coarse - fine| of the corrected quantity and half of the size of the
# correction (|r_2| + 2|I_1 r_1| + r_1^2 for variances). Evaluation: Euler and Maclaurin in double precision (bsy_zdh.py, relative error 1e-13 against mpmath at 30 digits, checked in the log
# at six points; second parameter set nf = 0.45, M = 36 and L = 37 levels on 600 random nodes of every row: the largest differences of X and Y are logged).
# Output: data/bsy_cm.json, the log data/bsy_cm.log. Usage: python3 -I bsy_cm.py [--quick] [--out FILE] [--procs P]; --quick (about 25 s, one process) writes to /tmp unless
# --out is given; the full case (about 5e6 evaluations of f, up to height 1.4e4) runs on the cluster:
#   bash misc/tools/hopper_run.sh -c 16 -t 150 -m 32G -g bsy_cm.json -g bsy_cm.log code/bsy_cm.py code/bsy_zdh.py code/bsy_levy.py data/zeros_1700.txt -- --out bsy_cm.json
# Option --rows "zeta:0.5:5000" (kind:theta:t0, comma separated) computes other rows than the six of the table (task N6 of misc/reviews/nick-b4/reconcile-M3.md: the row t0 = 5000 of
# Table bsy:tab:cm); then --out is required (the default file data/bsy_cm.json is never overwritten), the ordinates of zeta in the window [t0 - W2, t0 + W2] are computed with
# Arb (python-flint acb.zeta_zeros; mpmath.zetazero where flint is missing) when data/zeros_1700.txt does not reach it, and the evaluator is also compared with mpmath at the
# heights t0 and t0 + W3. The six rows of the table are unchanged by the option (the default run and its --quick check give the same numbers). The row t0 = 5000 on the cluster:
#   bash misc/tools/hopper_run.sh -c 16 -t 300 -m 32G -n bsy_cm_5000 -g bsy_cm_5000.json -g bsy_cm_5000.log code/bsy_cm.py code/bsy_zdh.py code/bsy_levy.py data/zeros_1700.txt -- --rows zeta:0.5:5000 --out bsy_cm_5000.json
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse
import json
import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_zdh as Z

HTS = (0.125, 0.0625)
ROWS = [("zeta", 0.5, 0.0), ("zeta", 0.5, 85.7), ("zeta", 0.5, 1000.0), ("zeta", 0.75, 85.7), ("DH", 0.5, 85.7), ("DH", 0.75, 85.7)]


def ts_rule(h):
    K = int(math.ceil(math.asinh(2 / math.pi * math.atanh(1 - 1e-13)) / h))
    k = np.arange(-K, K + 1)
    u = math.pi / 2 * np.sinh(k * h)
    return np.tanh(u), h * math.pi / 2 * np.cosh(k * h) / np.cosh(u) ** 2


def nodes(t0, a, W2, W3, step3, breaks, hts, h3scale):
    """Panels between the breakpoints (zeros, jump heights) on [t0-W2, t0+W2] with a tanh-sinh rule of step hts (log singularities at the ends are integrated to
    1e-12), weights = Poisson density times rule weight; beyond W2 up to W3 the midpoint rule of step step3*h3scale with the exact mass of each cell."""
    lo, hi = t0 - W2, t0 + W2
    br = np.unique(np.concatenate([[lo, hi], breaks[(breaks > lo) & (breaks < hi)]]))
    x, wx = ts_rule(hts)
    ts, ws = [], []
    for p, q in zip(br[:-1], br[1:]):
        if q - p < 1e-9:
            continue
        t = (p + q) / 2 + (q - p) / 2 * x
        d = t - t0
        ts.append(t)
        ws.append(wx * (q - p) / 2 * a / (math.pi * (a * a + d * d)))
    n = int(round((W3 - W2) / (step3 * h3scale)))
    h = (W3 - W2) / n
    d = W2 + (np.arange(n) + 0.5) * h
    w = (np.arctan((d + h / 2) / a) - np.arctan((d - h / 2) / a)) / math.pi
    ts += [t0 - d, t0 + d]
    ws += [w, w]
    return np.concatenate(ts), np.concatenate(ws)


def work(args):
    kind, theta, tau, L, kw = args
    tau = np.where(np.abs(tau) < 1e-6, 1e-6, tau)  # a node exactly at tau = 0 hits s = 1 (removable for f, a cancellation in the Hurwitz sums)
    k = "zeta" if kind == "zeta" else "dh"
    Y, h = Z.arg_ray(k, theta, tau, theta, L=L, **kw)
    return np.log(np.abs(h)), Y


def evaluate(pool, kind, theta, tau, L=24, kw=None, chunk=500):
    kw = kw or {}
    order = np.argsort(np.abs(tau))
    jobs = [(kind, theta, tau[order[i:i + chunk]], L, kw) for i in range(0, len(tau), chunk)]
    res = list(pool.imap(work, jobs)) if pool else [work(j) for j in jobs]
    X = np.concatenate([r[0] for r in res])
    Y = np.concatenate([r[1] for r in res])
    Xo = np.empty_like(X)
    Yo = np.empty_like(Y)
    Xo[order] = X
    Yo[order] = Y
    return Xo, Yo


def breaks_for(kind, theta, t0, W2, zlist, dhz):
    """Breakpoints: zeta, theta = 1/2: the ordinates (both signs) and 0; DH, theta = 1/2: its line zeros and the heights of its zeros off the line (jumps of Y);
    theta > 1/2: the heights of the zeros with beta > theta (jumps of Y) and a uniform grid of step 1 (no singularities, keeps the panels short)."""
    if kind == "zeta" and theta == 0.5:
        return np.concatenate([zlist, -zlist, [0.0]])
    if kind == "DH" and theta == 0.5:
        return np.concatenate([dhz["line"], dhz["off"][:, 1], -dhz["line"], -dhz["off"][:, 1]])
    base = np.arange(math.floor(t0 - W2), math.ceil(t0 + W2) + 1.0, 1.0)
    if kind == "DH":
        base = np.concatenate([base, dhz["off"][dhz["off"][:, 0] > theta, 1]])
    return base


def extend_zeros(zlist, t0s, W2, log):
    """Ordinates of zeta in [t0 - W2 - 5, t0 + W2 + 5] for each t0 beyond the file: Arb (python-flint) in blocks from the index N(t0 - W2 - 5) + 1, mpmath.zetazero if flint is
    missing. The count N(T) is mpmath.nzeros; the first ordinate found is checked to lie above T and the previous one (index n0 - 1) below it."""
    import mpmath as mp
    out = [zlist]
    for t0 in t0s:
        lo, hi = t0 - W2 - 5.0, t0 + W2 + 5.0
        if hi <= zlist[-1]:
            continue
        n0 = int(mp.nzeros(lo)) + 1
        zs = []
        try:
            from flint import acb, ctx
            ctx.prec = 64
            n = n0
            while not zs or zs[-1] <= hi:
                zs += [float(z.imag.mid()) for z in acb.zeta_zeros(n, 400)]
                n += 400
            how = "Arb (python-flint, 64 bits)"
        except ImportError:
            mp.mp.dps = 15
            n = n0
            while not zs or zs[-1] <= hi:
                zs.append(float(mp.zetazero(n).imag))
                n += 1
            how = "mpmath.zetazero"
        zs = np.array([z for z in zs if z <= hi])
        log(f"  ordinates of zeta in [{lo:g}, {hi:g}]: {len(zs)} from index {n0} ({how}), first {zs[0]:.6f}, last {zs[-1]:.6f}; mpmath.nzeros(first) - nzeros(last) check: "
            f"{int(mp.nzeros(zs[-1] + 1e-3)) - int(mp.nzeros(lo)) } expected {len(zs)}")
        out.append(zs)
    z = np.unique(np.concatenate(out))
    return z


def row(pool, kind, theta, t0, W2, W3, step3, zlist, dhz, log):
    a = 1 - theta
    k = "zeta" if kind == "zeta" else "dh"
    if kind == "zeta" and t0 == 0.0 and theta == 0.5:
        x0, y0 = 0.0, 0.0
    else:
        Yr, h1 = Z.arg_ray(k, theta, [t0], 1.0)
        x0, y0 = float(np.log(abs(h1[0]))), float(Yr[0])
    out = {"kind": kind, "theta": theta, "t0": t0, "x0": x0, "y0": y0}
    brk = breaks_for(kind, theta, t0, W2, zlist, dhz)
    out["n_breaks"] = int(len(brk))
    for name, hts, sc in (("coarse", HTS[0], 1.0), ("fine", HTS[1], 0.5)):
        tau, w = nodes(t0, a, W2, W3, step3, brk, hts, sc)
        X, Y = evaluate(pool, kind, theta, tau)
        m_out = 2 / math.pi * math.atan(a / W3)
        tcheck = abs(w.sum() + m_out - 1)
        far = np.abs(tau - t0) > (W3 + W2) / 2
        r1, r2 = m_out * X[far].mean(), m_out * (X[far] ** 2).mean()
        s1, s2 = m_out * Y[far].mean(), m_out * (Y[far] ** 2).mean()
        rxy = m_out * (X[far] * Y[far]).mean()
        I = {"X1": (w * X).sum(), "X2": (w * X * X).sum(), "Y1": (w * Y).sum(), "Y2": (w * Y * Y).sum(), "XY": (w * X * Y).sum()}
        EX, EY = I["X1"] + r1, I["Y1"] + s1
        res = {"nodes": int(len(tau)), "mass_defect": tcheck, "m_out": m_out, "I": I, "r": [r1, r2, s1, s2, rxy],
               "EX-x0": EX - x0, "EY-y0": EY - y0, "VarX": I["X2"] + r2 - EX ** 2, "VarY": I["Y2"] + s2 - EY ** 2,
               "Cov": I["XY"] + rxy - EX * EY,
               "VarX_nocorr": I["X2"] - I["X1"] ** 2, "VarY_nocorr": I["Y2"] - I["Y1"] ** 2, "Cov_nocorr": I["XY"] - I["X1"] * I["Y1"],
               "corr": {"EX": abs(r1), "EY": abs(s1), "VarX": abs(r2) + 2 * abs(I["X1"] * r1) + r1 ** 2,
                        "VarY": abs(s2) + 2 * abs(I["Y1"] * s1) + s1 ** 2, "Cov": abs(rxy) + abs(I["X1"] * s1) + abs(I["Y1"] * r1) + abs(r1 * s1)}}
        out[name] = res
        log(f"  {kind} theta={theta} t0={t0} {name}: nodes {len(tau)}, mass defect {tcheck:.1e}, m_out {m_out:.2e}, E X - x0 {EX-x0:+.6f}, E Y - y0 {EY-y0:+.6f}, "
            f"Var X {res['VarX']:.6f}, Var Y {res['VarY']:.6f}, Cov {res['Cov']:+.6f} (uncorrected Var X {res['VarX_nocorr']:.6f}, Var Y {res['VarY_nocorr']:.6f})")
    out["x0y0"] = [x0, y0]
    err = {}
    for q in ("EX-x0", "EY-y0", "VarX", "VarY", "Cov"):
        key = {"EX-x0": "EX", "EY-y0": "EY"}.get(q, q)
        err[q] = max(abs(out["coarse"][q] - out["fine"][q]), 0.5 * out["fine"]["corr"][key])
    out["error"] = err
    log(f"  error column (max of mesh difference and half the tail correction): " + ", ".join(f"{q} {v:.1e}" for q, v in err.items()))
    # second parameter set on random nodes of the fine mesh
    rng = np.random.default_rng(1)
    tau, _ = nodes(t0, a, W2, W3, step3, brk, HTS[1], 0.5)
    sub = rng.choice(len(tau), size=min(600, len(tau)), replace=False)
    X1, Y1 = evaluate(None, kind, theta, tau[sub])
    X2, Y2 = evaluate(None, kind, theta, tau[sub], L=37, kw={"nf": 0.45, "M": 36})
    ok = (X1 > -18.0) & (X2 > -18.0)  # nodes within 1e-8 of a zero (|h| < 1e-8, weight below 1e-8) carry rounding noise in log|h|
    out["check_params"] = {"max_dX": float(np.abs(X1 - X2)[ok].max()), "max_dY": float(np.abs(Y1 - Y2)[ok].max()), "excluded_nodes": int((~ok).sum()),
                           "max_dX_excluded": float(np.abs(X1 - X2)[~ok].max()) if (~ok).any() else 0.0}
    log(f"  second parameters (nf 0.45, M 36, 37 levels) on {len(sub)} nodes ({int((~ok).sum())} nodes with |h| < 1e-8 left out, there |dX| up to {out['check_params']['max_dX_excluded']:.1e}): "
        f"max |dX| {out['check_params']['max_dX']:.2e}, max |dY| {out['check_params']['max_dY']:.2e}")
    return out


def main():
    global HTS
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--procs", type=int, default=0)
    ap.add_argument("--zeros", default=None)
    ap.add_argument("--rows", default=None)
    A = ap.parse_args()
    global ROWS
    if A.rows:
        if not A.out:
            ap.error("--rows needs --out (data/bsy_cm.json is the table's file)")
        ROWS = [(r.split(":")[0], float(r.split(":")[1]), float(r.split(":")[2])) for r in A.rows.split(",")]
    if A.quick:
        W2, W3, step3 = 100.0, 400.0, 2.0
        HTS = (0.25, 0.125)
        outp = A.out or "/tmp/bsy_cm_quick.json"
    else:
        W2, W3, step3 = 1000.0, 13000.0, 0.5
        outp = A.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "bsy_cm.json")
    logp = os.path.splitext(outp)[0] + ".log"
    lf = open(logp, "w")

    def log(m):
        print(m, flush=True)
        lf.write(m + "\n")
        lf.flush()

    t00 = time.time()
    log(f"bsy_cm.py {'quick' if A.quick else 'full'}: panels on |tau - t0| <= {W2:g} (tanh-sinh, steps {HTS[0]} and {HTS[1]}), midpoint cells to {W3:g} (steps {step3} and {step3 / 2})")
    here = os.path.dirname(os.path.abspath(__file__))
    zp = A.zeros or next(q for q in (os.path.join(here, "zeros_1700.txt"), os.path.join(here, "..", "data", "zeros_1700.txt")) if os.path.exists(q))
    zlist = np.array([float(x) for x in open(zp).read().split()])
    if A.rows:
        zlist = extend_zeros(zlist, [t0 for k, th, t0 in ROWS if k == "zeta"], W2, log)
    log(f"  ordinates of zeta: {len(zlist)} from {os.path.basename(zp)} (largest {zlist[-1]:.3f}); needed up to {1000 + W2 + 100:g}")
    import bsy_levy as LV
    Tz = 85.7 + W2 + 5.0
    if A.rows and not any(k == "DH" for k, th, t0 in ROWS):
        dhz = None
        Tz = 0.0
    if Tz > 0:
        dl, _ = LV.line_zeros(Tz, 0.005)
        dl2, _ = LV.line_zeros(Tz, 0.0071)
        do, _ = LV.offline_zeros(300.0)
        dhz = {"line": dl, "off": np.array([[z.real, z.imag] for z in do])}
        log(f"  DH zeros: {len(dl)} sign changes of Lambda(1/2 + it) on (0, {Tz:g}] (second grid: {len(dl2)}), {len(do)} zeros off the line below 300 (beta > 1/2): " + ", ".join(f"{z.real:.4f}+{z.imag:.3f}i" for z in do))
    pts = [0.5 + 85.7j, 0.75 + 400j, 0.6 + 13.2j]
    if A.rows:
        pts = [complex(0.5, t0) for k, th, t0 in ROWS] + [complex(0.5, t0 + W3) for k, th, t0 in ROWS]
    for s in pts:
        for kind in ("zeta", "dh"):
            v = Z.hfun(kind, 0.5, [s])[0]
            m = Z.mp_h(kind, 0.5, s)
            log(f"  evaluator vs mpmath (30 digits) {kind} at {s}: relative difference {abs(v - m) / abs(m):.1e}")
    procs = A.procs or (1 if A.quick else (os.cpu_count() or 1))
    pool = None
    if procs > 1:
        import multiprocessing as mpc
        pool = mpc.Pool(procs)
    res = []
    for kind, theta, t0 in ROWS:
        t1 = time.time()
        log(f"row {kind} theta={theta} t0={t0}")
        res.append(row(pool, kind, theta, t0, W2, W3, step3, zlist, dhz, log))
        log(f"  row time {time.time() - t1:.0f} s")
    json.dump({"W2": W2, "W3": W3, "step3": step3, "quick": A.quick, "rows": res}, open(outp, "w"), indent=1)
    log(f"done in {time.time() - t00:.0f} s; wrote {outp}")


if __name__ == "__main__":
    main()
