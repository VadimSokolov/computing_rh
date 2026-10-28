# bsy_mc.py: the two Monte Carlo columns of Table bsy:tab:defects (plain Monte Carlo and vertical likelihood Monte Carlo) for every row of the table, Chapter 35
# (Balazard, Saias and Yor meet Thorin), group BSY-D (final review 35b, item M1). Added for the book in October 2026, written from the mathematics of the text of
# Section bsy:sec:mc (the collaborator's scripts were not read, run or copied).
# Target. f(tau) = log|F(theta + i tau)|, F = zeta or the Davenport and Heilbronn function (DH, code/bsy_zdh.py), P = Cauchy law of scale a = 1 - theta centred at t0,
# tau = t0 + a tan(pi (u - 1/2)). The defect is D_theta(t0) = E_P[f] - log|h_theta(1 + i t0)| (h = zeta (s - 1)/(s - 2 theta + 1), h(1) = 1/(2 - 2 theta); h = F for DH).
# Plain: n iid exact Cauchy draws (n = 20000 for zeta, 4000 for DH, the caption), estimate mean f, standard error sd/sqrt(n).
# Vertical (Section bsy:sec:mc, "Vertical likelihood Monte Carlo"): P_w(dtau) proportional to P(dtau) w(f(tau)), w(f) = |f - c| + eps, c a pilot estimate (the mean of
# n/20 plain draws, separate from the plain run), eps = 0.05 sd(pilot) (the text says only "regularised": these two choices and the three below are the script's).
# Metropolis chain mixing, with probability 1/2 each, (i) the Cauchy independence proposal tau' ~ P (acceptance min(1, w(f(tau'))/w(f(tau)))) and (ii) a local Gaussian
# step tau' = tau + (a/2) N(0,1) (acceptance min(1, p(tau') w(f(tau'))/(p(tau) w(f(tau)))), p the Cauchy density); a mixture of two P_w reversible kernels. Start: an
# independence draw; burn-in n/20 steps; recorded steps n - n/20 - n/20 (the total of f evaluations, pilot and burn-in included, is n). Estimate: self normalised weights
# 1/w, mu_hat = sum f_i/w_i / sum 1/w_i. Standard error: batch means (B = 50 equal batches) of Y_i = (f_i - mu_hat)/w_i and Z_i = 1/w_i, se = sd_B(Y)/(sqrt(B) mean(Z)).
# Depth: the smallest f visited by each sampler (plain: its n draws; vertical: the recorded chain), and the number of draws with f < -1, -2 (vertical: unweighted
# visits of the chain, which is the point of the weighting).
# Seeds. Table run: np.random.SeedSequence([20261025, row, 0 for the plain run, 1 for the pilot, 2 for the chain]); calibration replicates r = 0..R-1: the same with the
# extra entry r + 1. All stored in the JSON. Error control: the evaluator at (nf, M) = (0.3, 28) is compared with (0.45, 36) on the 2000 first plain draws of every row
# (the maximum difference is logged); the second check is the calibration of the standard errors over R replicates (mean and sd of z = (estimate - quadrature)/se).
# Reference value of each row: the quadrature column of data/bsy_defects.json (W = 200 if present, else W = 60); the windowed quadrature is within its tail bound of the
# line integral, far below the Monte Carlo errors (so the comparison in standard errors is that with the exact defect, up to the DH prediction rows noted in the log).
# Output: data/bsy_mc.json and data/bsy_mc.log. Usage: python3 bsy_mc.py [--quick] [--out FILE] [--procs P] [--reps R]; --quick: n/10, no calibration, laptop (30 s);
# full on the cluster: bash misc/tools/hopper_run.sh -c 8 -t 120 -g bsy_mc.json -g bsy_mc.log code/bsy_mc.py data/bsy_defects.json bsy_zdh.py -- --out bsy_mc.json --reps 100
import argparse, json, math, os, sys, time
from pathlib import Path
import numpy as np

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))
import bsy_zdh as Z

def _block_scaled(s, a, N, M, deriv):
    """Euler and Maclaurin block of bsy_zdh._block with the products (s)_(2j-1) (N+a)^(-2j+1) updated together: the shared evaluator overflows (s)_(2j-1) for |s| > 4e5
    (nan at |tau| = 4.5e7, a Cauchy draw that occurs in about 1 replicate in 100); the sum is the same formula."""
    assert not deriv
    k = np.arange(N, dtype=float) + a
    lk = np.log(k)
    S = np.exp(-np.multiply.outer(s, lk)).sum(1)
    X = N + a; lX = math.log(X); Xs = np.exp(-s * lX)
    T = X * Xs / (s - 1) + Xs / 2
    pq = s / X
    for j in range(1, M + 1):
        T = T + Z._C[j] * pq * Xs
        pq = pq * (s + 2 * j - 1) * (s + 2 * j) / (X * X)
    return S + T, None

Z._block = _block_scaled

ROWS = [("zeta", 0.5, 0.0), ("zeta", 0.875, 0.0), ("zeta", 0.5, 85.7), ("zeta", 0.75, 85.7), ("zeta", 0.875, 85.7),
        ("dh", 0.5, 85.7), ("dh", 0.75, 85.7), ("dh", 0.75, 100.0), ("dh", 0.875, 85.7)]
NB = 50

def F(kind, theta, tau, **kw):
    s = theta + 1j * np.atleast_1d(np.asarray(tau, float))
    v = Z.zeta(s, **kw) if kind == "zeta" else Z.dh(s, **kw)
    return np.log(np.abs(v))

def logh1(kind, theta, t0):
    if kind == "zeta":
        if t0 == 0.0:
            return -math.log(2 - 2 * theta)
        return float(np.log(np.abs(Z.hfun("zeta", theta, np.array([1 + 1j * t0]))[0])))
    return float(F("dh", 1.0, t0)[0])   # |f(1 + i t0)|: theta = 1 gives s = 1 + i t0

def cauchy(t0, a, u):
    return t0 + a * np.tan(np.pi * (u - 0.5))

def dens(t0, a, t):
    return a / (math.pi * (a * a + (t - t0) ** 2))

def plain(kind, theta, t0, n, rng):
    a = 1 - theta
    tau = cauchy(t0, a, rng.random(n))
    f = F(kind, theta, tau)
    assert np.isfinite(f).all(), "non finite f in the plain run"
    return {"est": float(f.mean()), "se": float(f.std(ddof=1) / math.sqrt(n)), "min_f": float(f.min()), "n_below_-1": int((f < -1).sum()), "n_below_-2": int((f < -2).sum()), "sd": float(f.std(ddof=1))}

def vertical(kind, theta, t0, n, rng_pilot, rng_chain):
    a = 1 - theta
    npilot = n // 20; burn = n // 20; steps = n - npilot - burn
    fp = F(kind, theta, cauchy(t0, a, rng_pilot.random(npilot)))
    c = float(fp.mean()); eps = 0.05 * float(fp.std(ddof=1))
    w = lambda f: np.abs(f - c) + eps
    tot = burn + steps
    isind = rng_chain.random(tot + 1) < 0.5
    tau_ind = cauchy(t0, a, rng_chain.random(tot + 1))
    f_ind = F(kind, theta, tau_ind)                 # evaluated for all (only the independence ones are used; the others are discarded, not counted)
    xi = rng_chain.standard_normal(tot + 1) * (a / 2)
    ua = np.log(rng_chain.random(tot + 1))
    x = tau_ind[0]; fx = f_ind[0]; wx = float(w(fx)); px = dens(t0, a, x)
    assert np.isfinite(f_ind).all(), "non finite f at an independence proposal"
    fs = np.empty(steps); nacc = 0; nev = 0
    for i in range(1, tot + 1):
        if isind[i]:
            y = tau_ind[i]; fy = f_ind[i]; wy = float(w(fy))
            lr = math.log(wy) - math.log(wx); py = dens(t0, a, y)
        else:
            y = x + xi[i]; fy = float(F(kind, theta, y)[0]); assert math.isfinite(fy); wy = float(w(fy)); py = dens(t0, a, y)
            lr = math.log(py * wy) - math.log(px * wx)
        nev += 1
        if ua[i] < lr:
            x, fx, wx, px = y, fy, wy, py
            nacc += 1 if i > burn else 0
        if i > burn:
            fs[i - burn - 1] = fx
    ws = w(fs); iw = 1.0 / ws
    mu = float((fs * iw).sum() / iw.sum())
    Y = (fs - mu) * iw
    m = steps // NB
    Yb = Y[:m * NB].reshape(NB, m).mean(1); Zb = iw[:m * NB].reshape(NB, m).mean(1)
    se = float(Yb.std(ddof=1) / (math.sqrt(NB) * Zb.mean()))
    return {"est": mu, "se": se, "c": c, "eps": eps, "min_f": float(fs.min()), "n_below_-1": int((fs < -1).sum()), "n_below_-2": int((fs < -2).sum()),
            "acc_rate": nacc / steps, "unweighted_mean": float(fs.mean()), "steps": steps}

def one(args_):
    ri, rep, n, check = args_
    kind, theta, t0 = ROWS[ri]
    tag = 0 if rep is None else rep + 1
    ss = lambda k: np.random.default_rng(np.random.SeedSequence([20261025, ri, k, tag]))
    pl = plain(kind, theta, t0, n, ss(0))
    vt = vertical(kind, theta, t0, n, ss(1), ss(2))
    out = {"row": ri, "rep": rep, "plain": pl, "vertical": vt}
    if check:
        tau = cauchy(t0, 1 - theta, ss(0).random(min(n, 2000)))
        out["eval_check"] = float(np.max(np.abs(F(kind, theta, tau) - F(kind, theta, tau, nf=0.45, M=36))))
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--procs", type=int, default=0); ap.add_argument("--reps", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out) if args.out else (Path("/tmp/bsy_mc_quick.json") if args.quick else here.parent / "data" / "bsy_mc.json")
    qf = here.parent / "data" / "bsy_defects.json"
    if not qf.exists(): qf = here / "bsy_defects.json"
    if not qf.exists(): qf = Path("bsy_defects.json")
    quad = json.load(open(qf))["rows"]
    log = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
    def nof(kind): return (20000 if kind == "zeta" else 4000) // (10 if args.quick else 1)
    reps = 0 if args.quick else args.reps
    jobs = [(ri, None, nof(ROWS[ri][0]), True) for ri in range(len(ROWS))] + [(ri, r, nof(ROWS[ri][0]), False) for ri in range(len(ROWS)) for r in range(reps)]
    procs = args.procs or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1)
    t_start = time.time()
    if procs > 1:
        import multiprocessing as mpc
        with mpc.Pool(procs) as pool: res = pool.map(one, jobs, chunksize=1)
    else:
        res = [one(j) for j in jobs]
    P("compute time", round(time.time() - t_start, 1), "s on", procs, "processes; reps", reps, "; quick", args.quick)
    rows_out = []
    for ri, (kind, theta, t0) in enumerate(ROWS):
        q = next(r for r in quad if r["function"] == kind and r["theta"] == theta and r["t0"] == t0)
        ref = q["W200"]["D_W"] if "W200" in q else q["W60"]["D_W"]
        lh = logh1(kind, theta, t0)
        main = next(r for r in res if r["row"] == ri and r["rep"] is None)
        row = {"function": kind, "theta": theta, "t0": t0, "n": nof(kind), "quad_W200_or_W60": ref, "quad_prediction": q["prediction"], "log_h1": lh, "seed": [20261025, ri, "0 plain, 1 pilot, 2 chain", 0], "eval_check_nf0.3_vs_0.45": main["eval_check"]}
        for m in ("plain", "vertical"):
            r = main[m]; d = r["est"] - lh
            row[m] = {"D": d, "se": r["se"], "z_vs_quadrature": (d - ref) / r["se"], "z_vs_prediction": (d - q["prediction"]) / r["se"], **{k: v for k, v in r.items() if k not in ("est", "se")}}
        cal = [r for r in res if r["row"] == ri and r["rep"] is not None]
        if cal:
            for m in ("plain", "vertical"):
                e = np.array([r[m]["est"] - lh - ref for r in cal]); s = np.array([r[m]["se"] for r in cal]); zz = e / s
                row[m]["calibration"] = {"reps": len(cal), "mean_error": float(e.mean()), "rms_error": float(np.sqrt((e ** 2).mean())), "mean_se": float(s.mean()), "rms_over_mean_se": float(np.sqrt((e ** 2).mean()) / s.mean()),
                                         "z_mean": float(zz.mean()), "z_sd": float(zz.std(ddof=1)), "coverage_2se": float((np.abs(zz) <= 2).mean()), "min_f_median": float(np.median([r[m]["min_f"] for r in cal]))}
        rows_out.append(row)
        P(f"{kind} ({theta}, {t0}) n={row['n']} quad={ref:.4g} pred={q['prediction']:.4g} evalcheck={row['eval_check_nf0.3_vs_0.45']:.1e}")
        for m in ("plain", "vertical"):
            r = row[m]
            P(f"   {m:8s} D = {r['D']:+.4f} +- {r['se']:.4f}  z(quad) = {r['z_vs_quadrature']:+.2f}  min f = {r['min_f']:.2f}  n(f<-1) = {r['n_below_-1']}  n(f<-2) = {r['n_below_-2']}" + (f"  acc {r['acc_rate']:.2f}" if m == "vertical" else ""))
            if "calibration" in r:
                c = r["calibration"]; P(f"      calibration {c['reps']} reps: rms error {c['rms_error']:.4f}, mean se {c['mean_se']:.4f}, z mean {c['z_mean']:+.2f} sd {c['z_sd']:.2f}, coverage(2se) {c['coverage_2se']:.2f}, median min f {c['min_f_median']:.2f}")
    out.write_text(json.dumps({"rows": rows_out, "reps": reps, "seed_base": 20261025}, indent=1))
    Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
