# bsy_yks.py: Table bsy:tab:yks (root mean square errors of plain Monte Carlo, the sorted trapezoid rule of Yakowitz, Krimmel and Szidarovszky (YKS) and YKS with
# the Euler product control, for the Cauchy average of log|zeta(theta + i tau)|, theta = 1/2 and 7/8, t0 = 0, N = 500, 2000, 8000, 16000, 40 replications), its
# fitted slopes, the gain factors of the prose after the table, and the synthetic slopes of the text after Proposition bsy:prop:rate (cos(tau log 2) and the Selberg
# smoothed polynomial Re D_1000). Task C7 of reconcile-M.md (fresh run from the definitions in the text only). Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Definitions (from the chapter). P = Cauchy law centred at t0 with scale a = 1 - theta; tau(u) = t0 + a tan(pi (u - 1/2)); E_P[f] = int_0^1 f(tau(u)) du. Plain: the
# mean of f(tau(u_i)) over N iid uniforms. YKS: T_N = sum_{i=0}^N (u_(i+1) - u_(i)) (g(u_(i)) + g(u_(i+1)))/2 with u_(0) = 0, u_(N+1) = 1, g = f o tau, g(0) = g(1) = 0.
# Control: D_M(s) = sum_{n <= M} b_n n^-s, b_n = Lambda(n)/log n (1 - log n/log M) (Lambda the von Mangoldt function; the text writes log M, not log(M + 1)), the
# trapezoid rule is applied to g - Re D_M o tau - ell o tau with ell(tau) = log|(theta + 1 + i tau)/(theta - 1 + i tau)|, and the exact means E_P[Re D_M] = Re D_M(1 + i t0),
# E_P[ell] = log|(2 + i t0)/(2 - 2 theta + i t0)| (Lemma bsy:lem:cv) are added back; M = 100 at theta = 1/2 and M = 1000 at theta = 7/8. The true value E_P[f] is
# log|h_theta(1 + i t0)| (Theorem bsy:thm:bsy, no zero with beta > theta, i.e. RH at 1/2 and QRH(7/8)): 0 at theta = 1/2 and log 4 at theta = 7/8 for t0 = 0 (the true
# value differs from these by less than 7.4e-13 and 4.6e-14 by bsy_tailbound.py). Error = estimate minus this value; rms over the replications.
# Design. For each replication one array of 16000 uniforms is drawn and the sample for N is its first N points (an iid sample of size N; the errors for different N
# are correlated within a replication, each marginal is exact); 40 replications, seeds 1000 + r. Slopes: least squares of log rms on log N over the four values of N.
# The synthetic slopes use N = 500, 2000, 8000, 32000, 128000 with the same 40 replications: f = cos(tau log 2) (exact mean 2^-a cos(0) at t0 = 0) and f = Re D_1000
# (exact mean Re D_1000(1 + i t0)), at theta = 1/2 and 7/8, plain and YKS.
# log|zeta| is evaluated by python-flint acb.zeta at 53 bits (Hopper) or mpmath (quick, laptop). Error control: on 2000 random nodes of the first replication the value is
# recomputed at 106 bits and the largest difference is logged; the replications' spread gives the sampling error of every rms (relative 1/sqrt(2 R) = 11%).
# Column mean_ctrl (added 25 October 2026): the same control with coefficient one applied to the plain sample mean, mean(g - Re D_M o tau - ell o tau) + exact means - truth.
# Output: data/bsy_yks.json (raw replicate errors included) and the log bsy_yks.log. Usage: python3 bsy_yks.py [--quick] [--out FILE] [--procs P]; --quick: 4 replications,
# N = 100, 200, 400, 800, mpmath, a few seconds; full on the cluster: bash misc/tools/hopper_run.sh -c 8 -t 120 -g bsy_yks.json -g bsy_yks.log code/bsy_yks.py -- --out bsy_yks.json
import argparse, json, math, os, time
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--procs", type=int, default=0); ap.add_argument("--reps", type=int, default=0)
args = ap.parse_args()
here = Path(__file__).resolve().parent
out = Path(args.out) if args.out else (Path("/tmp/bsy_yks_quick.json") if args.quick else here.parent / "data" / "bsy_yks.json")
try:
    import flint
    from flint import acb, ctx
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False
    import mpmath as mp
    mp.mp.dps = 15

def logabs_zeta(sigma, taus, prec=53):
    taus = np.abs(np.asarray(taus, dtype=float))   # |zeta(sigma - i t)| = |zeta(sigma + i t)|
    res = np.empty(len(taus))
    if HAVE_FLINT:
        ctx.prec = prec
        for i, t in enumerate(taus):
            z = acb(sigma, float(t)).zeta()
            re = float(z.real.mid()); im = float(z.imag.mid())
            res[i] = 0.5*math.log(re*re + im*im)
    else:
        for i, t in enumerate(taus):
            z = complex(mp.zeta(mp.mpc(sigma, float(t))))
            res[i] = math.log(abs(z))
    return res

def primes_upto(n):
    sieve = np.ones(n + 1, bool); sieve[:2] = False
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]: sieve[p*p::p] = False
    return np.nonzero(sieve)[0]
def dirichlet_terms(M):
    ns, bs = [], []
    for p in primes_upto(M):
        q = int(p)
        while q <= M:
            ns.append(q); bs.append(math.log(p)/math.log(q)*(1 - math.log(q)/math.log(M))); q *= int(p)
    return np.array(ns, float), np.array(bs)
def reD(theta, taus, M):
    ns, bs = dirichlet_terms(M)
    taus = np.asarray(taus, float)
    out = np.zeros(len(taus))
    for lo in range(0, len(taus), 4000):
        t = taus[lo:lo+4000, None]
        out[lo:lo+4000] = (bs[None, :]*ns[None, :]**(-theta)*np.cos(t*np.log(ns)[None, :])).sum(axis=1)
    return out
def reD_exact_mean(t0, M):
    ns, bs = dirichlet_terms(M)
    return float((bs/ns*np.cos(t0*np.log(ns))).sum())
def ell(theta, taus):
    taus = np.asarray(taus, float)
    return 0.5*np.log(((theta + 1)**2 + taus**2)/((theta - 1)**2 + taus**2))
def ell_mean(theta, t0):
    return 0.5*math.log((4 + t0**2)/((2 - 2*theta)**2 + t0**2))

def yks_rule(u, g):
    o = np.argsort(u); us = np.concatenate(([0.0], u[o], [1.0])); gs = np.concatenate(([0.0], g[o], [0.0]))
    return float(np.sum(np.diff(us)*(gs[:-1] + gs[1:])/2))

THETAS = [(0.5, 100), (0.875, 1000)]
NS = [100, 200, 400, 800] if args.quick else [500, 2000, 8000, 16000]
NSYN = [100, 200, 400, 800] if args.quick else [500, 2000, 8000, 32000, 128000]
REPS = args.reps or (4 if args.quick else 40)   # --reps 400: robustness run (same seeds, the first 40 replications are those of the default run)
t0 = 0.0

def task(args_):
    kind, theta, M, r = args_
    a = 1 - theta
    nmax = max(NS) if kind == "zeta" else max(NSYN)
    rng = np.random.default_rng(1000 + r)
    u = rng.random(nmax)
    tau = t0 + a*np.tan(np.pi*(u - 0.5))
    if kind == "zeta":
        f = logabs_zeta(theta, tau)
        truth = 0.0 if theta == 0.5 else math.log(1/(2 - 2*theta))
        ctrl = reD(theta, tau, M) + ell(theta, tau)
        cmean = reD_exact_mean(t0, M) + ell_mean(theta, t0)
        sizes = NS
    elif kind == "cos":
        f = np.cos(tau*math.log(2)); truth = 2**(-a); ctrl = None; sizes = NSYN
    else:  # selberg polynomial
        f = reD(theta, tau, 1000); truth = reD_exact_mean(t0, 1000); ctrl = None; sizes = NSYN
    res = {}
    for N in sizes:
        uu = u[:N]
        plain = float(np.mean(f[:N])) - truth
        y = yks_rule(uu, f[:N]) - truth
        d = {"plain": plain, "yks": y}
        if ctrl is not None:
            d["yks_ctrl"] = yks_rule(uu, f[:N] - ctrl[:N]) + cmean - truth
            d["mean_ctrl"] = float(np.mean(f[:N] - ctrl[:N])) + cmean - truth   # control with coefficient one applied to the sample mean (final review 35b, M2)
        res[N] = d
    return (kind, theta, r, res)

if __name__ == "__main__":
    log = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
    P("flint" if HAVE_FLINT else "mpmath", "replications", REPS, "N", NS)
    jobs = [("zeta", th, M, r) for th, M in THETAS for r in range(REPS)] + [(k, th, 0, r) for k in ("cos", "sel") for th, M in THETAS for r in range(REPS)]
    procs = args.procs or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1)
    t_start = time.time()
    if procs > 1:
        import multiprocessing as mpc
        with mpc.Pool(procs) as pool: results = pool.map(task, jobs, chunksize=1)
    else:
        results = [task(j) for j in jobs]
    P("compute time", round(time.time() - t_start, 1), "s on", procs, "processes")
    # precision check
    rng = np.random.default_rng(1000); u = rng.random(2000); tau = 0.5*np.tan(np.pi*(u - 0.5))
    if HAVE_FLINT:
        d = np.max(np.abs(logabs_zeta(0.5, tau, 53) - logabs_zeta(0.5, tau, 106)))
        P("precision check: max |log|zeta| (53 bits) - (106 bits)| on 2000 nodes at theta = 1/2:", f"{d:.2e}")
    store = {}
    for kind, theta, r, res in results:
        for N, d in res.items():
            for m, v in d.items(): store.setdefault((kind, theta, m, N), []).append(v)
    summary = {}
    for kind, label in (("zeta", "zeta"), ("cos", "cos(tau log 2)"), ("sel", "Re D_1000")):
        for theta, M in THETAS:
            sizes = NS if kind == "zeta" else NSYN
            rms = {}
            for m in ("plain", "yks", "yks_ctrl", "mean_ctrl"):
                if (kind, theta, m, sizes[0]) not in store: continue
                rms[m] = [math.sqrt(np.mean(np.array(store[(kind, theta, m, N)])**2)) for N in sizes]
            sl = {m: float(np.polyfit(np.log(sizes), np.log(v), 1)[0]) for m, v in rms.items()}
            P(f"{label} theta={theta}: N={sizes}")
            for m, v in rms.items(): P(f"   {m:9s} rms " + " ".join(f"{x:.2e}" for x in v) + f"   slope {sl[m]:.3f}")
            if "yks_ctrl" in rms:
                P("   ratio plain/mean_ctrl " + " ".join(f"{p/y:.2f}" for p, y in zip(rms['plain'], rms['mean_ctrl'])))
                P("   ratio plain/YKS      " + " ".join(f"{p/y:.1f}" for p, y in zip(rms['plain'], rms['yks'])))
                P("   ratio plain/YKS+ctrl " + " ".join(f"{p/y:.1f}" for p, y in zip(rms['plain'], rms['yks_ctrl'])))
                P("   ratio YKS/ctrl       " + " ".join(f"{p/y:.1f}" for p, y in zip(rms['yks'], rms['yks_ctrl'])))
                P("   ratio plain/ctrl squared (equivalent plain samples) " + " ".join(f"{(p/y)**2:.0f}" for p, y in zip(rms['plain'], rms['yks_ctrl'])))
            summary[f"{kind}_{theta}"] = {"N": sizes, "rms": rms, "slope": sl}
    raw = {f"{k[0]}|{k[1]}|{k[2]}|{k[3]}": v for k, v in store.items()}
    out.write_text(json.dumps({"summary": summary, "raw_errors": raw, "reps": REPS}, indent=0))
    Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
