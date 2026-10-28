# bsy_isbound.py: the ceiling factor of Remark bsy:rem:weight (Chapter 35, Balazard, Saias and Yor meet Thorin, section "Vertical likelihood Monte Carlo"), group BSY-G, task N4 of
# misc/reviews/nick-b4/reconcile-M3.md. Claim supported: for independent draws no weight w > 0 gains more than the factor R = Var_P f / (E_P |f - mu|)^2 over plain Monte Carlo, for
# the integrands f = log|h_theta(theta + i tau)| of Table bsy:tab:defects (the nine rows of the table: zeta at (1/2, 0), (7/8, 0), (1/2, 85.7), (7/8, 85.7), (3/4, 85.7) and the
# Davenport and Heilbronn function (DH) at (1/2, 85.7), (3/4, 85.7), (7/8, 85.7), (3/4, 100)); the note's own figures (M3:558-565) are not used.
# Conventions as in bsy_cm.py: h_theta(s) = zeta(s)(s-1)/(s-2 theta+1) (h = f, entire, for DH), a = 1 - theta, P(tau) = a/(pi (a^2 + (tau - t0)^2)) the Cauchy law (the same P as
# bsy_cm.py and bsy_defects.py), f(tau) = log|h_theta(theta + i tau)|, mu = E_P f, Var_P f = E_P f^2 - mu^2, A = E_P |f - mu|, R = Var_P f / A^2 >= 1 (Cauchy and Schwarz).
# Method. The nodes, the weights (Poisson density times tanh-sinh rule weight between the zeros as breakpoints, midpoint cells with the exact Poisson mass beyond distance 1000 out to
# |tau - t0| <= 13000) and the evaluator are those of bsy_cm.py and bsy_zdh.py, imported, not copied. mu is computed first from the same nodes; then A with the weights and the same
# tail estimate (m_out times the unweighted mean of |f - mu| over the outer half of the midpoint region; the tail of f and f^2 as in bsy_cm.py). Two meshes (steps 0.125 and 0.0625
# of the tanh-sinh rule, midpoint cells 0.5 and 0.25); the values printed are from the finer one. Error control: the error of R is the larger of |R_coarse - R_fine| and half of
# |R with the tail correction - R without it| (fine mesh); a second parameter set of the evaluator (nf 0.45, M 36, 37 levels) on 600 random nodes of every fine mesh, the largest
# differences of f logged; the mass defect |sum of weights + m_out - 1| of every mesh is logged.
# Output: data/bsy_isbound.json and the log data/bsy_isbound.log. Usage: python3 -I bsy_isbound.py [--quick] [--out FILE] [--procs P]; --quick (a few seconds, window 100 and 400)
# writes to /tmp unless --out is given; the full case on the cluster:
#   bash misc/tools/hopper_run.sh -c 16 -t 120 -m 32G -n bsy_isbound -g bsy_isbound.json -g bsy_isbound.log code/bsy_isbound.py code/bsy_cm.py code/bsy_zdh.py code/bsy_levy.py data/zeros_1700.txt -- --out bsy_isbound.json
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
import bsy_cm as C

ROWS = [("zeta", 0.5, 0.0), ("zeta", 0.875, 0.0), ("zeta", 0.5, 85.7), ("zeta", 0.875, 85.7), ("zeta", 0.75, 85.7),
        ("DH", 0.5, 85.7), ("DH", 0.75, 85.7), ("DH", 0.875, 85.7), ("DH", 0.75, 100.0)]


def workX(args):
    kind, theta, tau, kw = args
    tau = np.where(np.abs(tau) < 1e-6, 1e-6, tau)
    return np.log(np.abs(Z.hfun("zeta" if kind == "zeta" else "dh", theta, theta + 1j * tau, **kw)))


def evalX(pool, kind, theta, tau, kw=None, chunk=500):
    kw = kw or {}
    order = np.argsort(np.abs(tau))
    jobs = [(kind, theta, tau[order[i:i + chunk]], kw) for i in range(0, len(tau), chunk)]
    res = list(pool.imap(workX, jobs)) if pool else [workX(j) for j in jobs]
    X = np.empty(len(tau))
    X[order] = np.concatenate(res)
    return X


def stats(w, X, tau, t0, W2, W3, a, tail=True):
    m_out = 2 / math.pi * math.atan(a / W3)
    far = np.abs(tau - t0) > (W3 + W2) / 2
    r1, r2 = m_out * X[far].mean(), m_out * (X[far] ** 2).mean()
    if not tail:
        r1 = r2 = 0.0
    mu = (w * X).sum() + r1
    var = (w * X * X).sum() + r2 - mu * mu
    rA = m_out * np.abs(X[far] - mu).mean() if tail else 0.0
    A = (w * np.abs(X - mu)).sum() + rA
    return {"mu": mu, "var": var, "A": A, "R": var / A ** 2, "tail_mu": r1, "tail_A": rA}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--procs", type=int, default=0)
    A = ap.parse_args()
    if A.quick:
        W2, W3, step3 = 100.0, 400.0, 2.0
        hts = (0.25, 0.125)
        outp = A.out or "/tmp/bsy_isbound_quick.json"
    else:
        W2, W3, step3 = 1000.0, 13000.0, 0.5
        hts = C.HTS
        outp = A.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "bsy_isbound.json")
    lf = open(os.path.splitext(outp)[0] + ".log", "w")

    def log(m):
        print(m, flush=True)
        lf.write(m + "\n")
        lf.flush()

    t00 = time.time()
    log(f"bsy_isbound.py {'quick' if A.quick else 'full'}: R = Var_P f / (E_P|f - mu|)^2, f = log|h_theta(theta + i tau)|; panels on |tau - t0| <= {W2:g} (tanh-sinh steps {hts}), midpoint cells to {W3:g}")
    here = os.path.dirname(os.path.abspath(__file__))
    zp = next(q for q in (os.path.join(here, "zeros_1700.txt"), os.path.join(here, "..", "data", "zeros_1700.txt")) if os.path.exists(q))
    zlist = np.array([float(x) for x in open(zp).read().split()])
    import bsy_levy as LV
    Tz = 100.0 + W2 + 5.0
    dl, _ = LV.line_zeros(Tz, 0.005)
    do, _ = LV.offline_zeros(300.0)
    dhz = {"line": dl, "off": np.array([[z.real, z.imag] for z in do])}
    log(f"  zeta ordinates {len(zlist)} (largest {zlist[-1]:.3f}); DH: {len(dl)} sign changes of Lambda(1/2 + it) on (0, {Tz:g}], {len(do)} zeros off the line below 300")
    procs = A.procs or (1 if A.quick else (os.cpu_count() or 1))
    pool = None
    if procs > 1:
        import multiprocessing as mpc
        pool = mpc.Pool(procs)
    out = []
    for kind, theta, t0 in ROWS:
        t1 = time.time()
        a = 1 - theta
        log(f"row {kind} theta={theta} t0={t0}")
        brk = C.breaks_for(kind, theta, t0, W2, zlist, dhz)
        r = {"kind": kind, "theta": theta, "t0": t0}
        for name, h, sc in (("coarse", hts[0], 1.0), ("fine", hts[1], 0.5)):
            tau, w = C.nodes(t0, a, W2, W3, step3, brk, h, sc)
            X = evalX(pool, kind, theta, tau)
            m_out = 2 / math.pi * math.atan(a / W3)
            s = stats(w, X, tau, t0, W2, W3, a)
            s0 = stats(w, X, tau, t0, W2, W3, a, tail=False)
            s.update({"nodes": int(len(tau)), "mass_defect": float(abs(w.sum() + m_out - 1)), "m_out": m_out, "R_notail": s0["R"], "mu_notail": s0["mu"]})
            r[name] = s
            log(f"  {name}: nodes {len(tau)}, mass defect {s['mass_defect']:.1e}, mu {s['mu']:+.6f}, Var_P f {s['var']:.6f}, E|f-mu| {s['A']:.6f}, R {s['R']:.5f} (without tail {s0['R']:.5f})")
        r["R"] = r["fine"]["R"]
        r["error"] = max(abs(r["coarse"]["R"] - r["fine"]["R"]), 0.5 * abs(r["fine"]["R"] - r["fine"]["R_notail"]))
        rng = np.random.default_rng(1)
        tau, _ = C.nodes(t0, a, W2, W3, step3, brk, hts[1], 0.5)
        sub = rng.choice(len(tau), size=min(600, len(tau)), replace=False)
        X1 = evalX(None, kind, theta, tau[sub])
        X2 = evalX(None, kind, theta, tau[sub], kw={"nf": 0.45, "M": 36})
        ok = (X1 > -18.0) & (X2 > -18.0)
        r["check_params"] = {"max_dX": float(np.abs(X1 - X2)[ok].max()), "excluded_nodes": int((~ok).sum())}
        log(f"  R = {r['R']:.5f} +- {r['error']:.1e}; second parameters (nf 0.45, M 36) on {len(sub)} nodes ({int((~ok).sum())} with |h| < 1e-8 left out): max |dX| {r['check_params']['max_dX']:.2e}; row time {time.time() - t1:.0f} s")
        out.append(r)
    json.dump({"W2": W2, "W3": W3, "step3": step3, "quick": A.quick, "rows": out}, open(outp, "w"), indent=1)
    log(f"done in {time.time() - t00:.0f} s; wrote {outp}")


if __name__ == "__main__":
    main()
