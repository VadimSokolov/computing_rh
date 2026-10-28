# bsy_defl_dh.py: Table bsy:tab:defl of Chapter 35 (ch/bsy2.tex, end of Section bsy:sec:2d), group BSY-C: the deflated shifted grid for the Davenport and Heilbronn (DH) function in the
# box of Table bsy:tab:2d, root mean square error of I(phi) over 10 random shifts as a function of the error eps of the deflating zero and of N = 4^k nodes, k = 5 .. 10 (N up to 2^20).
# Claims supported: Proposition bsy:prop:defl (two regimes: error of order eps N^(-1/2) up to a logarithm while the grid does not separate the two logarithmic singularities, of order N^(-1)
# afterwards; no floor), the row "none" (grid without deflation, rate N^(-1) log N of Proposition bsy:prop:rates2d(ii)), the constant of the text ("about 0.02 eps at N = 262144"), and
# the comparison with the synthetic test code/bsy_defl.py (data/bsy_defl.json) of reconcile-M D7. Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions and evaluator: those of code/bsy_riesz2d.py (imported: box, bump phi, Delta phi, the zero rho of data/cf_dhzeros.json, the Euler and Maclaurin evaluator of code/bsy_zdh.py).
# Method. For each shift (a pair u uniform in [0,1)^2, seed 20261026 + shift) and each N = m^2 the grid x_j = -1 + 2 (j + u)/m of [-1,1]^2 is evaluated once, log|f| at the nodes
# with Delta phi != 0; the estimate for each eps is (4 hs ht/N) sum_j (log|f| - log|s_j - rho_hat|) Delta phi + 2 pi phi(rho_hat), rho_hat = rho + eps (1 + i)/sqrt 2 (the shift of the
# estimate is in the (sigma, t) plane), and for "none" the same without deflation. The error is the estimate minus 2 pi phi(rho); the logged rms is the root mean square over the shifts, with
# the bias, the standard error of the rms (1/sqrt(2 R)), the local log-log slopes of the rms in N (consecutive pairs), rms/eps, rms/(eps N^(-1/2)) and rms N.
# Error control. The evaluator is run at (nf, M) = (0.3, 28) and (0.45, 36) on the nodes of shift 0 at N = 16384 and the difference of the estimates and of log|f| is logged; the
# exact value 2 pi phi(rho) uses the certified zero (15 digits certified, 22 digits stored): replacing rho by rho + 1e-15 changes it by less than 1e-15 (logged).
# Output: data/bsy_defl_dh.json and the log data/bsy_defl_dh.log. Usage: python3 bsy_defl_dh.py [--quick] [--out FILE] [--procs P]; --quick: N = 1024, 4096, 3 shifts, seconds, to
# misc/reviews/nick-b4/quick/BSY-C/; full run: bash misc/tools/hopper_run.sh -c 10 -t 120 -g bsy_defl_dh.json -g bsy_defl_dh.log code/bsy_defl_dh.py code/bsy_riesz2d.py code/bsy_zdh.py
# data/cf_dhzeros.json -- --out bsy_defl_dh.json
import argparse, json, math, os, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_riesz2d as R2

HERE = Path(__file__).resolve().parent
EPS = [0.0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2]


def one_shift(args):
    shift, Ns, check = args
    rng = np.random.default_rng(20261026 + shift)
    u = rng.random(2)
    exact = R2.exact_value("dh")
    area = 4.0 * R2.HS * R2.HT
    out = {}
    for N in Ns:
        m = int(round(math.sqrt(N)))
        g = -1 + 2.0 * (np.arange(m) + 0) / m
        X, Y = np.meshgrid(g + 2.0 * u[0] / m, g + 2.0 * u[1] / m, indexing="ij")
        x, y = X.ravel(), Y.ravel()
        L = R2.lap_phi(x, y)
        msk = L != 0
        s = (R2.SC + R2.HS * x[msk]) + 1j * (R2.TC + R2.HT * y[msk])
        lg = R2.logf("dh", s)
        Lm = L[msk]
        e = {"none": area * float((lg * Lm).sum()) / N - exact}
        for eps in EPS:
            rh = R2.RHO + eps * (1 + 1j) / math.sqrt(2)
            e[repr(eps)] = area * float(((lg - np.log(np.abs(s - rh))) * Lm).sum()) / N + R2.add_back(rh) - exact
        out[N] = e
        if check and shift == 0 and N == 16384:
            lg2 = R2.logf("dh", s, 0.45, 36)
            out["check"] = {"max_diff_logf": float(np.abs(lg - lg2).max()), "diff_est_eps0": area * float(((lg2 - lg) * Lm).sum()) / N}
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--procs", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out) if args.out else (Path("/tmp/bsy_defl_dh_quick.json") if args.quick else HERE.parent / "data" / "bsy_defl_dh.json")
    log = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
    Ns = [1024, 4096] if args.quick else [4 ** k for k in range(5, 11)]
    nsh = 3 if args.quick else 10
    P(f"bsy_defl_dh: N = {Ns}, shifts {nsh}, eps = {EPS}, rho = {R2.RHO!r}, exact 2 pi phi(rho) = {R2.exact_value('dh'):.15f}")
    dphi = 2 * math.pi * float(R2.phi_st(R2.RHO.real + 1e-15, R2.RHO.imag)) - R2.exact_value("dh")
    P("exact value change under rho + 1e-15:", dphi)
    t0 = time.time()
    jobs = [(k, Ns, True) for k in range(nsh)]
    procs = args.procs or (1 if args.quick else min(nsh, os.cpu_count() or 1))
    if procs > 1:
        import multiprocessing as mpc
        with mpc.Pool(procs) as pool:
            res = pool.map(one_shift, jobs, chunksize=1)
    else:
        res = [one_shift(j) for j in jobs]
    P(f"done in {time.time() - t0:.0f} s on {procs} processes")
    rows = ["none"] + [repr(e) for e in EPS]
    summ = {}
    for r in rows:
        summ[r] = {}
        for N in Ns:
            v = np.array([res[k][N][r] for k in range(nsh)])
            summ[r][N] = {"rms": float(np.sqrt((v ** 2).mean())), "bias": float(v.mean()), "errors": v.tolist()}
    for r in rows:
        P("eps" if r != "none" else "", r, "  ".join(f"N={N}: {summ[r][N]['rms']:.2e}" for N in Ns))
    slopes = {}
    for r in rows:
        sl = [math.log(summ[r][Ns[i + 1]]["rms"] / summ[r][Ns[i]]["rms"]) / math.log(Ns[i + 1] / Ns[i]) for i in range(len(Ns) - 1)]
        slopes[r] = sl
        P("local slopes", r, " ".join(f"{x:.2f}" for x in sl))
    for r in rows[1:]:
        e = float(r)
        if e > 0:
            P("rms/eps", r, " ".join(f"{summ[r][N]['rms'] / e:.3g}" for N in Ns), "| eps sqrt(N) =", " ".join(f"{e * math.sqrt(N):.3g}" for N in Ns))
    P("relative standard error of each rms (1/sqrt(2R)):", 1 / math.sqrt(2 * nsh))
    chk = res[0].get("check")
    P("check (nf, M) = (0.3, 28) vs (0.45, 36), shift 0, N = 16384:", chk)
    try:
        d = json.load(open(HERE.parent / "data" / "bsy_defl.json"))
        P("synthetic test (code/bsy_defl.py) keys:", list(d.keys())[:8])
    except Exception as ex:
        P("synthetic test file not read:", ex)
    out.write_text(json.dumps({"Ns": Ns, "shifts": nsh, "eps": EPS, "rho": [R2.RHO.real, R2.RHO.imag], "exact": R2.exact_value("dh"), "summary": summ, "slopes": slopes, "check": chk}))
    Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
