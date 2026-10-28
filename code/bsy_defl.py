# bsy_defl.py: the two regimes of approximate deflation, Chapter 35 Proposition bsy:prop:defl and Table bsy:tab:defl (task C6 of reconcile-M.md).
# Claims supported: with a zero estimate at distance eps the shifted grid rule has root mean square error of order eps N^(-1/2) (up to a logarithm) while
# eps sqrt(N) <~ 1, and of order N^(-1) once eps sqrt(N) >~ 1; no floor; the crossover is near N = eps^-2. Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. Torus [0,1)^2, rho = (1/2, 1/2), rho_hat = rho + eps (1 + i)/sqrt 2, residual r(s) = psi(s) (log|s - rho| - log|s - rho_hat|), psi(s) =
# exp(1 - 1/(1 - (|s - rho|/R)^2)) for |s - rho| < R = 0.4, else 0 (a C^infinity bump, radial about rho). Rule: Q(delta) = (m1 m2)^-1 sum_j r(j/m + delta), N = m1 m2 =
# 2^k, m1 = 2^ceil(k/2), m2 = 2^floor(k/2), delta uniform on [0, 1/m1) x [0, 1/m2) (for odd k a rectangular grid).
# Exact value. The circle average of log|s - rho_hat| over |s - rho| = r is log max(r, eps), so int r = 2 pi int_0^eps psi(r) r (log r - log eps) dr, evaluated by
# mp.quad (30 digits). This is the reference: the error of a replicate is Q(delta) - int r.
# Statistics. R(N) independent uniform shifts (R = 400 for N <= 2^14, 200 for 2^16, 100 for 2^18, 40 for 2^20, 20 for 2^22), seed 12345 + k; the logged quantities
# are the root mean square error (rms), the mean error over the shifts (bias; the rule is unbiased), the relative standard error 1/sqrt(2 R) of the rms, and for each
# eps the local log-log slope on the two regimes (N <= eps^-2/16 and N >= 16 eps^-2), the ratios rms/(eps N^-1/2) and rms N, and the ratio of rms to the bound of
# the proposition, B = eps N^-1/2 (1 + log^+(1/(eps sqrt N)))^(1/2) + eps log(2/eps) N^-1 (eps sqrt N <= 1), B = N^-1 otherwise (its max and min over N is the constant C).
# Error control: the reference integral is recomputed at 50 digits (difference logged), and every rms is recomputed on a second independent seed block for
# N <= 2^16 (the two rms agree within the logged statistical error). The DH extension of Table bsy:tab:defl (optional) is not run here.
# Output: data/bsy_defl.json and the log bsy_defl.log. Usage: python3 bsy_defl.py [--quick] [--out FILE]; --quick: eps = 1e-2, 1e-3, N = 2^10..2^14, few shifts, a few
# seconds; the full run is about 2 minutes on one core: bash misc/tools/hopper_run.sh -c 1 -t 60 -g bsy_defl.json -g bsy_defl.log code/bsy_defl.py -- --out bsy_defl.json
import argparse, json, math
from pathlib import Path
import numpy as np
import mpmath as mp

ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
here = Path(__file__).resolve().parent
out = Path(args.out) if args.out else (Path("/tmp/bsy_defl_quick.json") if args.quick else here.parent / "data" / "bsy_defl.json")
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)

R0 = 0.4
def psi_np(d2):
    r2 = d2/R0**2
    out = np.zeros_like(d2)
    m = r2 < 1
    out[m] = np.exp(1.0 - 1.0/(1.0 - r2[m]))
    return out

def exact_int(eps, dps):
    mp.mp.dps = dps
    e = mp.mpf(eps)
    psi = lambda r: mp.exp(1 - 1/(1 - (r/R0)**2))
    return 2*mp.pi*mp.quad(lambda r: psi(r)*r*(mp.log(r) - mp.log(e)), [0, e])

def rule(eps, k, rng):
    m1 = 2**((k + 1)//2); m2 = 2**(k//2)
    d1 = rng.random()/m1; d2_ = rng.random()/m2
    xs = np.arange(m1)/m1 + d1 - 0.5          # relative to rho
    ys = np.arange(m2)/m2 + d2_ - 0.5
    sh = eps/math.sqrt(2.0)
    tot = 0.0
    step = max(1, (1 << 20)//m2)
    for i in range(0, m1, step):
        X = xs[i:i+step, None]; Y = ys[None, :]
        da = X*X + Y*Y
        db = (X - sh)**2 + (Y - sh)**2
        ps = psi_np(da)
        mask = ps > 0
        val = np.zeros_like(da)
        val[mask] = 0.5*ps[mask]*(np.log(da[mask]) - np.log(db[mask]))
        tot += val.sum()
    return tot/(m1*m2)

epss = [1e-2, 1e-3] if args.quick else [1e-2, 1e-3, 1e-4]
ks = [10, 11, 12, 13, 14] if args.quick else list(range(10, 23))
def nrep(k):
    if args.quick: return 8
    return 400 if k <= 14 else 200 if k <= 16 else 100 if k <= 18 else 40 if k <= 20 else 20
res = {"eps": {}}
for eps in epss:
    I1 = exact_int(eps, 30); I2 = exact_int(eps, 50)
    P(f"eps {eps}: exact integral {mp.nstr(I1, 15)}  (dps 30 minus dps 50: {mp.nstr(I1 - I2, 3)})")
    Iex = float(I1)
    rows = []
    for k in ks:
        N = 2**k; Rn = nrep(k)
        errs = []
        for blk in (0, 1) if (k <= 16 and not args.quick) else (0,):
            rng = np.random.default_rng(12345 + k + 1000*blk)
            e = np.array([rule(eps, k, rng) - Iex for _ in range(Rn)])
            errs.append(e)
        allv = np.concatenate(errs)
        rms = math.sqrt(np.mean(allv**2)); bias = float(np.mean(allv))
        rms_blocks = [math.sqrt(np.mean(e**2)) for e in errs]
        x = eps*math.sqrt(N)
        B = eps*N**-0.5*(1 + max(0.0, math.log(1/x)))**0.5 + eps*math.log(2/eps)/N if x <= 1 else 1.0/N
        rows.append({"k": k, "N": N, "R": len(allv), "rms": rms, "bias": bias, "rms_blocks": rms_blocks, "eps_sqrtN": x, "ratio_to_bound": rms/B,
                     "rms/(eps N^-1/2)": rms/(eps*N**-0.5), "rms*N": rms*N})
        P(f"  eps {eps:g} N=2^{k} R={len(allv)} rms={rms:.4e} bias={bias:+.2e} (se of rms {100/math.sqrt(2*len(allv)):.0f}%) eps*sqrtN={x:.3g} rms/(eps N^-1/2)={rms/(eps*N**-0.5):.3f} rms*N={rms*N:.3f} rms/B={rms/B:.3f}")
    # slopes on the two regimes
    def slope(sel):
        if len(sel) < 2: return None
        a = np.polyfit([math.log(r["N"]) for r in sel], [math.log(r["rms"]) for r in sel], 1)[0]
        return float(a)
    r1 = [r for r in rows if r["N"] <= eps**-2/16]; r2 = [r for r in rows if r["N"] >= 16*eps**-2]
    s1 = slope(r1); s2 = slope(r2)
    P(f"  eps {eps:g}: crossover N = eps^-2 = {eps**-2:.3g} (2^{math.log2(eps**-2):.1f}); slope for N <= eps^-2/16: {s1}; slope for N >= 16 eps^-2: {s2}")
    rb = [r["ratio_to_bound"] for r in rows]
    P(f"  eps {eps:g}: rms/B between {min(rb):.3f} and {max(rb):.3f}")
    res["eps"][str(eps)] = {"exact": mp.nstr(I1, 15), "rows": rows, "slope_regime1": s1, "slope_regime2": s2, "rms_over_bound_min": min(rb), "rms_over_bound_max": max(rb)}
    for r in rows:
        if r["k"] == 18: P(f"  eps {eps:g}: rms/eps at N=262144: {r['rms']/eps:.4f}")
Path(out).write_text(json.dumps(res, indent=1)); Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
