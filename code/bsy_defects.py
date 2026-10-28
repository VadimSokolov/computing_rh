# bsy_defects.py: the quadrature column of Table bsy:tab:defects and the rows of the former Table tab:qrh of note L (Chapter 34): the local defect
# D_theta(t0) = int log|F(theta + i tau)| P_{1 - theta}(t0 - tau) d tau - log|h_theta(1 + i t0)| for zeta and for the Davenport and Heilbronn (DH) function, windows
# W = 60 and 200, with a rigorous-in-form bound of the neglected tail. Task C3 of reconcile-M.md (and row N24 of reconcile-L.md). Added for the book in October 2026,
# written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. P_a(x) = a/(pi (a^2 + x^2)), a = 1 - theta. F = zeta with h_theta(s) = zeta(s)(s - 1)/(s - 2 theta + 1), h_theta(1) = 1/(2 - 2 theta), or F = f, the DH
# function f(s) = 5^-s [zeta(s, 1/5) + kappa zeta(s, 2/5) - kappa zeta(s, 3/5) - zeta(s, 4/5)], kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1), with h = f (entire). The
# windowed value D_W = int_{|tau - t0| <= W} log|F| P d tau - log|h(1 + i t0)| (the neglected Cauchy mass is 1 - (2/pi) arctan(W/a) per window, logged).
# Quadrature. scipy.integrate.quad (QUADPACK, double precision, epsabs = epsrel = 1e-12, limit 200) on every piece between consecutive breakpoints; the breakpoints are the
# window ends, 0, and the ordinates +-gamma of the zeros on the critical line (zeta: the first 1700 of data/zeros_1700.txt; DH: all sign changes of the real function
# Lambda(1/2 + it) = (5/pi)^(s/2) Gamma((s + 1)/2) f(s), found on a grid of step 0.02 and bisected, plus the off line zeros, found as local minima of |Lambda| without a sign
# change and refined by Newton's method). The logged error is the sum of the QUADPACK error estimates; a second pass with epsabs = 1e-10 and limit 100 gives the
# second difference. The evaluation is double precision from python-flint acb (53 bits); for the row zeta (1/2, 85.7), W = 60 it is repeated at 106 bits.
# Tail. By Cauchy and Schwarz on blocks B_k = [k, k + 1] (the two partial blocks at the window ends are treated as blocks), |int_tail f P| <= sum_k (int_{B_k} f^2)^(1/2)
# (int_{B_k} P^2)^(1/2), the second factor in closed form. The first factors m_k^2 = int_{B_k} f^2 are computed (same quadrature) for every block up to height Kmax
# (zeta 1000, DH 500); beyond Kmax the bound uses max_k m_k times sum of (int_B P^2)^(1/2) over all further blocks (closed form, to 10^7); that last part is a
# heuristic extrapolation (the unit interval second moment of log|zeta| grows like (1/2) log log T) and is logged separately from the computed part.
# Predictions. zeta: 0 (QRH(theta)). DH: the exact sum of G_theta(1 + i t0, rho) = log|(1 + i t0 - 2 theta + conj rho)/(1 + i t0 - rho)| over the off line zeros found (all zeros
# with beta > theta, both signs of gamma) with the contribution of the zero near 0.8085 + 85.6993 i listed separately.
# Output: data/bsy_defects.json and the log bsy_defects.log. Usage: python3 bsy_defects.py [--quick] [--out FILE] [--procs P] [--zeros FILE]; --quick: three rows, W = 6 and 12,
# mpmath, about 25 seconds, to misc/reviews/nick-b4/quick/BSY-B/; full run on the cluster: bash misc/tools/hopper_run.sh -c 8 -t 120 -g bsy_defects.json -g bsy_defects.log
# code/bsy_defects.py data/zeros_1700.txt -- --out bsy_defects.json --zeros zeros_1700.txt
import argparse, json, math, os, time
from pathlib import Path
import numpy as np
from scipy.integrate import quad

ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None); ap.add_argument("--procs", type=int, default=0)
ap.add_argument("--zeros", default=None)
args = ap.parse_args()
here = Path(__file__).resolve().parent
out = Path(args.out) if args.out else (Path("/tmp/bsy_defects_quick.json") if args.quick else here.parent / "data" / "bsy_defects.json")
zfile = Path(args.zeros) if args.zeros else here.parent / "data" / "zeros_1700.txt"
try:
    from flint import acb, arb, ctx
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False
    import mpmath as mp
    mp.mp.dps = 20
PREC = [53]
KAPPA = ((math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1))

def zeta_c(s):
    if HAVE_FLINT:
        ctx.prec = PREC[0]; z = acb(s.real, s.imag).zeta(); return complex(float(z.real.mid()), float(z.imag.mid()))
    return complex(mp.zeta(mp.mpc(s.real, s.imag)))
def dh_c(s):
    if HAVE_FLINT:
        ctx.prec = PREC[0]
        ss = acb(s.real, s.imag)
        kap = (arb(10) - 2*arb(5).sqrt()).sqrt() - 2
        kap = kap/(arb(5).sqrt() - 1)
        v = ss.zeta(acb(arb(1)/5)) + kap*ss.zeta(acb(arb(2)/5)) - kap*ss.zeta(acb(arb(3)/5)) - ss.zeta(acb(arb(4)/5))
        v = v*(acb(5)**(-ss))
        return complex(float(v.real.mid()), float(v.imag.mid()))
    S = mp.mpc(s.real, s.imag); K = mp.mpf(KAPPA)
    v = (mp.zeta(S, mp.mpf(1)/5) + K*mp.zeta(S, mp.mpf(2)/5) - K*mp.zeta(S, mp.mpf(3)/5) - mp.zeta(S, mp.mpf(4)/5))*mp.mpf(5)**(-S)
    return complex(v)
def lam_real(t):   # Lambda(1/2 + it), real
    s = complex(0.5, t)
    if HAVE_FLINT:
        ctx.prec = PREC[0]
        ss = acb(0.5, t)
        v = dh_c(s)
        g = ((ss + 1)/2).gamma(); c = (acb(5)/acb.pi())**(ss/2)
        w = c*g
        return (complex(float(w.real.mid()), float(w.imag.mid()))*v).real
    return (complex((5/mp.pi)**(S/2)*mp.gamma((S + 1)/2)) * v).real if False else (complex(mp.mpf(5)**(mp.mpc(s)/2)/mp.pi**(mp.mpc(s)/2)*mp.gamma((mp.mpc(s) + 1)/2))*dh_c(s)).real
def logabs(kind, theta, tau):
    s = complex(theta, abs(tau))
    z = zeta_c(s) if kind == "zeta" else dh_c(s)
    return math.log(abs(z))
def logh(kind, theta, t0):
    if kind == "zeta":
        if t0 == 0: return math.log(1/(2 - 2*theta))
        s = complex(1, t0); return math.log(abs(zeta_c(s)*(s - 1)/(s - 2*theta + 1)))
    return math.log(abs(dh_c(complex(1, t0))))
def Pa(a, x): return a/(math.pi*(a*a + x*x))
def intP2(a, x0, x1):
    F = lambda x: a*a/math.pi**2*(x/(2*a*a*(a*a + x*x)) + math.atan(x/a)/(2*a**3))
    return F(x1) - F(x0)

def pieces(lo, hi, bps):
    pts = sorted(set([lo, hi] + [b for b in bps if lo < b < hi]))
    return list(zip(pts[:-1], pts[1:]))
def qint(fun, lo, hi, bps, tight=True):
    tot = 0.0; err = 0.0
    for a_, b_ in pieces(lo, hi, bps):
        v, e = quad(fun, a_, b_, epsabs=1e-12 if tight else 1e-10, epsrel=1e-12 if tight else 1e-10, limit=200 if tight else 100)
        tot += v; err += e
    return tot, err

ZETA_ZEROS = None; DH_ZEROS = []; DH_OFF = []
def breakpoints(kind):
    zs = ZETA_ZEROS if kind == "zeta" else DH_ZEROS
    return sorted(set([0.0] + list(zs) + [-z for z in zs]))

def m2_block(kind, theta, x0, x1):
    bps = breakpoints(kind)
    v, e = qint(lambda t: logabs(kind, theta, t)**2, x0, x1, bps, tight=False)
    return v

def run_row(job):
    kind, theta, t0, W, tight = job
    a = 1 - theta; bps = breakpoints(kind)
    t_start = time.time()
    I, err = qint(lambda t: logabs(kind, theta, t)*Pa(a, t0 - t), t0 - W, t0 + W, bps, tight)
    return {"kind": kind, "theta": theta, "t0": t0, "W": W, "tight": tight, "I": I, "err": err, "logh": logh(kind, theta, t0), "sec": time.time() - t_start}

def run_m2(job):
    kind, theta, K = job
    return (kind, theta, [m2_block(kind, theta, k, k + 1) for k in range(K)])

def find_dh_zeros(T, step):
    ts = np.arange(0.5*step, T, step)
    vals = np.array([lam_real(t) for t in ts])
    zs = []
    for i in range(len(ts) - 1):
        if vals[i]*vals[i + 1] < 0:
            a_, b_ = ts[i], ts[i + 1]; fa = vals[i]
            for _ in range(40):
                m_ = 0.5*(a_ + b_); fm = lam_real(m_)
                if fa*fm <= 0: b_ = m_
                else: a_, fa = m_, fm
            zs.append(0.5*(a_ + b_))
    # off line candidates: local minima of |Lambda| without a sign change nearby
    off = []
    av = np.abs(vals)
    for i in range(1, len(ts) - 1):
        if av[i] < av[i - 1] and av[i] < av[i + 1] and vals[i - 1]*vals[i] > 0 and vals[i]*vals[i + 1] > 0:
            off.append(ts[i])
    return zs, off

def newton_offline(t_guess):
    # Newton for f(s) = 0 near 0.5 + 0.3 + i t_guess (complex secant via mpmath-free finite differences)
    best = None
    for b0 in (0.8, 0.7, 0.9, 0.65):
        s = complex(b0, t_guess)
        for _ in range(60):
            f0 = dh_c(s); h = 1e-6
            d = (dh_c(s + h) - dh_c(s - h))/(2*h)
            if d == 0: break
            ds = f0/d; s = s - ds
            if abs(ds) < 1e-12: break
        if abs(dh_c(s)) < 1e-9 and abs(s.real - 0.5) > 1e-4 and abs(s.real - 0.5) < 1.5:
            best = s; break
    return best

def tail_bound(kind, theta, t0, W, a, m2v):
    """Cauchy-Schwarz bound of |int_{|tau - t0| > W} log|F| P d tau| over integer aligned unit blocks (computed up to height K) plus an extrapolation beyond K."""
    if m2v is None: return None
    K = len(m2v)
    F = lambda x: a*a/math.pi**2*(x/(2*a*a*(a*a + x*x)) + np.arctan(x/a)/(2*a**3))
    def finite(lo, hi):
        pts = [lo] + [float(k) for k in range(math.floor(lo) + 1, math.ceil(hi)) if lo < k < hi] + [hi]
        tot = 0.0
        for x0, x1 in zip(pts[:-1], pts[1:]):
            if x1 <= x0: continue
            h0, h1 = sorted((abs(x0), abs(x1)))
            full = abs(x1 - x0 - 1) < 1e-12 and abs(h0 - round(h0)) < 1e-12 or (abs(x1 - x0 - 1) < 1e-12 and abs(h1 - round(h1)) < 1e-12)
            mv = m2v[int(round(min(abs(x0), abs(x1))))] if (full and int(round(min(abs(x0), abs(x1)))) < K) else m2_block(kind, theta, h0, h1)
            tot += math.sqrt(max(mv, 0.0)*max(F(x1 - t0) - F(x0 - t0), 0.0))
        return tot
    eR, eL = t0 + W, t0 - W
    computed = finite(eR, float(K)) + finite(-float(K), eL)
    mmax = math.sqrt(max(m2v))
    ks = np.arange(K, 10**7, dtype=float)
    extra = 0.0
    for lo_arr, hi_arr in ((ks, ks + 1), (-ks - 1, -ks)):
        p2 = F(hi_arr - t0) - F(lo_arr - t0)
        extra += mmax*float(np.sum(np.sqrt(np.maximum(p2, 0))))
    return {"computed_blocks": computed, "extrapolated_beyond_Kmax": extra}

if __name__ == "__main__":
    log = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
    t_start = time.time()
    ZETA_ZEROS = [float(x) for x in zfile.read_text().split()[:1700]]
    quick = args.quick
    Tdh = 40 if quick else 700
    DH_ZEROS, off_t = find_dh_zeros(Tdh, 0.05 if quick else 0.02)
    P(f"DH: {len(DH_ZEROS)} sign changes of Lambda(1/2 + it) on (0, {Tdh}); {len(off_t)} candidate minima without sign change at", [round(x, 2) for x in off_t])
    for tg in off_t:
        z = newton_offline(tg)
        if z is not None and z.imag > 0 and all(abs(z - w) > 1e-6 for w in DH_OFF): DH_OFF.append(z)
    DH_OFF = [z for z in DH_OFF if z.real > 0.5]
    P("DH off line zeros with beta > 1/2 (Newton, |f| < 1e-9):", [f"{z.real:.6f}+{z.imag:.6f}i" for z in DH_OFF])
    if not quick:
        P("count check for 0 < gamma < 300: sign changes", sum(1 for z in DH_ZEROS if z < 300), "+ 2 x off line pairs", 2*sum(1 for z in DH_OFF if z.imag < 300), "=",
          sum(1 for z in DH_ZEROS if z < 300) + 2*sum(1 for z in DH_OFF if z.imag < 300), "; main term (T/2pi) log(5T/(2 pi e)) =", round(300/(2*math.pi)*math.log(5*300/(2*math.pi*math.e)), 2))
    # breakpoints: also the off line ordinates
    DH_ZEROS = sorted(DH_ZEROS + [z.imag for z in DH_OFF])
    Ws = [6, 12] if quick else [60, 200]
    rows = [("zeta", 0.5, 0.0), ("zeta", 0.5, 25.0), ("dh", 0.5, 5.0)] if quick else [("zeta", 0.5, 0.0), ("zeta", 0.875, 0.0), ("zeta", 0.5, 85.7), ("zeta", 0.875, 85.7), ("zeta", 0.75, 85.7),
            ("dh", 0.5, 85.7), ("dh", 0.75, 85.7), ("dh", 0.875, 85.7), ("dh", 0.75, 100.0)]
    jobs = [(k, th, t0, W, True) for (k, th, t0) in rows for W in Ws] + [(k, th, t0, Ws[0], False) for (k, th, t0) in rows]
    m2jobs = [] if quick else [(k, th, (1000 if k == "zeta" else 500)) for k, th in sorted(set((k, th) for k, th, _ in rows))]
    if quick: m2jobs = [("zeta", 0.5, 40)] if False else []
    procs = args.procs or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1)
    import multiprocessing as mpc
    if procs > 1 and not quick:
        with mpc.Pool(procs) as pool:
            r1 = pool.map_async(run_row, jobs, chunksize=1); r2 = pool.map_async(run_m2, m2jobs, chunksize=1)
            results = r1.get(); m2res = r2.get()
    else:
        results = [run_row(j) for j in jobs]; m2res = [run_m2(j) for j in m2jobs]
    m2 = {(k, th): v for k, th, v in m2res}
    # precision check of the evaluator
    if HAVE_FLINT:
        ts = np.random.default_rng(7).uniform(0, 150, 300)
        d = 0.0
        for kind in ("zeta", "dh"):
            PREC[0] = 53; a1 = np.array([logabs(kind, 0.75, t) for t in ts]); PREC[0] = 106; a2 = np.array([logabs(kind, 0.75, t) for t in ts]); PREC[0] = 53
            d = max(d, float(np.max(np.abs(a1 - a2))))
        P(f"evaluator check: max |log|F| at 53 bits - at 106 bits| over 300 random heights, both functions, theta = 3/4: {d:.2e}")
    # predictions for DH
    def dh_pred(theta, t0, skip85=False):
        tot = 0.0
        for z in DH_OFF:
            if z.real <= theta: continue
            if skip85 and abs(z.imag - 85.6993) < 0.01: continue
            for rho in (z, z.conjugate()):
                s = complex(1, t0)
                tot += math.log(abs((s - 2*theta + rho.conjugate())/(s - rho)))
        return tot
    summary = []
    byrow = {}
    for r in results: byrow.setdefault((r["kind"], r["theta"], r["t0"]), {})[(r["W"], r["tight"])] = r
    for (k, th, t0) in rows:
        a = 1 - th
        pred = 0.0 if k == "zeta" else dh_pred(th, t0)
        rec = {"function": k, "theta": th, "t0": t0, "prediction": pred}
        if k == "dh":
            rec["prediction_single_zero_0.8085"] = dh_pred(th, t0) - dh_pred(th, t0, True)
            rec["other_zeros_contribution"] = dh_pred(th, t0, True)
        for W in Ws:
            r = byrow[(k, th, t0)][(W, True)]
            D = r["I"] - r["logh"]
            mass = 1 - 2/math.pi*math.atan(W/a)
            # tail bound
            tb = tail_bound(k, th, t0, W, a, m2.get((k, th)))
            rec[f"W{W}"] = {"D_W": D, "quad_err": r["err"], "Cauchy_mass_outside": mass, "tail_bound": tb, "sec": r["sec"]}
        r2 = byrow[(k, th, t0)][(Ws[0], False)]
        rec["second_pass_difference_W%d" % Ws[0]] = (r2["I"] - r2["logh"]) - rec[f"W{Ws[0]}"]["D_W"]
        summary.append(rec)
        P(f"{k:4s} theta={th} t0={t0}: prediction {pred:.6g}" + (f" (single zero {rec['prediction_single_zero_0.8085']:.6g}, others {rec['other_zeros_contribution']:.3g})" if k == "dh" else ""))
        for W in Ws:
            x = rec[f"W{W}"]
            tbs = "" if x["tail_bound"] is None else f" tail bound (computed blocks) {x['tail_bound']['computed_blocks']:.3e} + (extrapolated) {x['tail_bound']['extrapolated_beyond_Kmax']:.3e}"
            P(f"    W={W}: D_W = {x['D_W']:.6e}  (quad error estimate {x['quad_err']:.1e}; Cauchy mass outside {x['Cauchy_mass_outside']:.4f}){tbs}")
        P(f"    second pass (epsabs 1e-10) minus first at W={Ws[0]}: {rec['second_pass_difference_W%d' % Ws[0]]:.2e}")
    P("total time", round(time.time() - t_start, 1), "s;", "flint" if HAVE_FLINT else "mpmath")
    Path(out).write_text(json.dumps({"rows": summary, "dh_off": [str(z) for z in DH_OFF], "dh_online_count_below_700": len(DH_ZEROS)}, indent=1))
    Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
