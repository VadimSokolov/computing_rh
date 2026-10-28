# bsy_riesz2d.py: Table bsy:tab:2d and Figures bsy:fig:maps and bsy:fig:rates2d of Chapter 35 (ch/bsy2.tex, Section bsy:sec:2d, "Two dimensions: the smoothed Riesz identity"), group BSY-C.
# Claims supported: the smoothed Riesz identity I(phi) = iint log|f| Delta phi dsigma dt = 2 pi (sum phi(rho)) of Theorem bsy:thm:riesz2d is estimated by plain Monte Carlo, Voronoi
# weights, Delaunay weights (uniform and jittered nodes), the randomly shifted product trapezoid (grid) rule and the deflated grid and Delaunay rules; root mean square errors over
# replications for zeta (exact value 0, no zero in the box) and for the Davenport and Heilbronn function (DH; exact value 2 pi phi(rho), rho the zero near 0.8085 + 85.6993 i); fitted
# slopes; the maps of the integrand (Figure bsy:fig:maps). Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. Box 0.6 <= sigma <= 1.4, 83.7 <= t <= 87.7, centre (1, 85.7), half widths hs = 0.4 and ht = 2; x = (sigma - 1)/hs, y = (t - 85.7)/ht in [-1, 1]^2. phi = b(x) b(y), b(x) =
# exp(-1/(1 - x^2)) for |x| < 1 (else 0), so Delta phi = b''(x) b(y)/hs^2 + b(x) b''(y)/ht^2 in the (sigma, t) plane, b'' = b (4x^2/u^4 - 2/u^2 - 8x^2/u^3), u = 1 - x^2; the integral is
# hs ht times the integral over [-1, 1]^2. The zero rho is the certified zero z85 of data/cf_dhzeros.json (T2). Deflation: log|f| is replaced by log|f| - log|s - rho_hat| and
# 2 pi phi(rho_hat) is added back (exact by Theorem bsy:thm:riesz2d); rho_hat = rho + eps (1 + i)/sqrt 2.
# Estimators (N nodes in [-1, 1]^2, all weights sum to the box area): plain = area * mean; Voronoi = weights are the areas of the Voronoi cells of the uniform nodes clipped to the
# box (cells of the nodes and of the mirror images of the nodes in a strip along each side, scipy.spatial.Voronoi; the sum of the areas is checked to equal 4); Delaunay = the
# integral of the piecewise linear interpolant on the Delaunay triangulation of the nodes and the four corners (value 0 at the corners, where Delta phi = 0), for uniform nodes and for
# jittered nodes (one uniform point in each cell of an m x m grid, m = sqrt N); grid = (area/N) sum of G over x_j = -1 + 2 (j + u)/m, u uniform in [0,1)^2 (random shift), G = the integrand
# supported in the box; deflated grid and deflated jittered Delaunay use the deflated integrand with the same nodes.
# Error control. The Euler and Maclaurin evaluator (code/bsy_zdh.py, imported) is run at (nf, M) = (0.3, 28) and again at (0.45, 36) on the nodes of one replicate: the maximum
# difference of log|f| and the difference of the grid estimate are logged; Delta phi is compared with a central finite difference at random points; the Voronoi areas sum to the box area.
# The seeds are 20261025 + replicate.
# Output: data/bsy_riesz2d.json, the log data/bsy_riesz2d.log. Usage: python3 bsy_riesz2d.py [--quick] [--out FILE] [--procs P] [--reps R]; --quick: N = 1024, 4096, 3 replicates, a few
# seconds, to misc/reviews/nick-b4/quick/BSY-C/; full run (N = 4^5 .. 4^9, 20 replicates) on the cluster: bash misc/tools/hopper_run.sh -c 8 -t 90 -g bsy_riesz2d.json -g bsy_riesz2d.log
# code/bsy_riesz2d.py -- --out bsy_riesz2d.json
import argparse, json, math, os, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial import Delaunay, Voronoi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_zdh as Z

HS, HT, SC, TC = 0.4, 2.0, 1.0, 85.7
HERE = Path(__file__).resolve().parent


def load_rho():
    for c in (HERE.parent / "data", HERE / "data", HERE):
        if (c / "cf_dhzeros.json").exists():
            d = json.load(open(c / "cf_dhzeros.json")); break
    z = [z for z in d["T2"]["zeros"] if z["tag"] == "z85"][0]
    return complex(float(z["center_re"]), float(z["center_im"]))


RHO = load_rho()


def bfun(x):
    out = np.zeros_like(x)
    m = np.abs(x) < 1
    out[m] = np.exp(-1.0 / (1.0 - x[m] ** 2))
    return out


def b2fun(x):
    out = np.zeros_like(x)
    m = np.abs(x) < 1
    xm = x[m]; u = 1.0 - xm ** 2
    out[m] = np.exp(-1.0 / u) * (4 * xm ** 2 / u ** 4 - 2 / u ** 2 - 8 * xm ** 2 / u ** 3)
    return out


def phi_st(sig, t):
    return bfun((np.asarray(sig) - SC) / HS) * bfun((np.asarray(t) - TC) / HT)


def lap_phi(x, y):
    return b2fun(x) * bfun(y) / HS ** 2 + bfun(x) * b2fun(y) / HT ** 2


def logf(kind, s, nf=0.3, M=28):
    f = Z.zeta(s, nf=nf, M=M) if kind == "zeta" else Z.dh(s, nf=nf, M=M)
    return np.log(np.abs(f))


def integrand(kind, x, y, rho_hat=None, nf=0.3, M=28):
    """G = log|f| Delta phi (minus log|s - rho_hat| Delta phi when deflating) at the points (x, y) of [-1,1]^2 (flat arrays)."""
    L = lap_phi(x, y)
    G = np.zeros_like(x)
    m = L != 0
    s = (SC + HS * x[m]) + 1j * (TC + HT * y[m])
    lg = logf(kind, s, nf, M)
    if rho_hat is not None:
        lg = lg - np.log(np.abs(s - rho_hat))
    G[m] = lg * L[m]
    return G


def exact_value(kind):
    return 0.0 if kind == "zeta" else 2 * math.pi * float(phi_st(RHO.real, RHO.imag))


def add_back(rho_hat):
    return 2 * math.pi * float(phi_st(rho_hat.real, rho_hat.imag))


def voronoi_areas(P, m):
    """areas of the Voronoi cells of the points P in [-1,1]^2 clipped to the box, by mirror images in a strip of width w along each side."""
    w = min(1.0, 24.0 / m)
    imgs = [P]
    for ax in (0, 1):
        for side in (-1.0, 1.0):
            sel = P[:, ax] * side > 1 - w
            Q = P[sel].copy(); Q[:, ax] = 2 * side - Q[:, ax]
            imgs.append(Q)
    A = np.vstack(imgs)
    vor = Voronoi(A)
    n = len(P)
    regs = [vor.regions[vor.point_region[i]] for i in range(n)]
    lens = np.fromiter((len(r) for r in regs), dtype=np.int64, count=n)
    flat = np.concatenate([np.asarray(r, dtype=np.int64) for r in regs])
    assert (flat >= 0).all(), "unbounded Voronoi cell"
    off = np.concatenate([[0], np.cumsum(lens)[:-1]])
    nxt = np.arange(len(flat)) + 1
    nxt[off + lens - 1] = off
    X = vor.vertices[flat]; Y = vor.vertices[flat[nxt]]
    cross = X[:, 0] * Y[:, 1] - X[:, 1] * Y[:, 0]
    return np.abs(0.5 * np.add.reduceat(cross, off))


def delaunay_integral(P, G):
    corners = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]], dtype=float)
    A = np.vstack([P, corners]); GA = np.concatenate([G, np.zeros(4)])
    T = Delaunay(A).simplices
    a, b, c = A[T[:, 0]], A[T[:, 1]], A[T[:, 2]]
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))
    return float((area * (GA[T[:, 0]] + GA[T[:, 1]] + GA[T[:, 2]]) / 3.0).sum())


def nodes(rng, N):
    m = int(round(math.sqrt(N)))
    U = rng.uniform(-1, 1, (N, 2))
    j = np.stack(np.meshgrid(np.arange(m), np.arange(m), indexing="ij"), -1).reshape(-1, 2)
    J = -1 + 2.0 * (j + rng.random((N, 2))) / m
    sh = rng.random(2)
    Gd = -1 + 2.0 * (j + sh) / m
    return m, U, J, Gd


def replicate(args):
    rep, Ns, kinds, eps_defl, check = args
    out = {}
    for N in Ns:
        rng = np.random.default_rng(20261025 + 1000 * rep + N)
        m, U, J, Gd = nodes(rng, N)
        area = 4.0 * HS * HT
        wv = voronoi_areas(U, m)
        for kind in kinds:
            ex = exact_value(kind)
            r = {}
            GU = integrand(kind, U[:, 0], U[:, 1]); GJ = integrand(kind, J[:, 0], J[:, 1]); GG = integrand(kind, Gd[:, 0], Gd[:, 1])
            r["plain"] = area * GU.mean()
            r["voronoi"] = HS * HT * float((wv * GU).sum())
            r["delaunay"] = HS * HT * delaunay_integral(U, GU)
            r["jittered"] = HS * HT * delaunay_integral(J, GJ)
            r["grid"] = area * GG.mean()
            if kind == "dh":
                for eps, tag in eps_defl:
                    rh = RHO + eps * (1 + 1j) / math.sqrt(2)
                    GGd = integrand(kind, Gd[:, 0], Gd[:, 1], rh)
                    r["grid_defl_" + tag] = area * GGd.mean() + add_back(rh)
                rh0 = RHO
                GJd = integrand(kind, J[:, 0], J[:, 1], rh0)
                r["jittered_defl_0"] = HS * HT * delaunay_integral(J, GJd) + add_back(rh0)
            out[f"{kind}|{N}"] = {k: v - ex for k, v in r.items()}
            if check and rep == 0:
                GG2 = integrand(kind, Gd[:, 0], Gd[:, 1], None, 0.45, 36)
                s = (SC + HS * Gd[:, 0]) + 1j * (TC + HT * Gd[:, 1])
                dl = np.abs(logf(kind, s) - logf(kind, s, 0.45, 36)).max()
                out[f"check|{kind}|{N}"] = {"max_diff_logf": float(dl), "grid_diff": float(area * (GG2.mean() - GG.mean())), "voronoi_area_sum": float(wv.sum())}
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
    ap.add_argument("--procs", type=int, default=0); ap.add_argument("--reps", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out) if args.out else (Path("/tmp/bsy_riesz2d_quick.json") if args.quick else HERE.parent / "data" / "bsy_riesz2d.json")
    log = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
    Ns = [1024, 4096] if args.quick else [1024, 4096, 16384, 65536, 262144]
    R = args.reps or (3 if args.quick else 20)
    eps_defl = [(1.4e-4, "1.4e-4"), (0.0, "0")]
    P(f"bsy_riesz2d: N = {Ns}, replicates {R}, rho = {RHO!r}, exact DH value 2 pi phi(rho) = {exact_value('dh'):.12f}, box {SC-HS}..{SC+HS} x {TC-HT}..{TC+HT}")
    # Delta phi against a central finite difference
    rg = np.random.default_rng(1)
    xs, ys = rg.uniform(-0.9, 0.9, 5), rg.uniform(-0.9, 0.9, 5)
    d = 1e-4
    ph = lambda x, y: bfun(x) * bfun(y)
    fd = (ph(xs + d, ys) + ph(xs - d, ys)) / HS ** 2 * 0 + ((ph(xs + d, ys) - 2 * ph(xs, ys) + ph(xs - d, ys)) / (d * HS) ** 2 + (ph(xs, ys + d) - 2 * ph(xs, ys) + ph(xs, ys - d)) / (d * HT) ** 2)
    P("Delta phi vs finite difference, max relative difference:", float(np.max(np.abs(fd - lap_phi(xs, ys)) / np.abs(lap_phi(xs, ys)))))
    t0 = time.time()
    jobs = [(r, Ns, ["zeta", "dh"], eps_defl, True) for r in range(R)]
    procs = args.procs or (1 if args.quick else min(R, os.cpu_count() or 1))
    if procs > 1:
        import multiprocessing as mpc
        with mpc.Pool(procs) as pool:
            res = pool.map(replicate, jobs, chunksize=1)
    else:
        res = [replicate(j) for j in jobs]
    P(f"replicates done in {time.time() - t0:.0f} s on {procs} processes")
    errors = {}; checks = {}
    for r in res:
        for k, v in r.items():
            if k.startswith("check|"):
                checks[k] = v
            else:
                for est, e in v.items():
                    errors.setdefault(f"{k}|{est}", []).append(e)
    summary = {}
    for k, v in errors.items():
        v = np.array(v)
        summary[k] = {"rms": float(np.sqrt((v ** 2).mean())), "bias": float(v.mean()), "reps": len(v)}
    # slopes for N >= 4096
    slopes = {}
    ests = sorted({k.split("|")[2] for k in summary})
    for kind in ("zeta", "dh"):
        for est in ests:
            pts = [(N, summary[f"{kind}|{N}|{est}"]["rms"]) for N in Ns if N >= 4096 and f"{kind}|{N}|{est}" in summary]
            if len(pts) >= 2:
                slopes[f"{kind}|{est}"] = float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0])
    for kind in ("zeta", "dh"):
        for N in Ns:
            P(kind, N, "  ".join(f"{e}={summary[f'{kind}|{N}|{e}']['rms']:.2e}" for e in ests if f"{kind}|{N}|{e}" in summary))
    for k, v in slopes.items():
        P("slope (N >= 4096)", k, f"{v:.2f}")
    for k, v in checks.items():
        P(k, v)
    # maps
    nm = 41 if args.quick else 121
    gx = np.linspace(-1, 1, nm); X, Y = np.meshgrid(gx, gx, indexing="ij")
    maps = {}
    for kind in ("zeta", "dh"):
        G = integrand(kind, X.ravel(), Y.ravel()).reshape(nm, nm)
        maps[kind] = [[float(f"{v:.5g}") for v in row] for row in G]
        P(kind, "map: integral by trapezoid on the map grid (a sanity check only)", float(HS * HT * np.trapezoid(np.trapezoid(G, gx, axis=1), gx)), "exact", exact_value(kind))
    out.write_text(json.dumps({"Ns": Ns, "reps": R, "rho": [RHO.real, RHO.imag], "exact": {"zeta": 0.0, "dh": exact_value("dh")}, "box": {"sigma": [SC - HS, SC + HS], "t": [TC - HT, TC + HT]},
                               "eps_defl": [e for e, _ in eps_defl], "summary": summary, "slopes": slopes, "checks": checks, "map_x": gx.tolist(), "maps": maps, "errors": errors}))
    Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
