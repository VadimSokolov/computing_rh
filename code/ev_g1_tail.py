# ev_g1_tail.py: numerical claims of Nick Polson's extreme value hitting times chapter (E; Chapter 25 after the merge) about the hitting time T of the
# reciprocal Xi string: the Talbot estimate of the principal eigenvalue (E:155, prop:tail), the sup distances of the shifted exponential B splines to the density
# of X = sum_k E_k/gamma_k^2 (E:186, prop:xibs) and the zeros of the Dyck chain Dirichlet polynomials F_n(s) = E[Y_n^s] (E:80 to E:82, rows 5 and 6 of
# tab:approx at E:96 and E:97). Added for the book in October 2026 by computation group G1, written from the mathematics (the collaborator's scripts were not
# read, run or copied).
# Conventions. L(s) = E exp(-sT) = xi(1/2)/xi(1/2 + sqrt s), xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, T = sum_k E_k/gamma_k^2. Density f of T and tail
# S(x) = P(T > x) by inverse Laplace transform of L(s) and (1 - L(s))/s with the fixed Talbot rule of Abate and Valko (delta_0 = 2M/(5t), delta_k = (2 k pi/5)
# (cot(k pi/M) + i), f(t) = (2/(5t)) sum_{k<M} Re(gamma_k F(delta_k/t))), in ball arithmetic (python-flint) at 3 bits per decimal digit of margin. The rate estimate
# is the square root of the hazard f/S at x (the chapter does not define its estimator: this is the natural one, since S ~ c_1 exp(-gamma_1^2 x) gives f/S -> gamma_1^2).
# Deflation: S_1 = S - c_1 exp(-lam_1 x) with lam_1 and c_1 = S(x_a) exp(lam_1 x_a) estimated from the window, then the hazard of S_1 for the second atom, and so on (the
# chapter's windows for the deflation are not given: the script scans x and reports the x at which the chapter's printed values occur, which is not a reproduction).
# B splines: K-atom spline with knots lam_k = gamma_k^2, shifted by r_K = mu_1 - sum_{k<=K} lam_k^{-1}, mu_1 = xi''(1/2)/(2 xi(1/2)) (Cauchy integral): density
# sum_{k<=K} C_k lam_k exp(-lam_k (x - r_K)) on x > r_K with C_k = prod_{j != k} (1 - lam_k/lam_j)^{-1}; the sup over a grid of x in [0, 0.4] (spacing 5e-4, plus r_K and
# the 40 nearest points around it on both sides) of |spline - f|. Dyck chain: exact distribution of the height H_n of a uniform Dyck path of length 2n by
# the reflection formula (de Bruijn, Knuth and Rice 1972): #{height <= h} = sum_{k in Z} [C(2n, n + k(h+2)) - C(2n, n + k(h+2) - 1)], Y_n = H_n/sqrt(pi n),
# F_n(s) = sum_h P(H_n = h) (h/sqrt(pi n))^s; the zero near 1/2 + i gamma_1 by Newton's method in ball arithmetic at 400 bits; the extrapolation is done in two ways
# (powers of n^(-1/2) and of n^(-1)) because the chapter does not state its method.
# Output: data/ev_g1_tail.json and data/ev_g1_tail.log (--out FILE).
# Usage: python3 ev_g1_tail.py [--quick] [--out FILE] [--part talbot|bspline|dyck|all]; each part runs in seconds to about 25 seconds on a laptop; cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g ev_g1_tail.json -g ev_g1_tail.log code/ev_g1_tail.py data/zeros_1700.txt -- --out ev_g1_tail.json
import argparse, json, math, os, sys, time
from flint import arb, acb, ctx

START = time.time()
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


def xi(s):
    pi = arb.pi()
    return s * (s - 1) / 2 * ((-s / 2) * pi.log()).exp() * (s / 2).gamma() * s.zeta()


def Lfun(s, x0):
    return x0 / xi(acb(0.5) + s.sqrt())


def talbot(F, t, M):
    """f(t) = (2/(5t)) sum_{k<M} Re(gamma_k F(delta_k/t))."""
    pi = arb.pi()
    t = arb(t)
    tot = arb(0)
    d0 = arb(2 * M) / 5
    tot += (acb(d0).exp() / 2 * F(acb(d0) / t)).real
    for k in range(1, M):
        th = pi * k / M
        cot = th.cos() / th.sin()
        d = acb(2 * k * pi / 5 * cot, 2 * k * pi / 5)
        g = acb(1, th * (1 + cot * cot) - cot) * d.exp()
        tot += (g * F(d / t)).real
    return 2 * tot / (5 * t)


def load_zeros(path):
    return [arb(l.split()[0]) for l in open(path) if l.strip()]


def run_talbot(args, gam, out):
    pi = arb.pi()
    res = {}
    for (M, bits) in ((100, 700), (130, 900)) if not args.quick else ((60, 500), (80, 600)):
        ctx.prec = bits
        x0 = xi(acb(0.5)).real
        FL = lambda s: Lfun(s, x0)
        FS = lambda s: (1 - Lfun(s, x0)) / s
        g1 = gam[0]
        rows = []
        for x in ("0.18", "0.19", "0.20"):
            xv = arb(x)
            f = talbot(FL, xv, M)
            S = talbot(FS, xv, M)
            gh = (f / S).sqrt()
            rows.append((x, gh))
            log(f"[talbot] M={M} bits={bits} x={x}: sqrt(f/S) = {gh.str(30, radius=False)}  gamma_1 = {g1.str(30, radius=False)}  diff = {(gh - g1).str(3)}  (S = {S.str(10, radius=False)})")
        res[f"M{M}"] = [(x, gh.str(35, radius=False)) for x, gh in rows]
    ctx.prec = 700
    # deflation scan
    M = 100
    x0 = xi(acb(0.5)).real
    FL = lambda s: Lfun(s, x0)
    FS = lambda s: (1 - Lfun(s, x0)) / s
    xa = arb("0.19")
    S_a = talbot(FS, xa, M)
    lam1 = (talbot(FL, xa, M) / S_a)
    c1 = S_a * (lam1 * xa).exp()
    log(f"[deflate] lambda_1 estimate {lam1.str(25, radius=False)} (gamma_1^2 = {(gam[0] ** 2).str(25, radius=False)}), c_1 estimate {c1.str(15, radius=False)}")
    ests = []
    for xs in ("0.02", "0.03", "0.04", "0.05", "0.06", "0.08", "0.10", "0.12"):
        x = arb(xs)
        S = talbot(FS, x, M)
        f = talbot(FL, x, M)
        S1 = S - c1 * (-lam1 * x).exp()
        f1 = f - c1 * lam1 * (-lam1 * x).exp()
        g2 = (f1 / S1).sqrt()
        ests.append((xs, g2))
        log(f"[deflate] x={xs}: second atom estimate sqrt(f1/S1) = {g2.str(12, radius=False)} (gamma_2 = {gam[1].str(12, radius=False)}; chapter prints 21.0196)")
    out["talbot"] = res
    out["deflate"] = [(x, g.str(15, radius=False)) for x, g in ests]


def run_bspline(args, gam, out):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from ev_g1_string import xi_dd_half
    ctx.prec = 300
    mu1, x0 = xi_dd_half()
    M = 40 if not args.quick else 30
    ctx.prec = 400
    FL = lambda s: Lfun(s, x0)
    xs = [arb(5) / 10000 * j for j in range(6, 801)]
    if args.quick:
        xs = [arb(2) / 1000 * j for j in range(2, 101)]
    # density of X by Talbot at grid points and at the knots r_K +- small
    lam = [g * g for g in gam]
    Ks = (1, 2, 5, 20, 40)
    rK = {K: mu1 - sum((1 / lam[k] for k in range(K)), arb(0)) for K in Ks}
    pts = list(xs)
    for K in Ks:
        for j in range(-20, 21):
            pts.append(rK[K] + arb(j) / 4000)
    pts = [p for p in pts if p > arb(3) / 1000]
    t0 = time.time()
    dens = {}
    for p in pts:
        key = p.str(25, radius=False)
        dens[key] = talbot(FL, p, M)
    log(f"[bspline] Talbot (M={M}, 300 bits) density of X at {len(pts)} points in {time.time() - t0:.1f}s; mu_1 = {mu1.str(15, radius=False)}")
    pk_x, peak = max(dens.items(), key=lambda kv: float(kv[1].mid()))
    log(f"[bspline] peak density on the grid: {peak.str(8, radius=False)} at x = {pk_x[:8]}")
    res = {}
    for K in Ks:
        C = []
        for k in range(K):
            c = arb(1)
            for j in range(K):
                if j != k:
                    c *= 1 / (1 - lam[k] / lam[j])
            C.append(c)

        def spline(x):
            if x < rK[K]:
                return arb(0)
            return sum((C[k] * lam[k] * (-lam[k] * (x - rK[K])).exp() for k in range(K)), arb(0))
        best, bx = arb(0), None
        for p in pts:
            d = abs(spline(p) - dens[p.str(25, radius=False)])
            if d > best:
                best, bx = d, p
        bu = arb(0)
        for p in xs:
            if p.str(25, radius=False) not in dens:
                continue
            d = abs(spline(p) - dens[p.str(25, radius=False)])
            if d > bu:
                bu = d
        log(f"[bspline] K={K}: sup over the uniform grid only (spacing 5e-4, the knot r_K itself excluded) = {bu.str(5, radius=False)}")
        # value just right of the knot (jump of K = 1)
        jump = spline(rK[K]) - dens[(rK[K]).str(25, radius=False)] if rK[K] > 0 and (rK[K]).str(25, radius=False) in dens else None
        log(f"[bspline] K={K}: r_K = {rK[K].str(8, radius=False)}; sup |spline - f| on the grid = {best.str(5, radius=False)} at x = {bx.str(6, radius=False)} (peak of f: {peak.str(5, radius=False)})")
        res[f"K{K}"] = dict(sup_uniform=bu.str(8, radius=False), sup=best.str(8, radius=False), at=bx.str(8, radius=False), rK=rK[K].str(10, radius=False))
    out["bspline"] = res


def dyck_zero(n, s0, bits=400):
    ctx.prec = bits
    N2 = 2 * n
    row = [1]
    for j in range(N2):
        row.append(row[-1] * (N2 - j) // (j + 1))
    cat = row[n] // (n + 1)

    def C(j):
        return row[j] if 0 <= j <= N2 else 0

    def cnt_le(h):
        if h < 0:
            return 0
        w = h + 2
        kmax = (n + 1) // w + 2
        return sum(C(n + k * w) - C(n + k * w - 1) for k in range(-kmax, kmax + 1))
    pi = arb.pi()
    ph = []
    prev = 0
    hs = []
    hmax = int(14 * math.sqrt(n)) + 10
    for h in range(1, min(hmax, n) + 1):
        c = cnt_le(h)
        ph.append(arb(c - prev) / arb(cat))
        hs.append(arb(h) / (pi * n).sqrt())
        prev = c
    tot_p = sum(ph, arb(0))
    ln_y = [y.log() for y in hs]

    def F(s):
        return sum((p * (s * l).exp() for p, l in zip(ph, ln_y)), acb(0))

    def dF(s):
        return sum((p * l * (s * l).exp() for p, l in zip(ph, ln_y)), acb(0))
    if s0 is None:
        best = None
        for re in [0.05 * i for i in range(0, 21)]:
            for im in [12.0 + 0.1 * j for j in range(0, 41)]:
                v = abs(F(acb(re, im))).mid()
                if best is None or float(v) < best[0]:
                    best = (float(v), re, im)
        s0 = acb(best[1], best[2])
        for _ in range(5):
            pass
    s = acb(s0)
    for _ in range(30):
        step = F(s) / dF(s)
        s = s - step
        s = acb(s.real.mid(), s.imag.mid())
        if abs(step).mid() < arb(10) ** -60:
            break
    return s, tot_p, abs(step)


def run_dyck(args, out):
    ns = [100, 200, 400, 800, 1600, 3200, 6400, 12800] if not args.quick else [100, 400, 1600]
    zs = []
    s0 = None
    t0 = time.time()
    for n in ns:
        s, tp, st = dyck_zero(n, s0)
        s0 = s
        zs.append((n, s))
        log(f"[dyck] n={n}: sum of P(H=h) = {tp.str(12, radius=False)}; zero of F_n near 1/2 + 14.13i: {s.real.str(10, radius=False)} + {s.imag.str(10, radius=False)} i  (Newton step {st.str(2)}; {time.time() - t0:.1f}s)")
    out["dyck"] = [(n, s.real.str(15, radius=False), s.imag.str(15, radius=False)) for n, s in zs]
    # extrapolation in n^(-1/2) and n^(-1) from the last three doubled n
    if len(zs) >= 3:
        for name, p in (("n^(-1/2)", 0.5), ("n^(-1)", 1.0)):
            a = [z for n, z in zs[-3:]]
            q = arb(2) ** p
            # Richardson with two orders: error ~ c1 n^-p + c2 n^-2p
            r1 = [(q * a[1] - a[0]) / (q - 1), (q * a[2] - a[1]) / (q - 1)]
            r2 = (q * q * r1[1] - r1[0]) / (q * q - 1)
            log(f"[dyck] Richardson in powers of {name} from the last three n ({[n for n, z in zs[-3:]]}): one step {r1[1].real.str(8, radius=False)} + {r1[1].imag.str(8, radius=False)} i, two steps {r2.real.str(8, radius=False)} + {r2.imag.str(8, radius=False)} i")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--part", default="all")
    ap.add_argument("--zeros", default=None)
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    zp = args.zeros or next(p for p in (os.path.join(here, "..", "data", "zeros_1700.txt"), "zeros_1700.txt") if os.path.exists(p))
    ctx.prec = 700
    gam = load_zeros(zp)
    out = {"quick": args.quick}
    if args.part in ("talbot", "all"):
        run_talbot(args, gam, out)
    if args.part in ("bspline", "all"):
        run_bspline(args, gam, out)
    if args.part in ("dyck", "all"):
        run_dyck(args, out)
    out["seconds"] = time.time() - START
    outp = args.out or ("/tmp/ev_g1_tail_quick.json" if args.quick else os.path.join(here, "..", "data", "ev_g1_tail.json"))
    json.dump(out, open(outp, "w"), indent=1, default=str)
    open(os.path.splitext(outp)[0] + ".log", "w").write("\n".join(LOG) + "\n")
    log(f"done in {time.time() - START:.1f}s; wrote {outp}")


if __name__ == "__main__":
    main()
