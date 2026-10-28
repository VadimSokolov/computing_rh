# bsy_cuts.py: check of Proposition bsy:prop:cuts (complex moments through the cuts, Chapter 35, Balazard, Saias and Yor meet Thorin, subsection "The defect of the complex
# moments") on the Davenport and Heilbronn function f (h = f), group BSY-G, task N3 of misc/reviews/nick-b4/reconcile-M3.md. Claims supported: the jump of L across every cut is 2 pi i;
# at theta = 1/2 and 3/4, t0 = 85.7: D = E Z - z0 = sum over the cuts of int_C K ds, E Z^2 - z0^2 = 2 sum int_C (L_- + pi i) K ds and the balance defect
# Var X - Var Y + 2 i Cov(X, Y) = (E Z^2 - z0^2) - 2 z0 D - D^2, against the DH rows of Table bsy:tab:cm (data/bsy_cm.json, quadrature of the moments of Z, an independent computation).
# Conventions of the Proposition: z = 1 + i t0, z* = 2 theta - 1 + i t0, K(s) = 1/(s - z*) - 1/(s - z), cuts C = [theta + i gamma, beta + i gamma] for the (simple) zeros rho = beta + i gamma
# of f with beta > theta, below height 700, and their conjugates (f has real Dirichlet coefficients); L = log f continued along horizontal paths from Re s = +infinity (L -> 0),
# L_- the boundary value from below on the cut; the value above is L_- + 2 pi i. Zeros: the ten zeros with beta > 1/2 of data/bsy_defects.json (below height 700, found there by the scan of
# bsy_defects.py), each refined here by the secant method in mpmath at 30 digits from f in its definition by Hurwitz zeta functions (residual |f(rho)| logged); none with beta > 1 is known below 700
# (the table of defects scans sigma in [0.52, 1.40]); a zero with beta >= 1 on a path would show as a small minimum of |f|, which is logged (minimum of |f| on every path).
# Method. f'/f from bsy_zdh.py (Euler and Maclaurin, double precision), imported, not copied. L at b + R + i gamma (b = beta, R = 0.03) is Log f(sigma0 + i gamma) minus the integral of f'/f
# from b + R to sigma0 = 45 (|Log f| < 1e-14 there; composite Gauss and Legendre, panels of width 1). With g(s) = f'/f(s) - 1/(s - rho) analytic near rho, g on [b - R, b + R] is
# interpolated by a Chebyshev polynomial of degree 23 from 24 nodes (the nodes stay at least 2e-3 from rho, so the cancellation in g is harmless; rho is known to 1e-28), and then
#   L_-(x) + pi i = ell(x) + log(b - x), ell(x) = L(b + R) - log R - int_x^(b+R) g (the continuation along a path below rho: arg(s - rho) goes from 0 to -pi),
# which is the formula checked here. Left of b - R, ell is continued by Gauss and Legendre (panels of width 0.05, 16 nodes). The cut integrals use u = b - x and geometric panels
# [u_max 2^-(j+1), u_max 2^-j] (16 nodes, down to u = 1e-22; the neglected piece is below 1e-20) for the log factor, plain panels beyond R. D is also computed in closed form,
# the sum of log|(rho - z*)/(rho - z)| + i [Arg(...)] differences along the straight cut (no pole of K on a cut). z0 = L(1 + i t0) by the same path integral.
# Checks. (1) Rectangle: for three points x on every cut, L(x + i(gamma + eta)) - L(x + i(gamma - eta)) - int_(x - i eta)^(x + i eta) f'/f ds = 2 pi i (argument principle, eta = 0.005,
# direct line integrals, no g, no interpolation); (2) the formula for L_- at the same points against the direct value on the lower line plus the vertical integral up to the cut;
# (3) everything repeated with other parameters (evaluator nf 0.45, M 36; R = 0.02; degree 19; panels 0.04; 20 nodes; sigma0 = 50) and the differences logged; (4) z0 against the arg_ray
# value of bsy_cm.py (data/bsy_cm.json x0, y0). Output: data/bsy_cuts.json and the log data/bsy_cuts.log. Usage: python3 -I bsy_cuts.py [--quick] [--out FILE]; --quick (3 zeros, about
# 25 s) writes to /tmp unless --out is given; the full case is a few minutes (run it on the cluster if it exceeds 30 s locally):
#   bash misc/tools/hopper_run.sh -c 2 -t 30 -m 4G -n bsy_cuts -g bsy_cuts.json -g bsy_cuts.log code/bsy_cuts.py code/bsy_zdh.py data/bsy_defects.json data/bsy_cm.json -- --out bsy_cuts.json
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse
import json
import math
import os
import sys
import time
import numpy as np
import mpmath as mp
from numpy.polynomial import chebyshev as Ch
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsy_zdh as Z

T0 = 85.7
PI = math.pi


def find(name):
    here = os.path.dirname(os.path.abspath(__file__))
    return next(q for q in (os.path.join(here, name), os.path.join(here, "..", "data", name)) if os.path.exists(q))


def fmp(s):
    k = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    h = lambda a: mp.zeta(s, mp.mpf(a) / 5)
    return mp.mpf(5) ** (-s) * (h(1) - h(4) + k * (h(2) - h(3)))


def scan_zeros(t_lo, t_hi, kw=None):
    """Zeros of f with 1/2 < beta < 1.4 and t_lo < gamma <= t_hi: grid (step 0.05 in t, 0.02 in sigma), local minima of |f| below 0.5, Newton in double precision."""
    kw = kw or {}
    ts = np.arange(t_lo, t_hi, 0.05)
    sg = np.arange(0.52, 1.40 + 1e-9, 0.02)
    F = np.empty((len(ts), len(sg)))
    for j, sj in enumerate(sg):
        F[:, j] = np.abs(Z.dh(sj + 1j * ts, **kw))
    c = (F[1:-1, 1:-1] < F[:-2, 1:-1]) & (F[1:-1, 1:-1] < F[2:, 1:-1]) & (F[1:-1, 1:-1] < F[1:-1, :-2]) & (F[1:-1, 1:-1] < F[1:-1, 2:]) & (F[1:-1, 1:-1] < 0.5)
    ii, jj = np.nonzero(c)
    found = []
    for z in sg[jj + 1] + 1j * ts[ii + 1]:
        s_ = np.array([z])
        for _ in range(30):
            f, d = Z.dh(s_, deriv=True)
            s_ = s_ - f / d
        f, d = Z.dh(s_, deriv=True)
        if abs(f[0]) < 1e-11 and 0.5 + 1e-6 < s_[0].real < 1.45 and t_lo < s_[0].imag <= t_hi and not any(abs(s_[0] - q) < 1e-6 for q in found):
            found.append(complex(s_[0]))
    return sorted(found, key=lambda z: z.imag), int(c.sum())


def refine(seed, dps):
    mp.mp.dps = dps
    r = mp.findroot(fmp, (mp.mpc(seed) , mp.mpc(seed) + mp.mpf("1e-6")), solver="secant", tol=mp.mpf(10) ** (-dps + 4), maxsteps=40)
    return r, abs(fmp(r))


class Run:
    """One parameter set of the whole computation."""

    def __init__(self, kw, R, deg_n, wid, ngl, sigma0, log):
        self.kw, self.R, self.n, self.wid, self.ngl, self.s0, self.log = kw, R, deg_n, wid, ngl, sigma0, log
        self.xg, self.wg = leggauss(ngl)
        self.minabs = 1e9

    def fp(self, s):
        f, d = Z.dh(np.atleast_1d(np.asarray(s, dtype=complex)), deriv=True, **self.kw)
        self.minabs = min(self.minabs, float(np.abs(f).min()))
        return d / f

    def gl(self, a, b):
        return (a + b) / 2 + (b - a) / 2 * self.xg, (b - a) / 2 * self.wg

    def Lline(self, tau, x, d0=0.15, center=None, eta=None):
        """L(x + i tau) along the line Im s = tau from sigma0: Log f(sigma0 + i tau) - int_x^sigma0 f'/f. Panels: doubling from x (first width d0, so a pole at distance
        d0 left of x or more is resolved), unit panels, and, if a zero is at real part center and distance eta from the line, panel edges center +- eta 2^k."""
        f0, _ = Z.dh(np.array([self.s0 + 1j * tau]), deriv=True, **self.kw)
        L = np.log(f0[0])
        e = [x + d0 * (2.0 ** np.arange(0, 12) - 1), np.arange(x, self.s0, 1.0), [self.s0]]
        if center is not None:
            k = 2.0 ** np.arange(-1, 11)
            e += [center + eta * k, center - eta * k, [center]]
        edges = np.unique(np.concatenate(e))
        edges = edges[(edges >= x) & (edges <= self.s0)]
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a < 1e-14:
                continue
            xs, ws = self.gl(a, b)
            L = L - (ws * self.fp(xs + 1j * tau)).sum()
        return L

    def intg(self, g, a, b):
        xs, ws = self.gl(a, b)
        return (ws * g(xs)).sum()

    def zero(self, rho_mp):
        """Data of one zero (upper half plane): g, the Chebyshev antiderivative and ell."""
        b_hi = float(rho_mp.real)
        b_lo = float(rho_mp.real - b_hi)
        g_hi = float(rho_mp.imag)
        g_lo = float(rho_mp.imag - g_hi)
        R = self.R

        def g(xs):
            s = xs + 1j * g_hi
            return self.fp(s) - 1 / ((xs - b_hi - b_lo) - 1j * g_lo)

        Lb = self.Lline(g_hi, b_hi + R, d0=R, center=b_hi, eta=R)
        k = np.arange(self.n)
        t = np.cos(PI * (k + 0.5) / self.n)
        c = Ch.chebfit(t, g(b_hi + R * t), self.n - 1)
        ci = Ch.chebint(c, lbnd=1)
        # second fit of one degree more as a check of the interpolation
        t2 = np.cos(PI * (np.arange(self.n + 4) + 0.5) / (self.n + 4))
        c2 = Ch.chebfit(t2, g(b_hi + R * t2), self.n + 3)
        ci2 = Ch.chebint(c2, lbnd=1)
        tt = np.linspace(-1, 1, 7)
        interp_diff = float(np.abs(R * (Ch.chebval(tt, ci) - Ch.chebval(tt, ci2))).max())

        def ell_cheb(x):  # x in [b - R, b + R]
            return Lb + R * Ch.chebval((x - b_hi) / R, ci) - math.log(R)

        zd = {"b": b_hi, "gamma": g_hi, "g": g, "ell_cheb": ell_cheb, "Lb": Lb, "interp_diff": interp_diff}

        def ell(x):  # any x <= b + R, scalar
            if x >= b_hi - R:
                return ell_cheb(x)
            e = ell_cheb(b_hi - R)
            hi = b_hi - R
            while hi > x + 1e-15:
                lo = max(x, hi - self.wid)
                e = e - self.intg(g, lo, hi)
                hi = lo
            return e

        zd["ell"] = ell
        return zd

    def rule(self, zd, theta):
        """Nodes and weights x (real part), u = b - x, ell(x) of the cut [theta, b] (upper half plane)."""
        b, R = zd["b"], self.R
        umax = min(R, b - theta)
        X, W, Ell, U = [], [], [], []
        edges = umax * 0.5 ** np.arange(0, 74)
        for hi, lo in zip(edges[:-1], edges[1:]):
            us, ws = self.gl(lo, hi)
            X.append(b - us)
            U.append(us)
            W.append(ws)
            Ell.append(zd["ell_cheb"](b - us))
        if b - theta > R:
            ehi = b - R
            ell_hi = zd["ell_cheb"](ehi)
            while ehi > theta + 1e-15:
                elo = max(theta, ehi - self.wid)
                xs, ws = self.gl(elo, ehi)
                e = np.array([ell_hi - self.intg(zd["g"], x, ehi) for x in xs])
                X.append(xs)
                U.append(b - xs)
                W.append(ws)
                Ell.append(e)
                ell_hi = ell_hi - self.intg(zd["g"], elo, ehi)
                ehi = elo
        return np.concatenate(X), np.concatenate(W), np.concatenate(Ell), np.concatenate(U)


def Kf(s, theta):
    zs = 2 * theta - 1 + 1j * T0
    z = 1 + 1j * T0
    return 1 / (s - zs) - 1 / (s - z)


def cutsum(run, zeros, theta, z0):
    """D (quadrature), D (closed form), E Z^2 - z0^2 and the balance defect for the cuts with beta > theta."""
    Dq = 0j
    Dc = 0j
    S1 = 0j
    ncut = 0
    zs = 2 * theta - 1 + 1j * T0
    z = 1 + 1j * T0
    for rho, zd in zeros:
        if zd["b"] <= theta:
            continue
        ncut += 1
        X, W, Ell, u = run.rule(zd, theta)
        br = Ell + np.log(u)  # L_- + pi i
        for sgn in (1, -1):
            s = X + sgn * 1j * zd["gamma"]
            Ks = Kf(s, theta)
            bb = br if sgn == 1 else np.conj(br)
            Dq += (W * Ks).sum()
            S1 += (W * bb * Ks).sum()
            r = complex(zd["b"], sgn * zd["gamma"])
            l = complex(theta, sgn * zd["gamma"])
            Dc += math.log(abs((r - zs) / (r - z)) / abs((l - zs) / (l - z))) + 1j * (np.angle((r - zs) / (l - zs)) - np.angle((r - z) / (l - z)))
    EZ2 = 2 * S1
    bal = EZ2 - 2 * z0 * Dq - Dq ** 2
    return {"cuts": ncut, "D_quad": Dq, "D_closed": Dc, "EZ2_minus_z0sq": EZ2, "balance": bal}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--extend", type=float, default=0.0)
    A = ap.parse_args()
    if A.extend and not A.out:
        ap.error("--extend needs --out")
    outp = A.out or ("/tmp/bsy_cuts_quick.json" if A.quick else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "bsy_cuts.json"))
    lf = open(os.path.splitext(outp)[0] + ".log", "w")

    def log(m):
        print(m, flush=True)
        lf.write(m + "\n")
        lf.flush()

    t00 = time.time()
    dps = 20 if A.quick else 30
    seeds = [complex(x) for x in json.load(open(find("bsy_defects.json")))["dh_off"]]
    if A.quick:
        seeds = seeds[:3]
    if A.extend:
        ext, nc = scan_zeros(700.0, A.extend)
        log(f"extension: zeros of f with 1/2 < beta < 1.4 and 700 < gamma <= {A.extend:g}: {len(ext)} found from {nc} candidate minima: " + ", ".join(f"{z.real:.4f}+{z.imag:.3f}i" for z in ext))
        seeds = seeds + ext
    log(f"bsy_cuts.py {'quick' if A.quick else 'full'}: DH, t0 = {T0}; {len(seeds)} zeros with beta > 1/2 (data/bsy_defects.json), refined in mpmath at {dps} digits")
    rhos = []
    for sd in seeds:
        r, res = refine(sd, dps)
        rhos.append(r)
        log(f"  rho = {mp.nstr(r.real, 20)} + {mp.nstr(r.imag, 22)} i, |f(rho)| = {mp.nstr(res, 2)}; seed difference {abs(complex(r) - sd):.1e}")
    runs = [Run({}, 0.03, 24, 0.05, 16, 45.0, log), Run({"nf": 0.45, "M": 36}, 0.02, 20, 0.04, 20, 50.0, log)]
    out = {"t0": T0, "zeros": [[mp.nstr(r.real, 25), mp.nstr(r.imag, 25)] for r in rhos], "runs": []}
    ref = json.load(open(find("bsy_cm.json")))["rows"]
    for ri, run in enumerate(runs):
        t1 = time.time()
        log(f"run {ri + 1}: evaluator {run.kw or 'default'}, R {run.R}, Chebyshev nodes {run.n}, panel width {run.wid}, Gauss nodes {run.ngl}, sigma0 {run.s0}")
        zeros = [(r, run.zero(r)) for r in rhos]
        for r, zd in zeros:
            log(f"  zero {zd['b']:.6f} + {zd['gamma']:.6f} i: L(b + R + i gamma) = {zd['Lb']:.12f}, degree test (n vs n + 4) {zd['interp_diff']:.1e}")
        res = {"theta": {}}
        # z0 and its comparison with the arg_ray value of bsy_cm.py
        z0 = run.Lline(T0, 1.0)
        rowref = {0.5: ref[4], 0.75: ref[5]}
        log(f"  z0 = L(1 + i t0) = {z0.real:.10f} {z0.imag:+.10f} i (bsy_cm.py arg_ray: {ref[4]['x0']:.10f} {ref[4]['y0']:+.10f} i; difference {abs(z0 - complex(ref[4]['x0'], ref[4]['y0'])):.1e})")
        res["z0"] = [z0.real, z0.imag]
        for theta, H in [(th, hh) for th in (0.5, 0.75) for hh in ([None] if not A.extend else [700.0, A.extend])]:
            zsel = zeros if H is None else [(r, zd) for r, zd in zeros if zd["gamma"] <= H]
            c = cutsum(run, zsel, theta, z0)
            rf = rowref[theta]["fine"]
            tab_D = complex(rf["EX-x0"], rf["EY-y0"])
            tab_B = complex(rf["VarX"] - rf["VarY"], 2 * rf["Cov"])
            key = str(theta) if H is None else f"{theta}_H{H:g}"
            res["theta"][key] = {k: ([v.real, v.imag] if isinstance(v, complex) else v) for k, v in c.items()}
            res["theta"][key]["table_D"] = [tab_D.real, tab_D.imag]
            res["theta"][key]["table_balance"] = [tab_B.real, tab_B.imag]
            log(f"  theta {theta}" + ("" if H is None else f", zeros below height {H:g}") + f": {c['cuts']} zeros with beta > theta (and conjugates)")
            log(f"    D = E Z - z0: quadrature {c['D_quad'].real:+.8f} {c['D_quad'].imag:+.8f} i, closed form {c['D_closed'].real:+.8f} {c['D_closed'].imag:+.8f} i (difference {abs(c['D_quad'] - c['D_closed']):.1e}); Table cm {tab_D.real:+.4f} {tab_D.imag:+.4f} i")
            log(f"    E Z^2 - z0^2 = {c['EZ2_minus_z0sq'].real:+.8f} {c['EZ2_minus_z0sq'].imag:+.8f} i")
            log(f"    balance defect Var X - Var Y + 2i Cov = {c['balance'].real:+.8f} {c['balance'].imag:+.8f} i; Table cm {tab_B.real:+.4f} {tab_B.imag:+.4f} i (difference {abs(c['balance'] - tab_B):.1e})")
        # rectangle check and the formula for L_- on theta = 1/2 cuts
        eta = 0.005
        worst_j, worst_f = 0.0, 0.0
        for r, zd in zeros:
            b, g = zd["b"], zd["gamma"]
            for u in (0.003, 0.02, 0.1):
                x = b - u
                if x <= 0.5:
                    continue
                Lup = run.Lline(g + eta, x, d0=eta, center=b, eta=eta)
                Llo = run.Lline(g - eta, x, d0=eta, center=b, eta=eta)
                vs, vw = run.gl(-eta, eta)
                vert = 1j * (vw * run.fp(x + 1j * (g + vs))).sum()
                J = Lup - Llo - vert
                worst_j = max(worst_j, abs(J - 2j * PI))
                hs, hw = run.gl(-eta, 0.0)
                low_to_cut = Llo + 1j * (hw * run.fp(x + 1j * (g + hs))).sum()
                form = zd["ell"](x) + math.log(u) - 1j * PI
                worst_f = max(worst_f, abs(form - low_to_cut))
        res["max_jump_error"] = worst_j
        res["max_Lminus_error"] = worst_f
        log(f"  rectangle check (jump 2 pi i, argument principle): max |J - 2 pi i| = {worst_j:.1e}; formula for L_- against the lower line continued to the cut: max difference {worst_f:.1e}; minimum of |f| on all paths {run.minabs:.2e}; run time {time.time() - t1:.0f} s")
        out["runs"].append(res)
    a, b = out["runs"]
    for th in (["0.5", "0.75"] if not A.extend else [f"{t}_H{h:g}" for t in (0.5, 0.75) for h in (700.0, A.extend)]):
        dD = abs(complex(*a["theta"][th]["D_quad"]) - complex(*b["theta"][th]["D_quad"]))
        dB = abs(complex(*a["theta"][th]["balance"]) - complex(*b["theta"][th]["balance"]))
        out.setdefault("second_run_difference", {})[th] = {"D": dD, "balance": dB}
        log(f"second parameter set against the first, theta {th}: |dD| {dD:.1e}, |d balance| {dB:.1e}")
    json.dump(out, open(outp, "w"), indent=1)
    log(f"done in {time.time() - t00:.0f} s; wrote {outp}")


if __name__ == "__main__":
    main()
