# ev_g1_string.py: the numerical claims of Nick Polson's extreme value hitting times chapter (E; Chapter 25 after the merge) about the string test and
# the shifted Bessel strings: the Hankel determinants D_n, D_n' and mu_1 of Section sec:string (E:224 to E:243, claim at E:241), the finite difference
# and Bessel zeros of the Liouville string (E:220, prop:liouville), the partial fractions of tab:shift (E:359 to E:374) and the eigenvalue lists of
# E:357 (zeros of Re K_{mu + i nu}(a)). Added for the book in October 2026 by computation group G1, written from the mathematics (the collaborator's
# scripts were not read, run or copied).
# Conventions. Phi_exact(u) = sum_n phi_n(u), phi_n = (2a^2 e^{9u} - 3a e^{5u}) exp(-a e^{4u}), a = pi n^2; folded N = 2: phi_1 + phi_2; Polya's kernel
# (4 pi^2 cosh 9u - 6 pi cosh 5u) exp(-2 pi cosh 4u); all on u >= 0. M_{2m} = int_0^inf u^{2m} Phi du, h(s) = sum_m 4^m M_{2m} s^m/((2m)! M_0) =
# H(2i sqrt s)/H(0); log h(s) = sum_j (-1)^{j+1} mu_j s^j / j; D_n = det[mu_{i+j+1}], D_n' = det[mu_{i+j+2}] (n x n, i, j = 0..n-1).
# For the shifted strings: K_{nu+1} = K_{nu-1} + (2 nu/x) K_nu, K_{i v + 1} = (i v K_{iv} - x K_{iv}')/x; Re K_{mu + i v}(a) = A_mu(v) K_{iv}(a) + B_mu(v) a K_{iv}'(a)
# with real polynomials, even in v, so polynomials in lambda = v^2; beta_mu(lambda) = -A_mu/B_mu = alpha - m lambda - sum_j sigma_j/(c_j - lambda) (polynomial
# division and residues at the roots of B_mu); the extension of prop:ext is positive in the chapter's sense when m >= 0 and all sigma_j > 0.
# Method. Moments by Gauss-Legendre quadrature on [0, U] in ball arithmetic (python-flint, arb.legendre_p_root) with two node counts; every D_n is computed
# for both and the relative change is logged (a quadrature error is not certified: the status is 'computed'). Real zeros of K_{iz}(2 pi) and of
# Re K_{mu+i nu}(a) from sign changes with certified signs on a grid, refined by bisection. Finite differences: second order on (0, 6), Dirichlet at both ends,
# a = 2 pi, numpy/scipy tridiagonal eigenvalues, for several mesh sizes (the chapter does not give its mesh).
# Output: data/ev_g1_string.json and data/ev_g1_string.log (--out FILE).
# Usage: python3 ev_g1_string.py [--quick] [--out FILE] [--part hankel|liouville|shift|all]; quick: 6 determinants and 150 digits in a second or two;
# full case in about 20 seconds on a laptop; cluster: bash misc/tools/hopper_run.sh -c 1 -t 30 -g ev_g1_string.json -g ev_g1_string.log code/ev_g1_string.py -- --out ev_g1_string.json
import argparse, json, math, os, sys, time
from flint import arb, acb, ctx, arb_poly, acb_poly, arb_mat

START = time.time()
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


def pi_():
    return arb.pi()


# ---------------- Hankel ----------------
def kernel(name):
    pi = pi_()
    if name == "exact":
        def Phi(u):
            tot = arb(0)
            n = 1
            e4 = (4 * u).exp()
            e5 = (5 * u).exp()
            e9 = (9 * u).exp()
            while True:
                a = pi * n * n
                t = (2 * a * a * e9 - 3 * a * e5) * (-a * e4).exp()
                tot += t
                if n > 2 and abs(t).upper() < arb(2) ** (-ctx.prec - 30) and float(a * e4) > 40:
                    break
                n += 1
            return tot
    elif name == "folded2":
        def Phi(u):
            tot = arb(0)
            e4 = (4 * u).exp()
            e5 = (5 * u).exp()
            e9 = (9 * u).exp()
            for n in (1, 2):
                a = pi * n * n
                tot += (2 * a * a * e9 - 3 * a * e5) * (-a * e4).exp()
            return tot
    else:
        def Phi(u):
            c9 = ((9 * u).exp() + (-9 * u).exp()) / 2
            c5 = ((5 * u).exp() + (-5 * u).exp()) / 2
            c4 = ((4 * u).exp() + (-4 * u).exp()) / 2
            return (4 * pi * pi * c9 - 6 * pi * c5) * (-2 * pi * c4).exp()
    return Phi


def moments(Phi, M, nodes, U):
    xs = []
    for k in range(nodes):
        r, w = arb.legendre_p_root(nodes, k, weight=True)
        xs.append((U * (r + 1) / 2, w * U / 2))
    vals = [(u * u, w * Phi(u)) for u, w in xs]
    out = [arb(0)] * (M + 1)
    for u2, wf in vals:
        pw = wf
        for m in range(M + 1):
            out[m] = out[m] + pw
            pw = pw * u2
    return out


def mus_from_moments(Ms, J):
    c = [Ms[m] * arb(4) ** m / arb(math.factorial(2 * m)) / Ms[0] for m in range(J + 1)]
    l = [arb(0)] * (J + 1)
    for j in range(1, J + 1):
        # j c_j = sum_{i=1}^{j} i l_i c_{j-i}
        s = j * c[j]
        for i in range(1, j):
            s -= i * l[i] * c[j - i]
        l[j] = s / j      # c_0 = 1
    mu = [None] + [(-1) ** (j + 1) * j * l[j] for j in range(1, J + 1)]
    return mu


def hankel(mu, nmax):
    Ds, Dps = [], []
    for n in range(1, nmax + 1):
        Ds.append(arb_mat(n, n, [mu[i + j + 1] for i in range(n) for j in range(n)]).det())
        Dps.append(arb_mat(n, n, [mu[i + j + 2] for i in range(n) for j in range(n)]).det())
    return Ds, Dps


def xi_dd_half():
    """xi''(1/2)/(2 xi(1/2)) by a Cauchy integral on |s - 1/2| = 1 (128 points); xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2."""
    pi = pi_()

    def xi(s):
        return s * (s - 1) / 2 * ((-s / 2) * pi.log()).exp() * (s / 2).gamma() * s.zeta()
    n = 128
    acc = acb(0)
    for k in range(n):
        th = 2 * pi * k / n
        z = acb(th.cos(), th.sin())
        acc += xi(acb(0.5) + z) / (z * z)
    d2 = 2 * acc / n          # f''(1/2) = 2!/(2 pi i) oint f/(z-z0)^3 dz = 2 * mean(f z^{-2})
    return (d2 / (2 * xi(acb(0.5)))).real, xi(acb(0.5)).real


def run_hankel(args, out):
    nmax = 6 if args.quick else 11
    J = 2 * nmax + 1
    Mm = J
    bits = 700 if args.quick else 1700
    nodes_list = (160, 240) if args.quick else (420, 560)
    U = 2.2
    res = {}
    for name in ((args.kernel,) if args.kernel else ("exact", "folded2", "polya")):
        t0 = time.time()
        runs = []
        for nodes in nodes_list:
            ctx.prec = bits
            Phi = kernel(name)
            Ms = moments(Phi, Mm, nodes, arb(U))
            mu = mus_from_moments(Ms, J)
            Ds, Dps = hankel(mu, nmax)
            runs.append((Ms, mu, Ds, Dps))
        ctx.prec = 300
        (Ms, mu, Ds, Dps), (Ms2, mu2, Ds2, Dps2) = runs
        relch = max(float(abs((Ds[n] - Ds2[n]) / Ds[n]).upper()) for n in range(nmax))
        relch2 = max(float(abs((Dps[n] - Dps2[n]) / Dps[n]).upper()) for n in range(nmax))
        pos = all(d > 0 for d in Ds) and all(d > 0 for d in Dps)
        pos2 = all(d > 0 for d in Ds2) and all(d > 0 for d in Dps2)
        log(f"[hankel] {name}: {nmax} determinants, {bits} bits (~{int(bits * 0.30103)} digits), nodes {nodes_list}: all D_n and D_n' positive (ball excludes 0): run1 {pos}, run2 {pos2}; max relative change of D_n between the node counts {relch:.2e}, of D_n' {relch2:.2e}; time {time.time() - t0:.1f}s")
        log(f"[hankel] {name}: mu_1 = {mu[1].str(15, radius=False)}  mu_2 = {mu[2].str(12, radius=False)}  D_1 = {Ds[0].str(8, radius=False)}  D_{nmax} = {Ds[-1].str(6, radius=False)}  D'_{nmax} = {Dps[-1].str(6, radius=False)}")
        # the chapter's precision: 160 digits (about 532 bits); are the signs decided by the balls?
        ctx.prec = 600
        Phi = kernel(name)
        Msl = moments(Phi, Mm, nodes_list[1], arb(U))
        mul = mus_from_moments(Msl, J)
        Dl, Dpl = hankel(mul, nmax)
        okD = [bool(d > 0) for d in Dl]
        okDp = [bool(d > 0) for d in Dpl]
        log(f"[hankel] {name}: at 600 bits (~180 digits, all intermediate quantities in balls) D_n certified positive for n = " + ",".join(str(i + 1) for i, o in enumerate(okD) if o) + "; D_n' certified positive for n = " + ",".join(str(i + 1) for i, o in enumerate(okDp) if o) + f" (of 1..{nmax}); first undecided D_n: " + str(next((i + 1 for i, o in enumerate(okD) if not o), None)) + ", D_n': " + str(next((i + 1 for i, o in enumerate(okDp) if not o), None)))
        ctx.prec = 300
        res[name] = dict(mu1=mu[1].str(20, radius=False), all_positive=bool(pos and pos2), relchange=max(relch, relch2), nmax=nmax, Dn=[d.str(8, radius=False) for d in Ds], Dpn=[d.str(8, radius=False) for d in Dps])
        if name == "exact":
            ctx.prec = 300
            val, x0 = xi_dd_half()
            log(f"[hankel] exact kernel: xi''(1/2)/(2 xi(1/2)) by Cauchy integral = {val.str(15, radius=False)} (mu_1 above), M_0 = {Ms[0].str(15, radius=False)} vs xi(1/2)/8 = {(x0 / 8).str(15, radius=False)}")
            res[name]["mu1_xi"] = val.str(20, radius=False)
            # digits needed: smallest D_n over the run, sizes
            log("[hankel] exact kernel D_n sizes: " + ", ".join(f"{float(d.mid()):.2e}" for d in Ds))
    out["hankel"] = res


# ---------------- Liouville ----------------
def run_liouville(args, out):
    import numpy as np
    from scipy.linalg import eigh_tridiagonal
    ctx.prec = 200
    pi = pi_()
    x = 2 * pi

    def K(z):
        return x.bessel_k(acb(0, 1) * acb(z)) if False else acb(x).bessel_k(acb(0, 1) * acb(z))
    fun = lambda z: K(arb(z)).real
    zs = []
    prev = None
    for j in range(1, 421):
        z = 0.05 * j
        v = fun(z)
        if v.contains(0):
            continue
        sg = v > 0
        if prev is not None and sg != prev[1]:
            lo, hi = prev[0], z
            slo = fun(lo) > 0
            for _ in range(60):
                mid = (lo + hi) / 2
                vm = fun(mid)
                if vm.contains(0):
                    break
                if (vm > 0) == slo:
                    lo = mid
                else:
                    hi = mid
            zs.append((lo + hi) / 2)
        prev = (z, sg)
    log("[liouville] zeros of K_{iz}(2 pi) (certified signs, bisection): " + ", ".join(f"{z:.8f}" for z in zs[:6]) + " ...; number below 18: " + str(sum(1 for z in zs if z < 18)))
    res = {"zeros": zs[:8]}
    a = float(2 * pi)
    fd = {}
    for n in ((2000, 8000, 32000) if args.quick else (2000, 8000, 32000, 64000, 100000, 128000, 200000)):
        h = 6.0 / n
        t = h * np.arange(1, n)
        d = 2.0 / h ** 2 + a * a * np.exp(2 * t)
        e = -np.ones(n - 2) / h ** 2
        w = eigh_tridiagonal(d, e, select="i", select_range=(0, 4), eigvals_only=True)
        fd[n] = np.sqrt(w)
        log(f"[liouville] finite differences (second order, (0,6), n = {n}): sqrt(lambda_1..5) = " + ", ".join(f"{v:.6f}" for v in fd[n]))
    # Richardson from the two finest
    ns = sorted(fd)
    r = (4 * fd[ns[-1]] - fd[ns[-2]]) / 3 if ns[-1] == 2 * ns[-2] else None
    ex = np.array(zs[:5])
    log("[liouville] differences (FD n=" + str(ns[-1]) + " minus zero of K): " + ", ".join(f"{v:.2e}" for v in (fd[ns[-1]] - ex)))
    out["liouville"] = dict(zeros=zs[:8], fd={str(k): list(map(float, v)) for k, v in fd.items()})


# ---------------- tab:shift and E:357 ----------------
def AB(mu, a):
    """(A, B) as acb_poly in v with K_{iv+mu} = A K_{iv} + B (a K_{iv}')."""
    one = acb(1)
    I = acb(0, 1)
    A0, B0 = acb_poly([1]), acb_poly([0])
    A1, B1 = acb_poly([0, I / a]), acb_poly([-1 / a])
    if mu == 0:
        return A0, B0
    prev, cur = (A0, B0), (A1, B1)
    for j in range(1, mu):
        f = acb_poly([acb(2 * j) / a, 2 * I / a])
        nxt = (prev[0] + f * cur[0], prev[1] + f * cur[1])
        prev, cur = cur, nxt
    return cur


def real_lambda_poly(P):
    cs = P.coeffs()
    re = [c.real for c in cs]
    odd = max([abs(c.imag).upper() for c in cs[0::2]] + [abs(c.real).upper() for c in cs[1::2]] + [arb(0)])
    return arb_poly([re[k] for k in range(0, len(re), 2)]), odd


def run_shift(args, out):
    ctx.prec = 300
    pi = pi_()
    res = {}
    for a, aname in ((2 * pi, "2pi"), (arb(2), "2")):
        for mu in range(2, 9):
            A, B = AB(mu, a)
            Al, oa = real_lambda_poly(A)
            Bl, ob = real_lambda_poly(B)
            q, r = divmod(-Al, Bl)
            qc = q.coeffs() + [arb(0), arb(0)]
            alpha, m = qc[0], -qc[1]
            roots = [rt for rt in Bl.complex_roots()] if Bl.degree() > 0 else []
            poles = []
            dB = Bl.derivative()
            for rt in roots:
                c = rt[0] if isinstance(rt, tuple) else rt
                cr = c.real
                sig = (r(cr) / dB(cr))
                poles.append((float(cr.mid()), float(sig.mid())))
            poles.sort()
            positive = (float(m.mid()) >= -1e-12) and all(s > 0 for c, s in poles)
            log(f"[shift] a={aname} mu={mu}: odd part check {max(oa, ob).str(2)}; alpha = {float(alpha.mid()):.5g}, m = {float(m.mid()):.6g}, poles c_j = " + ", ".join(f"{c:.4g}" for c, s in poles) + "; sigma_j = " + ", ".join(f"{s:.4g}" for c, s in poles) + f"; positive (m >= 0 and all sigma_j > 0): {positive}")
            res[f"{aname}_{mu}"] = dict(m=float(m.mid()), alpha=float(alpha.mid()), poles=poles, positive=positive)
    out["shift"] = res
    # zeros of Re K_{mu + i nu}(a)
    zr = {}
    for mu, a, nm, count in ((3, arb(2), "2", 5), (4, 2 * pi, "2pi", 5), (6, 2 * pi, "2pi", 5), (3, 2 * pi, "2pi", 3)):
        ac = acb(a)
        fun = lambda v: ac.bessel_k(acb(mu) + acb(0, 1) * acb(v)).real
        zs, prev = [], None
        for j in range(1, 800):
            v = 0.02 * j
            val = fun(v)
            if val.contains(0):
                continue
            sg = val > 0
            if prev is not None and sg != prev[1]:
                lo, hi = prev[0], v
                slo = fun(lo) > 0
                for _ in range(50):
                    mid = (lo + hi) / 2
                    vm = fun(mid)
                    if vm.contains(0):
                        break
                    if (vm > 0) == slo:
                        lo = mid
                    else:
                        hi = mid
                zs.append((lo + hi) / 2)
                if len(zs) >= count:
                    break
            prev = (v, sg)
        log(f"[shift] zeros of Re K_{{{mu}+i nu}}({nm}) (first {count}): " + ", ".join(f"{z:.5f}" for z in zs))
        zr[f"mu{mu}_a{nm}"] = zs
    out["shift_zeros"] = zr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--part", default="all")
    ap.add_argument("--kernel", default=None)
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = {"quick": args.quick}
    if args.part in ("hankel", "all"):
        run_hankel(args, out)
    if args.part in ("liouville", "all"):
        run_liouville(args, out)
    if args.part in ("shift", "all"):
        run_shift(args, out)
    out["seconds"] = time.time() - START
    outp = args.out or ("/tmp/ev_g1_string_quick.json" if args.quick else os.path.join(here, "..", "data", "ev_g1_string.json"))
    json.dump(out, open(outp, "w"), indent=1, default=str)
    open(os.path.splitext(outp)[0] + ".log", "w").write("\n".join(LOG) + "\n")
    log(f"done in {time.time() - START:.1f}s; wrote {outp}")


if __name__ == "__main__":
    main()
