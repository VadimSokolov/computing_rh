# gp_skeleton.py: the numerical claims N1.1 to N1.5 of the zeros of the skeleton Xi^o(t) = Re K_{1/4+it/2}(2 pi) (Proposition pol:prop:skeleton, the
# paragraph "Numerically, Xi^o has ... zeros in (0,100]" after its proof; notes G and G2 of the merge): the number of zeros t^o_n in (0,100], the
# first zero, the range of theta(t^o_n)/pi - (n - 1/2), and the displacement t^o_n - gamma_(n+1) against the zeros of zeta (mean absolute value,
# Pearson correlation with pi S(gamma_(n+1))/theta'(gamma_(n+1))). Added for the book in October 2026, written from the mathematics (the
# collaborator's scripts were not read, run or copied).
# Conventions. t is the book's ordinate (the same unit as the ordinates gamma_n of zeta; K_nu(2 pi) at nu = 1/4 + i t/2). theta(t) = Im log
# Gamma(1/4 + i t/2) - (t/2) log pi is the Riemann Siegel theta function, N(T) = theta(T)/pi + 1 + S(T), and at an ordinate gamma_(n+1) (the
# (n+1)st zero, counted by ordinate) the average of the one sided limits gives S(gamma_(n+1)) = n - 1/2 - theta(gamma_(n+1))/pi. The zeros of the
# skeleton are indexed t^o_1 < t^o_2 < ...; t^o_n is paired with gamma_(n+1), and gamma_1 has no partner. All times are in t units.
# Method. E(t) = Re K_nu(2 pi) is evaluated in ball arithmetic (python-flint acb.bessel_k) at --prec bits.
# (1) Certified zeros. The interval (0, T] is split at the points j/20 (balls of width below 1e-37). At every point the ball E must exclude 0,
# and every cell on which the sign changes contains a zero (continuity): these are certified existence statements, so the number of sign
# change cells is a certified lower bound for the number of zeros. That a sign change cell holds exactly one zero and the other cells none is
# not certified by ball arithmetic here (the ball enclosures of E over a cell, or on a disc of the complex plane for a Cauchy estimate of E',
# are too crude: they are as wide as ten times the amplitude of E). It follows from Proposition pol:prop:skeleton(i) (Xi^o = |K| cos psi with
# psi strictly increasing) once the increments of psi over a cell are below pi, which step (3) checks on the grid (computed; an extra
# increment of 2 pi in one cell would need psi' > 125, and the observed maximum increment per cell is logged: about 0.1, that is psi' about 2 at the start).
# (2) Each zero is refined by bisection on the sign of the ball E to width 2e-11 (the midpoint is t^o_n, correct to 1e-10).
# (3) The continuous argument psi(t) of K_nu(2 pi) (Xi^o = |K| cos psi) from the principal arguments at the grid points (increments in
# (-pi, pi], logged between which values), and floor(psi(T)/pi + 1/2), compared with the number of sign change cells (consistency, computed).
# (4) theta, theta' in ball arithmetic at the zeros t^o_n and at the gamma_(n+1) of data/zeros_1700.txt (68 digits, read as 60 digits); N1.3,
# N1.4, N1.5 are computed from these (statistics in mpmath at 40 digits from the ball midpoints) and are COMPUTED, not certified, because the
# gamma_(n+1) come from the file and the exactness of the count is not certified (step 1).
# (5) Second run at --prec2 bits (default 1.5 times --prec): the maximum difference of the t^o_n is logged.
# Output: data/gp_skeleton.json (numbers as strings) and the log data/gp_skeleton.log next to it (or --out FILE).
# Usage: python3 gp_skeleton.py [zeros file] [--T 102] [--quick] [--out FILE]; --quick is T = 30 (about 1 second); the full case runs in about 10
# seconds on the laptop with the python-flint venv. On the cluster: bash misc/tools/hopper_run.sh -c 1 -t 10 -g gp_skeleton.json -g gp_skeleton.log
# code/gp_skeleton.py data/zeros_1700.txt -- --out gp_skeleton.json
import argparse
import json
import os
import sys
import time

import mpmath as mp
from flint import acb, arb, ctx

LOG = []


def log(s):
    print(s)
    sys.stdout.flush()
    LOG.append(s)


def Kball(t):
    """K_{1/4+it/2}(2 pi) for t an arb ball."""
    nu = acb(arb(1) / 4) + acb(0, 1) * acb(t) / 2
    return acb(2 * arb.pi()).bessel_k(nu)


def Eball(t):
    """Xi^o(t) = Re K_{1/4+it/2}(2 pi), t a real ball."""
    return Kball(t).real


def sign_of(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def bisect_zero(a, b, width=arb(2) / 10 ** 11):
    """Bisection on the sign of the ball E; a, b arb points with opposite signs."""
    sa = sign_of(Eball(a))
    assert sa != 0
    while (b - a) > width:
        m = arb(((a + b) / 2).mid())
        s = sign_of(Eball(m))
        assert s != 0, "ball contains 0 at bisection midpoint"
        if s == sa:
            a = m
        else:
            b = m
    return arb(((a + b) / 2).mid())


def zeros_cells(T, prec, nper=20):
    ctx.prec = prec
    ncell = int(round(T * nper))
    sg = []
    for i in range(ncell + 1):
        e = Eball(arb(i) / nper)
        s = sign_of(e)
        if s == 0:
            raise RuntimeError("E ball contains 0 at grid point %d/20: shift the grid" % i)
        sg.append(s)
    zs = []
    for i in range(ncell):
        if sg[i] != sg[i + 1]:
            zs.append(bisect_zero(arb(i) / nper, arb(i + 1) / nper))
    return zs


def psi_count(T, prec, nper=20):
    """Continuous argument of K at the grid points; floor(psi(T)/pi + 1/2), the largest and smallest increment."""
    ctx.prec = prec
    ncell = int(round(T * nper))
    mp.mp.dps = 30
    psi = mp.mpf(0)
    maxinc, mininc = mp.mpf(-10), mp.mpf(10)
    k0 = Kball(arb(0))
    assert k0.real > 0 and abs(k0.imag) < arb(10) ** -20
    prevarg = mp.mpf(0)
    for i in range(1, ncell + 1):
        k = Kball(arb(i) / nper)
        arg = mp.atan2(mp.mpf(k.imag.mid().str(30, radius=False)), mp.mpf(k.real.mid().str(30, radius=False)))
        inc = arg - prevarg
        while inc <= -mp.pi:
            inc += 2 * mp.pi
        while inc > mp.pi:
            inc -= 2 * mp.pi
        maxinc, mininc = max(maxinc, inc), min(mininc, inc)
        psi += inc
        prevarg = arg
    return int(mp.floor(psi / mp.pi + mp.mpf(1) / 2)), psi, maxinc, mininc


def theta_ball(t):
    z = acb(arb(1) / 4) + acb(0, 1) * acb(t) / 2
    return z.lgamma().imag - (arb(t) / 2) * arb.pi().log()


def theta_prime_ball(t):
    z = acb(arb(1) / 4) + acb(0, 1) * acb(t) / 2
    return z.digamma().real / 2 - arb.pi().log() / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("zfile", nargs="?", default=None)
    ap.add_argument("--T", type=float, default=102.0)
    ap.add_argument("--Nz", type=float, default=100.0)
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--prec2", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    zfile = a.zfile or os.path.join(root, "data", "zeros_1700.txt")
    T = 30.0 if a.quick else a.T
    Nz = 30.0 if a.quick else a.Nz
    out = a.out or os.path.join(root, "data", "gp_skeleton_quick.json" if a.quick else "gp_skeleton.json")
    prec2 = a.prec2 or int(1.5 * a.prec)
    res = {"T": T, "Nz": Nz, "prec": a.prec, "prec2": prec2}
    t0 = time.time()
    zeros = zeros_cells(T, a.prec)
    n_all = len(zeros)
    log("T = %g, %d bits: %d sign change cells of width 1/20 on (0, %g], each with a zero located by bisection to 2e-11 (time %.1f s)" % (T, a.prec, n_all, T, time.time() - t0))
    zeros2 = zeros_cells(T, prec2)
    if len(zeros2) == n_all:
        dmax = max(abs(x - y) for x, y in zip(zeros, zeros2))
        log("second run at %d bits: same count %d, maximum difference of the zeros %s" % (prec2, n_all, dmax.str(3)))
        res["max_diff_precisions"] = dmax.str(3)
    else:
        log("SECOND RUN DISAGREES in the count: %d against %d" % (n_all, len(zeros2)))
    ctx.prec = a.prec
    kpsi, psi, mx, mn = psi_count(T, a.prec)
    log("continuous argument: psi(%g) = %s, floor(psi/pi + 1/2) = %d against the count %d; increments per cell between %s and %s (all in (0, 1/2): %s, so psi' < %.2f on the grid)" % (T, mp.nstr(psi, 12), kpsi, n_all, mp.nstr(mn, 4), mp.nstr(mx, 4), bool(mn > 0 and mx < 0.5), float(mx * 20)))
    res["psi_T"] = mp.nstr(psi, 15)
    res["psi_count"] = kpsi
    mp.mp.dps = 40
    tvals = [mp.mpf(z.mid().str(35, radius=False)) for z in zeros]
    nN = sum(1 for v in tvals if v <= Nz)
    log("N1.1: zeros of Xi^o in (0, %g]: %d (in (0, %g]: %d)" % (Nz, nN, T, n_all))
    res["N1_1"] = nN
    log("N1.2: t^o_1 = %s (t^o_2 = %s, t^o_3 = %s, t^o_%d = %s)" % (mp.nstr(tvals[0], 12), mp.nstr(tvals[1], 12), mp.nstr(tvals[2], 12), n_all, mp.nstr(tvals[-1], 12)))
    res["N1_2"] = mp.nstr(tvals[0], 12)
    res["zeros"] = [mp.nstr(v, 14) for v in tvals]
    zs = []
    with open(zfile) as f:
        for line in f:
            if line.strip():
                zs.append(line.strip())
            if len(zs) >= n_all + 2:
                break
    gam = [mp.mpf(s[:60]) for s in zs]
    r = [mp.mpf(theta_ball(z).mid().str(35, radius=False)) / mp.pi - (n - mp.mpf(1) / 2) for n, z in enumerate(zeros, 1)]
    lo, hi = min(r[:nN]), max(r[:nN])
    log("N1.3: theta(t^o_n)/pi - (n - 1/2) for n <= %d: min %s (n = %d), max %s (n = %d); n = 1: %s, n = 2: %s, n = %d: %s" % (nN, mp.nstr(lo, 6), 1 + r[:nN].index(lo), mp.nstr(hi, 6), 1 + r[:nN].index(hi), mp.nstr(r[0], 6), mp.nstr(r[1], 6), n_all, mp.nstr(r[n_all - 1], 6)))
    res["N1_3"] = [mp.nstr(lo, 8), mp.nstr(hi, 8)]
    res["theta_resid"] = [mp.nstr(v, 8) for v in r]
    d, xs, gp = [], [], []
    for n in range(1, n_all + 1):
        g = gam[n]  # gamma_(n+1)
        gs = mp.nstr(g, 40)
        th = mp.mpf(theta_ball(arb(gs)).mid().str(35, radius=False))
        S = n - mp.mpf(1) / 2 - th / mp.pi
        thp = mp.mpf(theta_prime_ball(arb(gs)).mid().str(35, radius=False))
        d.append(tvals[n - 1] - g)
        xs.append(mp.pi * S / thp)
        gp.append(S)

    def stats(m):
        dd, xx = d[:m], xs[:m]
        mad = sum(abs(v) for v in dd) / m
        md, mxx = sum(dd) / m, sum(xx) / m
        cov = sum((u - md) * (v - mxx) for u, v in zip(dd, xx))
        vd = sum((u - md) ** 2 for u in dd)
        vx = sum((v - mxx) ** 2 for v in xx)
        return mad, cov / mp.sqrt(vd * vx), md
    for m in sorted({nN, n_all}):
        mad, cor, md = stats(m)
        log("n <= %d: mean |t^o_n - gamma_(n+1)| = %s, mean signed = %s, Pearson correlation with pi S/theta' = %s" % (m, mp.nstr(mad, 10), mp.nstr(md, 10), mp.nstr(cor, 12)))
        res["N1_45_n%d" % m] = [mp.nstr(mad, 12), mp.nstr(cor, 14), mp.nstr(md, 12)]
    res["N1_4"], res["N1_5"] = res["N1_45_n%d" % nN][0], res["N1_45_n%d" % nN][1]
    res["disp"] = [mp.nstr(v, 10) for v in d]
    res["Spred"] = [mp.nstr(v, 10) for v in xs]
    res["S_at_gamma"] = [mp.nstr(v, 8) for v in gp]
    log("displacement t^o_n - gamma_(n+1): smallest %s, largest %s (n <= %d); range of S(gamma_(n+1)): %s to %s" % (mp.nstr(min(d[:nN]), 6), mp.nstr(max(d[:nN]), 6), nN, mp.nstr(min(gp[:nN]), 5), mp.nstr(max(gp[:nN]), 5)))
    json.dump(res, open(out, "w"), indent=1)
    open(out.replace(".json", ".log"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
