#!/usr/bin/env python3
# cx_min.py: ball arithmetic certificates of the global minimum of the Thorin density u_theta (Chapter "Complexity", thm:globalmin,
# thm:uniform, thm:lbl, prop:slope and lem:pairing) with exact dyadic cell endpoints and outward rounding; zero counting by HSW22.
# Added for the book in October 2026 (group CX2), written from the mathematics (the collaborator's scripts were not read, run or copied).
#
# Conventions. xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, L = xi'/xi, u_theta(v) = Re L(theta + i v)/pi (Thorin density of U_theta),
# P_a(x) = a/(pi (a^2 + x^2)). Ordinates: the first 1000 zeros of data/zeros_1700.txt (68 digits). Each one is certified by a sign change
# of the real function Xi(t) = xi(1/2 + i t) at gamma -/+ 1e-25 (Arb, 160 bits) and enters every sum as the exact enclosure [gamma - 1e-25,
# gamma + 1e-25] (mode "tight"), or, to reproduce the chapter, as the grid interval of width 0.02 that contains it (mode "ch"). That the
# 1000 ordinates are all the zeros below 1419.9 (N(1419.9) = 1000, all on the line) is thm:certified, certified by group CX1, not here.
# Parts (--part). ords: certification, sums of gamma^-2, far sums by partial summation against HSW22, window counts, closed form
# constants. W: W(v) = pi^-1 sum_rho 1/(1+(v-gamma)^2) >= w* on [13.5, 1418.9] and its least certified bound (thm:lbl). S: S(v) > c*
# on the cells of the two cases of thm:uniform. B: inf of B_c over [13.5, 1418.9] for c = 1, 1/2, 1/4, 0 and V_c(0) (prop:slope).
# U1, U78: u_theta(v) > u_theta(0) on [13.5, 1420] by lem:cell (Taylor with the sup of |L''| over the cell as a ball) at theta = 1, 7/8.
# Cells. Dyadic endpoints (13.5 + k 3/16 for W, S, B; 13.5 + k/8 for U), bisected where the bound fails; a cell bound is evaluated with
# every term at its worst point of the cell and of the ordinate enclosure, so it is a lower bound for every point of the cell.
# The least certified bound is refined by bisection of the cells whose bound is below (the least midpoint value - tol) down to radius
# 2^-14 (2^-13 for U); the least midpoint value (an upper bound for the infimum) is logged, so the least bound is sandwiched.
# Error control. Every number is a ball; the leaf cell with the least bound is re-evaluated at 256 bits and the difference logged.
# Options: --quick (a short range and 40 ordinates, seconds), --out FILE, --procs N.
import os, sys, json, time, argparse
import math
import multiprocessing
from multiprocessing import Pool
from flint import arb, acb, acb_series, ctx

PREC = 128
ctx.prec = PREC
ctx.cap = 6
CFG = {}
LO = []
HI = []
T0 = time.time()


def log(*a):
    print("[%7.1f s]" % (time.time() - T0), *a, flush=True)


def cst(s):
    return arb(s).upper()


def fl(x):
    return float(x.mid())


def upf(x):
    f = float(x.mid())
    if arb(f) < x:
        f = math.nextafter(f, float('inf'))
    return f


# ---------------------------------------------------------------- ordinates
def zeros_file():
    for p in ["zeros_1700.txt", "data/zeros_1700.txt", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "zeros_1700.txt")]:
        if os.path.exists(p):
            return p
    raise SystemExit("zeros_1700.txt not found")


def load_zeros(n):
    ctx.prec = 320
    xs = [arb(t) for t in open(zeros_file()).read().split()[:n]]
    ctx.prec = PREC
    return xs


def xi(s):
    return (s * (s - 1) / 2) * ((-(s / 2) * arb.pi().log()).exp()) * (s / 2).gamma() * s.zeta()


def certify(gs, delta="1e-25"):
    ctx.prec = 160
    d = arb(delta)
    half = arb(1) / 2
    good = 0
    enc = []
    for g in gs:
        r0 = xi(acb(half, g - d)).real
        r1 = xi(acb(half, g + d)).real
        if (r0 > 0 and r1 < 0) or (r0 < 0 and r1 > 0):
            good += 1
        enc.append(((g - d).lower(), (g + d).upper()))
    ctx.prec = PREC
    return good, enc


def chapter_enc(gs):
    out = []
    for g in gs:
        k = int(float((g * 50).floor().mid()))  # g in [k/50, (k+1)/50)
        lo = (arb(k) / 50).lower()
        hi = (arb(k + 1) / 50).upper()
        assert lo <= g and g <= hi
        out.append((lo, hi))
    return out


# ---------------------------------------------------------------- HSW22 (the book's version, ch/basepoint.tex)
X0 = arb(30610046000)


def Mmain(T):
    return T / (2 * arb.pi()) * (T / (2 * arb.pi() * arb(1).exp())).log() + arb(7) / 8


def E_small(T):  # |N(T) - Mmain(T)| <= 2.5167 + 1/(50T) for e <= T <= X0
    return cst("2.5167") + 1 / (50 * T)


def E_gen(T):  # for every T >= e
    return cst("0.1095") * T.log() + cst("0.2042") * T.log().log() + cst("3.0305") + 1 / (50 * T)


def E_tru(T):  # Trudgian 2014 (incomplete proof), for comparison only
    return cst("0.112") * T.log() + cst("0.278") * T.log().log() + cst("2.51") + cst("0.2") / T


def E_task(T):  # the constants quoted in the task statement, for comparison only
    return cst("0.1038") * T.log() + cst("0.2573") * T.log().log() + cst("9.3675")


def count_lb(v, E, H=8):
    return (Mmain(v + H) - Mmain(v - H) - E(v + H) - E(v - H)).lower()


def far_pos(Tz, N0):
    # sum_{gamma > Tz} gamma^-2 <= int_{Tz}^inf 2 t^-3 (N(t) - N0) dt with N(t) <= Mmain + E, E = E_small up to X0, E_gen beyond
    k = 1 + (2 * arb.pi()).log()
    pi = arb.pi()

    def G(t, c0):
        return -(t.log() - k + 1) / (t * pi) - c0 / t ** 2 - 1 / (75 * t ** 3)

    c0 = arb(7) / 8 + cst("2.5167") - N0
    F1 = G(X0, c0) - G(Tz, c0)
    c1 = arb(7) / 8 - N0
    Mt = (X0.log() - k + 1) / (X0 * pi) + c1 / X0 ** 2
    a = cst("0.3137")
    b = cst("3.0305")
    Et = (a * X0.log() + b) / X0 ** 2 + a / (2 * X0 ** 2) + 1 / (75 * X0 ** 3)
    return (F1 + Mt + Et).upper()


def far_trudgian(Tz, N0):
    import mpmath as mp
    mp.mp.dps = 30
    Tz = mp.mpf(Tz)
    f = lambda t: 2 * t ** -3 * (t / (2 * mp.pi) * mp.log(t / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8 + 0.112 * mp.log(t) + 0.278 * mp.log(mp.log(t)) + 2.51 + 0.2 / t - N0)
    return mp.quad(f, [Tz, 10 * Tz, 1e3 * Tz, 1e6 * Tz, 1e9 * Tz, 1e12 * Tz, mp.inf])


# ---------------------------------------------------------------- cell machinery (W, S, B)
def far_dists(a, b):
    dp = []
    dn = []
    for lo, hi in zip(LO, HI):
        d1 = (b - lo).upper()
        d2 = (hi - a).upper()
        dp.append(d1 if d1 > d2 else d2)
        dn.append((b + hi).upper())
    return dp, dn


def lb_W(a, b):
    dp, dn = far_dists(a, b)
    s = arb(0)
    for d in (dp if CFG.get("pos") else dp + dn):
        s += 1 / (1 + d * d)
    return (s / arb.pi()).lower()


def lb_B(a, b):
    c2 = arb(CFG["c"]) ** 2 / 4
    dp, dn = far_dists(a, b)
    s = arb(0)
    for d in dp + dn:
        s += 1 / (c2 + d * d)
    return (s - CFG["V0up"] - CFG["FAR2"]).lower()


def mid_W(m):
    s = arb(0)
    for lo, hi in zip(LO, HI):
        g = lo.union(hi)
        s += 1 / (1 + (m - g) ** 2)
        if not CFG.get("pos"):
            s += 1 / (1 + (m + g) ** 2)
    return (s / arb.pi()).upper()


def mid_B(m):
    c2 = arb(CFG["c"]) ** 2 / 4
    s = arb(0)
    for lo, hi in zip(LO, HI):
        g = lo.union(hi)
        s += 1 / (c2 + (m - g) ** 2) + 1 / (c2 + (m + g) ** 2)
    return (s - CFG["V0lo"] - CFG["FAR2lo"]).upper()


def refine_cell(task):
    a, b, target, rmin = task
    target = arb(target)
    lbf = lb_W if CFG["kind"] == "W" else lb_B
    out = []
    stack = [(a, b)]
    while stack:
        a, b = stack.pop()
        lb = lbf(arb(a), arb(b))
        r = (b - a) / 2
        if lb >= target or r <= rmin:
            out.append((a, b, fl(lb), bool(lb >= target)))
        else:
            m = (a + b) / 2
            stack.append((a, m))
            stack.append((m, b))
    return out


def mid_task(m):
    return (m, fl((mid_W if CFG["kind"] == "W" else mid_B)(arb(m))))


def run_min(pool, cells, target1, tol, rmin, label):
    kind = CFG["kind"]
    leaves = []
    for res in pool.imap_unordered(refine_cell, [(a, b, upf(target1), 2.0 ** -20) for a, b in cells], chunksize=4):
        leaves += res
    nfail = sum(1 for x in leaves if not x[3])
    lbmin = min(x[2] for x in leaves)
    log(label, "pass 1: leaves", len(leaves), "failing", nfail, "least leaf bound", lbmin)
    near = [x for x in leaves if x[2] <= lbmin + 0.02]
    mids = pool.map(mid_task, [(x[0] + x[1]) / 2 for x in near], chunksize=2)
    mstar = min(mids, key=lambda t: t[1])
    log(label, "least midpoint value (upper bound of the infimum)", mstar[1], "at", mstar[0], "from", len(near), "midpoints")
    tgt = mstar[1] - tol
    todo = [(x[0], x[1]) for x in leaves if x[2] < mstar[1] - tol]
    leaves2 = []
    for res in pool.imap_unordered(refine_cell, [(a, b, tgt, rmin) for a, b in todo], chunksize=1):
        leaves2 += res
    rest = [x for x in leaves if x[2] >= mstar[1] - tol]
    allv = leaves2 + rest
    least = min(allv, key=lambda x: x[2])
    log(label, "pass 2: refined", len(todo), "cells into", len(leaves2), "leaves; certified least bound", least[2], "in cell", (least[0], least[1]))
    return dict(leaves1=len(leaves), failing1=nfail, least1=lbmin, mid_least_upper=mstar[1], mid_at=mstar[0], leaves2=len(leaves2),
                least=least[2], least_cell=[least[0], least[1]]), least


def recheck_lb(cell, kind):
    ctx.prec = 256
    try:
        v = (lb_W if kind == "W" else lb_B)(arb(cell[0]), arb(cell[1]))
        return fl(v)
    finally:
        ctx.prec = PREC


def make_cells(a0, a1, w):
    out = []
    a = a0
    while a < a1:
        out.append((a, a + w))
        a += w
    return out


def set_ords(enc):
    global LO, HI
    LO = [x[0] for x in enc]
    HI = [x[1] for x in enc]


def coarse_task(k):
    # the chapter's cell k: [13.5 + 0.2 k, 13.7 + 0.2 k] (radius 0.1), exact rational endpoints enclosed by balls
    a = (arb(135 + 2 * k) / 10).lower()
    b = (arb(137 + 2 * k) / 10).upper()
    CFG["pos"] = True
    p = lb_W(a, b)
    CFG["pos"] = False
    q = lb_W(a, b)
    return (k, fl(p), fl(q))


def coarse_grid_floats(k):
    return (math.nextafter((135 + 2 * k) / 10, -float("inf")), math.nextafter((137 + 2 * k) / 10, float("inf")))


# ---------------------------------------------------------------- S cells (thm:uniform)
def S_cell(task):
    a0, b0 = task
    out = []
    stack = [(a0, b0)]
    q = arb(1) / 4
    while stack:
        a, b = stack.pop()
        Aa = arb(a)
        Bb = arb(b)
        dp, dn = far_dists(Aa, Bb)
        done = False
        if any(d <= CFG["d105"] for d in dp):
            out.append((a, b, 1, 0.0))
            done = True
        else:
            near = None
            for lo, hi in zip(LO, HI):
                if hi < Aa:
                    n1 = (Aa - hi).lower()
                elif lo > Bb:
                    n1 = (lo - Bb).lower()
                else:
                    n1 = arb(0)
                n2 = (Aa + lo).lower()
                n = n1 if n1 < n2 else n2
                near = n if (near is None or n < near) else near
            if near >= CFG["d095"]:
                s = arb(0)
                for d in dp + dn:
                    d2 = d * d
                    s += (d2 - q) / (d2 + q) ** 2
                s = (s / arb.pi()).lower()
                if s > CFG["cstar_up"]:
                    out.append((a, b, 2, fl(s)))
                    done = True
        if not done:
            if (b - a) / 2 <= 2.0 ** -20:
                out.append((a, b, 0, 0.0))
            else:
                m = (a + b) / 2
                stack.append((a, m))
                stack.append((m, b))
    return out


# ---------------------------------------------------------------- u cells (thm:globalmin, lem:cell)
def Lseries(s):
    t = acb_series([s, 1])
    Z = t.zeta()
    G = (t / 2).lgamma().derivative()
    Ls = 1 / t + 1 / (t - 1) - arb.pi().log() / 2 + G + Z.derivative() / Z
    return Ls[0], Ls[1], 2 * Ls[2]


def u_cell_bound(m, r, theta, ball):
    pi = arb.pi()
    L, Lp, _ = Lseries(acb(theta, m))
    um = L.real / pi
    up = -Lp.imag / pi
    _, _, Lpp = Lseries(acb(theta, ball))
    lb = um - r * abs(up) - r * r / (2 * pi) * abs(Lpp).upper()
    return lb.lower(), um


def u_task(task):
    a, b, theta, target, rmin = task
    target = arb(target)
    th = arb(theta)
    out = []
    stack = [(a, b)]
    while stack:
        a, b = stack.pop()
        m = (a + b) / 2
        r = (b - a) / 2
        try:
            lb, um = u_cell_bound(arb(m), arb(r), th, arb(a).union(arb(b)))
            ok = bool(lb >= target)
            fin = lb.is_finite()
        except Exception:
            lb, um, ok, fin = arb(0), arb(0), False, False
        if (ok and fin) or r <= rmin:
            out.append((a, b, fl(lb) if fin else float("nan"), fl(um) if fin else float("nan"), ok and fin))
        else:
            stack.append((a, m))
            stack.append((m, b))
    return out


def run_u(pool, theta, a0, a1, w, tol, rmin, label, fixed=False):
    th = arb(theta)
    if theta == 1:
        pi = arb.pi()
        u0 = (1 + arb.const_euler() / 2 - (4 * pi).log() / 2) / pi
    else:
        L, _, _ = Lseries(acb(th, 0))
        u0 = L.real / arb.pi()
    u0up = u0.upper()
    log(label, "u_theta(0) =", u0.str(15), " upper end", fl(u0up))
    cells = make_cells(a0, a1, w)
    if fixed:  # the chapter's grid: radius 0.004, cells [13.5 + 0.008 k, +0.008] (floats, an emulation, not a certificate), no refinement
        cells = [(13.5 + 0.008 * k, 13.5 + 0.008 * (k + 1)) for k in range(int((a1 - 13.5) / 0.008 - 1e-9) + 1)]
    leaves = []
    for res in pool.imap_unordered(u_task, [(a, b, theta, upf(u0up), 2.0 ** -12 if not fixed else 2.0 ** -18) for a, b in cells], chunksize=2):
        leaves += res
    nbad = sum(1 for x in leaves if not x[4])
    lbmin = min(leaves, key=lambda x: x[2])
    mid_min = min(leaves, key=lambda x: x[3])
    log(label, "pass 1: initial cells", len(cells), "leaves", len(leaves), "failing", nbad, "least bound", lbmin[2], "at", (lbmin[0] + lbmin[1]) / 2,
        "least midpoint value", mid_min[3], "at", (mid_min[0] + mid_min[1]) / 2)
    if fixed:
        return dict(theta=theta, u0=u0.str(15), cells0=len(cells), leaves=len(leaves), failing=nbad, least=lbmin[2],
                    least_at=(lbmin[0] + lbmin[1]) / 2, least_cell=[lbmin[0], lbmin[1]])
    mstar = mid_min[3]
    tgt = mstar - tol
    todo = [(x[0], x[1]) for x in leaves if x[2] < mstar - tol]
    leaves2 = []
    for res in pool.imap_unordered(u_task, [(a, b, theta, tgt, rmin) for a, b in todo], chunksize=1):
        leaves2 += res
    rest = [x for x in leaves if not (x[2] < mstar - tol)]
    allv = leaves2 + rest
    least = min(allv, key=lambda x: x[2])
    log(label, "pass 2: refined", len(todo), "cells into", len(leaves2), "leaves; certified least bound", least[2], "at", (least[0] + least[1]) / 2)
    # second precision
    ctx.prec = 256
    try:
        lb2, _ = u_cell_bound(arb((least[0] + least[1]) / 2), arb((least[1] - least[0]) / 2), arb(theta), arb(least[0]).union(arb(least[1])))
    finally:
        ctx.prec = PREC
    log(label, "least cell re-evaluated at 256 bits:", fl(lb2), "difference", abs(fl(lb2) - least[2]))
    return dict(theta=theta, u0=u0.str(15), u0_upper=fl(u0up), cells0=len(cells), leaves1=len(leaves), failing1=nbad,
                least1=lbmin[2], least1_at=(lbmin[0] + lbmin[1]) / 2, mid_least=mid_min[3], mid_least_at=(mid_min[0] + mid_min[1]) / 2,
                leaves2=len(leaves2), least=least[2], least_at=(least[0] + least[1]) / 2, least_r=(least[1] - least[0]) / 2,
                least_256=fl(lb2), failing2=sum(1 for x in leaves2 if not x[4]),
                margin_over_u0=least[2] - fl(u0up))


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True, choices=["ords", "W", "S", "B", "U1", "U78", "U1c", "U78c"])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--procs", type=int, default=0)
    ap.add_argument("--vmax", type=float, default=0.0, help="stop the cell range here (partial run, all ordinates)")
    args = ap.parse_args()
    multiprocessing.set_start_method('fork')
    nproc = args.procs or (len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count())
    nord = 40 if args.quick else 1000
    vmax = args.vmax or (40.0 if args.quick else 1418.9)
    gs = load_zeros(nord)
    res = dict(part=args.part, quick=args.quick, n_ordinates=nord)
    good, enc = certify(gs)
    log("certified sign changes", good, "of", nord)
    res["certified_sign_changes"] = good
    res["gamma_last"] = gs[-1].str(20, radius=False)
    chenc = chapter_enc(gs)
    # constants from ords (needed by every part)
    Tz = arb("1419.9")
    N0 = arb(1000)
    Fp = far_pos(Tz, N0)
    Fp_nocount = far_pos(Tz, (Mmain(Tz) - E_small(Tz)).lower())
    FAR = 2 * Fp  # both signs
    sums = {}
    for mode, en in (("tight", enc), ("ch", chenc)):
        lo = [x[0] for x in en]
        sums[mode] = {c: (sum(2 / (arb(c) ** 2 / 4 + l * l) for l in lo)).upper() for c in (1, 0.5, 0.25, 0)}
    if args.part == "ords":
        res["far_sum_pos_HSW22_N1000"] = fl(Fp)
        res["far_sum_both_HSW22_N1000"] = fl(FAR)
        res["far_sum_both_HSW22_Ncount_lower"] = fl(2 * Fp_nocount)
        try:
            res["far_sum_both_Trudgian_N1000_mpmath"] = float(2 * far_trudgian(1419.9, 1000))
        except Exception as e:
            res["far_sum_both_Trudgian_N1000_mpmath"] = str(e)
        for mode in sums:
            res["sum_inv_gamma2_both_upper_" + mode] = fl(sums[mode][0])
            res["V_c0_upper_" + mode] = {str(c): fl(v) for c, v in sums[mode].items()}
            cs = (sums[mode][0] + FAR) / arb.pi()
            res["cstar_upper_" + mode] = fl(cs.upper())
        # window counts
        v = arb(14189) / 10
        res["count_HSW22_book_1418.9"] = fl(count_lb(v, E_small))
        res["count_HSW22_book_1420"] = fl(count_lb(arb(1420), E_small))
        res["count_HSW22_book_at_X0-8"] = fl(count_lb(X0 - 8, E_small))
        res["count_HSW22_gen_at_X0-8"] = fl(count_lb(X0 - 8, E_gen))
        res["count_Trudgian_1418.9"] = fl(count_lb(v, E_tru))
        res["count_Trudgian_1420"] = fl(count_lb(arb(1420), E_tru))
        res["count_task_constants_1418.9"] = fl(count_lb(v, E_task))
        res["count_HSW22_gen_1418.9"] = fl(count_lb(v, E_gen))
        # derivative of the lower bound is positive: D'(v) - E'(v+8) - E'(v-8) on [1418.9, X0-8] and beyond
        vs = [arb("1418.9"), arb(10) ** 4, arb(10) ** 6, arb(10) ** 8, arb(10) ** 10]
        res["monotone_check_small"] = [fl(((((v + 8) / (v - 8)).log() / (2 * arb.pi())) - (1 / 50) * (1 / (v - 8) ** 2 + 1 / (v + 8) ** 2))) for v in vs]
        res["monotone_check_gen"] = [fl(((((v + 8) / (v - 8)).log() / (2 * arb.pi())) - (cst("0.1095") + cst("0.2042") / (v - 8).log()) / (v - 8) - (cst("0.1095") + cst("0.2042") / (v + 8).log()) / (v + 8))) for v in [X0, arb(10) ** 12, arb(10) ** 20]]
        pi = arb.pi()
        P = lambda a, x: a / (pi * (a * a + x * x))
        u1_0 = (1 + arb.const_euler() / 2 - (4 * pi).log() / 2) / pi
        L78, _, _ = Lseries(acb(arb(7) / 8, 0))
        u78_0 = L78.real / pi
        res["u1_0"] = u1_0.str(15, radius=False)
        res["u78_0"] = u78_0.str(15, radius=False)
        res["P1_8"] = fl(P(arb(1), 8))
        res["P34_8"] = fl(P(arb(3) / 4, 8))
        res["P12_8"] = fl(P(arb(1) / 2, 8))
        for n in (7, 9):
            res["globalmin_theta1_n%d" % n] = fl(n * P(arb(1), 8) / 2)
            res["globalmin_theta78_n%d" % n] = fl(n * P(arb(3) / 4, 8) / 2)
            res["rhmin_n%d" % n] = fl(n * P(arb(1) / 2, 8))
        res["n_needed_theta1"] = fl(2 * u1_0 / P(arb(1), 8))
        res["n_needed_theta78"] = fl(2 * u78_0 / P(arb(3) / 4, 8))
        res["P38_105"] = fl(P(arb(3) / 8, arb(105) / 100))
        res["7_over_65pi"] = fl(7 / (65 * pi))
        res["9_over_65pi"] = fl(9 / (65 * pi))
        # prop:slope counting range: min over c in [0,1] of n/(c^2+64) - V_c(0) - 2*FAR (decided form: far sum not dropped)
        for mode in sums:
            for n in (7, 9):
                best = None
                for j in range(0, 65):
                    c = arb(j) / 64
                    V0 = (sum(2 / (c * c / 4 + l * l) for l in [x[0] for x in (enc if mode == "tight" else chenc)])).upper()
                    val = n / (c * c + 64) - V0 - 2 * FAR
                    best = val if best is None or val.lower() < best.lower() else best
                res["slope_far_range_%s_n%d_min_over_c_grid" % (mode, n)] = fl(best.lower())
        log(json.dumps(res, indent=1))
    else:
        modes = [("tight", enc), ("ch", chenc)]
        if args.part == "W":  # also the positive ordinates only (the chapter's 0.0617 is of this kind)
            modes += [("tight_pos", enc), ("ch_pos", chenc)]
        for mode, en in modes:
            set_ords(en)
            label = "%s/%s" % (args.part, mode)
            if args.part in ("U1", "U78", "U1c", "U78c"):
                if mode == "ch":
                    continue
                a1 = (args.vmax or 40.0) if (args.quick or args.vmax) else 1420.0
                with Pool(nproc) as pool:
                    res[mode] = run_u(pool, 1 if args.part.startswith("U1") else 0.875, 13.5, a1, 0.125, 1e-5, 2.0 ** -13, label, fixed=args.part.endswith("c"))
                continue
            cells = make_cells(13.5, vmax, 0.1875)
            if args.part == "W":
                CFG.update(kind="W", pos=mode.endswith("_pos"))
                cstar_ours = (sums[mode.split("_")[0]][0] + FAR) / arb.pi()
                with Pool(nproc) as pool:
                    d, least = run_min(pool, cells, cst("0.029445"), 1e-5, 2.0 ** -14, label)
                d["w_star_chapter"] = 0.029445
                d["w_star_ours_2cstar"] = fl(2 * cstar_ours.upper())
                d["recheck_256"] = recheck_lb(least, "W")
                d["bracket_mid_range_pi(least-wstar)"] = float(3.141592653589793 * (d["least"] - 0.029445))
                if mode in ("ch", "ch_pos"):
                    ncoarse = int((vmax - 13.5) / 0.2 - 1e-9) + 1
                    with Pool(nproc) as pool:
                        cr = pool.map(coarse_task, range(ncoarse), chunksize=8)
                    if mode == "ch":
                        d["coarse_cells_r0.1"] = ncoarse
                        d["coarse_least_pos_only"] = min(x[1] for x in cr)
                        d["coarse_least_pos_only_centre"] = 13.6 + 0.2 * min(cr, key=lambda x: x[1])[0]
                        d["coarse_least_both_signs"] = min(x[2] for x in cr)
                        d["coarse_least_both_signs_centre"] = 13.6 + 0.2 * min(cr, key=lambda x: x[2])[0]
                        d["coarse_pi_least_pos_minus_wstar"] = float(3.141592653589793 * (d["coarse_least_pos_only"] - 0.029445))
                        log(label, "chapter grid (radius 0.1, positive ordinates, ch enclosures):", ncoarse, "cells, least", d["coarse_least_pos_only"], "centre", d["coarse_least_pos_only_centre"], "; both signs", d["coarse_least_both_signs"])
                res[mode] = d
            elif args.part == "S":
                CFG.update(d105=(arb(21) / 20).lower(), d095=(arb(19) / 20).upper(), cstar_up=cst("0.0147225"))
                leaves = []
                with Pool(nproc) as pool:
                    for r in pool.imap_unordered(S_cell, cells, chunksize=4):
                        leaves += r
                k1 = sum(1 for x in leaves if x[2] == 1)
                k2 = [x for x in leaves if x[2] == 2]
                k0 = sum(1 for x in leaves if x[2] == 0)
                cs = ((sums[mode][0] + FAR) / arb.pi()).upper()
                res[mode] = dict(initial_cells=len(cells), leaves=len(leaves), first_kind=k1, second_kind=len(k2), failing=k0,
                                 least_S_second_kind=min(x[3] for x in k2) if k2 else None, cstar_chapter=0.0147225, cstar_ours=fl(cs),
                                 margin_over_cstar_chapter=(min(x[3] for x in k2) - 0.0147225) if k2 else None)
                log(label, json.dumps(res[mode]))
                ncoarse = int((vmax - 13.5) / 0.2 - 1e-9) + 1
                leaves = []
                with Pool(nproc) as pool:
                    for r in pool.imap_unordered(S_cell, [coarse_grid_floats(k) for k in range(ncoarse)], chunksize=4):
                        leaves += r
                k2 = [x for x in leaves if x[2] == 2]
                res[mode + "_r0.1"] = dict(initial_cells=ncoarse, leaves=len(leaves), first_kind=sum(1 for x in leaves if x[2] == 1), second_kind=len(k2),
                                           failing=sum(1 for x in leaves if x[2] == 0), least_S_second_kind=min(x[3] for x in k2) if k2 else None)
                log(label, "grid of radius 0.1:", json.dumps(res[mode + "_r0.1"]))
            elif args.part == "B":
                res[mode] = {}
                for c in (1, 0.5, 0.25, 0):
                    V0up = sums[mode][c]
                    V0lo = (sum(2 / (arb(c) ** 2 / 4 + h * h) for h in HI)).lower()
                    CFG.update(kind="B", c=c, V0up=V0up, V0lo=V0lo, FAR2=2 * FAR, FAR2lo=2 * FAR)
                    with Pool(nproc) as pool:
                        d, least = run_min(pool, cells, arb(0.17), 1e-6, 2.0 ** -14, label + "/c=%s" % c)
                    d["V_c0_upper"] = fl(V0up)
                    d["FAR2_used"] = fl(2 * FAR)
                    d["recheck_256"] = recheck_lb(least, "B")
                    res[mode][str(c)] = d
    outp = args.out or ("cx_c2_%s%s.json" % (args.part, "_quick" if args.quick else ""))
    json.dump(res, open(outp, "w"), indent=1)
    log("wrote", outp)


if __name__ == "__main__":
    main()
