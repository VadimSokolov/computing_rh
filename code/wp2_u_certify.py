# wp2_u_certify.py: work package 2 (the Thorin level) of misc/reviews/nick-runs/vadim_runs/rh_compute.tex.
# Extends code/cx_u_certify.py, part (i), to any level alpha in (1/2, 1) and any height range [t0, t1] (integers).
# Certified claim (Lemma cell of cx_u_certify.py, ball arithmetic): on every accepted cell {alpha + i w : w in [a, b]} the ball of zeta(s0)
# excludes 0 (so xi has no zero in the cell) and u_alpha(w) >= LB(cell) := u(m) - r |u'(m)| - r^2/(2 pi) sup|L''|, with u = Re(xi'/xi)/pi, sup|L''| <= 2|c_2|
# over the ball cell. Cells have exact dyadic endpoints [i 2^-j, (i+1) 2^-j], start at j = j0, are bisected until accepted (LB > 0), up to j = JMAX.
# A cell with LB > 0 certifies u_alpha > 0 on it, which is the absence of zeros of xi on that segment of the line Re s = alpha; all cells together
# certify u_alpha > 0 on [t0, t1] (u is even, so negative heights follow).
# Additional records, by height bin (--bins): the least cell bound (certified), the least midpoint value u(m) over accepted cells (computed in the same
# ball arithmetic but a midpoint, not a bound: an upper bound for the infimum of u on the bin), the cell count and maximal depth, and
# the largest right end of a cell with LB <= u_alpha(0) (so u_alpha(v) > u_alpha(0) is certified for all v above it up to t1).
# u_alpha(0) is exact in ball arithmetic. The constant factor of xi is irrelevant (only xi'/xi enters).
# Usage: python3 wp2_u_certify.py --alpha 3/4 --t0 0 --t1 100 --prec 96 --procs 8 --j0 7 --bins 0,10,100 --out wp2_u_a075_0_100.json
import argparse, json, math, time
from fractions import Fraction
from multiprocessing import Pool
from flint import arb, acb, acb_series, ctx

JMAX = 30
U = None
LEMMA = "c2"


def setprec(p):
    global U
    ctx.prec = p
    U = arb(0, 1)


def L_coeffs(s0, n):
    x = acb_series([s0, 1], prec=n + 1)
    Z = x.zeta()
    if not (Z.coeffs()[0].abs_lower() > 0):
        return None
    g = (x / 2).lgamma()
    L = x.inv() + (x - 1).inv() + g.derivative() + Z.derivative() / Z - arb.pi().log() / 2
    return L.coeffs()


def cell_lb(theta, i, j):
    mid = arb(2 * i + 1) / (2 ** (j + 1))
    r = arb(1) / (2 ** (j + 1))
    if LEMMA == "c2":
        cb = L_coeffs(acb(theta, mid + r * U), 3)   # the ball first: a cell whose ball of zeta contains 0 is bisected after one evaluation
        if cb is None:
            return None, None
        cm = L_coeffs(acb(theta, mid), 2)
        if cm is None:
            return None, None
        sup2 = 2 * cb[2].abs_upper()
    else:
        # Lemma cell with a sharper bound of sup |L''| (this file, not cx_u_certify.py): with c_k(m) the Taylor coefficients of L at the exact midpoint
        # and c_4(ball) the fourth coefficient of the series at the ball cell (it encloses c_4(xi) for every xi in the cell),
        # |L''(m+x)| <= 2|c_2(m)| + 6|c_3(m)| |x| + 12 |c_4(xi)| |x|^2 (Taylor of L'' to first order with Lagrange remainder), |x| <= r.
        cm = L_coeffs(acb(theta, mid), 4)
        if cm is None:
            return None, None
        cb = L_coeffs(acb(theta, mid + r * U), 5)
        if cb is None:
            return None, None
        sup2 = 2 * cm[2].abs_upper() + 6 * cm[3].abs_upper() * r + 12 * cb[4].abs_upper() * r * r
    pi = arb.pi()
    u0 = cm[0].real / pi
    u1 = -cm[1].imag / pi
    return u0 - r * abs(u1) - r * r / (2 * pi) * sup2, u0


def process(theta, stack, leaves, fails, st):
    """Bisect the cells on the stack until each is accepted (LB > 0, ball of zeta excludes 0) or reaches JMAX; accepted cells go to leaves."""
    while stack:
        i, j = stack.pop()
        lb, um = cell_lb(theta, i, j)
        if lb is not None and lb > 0:
            leaves.append((float(lb.lower().mid()), i, j, float(um.mid())))
        else:
            if j >= JMAX:
                fails.append((i, j))
                continue
            st[0] += 1
            stack.append((2 * i + 1, j + 1))
            stack.append((2 * i, j + 1))


def cell_task(a):
    an, ad, i0, i1, j0, prec, thr, lem, q = a
    global LEMMA
    LEMMA = lem
    setprec(prec)
    theta = arb(an) / ad
    thr = float(thr)
    leaves = []
    fails = []
    st = [0]
    process(theta, [(i, j0) for i in range(i1 - 1, i0 - 1, -1)], leaves, fails, st)
    # refinement of the least bounds: bisect the leaf with the least LB until it is within the fraction q of the least midpoint value of the task,
    # so that the least certified bound is a tight lower bound of the infimum of u (a leaf is a cell with LB > 0; bisecting keeps all bounds valid)
    if leaves and q > 0:
        for _ in range(400):
            target = q * min(l[3] for l in leaves)
            k = min(range(len(leaves)), key=lambda t: leaves[t][0])
            if leaves[k][0] >= target or leaves[k][2] >= JMAX:
                break
            _, i, j, _ = leaves.pop(k)
            process(theta, [(2 * i + 1, j + 1), (2 * i, j + 1)], leaves, fails, st)
    n = len(leaves)
    least = leastmid = endthr = None
    maxj = j0
    for v, i, j, vm in leaves:
        h = (i + 0.5) / 2 ** j
        if least is None or v < least[0]:
            least = (v, h)
        if leastmid is None or vm < leastmid[0]:
            leastmid = (vm, h)
        if v <= thr:
            e = (i + 1) / 2 ** j
            if endthr is None or e > endthr:
                endthr = e
        maxj = max(maxj, j)
    return n, least, leastmid, maxj, fails, st[0], endthr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="3/4")
    ap.add_argument("--t0", type=int, default=0)
    ap.add_argument("--t1", type=int, default=100)
    ap.add_argument("--prec", type=int, default=96)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--j0", type=int, default=7)
    ap.add_argument("--lemma", default="c2", choices=["c2", "c4"])
    ap.add_argument("--q", type=float, default=0.98, help="refine the least bounds to this fraction of the least midpoint value (0 = no refinement)")
    ap.add_argument("--chunk", type=int, default=64)
    ap.add_argument("--bins", default="0,10,100,1000,10000,100000")
    ap.add_argument("--out", required=True)
    A = ap.parse_args()
    fr = Fraction(A.alpha)
    an, ad = fr.numerator, fr.denominator
    setprec(A.prec)
    logf = open(A.out.replace(".json", ".log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    edges = [int(x) for x in A.bins.split(",")]
    edges = [e for e in edges if A.t0 < e < A.t1]
    edges = [A.t0] + edges + [A.t1]
    c0 = L_coeffs(acb(arb(an) / ad, 0), 1)
    u0 = c0[0].real / arb.pi()
    thr = float(u0.upper())
    log("wp2_u_certify: alpha", A.alpha, "range [%d, %d]" % (A.t0, A.t1), "prec", A.prec, "j0", A.j0, "procs", A.procs, "u(0) =", u0.str(15, radius=True),
        "pi u(0) = b_alpha =", (u0 * arb.pi()).str(15, radius=True))
    pool = Pool(A.procs)
    t0 = time.time()
    res = {"alpha": A.alpha, "prec": A.prec, "j0": A.j0, "lemma": A.lemma, "q": A.q, "u0": u0.str(15, radius=True), "u0_float": thr, "bins": []}
    for lo, hi in zip(edges, edges[1:]):
        c_lo, c_hi = lo << A.j0, hi << A.j0
        tasks = [(an, ad, i, min(i + A.chunk, c_hi), A.j0, A.prec, thr, A.lemma, A.q) for i in range(c_lo, c_hi, A.chunk)]
        ncells = 0
        nb = 0
        maxj = 0
        fails = []
        least = None
        leastmid = None
        endthr = None
        for n, ls, lm, mj, f, b, e in pool.imap_unordered(cell_task, tasks):
            ncells += n
            nb += b
            maxj = max(maxj, mj)
            fails += f
            if ls is not None and (least is None or ls[0] < least[0]):
                least = ls
            if lm is not None and (leastmid is None or lm[0] < leastmid[0]):
                leastmid = lm
            if e is not None and (endthr is None or e > endthr):
                endthr = e
        d = {"lo": lo, "hi": hi, "cells": ncells, "bisections": nb, "maxj": maxj, "fails": len(fails), "fail_list": [(int(i), int(j)) for i, j in fails[:20]],
             "least_lb": least[0] if least else None, "least_lb_at": least[1] if least else None,
             "least_mid": leastmid[0] if leastmid else None, "least_mid_at": leastmid[1] if leastmid else None,
             "last_end_lb_le_u0": endthr, "seconds": time.time() - t0}
        res["bins"].append(d)
        log("bin [%d,%d]: cells %d bisections %d maxj %d fails %d | least LB %.6e at v=%.5f | least midpoint u %.6e at v=%.5f | last end of cell with LB<=u(0): %s | t=%.0fs"
            % (lo, hi, ncells, nb, maxj, len(fails), least[0] if least else float("nan"), least[1] if least else float("nan"),
               leastmid[0] if leastmid else float("nan"), leastmid[1] if leastmid else float("nan"), endthr, time.time() - t0))
        json.dump(res, open(A.out, "w"), indent=1)
    log("done")


if __name__ == "__main__":
    main()
