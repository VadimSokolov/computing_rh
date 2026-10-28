# cx_u_certify.py: certificate of Theorem thm:certified of the chapter "complexity" (Nick Polson's notes, source2/lvjbook/chapter_complexity.tex),
# redone in Arb from the mathematics (the collaborator's scripts were not read, run or copied). Task C1 of misc/reviews/nick-b4/reconcile-I.md.
# Claims reproduced. (i) u_theta(v) = pi^{-1} Re (xi'/xi)(theta + i v) > 0 on [0, 1420] for theta = 3/4 and 7/8, with the least cell lower bounds
# (chapter 1.1e-5 and 2.7e-4) and u_theta(0); (ii) N(T*) = 1000 at T* = 14199/10 exactly; (iii) 2 N(3/4, T*) = 0 and M_{3/4}(T*) (chapter
# 1000.0195 +- 2.6e-5); (iv) 1000 sign changes of Xi(v) = xi(1/2 + i v) at v = k/50, k = 0..70995.
# Conventions. xi(s) = s (s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s); only xi'/xi and signs enter, so the constant factor is irrelevant. L = xi'/xi.
# Method (i), Lemma cell: for a cell {theta + i w : |w - m| <= r} on which xi has no zero, u(m+t) >= u(m) - r |u'(m)| - r^2/(2 pi) sup |L''|, with
# u'(m) = -pi^{-1} Im L'(m). Cell endpoints are exact dyadics (cells [i 2^-j, (i+1) 2^-j], initial j = 7, i < 181760, so the last cell ends at
# exactly 1420; bisection to j <= 26). L and L' are evaluated at the exact midpoint, and sup |L''| <= 2 |c_2| over the whole cell, c_2 the x^2 coefficient
# of the acb_series expansion of L at the ball s0 = theta + i [m-r, m+r] (the ball is exactly the cell), so c_2 encloses L''(s)/2 for every s in the cell.
# L = 1/s + 1/(s-1) - (log pi)/2 + d/dx lgamma((s0+x)/2) + zeta'/zeta, all as series in x. A cell is accepted only if the lower bound is certainly
# positive (which also needs the ball of zeta(s0) to exclude 0, so the cell has no zero of xi); otherwise it is bisected. u is even, so [0,1420] suffices.
# Method (ii), (iii), Lemma segment: along a segment [a,b] of length d, if d sup |xi'| < |xi(a)| (sup over a ball covering the segment, xi' from the
# acb_series of xi at that ball), the argument change is Arg(xi(b)/xi(a)) in (-pi/2, pi/2); otherwise the segment is bisected. All endpoints are exact
# dyadic except the corner T* = 14199/10, which is a ball containing it (every path through the corner uses the same ball object, so the
# contours close; all values are enclosures of the exact values at the exact corner). pi N(T*) = Delta arg along 2 -> 2 + iT* -> 1/2 + iT* (right
# line R, top 2 -> 3/4 TopA, top 3/4 -> 1/2 TopB; xi real positive on the real axis s > 1/2 and real on Re s = 1/2, by symmetry the full
# contour gives 4 times this; standard Backlund). 2 pi N(3/4,T) = R + TopA - Left (Left = Delta arg on Re s = 3/4 from 0 to T), the bottom side
# [3/4, 2] contributing 0 because xi > 0 on (1/2, infinity) (standard fact, not certified here; zeta < 0 on (0,1), > 0 on (1, infinity)). M_{3/4}(T) = Left / pi.
# Method (iv): Xi(k/50) at the ball arb(k)/50 (xi(1/2 + i v) is real, so the real part of the enclosure holds its sign); sign changes counted.
# Error control: the whole run is made at two precisions (--prec 96 and 160); the least bounds, the enclosures and the counts are compared in the checks file.
# Usage: python3 cx_u_certify.py [--quick] [--prec P] [--procs N] [--out FILE]; --quick uses T* = 81/2 (6 zeros) and finishes in a minute on Hopper (writes to
# /tmp unless --out is given); the full run (T* = 14199/10) is a Hopper job: bash misc/tools/hopper_run.sh -c 24 -t 180 -g cx_u_certify_96.json
# -g cx_u_certify_96.log code/cx_u_certify.py -- --prec 96 --procs 24 --out cx_u_certify_96.json
import argparse, json, sys, time, math
from fractions import Fraction
from multiprocessing import Pool
from flint import arb, acb, acb_series, ctx

PREC = 96
JMAX = 26
U = None  # arb(0,1), the unit radius ball


def setprec(p):
    global U
    ctx.prec = p
    U = arb(0, 1)


def pack(x):
    m, e = x.mid().man_exp()
    rm, re_ = x.rad().man_exp()
    return (int(m), int(e), int(rm), int(re_))


def unpack(t):
    m, e, rm, re_ = t
    old = ctx.prec
    ctx.prec = max(old, 2 * (abs(m).bit_length() + 8))
    mid = arb(m) * arb(2) ** e
    rad = arb(rm) * arb(2) ** re_
    x = mid + rad * arb(0, 1)
    ctx.prec = old
    return x


def L_coeffs(s0, n):
    """First n coefficients of L(s0 + x) = (xi'/xi)(s0 + x) as balls; None if the ball of zeta(s0) contains 0."""
    x = acb_series([s0, 1], prec=n + 1)
    Z = x.zeta()
    if not (Z.coeffs()[0].abs_lower() > 0):
        return None
    g = (x / 2).lgamma()
    L = x.inv() + (x - 1).inv() + g.derivative() + Z.derivative() / Z - arb.pi().log() / 2
    return L.coeffs()


def xi_val(s):
    """xi(s) (one ball) and, if wanted, xi'(s), via the series of xi at s."""
    return xi_coeffs(s, 1)[0]


def xi_coeffs(s0, n):
    x = acb_series([s0, 1], prec=n)
    pre = x * (x - 1) / 2
    ex = (-(x / 2) * arb.pi().log()).exp()
    G = (x / 2).gamma()
    Z = x.zeta()
    return (pre * ex * G * Z).coeffs()


# ---------------- part (i): cells ----------------
def cell_lb(theta, i, j):
    mid = arb(2 * i + 1) / (2 ** (j + 1))
    r = arb(1) / (2 ** (j + 1))
    cm = L_coeffs(acb(theta, mid), 2)
    if cm is None:
        return None
    cb = L_coeffs(acb(theta, mid + r * U), 3)
    if cb is None:
        return None
    pi = arb.pi()
    u0 = cm[0].real / pi
    u1 = -cm[1].imag / pi
    sup2 = 2 * cb[2].abs_upper()
    return u0 - r * abs(u1) - r * r / (2 * pi) * sup2


def cell_task(a):
    th_num, th_den, i0, i1, j0, prec = a
    setprec(prec)
    theta = arb(th_num) / th_den
    stack = [(i, j0) for i in range(i1 - 1, i0 - 1, -1)]
    n = 0
    least = None
    maxj = j0
    fails = []
    nb = 0
    while stack:
        i, j = stack.pop()
        lb = cell_lb(theta, i, j)
        if lb is not None and lb > 0:
            n += 1
            v = float(lb.lower().mid())
            if least is None or v < least[0]:
                least = (v, i, j, pack(lb.lower()))
            maxj = max(maxj, j)
        else:
            if j >= JMAX:
                fails.append((i, j))
                continue
            nb += 1
            stack.append((2 * i + 1, j + 1))
            stack.append((2 * i, j + 1))
    return n, least, maxj, fails, nb


# ---------------- parts (ii), (iii): argument change ----------------
def arg_path(a, b, vert, depth=0):
    """Delta arg xi along the segment a -> b (acb endpoints; vert: same real part), Lemma segment with bisection. Returns (arb, nseg)."""
    xa = xi_val(a)
    xb = xi_val(b)
    c = (a + b) / 2
    if vert:
        d = (b.imag - a.imag).abs_upper()
        ball = acb(c.real, c.imag + (d / 2) * U)
    else:
        d = (b.real - a.real).abs_upper()
        ball = acb(c.real + (d / 2) * U, c.imag)
    sup = xi_coeffs(ball, 2)[1].abs_upper()
    if d * sup < xa.abs_lower():
        return (xb / xa).arg(), 1
    if depth > 40:
        raise RuntimeError("segment bisection too deep")
    s1, n1 = arg_path(a, c, vert, depth + 1)
    s2, n2 = arg_path(c, b, vert, depth + 1)
    return s1 + s2, n1 + n2


def P_ball(Pf):
    return arb(Pf.numerator) / Pf.denominator


def vert_task(a):
    sig_num, sig_den, j, last, Pnum, Pden, prec = a
    setprec(prec)
    sig = arb(sig_num) / sig_den
    Pf = Fraction(Pnum, Pden)
    ya = arb(j) / 8
    yb = P_ball(Pf) if last else arb(j + 1) / 8
    d, n = arg_path(acb(sig, ya), acb(sig, yb), True)
    return pack(d), n


def horiz_task(a):
    start, i, nlen, Pnum, Pden, prec = a
    setprec(prec)
    P = P_ball(Fraction(Pnum, Pden))
    xa = arb(start - i) / 64
    xb = arb(start - i - 1) / 64
    d, n = arg_path(acb(xa, P), acb(xb, P), False)
    return pack(d), n


def horiz_line(pool, start, npieces, Pf, prec):
    tasks = [(start, i, npieces, Pf.numerator, Pf.denominator, prec) for i in range(npieces)]
    tot = arb(0)
    nseg = 0
    for t, n in pool.imap_unordered(horiz_task, tasks):
        tot += unpack(t)
        nseg += n
    return tot, nseg


def vert_line(pool, sig_num, sig_den, Pf, prec):
    n_full = math.floor(8 * Pf)
    tasks = [(sig_num, sig_den, j, False, Pf.numerator, Pf.denominator, prec) for j in range(n_full)]
    if Fraction(n_full, 8) < Pf:
        tasks.append((sig_num, sig_den, n_full, True, Pf.numerator, Pf.denominator, prec))
    tot = arb(0)
    nseg = 0
    for t, n in pool.imap_unordered(vert_task, tasks, chunksize=8):
        tot += unpack(t)
        nseg += n
    return tot, nseg


# ---------------- part (iv): signs ----------------
def sign_task(a):
    k0, k1, prec = a
    setprec(prec)
    out = []
    half = arb(1) / 2
    for k in range(k0, k1):
        v = arb(k) / 50
        x = xi_val(acb(half, v))
        re = x.real
        out.append("+" if re > 0 else ("-" if re < 0 else "?"))
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--prec", type=int, default=96)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--out", default=None)
    A = ap.parse_args()
    out = A.out or ("/tmp/cx_u_certify_quick.json" if A.quick else "cx_u_certify.json")
    logf = open(out.replace(".json", ".log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        logf.write(s + "\n")
        logf.flush()

    Pf = Fraction(81, 2) if A.quick else Fraction(14199, 10)
    prec = A.prec
    setprec(prec)
    res = {"prec": prec, "T": str(Pf), "quick": A.quick}
    t0 = time.time()
    pool = Pool(A.procs)
    ncell0 = math.ceil(Pf * 128) if A.quick else 1420 * 128
    top = Fraction(ncell0, 128) if A.quick else Fraction(1420)
    log("cx_u_certify: prec", prec, "T*", Pf, "cells cover [0,%s]" % top, "procs", A.procs)
    # (i)
    for (tn, td, name) in ((3, 4, "3/4"), (7, 8, "7/8")):
        tasks = [(tn, td, i, min(i + 64, ncell0), 7, prec) for i in range(0, ncell0, 64)]
        ncells = 0
        least = None
        maxj = 0
        fails = []
        nb = 0
        for n, lst, mj, f, b in pool.imap_unordered(cell_task, tasks):
            ncells += n
            nb += b
            maxj = max(maxj, mj)
            fails += f
            if lst is not None and (least is None or lst[0] < least[0]):
                least = lst
        i, j = least[1], least[2]
        loc = float(Fraction(2 * i + 1, 2 ** (j + 1)))
        log("(i) theta=%s: cells %d (initial %d, bisections %d), max depth j=%d, fails %d, least cell lower bound %.4e at v=%.6f (radius 2^-%d)"
            % (name, ncells, ncell0, nb, maxj, len(fails), least[0], loc, j + 1))
        s0 = acb(arb(tn) / td, 0)
        c = L_coeffs(s0, 1)
        u0 = c[0].real / arb.pi()
        log("    u_theta(0) =", u0.str(15, radius=True))
        res["cells_" + name] = {"cells": ncells, "initial": ncell0, "bisections": nb, "maxj": maxj, "fails": len(fails),
                                "least": least[0], "least_at_v": loc, "least_radius_exp": -(j + 1), "u0": u0.str(15, radius=True)}
        log("    time %.0f s" % (time.time() - t0))
    # (ii), (iii)
    R, nR = vert_line(pool, 2, 1, Pf, prec)
    Lf, nL = vert_line(pool, 3, 4, Pf, prec)
    log("right line Re s=2: segments", nR, " left line Re s=3/4: segments", nL, " time %.0f s" % (time.time() - t0))
    P = P_ball(Pf)
    ta, na = horiz_line(pool, 128, 80, Pf, prec)   # 2 -> 3/4 in pieces of length 1/64
    tb, nb_ = horiz_line(pool, 48, 16, Pf, prec)   # 3/4 -> 1/2
    log("top segments:", na, nb_)
    pi = arb.pi()
    N = (R + ta + tb) / pi
    D = (R + ta - Lf) / pi
    M = Lf / pi
    log("(ii) N(T*) =", N.str(20, radius=True))
    log("(iii) 2N(3/4,T*) =", D.str(20, radius=True), " contains 0:", D.contains(0))
    log("(iii) M_3/4(T*) =", M.str(20, radius=True))
    res.update({"N": N.str(20, radius=True), "twoN34": D.str(20, radius=True), "M34": M.str(20, radius=True),
                "segments": {"right": nR, "left34": nL, "top": [na, nb_]}, "twoN34_contains_0": bool(D.contains(0)),
                "N_contains_1000": bool(N.contains(1000)),
                "N_rad": float(N.rad()), "twoN34_rad": float(D.rad()), "M34_rad": float(M.rad())})
    # (iv)
    nk = int(Pf * 50)
    tasks = [(k, min(k + 500, nk + 1), prec) for k in range(0, nk + 1, 500)]
    sg = "".join(pool.imap(sign_task, tasks))
    und = sg.count("?")
    ch = sum(1 for p, q in zip(sg, sg[1:]) if "?" not in (p, q) and p != q)
    log("(iv) points %d, undetermined %d, sign changes %d, first sign %s" % (len(sg), und, ch, sg[0]))
    res["signs"] = {"points": len(sg), "undetermined": und, "changes": ch}
    log("total time %.0f s" % (time.time() - t0))
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
