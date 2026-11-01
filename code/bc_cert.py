# bc_cert.py: a certified version of Proposition barrier:prop:check of ch/nc_barrier.tex (Section barrier:sec:check, Table barrier:tab:check),
# which code/nc_barrier_check.py computes in floating point. Ball arithmetic (python-flint, Arb) throughout; nothing is tuned.
#
# Setting. xi_t(s) = 4 int_R e^{t u^2} Phi(u) e^{2 w u} du, w = s - 1/2 = x + i theta, Phi = Phi_dBN (code/phase_rig.py, Section cp:sec:march).
# B = [0.6, 1] x [497.4, 497.9], t in [0, 0.2]. G_k(s) = d_t^k xi_t(s) at the base time t_b: the same integral with u^{2k} e^{t_b u^2}.
# What is proved (every inequality below is evaluated in ball arithmetic, upper and lower bounds are the exact endpoints of balls):
#  (A) the evaluator. S_k(s) = 8 h sum_{j=0}^{J} omega_j u_j^{2k} e^{t_b u_j^2} Phi(u_j) cosh(2 w u_j), u_j = j h, h = 1/512, J = 1536 (U = 3), and
#      G_k - S_k = E_k with |E_k| <= Etrap + Ecut, Etrap = 2 cosh(2 a theta) Kt_k(a, |x|) / (e^{2 pi a / h} - 1) by Poisson summation on the strip
#      |Im u| < a = 0.30 (f_k analytic there; Kt_k(eta, x) = 8 int_0^oo |u + i eta|^{2k} e^{t_b (u^2 - eta^2)} |Phi(u + i eta)| cosh(2 x u) du),
#      Ecut = 16 G(U)/lambda for the nodes beyond U and the integral beyond U (|Phi(u)| <= 60 e^{9u} e^{-pi e^{4u}}, u >= 3, G decreasing, -G'/G >= lambda).
#      Kt_k is bounded above by interval cells (1536 cells of [0, 3], |Phi(v + i eta)| by phase_rig.phi, the tail beyond 3 by the kernel bound
#      |Phi(u + i eta)| <= 60 e^{9u} exp(-pi cos(4 eta) e^{4u}) and the same lambda method).
#  (B) the tilted bound (cp:lem:tilt, for the derivatives G_k): |G_k(z)| <= Kt_k(eta, |x|) e^{-2 eta Im z}, Im z >= 0, eta = 0.38 < pi/8.
#  (C) Taylor model of G_k on a disc |z - c| <= rho: p_k(z) = sum_{j<=N} c_{kj} (z - c)^j (c_{kj} the exact Taylor coefficients of S_k at c, balls)
#      and Rem_k = Emax + (M + Emax) q^{N+1}/(1 - q), q = rho / r, r = 1/2, M = Kt_k(0.38, |x_c| + r) e^{-2 eta (theta_c - r)}, Emax = sup of the bound of
#      (A) on the disc of radius r (Cauchy's estimate for the coefficients of S_k beyond N, whose sup norm is at most M + Emax). For the derivative in s:
#      Rem'_k = Emax/(r - rho) + (M + Emax)/r * sum_{j>N} j q^{j-1}.
#  (D) the Taylor tail in t: |sum_{k>8} D^k/k! G_k(s)| <= 8 e^{-2 eta 497.4} int e^{t_b (u^2 - eta^2)} |Phi(u + i eta)| cosh(u) r^9 e^r / 9! du,
#      r = D |u + i eta|^2 (s in B: x <= 1/2, theta >= 497.4); for the s-derivative the integrand carries the extra factor 2 |u + i eta|.
#  (E) boundary: cells are segments, discs about their midpoints. On a cell with ell = |c_00| - sum_{j>=1} |c_0j| rho^j - Rem_0 > 0 the image of the cell lies in the open
#      disc about c_00 that excludes 0 (the analogue of cp:lem:grid), so the change of arg along it is Arg(xi(end)/c_00) + Arg(c_00/xi(start)); the sum over the
#      cells is a ball that must lie in (-pi, pi) for winding number 0 (for the control, in (pi, 3 pi) for winding number 1).
#      If moreover sum_{k=1}^8 D^k/k! (sum_j |c_kj| rho^j + Rem_k) + Tail(D) < ell on every cell of the boundary, then |xi_t(s) - xi_{t_b}(s)| < |xi_{t_b}(s)| for every
#      s on the boundary and every t in [t_b, t_b + D]; with t_b = 0, D = 0.2, xi_t has no zero on the boundary for t in [0, 0.2], so the winding number is the same for
#      every t (continuity), and the winding number at t = 0 is certified by the same cells.
#  (F) density: v_t = Re(xi_t'/xi_t) > 0 on B x [0, 0.2]. Cells: squares of side 0.1 (discs of radius 0.0708) times time intervals [t_b, t_b + 0.02]. At the base time,
#      L = xi'/xi has a Taylor model Lambda(z) = sum_{j<=14} l_j z^j (power series division), with  |L - Lambda| <= (sum_j |d_j| rho^j + Rem'_0 + Rem_0 sum |l_j| rho^j) / ell,
#      d = p_0' - p_0 Lambda exactly as a polynomial; hence Re L >= Re l_0 - sum_{j>=1} |l_j| rho^j - that. In time, with A = sum_{k=1}^8 D^k/k! (abs of p_k + Rem_k) + Tail and
#      B likewise for the derivative, |L_t - L_{t_b}| <= (B + max|L| A) / (ell - A). The minimum over all cells and base times is a certified lower bound of v on B x [0, 0.2].
# Not proved here: nothing in the chapter's other claims (winding at the 21 grid times and the 35 points are floating point, replaced by the statements for all t).
# Usage: python3 bc_cert.py [--quick] [--workers N] [--out FILE]
import argparse
import json
import math
import multiprocessing as mpr
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import acb, arb, ctx  # noqa: E402
import phase_rig as pr  # noqa: E402

PREC = 1000
KPREC = 128
H = arb(1) / 512
J = 1536
NODES = J + 1
N = 24
MQ = 14
KMAX = 8
R = arb(1) / 2
ETA_S = "0.38"
A_S = "0.30"
NCELL_K = 3 * 2 ** 9
T0 = None
T1 = None
TH_MIN = 497.4
HALF = arb(1) / 2

PHI = None
KB = None


def upper(x):
    return x.upper()


def lower(x):
    return x.lower()


def abs_up(z):
    return abs(z).upper()


def gtail_gen(U, A, B, m, mult, c, cst=None):
    """mult * cst * 60 G(U) / lambda >= mult * cst * int_U^oo 60 exp(A u^2 + B u + m ln(u + 1/2) - pi c e^{4u}) du, U >= 3, c <= cos(4 eta) (c = 1 on the real axis)."""
    pi = arb.pi()
    e4U = (4 * U).exp()
    lam = 4 * pi * c * e4U - 2 * A * U - B - m / (U + HALF)
    if not (U >= 3 and lam > 0 and pi * c * e4U > 1000):
        raise RuntimeError("tail hypotheses fail")
    G = (A * U * U + B * U + m * (U + HALF).log() - pi * c * e4U).exp()
    out = mult * 60 * G / lam
    if cst is not None:
        out = out * cst
    return out.upper()


class KBounds:
    """Upper bounds Kt_k(eta, x; t_b) and the Taylor tails, with the cells of |Phi(v + i eta)| cached."""

    def __init__(self):
        self.cells = {}
        self.cache = {}

    def phicells(self, eta_s):
        if eta_s not in self.cells:
            old = ctx.prec
            ctx.prec = KPREC
            eta = arb(eta_s)
            w = arb(1) / 2 ** 9
            out = []
            for i in range(NCELL_K):
                lo = w * i
                ur = lo.union(lo + w)
                ph = abs(pr.phi(acb(ur, eta), KPREC, True)).upper()
                if not (ph == ph):
                    raise RuntimeError("indeterminate Phi cell")
                out.append((ur, ph, (lo + w)))
            self.cells[eta_s] = out
            ctx.prec = old
        return self.cells[eta_s]

    def Kt(self, eta_s, k, tb_s, x_s):
        key = ("K", eta_s, k, tb_s, x_s)
        if key in self.cache:
            return self.cache[key]
        old = ctx.prec
        ctx.prec = KPREC
        eta, tb, x = arb(eta_s), arb(tb_s), arb(x_s)
        w = arb(1) / 2 ** 9
        tot = arb(0)
        for ur, ph, hi in self.phicells(eta_s):
            F = (ur * ur + eta * eta) ** k * (tb * (ur * ur - eta * eta)).exp() * ph
            tot += F.upper() * (2 * x * hi).cosh().upper()
        c = (4 * eta).cos().lower()
        K = (8 * w * tot).upper() + gtail_gen(arb(3), tb, 9 + 2 * x, 2 * k, 8, c)
        K = K.upper()
        ctx.prec = old
        self.cache[key] = K
        return K

    def tail_t(self, tb_s, D_s, pw):
        """Upper bound for sum_{k>8} D^k/k! |G_k(s)| (pw = 0) or |G_k'(s)| (pw = 1) over s in B, G_k at the base time tb."""
        key = ("T", tb_s, D_s, pw)
        if key in self.cache:
            return self.cache[key]
        old = ctx.prec
        ctx.prec = KPREC
        eta, tb, D = arb(ETA_S), arb(tb_s), arb(D_s)
        w = arb(1) / 2 ** 9
        tot = arb(0)
        for ur, ph, hi in self.phicells(ETA_S):
            v2 = ur * ur + eta * eta
            r = D * v2
            F = (tb * (ur * ur - eta * eta)).exp() * ph * (r ** 9) * r.exp() / 362880
            if pw == 1:
                F = F * v2.sqrt()
            tot += F.upper() * hi.cosh().upper()
        c = (4 * eta).cos().lower()
        mult = 8 * (2 if pw == 1 else 1)
        T = (mult * w * tot).upper() + gtail_gen(arb(3), tb + 2 * D, 10 + 2 * D, pw, mult, c, cst=(D / 2).exp())
        T = (T * (-2 * eta * TH_MIN).exp()).upper()
        ctx.prec = old
        self.cache[key] = T
        return T


def init_proc():
    global PHI, KB, T0, T1
    ctx.prec = PREC
    T0, T1 = arb("497.4"), arb("497.9")
    PHI = [pr.phi(H * j, PREC, False) for j in range(NODES)]
    KB = KBounds()


class Model:
    pass


def xm_string(sig_c):
    v = abs(float(sig_c) - 0.5) + 0.5
    return "%.3f" % ((math.ceil(v * 1000) + 1) / 1000)


def cell_model(tb_s, C, rho, kmax=KMAX):
    """Taylor model of G_k (k = 0..kmax) at the base time tb on the disc of radius rho about the acb C."""
    ctx.prec = PREC
    tb = arb(tb_s)
    m = Model()
    m.C, m.rho, m.tb_s = C, rho, tb_s
    wc = C - HALF
    u = [H * j for j in range(NODES)]
    ew = [(2 * wc * ui).exp() for ui in u]
    ewi = [1 / e for e in ew]
    om = [HALF if j == 0 else arb(1) for j in range(NODES)]
    ph_t = [(tb * ui * ui).exp() * PHI[j] for j, ui in enumerate(u)]
    # weights W_k,i = 8 h omega_i u_i^{2k} e^{tb u_i^2} Phi(u_i)
    W = []
    for k in range(kmax + 1):
        W.append([8 * H * om[j] * (u[j] ** (2 * k)) * ph_t[j] for j in range(NODES)])
    pw = [arb(1)] * NODES
    coef = [[None] * (N + 1) for _ in range(kmax + 1)]
    for j in range(N + 1):
        if j > 0:
            pw = [pw[i] * (2 * u[i]) / j for i in range(NODES)]
        sg = 1 if j % 2 == 0 else -1
        Z = [pw[i] * (ew[i] + sg * ewi[i]) / 2 for i in range(NODES)]
        for k in range(kmax + 1):
            Wk = W[k]
            s = acb(0)
            for i in range(NODES):
                s += Wk[i] * Z[i]
            coef[k][j] = s
    m.coef = coef
    # remainders
    sig_c = float(C.real.mid())
    th_c = float(C.imag.mid())
    xs = xm_string(sig_c)
    xm = arb(xs)
    th_lo = C.imag - R
    th_hi = C.imag + R
    ctx.prec = KPREC
    pi = arb.pi()
    q = (rho / R).upper()
    tg = (q ** (N + 1) / (1 - q)).upper()
    tgp = (((N + 1) * q ** N * (1 - q) + q ** (N + 1)) / (1 - q) ** 2).upper()
    expo = ((2 * pi * arb(A_S) / H).exp() - 1).lower()
    m.rem, m.remp, m.M, m.Emax = [], [], [], []
    for k in range(kmax + 1):
        Kc = KB.Kt(ETA_S, k, tb_s, xs)
        M = (Kc * (-2 * arb(ETA_S) * th_lo).exp()).upper()
        Ks = KB.Kt(A_S, k, tb_s, xs)
        Etr = (2 * (2 * arb(A_S) * th_hi).cosh() * Ks / expo).upper()
        Ecut = gtail_gen(arb(3), arb(tb_s), 9 + 2 * xm, 2 * k, 16, arb(1))
        E = (Etr + Ecut).upper()
        rem = (E + (M + E) * tg).upper()
        remp = (E / (R - rho) + (M + E) / R * tgp).upper()
        m.rem.append(rem)
        m.remp.append(remp)
        m.M.append(M)
        m.Emax.append(E)
    ctx.prec = PREC
    return m


def abs_poly_sum(cs, rho, start=0, dermul=False):
    """Upper bound of sum_{j>=start} |c_j| rho^j (or j |c_j| rho^{j-1} when dermul)."""
    tot = arb(0)
    rp = rho ** start if not dermul else (rho ** (start - 1) if start >= 1 else arb(0))
    for j in range(start, len(cs)):
        if dermul:
            tot += (abs_up(cs[j]) * j * (rho ** (j - 1) if j >= 1 else arb(0))).upper()
        else:
            tot += (abs_up(cs[j]) * rho ** j).upper()
    return tot.upper()


def ell_of(m):
    c0 = m.coef[0]
    return (abs(c0[0]).lower() - abs_poly_sum(c0, m.rho, 1) - m.rem[0]).lower()


def a_dev(m, D_s, tail_s, deriv=False):
    """sum_{k=1}^8 D^k/k! (sup of |G_k| (or |G_k'|) on the disc by the model) + tail."""
    D = arb(D_s)
    tot = arb(0)
    for k in range(1, KMAX + 1):
        if deriv:
            sup = abs_poly_sum(m.coef[k], m.rho, 1, True) + m.remp[k]
        else:
            sup = abs_poly_sum(m.coef[k], m.rho, 0) + m.rem[k]
        tot += (D ** k / math.factorial(k) * sup).upper()
    return (tot + tail_s).upper()


def poly_mul(a, b):
    out = [acb(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = out[i + j] + x * y
    return out


def wind_piece(m, S, E, ctrl):
    """Ball for the argument change along the segment from S to E (inside the disc of m), and ell; None if the cell condition fails."""
    c0 = list(m.coef[0])
    rem = m.rem[0]
    if ctrl:
        w0 = acb(arb("0.25"), arb("497.65"))
        wc = m.C - HALF
        f1 = [wc * wc - w0 * w0, 2 * wc, acb(1)]
        f2 = [wc * wc - w0.conjugate() ** 2, 2 * wc, acb(1)]
        Q = poly_mul(f1, f2)
        c0 = poly_mul(c0, Q)
        Qmax = abs_poly_sum(Q, m.rho, 0)
        rem = (rem * Qmax).upper()
    A = c0[0]
    ell = (abs(A).lower() - abs_poly_sum(c0, m.rho, 1) - rem).lower()
    if not (ell > 0):
        return None, ell

    def val(P):
        z = P - m.C
        s = acb(0)
        for j in range(len(c0) - 1, -1, -1):
            s = s * z + c0[j]
        return s + acb(arb(0, rem), arb(0, rem))

    vs, ve = val(S), val(E)
    ang = (ve / A).arg() + (A / vs).arg()
    return ang, ell


def boundary_cells(depth=0):
    pts = []
    s0, s1 = arb("0.6"), arb("1")
    n = 4
    # counterclockwise: bottom (0.6 -> 1 at T0), right (T0 -> T1 at 1), top (1 -> 0.6 at T1), left (T1 -> T0 at 0.6)
    bottom = [(acb(s0 + (s1 - s0) * i / n, T0), acb(s0 + (s1 - s0) * (i + 1) / n, T0)) for i in range(n)]
    right = [(acb(s1, T0 + (T1 - T0) * i / 5), acb(s1, T0 + (T1 - T0) * (i + 1) / 5)) for i in range(5)]
    top = [(acb(s1 - (s1 - s0) * i / n, T1), acb(s1 - (s1 - s0) * (i + 1) / n, T1)) for i in range(n)]
    left = [(acb(s0, T1 - (T1 - T0) * i / 5), acb(s0, T1 - (T1 - T0) * (i + 1) / 5)) for i in range(5)]
    return bottom + right + top + left


def boundary_job(args):
    """tb, list of D (strings), ctrl flag: the cells of the boundary with ell, the winding sum and the Taylor criterion for each D."""
    tb_s, Ds, ctrl = args
    t_start = time.time()
    ctx.prec = PREC
    D_tail = {D: KB.tail_t(tb_s, D, 0) for D in Ds}
    total = arb(0)
    cells = []
    stack = [(S, E, 0) for S, E in reversed(boundary_cells())]
    out = []
    while stack:
        S, E, d = stack.pop()
        C = (S + E) / 2
        rho = (abs(E - S) / 2).upper()
        m = cell_model(tb_s, C, rho, kmax=0 if ctrl else KMAX)
        ang, ell = wind_piece(m, S, E, ctrl)
        row = dict(S=[S.real.str(8, radius=False), S.imag.str(10, radius=False)], E=[E.real.str(8, radius=False), E.imag.str(10, radius=False)],
                   rho=float(rho.mid()), ell=ell.str(6, radius=False), ell_f=float(ell.mid()))
        ok = ang is not None
        if ok and not ctrl:
            for D in Ds:
                A = a_dev(m, D, D_tail[D])
                row["ratio_D" + D] = float((A / ell).upper().mid())
                row["ratio_floor_D" + D] = float((a_dev(m, D, arb(0)) / ell).upper().mid())
                ok = ok and bool(A < ell)
        if not ok:
            if d >= 4:
                return dict(tb=tb_s, ok=False, failed_cell=row)
            mid = (S + E) / 2
            stack.append((mid, E, d + 1))
            stack.append((S, mid, d + 1))
            continue
        total += ang
        row["arg"] = ang.str(8)
        out.append(row)
    pi = arb.pi()
    res = dict(tb=tb_s, ok=True, ctrl=ctrl, ncells=len(out), cells=out, argsum=total.str(12), winding=(total / (2 * pi)).str(8),
               seconds=round(time.time() - t_start, 1))
    if ctrl:
        res["winding_is_one"] = bool(total > pi and total < 3 * pi)
    else:
        res["winding_is_zero"] = bool(abs(total).upper() < pi)
        res["min_ell"] = min(c["ell_f"] for c in out)
        for D in Ds:
            res["max_ratio_D" + D] = max(c["ratio_D" + D] for c in out)
            res["max_ratio_floor_D" + D] = max(c["ratio_floor_D" + D] for c in out)
            res["tail_D" + D] = D_tail[D].str(6, radius=False)
    return res


def interior_job(args):
    """Cells (4 x 5 squares of side 0.1) of B at the base time tb, step D: certified lower bound of v on B x [tb, tb + D] (D = '0': at tb only)."""
    tb_s, D_s = args
    t_start = time.time()
    ctx.prec = PREC
    pi = arb.pi()
    Dz = D_s != "0"
    tailA = KB.tail_t(tb_s, D_s, 0) if Dz else arb(0)
    tailB = KB.tail_t(tb_s, D_s, 1) if Dz else arb(0)
    rho = arb("0.0708")
    best = None
    bestb = None
    cells = []
    for i in range(4):
        sg = arb("0.65") + arb("0.1") * i
        for jj in range(5):
            th = arb("497.45") + arb("0.1") * jj
            C = acb(sg, th)
            m = cell_model(tb_s, C, rho)
            ell = ell_of(m)
            if not (ell > 0):
                return dict(tb=tb_s, D=D_s, ok=False, cell=[float(sg.mid()), float(th.mid())], ell=ell.str(5))
            Pa = list(m.coef[0])
            Pb = [(j + 1) * Pa[j + 1] for j in range(N)]
            # Taylor series of L = Pb / Pa to order MQ
            lj = []
            for j in range(MQ + 1):
                s = Pb[j] if j < len(Pb) else acb(0)
                for ii in range(1, j + 1):
                    s = s - Pa[ii] * lj[j - ii]
                lj.append(s / Pa[0])
            prod = poly_mul(Pa, lj)
            dpoly = [(Pb[j] if j < len(Pb) else acb(0)) - prod[j] for j in range(len(prod))]
            sumd = abs_poly_sum(dpoly, rho, 0)
            suml = abs_poly_sum(lj, rho, 1)
            remL = ((sumd + m.remp[0] + (m.rem[0] * (abs_up(lj[0]) + suml))) / ell).upper()
            vb = (lj[0].real - suml - remL).lower()
            Lmax = (abs_up(lj[0]) + suml + remL).upper()
            if Dz:
                A = a_dev(m, D_s, tailA)
                Bd = a_dev(m, D_s, tailB, deriv=True)
                if not (ell > A):
                    return dict(tb=tb_s, D=D_s, ok=False, cell=[float(sg.mid()), float(th.mid())], ell=ell.str(5), A=A.str(5))
                dL = ((Bd + Lmax * A) / (ell - A)).upper()
            else:
                dL = arb(0)
            vl = (vb - dL).lower()
            cells.append(dict(c=[float(sg.mid()), float(th.mid())], v_at_tb=float(vb.mid()), dL=float(dL.mid()), vlow=float(vl.mid()),
                              ell=float(ell.mid()), remL=float(remL.mid())))
            if best is None or vl < best:
                best = vl
            if bestb is None or vb < bestb:
                bestb = vb
    return dict(tb=tb_s, D=D_s, ok=True, vlow=best.str(10, radius=False), vlow_at_tb=bestb.str(10, radius=False),
                vlow_f=float(best.mid()), vlow_at_tb_f=float(bestb.mid()), cells=cells, seconds=round(time.time() - t_start, 1))


def set_h(den):
    global H, J, NODES, PHI
    H = arb(1) / den
    J = 3 * den
    NODES = J + 1
    ctx.prec = PREC
    PHI = [pr.phi(H * j, PREC, False) for j in range(NODES)]


def selftest(_):
    """Independent checks (floating point mpmath at 80 digits, and a second step h = 1/1024): the cell model about (0.8, 497.65), the corners of B,
    the values of v_t at the points where the chapter reports the minimum, and the heat equation d_t xi = xi''/4."""
    import mpmath as mp
    mp.mp.dps = 80
    out = {}
    ctx.prec = PREC

    def xi_mp(s):
        return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

    def mpv(x):
        return mp.mpf(x.mid().str(75, radius=False))

    C = acb(arb("0.8"), arb("497.65"))
    m = cell_model("0", C, arb("0.0708"))
    out.update(rho=m.rho.str(5, radius=False), ell_lb=ell_of(m).str(6, radius=False), rem0=m.rem[0].str(5, radius=False), M0=m.M[0].str(5, radius=False),
               Emax0=m.Emax[0].str(5, radius=False), rem8=m.rem[8].str(5, radius=False), M8=m.M[8].str(5, radius=False))
    diffs = []
    for zr, zi in (("0", "0"), ("0.03", "0.02"), ("0", "-0.05"), ("0.06", "-0.01")):
        zz = acb(arb(zr), arb(zi))
        pz = acb(0)
        for j in range(N, -1, -1):
            pz = pz * zz + m.coef[0][j]
        ex = xi_mp(mp.mpc(mp.mpf("0.8") + mp.mpf(zr), mp.mpf("497.65") + mp.mpf(zi)))
        pv = mp.mpc(mpv(pz.real), mpv(pz.imag))
        diffs.append(float(abs(pv - ex) / abs(ex)))
    out["cell_rel_diff_vs_mpmath"] = diffs
    hd = []
    for j in range(0, 6):
        d = m.coef[1][j] - (j + 1) * (j + 2) * m.coef[0][j + 2] / 4
        hd.append(float((abs(d).upper() / abs(m.coef[1][j]).lower()).mid()))
    out["heat_rel"] = hd
    # corners at t = 0 and the point values of v_t, for two steps h
    pts = {"corners": [("0.6", "497.4"), ("1", "497.4"), ("1", "497.9"), ("0.6", "497.9")], "v": [("0.6", "497.45"), ("0.6", "497.5")]}
    rows = {}
    for den in (512, 1024):
        if den != 512:
            set_h(den)
        r = {}
        for sg, th in pts["corners"]:
            mm = cell_model("0", acb(arb(sg), arb(th)), arb("1e-30"))
            c00, c01 = mm.coef[0][0], mm.coef[0][1]
            ex = xi_mp(mp.mpc(mp.mpf(sg), mp.mpf(th)))
            exd = mp.diff(lambda z: xi_mp(z), mp.mpc(mp.mpf(sg), mp.mpf(th)))
            r["xi(%s+%si)" % (sg, th)] = dict(abs=c00.real.str(1, radius=False) and abs(c00).str(12), rel_xi_vs_mpmath=float(abs(mp.mpc(mpv(c00.real), mpv(c00.imag)) - ex) / abs(ex)),
                                               rel_dxi_vs_mpmath=float(abs(mp.mpc(mpv(c01.real), mpv(c01.imag)) - exd) / abs(exd)), radius_rel=float((abs(c00).rad() / abs(c00).lower()).mid()))
        for t_s in ("0", "0.17", "0.18", "0.2"):
            for sg, th in pts["v"]:
                mm = cell_model(t_s, acb(arb(sg), arb(th)), arb("1e-30"))
                v = (mm.coef[0][1] / mm.coef[0][0]).real
                r["v(t=%s;%s,%s)" % (t_s, sg, th)] = v.str(15)
        rows["h=1/%d" % den] = r
    out["checks"] = rows
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t_start = time.time()
    here = os.path.dirname(os.path.abspath(__file__))
    res = {}
    if a.selftest:
        init_proc()
        res["selftest"] = selftest(0)
        print(json.dumps(res, indent=1))
        dd = os.path.join(here, "..", "data")
        json.dump(res, open(os.path.join(dd if os.path.isdir(dd) else here, "bc_selftest.json"), "w"), indent=1)
        return
    ctx.prec = PREC
    jobsB = [("0", ["0.2", "0.01"], False), ("0", [], True), ("0.19", ["0.01"], False)]
    ts = [("%.2f" % (0.02 * i), "0.02") for i in range(10)] + [("0.20", "0")]
    ts += [("0.00", "0")]
    if a.quick:
        jobsB, ts = jobsB[:1], ts[:1]
    with mpr.Pool(a.workers, initializer=init_proc) as pool:
        rB = pool.map(boundary_job, jobsB, chunksize=1)
        rI = pool.map(interior_job, ts, chunksize=1)
    res["boundary"] = rB
    res["interior"] = rI
    okI = [r for r in rI if r["ok"] and r["D"] != "0"]
    if okI:
        res["vlow_all_t"] = min(r["vlow_f"] for r in okI)
        res["vlow_all_t_ball"] = min(okI, key=lambda r: r["vlow_f"])["vlow"]
    res["seconds"] = round(time.time() - t_start, 1)
    out = a.out or os.path.join(here, "..", "data", "bc_cert%s.json" % ("_quick" if a.quick else ""))
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("boundary", "interior")}, indent=1))
    for r in rB:
        print({k: v for k, v in r.items() if k != "cells"})
    for r in rI:
        print({k: v for k, v in r.items() if k != "cells"})


if __name__ == "__main__":
    main()
