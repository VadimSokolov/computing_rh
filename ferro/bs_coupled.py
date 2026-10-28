# bs_coupled.py: the numerical claims N3.1 to N3.6 of two coupled Bessel spins (Chapter 13, Section fm:sec:bessel, the paragraph "Numerically, for the
# spins (a, l) = (2 pi, 1) and (2 pi, 1.3)" after Proposition fm:prop:levelrep; note B of the merge): the number of zeros in (0, 45) of the transform
# of X_1 + X_2 for J = 0, 1/8, 1/2; how well Proposition fm:prop:levelrep(i) predicts the zeros at J = 1/8; and, for two spins (2 pi, 1), the
# separation of the pairs into which the double zeros of Z_0 split at J = 3/40 and 1/4. Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. q_i(t) proportional to exp(-a_i cosh(l_i t)) on R (m = 1), i = 1, 2; (X_1, X_2) has density proportional to exp(J t_1 t_2) q_1(t_1) q_2(t_2),
# and the transform of X_1 + X_2 is, up to a constant, Z_J(theta) = int int exp(J t_1 t_2) q_1 q_2 exp(i theta (t_1 + t_2)) dt_1 dt_2 = int q_1(u) e^(i theta u)
# (2/l_2) K_{(J u + i theta)/l_2}(a_2) du (the inner integral is the Bessel integral), real and even in theta. phi_i(theta) = int q_i e^(i theta t) dt is
# (2/l_i) K_{i theta/l_i}(a_i) up to the constant, and Z_J = phi_1 phi_2 - J phi_1' phi_2' + O(J^2) (theta is the variable of the book's transform; the
# reviewers' z is 2 theta). Zeros are real positive theta; counts are of zeros in the open interval (0, 45).
# Method. The integrand is entire and decays like exp(-a cos(l eta) cosh(l x)) along Im t = eta for |eta| < pi/(2 l), so the contour may be moved to
# Im t_i = eta_i = f pi/(2 l_i), which removes the cancellation (Z_0 is of order exp(-S theta), S = (pi/2)(1/l_1 + 1/l_2) = 2.78, about 1e-55 at
# theta = 45, while the integrand is O(1)). The price: with g = 1 - f, (1) the loss in double precision is exp(g S theta) and (2) the trapezoid sum of step h
# on the shifted lines (Z = sum_(i,j) W1_i E_ij W2_j, W_i = h q(u) e^(i theta u), E_ij = exp(J u_i v_j), u_i = x_i + i eta_1, v_j = y_j + i eta_2,
# |x| <= X_1, |y| <= X_2) carries the aliasing terms of the Poisson summation, which are the double transform at theta - 2 pi k/h multiplied by
# exp(2 pi (k eta_1 + m eta_2)/h); relative to Z_J they are of size exp(S (2 theta - 2 pi g/h)) (a first attempt with g = 0.03, h = 0.02 showed it:
# relative error 7e-4 at theta = 1). So the contour and the step depend on the band of theta, with g theta_hi S <= 14 and h <= pi g/(theta_hi + 14):
# bands [0,16): g = 0.30, h = 0.03; [16,26): g = 0.19, h = 0.014; [26,36): g = 0.14, h = 0.008; [36,46.5): g = 0.11, h = 0.0055 (second settings: g * 0.85,
# h * 0.7). The remaining loss is below 1e-9 in relative accuracy, and the aliasing below 1e-17. The real part is taken (the imaginary part is
# rounding, logged). Zeros: sign changes of Z_J on a grid of step 0.02 in theta on [0.1, 46], refined by bisection (45 steps). Count check by the argument
# principle on the half rectangle (0, 45) x (-0.3, 0.3) (Z_J is even, real and positive on the imaginary axis: its number of zeros inside is
# (2 D_bottom + D_right)/(2 pi), with the edges closed exactly at the corner, as in gp_flow.py); only valid if no zero lies within 0.05 of theta = 45,
# which is logged. Prediction of Proposition fm:prop:levelrep(i) at J: for a zero alpha of phi_1, alpha + J phi_2'(alpha)/phi_2(alpha); for a zero of
# phi_2, alpha + J phi_1'(alpha)/phi_1(alpha); the error is the distance to the nearest zero of Z_J, and the distance of alpha to the nearest zero of the
# other family is logged beside it (the proposition is asymptotic as J -> 0, its error is O(J^2) with a constant that grows when two zeros of different
# families are close). Equal spins: the zeros of Z_J in (0, 45) are paired consecutively (each double zero of Z_0 = phi_1^2 splits in two).
# Checks (logged): (a) the banded sum against the untilted sum (f = 0) at theta = 1 and 2; (b) Z_J against phi_1 phi_2 - J phi_1' phi_2' + J^2/2 phi_1''
# phi_2'' at theta = 3 for J = 0.01 and 0.02 (the error must scale like J^3); (c) the whole computation is repeated with the second settings and the
# maximum difference of the zeros is logged.
# Output: ferro/bs_coupled.json and ferro/bs_coupled.log (or --out FILE). Case A: spins (2 pi, 1), (2 pi, 1.3) at J = 0, 1/8, 1/2 with the predictions;
# case B: spins (2 pi, 1) twice at J = 3/40, 1/4. Everything is COMPUTED (double precision), not certified.
# Usage: timeout 30 python3 -I ferro/bs_coupled.py --case A|B|all [--quick] [--out FILE]; --quick uses theta < 25 and coarser steps (a few seconds). On the
# cluster: bash misc/tools/hopper_run.sh -c 4 -t 20 -g bs_coupled.json -g bs_coupled.log ferro/bs_coupled.py -- --case all --out bs_coupled.json
import argparse
import json
import os
import sys
import time

import numpy as np

LOG = []
A2PI = 2 * np.pi


def log(s):
    print(s)
    sys.stdout.flush()
    LOG.append(s)


class Spin:
    def __init__(self, a, ell, f, h, X=None):
        self.a, self.ell, self.f, self.h = a, ell, f, h
        self.eta = f * np.pi / (2 * ell)
        c = np.cos(ell * self.eta)
        if X is None:
            # a cos(l eta) cosh(l x) >= 80 + 0.5 x^2 + ...: solve numerically on a grid
            xs = np.arange(0.5, 12, 0.05)
            ok = a * c * np.cosh(ell * xs) > 90 + 0.5 * xs ** 2
            X = float(xs[np.argmax(ok)]) + 0.5
        self.X = X
        n = int(round(X / h))
        self.x = h * np.arange(-n, n + 1)
        self.u = self.x + 1j * self.eta
        self.q = np.exp(-a * np.cosh(ell * self.u)) * h

    def W(self, th, deriv=0):
        """rows: theta values (complex allowed); q(u) (i u)^deriv e^(i theta u) h."""
        th = np.asarray(th, dtype=complex)
        w = np.exp(1j * th[:, None] * self.u[None, :]) * self.q[None, :]
        if deriv:
            w = w * (1j * self.u[None, :]) ** deriv
        return w

    def phi(self, th, deriv=0):
        out = np.empty(len(th), dtype=complex)
        for i in range(0, len(th), 2000):
            out[i:i + 2000] = self.W(th[i:i + 2000], deriv).sum(axis=1)
        return out


BANDS = [(0.0, 16.0, 0.30, 0.03), (16.0, 26.0, 0.19, 0.014), (26.0, 36.0, 0.14, 0.008), (36.0, 1e9, 0.11, 0.0055)]


def bands(scale_g=1.0, scale_h=1.0):
    return [(lo, hi, 1 - g * scale_g, h * scale_h) for (lo, hi, g, h) in BANDS]


class Pair:
    """Z_J for the spins (a1, l1), (a2, l2); the contour and step depend on the band of Re theta."""

    def __init__(self, sp1, sp2, J, bnd):
        self.J, self.bnd = J, bnd
        self.sp = []
        for (lo, hi, f, h) in bnd:
            s1 = Spin(sp1[0], sp1[1], f, h)
            s2 = Spin(sp2[0], sp2[1], f, h)
            E = None if J == 0 else np.exp(J * s1.u[:, None] * s2.u[None, :])
            self.sp.append((lo, hi, s1, s2, E))

    def Z(self, th):
        th = np.asarray(th, dtype=complex)
        out = np.zeros(len(th), dtype=complex)
        re = np.abs(th.real)
        for (lo, hi, s1, s2, E) in self.sp:
            idx = np.nonzero((re >= lo) & (re < hi))[0]
            for i in range(0, len(idx), 300):
                ii = idx[i:i + 300]
                w1 = s1.W(th[ii])
                w2 = s2.W(th[ii])
                out[ii] = w1.sum(axis=1) * w2.sum(axis=1) if E is None else ((w1 @ E) * w2).sum(axis=1)
        return out


class Phi:
    """phi(theta) = int q e^(i theta t) dt, derivatives in theta, with banded contours."""

    def __init__(self, sp, bnd):
        self.sp = [(lo, hi, Spin(sp[0], sp[1], f, h)) for (lo, hi, f, h) in bnd]

    def __call__(self, th, deriv=0):
        th = np.asarray(th, dtype=complex)
        out = np.zeros(len(th), dtype=complex)
        re = np.abs(th.real)
        for (lo, hi, s) in self.sp:
            idx = np.nonzero((re >= lo) & (re < hi))[0]
            for i in range(0, len(idx), 1000):
                ii = idx[i:i + 1000]
                out[ii] = s.W(th[ii], deriv).sum(axis=1)
        return out


def scan_zeros(fun, lo, hi, grid=0.02, steps=45):
    th = np.arange(lo, hi + grid, grid)
    v = fun(th).real
    s = np.sign(v)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    a, b, fa = th[idx].copy(), th[idx + 1].copy(), v[idx]
    for _ in range(steps):
        m = 0.5 * (a + b)
        fm = fun(m).real
        left = np.sign(fm) == np.sign(fa)
        a = np.where(left, m, a)
        fa = np.where(left, fm, fa)
        b = np.where(left, b, m)
    return 0.5 * (a + b)


def arg_count(fun, Tp, y=0.3, grid=0.05):
    sb = np.linspace(0.0, Tp, int(np.ceil(Tp / grid)) + 1)
    hb = fun(sb - 1j * y)
    ab = np.angle(hb)
    ib = (np.diff(ab) + np.pi) % (2 * np.pi) - np.pi
    sr = np.linspace(-y, y, int(round(2 * y / 0.01)) + 1)
    hr = fun(Tp + 1j * sr)
    ar = np.angle(hr)
    ir = (np.diff(ar) + np.pi) % (2 * np.pi) - np.pi
    junc = (ar[0] - ab[-1] + np.pi) % (2 * np.pi) - np.pi
    return (2 * ib.sum() + junc + ir.sum()) / (2 * np.pi), max(np.abs(ib).max(), np.abs(ir).max())


def counts(zs, T):
    return int(np.sum((zs > 0) & (zs < T)))


def caseA(bnd, T, res, tag):
    sp1, sp2 = (A2PI, 1.0), (A2PI, 1.3)
    p1, p2 = Phi(sp1, bnd), Phi(sp2, bnd)
    zs1 = scan_zeros(lambda t: p1(t), 0.1, T + 1)
    zs2 = scan_zeros(lambda t: p2(t), 0.1, T + 1)
    out = {"phi1": zs1.tolist(), "phi2": zs2.tolist()}
    log("[%s] zeros of phi_1 (2 pi, 1) in (0, %g): %d, of phi_2 (2 pi, 1.3): %d" % (tag, T, counts(zs1, T), counts(zs2, T)))
    zJ = {}
    for J in (0.0, 1.0 / 8, 0.5):
        P = Pair(sp1, sp2, J, bnd)
        zz = scan_zeros(lambda t: P.Z(t), 0.1, T + 1)
        c = counts(zz, T)
        nearT = float(np.min(np.abs(zz - T)))
        ac, mi = arg_count(lambda t: P.Z(t), T)
        log("[%s] J = %.4f: zeros of Z_J in (0, %g): %d by sign changes; argument principle on (0, %g) x (-0.3, 0.3): %.6f (largest increment %.3f); nearest zero to %g at distance %.3f" % (tag, J, T, c, T, ac, mi, T, nearT))
        zJ[J] = zz
        out["J=%.4f" % J] = {"count": c, "argcount": ac, "zeros": zz.tolist()}
    J = 1.0 / 8
    pred = []
    for (zf, pa, pb, nm) in ((zs1, p1, p2, "phi1"), (zs2, p2, p1, "phi2")):
        zin = zf[zf < T]
        shift = J * pb(zin, 1) / pb(zin)
        for al, sh in zip(zin, shift):
            pred.append((al, al + sh.real, nm, abs(sh.imag)))
    zz = zJ[J]
    errs = []
    for al, p, nm, im in pred:
        e = float(np.min(np.abs(zz - p)))
        other = zs2 if nm == "phi1" else zs1
        errs.append((e, al, p, nm, float(np.min(np.abs(other - al)))))
    errs.sort(reverse=True)
    es = np.array([e[0] for e in errs])
    shifts = np.array([abs(p - al) for al, p, nm, im in pred])
    log("[%s] first order predictions at J = 1/8 for the %d zeros of phi_1 phi_2 in (0, %g): largest error %.5f (alpha = %.4f of %s, predicted %.4f, other family %.4f away), median %.5f, 90th percentile %.5f; largest predicted shift %.4f; imaginary parts of the shifts below %.1e" % (tag, len(pred), T, errs[0][0], errs[0][1], errs[0][3], errs[0][2], errs[0][4], np.median(es), np.percentile(es, 90), shifts.max(), max(p[3] for p in pred)))
    log("[%s] the five largest errors (error, alpha, predicted, family, distance to the other family): %s" % (tag, "; ".join("%.4f, %.3f, %.3f, %s, %.3f" % (e[0], e[1], e[2], e[3], e[4]) for e in errs[:5])))
    for dist in (0.5, 1.0):
        far = [e for e in errs if e[4] > dist]
        log("[%s] among the %d zeros at distance above %g from the other family: largest error %.5f, median %.5f" % (tag, len(far), dist, far[0][0], np.median([e[0] for e in far])))
        out["pred_far%g" % dist] = {"n": len(far), "max_err": far[0][0], "median_err": float(np.median([e[0] for e in far]))}
    out["pred"] = {"max_err": errs[0][0], "median_err": float(np.median(es)), "n": len(pred), "worst": errs[:5]}
    res[tag] = out
    return out


def caseB(bnd, T, res, tag):
    sp = (A2PI, 1.0)
    p1 = Phi(sp, bnd)
    out = {}
    zphi = scan_zeros(lambda t: p1(t), 0.1, T + 1)
    nphi = counts(zphi, T)
    for J in (3.0 / 40, 0.25):
        P = Pair(sp, sp, J, bnd)
        zz = scan_zeros(lambda t: P.Z(t), 0.1, T + 1)
        c = counts(zz, T)
        ac, mi = arg_count(lambda t: P.Z(t), T)
        zin = zz[:2 * nphi]
        pairs = zin.reshape(-1, 2)
        sep = pairs[:, 1] - pairs[:, 0]
        mid = pairs.mean(axis=1)
        dev = np.abs(mid - zphi[:nphi]).max()
        sq = 2 * np.sqrt(J)
        inside = int(np.sum(pairs[:, 1] < T))
        log("[%s] equal spins (2 pi, 1), J = %.4f: %d zeros of Z_J in (0, %g) by sign changes, argument principle %.6f; %d zeros of phi in (0, %g), so %d pairs (%d with both members below %g); separation of the pairs: mean %.5f, min %.5f (pair %d, at %.3f), max %.5f (pair %d, at %.3f), first %.5f, last %.5f; 2 sqrt(J) = %.5f; mean/(2 sqrt J) = %.5f; largest |midpoint - zero of phi| = %.4f" % (tag, J, c, T, ac, nphi, T, len(pairs), inside, T, sep.mean(), sep.min(), 1 + int(np.argmin(sep)), mid[np.argmin(sep)], sep.max(), 1 + int(np.argmax(sep)), mid[np.argmax(sep)], sep[0], sep[-1], sq, sep.mean() / sq, dev))
        out["J=%.4f" % J] = {"count": c, "argcount": ac, "npairs": len(pairs), "sep": sep.tolist(), "mid": mid.tolist(), "two_sqrt_J": sq, "mean": float(sep.mean()), "min": float(sep.min()), "max": float(sep.max()), "zeros": zz.tolist()}
    res[tag] = out
    return out


def checks(bnd):
    sp1, sp2 = (A2PI, 1.0), (A2PI, 1.3)
    nb = [(0.0, 1e9, 1.0, 0.02)]
    for J in (0.5,):
        P = Pair(sp1, sp2, J, bnd)
        Q = Pair(sp1, sp2, J, nb)
        Q.sp = [(lo, hi, Spin(A2PI, 1.0, 0.0, 0.02, X=7.0), Spin(A2PI, 1.3, 0.0, 0.02, X=5.5), E) for (lo, hi, s1, s2, E) in Q.sp]
        s1, s2 = Q.sp[0][2], Q.sp[0][3]
        Q.sp[0] = (0.0, 1e9, s1, s2, np.exp(J * s1.u[:, None] * s2.u[None, :]))
        for th in (1.0, 2.0, 10.0):
            a = P.Z(np.array([th]))[0]
            b = Q.Z(np.array([th]))[0]
            log("check (a): J = %g, theta = %g: banded %.12e, untilted %.12e, relative difference %.1e, imaginary/real of the banded %.1e" % (J, th, a.real, b.real, abs(a - b) / abs(b), abs(a.imag / a.real)))
    th = np.array([3.0])
    p1, p2 = Phi(sp1, bnd), Phi(sp2, bnd)
    errs = []
    for J in (0.01, 0.02):
        z = Pair(sp1, sp2, J, bnd).Z(th)[0]
        ex = p1(th)[0] * p2(th)[0] - J * p1(th, 1)[0] * p2(th, 1)[0] + J ** 2 / 2 * p1(th, 2)[0] * p2(th, 2)[0]
        errs.append(abs(z - ex))
    log("check (b): theta = 3, |Z_J - three term expansion| = %.3e (J = 0.01), %.3e (J = 0.02), ratio %.2f (J^3 scaling gives 8)" % (errs[0], errs[1], errs[1] / errs[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="all")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--T", type=float, default=45.0)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    T = 25.0 if a.quick else a.T
    out = a.out or os.path.join(here, "bs_coupled_quick.json" if a.quick else "bs_coupled.json")
    res = {"T": T, "bands": BANDS, "second": "g * 0.85, h * 0.7"}
    t0 = time.time()
    checks(bands())
    for case, fn in (("A", caseA), ("B", caseB)):
        if a.case not in (case, "all"):
            continue
        r1, r2 = {}, {}
        fn(bands(), T, r1, case + "1")
        fn(bands(0.85, 0.7), T, r2, case + "2")
        res[case] = r1[case + "1"]
        res[case + "_second"] = r2[case + "2"]
        # difference of zeros
        dm = 0.0
        for k, v in r1[case + "1"].items():
            if isinstance(v, dict) and "zeros" in v:
                v2 = r2[case + "2"][k]["zeros"]
                if len(v2) == len(v["zeros"]):
                    dm = max(dm, float(np.abs(np.array(v["zeros"]) - np.array(v2)).max()))
                else:
                    log("[%s] ZERO COUNT DIFFERS between settings for %s: %d against %d (zeros up to %g + 1)" % (case, k, len(v["zeros"]), len(v2), T))
        log("[%s] maximum difference of the zeros of Z_J between the two settings: %.2e" % (case, dm))
        res[case]["maxdiff_settings"] = dm
    log("time %.1f s" % (time.time() - t0))
    json.dump(res, open(out, "w"), indent=1)
    open(out.replace(".json", ".log"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
