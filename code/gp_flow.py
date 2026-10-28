# gp_flow.py: the numerical claims N2.1 to N2.9 of the forward flow (Chapter 21, Section ch:dbn, the paragraph after "forward in time the zeros repel
# each other", note G3 of the merge): the first nz = 120 positive zeros of H_t(2 theta), t = 0, 0.05, 0.1, 0.2, unfolded by x_j = theta(theta_j)/pi
# (theta the Riemann Siegel theta function), the population variance of the gaps x_(j+1) - x_j and their minimum. Added for the book in October
# 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. H_t(z) = int_0^infty exp(t u^2) Phi_dBN(u) cos(z u) du with Phi_dBN(u) = sum_(n >= 1) (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u)) exp(-pi n^2
# e^(4u)) (Section sec:ch21:flow); a zero z of H_t is the ordinate theta = z/2, and the unit of theta is the unit of the ordinates gamma_n of zeta
# (H_0(2 theta) = Xi(theta)/8). The flow time t is the book's. The zeros theta_j are the positive ones, j = 1, 2, ...; x_j = vartheta(theta_j)/pi is
# dimensionless (mean gap 1), vartheta(theta) = Im log Gamma(1/4 + i theta/2) - (theta/2) log pi. The variance is the population variance
# (divisor n - 1 gaps, not n - 2) of the nz - 1 = 119 gaps about their mean; the minimum is the smallest gap.
# Method. The integrand is entire, so the real line may be replaced by the line Im u = eta, and the trapezoid sum of step h on it is
#   H_t(z) = (h/2) e^(-z eta) [A_0 + sum_(j >= 1) (A_j e^(i z x_j) + conj(A_j) e^(-i z x_j))],   A_j = e^(t u_j^2) Phi_dBN(u_j),  u_j = j h + i eta,
# (Phi_dBN even and real on the real axis gives A_(-j) = conj(A_j)); for real z this is (h/2) e^(-z eta)[A_0 + 2 Re sum A_j e^(i z x_j)]. The shift
# removes the cancellation of the real line integral (H_0(z) is of order e^(-pi z/8) while the integrand is O(1)): with eta close to pi/8 the loss is
# e^((pi/8 - eta) z), about 1e5 in double precision at z = 524. Phi_dBN(u_j) decays like exp(-pi e^(4 x) cos(4 eta)), so |x| <= xmax = 2.6 and
# n <= nmax = 40 (the last terms are logged). By Poisson summation the sum equals H_t(z) plus H_t(z +- 2 pi k/h), which is below 1e-100 relative
# for h <= 0.0025 at z <= 524. The two settings are (eta, h) = (0.38, 0.00125) and (0.37, 0.0025); the eta = 0.385, h = 0.005 first tried had a relative error of 1e-1 at theta = 5 (checked against the Riemann xi function at 40 digits, which
# agrees with the final settings to 1e-13 at theta <= 250), and eta = 0.2 loses all accuracy by theta = 100. All in double precision (numpy), so the zeros are COMPUTED, not certified.
# Zeros: sign changes of H_t(2 theta) on the grid of step 0.04 in theta, refined by bisection (vectorised, 45 steps, final width 1e-14).
# Count check by the argument principle on the half rectangle (0, T') x (-1/2, 1/2) in theta (T' the midpoint of theta_nz and theta_(nz+1)):
# H_t(2 theta) is even, real and positive on the imaginary axis, so the number of zeros inside is (2 D_bottom + D_right)/(2 pi), D the change of the
# continuous argument along the bottom edge (step 0.05) and the right edge (step 0.01); the increments must stay below pi/2 (logged). It must equal nz.
# Error control: the whole computation is repeated with a second (eta, h) and the maximum difference of the zeros theta_j is logged; at t = 0 the
# zeros are compared with the ordinates of data/zeros_1700.txt (60 digits read, converted to double).
# Output: data/gp_flow.json and the log data/gp_flow.log (or --out FILE). --quick: t = 0 and 0.2, nz = 30 (seconds).
# Usage: timeout 30 python3 -I code/gp_flow.py [--t 0 0.05 0.1 0.2] [--nz 120] [--quick] [--out FILE]. The default run takes about 10 seconds on the
# laptop; on the cluster: bash misc/tools/hopper_run.sh -c 1 -t 10 -g gp_flow.json -g gp_flow.log code/gp_flow.py data/zeros_1700.txt -- --out gp_flow.json
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.special import loggamma

LOG = []
TWOPI = 2 * np.pi


def log(s):
    print(s)
    sys.stdout.flush()
    LOG.append(s)


def vartheta(th):
    th = np.asarray(th, dtype=float)
    return loggamma(0.25 + 0.5j * th).imag - 0.5 * th * np.log(np.pi)


class Evaluator:
    """H_t(z) = (h/2) e^(-z eta) [A_0 + sum_j (A_j e^(i z x_j) + conj(A_j) e^(-i z x_j))] with the A_j of the contour Im u = eta."""

    def __init__(self, t, eta, h, xmax=2.6, nmax=40):
        self.t, self.eta, self.h = t, eta, h
        x = h * np.arange(0, int(round(xmax / h)) + 1)
        u = x + 1j * eta
        n = np.arange(1, nmax + 1, dtype=float)[:, None]
        w = np.exp(4 * u)[None, :]
        e9 = np.exp(9 * u)[None, :]
        e5 = np.exp(5 * u)[None, :]
        terms = (2 * np.pi ** 2 * n ** 4 * e9 - 3 * np.pi * n ** 2 * e5) * np.exp(-np.pi * n ** 2 * w)
        self.last = np.abs(terms[-1]).max()
        self.phi = terms.sum(axis=0)
        self.A = np.exp(t * u ** 2) * self.phi
        self.x = x
        self.xmax = xmax
        self.A[0] *= 0.5  # weight of j = 0 (counted twice below)

    def H(self, z):
        """z: complex or real 1d array. Returns H_t(z) up to the constant factor h/2 (positive)."""
        z = np.asarray(z, dtype=complex)
        out = np.empty(z.shape, dtype=complex)
        step = 1500
        for i in range(0, len(z), step):
            zz = z[i:i + step, None]
            e = np.exp(1j * zz * self.x[None, :])
            em = np.conj(e) if not np.any(zz.imag) else np.exp(-1j * zz * self.x[None, :])
            s = (e * self.A[None, :]).sum(axis=1) + (em * np.conj(self.A)[None, :]).sum(axis=1)
            out[i:i + step] = s * np.exp(-zz[:, 0] * self.eta) * (self.h / 2)
        return out


def find_zeros(ev, thmax, grid=0.04):
    th = np.arange(0.0, thmax + grid, grid)
    f = ev.H(2 * th).real
    s = np.sign(f)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    a, b = th[idx].copy(), th[idx + 1].copy()
    fa = f[idx]
    for _ in range(45):
        m = 0.5 * (a + b)
        fm = ev.H(2 * m).real
        left = np.sign(fm) == np.sign(fa)
        a = np.where(left, m, a)
        fa = np.where(left, fm, fa)
        b = np.where(left, b, m)
    return 0.5 * (a + b), f.max() / max(1e-300, 1), th


def arg_count(ev, Tp, y=0.5, grid=0.05):
    """(2 D_bottom + D_right)/(2 pi) for the half rectangle (0, Tp) x (-y, y) in theta."""
    sb = np.linspace(0.0, Tp, int(np.ceil(Tp / grid)) + 1)
    zb = 2 * (sb - 1j * y)
    hb = ev.H(zb)
    ab = np.angle(hb)
    ib = np.diff(ab)
    ib = (ib + np.pi) % (2 * np.pi) - np.pi
    sr = np.linspace(-y, y, int(round(2 * y / 0.01)) + 1)
    zr = 2 * (Tp + 1j * sr)
    hr = ev.H(zr)
    ar = np.angle(hr)
    ir = np.diff(ar)
    ir = (ir + np.pi) % (2 * np.pi) - np.pi
    # the start point has argument 0 (H real positive at theta = -0.5 i) and the bottom end joins the right edge
    # (the junction increment is the principal difference of the two arguments at the same point)
    junc = (ar[0] - ab[-1] + np.pi) % (2 * np.pi) - np.pi
    Db = ib.sum()
    Dr = junc + ir.sum()
    maxinc = max(np.abs(ib).max(), np.abs(ir).max())
    return (2 * Db + Dr) / (2 * np.pi), maxinc, abs(hb[0].imag / hb[0].real)


def run_t(t, nz, eta, h, gam, thmax_pad=1.0):
    ev = Evaluator(t, eta, h)
    thmax = gam[nz + 1] + thmax_pad
    zs, _, th = find_zeros(ev, thmax)
    return ev, zs


def stats(zs, nz):
    x = vartheta(zs[:nz]) / np.pi
    g = np.diff(x)
    return float(np.var(g)), float(g.min()), int(np.argmin(g)) + 1, float(g.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("zfile", nargs="?", default=None)
    ap.add_argument("--t", type=float, nargs="*", default=[0.0, 0.05, 0.1, 0.2])
    ap.add_argument("--nz", type=int, default=120)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    zfile = a.zfile or os.path.join(root, "data", "zeros_1700.txt")
    ts = [0.0, 0.2] if a.quick else a.t
    nz = 30 if a.quick else a.nz
    out = a.out or os.path.join(root, "data", "gp_flow_quick.json" if a.quick else "gp_flow.json")
    gam = []
    with open(zfile) as f:
        for line in f:
            if line.strip():
                gam.append(float(line.strip()[:40]))
            if len(gam) >= nz + 5:
                break
    gam = np.array(gam)
    res = {"nz": nz, "t": ts, "settings": {"eta1": 0.38, "h1": 0.00125, "eta2": 0.37, "h2": 0.0025}}
    t0 = time.time()
    for t in ts:
        ev1, z1 = run_t(t, nz, 0.38, 0.00125, gam)
        ev2, z2 = run_t(t, nz, 0.37, 0.0025, gam)
        n1, n2 = len(z1), len(z2)
        k = min(n1, n2)
        dmax = float(np.abs(z1[:k] - z2[:k]).max())
        log("t = %g: %d sign changes up to theta = %.3f (settings 1), %d (settings 2); max |difference of zeros| over %d = %.2e; last Phi term %.1e" % (t, n1, gam[nz + 1] + 1.0, n2, k, dmax, ev1.last))
        assert n1 >= nz + 1 and n1 == n2, "zero count mismatch"
        Tp = 0.5 * (z1[nz - 1] + z1[nz])
        c1, mi1, im1 = arg_count(ev2, Tp)
        c2, mi2, im2 = arg_count(ev2, Tp, y=0.25)
        log("t = %g: argument principle on (0, %.3f) x (-1/2, 1/2): %.6f (half height 1/2) and %.6f (half height 1/4) zeros (settings 2) against %d sign changes; largest argument increment %.3f, imaginary part at the corner %.1e" % (t, Tp, c1, c2, nz, max(mi1, mi2), im1))
        v, mn, jm, mean = stats(z1, nz)
        v2, mn2, jm2, _ = stats(z2, nz)
        log("t = %g: unfolded gaps of the first %d zeros (%d gaps): mean %.6f, variance %.10f, minimum %.10f at j = %d (second settings: %.10f, %.10f, j = %d; difference %.1e, %.1e)" % (t, nz, nz - 1, mean, v, mn, jm, v2, mn2, jm2, abs(v - v2), abs(mn - mn2)))
        entry = {"var": v, "min": mn, "jmin": jm, "mean": mean, "var2": v2, "min2": mn2, "argcount": [c1, c2], "maxdiff_zeros": dmax, "zeros": [float(x) for x in z1[:nz]]}
        if t == 0.0:
            d0 = np.abs(z1[:nz] - gam[:nz]).max()
            log("t = 0: max |theta_j - gamma_j| against data/zeros_1700.txt over j <= %d: %.2e" % (nz, d0))
            entry["max_diff_book"] = float(d0)
        for m in (80, 100):
            if m <= nz:
                vm, mnm, jmm, _ = stats(z1, m)
                log("t = %g, first %d zeros: variance %.10f, minimum %.10f at j = %d" % (t, m, vm, mnm, jmm))
                entry["first%d" % m] = [vm, mnm, jmm]
        res["t=%g" % t] = entry
    log("time %.1f s" % (time.time() - t0))
    json.dump(res, open(out, "w"), indent=1)
    open(out.replace(".json", ".log"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
