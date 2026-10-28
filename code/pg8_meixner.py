# pg8_meixner.py: the real and nonreal zeros of q_phi * Xi for phi = 0, pi/9, pi/4, pi/2 (Chapter 22, Section mi:sec:mixture, Proposition
# mi:prop:meixner and the three computations of the last paragraph; placeholders [N8.3a] to [N8.3e] of reconcile-X8.md). Added for the book in
# October 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Object. F_phi(z) = int Phi(u) (1 + cos phi)/(cosh u + cos phi) e^{izu} du = (q_phi * Xi)(z) = E H_{-varsigma_phi}(z) (Proposition
# mi:prop:meixner (ii)); phi = pi/2 is g_1 of Theorem mi:thm:levyarea, whose 51 real zeros in (0, 250) (mi4.tex, "51 and 51 at b = 1") and the
# 14 real zeros of phi = pi/9 below 65 are the required cross checks. The count for Xi itself is 108 (mi4.tex).
# Method. The contour method of pg8_constants.py (line Im u = 0.77, scaled by e^{0.77 z}; H = 1000 nodes per unit, v <= 4.2, 160 bits): real
# zeros by sign changes of F_s(t) on a grid of step 0.01 (double precision) refined by Newton in ball arithmetic and certified as sign
# changes of F_s at root -+ 1e-12; the argument principle on [0.5, x1] x [-3, 3] (x1 = 250 moved off the zeros) with a Lipschitz test on every
# step (track() in pg8_constants.py); the nonreal zeros from grid minima of |F_s| plus Newton; the first nonreal zero for phi = pi/9 enclosed by
# interval Newton (Krawczyk) with an allowance for the quadrature error. Second quadrature: H = 1250, v <= 4.6 (a different step and cut), the
# difference is logged. Status: the ball arithmetic bounds the rounding error; the quadrature error is estimated by the second run and not
# bounded, so every statement is "computed", and "certified" only for the two things the ball arithmetic does prove for the discretised F_s
# (nonvanishing on the boundary and the winding number, the sign changes, the enclosure of the zero).
# Output: data/pg8_meixner.json and data/pg8_meixner.log. Usage: python3 pg8_meixner.py [--quick] [--out FILE] [--cpus N] [--tmax T].
# --quick (H = 250, v <= 3.6, 64 bits, scan and counts only, the argument principle on [0.5, 30]) takes 20 s on the laptop.
# Hopper: bash misc/tools/hopper_run.sh -c 8 -t 240 -n pg8_meixner -g pg8_meixner.json -g pg8_meixner.log code/pg8_meixner.py code/pg8_constants.py data/zeros_1700.txt -- --out pg8_meixner.json
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
import numpy as np
import pg8_constants as C
from flint import acb, arb, ctx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--H", type=int, default=None, help="nodes per unit (default 1000; for tests only)")
    ap.add_argument("--V", type=float, default=None, help="cut in v (default 4.2)")
    ap.add_argument("--prec", type=int, default=None, help="bits (default 160)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--tmax", type=float, default=250.0)
    ap.add_argument("--zeros", default=None)
    ap.add_argument("--only", default=None, help="comma separated list of phi names, e.g. pi/9")
    ap.add_argument("--x1", type=float, default=None, help="right edge of the argument principle box (quick default 30)")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or (os.path.join(here, "..", "data", "pg8_meixner.json") if not a.quick else "/tmp/pg8_meixner.json")
    log = C.Log(os.path.splitext(out)[0] + ".log")
    if a.quick:
        H, V, H2, V2, prec, box_x1 = 250, 3.6, 300, 3.9, 64, (a.x1 or 30.0)
    else:
        H, V, prec = (a.H or 1000), (a.V or 4.2), (a.prec or 160)
        H2, V2 = int(H * 1.25), V + 0.4
        box_x1 = a.tmax
    ctx.prec = prec
    pi = arb.pi()
    phis = [("0", arb(0)), ("pi/9", pi / 9), ("pi/4", pi / 4), ("pi/2", pi / 2)]
    zpath = a.zeros or next((p for p in [os.path.join(here, "..", "data", "zeros_1700.txt"), "zeros_1700.txt"] if os.path.exists(p)), None)
    xi_zeros = [float(z) for z in C.load_ordinates(zpath)] if zpath else []
    results = {"parameters": {"H": H, "V": V, "H2": H2, "V2": V2, "prec": prec, "tmax": a.tmax, "tau": 0.77, "quick": a.quick}}
    t0 = time.time()
    for name, phi in phis:
        if a.only and name not in a.only.split(","):
            continue
        qa = C.Quad(C.kern_meixner(phi), H, V)
        qb = None if a.quick else C.Quad(C.kern_meixner(phi), H2, V2)
        log("== phi = %s: nodes built (%.1f s)" % (name, time.time() - t0))
        r = C.analyse("phi=" + name, qa, qb, a.tmax, (0.5, box_x1, 3.0), a.cpus, log, quick=a.quick)
        rz = [float(x) for x in r["roots"]]
        n_all = sum(1 for x in rz if x <= a.tmax)
        r["real_zeros_0_tmax"] = n_all
        r["real_zeros_below_65"] = sum(1 for x in rz if x <= 65)
        # first deficit of two against Xi
        if rz and xi_zeros:
            first = None
            for k, g in enumerate(xi_zeros, start=1):
                if g > a.tmax:
                    break
                nphi = sum(1 for x in rz if x < g)
                if k - nphi >= 2:
                    first = (k, g)
                    break
            r["first_deficit_of_two_at_xi_zero"] = None if first is None else {"index": first[0], "ordinate": "%.10g" % first[1]}
            log("[phi=%s] real zeros in (0,%g]: %d; below 65: %d; Xi has %d below %g; first Xi zero at which the deficit reaches 2: %s"
                % (name, a.tmax, n_all, r["real_zeros_below_65"], sum(1 for g in xi_zeros if g <= a.tmax), a.tmax, first))
        if r.get("nonreal_zeros_upper"):
            zs = r["nonreal_zeros_upper"]
            r["smallest_real_part_of_a_nonreal_zero"] = zs[0][0]
            log("[phi=%s] all zeros with real part below %s are real (smallest real part of a nonreal zero); first nonreal zero %s + %s i" % (name, zs[0][0], zs[0][0], zs[0][1]))
            if name == "pi/9" and not a.quick:
                z1 = acb(arb(zs[0][0]), arb(zs[0][1]))
                z1 = C.newton_c(qa, z1, tol=1e-34)
                # quadrature allowance: difference of the two quadratures at z1
                qerr = float(abs(qa.eval(z1) - qb.eval(z1)).mid()) * 10 + 1e-40
                ok, rad, Nx = C.enclose(qa, z1, qerr)
                r["first_nonreal_zero"] = {"re": C.S(z1.real, 30), "im": C.S(z1.imag, 30), "enclosure_radius": "%g" % rad, "krawczyk_contains": ok, "quadrature_allowance": "%.3g" % qerr}
                log("[phi=pi/9] first nonreal zero z = %s + %s i; interval Newton in the box +-%g (allowance %.3g for the quadrature error): %s"
                    % (C.S(z1.real, 30), C.S(z1.imag, 30), rad, qerr, "encloses exactly one simple zero" if ok else "FAILED"))
        results["phi=" + name] = r
        with open(out, "w") as fh:
            json.dump(results, fh, indent=1)
    if "phi=pi/2" in results and "phi=pi/9" in results:
        log("cross checks: phi = pi/2 real zeros in (0,%g]: %s (required 51); phi = pi/9 below 65: %s (required 14); total time %.1f s"
            % (a.tmax, results["phi=pi/2"]["real_zeros_0_tmax"], results["phi=pi/9"]["real_zeros_below_65"], time.time() - t0))
    with open(out, "w") as fh:
        json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
