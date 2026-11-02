"""Independent computations for report.md. No input from the manuscript's code.

Run: PYTHONPATH=./pylib python3 referee_checks.py [core|density|zeros]
The core and density modes use Arb complex balls through python-flint.
The zeros mode separates finite ball sums from an unbounded model tail.
"""

from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import flint
from flint import acb, acb_series, arb, ctx
import mpmath as mp

ROOT = Path(__file__).resolve().parent
ctx.prec = 256
ctx.cap = 20
mp.mp.dps = 85


def xi_series(s, order=8):
    """Removable singularity at s=1 is evaluated using deflated zeta."""
    q = acb_series([s, 1], prec=order)
    return (q / 2 * (-q * arb.pi().log() / 2).exp()
            * (q / 2).gamma() * (1 + (q - 1) * q.zeta(deflate=True)))


def xi(s):
    return xi_series(s, 1)[0]


def tau_ball():
    return ((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)


def dh_series(s, order=8, completed=True):
    q = acb_series([s, 1], prec=order)
    tau = tau_ball()
    coeffs = [1, tau, -tau, -1]
    hurwitz = sum(c * q.zeta(a=arb(k) / 5, deflate=True)
                  for k, c in enumerate(coeffs, 1))
    d = (-q * arb(5).log()).exp() * hurwitz
    if completed:
        d *= (((q + 1) / 2) * (arb(5) / arb.pi()).log()).exp()
        d *= ((q + 1) / 2).gamma()
    return d


def dh(s, completed=True):
    return dh_series(s, 1, completed)[0]


def mp_xi(s):
    if s == 0 or s == 1:
        return mp.mpf("0.5")
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def mp_dh(s):
    tau = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    return mp.fsum(c * mp.zeta(s, mp.mpf(k) / 5)
                   for k, c in enumerate([1, tau, -tau, -1], 1)) / 5 ** s


def text_ball(x, digits=38):
    return x.str(digits)


def metadata():
    return {"python": platform.python_version(), "python_flint": flint.__version__,
            "mpmath": mp.__version__, "precision_bits": ctx.prec,
            "note_sha256": hashlib.sha256((ROOT / "note.tex").read_bytes()).hexdigest()}


def save(name, data):
    (ROOT / "results").mkdir(exist_ok=True)
    path = ROOT / "results" / (name + ".json")
    path.write_text(json.dumps(data, indent=2) + "\n")
    print(path.name, flush=True)


def core():
    out = {"metadata": metadata(), "zeta": {}, "dh": {}}
    z = out["zeta"]
    z["xi_half"] = text_ball(xi(acb(arb(1) / 2)))
    logs = xi_series(acb(arb(1) / 2), 18).log()
    z["cumulants_X"] = {str(k): text_ball(logs[k] * math.factorial(k))
                          for k in range(2, 17, 2)}
    z["cumulants_S"] = {str(k): text_ball(logs[k] * math.factorial(k) * arb(k).zeta())
                          for k in (2, 4)}
    z["magnetisations"] = {}
    for a in ("1", "1.5", "2", "3"):
        value = xi_series(acb(a), 3).log()[1]
        assert value.real > 0 and value.imag.contains(0)
        z["magnetisations"][a] = text_ball(value)
    first = acb.zeta_zero(1)
    z["zero_1"] = text_ball(first)
    z["N_14"] = text_ball(arb(14).zeta_nzeros())
    z["N_15"] = text_ball(arb(15).zeta_nzeros())
    z["N_45"] = text_ball(arb(45).zeta_nzeros())
    assert arb(14).zeta_nzeros() == 0 and arb(15).zeta_nzeros() == 1
    z["first_nine_ordinates"] = [text_ball(r.imag) for r in acb.zeta_zeros(1, 9)]
    z["multiples"] = {}
    for k in (2, 3):
        v = xi(acb(arb(1) / 2, k * first.imag))
        assert not v.contains(0)
        z["multiples"][str(k)] = {"ordinate": text_ball(k * first.imag),
                                  "xi_value": text_ball(v),
                                  "abs_lower": text_ball(abs(v).lower())}
    z["mpmath_crosscheck"] = {
        "cumulant2_S": mp.nstr(mp.diff(lambda x: mp.log(mp_xi(x)), mp.mpf("0.5"), 2) * mp.zeta(2), 45),
        "cumulant4_S": mp.nstr(mp.diff(lambda x: mp.log(mp_xi(x)), mp.mpf("0.5"), 4) * mp.zeta(4), 45)}
    z["zeta2"] = text_ball(arb(2).zeta())
    z["zeta4"] = text_ball(arb(4).zeta())
    z["variance_fraction"] = text_ball(1 - 6 / arb.pi() ** 2)
    z["note_corrected_variance_minus_derivative"] = text_ball(
        arb("0.0760123845") - (logs[2] * 2 * arb(2).zeta()).real)

    d = out["dh"]
    d["tau"] = text_ball(tau_ball())
    root = mp.findroot(mp_dh, (mp.mpc("0.81", "85.70"), mp.mpc("0.80", "85.69")),
                       tol=mp.mpf("1e-80"), maxsteps=30)
    center_re, center_im = mp.nstr(root.real, 70), mp.nstr(root.imag, 70)
    center = acb(center_re, center_im)
    radius = arb("1e-25")
    box = acb(arb(center_re, "1e-25"), arb(center_im, "1e-25"))
    at_center = dh_series(center, 3, completed=False)
    at_box = dh_series(box, 3, completed=False)
    residual = abs(at_center[0]).upper()
    derivative_lower = abs(at_center[1]).lower()
    half_second_upper = abs(at_box[2]).upper()
    left = residual + half_second_upper * radius ** 2
    right = derivative_lower * radius
    assert left < right
    assert box.real > arb("0.8")
    d["root_certificate"] = {
        "center_real": center_re, "center_imag": center_im,
        "radius": text_ball(radius), "function": "uncompleted D(s)",
        "residual_abs_upper": text_ball(residual),
        "derivative_abs_lower": text_ball(derivative_lower),
        "half_second_abs_upper_on_box": text_ball(half_second_upper),
        "rouche_left_upper": text_ball(left.upper()),
        "rouche_right_lower": text_ball(right.lower()),
        "proved_inequality": "|D(c)| + sup|D''| r^2/2 < |D'(c)| r"}
    d["functional_equation_checks"] = {}
    for re, im in (("0.7", "12.3"), ("2", "3"), ("0.5", "85.7")):
        s = acb(re, im)
        difference = dh(s) - dh(1 - s)
        assert difference.contains(0)
        d["functional_equation_checks"][re + "+i" + im] = text_ball(difference)
    dh_logs = dh_series(acb(arb(1) / 2), 18).log()
    d["completed_value_half"] = text_ball(dh(acb(arb(1) / 2)))
    d["cumulants_X"] = {str(k): text_ball(dh_logs[k] * math.factorial(k))
                          for k in range(2, 17, 2)}
    for k in range(2, 17, 2):
        assert ((-1) ** (k // 2 - 1) * dh_logs[k]).real > 0
    d["cumulants_S"] = {str(k): text_ball(dh_logs[k] * math.factorial(k) * arb(k).zeta())
                          for k in (2, 4)}
    for a in ("1", "2"):
        d["magnetisation_" + a] = text_ball(dh_series(acb(a), 3).log()[1])
    d["zero_free_majorant_at_2"] = text_ball(arb(2).zeta() - 1)
    # A density boundary value just to the left of the off-line zero.
    s = acb("0.8", center_im)
    value = dh_series(s, 2).log()[1].real / arb.pi()
    assert value < 0
    d["boundary_Re_logderivative_at_0.8_gamma_over_pi"] = text_ball(value)
    d["elementary_zero_width_at_alpha2"] = text_ball(2 - box.real)
    save("core", out)


def density():
    """Certify f_W(x) using a contour shifted vertically by 1.

    Integral over [0,1]: -Im integral exp(-xt) xi(a)/xi(a+it) dt / pi.
    The omitted horizontal integral has modulus <=54 exp(-x)/pi.
    The proof of this explicit bound is supplied in report.md.
    """
    ctx.prec = 192
    first = acb.zeta_zero(1)
    assert arb(14).zeta_nzeros() == 0
    assert first.imag > 14 and first.imag < 15
    assert xi_series(acb(1), 2).log()[1].real < arb("0.024")
    out = {"metadata": metadata(), "method": "vertical contour plus rigorous horizontal remainder",
           "rows": {}}
    for alpha in ("1", "1.5", "2", "3"):
        a = acb(alpha)
        xa = xi(a)
        row = {"magnetisation": text_ball(xi_series(a, 2).log()[1])}
        for x in (50, 100, 200):
            def integrand(t, analytic):
                return (-x * t).exp() * xa / xi(a + acb(0, 1) * t)
            integral = acb.integral(integrand, 0, 1,
                                    rel_tol=arb("1e-30"), abs_tol=arb("1e-35"),
                                    eval_limit=20000, depth_limit=30)
            error = 54 * arb(-x).exp() / arb.pi()
            f = -integral.imag / arb.pi() + arb(0, error.upper())
            assert f > 0 and f.is_finite()
            row[str(x)] = {"density": text_ball(f),
                           "pi_x2_density": text_ball(arb.pi() * x * x * f),
                           "horizontal_error_upper": text_ball(error.upper())}
            print(alpha, x, text_ball(arb.pi() * x * x * f, 20), flush=True)
        out["rows"][alpha] = row
        save("density", out)


def zeros():
    """The finite 2000-zero sum is a ball enclosure; the density tail is a model."""
    ctx.prec = 100
    out = {"metadata": metadata(), "N": 2000}
    s2 = arb(0)
    s4 = arb(0)
    last = None
    for start in range(1, 2001, 100):
        block = acb.zeta_zeros(start, 100)
        for r in block:
            assert r.real == arb(1) / 2
            t = r.imag
            s2 += 1 / t ** 2
            s4 += 1 / t ** 4
            last = t
        print("zeros through", start + 99, flush=True)
    # Integrate the smooth counting density (2 pi)^(-1) log(t/(2 pi)).
    tail2 = (last / (2 * arb.pi())).log() / last + 1 / last
    tail2 /= 2 * arb.pi()
    tail4 = ((last / (2 * arb.pi())).log() + arb(1) / 3) / (6 * arb.pi() * last ** 3)
    out.update({"last_ordinate": text_ball(last),
                "N_2515": text_ball(arb(2515).zeta_nzeros()),
                "N_2516": text_ball(arb(2516).zeta_nzeros()),
                "sum_gamma_inverse2": text_ball(s2), "sum_gamma_inverse4": text_ball(s4),
                "finite_uniform_variance": text_ball(arb.pi() ** 2 / 3 * s2),
                "finite_uniform_cumulant4": text_ball(-2 * arb.pi() ** 4 / 15 * s4),
                "model_tail_inverse2": text_ball(tail2), "model_tail_inverse4": text_ball(tail4),
                "model_corrected_variance": text_ball(arb.pi() ** 2 / 3 * (s2 + tail2)),
                "model_corrected_cumulant4": text_ball(-2 * arb.pi() ** 4 / 15 * (s4 + tail4)),
                "tail_status": "heuristic: no bound on omitted zero-counting error"})
    assert arb(2515).zeta_nzeros() == 1999 and arb(2516).zeta_nzeros() == 2000
    save("zeros", out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["core", "density", "zeros"], default="core", nargs="?")
    args = parser.parse_args()
    started = time.monotonic()
    globals()[args.mode]()
    print("elapsed_seconds", round(time.monotonic() - started, 3), flush=True)
