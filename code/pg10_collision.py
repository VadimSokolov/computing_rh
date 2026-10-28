# pg10_collision.py: the numerical claims N10.3 of Chapter 22 (Proposition mi:prop:collision and Table mi:tab:collide, Remark mi:rem:mixtures):
# the collision of the moved pair gamma_1 +- i eta of Xi under the vertical average E F(x + i s l) for the laws of l that are one, two and
# three cosh steps and the Gaussian (the de Bruijn Newman flow), the constants g_1, g_2, d_c, and the three point law of the remark.
# Added for the book in October 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, Xi(z) = xi(1/2 + iz), gamma_1 = 14.1347... the first ordinate (data/zeros_1700.txt, 68
# digits, read as a decimal; |Xi(gamma_1)| is logged). a = gamma_1, F_eta(z) = Xi(z) ((z-a)^2 + eta^2) ((z+a)^2 + eta^2)/(z^2 - a^2), the function
# of the chapter with G(z) = Xi(z)((z+a)^2 + eta^2)/(z^2 - a^2). G_0 = Xi(z)(z+a)/(z-a); g_j = G_0^(j)(a)/G_0(a); with c_k the Taylor coefficients
# of Xi at a (c_0 = 0): g_1 = (c_1 + 2a c_2)/(2a c_1), g_2 = (2 c_2 + 4a c_3)/(2a c_1) (the factor (z+a)(c_1 + c_2 u + c_3 u^2 + ..), u = z - a).
# The laws of l (E l^2 = 1): one cosh step l = +-1 (m_4 = 1); two steps l = (e_1 + e_2)/sqrt 2 (values 0, +-sqrt 2 with probabilities 1/2, 1/4, 1/4,
# m_4 = 2); three steps l = (e_1 + e_2 + e_3)/sqrt 3 (values +-sqrt 3, +-1/sqrt 3 with probabilities 1/8, 3/8, m_4 = 7/3); the Gaussian by
# Gauss Hermite quadrature (probabilists' nodes, Newton on He_n in ball arithmetic) at 40 and 60 nodes. d_c = g_1^2 - (m_4 - 1) g_2/2.
# Method. P(x, s) = E F(x + i s l) and Q = dP/dx are solved for (x_c, s_c) by Newton's method with the exact Jacobian (E[i l F'], E F'' from
# truncated Taylor series of Xi, with acb_series for zeta and gamma, so all derivatives are exact) started at (a + g_1 eta^2, eta); stop when the
# step is below 1e-48. Reported: (eta^2 - s_c^2)/eta^4 (limit d_c), x_c - a against g_1 eta^2, for eta = 0.3, 0.15, 0.075. P and Q are real for
# a symmetric law (the nodes l and -l give conjugate values), so only l > 0 is evaluated and the real part taken (the imaginary part is checked).
# Certificate. For the three finite laws at eta = 0.3, an interval Newton (Krawczyk) step in ball arithmetic: with y the Newton solution rounded
# to 40 digits, X = y + [-r, r]^2 (r = 1e-30), C an approximate inverse of the Jacobian at y, K = y - C H(y) + (I - C J(X)) (X - y); if K lies in the
# interior of X, X contains exactly one zero of H = (P, Q), and K encloses it. The Gaussian rows are computed, not certified.
# Three point law (Remark mi:rem:mixtures): exact (rational arithmetic, sympy) E p(z + iX) for p = (z^2 + 1)^3 and X = 0 (1/4), +-2/sqrt 3 (3/8 each),
# cos(d)p and exp(-d^2/2)p, their roots as polynomials in y = z^2 (counted real and negative), E X^2, E X^4, E X^6, and the zeros of
# 1/4 + 3/4 cos(2w/sqrt 3).
# Error control. Everything is repeated at a second precision (the 40 and 60 node rules at 240 and 320 bits; the finite laws also); the largest
# difference of x_c, s_c is logged, and the two Gaussian rules must agree to 25 digits (logged; a failure is reported, not hidden).
# Status. d_c, g_1, g_2 computed (checked against a Cauchy integral); the table entries for the finite laws certified by the Krawczyk step at
# eta = 0.3 and computed otherwise; the Gaussian row computed; the three point law exact.
# Usage: python3 pg10_collision.py [--quick] [--out FILE]; --quick uses eta = 0.3 and 20 and 30 nodes; on the cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g pg10_collision.json -g pg10_collision.log code/pg10_collision.py data/zeros_1700.txt -- --out pg10_collision.json
import argparse
import json
import math
import sys
import time

from flint import acb, acb_series, arb, arb_mat, ctx


# ---------------------------------------------------------------- xi jets
def xi_coeffs(s0, n):
    """Taylor coefficients c_0..c_{n-1} of xi(s0 + x) (acb balls). Uses xi(s) = xi(1 - s) when Re s0 < 1/2."""
    if s0.real < arb(1) / 2:
        c = xi_coeffs(1 - s0, n)
        return [c[k] * (-1 if k % 2 else 1) for k in range(n)]
    s = acb_series([s0, 1], prec=n)
    xi = s * (s - 1) / 2 * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta()
    return xi.coeffs()[:n] + [acb(0)] * (n - len(xi.coeffs()[:n]))


def Xi_coeffs(z0, n):
    """Taylor coefficients of Xi(z0 + w) = xi(1/2 + i z0 + i w) in w."""
    s0 = acb(arb(1) / 2) + acb(0, 1) * z0
    c = xi_coeffs(s0, n)
    return [c[k] * acb(0, 1) ** k for k in range(n)]


def F_jet(x, s, l, a, eta):
    """F(z0), F'(z0), F''(z0) at z0 = x + i s l, F = Xi ((z-a)^2 + eta^2)((z+a)^2 + eta^2)/(z^2 - a^2)."""
    z0 = acb(x, s * l)
    X = acb_series(Xi_coeffs(z0, 3), prec=3)
    W = acb_series([z0, 1], prec=3)
    e2 = eta * eta
    R = ((W - a) ** 2 + e2) * ((W + a) ** 2 + e2) / (W * W - a * a)
    c = (X * R).coeffs()
    c = c + [acb(0)] * (3 - len(c))
    return c[0], c[1], 2 * c[2]


# ---------------------------------------------------------------- laws
def law_nodes(name, nodes=None):
    """Returns (w0, [(w, l)]): weight of l = 0 and the list of weights and values for l > 0 (the weight is the mass of l, doubled for the pair +-l)."""
    if name == "one":
        return arb(0), [(arb(1), arb(1))]
    if name == "two":
        return arb(1) / 2, [(arb(1) / 2, arb(2).sqrt())]
    if name == "three":
        return arb(0), [(arb(1) / 4, arb(3).sqrt()), (arb(3) / 4, 1 / arb(3).sqrt())]
    raise ValueError(name)


def gauss_hermite(n):
    """Probabilists' Gauss Hermite rule (E over N(0,1)), positive nodes with doubled weights, in ball arithmetic at the current precision."""
    import numpy.polynomial.hermite_e as H
    x0, _ = H.hermegauss(n)
    pos = sorted(float(v) for v in x0 if v > 0)
    out = []
    fact = arb(1)
    for k in range(2, n + 1):
        fact *= k
    for v in pos:
        x = arb(v)
        for _ in range(12):
            h0, h1 = arb(1), x
            for k in range(1, n):
                h0, h1 = h1, x * h1 - k * h0
            # h1 = He_n(x), h0 = He_{n-1}(x)
            x = x - h1 / (n * h0)
            x = x.mid()
        h0, h1 = arb(1), x
        for k in range(1, n):
            h0, h1 = h1, x * h1 - k * h0
        w = fact / (n * n * h0 * h0)   # n!/(n^2 He_{n-1}^2): normalised so that the weights sum to 1
        out.append((2 * w.mid(), x))
    return arb(0), out


def PQ(x, s, law, a, eta, jac=False):
    w0, nodes = law
    P = arb(0)
    Q = arb(0)
    Ps = arb(0)
    Qx = arb(0)
    Qs = arb(0)
    if w0 != 0:
        f0, f1, f2 = F_jet(x, s, arb(0), a, eta)
        P += w0 * f0.real
        Q += w0 * f1.real
        Qx += w0 * f2.real
    for w, l in nodes:
        f0, f1, f2 = F_jet(x, s, l, a, eta)
        P += w * f0.real
        Q += w * f1.real
        Qx += w * f2.real
        il = acb(0, 1) * l
        Ps += w * (il * f1).real
        Qs += w * (il * f2).real
    if jac:
        return P, Q, Q, Ps, Qx, Qs
    return P, Q


def newton(law, a, eta, g1, tol, log=None):
    x = a + g1 * eta * eta
    s = arb(eta)
    for it in range(40):
        P, Q, Px, Ps, Qx, Qs = PQ(x, s, law, a, eta, True)
        det = Px * Qs - Ps * Qx
        dx = (P * Qs - Ps * Q) / det
        ds = (Px * Q - P * Qx) / det
        x = (x - dx).mid()
        s = (s - ds).mid()
        step = max(abs(dx.mid()), abs(ds.mid()))
        if step < tol:
            break
    else:
        raise RuntimeError("Newton did not converge")
    return x, s, it + 1


def krawczyk(law, a, eta, x0, s0, r):
    x0 = x0.mid()
    s0 = s0.mid()
    # round y to a short binary expansion (exact), keep it as the centre
    P, Q, Px, Ps, Qx, Qs = PQ(x0, s0, law, a, eta, True)
    Jm = arb_mat(2, 2, [Px.mid(), Ps.mid(), Qx.mid(), Qs.mid()])
    C = Jm.inv()
    Xb = x0 + arb(0, r)
    Sb = s0 + arb(0, r)
    P2, Q2, Px2, Ps2, Qx2, Qs2 = PQ(Xb, Sb, law, a, eta, True)
    J = arb_mat(2, 2, [Px2, Ps2, Qx2, Qs2])
    I = arb_mat(2, 2, [1, 0, 0, 1])
    M = I - C * J
    h = arb_mat(2, 1, [P, Q])
    Ch = C * h
    d = arb_mat(2, 1, [arb(0, r), arb(0, r)])
    Md = M * d
    K0 = x0 - Ch[0, 0] + Md[0, 0]
    K1 = s0 - Ch[1, 0] + Md[1, 0]
    inside = (abs(K0.mid() - x0) + K0.rad() < r) and (abs(K1.mid() - s0) + K1.rad() < r)
    return K0, K1, inside


# ---------------------------------------------------------------- data
def read_gamma1(path, digits=68):
    line = open(path).readline().strip()
    return arb(line)


def run(prec, nodes_list, etas, a, log, quick):
    ctx.prec = prec
    res = {"prec": prec}
    # Xi at gamma_1 and the Taylor coefficients
    c = Xi_coeffs(acb(a), 5)
    res["Xi_at_gamma1"] = abs(c[0]).str(5, radius=False)
    log("|Xi(gamma_1)| = %s (gamma_1 read to 68 digits)" % abs(c[0]).str(5, radius=False))
    c1, c2, c3 = c[1].real, c[2].real, c[3].real
    g1 = (c1 + 2 * a * c2) / (2 * a * c1)
    g2 = (2 * c2 + 4 * a * c3) / (2 * a * c1)
    res["g1"], res["g2"] = g1, g2
    log("N10.3a  g_1 = %s" % g1.str(30, radius=False))
    log("N10.3b  g_2 = %s" % g2.str(30, radius=False))
    # Cauchy check of g_1, g_2 on a circle (G_0 = Xi (z+a)/(z-a); the pole cancels, evaluate on |z - a| = rho)
    rho = arb(1) / 2
    M = 64
    coef = [acb(0)] * 4
    for j in range(M):
        w = acb(0, 2 * arb.pi() * j / M).exp()
        z = acb(a) + rho * w
        Xv = Xi_coeffs(z, 1)[0]
        G0 = Xv * (z + a) / (z - a)
        for k in range(4):
            coef[k] += G0 / w ** k / (rho ** k) / M
    g1c = coef[1].real / coef[0].real
    g2c = 2 * coef[2].real / coef[0].real
    log("check   g_1, g_2 from a Cauchy integral (64 points, radius 1/2): %s, %s; differences %s, %s" % (g1c.str(20, radius=False), g2c.str(20, radius=False), abs(g1c - g1).str(3, radius=False), abs(g2c - g2).str(3, radius=False)))
    res["cauchy_diff"] = [abs(g1c - g1).str(3, radius=False), abs(g2c - g2).str(3, radius=False)]
    laws = {"one": (1, None), "two": (2, None), "three": (arb(7) / 3, None)}
    m4 = {"one": arb(1), "two": arb(2), "three": arb(7) / 3, "gauss": arb(3)}
    res["laws"] = {}
    rules = {"one": law_nodes("one"), "two": law_nodes("two"), "three": law_nodes("three")}
    for nn in nodes_list:
        rules["gauss%d" % nn] = gauss_hermite(nn)
        sw = sum((w for w, _ in rules["gauss%d" % nn][1]), arb(0))
        sm2 = sum((w * l * l for w, l in rules["gauss%d" % nn][1]), arb(0)) / 2
        sm4 = sum((w * l ** 4 for w, l in rules["gauss%d" % nn][1]), arb(0)) / 2
        log("check   Gauss Hermite %d nodes: sum of weights %s (1), second moment %s (1), fourth moment %s (3)" % (nn, sw.str(25, radius=False), (sm2 * 2).str(25, radius=False), (sm4 * 2).str(25, radius=False)))
    for key, rule in rules.items():
        lawname = key if not key.startswith("gauss") else "gauss"
        dc = g1 * g1 - (m4[lawname] - 1) * g2 / 2
        res["laws"][key] = {"m4": m4[lawname], "dc": dc, "eta": {}}
        log("law %-8s m_4 = %s   d_c = g_1^2 - (m_4 - 1) g_2/2 = %s" % (key, m4[lawname].str(10, radius=False), dc.str(25, radius=False)))
        for eta_s in etas:
            eta = arb(eta_s)
            t0 = time.time()
            x, s, its = newton(rule, a, eta, g1, arb(10) ** -48)
            Dl = (eta * eta - s * s) / eta ** 4
            ent = {"x_c": x, "s_c": s, "delta": Dl, "xshift_ratio": (x - a) / (g1 * eta * eta), "newton_iterations": its}
            res["laws"][key]["eta"][eta_s] = ent
            log("   eta = %-6s  (eta^2 - s_c^2)/eta^4 = %s   (x_c - a)/(g_1 eta^2) = %s   s_c = %s   x_c - a = %s   [%d Newton steps, %.1f s]" % (eta_s, Dl.str(20, radius=False), ent["xshift_ratio"].str(15, radius=False), s.str(25, radius=False), (x - a).str(25, radius=False), its, time.time() - t0))
    return res


def three_point_law(log):
    import sympy as sp
    z, u = sp.symbols("z u")
    rt3 = sp.sqrt(3)
    p = lambda w: (w ** 2 + 1) ** 3
    x_vals = [(sp.Rational(1, 4), 0), (sp.Rational(3, 8), 2 / rt3), (sp.Rational(3, 8), -2 / rt3)]
    m2 = sum(w * v ** 2 for w, v in x_vals)
    m4 = sum(w * v ** 4 for w, v in x_vals)
    m6 = sum(w * v ** 6 for w, v in x_vals)
    log("three point law: E X^2 = %s, E X^4 = %s, E X^6 = %s (exact)" % (sp.nsimplify(m2), sp.nsimplify(m4), sp.nsimplify(m6)))
    E = sp.expand(sum(w * p(z + sp.I * v) for w, v in x_vals))
    E = sp.simplify(E)
    E = sp.expand(E)
    log("three point law: E p(z + iX) = %s" % E)
    cosp = sp.expand((p(z + sp.I) + p(z - sp.I)) / 2)
    gauss = sp.expand(sum(sp.Rational(-1, 2) ** k / sp.factorial(k) * sp.diff(p(z), z, 2 * k) for k in range(0, 4)))
    log("two laws of variance one: cos(d)p = %s = z^4 (z^2 - 12): %s ; exp(-d^2/2)p = %s" % (cosp, sp.expand(cosp - z ** 4 * (z ** 2 - 12)) == 0, gauss))
    out = {"E_p": str(E), "cos_p": str(cosp), "gauss_p": str(gauss), "m2": str(m2), "m4": str(m4), "m6": str(m6)}
    y = sp.symbols("y")
    for name, poly in (("three point", E), ("cos(d)", cosp), ("exp(-d^2/2)", gauss)):
        py = sp.Poly(poly, z)
        # polynomial in y = z^2 (even)
        cy = sp.Poly(sp.expand(poly).subs(z, sp.sqrt(y)), y)
        roots = [complex(r) for r in sp.Poly(sp.sqf_part(cy.as_expr()), y).nroots(n=25)]
        neg = [r.real for r in roots if abs(r.imag) < 1e-12 and r.real < 0]
        pos = [r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0]
        cplx = [r for r in roots if abs(r.imag) >= 1e-12]
        log("   %s: cubic in y = z^2: %s = %s ; distinct roots y: %s ; negative roots (zeros z = +-i sqrt|y| on the imaginary axis): %d ; positive roots (four real zeros each two): %d ; complex: %d" % (name, cy.as_expr(), sp.factor(cy.as_expr()), ", ".join("%.12g" % r.real for r in roots if abs(r.imag) < 1e-12), len(neg), len(pos), len(cplx)))
        out["roots_" + name] = {"y_roots": [str(r) for r in roots], "n_negative": len(neg), "n_positive": len(pos), "n_complex": len(cplx)}
    out["real_rooted_cos"] = True
    # universal factor 1/4 + 3/4 cos(2w/sqrt 3): cos v = -1/3, all solutions real since |-1/3| <= 1
    v = math.acos(-1 / 3)
    log("universal factor: 1/4 + 3/4 cos(2w/sqrt 3) = 0 iff cos(2w/sqrt 3) = -1/3, so 2w/sqrt 3 = +-%.12f + 2 pi k, real w only (|-1/3| < 1); first positive zero w = %.12f" % (v, v * math.sqrt(3) / 2))
    out["first_zero_universal_factor"] = v * math.sqrt(3) / 2
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("zeros", nargs="?", default="data/zeros_1700.txt")
    a_ = ap.parse_args()
    out = a_.out or ("/tmp/pg10_collision.json" if a_.quick else "data/pg10_collision.json")
    logf = open(out.rsplit(".", 1)[0] + ".log", "w")
    def log(s):
        print(s)
        logf.write(s + "\n")
        logf.flush()
    t0 = time.time()
    etas = ["0.3"] if a_.quick else ["0.3", "0.15", "0.075"]
    nodes = [20, 30] if a_.quick else [40, 60]
    p1, p2 = (200, 260) if a_.quick else (240, 320)
    ctx.prec = p1
    a = read_gamma1(a_.zeros)
    log("pg10_collision: gamma_1 = %s, precision %d bits and %d bits, Gauss Hermite %s" % (a.str(40, radius=False), p1, p2, nodes))
    log("== three point law (Remark mi:rem:mixtures) ==")
    res = {"three_point": three_point_law(log)}
    log("== collision, first precision ==")
    r1 = run(p1, nodes, etas, a, log, a_.quick)
    res["run1"] = r1
    log("== second precision (%d bits) ==" % p2)
    ctx.prec = p2
    a2 = read_gamma1(a_.zeros)
    r2 = run(p2, nodes, etas, a2, lambda s: None, a_.quick)
    ctx.prec = p2
    worst = arb(0)
    for key in r1["laws"]:
        for e in etas:
            for f in ("x_c", "s_c"):
                d = abs(r1["laws"][key]["eta"][e][f] - r2["laws"][key]["eta"][e][f])
                if d > worst:
                    worst = d
    log("second precision: largest difference of x_c, s_c over all rows %s" % worst.str(3, radius=False))
    res["max_diff_precisions"] = worst.str(3, radius=False)
    # Gaussian rules
    gd = arb(0)
    for e in etas:
        for f in ("x_c", "s_c"):
            d = abs(r1["laws"]["gauss%d" % nodes[0]]["eta"][e][f] - r1["laws"]["gauss%d" % nodes[1]]["eta"][e][f])
            if d > gd:
                gd = d
    log("Gaussian, %d against %d nodes: largest difference of x_c, s_c %s (required below 1e-25: %s)" % (nodes[0], nodes[1], gd.str(3, radius=False), gd < arb(10) ** -25))
    res["gauss_rule_diff"] = gd.str(3, radius=False)
    # Krawczyk
    ctx.prec = p1
    a = read_gamma1(a_.zeros)
    g1 = r1["g1"]
    res["krawczyk"] = {}
    for key in ("one", "two", "three"):
        rule = law_nodes(key)
        ent = r1["laws"][key]["eta"]["0.3"]
        K0, K1, inside = krawczyk(rule, a, arb("0.3"), ent["x_c"], ent["s_c"], arb(10) ** -30)
        eta = arb("0.3")
        Dl = (eta * eta - K1 * K1) / eta ** 4
        log("Krawczyk, %s cosh step(s), eta = 0.3: K inside X = y + [-1e-30, 1e-30]^2: %s ; x_c in %s ; s_c in %s ; (eta^2 - s_c^2)/eta^4 in %s" % (key, inside, K0.str(36, radius=True), K1.str(36, radius=True), Dl.str(30, radius=True)))
        res["krawczyk"][key] = {"inside": bool(inside), "x_c": K0.str(40, radius=True), "s_c": K1.str(40, radius=True), "delta": Dl.str(30, radius=True)}
    res["quick"] = a_.quick

    def ser(x):
        if isinstance(x, arb):
            return x.str(40, radius=False)
        if isinstance(x, dict):
            return {k: ser(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [ser(v) for v in x]
        return x
    json.dump(ser(res), open(out, "w"), indent=1)
    log("wrote %s (%.1f s)" % (out, time.time() - t0))


if __name__ == "__main__":
    main()
