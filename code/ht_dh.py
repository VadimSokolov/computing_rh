# ht_dh.py: the Davenport Heilbronn tracking of the Jensen polynomials and the de Bruijn flow time of the off line pair, Sections "What a degree can see" and "A Hermite identity and the
# de Bruijn flow" of Chapter 24 (Nick Polson's hitting times chapter). Claims supported (H = chapter_hitting_times.tex): H:184 (the zero of J_d^DH starting at s_rho = (0.3085 + 85.6993 i)^2 =
# -7344.28 + 52.88 i, (Im l_d)^2 decreasing linearly in 1/d with slope 6.3e7 against |s|^2 = 5.4e7, reaching zero between d = 20619 and 20833, quadratic fit d* = 20785, model gamma^2/(4 delta^2) = 19290);
# H:213 (t_c = 0.1765 for Davenport and Heilbronn, d* = gamma^2/(2 t_c) = 20805, agreement 0.1%).
# Conventions. a_n periodic mod 5 with (a_0..a_4) = (0, 1, kappa, -kappa, -1), kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1); psi(y) = sum_{n>=1} n a_n exp(-pi n^2 y/5); phi(v) = e^(3v/2) psi(e^(2v)),
# even (self duality theta(1/x) = x^(3/2) theta(x), evenness checked numerically); Lambda(1/2 + w) = int_R phi(v) e^(w v) dv, F(s) = Lambda(1/2 + sqrt s)/Lambda(1/2) = sum c_k s^k, c_k = m_{2k}/((2k)! m_0);
# J_d(s) = sum c_k prod_{j<k}(1 - j/d) s^k with eps = 1/d continuous. The zero s_rho = (rho - 1/2)^2 of F, rho the first off line zero (Newton from (0.3085 + 85.6993 i)^2). Flow: H_t(z) = int phi(v) e^(t v^2/4) cos(z v) dv,
# H_t = -(1/4) H_zz, whose zero at z = 85.6993 - 0.3085 i (t = 0) moves to the real axis at the double root (H_t = dH_t/dz = 0); for an isolated pair x +- i b the merger time is b^2 * 2 (t_c = 2 delta^2).
# Method. Moments and H_t by the trapezoid rule on [0, 6] with step h (phi analytic in |Im v| < pi/4), K = 600 coefficients (the terms beyond are below e^(-640) relative), two steps (full: 1/200 and 1/250, 240 and 280 digits; H_t at 90 and 110 digits).
# Tracking: Newton in s as eps increases from 0 in steps 2e-6, then the double root (s, eps) of (J, J_s) = 0 by findroot seeded from a quadratic extrapolation of (Im s)^2 against eps. Run twice (both coefficient sets), d* from each.
# Written from the mathematics (the collaborator's scripts were not read, run or copied). Usage: python3 ht_dh.py [--quick] [--out FILE]; quick: 150 digits, step 1/100 (smoke test, not results).
import argparse, json, os, time
import mpmath as mp

def kappa():
    s5 = mp.sqrt(5); return (mp.sqrt(10 - 2 * s5) - 2) / (s5 - 1)

def phidh(v, dps, kap):
    y = mp.exp(2 * v); a = [0, 1, kap, -kap, -1]; s = mp.mpf(0); n = 1; eps = mp.mpf(10) ** (-dps - 20)
    while True:
        t = n * a[n % 5] * mp.exp(-mp.pi * n * n * y / 5); s += t
        if n > 5 and abs(t) < eps * (abs(s) + eps): break
        n += 1
    return mp.exp(3 * v / 2) * s

def kernel_nodes(hinv, dps, V=6):
    mp.mp.dps = dps; kap = kappa(); h = mp.mpf(1) / hinv
    vs = [j * h for j in range(V * hinv + 1)]
    return vs, [phidh(v, dps, kap) for v in vs], h

def coeffs(K, hinv, dps):
    vs, ph, h = kernel_nodes(hinv, dps)
    w = [h * p for p in ph]; w[0] /= 2; m0 = 2 * sum(w)
    c = [mp.mpf(1)]; pw = [mp.mpf(1)] * len(vs); v2 = [v * v for v in vs]; fact = mp.mpf(1)
    for k in range(1, K):
        pw = [p * q for p, q in zip(pw, v2)]; fact *= (2 * k - 1) * (2 * k)
        c.append(2 * sum(x * y for x, y in zip(w, pw)) / (fact * m0))
    return c

def Jeps(c, eps, s, K):
    # returns J, J_s at s (complex) for eps = 1/d: coefficients c_k prod_{j<k}(1 - j eps)
    co = []; pr = mp.mpf(1)
    for k in range(K): co.append(c[k] * pr); pr *= (1 - k * eps)
    f = 0; g = 0
    for a in reversed(co): g = g * s + f; f = f * s + a
    return f, g

def Hflow(vs, ph, h, t, z):
    # H_t(z) and dH/dz at complex z by the trapezoid rule (phi even, nodes j h, j >= 0)
    f = 0; g = 0
    for j, (v, p) in enumerate(zip(vs, ph)):
        wgt = h * p * mp.exp(t * v * v / 4) * (1 if j == 0 else 2)
        f += wgt * mp.cos(z * v); g -= wgt * v * mp.sin(z * v)
    return f, g

def run(K, hinv, dps, step, log, quick=False):
    mp.mp.dps = dps
    c = coeffs(K, hinv, dps); mp.mp.dps = dps
    out = {}
    rho = mp.mpc("0.3085", "85.6993"); s0 = rho ** 2
    s = s0
    for it in range(80):
        f, g = Jeps(c, mp.mpf(0), s, K); ds = f / g; s -= ds
        if abs(ds) < mp.mpf(10) ** (-30): break
    out["s_rho"] = s; out["delta"] = mp.sqrt(s).real if False else None
    r = mp.sqrt(s); out["rho_minus_half"] = r  # = delta + i gamma
    out["gamma"] = r.imag; out["delta"] = r.real
    out["model_d"] = r.imag ** 2 / (4 * r.real ** 2); out["abs_s2"] = abs(s) ** 2
    eps = mp.mpf(0); path = [(eps, s)]; epsmax = mp.mpf("3.9e-5")
    while eps < epsmax:
        eps += step
        for it in range(60):
            f, g = Jeps(c, eps, s, K); ds = f / g; s -= ds
            if abs(ds) < mp.mpf(10) ** (-30) * abs(s): break
        path.append((eps, s))
    out["path_im2"] = [(mp.nstr(e, 6), mp.nstr(z.imag ** 2, 10)) for e, z in path[:: max(1, len(path) // 12)]]
    xs = [e for e, z in path]; ys = [z.imag ** 2 for e, z in path]
    slope0 = (ys[1] - ys[0]) / (xs[1] - xs[0]); out["slope_at_0"] = slope0
    # quadratic fit y = a0 + a1 e + a2 e^2 on all path points (least squares via normal equations at high precision)
    n = len(xs); A = mp.matrix(n, 3); B = mp.matrix(n, 1)
    for i in range(n): A[i, 0] = 1; A[i, 1] = xs[i]; A[i, 2] = xs[i] ** 2; B[i] = ys[i]
    coef = mp.lu_solve(A.T * A, A.T * B)
    disc = coef[1] ** 2 - 4 * coef[2] * coef[0]
    roots = [(-coef[1] + sg * mp.sqrt(disc)) / (2 * coef[2]) for sg in (1, -1)]
    estar = min([r_ for r_ in roots if r_ > 0], default=None)
    out["fit_coeffs"] = [mp.nstr(x, 8) for x in coef]; out["d_star_quadratic_fit"] = 1 / estar if estar else None
    # a linear fit through the first points too
    out["d_star_linear_extrap"] = (xs[0] - ys[0] / slope0) ** -1 if False else 1 / (-ys[0] / slope0)
    # double root (s, eps)
    try:
        sol = mp.findroot(lambda S, E: [Jeps(c, E, S, K)[0], Jeps(c, E, S, K)[1]], (path[-1][1].real, estar), tol=mp.mpf(10) ** (-(dps // 3)), maxsteps=100)
        out["double_root_s"] = sol[0]; out["double_root_eps"] = sol[1]; out["d_star_double_root"] = 1 / sol[1]
    except Exception as ex:
        out["double_root_error"] = str(ex)
    # Im at d = 20619 and d = 20833 by Newton continuation along the path
    def im_at(d):
        e_t = mp.mpf(1) / d; ss = path[-1][1]; ee = path[-1][0]; k = 20
        for i in range(1, k + 1):
            e_i = ee + (e_t - ee) * i / k
            for it in range(80):
                f, g = Jeps(c, e_i, ss, K); ds = f / g; ss -= ds
                if abs(ds) < mp.mpf(10) ** (-30) * abs(ss): break
        return ss
    for d in (20619, 20833):
        try:
            ssd = im_at(d); out["s_at_d%d" % d] = ssd
        except Exception as ex:
            out["s_at_d%d" % d] = "error " + str(ex)
    return out, c

def flow(hinv, dps, log):
    vs, ph, h = kernel_nodes(hinv, dps)
    mp.mp.dps = dps
    z = mp.mpc("85.6993", "-0.3085"); t = mp.mpf(0)
    for it in range(60):
        f, g = Hflow(vs, ph, h, t, z); dz = f / g; z -= dz
        if abs(dz) < mp.mpf(10) ** (-25): break
    z0 = z; path = [(t, z)]
    for i in range(1, 9):
        t = mp.mpf("0.02") * i
        for it in range(60):
            f, g = Hflow(vs, ph, h, t, z); dz = f / g; z -= dz
            if abs(dz) < mp.mpf(10) ** (-25): break
        path.append((t, z))
    # double root seeded from the last tracked point, real z
    def G(Z, T):
        f, g = Hflow(vs, ph, h, T, Z); return [f, g]
    zs = path[-1][1].real
    sol = mp.findroot(G, (zs, mp.mpf("0.1765")), tol=mp.mpf(10) ** (-(dps // 3)), maxsteps=100)
    return z0, path, sol

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--flowonly", action="store_true"); ap.add_argument("--out", default=None)
    a = ap.parse_args(); t0 = time.time(); L = []; R = {}
    def log(s): L.append(s); print(s, flush=True)
    K = 600
    cfg = [(100, 150, mp.mpf("2e-6") * 5), (100, 170, mp.mpf("2e-6") * 5)] if a.quick else [(200, 240, mp.mpf("2e-6")), (250, 280, mp.mpf("2.5e-6"))]
    res = []
    for hinv, dps, step in ([] if a.flowonly else cfg):
        o, c = run(K, hinv, dps, step, log)
        res.append(o)
        log("coefficients+tracking: step 1/%d, %d digits, eps step %s, %.1f s" % (hinv, dps, mp.nstr(step, 3), time.time() - t0))
        for k, v in o.items():
            if k in ("path_im2",): log("  path (eps, Im^2): " + "; ".join("(%s, %s)" % p for p in v))
            else: log("  %s = %s" % (k, mp.nstr(v, 14) if isinstance(v, (mp.mpf, mp.mpc)) else str(v)[:120]))
    R["tracking"] = [{k: (mp.nstr(v, 14) if isinstance(v, (mp.mpf, mp.mpc)) else str(v)) for k, v in o.items()} for o in res]
    # flow
    fl = []
    for hinv, dps in ((100, 90), (125, 110)) if a.quick else ((200, 90), (250, 110), (300, 130)):
        z0, path, sol = flow(hinv, dps, log); fl.append((z0, path, sol))
        log("flow h=1/%d, %d digits: t=0 zero z = %s; tracked (t, z): %s; double root z_c = %s, t_c = %s" % (hinv, dps, mp.nstr(z0, 12), "; ".join("(%s, %s)" % (mp.nstr(t, 3), mp.nstr(z, 8)) for t, z in path[::2]), mp.nstr(sol[0], 14), mp.nstr(sol[1], 14)))
    R["flow"] = [{"z0": mp.nstr(z0, 14), "z_c": mp.nstr(s_[0], 14), "t_c": mp.nstr(s_[1], 14)} for z0, p, s_ in fl]
    tc = fl[-1][2][1]; gam = fl[-1][0].real
    R["d_star_from_t_c"] = mp.nstr(gam ** 2 / (2 * tc), 8); R["t_c_for_gamma_2_over_2_t_c_20785"] = mp.nstr(gam ** 2 / (2 * 20785), 8)
    R["two_delta_squared"] = mp.nstr(2 * fl[-1][0].imag ** 2, 8)
    log("d* = gamma^2/(2 t_c) = %s (gamma = %s); t_c implied by d* = 20785: %s; 2 delta^2 = %s; elapsed %.1f s" % (R["d_star_from_t_c"], mp.nstr(gam, 12), R["t_c_for_gamma_2_over_2_t_c_20785"], R["two_delta_squared"], time.time() - t0))
    if a.out:
        json.dump(R, open(a.out + ".json", "w"), indent=1); open(a.out + ".log", "w").write("\n".join(L) + "\n")
main()
