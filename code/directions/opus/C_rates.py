#!/usr/bin/env python3
"""Track C, script 6. Rates and finite prime laws.

(a) TV(f_K, p) against (1/4) tau_K ||p''''||_1, tau_K = sum_{k>K} gamma_k^{-4}
(b) the finite integer law: Phi_M(|t|) = terms n <= M of the theta series, symmetrised. Its total
    variation distance to p, the corner term, and the height where its real zeros stop tracking gamma_k
(c) Gaussian damping: TV(p_eta, p) for the laws e^{-eta x^2} p of Rule 3

Run: timeout 500 python3 -I C_rates.py
Status: computed in mpmath and double precision, not certified.
"""
import numpy as np
import mpmath as mp

XI0 = 0.4971207781883141
VARX = 0.0462099862308379
SUM4 = 3.71725992852697e-5
mp.mp.dps = 50


def p_density(x):
    x = np.abs(np.asarray(x, dtype=float))
    s = np.zeros_like(x)
    for n in range(1, 8):
        y = np.pi * n**2 * np.exp(2 * x)
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * x) - 3 * np.pi * n**2 * np.exp(2.5 * x)) * np.exp(-y)
    return 2 * s / XI0


def part_a():
    print("(a) f_K against the fourth derivative law")
    zeros = np.array([float(mp.zetazero(k).imag) for k in range(1, 321)])
    xs = np.arange(0, 2.4001, 0.0005)
    pv = p_density(xs)
    th = np.arange(0, 140.0001, 0.01)
    w = np.ones_like(th) * 0.01
    w[0] = w[-1] = 0.005
    # phi_X by cosine transform of p (double precision is enough here)
    xq = np.arange(0, 3.0001, 0.0005)
    pq = p_density(xq)
    phiX = np.array([2 * np.sum(pq * np.cos(t * xq)) * 0.0005 - pq[0] * 0.0005 for t in th[th <= 45]])
    th45 = th[th <= 45]
    w45 = w[: len(th45)]
    p4 = np.array([np.sum(w45 * th45**4 * phiX * np.cos(x * th45)) / np.pi for x in xs[::10]])
    norm4 = 2 * (np.sum(np.abs(p4)) - 0.5 * abs(p4[0]) - 0.5 * abs(p4[-1])) * 0.005
    print(f"  ||p''''||_1 = {norm4:.3f}")
    for K in (16, 50, 100, 150, 200, 320):
        g = zeros[:K]
        betaK = 0.5 * VARX - np.sum(g**-2.0)
        tau = SUM4 - np.sum(g**-4.0)
        logabs = -betaK * th**2
        sign = np.ones_like(th)
        for z in g:
            f = 1 - th**2 / z**2
            logabs += np.log(np.maximum(np.abs(f), 1e-300))
            sign *= np.sign(f)
        phiK = sign * np.exp(logabs)
        fK = np.array([np.sum(w * phiK * np.cos(x * th)) / np.pi for x in xs[::2]])
        d = np.abs(fK - pv[::2])
        tv = (np.sum(d) - 0.5 * d[0] - 0.5 * d[-1]) * 0.001
        print(f"  K={K:4d} gamma_K={g[-1]:9.3f} tau_K={tau:.4e} TV(f_K,p)={tv:.4e}  (1/4) tau_K ||p''''||={0.25*tau*norm4:.4e} ratio={tv/(0.25*tau*norm4):.4f}")


def xi(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def part_b():
    print("(b) finite integer laws: Phi_M(|t|) keeps n <= M; transform = Xi(th) - 2 int_0^inf R_M cos(th t) dt, up to the norm")
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)  # 48 nodes
    panels = 40
    h = mp.mpf("1.4") / panels
    ts, ws = [], []
    for i in range(panels):
        mid = (i + mp.mpf(1) / 2) * h
        for x, wg in nodes:
            ts.append(mid + x * h / 2)
            ws.append(wg * h / 2)
    term = lambda n, t: (2 * mp.pi**2 * n**4 * mp.exp(4.5 * t) - 3 * mp.pi * n**2 * mp.exp(2.5 * t)) * mp.exp(-mp.pi * n * n * mp.exp(2 * t))
    dterm = lambda n, t: mp.diff(lambda u: term(n, u), t)
    zeros = [mp.zetazero(k).imag for k in range(1, 60)]
    xi0 = xi(mp.mpf(1) / 2)
    for M in (1, 2, 3, 4):
        RM = [2 * sum(term(n, t) for n in range(M + 1, M + 12)) for t in ts]  # 2 * sum over n>M, as in Phi
        intR = 2 * mp.fsum(wg * r for wg, r in zip(ws, RM))  # mass of the difference on the line
        norm = xi0 - intR
        # corner: the derivative of R_M at 0
        dR0 = 2 * sum(dterm(n, mp.mpf(0)) for n in range(M + 1, M + 8))
        # TV(p^(M), p) = int_0^inf |Phi_M(t)/norm - Phi(t)/xi0| dt (both even, whole line mass 1)
        Phi_vals = [2 * sum(term(n, t) for n in range(1, M + 12)) for t in ts]
        PhiM_vals = [2 * sum(term(n, t) for n in range(1, M + 1)) for t in ts]
        tv = mp.fsum(wg * abs(a / norm - b / xi0) for wg, a, b in zip(ws, PhiM_vals, Phi_vals))
        # real zeros of the transform and where they stop. The remainder integral needs relative
        # accuracy only, so double precision on the scaled remainder is enough
        tmax = 4 * (M + 1) ** 2 + 30
        grid = [mp.mpf(j) / 10 for j in range(5, int(10 * tmax))]
        t_np = np.array([float(t) for t in ts])
        w_np = np.array([float(wg) for wg in ws])
        scale = max(abs(r) for r in RM)
        r_np = np.array([float(r / scale) for r in RM])
        vals = []
        for th in grid:
            c = 2 * scale * mp.mpf(float(np.sum(w_np * r_np * np.cos(float(th) * t_np))))
            vals.append((xi(mp.mpf(1) / 2 + 1j * th).real - c) / norm)
        sc = [grid[i] for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0]
        # crossing of the envelope: 2 |dR0| / th^2 = 2.2519 * xi0 * th^(7/4) e^(-pi th/4)
        cross = mp.findroot(lambda th: mp.log(2 * abs(dR0)) - 2 * mp.log(th) - mp.log(mp.mpf("2.2519") * xi0) - 1.75 * mp.log(th) + mp.pi * th / 4, 4 * (M + 1) ** 2)
        nearest = [mp.nstr(min((abs(z - s) for z in zeros)), 2) for s in sc]
        print(f"  M={M}: TV(p^(M),p)={mp.nstr(tv,3)}; R_M'(0)={mp.nstr(dR0,3)}; envelope crossing {mp.nstr(cross,5)} (4(M+1)^2={4*(M+1)**2}); last real sign change on (0.5,{tmax}): {mp.nstr(sc[-1],5) if sc else None}; count {len(sc)}; zeros gamma below the last: {sum(1 for z in zeros if sc and z <= sc[-1])}")
        print(f"       distance of each sign change to the nearest gamma: {nearest[:22]}")


def part_c():
    print("(c) Gaussian damping of Rule 3: p_eta = e^{-eta x^2} p / E e^{-eta X^2}; lambda = -4 eta in the flow")
    xs = np.arange(0, 3.0001, 0.0005)
    pv = p_density(xs)
    mean_abs = 2 * np.sum(pv * np.abs(xs**2 - VARX)) * 0.0005
    print(f"  E|X^2 - Var X| = {mean_abs:.6f}, so TV(p_eta,p) ~ eta * {mean_abs/2:.6f} for small eta")
    for eta in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6):
        q = np.exp(-eta * xs**2) * pv
        q = q / (2 * np.sum(q) * 0.0005 - q[0] * 0.0005)
        print(f"  eta={eta:.0e}: TV={np.sum(np.abs(q-pv))*0.0005:.4e}")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
