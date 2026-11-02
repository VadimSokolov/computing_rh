#!/usr/bin/env python3
"""Track C, script 1. Cumulants of the Polya law from the theta series only (no zeros),
then the zero free Hankel tests for sums of scaled uniforms and of scaled fair signs.

Run: timeout 300 python3 -I C_cumulants.py
Status of the output: computed in mpmath at high precision, not certified ball arithmetic.
"""
import json
import os
import mpmath as mp

mp.mp.dps = 160


def Phi(t):
    """Polya kernel for t >= 0, with the theta series truncated at n = 12."""
    t = mp.mpf(t)
    s = mp.mpf(0)
    for n in range(1, 13):
        y = mp.pi * n * n * mp.exp(2 * t)
        s += (2 * mp.pi**2 * n**4 * mp.exp(4.5 * t) - 3 * mp.pi * n**2 * mp.exp(2.5 * t)) * mp.exp(-y)
    return 2 * s


def gauss_panels(a, b, npanels, order):
    """Composite Gauss Legendre nodes and weights on [a, b]."""
    xs, ws = [], []
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(order, mp.mp.prec)
    h = (mp.mpf(b) - a) / npanels
    for i in range(npanels):
        mid = a + (i + mp.mpf(1) / 2) * h
        for x, w in nodes:
            xs.append(mid + x * h / 2)
            ws.append(w * h / 2)
    return xs, ws


def moments(jmax):
    """E X^(2j), j = 0..jmax, normalised by the integral of Phi on the line."""
    xs, ws = gauss_panels(mp.mpf(0), mp.mpf(3), 60, 5)  # 60 panels, degree 3*2^4 rule
    vals = [Phi(x) for x in xs]
    mom = []
    for j in range(jmax + 1):
        mom.append(2 * mp.fsum(w * v * x ** (2 * j) for x, w, v in zip(xs, ws, vals)))
    return mom


def cumulants(mom):
    """Even cumulants from even moments of a symmetric law, by the standard recursion."""
    n = 2 * (len(mom) - 1)
    m = [mp.mpf(0)] * (n + 1)
    for j, v in enumerate(mom):
        m[2 * j] = v / mom[0]
    k = [mp.mpf(0)] * (n + 1)
    for r in range(1, n + 1):
        k[r] = m[r] - mp.fsum(mp.binomial(r - 1, i - 1) * k[i] * m[r - i] for i in range(1, r))
    return k


def hankel_dets(seq, nmax):
    """det (seq[i+j]) for sizes 1..nmax, after scaling by the first ratio."""
    out = []
    for n in range(1, nmax + 1):
        if 2 * n - 1 > len(seq):
            break
        A = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                A[i, j] = seq[i + j]
        out.append(mp.det(A))
    return out


def main():
    jmax = 26
    mom = moments(jmax)
    xi_half = mom[0]
    k = cumulants(mom)
    print("xi(1/2) from the quadrature :", mp.nstr(xi_half, 25))
    print("direct control              :", mp.nstr(mp.mpf("0.4971207781883141"), 16))
    print("kappa_2,4,6 =", mp.nstr(k[2], 15), mp.nstr(k[4], 15), mp.nstr(k[6], 15))
    # power sums of the zeros, no hypothesis (complex zeros allowed), m >= 1
    p = {}
    for m in range(1, jmax + 1):
        p[m] = (-1) ** (m + 1) * m * k[2 * m] / mp.factorial(2 * m)
    print("Sum gamma^-2, gamma^-4 =", mp.nstr(p[1], 15), mp.nstr(p[2], 15))
    # uniform test: q_m = pi^(2m) p_m / zeta(2m), m >= 2, must be sum c_k^(2m)
    q = [mp.pi ** (2 * m) * p[m] / mp.zeta(2 * m) for m in range(2, jmax + 1)]
    # fair sign test: r_m = (pi/2)^(2m) p_m / ((1-2^(-2m)) zeta(2m)), must be sum a_k^(2m)
    r = [(mp.pi / 2) ** (2 * m) * p[m] / ((1 - mp.mpf(2) ** (-2 * m)) * mp.zeta(2 * m)) for m in range(2, jmax + 1)]
    # normalise to keep the determinants readable: divide term j by x^(j+2), x the limit ratio
    xu = q[-1] / q[-2]
    xs = r[-1] / r[-2]
    print("c_1^2 implied by the ratios :", mp.nstr(xu, 15), " so pi/c_1 =", mp.nstr(mp.pi / mp.sqrt(xu), 15))
    print("a_1^2 implied by the ratios :", mp.nstr(xs, 15), " so pi/(2 a_1) =", mp.nstr(mp.pi / (2 * mp.sqrt(xs)), 15))
    qn = [q[j] / xu ** (j + 2) for j in range(len(q))]
    rn = [r[j] / xs ** (j + 2) for j in range(len(r))]
    print("mass of the top atom (uniform normalisation):", mp.nstr(qn[-1], 15))
    print("mass of the top atom (sign normalisation)   :", mp.nstr(rn[-1], 15))
    Hu = hankel_dets(qn, 12)
    Hs = hankel_dets(rn, 12)
    print("\nHankel determinants, uniform test (no zeros used):")
    for n, d in enumerate(Hu, 1):
        print(f"  n={n:2d}  H_n = {mp.nstr(d, 8)}")
    print("Hankel determinants, fair sign test (no zeros used):")
    for n, d in enumerate(Hs, 1):
        print(f"  n={n:2d}  H_n = {mp.nstr(d, 8)}")
    # control: the same determinants from the first zeros and the signed atoms
    zs = [mp.zetazero(i).imag for i in range(1, 41)]
    print("\ncontrol with zeros: gamma_1 =", mp.nstr(zs[0], 15))
    ctrl = [sum(z ** (-2 * m) for z in zs) for m in range(2, 12)]
    print("relative error of p_m from the first 40 zeros, m=2..11:")
    print("  ", [mp.nstr((ctrl[m - 2] - p[m]) / p[m], 3) for m in range(2, 12)])
    print("|zeta(1/2+2i gamma_1)| =", mp.nstr(abs(mp.zeta(mp.mpc(0.5, 2 * zs[0]))), 12))
    print("|zeta(1/2+3i gamma_1)| =", mp.nstr(abs(mp.zeta(mp.mpc(0.5, 3 * zs[0]))), 12))
    out = {"kappa": {str(2 * m): mp.nstr(k[2 * m], 60) for m in range(1, jmax + 1)}, "xi_half": mp.nstr(xi_half, 60)}
    with open("/tmp/C_cumulants.json", "w") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main()
