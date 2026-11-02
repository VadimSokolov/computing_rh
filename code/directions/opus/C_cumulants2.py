#!/usr/bin/env python3
"""Track C, script 1b. Larger Hankel determinants for the fair sign test, and controls on true sums.
Run: timeout 500 python3 -I C_cumulants2.py   (computed in mpmath, not certified)
"""
import mpmath as mp
mp.mp.dps = 700


def Phi(t):
    s = mp.mpf(0)
    e2 = mp.exp(2 * t)
    a, b = mp.exp(4.5 * t), mp.exp(2.5 * t)
    for n in range(1, 36):
        s += (2 * mp.pi**2 * n**4 * a - 3 * mp.pi * n**2 * b) * mp.exp(-mp.pi * n * n * e2)
    return 2 * s


def moments(jmax):
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)  # 96 nodes
    xs, ws = [], []
    h = mp.mpf(3) / 60
    for i in range(60):
        mid = (i + mp.mpf(1) / 2) * h
        for x, w in nodes:
            xs.append(mid + x * h / 2)
            ws.append(w * h / 2)
    vals = [Phi(x) for x in xs]
    mom = []
    x2 = [x * x for x in xs]
    pw = [w * v for w, v in zip(ws, vals)]
    for j in range(jmax + 1):
        mom.append(2 * mp.fsum(pw))
        pw = [u * y for u, y in zip(pw, x2)]
    return mom


def cumulants(mom):
    n = 2 * (len(mom) - 1)
    m = [mp.mpf(0)] * (n + 1)
    for j, v in enumerate(mom):
        m[2 * j] = v / mom[0]
    k = [mp.mpf(0)] * (n + 1)
    binom = [[mp.mpf(0)] * (n + 1) for _ in range(n + 1)]
    for a in range(n + 1):
        binom[a][0] = mp.mpf(1)
        for b in range(1, a + 1):
            binom[a][b] = binom[a - 1][b - 1] + (binom[a - 1][b] if b <= a - 1 else 0)
    for r in range(2, n + 1, 2):
        k[r] = m[r] - mp.fsum(binom[r - 1][i - 1] * k[i] * m[r - i] for i in range(2, r, 2))
    return k


def dets(seq, nmax):
    out = []
    for n in range(1, nmax + 1):
        A = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                A[i, j] = seq[i + j]
        out.append(mp.det(A))
    return out


def main():
    jmax = 62
    mom = moments(jmax)
    k = cumulants(mom)
    print("xi(1/2) =", mp.nstr(mom[0], 22), " kappa_4 =", mp.nstr(k[4], 15))
    p = {m: (-1) ** (m + 1) * m * k[2 * m] / mp.factorial(2 * m) for m in range(1, jmax + 1)}
    r = [(mp.pi / 2) ** (2 * m) * p[m] / ((1 - mp.mpf(2) ** (-2 * m)) * mp.zeta(2 * m)) for m in range(2, jmax + 1)]
    x = r[-1] / r[-2]
    rn = [r[j] / x ** (j + 2) for j in range(len(r))]
    q = [mp.pi ** (2 * m) * p[m] / mp.zeta(2 * m) for m in range(2, jmax + 1)]
    xq = q[-1] / q[-2]
    qn = [q[j] / xq ** (j + 2) for j in range(len(q))]
    Hs = dets(rn, 30)
    Hu = dets(qn, 12)
    print("fair sign test: signs of H_n, n=1..30:", "".join("+" if d > 0 else "-" for d in Hs))
    print("first negative at n =", next((i + 1 for i, d in enumerate(Hs) if d < 0), None))
    print("uniform test  : signs of H_n, n=1..12:", "".join("+" if d > 0 else "-" for d in Hu))
    print("H_7, uniform:", mp.nstr(Hu[6], 8))
    # stability against the precision of the input: Hankel size versus digits lost
    print("smallest |H_n| fair sign:", mp.nstr(min(abs(d) for d in Hs), 5))
    # control: a true sum of uniforms, harmonic coefficients, must give positive determinants
    b = mp.mpf("1.25839697")
    ctrl = [mp.mpf(1) / 2 ** (2 * m) * mp.zeta(2 * m, 1 + b) for m in range(2, 26)]
    x_c = ctrl[-1] / ctrl[-2]
    ctrln = [ctrl[j] / x_c ** (j + 2) for j in range(len(ctrl))]
    print("control, harmonic uniform law, signs n=1..12:", "".join("+" if d > 0 else "-" for d in dets(ctrln, 12)))


if __name__ == "__main__":
    main()
