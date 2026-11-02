#!/usr/bin/env python3
"""Track C, proposal P1: ball arithmetic certificates for the numbers of C.md. NOT RUN on the laptop (no python-flint).
Hopper: one core, 4 GB, under one hour. Written against python-flint 0.6; expect small API fixes.
Every printed inequality follows from the enclosure printed with it, nothing else is claimed.
"""
from flint import arb, acb, arb_mat, arb_series, ctx


def xi(s):
    return s * (s - 1) / 2 * acb.pi() ** (-s / 2) * (s / 2).gamma() * s.zeta()


def phi(theta):
    """phi_X(theta) = xi(1/2 + i theta)/xi(1/2), for a real ball theta."""
    half = acb(1) / 2
    return xi(half + acb(0, 1) * acb(theta)) / xi(half)


def first_zero():
    """Isolate gamma_1 by bisection on the sign of the real function phi_X."""
    lo, hi = arb("14.1"), arb("14.2")
    for _ in range(40):
        mid = (lo + hi) / 2
        if phi(mid).real > 0:
            lo = mid
        else:
            hi = mid
    return lo.union(hi)


def part_floor():
    ctx.prec = 256
    g = first_zero()
    print("gamma_1 in", g.str(25))
    dphi = (phi(g + arb("1e-6")) - phi(g - arb("1e-6"))) / arb("2e-6")
    print("phi_X'(gamma_1), difference quotient only, for orientation:", dphi.str(10))
    for k, name in ((2, "2 gamma_1"), (3, "3 gamma_1")):
        s = acb(1) / 2 + acb(0, 1) * (k * g)
        print(f"|zeta(1/2 + i {name})| in", abs(s.zeta()).str(15), " |phi_X| in", abs(phi(k * g)).str(10))
    ball = g + arb(0, "2.4e-4")
    two_over_pi = 2 / arb.pi()
    main = phi(2 * ball).real
    corr = phi(6 * ball).real / 3
    print("uniform floor: (2/pi)|phi_X(2t) - phi_X(6t)/3| on the ball is at least", (two_over_pi * abs(main - corr)).lower().str(8))
    print("remaining terms j >= 2: phi_X((4j+2)t)/(2j+1), add the amplitude bound, |phi_X(10 t)| in", abs(phi(10 * ball)).str(5))
    main3 = phi(3 * ball).real
    corr3 = phi(9 * ball).real / 3
    print("sign floor: at least", (two_over_pi * abs(main3 - corr3)).lower().str(8))


def part_lobes():
    ctx.prec = 256
    tops = [15.5857, 22.098, 26.2722, 31.2318, 34.1933, 38.4982, 41.7367, 44.542, 48.623, 50.839, 53.969, 62.1099, 88.0, 123.6]
    for t in tops:
        print(f"|phi_X({t})| in", abs(phi(arb(str(t)))).str(10))


def part_hankel():
    ctx.prec = 20000  # 6000 digits
    ctx.cap = 72
    s = arb_series([arb(1) / 2, 1])
    series = s * (s - 1) / 2 * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta()
    log_xi = series.log()
    fact = arb(1)
    kappa = {}
    for r in range(1, 71):
        fact = fact * r
        if r % 2 == 0:
            kappa[r] = log_xi[r] * fact
    print("kappa_2", kappa[2].str(20), "kappa_4", kappa[4].str(20))
    fac = arb(1)
    p = {}
    f = arb(1)
    for m in range(1, 36):
        f = (arb(2 * m).gamma() * (2 * m))  # (2m)!
        p[m] = (-1) ** (m + 1) * m * kappa[2 * m] / f
    zeta2 = {m: arb(2 * m).zeta() for m in range(2, 36)}
    q = {m: arb.pi() ** (2 * m) * p[m] / zeta2[m] for m in range(2, 36)}
    r = {m: (arb.pi() / 2) ** (2 * m) * p[m] / ((1 - arb(2) ** (-2 * m)) * zeta2[m]) for m in range(2, 36)}

    def dets(seq, nmax):
        out = []
        for n in range(1, nmax + 1):
            entries = [seq[2 + i + j] for i in range(n) for j in range(n)]
            out.append(arb_mat(n, n, entries).det())
        return out

    Hu = dets(q, 16)
    Hs = dets(r, 16)
    Hp = dets(p, 16)
    print("uniform test H_n signs, ball excludes 0?:", [(n + 1, d > 0, d < 0) for n, d in enumerate(Hu)])
    print("H_7 uniform in", Hu[6].str(8))
    print("sign test H_16 in", Hs[15].str(8))
    print("p_m test positive for n <= 16:", all(d > 0 for d in Hp))


def part_dh():
    ctx.prec = 256
    k = ((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
    coeff = [arb(1), k, -k, arb(-1)]

    def f(s):
        return acb(5) ** (-s) * sum((acb(c) * s.zeta(acb(r + 1) / 5) for r, c in enumerate(coeff)), acb(0))

    def lam(s):
        return (acb(5) / acb.pi()) ** (s / 2) * ((s + 1) / 2).gamma() * f(s)

    half = acb(1) / 2
    z0 = acb("0.808517182457", "85.6993484854")
    print("|f(z0)| in", abs(f(z0)).str(5))
    zc = z0 - half
    radius = arb("0.3")
    arcs = 2000
    best = None
    for j in range(arcs):
        psi = (arb(j) + arb("0.5") + arb(0, "0.5")) * 2 * arb.pi() / arcs
        z = zc + radius * acb(psi.cos(), psi.sin())
        v = abs(lam(half + z) / lam(half)).lower()
        best = v if best is None or v < best else best
    print("lower bound of |L_DH| on the circle of radius 0.3:", best.str(6))


if __name__ == "__main__":
    part_floor()
    part_lobes()
    part_dh()
    part_hankel()
