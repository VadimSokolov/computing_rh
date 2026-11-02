#!/usr/bin/env python3
"""Track A, proposal for Hopper. Needs python-flint (the book's Hopper venv).

  saddle: for each gap of a list of zeta zeros, ball enclosure of the saddle of the
          Green potential Sigma(1/2+it) = sum log coth(pi |t-gamma|/2).
  recon:  Identity B, xi(s0) = integral Re xi(1+iv)/sin(pi(s0-iv)) dv, by the trapezoid
          rule, with noise, against direct xi(s0) and the bound eps*exp(Sigma).

Usage:
  python A_hopper_cert.py saddle --zeros FILE --lo A --hi B --rad 1e-9 [--L 15]
  python A_hopper_cert.py recon --H 1e4 --a 0.1 --eps 1e-6 [--h 0.05] [--halfw 30]

Certified means ball arithmetic conditional on the inputs: the zero list is right to the
radius given, complete, and on the line. The count constants below are from the literature
(|S(T)| <= 0.112 log T + 0.278 log log T + 2.51, Trudgian 2014, weaker than HSW22 and enough
for an upper bound) and must be checked before any use.
"""
import argparse, bisect, math, random
from flint import arb, acb, ctx


def count_bound(u):
    """Upper bound for the number of zeros in [u, u+1], u >= 3."""
    l = math.log(u + 1)
    return l / (2 * math.pi) + 2 * (0.112 * l + 0.278 * math.log(l) + 2.51) + 1


def dsigma_float(t, near):
    d = 0.0
    for g in near:
        x = t - g
        d -= math.copysign(math.pi / math.sinh(math.pi * abs(x)), x)
    return d


def saddle(args):
    ctx.prec = args.prec
    gs = [float(x) for x in open(args.zeros) if x.strip()]
    gs.sort()
    i0 = bisect.bisect_left(gs, args.lo)
    i1 = bisect.bisect_right(gs, args.hi)
    L = args.L
    pi = arb.pi()
    worst = None
    print("left right sigma_upper sigma_lower exp(-sigma_upper)")
    for j in range(i0, i1 - 1):
        a, b = gs[j], gs[j + 1]
        m0 = bisect.bisect_left(gs, a - L)
        m1 = bisect.bisect_right(gs, b + L)
        near = gs[m0:m1]
        if gs[0] > a - L or gs[-1] < b + L:
            raise SystemExit("zero list does not cover the margin L")
        # float location of the saddle (only selects t*, the balls below are what counts)
        lo, hi = a + 1e-12, b - 1e-12
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if dsigma_float(mid, near) < 0:
                lo = mid
            else:
                hi = mid
        tf = 0.5 * (lo + hi)
        ts = arb(tf)
        sig = arb(0)
        der = arb(0)
        for g in near:
            z = arb(g, args.rad)
            ax = ts - z if g < tf else z - ts
            sign = 1 if g < tf else -1
            sig += -(pi * ax / 2).tanh().log()
            der += (-pi / (pi * ax).sinh()) * sign
        # tail: log coth(pi x/2) = 2 artanh(exp(-pi x)) <= 2.01 exp(-pi x) for x >= L >= 2,
        # zeros per unit interval <= B + k at distance L + k, both sides
        q = math.exp(-math.pi)
        B = count_bound(tf + L + 1)
        tail = 2 * 2.01 * math.exp(-math.pi * L) * (B + 2) / (1 - q)
        dtail = math.pi * tail
        S = sig + arb(tail / 2, tail / 2)  # encloses Sigma(t*)
        up = S.upper()
        low = (sig - (der.abs_upper() + arb(dtail)) * arb(b - a)).lower()  # convexity on the gap
        e = (-up).exp().lower()  # lower bound for exp(-saddle)
        print(f"{a:.6f} {b:.6f} {up.str(8, radius=False)} {low.str(8, radius=False)} {e.str(8, radius=False)}")
        if worst is None or float(e.str(12, radius=False)) < worst[0]:
            worst = (float(e.str(12, radius=False)), a, b)
    if worst:
        print("window eps_iso lower bound:", worst[0], "at gap", worst[1], worst[2])


def xi(s):
    pi = acb.pi()
    return s * (s - 1) / 2 * (-(s / 2) * pi.log()).exp() * (s / 2).gamma() * s.zeta()


def recon(args):
    ctx.prec = args.prec
    random.seed(args.seed)
    t = args.H + 0.123
    s0 = acb(arb(0.5 + args.a), arb(t))
    direct = xi(s0)
    k = int(args.halfw / args.h)
    pi = acb.pi()
    s_clean, s_noisy = acb(0), acb(0)
    for j in range(-k, k + 1):
        v = t + j * args.h
        x = xi(acb(arb(1), arb(v)))
        A = x.real
        K = 1 / (pi * (s0 - acb(0, arb(v)))).sin()
        s_clean += acb(A) * K
        noise = args.eps * (2 * random.random() - 1)
        s_noisy += acb(A + arb(noise) * x.abs_upper()) * K  # noise relative to |xi(1+iv)|
    s_clean *= args.h
    s_noisy *= args.h
    print("H", args.H, "a", args.a, "|direct - clean|/|direct|", ((direct - s_clean).abs_upper() / direct.abs_lower()).str(5))
    print("|direct - noisy|/|direct|", ((direct - s_noisy).abs_upper() / direct.abs_lower()).str(5), "eps", args.eps)
    print("|xi(1+it)|/|direct| (e^Sigma scale)", (xi(acb(arb(1), arb(t))).abs_upper() / direct.abs_lower()).str(5))


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("saddle")
    a.add_argument("--zeros", required=True)
    a.add_argument("--lo", type=float, required=True)
    a.add_argument("--hi", type=float, required=True)
    a.add_argument("--rad", type=float, default=1e-9)
    a.add_argument("--L", type=float, default=15.0)
    a.add_argument("--prec", type=int, default=128)
    b = sub.add_parser("recon")
    b.add_argument("--H", type=float, required=True)
    b.add_argument("--a", type=float, default=0.1)
    b.add_argument("--eps", type=float, default=1e-6)
    b.add_argument("--h", type=float, default=0.05)
    b.add_argument("--halfw", type=float, default=30.0)
    b.add_argument("--prec", type=int, default=128)
    b.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    {"saddle": saddle, "recon": recon}[args.cmd](args)


if __name__ == "__main__":
    main()
