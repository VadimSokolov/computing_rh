#!/usr/bin/env python3
"""PROPOSAL, NEVER EXECUTED: python-flint is on Hopper, not on the laptop. Check the API calls on the first run.

Certify one zero of f(s) = zeta(s) + eps * r^s * L_c(s) in ball arithmetic, with c the closed form of
Theorem 1(a) of B.md. In Re s > 1/2 a zero of f is a zero of F_eps = xi + eps*M, as the factor
s(s-1) pi^(-s/2) Gamma(s/2)/2 has no zero there.
Usage on Hopper: python3 -I hopper_certify_fake.py N eps sigma0 t0 halfwidth prec
Examples:
  python3 -I hopper_certify_fake.py 2 1 1.428166777076261 88.15330268995843 1e-6 512
  python3 -I hopper_certify_fake.py 4 1 1.075265148537528 357.5759059228842 1e-6 512
Output: the largest ball residual of self duality, and a ball for the winding number along the square.
A ball containing exactly one integer certifies the number of zeros in the square.
"""
import sys
from flint import acb, arb, ctx

N = int(sys.argv[1]); eps = arb(sys.argv[2]); sigma0, t0, hw = sys.argv[3], sys.argv[4], sys.argv[5]
ctx.prec = int(sys.argv[6])
r = 4 * N; q = r * r; K = N * r; half = q // 2
pi = arb.pi()
y = [(2 * pi * j / q).cos() for j in range(K)]
def W(x):
    p = arb(1)
    for yj in y:
        p *= x - yj
    return p
W1 = W(arb(-1))
c = [arb(0)] * q
for n in range(K, half + 1):
    v = W((2 * pi * n / q).cos()) / W1
    if n % 2:
        v = -v
    c[n] = v; c[(q - n) % q] = v
# self duality in balls: every residual ball must contain 0 and be small
sq = arb(q).sqrt(); worst = arb(0)
for m in range(half + 1):
    res = sum((c[b] * (2 * pi * m * b / q).cos() for b in range(q) if b >= K and b <= q - K), arb(0)) / sq - c[m]
    if not res.contains(0):
        raise SystemExit(f"self duality fails at m={m}: {res}")
    if abs(res).upper() > worst:
        worst = abs(res).upper()
print("largest self duality residual radius:", worst)

def f(s):
    tot = acb(0)
    for b in range(K, q - K + 1):        # a = b/q; c_b = 0 outside this range
        tot += c[b] * s.zeta(acb(b) / q)
    return s.zeta() + eps * (acb(r) / q) ** s * tot

s0 = acb(arb(sigma0), arb(t0)); h = arb(hw)
corners = [s0 + acb(-h, -h), s0 + acb(h, -h), s0 + acb(h, h), s0 + acb(-h, h)]
def winding(M_seg=64):
    total = arb(0)
    for k in range(4):
        a, b = corners[k], corners[(k + 1) % 4]
        prev = f(a)
        for j in range(1, M_seg + 1):
            cur = f(a + (b - a) * j / M_seg)
            ratio = cur / prev
            if not ratio.real > 0:
                raise SystemExit("subdivide: the ratio is not certainly in the right half plane")
            total += ratio.log().imag
            prev = cur
    return total / (2 * pi)
print("winding number ball:", winding())
print("square lies in Re s > 1:", (s0.real - h) > 1)
