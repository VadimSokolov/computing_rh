# px_const.py: PX.1 of the reconciliation of the notes on the pole cancellation and the excursion (Chapter 22, Remark mi:rem:laws, and
# Chapter ht, Remark ht:rem:chung): the constants G(0) = 2 sum_n (pi n^2)^(-1/4) Gamma(1/4, pi n^2) (PX.1a), eps = G(0)/4 compared with
# 1 - 2 xi(1/2) (PX.1b), and P(Y > 1) = 2 sum_n (2 pi n^2 - 1) exp(-pi n^2) compared with 1/2 - psi(1) and with a quadrature of
# 2 Phi(u) exp(-u/2) on u > 0 (PX.1c). Y = sqrt(2/pi) M, M the maximum of a standard Brownian excursion, psi(x) = sum exp(-pi n^2 x),
# Phi the kernel of eq:ch2:Phi. Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Method. The series are summed in ball arithmetic (python-flint acb_hypgeom_gamma_upper, 400 bits) to n = 12 and n = 14; the tail beyond
# n = 12 is bounded by 3 exp(-pi 169) (Gamma(1/4, a) <= 1.5 a^(-3/4) exp(-a) for a >= 3, so the term n is at most 3 exp(-pi n^2) / (pi n^2));
# xi(1/2) = (1/2)(1/2)(-1/2) pi^(-1/4) Gamma(1/4) zeta(1/2) is evaluated in balls. The quadrature is mpmath at 30 digits (computed).
# Status of the values: certified (balls), except the quadrature (computed).
# Usage: python3 px_const.py [--quick] [--out FILE]
import argparse, json, sys, time
from flint import acb, arb, ctx
import mpmath

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default='px_const.json')
args = ap.parse_args()
ctx.prec = 400
out = {}
lines = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def G0(N):
    # mpmath at 60 digits (python-flint's gamma_upper loses all accuracy for pi n^2 > about 400, so the large n terms are not taken from it)
    mpmath.mp.dps = 60
    return 2 * sum((mpmath.pi * n * n) ** (-mpmath.mpf(1) / 4) * mpmath.gammainc(mpmath.mpf(1) / 4, mpmath.pi * n * n) for n in range(1, N + 1))


def S(x, nd=30):
    return x.str(nd, radius=False)


g12, g14 = G0(12), G0(14)
P('G(0) series (mpmath, 60 digits) N=12', mpmath.nstr(g12, 45), 'N=14', mpmath.nstr(g14, 45), 'diff', mpmath.nstr(g14 - g12, 5))
# flint check of the first terms (n <= 4)
fl = arb(0)
for n in range(1, 5):
    a = arb.pi() * n * n
    fl += 2 * (a ** (arb(-1) / 4)) * acb(a).gamma_upper(acb(arb(1) / 4)).real
mp4 = 2 * sum((mpmath.pi * n * n) ** (-mpmath.mpf(1) / 4) * mpmath.gammainc(mpmath.mpf(1) / 4, mpmath.pi * n * n) for n in range(1, 5))
P('first four terms: flint ball', S(fl, 40), 'mpmath', mpmath.nstr(mp4, 40))
# G(0) = 4 (1 - 2 xi(1/2)) in balls (the identity of Remark mi:rem:laws, which is the value of the series)
s = acb(arb(1) / 2)
xi = acb(arb(1) / 2) * s * (s - 1) * (acb(arb.pi()) ** (-s / 2)) * (s / 2).gamma() * s.zeta()
G = 4 * (1 - 2 * xi.real)
P('PX.1a G(0) = 4(1 - 2 xi(1/2)) in balls =', S(G, 30), 'radius', G.rad().str(3), '; series minus identity', mpmath.nstr(g12 - mpmath.mpf(S(G, 55)), 5))
eps = G / 4
# xi(1/2)
s = acb(arb(1) / 2)
xi = acb(arb(1) / 2) * s * (s - 1) * (acb(arb.pi()) ** (-s / 2)) * (s / 2).gamma() * s.zeta()
one_m = 1 - 2 * xi.real
P('PX.1b eps = G(0)/4 =', S(eps, 25), '; 1 - 2 xi(1/2) =', S(one_m, 25), '; difference', (eps - one_m).str(5), '(ball contains 0:', (eps - one_m).contains(0), ')')
psi1 = arb(0)
for n in range(1, 8):
    psi1 += (-arb.pi() * n * n).exp()
psi_tail = 2 * (-arb.pi() * 64).exp()
psi1 = arb(psi1.mid(), psi1.rad() + psi_tail.upper())
ser = arb(0)
for n in range(1, 8):
    a = arb.pi() * n * n
    ser += 2 * (2 * a - 1) * (-a).exp()
ser = arb(ser.mid(), ser.rad() + (8 * (-arb.pi() * 64).exp()).upper() + ser.rad())
half = arb(1) / 2 - psi1
P('PX.1c P(Y>1): series', S(ser, 25), '; 1/2 - psi(1)', S(half, 25), '; difference', (ser - half).str(5), '; psi(1) =', S(psi1, 25))
# quadrature of 2 Phi(u) exp(-u/2), u > 0
mpmath.mp.dps = 30


def Phi(u):
    tot = mpmath.mpf(0)
    e2 = mpmath.exp(2 * u)
    for n in range(1, 13):
        tot += (2 * mpmath.pi ** 2 * n ** 4 * mpmath.exp(4.5 * u) - 3 * mpmath.pi * n * n * mpmath.exp(2.5 * u)) * mpmath.exp(-mpmath.pi * n * n * e2)
    return 2 * tot


quad = mpmath.quad(lambda u: 2 * Phi(u) * mpmath.exp(-u / 2), [0, 0.5, 1, 2, 4])
P('quadrature of 2 Phi(u) exp(-u/2) on u>0 (mpmath, 30 digits):', mpmath.nstr(quad, 25), '; difference from the series', mpmath.nstr(quad - mpmath.mpf(S(ser, 35)), 5))
# extra sanity
P('E Y^2 check not run; E(Y-1)^+ = 2 psi(1) =', S(2 * psi1, 20))
out = {'G0': S(G, 30), 'G0_radius': G.rad().str(3), 'eps': S(eps, 30), 'one_minus_2xi': S(one_m, 30), 'eps_minus_oneminus2xi': (eps - one_m).str(5),
       'psi1': S(psi1, 30), 'PYgt1_series': S(ser, 30), 'half_minus_psi1': S(half, 30), 'series_minus_half': (ser - half).str(5),
       'quadrature': mpmath.nstr(quad, 28), 'quad_minus_series': mpmath.nstr(quad - mpmath.mpf(S(ser, 35)), 5), 'log': lines}
json.dump(out, open(args.out, 'w'), indent=1)
