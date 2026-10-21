# polya_mellin_z1.py: the zero z_1 = s_1 + 1/2 of the Mellin transform M of the Riemann kernel, Chapter 16 (Polya and the Zeros of Entire
# Functions), Section pol:sec:mellin. With Phi(u) = 2 sum_n (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)) and
# M(z) = int_0^oo Phi(u) u^(2z-1) du (meromorphic continuation to Re z < 0, poles at z = 0, -1, -2, ...), Theorem thm:ch13:mellin says that
# M(s+1/2) has a simple zero near s_1 = -8.494182256992559 + 2.937767068567584 i (an interval certificate of the first author's, proved there).
# This script checks it from the mathematics, at high precision: it evaluates M and M' at z_0 = s_1 + 1/2, runs Newton's method from z_0 to the
# zero z_1 of M, and records |z_1 - z_0|, |M(z_0)|, |M'(z_1)|, the minimum of |M| on the circle |z - z_1| = 1/2 and the winding number of M
# around that circle (the number of zeros of M in the disc, since M has no pole there: the poles are on the real axis).
# Method. M(z) = sum_j p_j d^(2z+2j)/(2z+2j) + int_d^2.1 Phi(u) u^(2z-1) du with d = 1/4: Phi is even and analytic in |Im u| < pi/4, so
# Phi(u) = sum_j p_j u^(2j) with radius of convergence pi/4, and the piece on (0, d) is integrated termwise (J = 90 terms; the last term is
# below 1e-60); p_j come from the Cauchy integral on |u| = 1/2 with N points (aliasing about (0.5/(pi/4))^N). The rest is Gauss Legendre of 48
# nodes on each of seven panels of ratio at most 3/2 (the integrand u^(2z-1) is singular only at u = 0, so the rule converges geometrically);
# Phi(2.1) is below 1e-90. At Re z = -8 the two parts are of size 1e9 and M is of size 1e-14, so the working precision is 60 digits (40 digits
# lose 24). Error control. A second layout is evaluated at z_0 and z_1: d = 0.35 (J = 140 terms), other panels, and tanh sinh quadrature
# (mpmath quad) in place of Gauss Legendre; the relative difference of M and M' between the layouts is recorded. Nothing here is a proof:
# the computation is floating point at the precision stated (60 digits), and M' and the winding number are computed at 256 points.
# Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were neither read nor run).
# Output: data/polya_mellin_z1.json and data/polya_mellin_z1.log next to --out (the log has the same name with the extension .log).
# Usage: python3 polya_mellin_z1.py [--quick] [--out FILE]; --quick uses 40 digits, N = 192 and 64 circle points (about 10 s on one core) and
# writes to /tmp unless --out is given; the full case is for the Hopper cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g polya_mellin_z1.json -g polya_mellin_z1.log code/polya_mellin_z1.py -- --out polya_mellin_z1.json
import argparse
import atexit
import json
import time
from pathlib import Path

import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

T0 = time.time()
LOGF = None

# s_1 of Theorem thm:ch13:mellin (as printed there, 16 digits) and z_0 = s_1 + 1/2
S1 = ('-8.494182256992559', '2.937767068567584')


def log(*a):
    line = '[%7.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(line, flush=True)
    if LOGF is not None:
        LOGF.write(line + '\n')
        LOGF.flush()


def fs(x, n=12):
    return mp.nstr(mp.mpf(x), n)


def fz(z, n=12):
    z = mp.mpmathify(z)
    return '%s%s%si' % (mp.nstr(z.real, n), '+' if z.imag >= 0 else '-', mp.nstr(abs(z.imag), n))


def phi_mp(u):
    u = mp.mpmathify(u)
    e2 = mp.exp(2 * u)
    s = 0
    for n in range(1, 400):
        term = (2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * u / 2)) * mp.exp(-mp.pi * n ** 2 * e2)
        s += term
        if n > 2 and abs(term) <= mp.mpf(10) ** (-mp.mp.dps - 8) * abs(s):
            break
    return 2 * s


def taylor_p(J, N, r=mp.mpf(1) / 2):
    """p_j, j < J, with Phi(u) = sum p_j u^(2j), from the Cauchy integral on |u| = r with N points (N a multiple of 4)."""
    vals = [phi_mp(r * mp.expj(2 * mp.pi * m / N)) for m in range(N)]
    out = []
    for j in range(J):
        k = 2 * j
        s = sum(vals[m] * mp.expj(-2 * mp.pi * m * k / N) for m in range(N)) / N
        out.append(mp.re(s) / r ** k)
    return out


class Layout:
    """M and M' by the termwise formula on (0, d) and a quadrature rule on (d, 2.1)."""

    def __init__(self, d, p, panels, degree, rule):
        self.d = mp.mpf(d)
        self.p = p
        self.rule = rule
        self.panels = [mp.mpf(x) for x in panels]
        if rule == 'gauss':
            xs = GaussLegendre(mp.mp).calc_nodes(degree, mp.mp.prec)
            self.nodes = []
            for a, b in zip(self.panels[:-1], self.panels[1:]):
                for x, w in xs:
                    u = (a + b) / 2 + (b - a) / 2 * x
                    self.nodes.append((u, w * (b - a) / 2, phi_mp(u), mp.log(u)))
        self.n_terms = len(p)

    def termwise(self, z):
        d = self.d
        lz = mp.log(d)
        base = mp.power(d, 2 * z)
        a = mp.mpf(0)
        da = mp.mpf(0)
        dj = mp.mpf(1)
        for j, pj in enumerate(self.p):
            c = 2 * z + 2 * j
            t = pj * base * dj / c
            a += t
            da += t * (2 * lz - 2 / c)
            dj *= d * d
        return a, da

    def tail(self, z):
        if self.rule == 'gauss':
            q = mp.mpf(0)
            dq = mp.mpf(0)
            for u, w, ph, lu in self.nodes:
                t = w * ph * mp.exp((2 * z - 1) * lu)
                q += t
                dq += t * 2 * lu
            return q, dq
        f = lambda u: phi_mp(u) * mp.power(u, 2 * z - 1)
        g = lambda u: phi_mp(u) * mp.power(u, 2 * z - 1) * 2 * mp.log(u)
        q = mp.quad(f, self.panels)
        dq = mp.quad(g, self.panels)
        return q, dq

    def M(self, z):
        a, da = self.termwise(z)
        q, dq = self.tail(z)
        return a + q, da + dq

    def M0(self, z):
        a, _ = self.termwise(z)
        if self.rule == 'gauss':
            q = mp.mpf(0)
            for u, w, ph, lu in self.nodes:
                q += w * ph * mp.exp((2 * z - 1) * lu)
        else:
            q = mp.quad(lambda u: phi_mp(u) * mp.power(u, 2 * z - 1), self.panels)
        return a + q


def newton(L, z, iters=30):
    for it in range(iters):
        m, dm = L.M(z)
        step = m / dm
        z = z - step
        if abs(step) < mp.mpf(10) ** (-(mp.mp.dps - 12)):
            break
    return z, it + 1


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    out_path = Path(args.out) if args.out else (Path('/tmp/polya_mellin_z1.json') if args.quick else Path(__file__).resolve().parent.parent / 'data' / 'polya_mellin_z1.json')
    LOGF = open(out_path.with_suffix('.log'), 'w')
    atexit.register(LOGF.close)
    dps = 40 if args.quick else 60
    N = 192 if args.quick else 512
    ncirc = 64 if args.quick else 256
    mp.mp.dps = dps
    R = {'script': 'code/polya_mellin_z1.py', 'mode': 'quick' if args.quick else 'full', 'mpmath': mp.__version__,
         'parameters': {'dps': dps, 'cauchy_points': N, 'circle_points': ncirc, 'delta_A': '1/4', 'terms_A': 90, 'delta_B': '0.35', 'terms_B': 140}}
    z0 = mp.mpc(mp.mpf(S1[0]) + mp.mpf(1) / 2, mp.mpf(S1[1]))
    R['z0_from_theorem'] = fz(z0, 16)
    log('z_0 = s_1 + 1/2 =', fz(z0, 16), 'working precision', dps, 'digits')

    # layout A: d = 1/4, Gauss Legendre on seven panels of ratio at most 3/2 (48 nodes each)
    pA = taylor_p(90, N)
    log('Taylor coefficients (90) done; p_0 =', fs(pA[0], 15), 'p_1 =', fs(pA[1], 15), 'last term at d = 1/4:', fs(abs(pA[89]) * mp.mpf(0.25) ** 178, 3))
    R['p0'] = fs(pA[0], 15)
    R['p1'] = fs(pA[1], 15)
    R['last_term_A'] = fs(abs(pA[89]) * mp.mpf(0.25) ** 178, 3)
    pan = [mp.mpf(1) / 4]
    while pan[-1] < mp.mpf(21) / 10:
        pan.append(min(pan[-1] * mp.mpf(3) / 2, mp.mpf(21) / 10))
    LA = Layout(mp.mpf(1) / 4, pA, pan, 5, 'gauss')
    log('layout A: %d nodes on %d panels' % (len(LA.nodes), len(pan) - 1))
    m0, dm0 = LA.M(z0)
    log('M(z_0) =', fz(m0, 6), ' |M(z_0)| =', fs(abs(m0), 6), ' |M\'(z_0)| =', fs(abs(dm0), 8))
    R['abs_M_at_z0'] = fs(abs(m0), 6)
    z1, its = newton(LA, z0)
    m1, dm1 = LA.M(z1)
    R['z1'] = fz(z1, 24)
    R['distance_z1_z0'] = fs(abs(z1 - z0), 4)
    R['newton_iterations'] = its
    R['abs_M_at_z1'] = fs(abs(m1), 4)
    R['abs_Mprime_at_z1'] = fs(abs(dm1), 12)
    log('Newton from z_0: z_1 =', fz(z1, 24), ' after', its, 'steps; |z_1 - z_0| =', fs(abs(z1 - z0), 4), ' |M(z_1)| =', fs(abs(m1), 4), ' |M\'(z_1)| =', fs(abs(dm1), 12))

    # the circle |z - z_1| = 1/2: minimum of |M| and winding number
    mn = None
    wind = mp.mpf(0)
    prev = None
    first = None
    for k in range(ncirc + 1):
        z = z1 + mp.mpf(1) / 2 * mp.expj(2 * mp.pi * (k % ncirc) / ncirc)
        m = LA.M0(z)
        a = abs(m)
        if mn is None or a < mn:
            mn = a
        if prev is not None:
            dang = mp.arg(m / prev)
            wind += dang
        prev = m
    R['circle_min_abs_M'] = fs(mn, 6)
    R['circle_winding_number'] = fs(wind / (2 * mp.pi), 8)
    log('circle |z - z_1| = 1/2 with %d points: min |M| =' % ncirc, fs(mn, 6), ' winding number =', fs(wind / (2 * mp.pi), 8))

    # layout B: d = 0.35, other panels, tanh sinh quadrature
    pB = taylor_p(140, N)
    LB = Layout(mp.mpf('0.35'), pB, ['0.35', '0.7', '1.4', '2.1'], None, 'quad')
    mB0, dmB0 = LB.M(z0)
    mB1, dmB1 = LB.M(z1)
    R['layout_B'] = {'abs_M_at_z0': fs(abs(mB0), 6), 'abs_M_at_z1': fs(abs(mB1), 4), 'abs_Mprime_at_z1': fs(abs(dmB1), 12),
                     'relative_difference_M_at_z0': fs(abs(mB0 - m0) / abs(m0), 3), 'relative_difference_Mprime_at_z1': fs(abs(dmB1 - dm1) / abs(dm1), 3)}
    log('layout B (d = 0.35, tanh sinh): |M(z_0)| =', fs(abs(mB0), 6), ' |M(z_1)| =', fs(abs(mB1), 4), ' |M\'(z_1)| =', fs(abs(dmB1), 12),
        ' relative differences:', R['layout_B']['relative_difference_M_at_z0'], R['layout_B']['relative_difference_Mprime_at_z1'])
    R['seconds'] = fs(time.time() - T0, 4)
    out_path.write_text(json.dumps(R, indent=1))
    log('written', out_path)


if __name__ == '__main__':
    main()
