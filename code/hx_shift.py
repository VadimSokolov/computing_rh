# code/hx_shift.py (October 2026, from misc/reviews/nick-b4/reconcile-P.md, task N1 and the lemma of R8, R12): the shifted Gauss Hermite evaluator of Lemma cp:lem:shift, written from the mathematics (the collaborator's scripts, notes and pickles were not read, run or copied).
# Definitions. xi(s) = g(s) zeta(s), g(s) = (1/2) s (s-1) pi^(-s/2) Gamma(s/2), l = g'/g = 1/s + 1/(s-1) - (1/2) log pi + (1/2) psi(s/2). For t >= 0: xi_t(s) = E xi(s + sqrt(t/2) Z) (Proposition xg:prop:gauss, Z standard normal), and by Lemma cp:lem:shift xi_t(s) = g(s) e^{t l(s)^2/4} zeta_t(s) with zeta_t(s) = E Psi(s, h), h = (t/2) l(s) + sqrt(t/2) Z, Psi(s,h) = xi(s+h)/(g(s) e^{h l(s)}). The expectation is a Gauss Hermite sum with K nodes: E f(Z) = sum_i c_i f(sqrt 2 x_i), c_i = w_i/sqrt(pi); so h_i = (t/2) l(s) + sqrt(t) x_i. Psi is evaluated as g(s+h) zeta(s+h)/(g(s) e^{h l}) from the closed forms (no logarithms, so no branch to follow); the ratio of the two Gamma factors is formed in extended exponent range arithmetic (mpmath, or Arb balls, whose exponents do not underflow: |xi_{0.2}(1/2+1000.3 i)| is about 3.7e-336, far below the double range), and zeta_t itself is of order one, so nothing is lost. The real function on the line is Z_t(theta) = e^{i(vartheta(theta) + (t/4) Im l(1/2+i theta)^2)} zeta_t(1/2+i theta), vartheta = Riemann Siegel theta (arg g(1/2+i theta) = vartheta + pi); it is real up to rounding, and its zeros are the zeros of xi_t(1/2+i theta).
# The library part (class Ev, gh_nodes, xi_trap) is imported by hx_scan.py, hx_far.py and hx_edge.py. Two backends with the same formulas: python-flint (Arb balls, midpoints used) when it is installed (Hopper), mpmath otherwise (laptop); the self test compares them.
# (1) The test of this script (task N1). At t = 0.2: the relative difference of zeta_t between K = 40 and K = 60 nodes (and K = 24 and 60) at theta = 1000.3 and 4000.3, each at 30, 45 and 70 working digits (the 70 digit run measures the quadrature error below the 30 digit rounding) ([N1a], [N1b]); the size of xi_{0.2}(1/2 + 1000.3 i); the comparison with the book's trapezoid evaluator of xi_t (Section sec:ch4:trapezoid with the weight e^{t u^2}, as in phase_rig.py, mpmath at 320 digits, step h = 1/800 and 1/1000 for the discretisation check (the kernel is analytic in |Im u| < pi/8, so the rule needs h < 2.47/(pi theta/4 + 560)), cut |u| <= 1.6) at theta = 50.3, 300.3, 600.3; xi_t(1/2) at t = 0 against xi(1/2); and the sign changes of Z_{0.2} at the four zeros sqrt(a_j(0.2)) of the table of xg_tilt.json (theta_j = z_j/2, at theta_j -+ 1e-6, K = 40 and 60).
# Error control. K = 24, 40, 60 and two working precisions are compared and the differences logged; the trapezoid rule is run with two steps.
# Output: hx_shift.json (values as strings) and hx_shift.log; python3 hx_shift.py [--quick] [--out FILE] [--backend flint|mpmath]. Quick: theta = 50.3 and 300.3, K = 24 and 40, trapezoid at 200 digits.
import argparse
import json
import math
import os
import sys
import time

import mpmath as mp

try:
    from flint import acb, arb, ctx
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False

T0 = time.time()
LOGF = None


def say(*a):
    msg = '[%7.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(msg, flush=True)
    if LOGF:
        with open(LOGF, 'a') as f:
            f.write(msg + '\n')


_GH = {}


def gh_nodes(K, dps=60):
    """Abscissae x_i and weights c_i = w_i/sqrt(pi) of the K point Gauss Hermite rule (physicists' weight e^{-x^2}), as mpf strings."""
    if (K, dps) in _GH:
        return _GH[(K, dps)]
    import numpy as np
    old = mp.mp.dps
    mp.mp.dps = dps + 15
    x0, _ = np.polynomial.hermite.hermgauss(K)
    out = []
    for xi in x0:
        x = mp.mpf(float(xi))
        for _ in range(7):
            p0, p1 = mp.mpf(1), 2 * x
            for n in range(1, K):
                p0, p1 = p1, 2 * x * p1 - 2 * n * p0
            x -= p1 / (2 * K * p0)
        p0, p1 = mp.mpf(1), 2 * x
        for n in range(1, K - 1):
            p0, p1 = p1, 2 * x * p1 - 2 * (n) * p0
        # now p1 = H_{K-1}(x)
        w = mp.mpf(2) ** (K - 1) * mp.factorial(K) * mp.sqrt(mp.pi) / (K * K * p1 ** 2)
        out.append((x, w / mp.sqrt(mp.pi)))
    s0 = sum(c for _, c in out)
    s2 = sum(c * 2 * x * x for x, c in out)
    s4 = sum(c * 4 * x ** 4 for x, c in out)
    assert abs(s0 - 1) < mp.mpf(10) ** (-dps) and abs(s2 - 1) < mp.mpf(10) ** (-dps) and abs(s4 - 3) < mp.mpf(10) ** (-dps + 3), (s0, s2, s4)
    res = [(mp.nstr(x, dps + 5, strip_zeros=False), mp.nstr(c, dps + 5, strip_zeros=False)) for x, c in out]
    mp.mp.dps = old
    _GH[(K, dps)] = res
    return res


class Ev:
    """Evaluator of zeta_t and Z_t at a working precision of dps digits, with Arb (flint) or mpmath."""

    def __init__(self, dps=30, backend=None):
        self.fl = HAVE_FLINT if backend is None else (backend == 'flint')
        if self.fl and not HAVE_FLINT:
            raise RuntimeError('python-flint is not installed')
        self.dps = dps
        self.cache = {}
        self.act()

    def act(self):
        if self.fl:
            ctx.prec = int(self.dps * 3.3219) + 24
        else:
            mp.mp.dps = self.dps + 4

    def real(self, x):
        return arb(str(x)) if self.fl else mp.mpf(str(x))

    def cx(self, re, im):
        return acb(self.real(re), self.real(im)) if self.fl else mp.mpc(self.real(re), self.real(im))

    def pi(self):
        return arb.pi() if self.fl else +mp.pi

    def g(self, s):
        if self.fl:
            return s * (s - 1) / 2 * (-s / 2 * self.pi().log()).exp() * (s / 2).gamma()
        return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2)

    def ell(self, s):
        if self.fl:
            return 1 / s + 1 / (s - 1) - self.pi().log() / 2 + (s / 2).digamma() / 2
        return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2

    def zeta(self, w):
        return w.zeta() if self.fl else mp.zeta(w)

    def exp(self, z):
        return z.exp() if self.fl else mp.exp(z)

    def nodes(self, K):
        key = (K, self.dps)
        if key not in self.cache:
            nd = gh_nodes(K, self.dps + 10)
            self.cache[key] = [(self.real(x), self.real(c)) for x, c in nd]
        return self.cache[key]

    def zt(self, s, t, K):
        """zeta_t(s) by the K node rule; s is a backend complex number."""
        self.act()
        if t == 0:
            return self.zeta(s)
        nd = self.nodes(K)
        tt = self.real(t)
        st = tt.sqrt() if self.fl else mp.sqrt(tt)
        ell = self.ell(s)
        gs = self.g(s)
        m = tt / 2 * ell
        acc = 0
        for x, c in nd:
            h = m + st * x
            w = s + h
            acc += c * self.g(w) * self.zeta(w) / (gs * self.exp(h * ell))
        return acc

    def theta_rs(self, theta):
        """Riemann Siegel theta at the real number theta (python number or string)."""
        if self.fl:
            th = self.real(theta)
            return acb(0.25, th / 2).lgamma().imag - th / 2 * self.pi().log()
        return mp.siegeltheta(self.real(theta))

    def Z(self, theta, t, K):
        """The real function Z_t(theta) (backend complex, imaginary part is rounding)."""
        self.act()
        s = self.cx('0.5', theta)
        z = self.zt(s, t, K)
        e2 = self.ell(s) ** 2
        ph = self.theta_rs(theta) + self.real(t) / 4 * e2.imag
        e = acb(0, ph).exp() if self.fl else mp.exp(1j * ph)
        return z * e

    def to_mp(self, z):
        """A backend number as an mpmath number at the current mpmath precision (midpoint)."""
        if self.fl:
            def f(a):
                return mp.mpf(a.mid().str(self.dps + 5, radius=False))
            return mp.mpc(f(z.real), f(z.imag)) if isinstance(z, acb) else f(z)
        return z

    def relrad(self, z):
        """Relative radius of the Arb ball (0 for mpmath)."""
        if not self.fl:
            return 0.0
        return float((z.real.rad() + z.imag.rad()) / abs(z).mid()) if isinstance(z, acb) else 0.0

    def to_c(self, z):
        return complex(self.to_mp(z))


# ---------------------------------------------------------------- the book's trapezoid evaluator of xi_t
def phi_dbn(u, dps):
    """Phi_dBN(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}), u >= 0."""
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    stop = mp.mpf(10) ** (-(dps + 8))
    tot = mp.mpf(0)
    n = 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        tot += term
        if n > 2 and abs(term) < stop * abs(tot) + mp.mpf(10) ** (-3 * dps):
            return tot
        n += 1


_TRAP = {}


def trap_setup(t, h_inv, umax, dps):
    key = (str(t), h_inv, str(umax), dps)
    if key not in _TRAP:
        old = mp.mp.dps
        mp.mp.dps = dps
        tt = mp.mpf(str(t))
        J = int(mp.ceil(mp.mpf(str(umax)) * h_inv))
        h = mp.mpf(1) / h_inv
        c = [mp.exp(tt * (j * h) ** 2) * phi_dbn(j * h, dps) for j in range(J + 1)]
        _TRAP[key] = (h, J, c)
        mp.mp.dps = old
    return _TRAP[key]


def xi_trap(theta, t, h_inv, umax, dps):
    """xi_t(1/2 + i theta) = 8 int_0^oo e^{t u^2} Phi_dBN(u) cos(2 theta u) du by the trapezoid rule (mpmath, dps digits)."""
    h, J, c = trap_setup(t, h_inv, umax, dps)
    old = mp.mp.dps
    mp.mp.dps = dps
    th = mp.mpf(str(theta))
    tot = c[0] / 2
    for j in range(1, J + 1):
        tot += c[j] * mp.cos(2 * th * j * h)
    res = 8 * h * tot
    mp.mp.dps = old
    return res


def xi_exact(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def reldiff(a, b):
    return abs(a - b) / abs(b)


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='hx_shift.json')
    ap.add_argument('--backend', default=None)
    ap.add_argument('--selftest', action='store_true', help='compare the two backends (needs python-flint)')
    a = ap.parse_args()
    LOGF = os.path.splitext(a.out)[0] + '.log'
    open(LOGF, 'w').close()
    t = '0.2'
    res = {'script': 'hx_shift.py', 'quick': a.quick, 'backend': ('flint' if (HAVE_FLINT if a.backend is None else a.backend == 'flint') else 'mpmath')}
    say('hx_shift.py, backend', res['backend'], 'quick' if a.quick else 'full')
    mp.mp.dps = 40
    if a.selftest:
        e1 = Ev(30, 'flint')
        e2 = Ev(30, 'mpmath')
        for th in ['14.7', '50.3', '300.3', '1000.3']:
            z1 = e1.to_mp(e1.zt(e1.cx('0.5', th), t, 40))
            z2 = e2.zt(e2.cx('0.5', th), t, 40)
            mp.mp.dps = 30
            say('selftest backends at theta', th, 'rel diff', mp.nstr(reldiff(z1, z2), 5))
    # (a) xi_t(1/2) at t = 0 against xi(1/2), and Psi at t = 0.2, theta = 14 against the trapezoid
    mp.mp.dps = 40
    x0 = xi_trap('0', '0', 200, '3', 60)
    xh = xi_exact(mp.mpf(1) / 2)
    say('trapezoid xi_0(1/2) =', mp.nstr(x0, 30), ' xi(1/2) =', mp.nstr(xh, 30), ' rel diff', mp.nstr(reldiff(x0, xh), 5))
    res['xi0_half_rel_diff'] = mp.nstr(reldiff(x0, xh), 5)
    # (b) K comparisons
    thetas = ['50.3', '300.3'] if a.quick else ['1000.3', '4000.3']
    Ks = (24, 40) if a.quick else (24, 40, 60)
    res['K_compare'] = {}
    for th in thetas:
        for dps in (30, 45, 70):
            ev = Ev(dps, a.backend)
            s = ev.cx('0.5', th)
            raw = {K: ev.zt(s, t, K) for K in Ks}
            rr = {K: ev.relrad(raw[K]) for K in Ks}
            val = {K: ev.to_mp(raw[K]) for K in Ks}
            mp.mp.dps = dps + 10
            val = {K: ev.to_mp(raw[K]) for K in Ks}
            kk = Ks[-1]
            d40 = reldiff(val[40], val[kk]) if 40 != kk else None
            d24 = reldiff(val[24], val[kk])
            mp.mp.dps = 20
            say('theta', th, 'dps', dps, 'Arb relative radii', {K: '%.1e' % rr[K] for K in Ks}, 'zeta_t =', mp.nstr(val[kk], 15), ' rel diff K=24 vs %d:' % kk, mp.nstr(d24, 4), ' K=40 vs %d:' % kk, (mp.nstr(d40, 4) if d40 is not None else '-'))
            res['K_compare']['theta=%s dps=%d' % (th, dps)] = {'zeta_t': mp.nstr(val[kk], 25), 'arb_rel_radius': {str(K): rr[K] for K in Ks}, 'rel_24_vs_%d' % kk: mp.nstr(d24, 6), 'rel_40_vs_%d' % kk: (mp.nstr(d40, 6) if d40 is not None else None)}
            if dps == 30:
                ev.act()
                mp.mp.dps = 30
                lg = ev.to_mp(ev.g(s))
                el = ev.to_mp(ev.ell(s))
                xi_mod = abs(lg) * mp.exp(mp.mpf(t) / 4 * (el * el).real) * abs(val[kk])
                say('   |xi_0.2(1/2+i theta)| =', mp.nstr(xi_mod, 6))
                res['K_compare']['theta=%s dps=30' % th]['abs_xi_t'] = mp.nstr(xi_mod, 6)
    # (c) the trapezoid evaluator
    res['trapezoid'] = {}
    tdps = 200 if a.quick else 320
    for th in (['50.3', '300.3'] if a.quick else ['50.3', '300.3', '600.3']):
        ev = Ev(30, a.backend)
        s = ev.cx('0.5', th)
        mp.mp.dps = 40
        zk = ev.to_mp(ev.zt(s, t, 60 if not a.quick else 40))
        mp.mp.dps = 30
        ls = mp.mpc(0.5, mp.mpf(th))
        pref = xi_exact_g = ls * (ls - 1) / 2 * mp.pi ** (-ls / 2) * mp.gamma(ls / 2)
        ell = 1 / ls + 1 / (ls - 1) - mp.log(mp.pi) / 2 + mp.digamma(ls / 2) / 2
        den = pref * mp.exp(mp.mpf(t) / 4 * ell * ell)
        hi1, hi2 = (150, 200) if a.quick else (800, 1000)
        um = '1.6'
        xa = xi_trap(th, t, hi1, um, tdps)
        xb = xi_trap(th, t, hi2, um, tdps)
        mp.mp.dps = 30
        za = xa / den
        zb = xb / den
        say('theta', th, ': trapezoid xi_t = %s (h=1/%d), %s (h=1/%d); zeta_t(trap) =' % (mp.nstr(xa, 12), hi1, mp.nstr(xb, 12), hi2), mp.nstr(zb, 15),
            ' GH =', mp.nstr(zk, 15), ' rel diff trap vs GH', mp.nstr(reldiff(zb, zk), 4), ' trap h vs h2', mp.nstr(reldiff(za, zb), 4))
        res['trapezoid'][th] = {'zeta_t_trap': mp.nstr(zb, 20), 'zeta_t_GH': mp.nstr(zk, 20), 'rel_trap_vs_GH': mp.nstr(reldiff(zb, zk), 6), 'rel_step_check': mp.nstr(reldiff(za, zb), 6), 'abs_xi_t': mp.nstr(abs(xb), 6)}
    # (d) sign changes at the four zeros of the xg_tilt table
    tab = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'xg_tilt.json'))) if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'xg_tilt.json')) else json.load(open('xg_tilt.json'))
    zs = tab['table']['0.2']['zeros_z']
    res['sign_changes'] = []
    ev = Ev(25, a.backend)
    for j, z in enumerate(zs):
        mp.mp.dps = 30
        th = mp.mpf(z) / 2
        row = {'j': j + 1, 'theta_j': mp.nstr(th, 15)}
        for K in ((40,) if a.quick else (40, 60)):
            zl = ev.to_c(ev.Z(mp.nstr(th - mp.mpf('1e-6'), 25), t, K))
            zr = ev.to_c(ev.Z(mp.nstr(th + mp.mpf('1e-6'), 25), t, K))
            row['K%d' % K] = {'Z_left': '%.6e' % zl.real, 'Z_right': '%.6e' % zr.real, 'max_imag_over_real': '%.2e' % (max(abs(zl.imag) / abs(zl.real), abs(zr.imag) / abs(zr.real))), 'sign_change': bool(zl.real * zr.real < 0)}
        say('zero', j + 1, 'theta', row['theta_j'], row)
        res['sign_changes'].append(row)
    res['seconds'] = round(time.time() - T0, 1)
    json.dump(res, open(a.out, 'w'), indent=1)
    say('wrote', a.out)


if __name__ == '__main__':
    main()
