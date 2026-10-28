# tlz_deform.py: TLZ.3 of reconcile-TLZ.md. First order motion of the zeros of D (the Mellin transform of the clock of the pair Y_eps + Y_eps') at eps = 0 for the
# density (1 + eps t^2)/(2 pi + (2/3) pi^3 eps) on [-pi, pi] (fm:res:deform).
# D_eps(w) = int_0^inf u^(w-1) (L_eps(sqrt u)^(-2) - 1) du for -1 < Re w < 0 (continuation of int u^(w-1) L^(-2) du); with u = e^(2t):  D_eps(w) = 2 int e^(2wt) (L_eps(e^t)^(-2) - 1) dt.
# L_eps = (L + eps M2)/(1 + eps pi^2/3), L(s) = sinh(pi s)/(pi s), M2(s) = (1/2pi) int t^2 e^(st) dt, Ldot = M2 - (pi^2/3) L.
# dD/deps at 0:  Ddot(w) = -4 int e^(2wt) Ldot(e^t)/L(e^t)^3 dt. The zero w_k(0) = -1/4 + i gamma_k/2 moves with dw_k/deps = -Ddot(w_k)/D0'(w_k),
# D0(w) = -2 pi (2w+1) zeta(-2w)/sin(pi w), D0' = 4 pi (2w+1) zeta'(-2w)/sin(pi w) at a zero.
# Quadrature: composite Gauss Legendre on t in [T, log 16] (30 nodes per panel of width 1/2, then 40 nodes per panel of width 1/4), dps 60; below T = -30 the integrand is
# summed from its Taylor series in u = e^(2t) (10 terms, closed form integrals). Checks: D_0 by quadrature against the closed form at -0.3+4i; the lower end T = -30 against -40;
# central difference of the zeros of the quadrature D_eps at eps = +-1e-3 and +-1e-4 (secant); the zeros at eps = 0.002.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import json, argparse, time
import mpmath as mp

J = 10
def smul(a, b):
    c = [mp.mpf(0)] * J
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < J: c[i + j] += x * y
    return c
def sinv(a):
    b = [mp.mpf(0)] * J; b[0] = 1 / a[0]
    for n in range(1, J):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1)) / a[0]
    return b

def series_pieces(eps):
    pi = mp.pi
    Ls = [pi ** (2 * j) / mp.factorial(2 * j + 1) for j in range(J)]
    M2 = [pi ** (2 * j + 2) / (mp.factorial(2 * j) * (2 * j + 3)) for j in range(J)]
    Ld = [M2[j] - (pi ** 2 / 3) * Ls[j] for j in range(J)]
    Le = [(Ls[j] + eps * M2[j]) / (1 + eps * pi ** 2 / 3) for j in range(J)]
    i2 = sinv(smul(Le, Le))               # 1/L_eps^2
    g = [i2[0] - 1] + i2[1:]
    h = smul(Ld, smul(sinv(Ls), smul(sinv(Ls), sinv(Ls))))   # Ldot/L^3
    return g, h

def LM(s):
    """L(s), M2(s) at s (mp), series for s < 0.5, closed form otherwise"""
    pi = mp.pi
    if s < 0.5:
        s2 = s * s
        L = sum(pi ** (2 * j) / mp.factorial(2 * j + 1) * s2 ** j for j in range(45))
        M = sum(pi ** (2 * j + 2) / (mp.factorial(2 * j) * (2 * j + 3)) * s2 ** j for j in range(45))
        return L, M
    sh, ch = mp.sinh(pi * s), mp.cosh(pi * s)
    return sh / (pi * s), (2 * pi ** 2 * sh / s - 4 * pi * ch / s ** 2 + 4 * sh / s ** 3) / (2 * pi)

class Quad:
    def __init__(self, T, width, n, tmax):
        self.T, self.tmax = T, tmax
        # Gauss Legendre nodes on [-1,1]: Newton on the Legendre polynomial from the numpy nodes
        import numpy as np
        x0, _ = np.polynomial.legendre.leggauss(n)
        X, W = [], []
        for x in x0:
            r = mp.findroot(lambda y: mp.legendre(n, y), mp.mpf(float(x)), tol=mp.mpf(10) ** (-mp.mp.dps + 5))
            dp = n * (r * mp.legendre(n, r) - mp.legendre(n - 1, r)) / (r * r - 1)
            X.append(r); W.append(2 / ((1 - r * r) * dp * dp))
        self.t, self.wt = [], []
        npan = int(mp.ceil((tmax - T) / width))
        a = mp.mpf(T); h = (mp.mpf(tmax) - a) / npan
        for p in range(npan):
            lo = a + p * h
            for x, w in zip(X, W):
                self.t.append(lo + h * (x + 1) / 2); self.wt.append(w * h / 2)
        # values of L, M2 at the nodes
        self.L, self.M = [], []
        for t in self.t:
            L, M = LM(mp.e ** t); self.L.append(L); self.M.append(M)
        self.cache = {}
    def gvals(self, eps):
        k = ('g', eps)
        if k not in self.cache:
            c = 1 + eps * mp.pi ** 2 / 3
            self.cache[k] = [1 / ((L + eps * M) / c) ** 2 - 1 for L, M in zip(self.L, self.M)]
        return self.cache[k]
    def hvals(self):
        if 'h' not in self.cache:
            self.cache['h'] = [(M - mp.pi ** 2 / 3 * L) / L ** 3 for L, M in zip(self.L, self.M)]
        return self.cache['h']
    def D(self, w, eps):
        g = self.gvals(eps)
        s = mp.fsum(wt * mp.exp(2 * w * t) * gv for t, wt, gv in zip(self.t, self.wt, g))
        gs, _ = series_pieces(eps)
        tail = mp.fsum(gs[j] * mp.exp((2 * w + 2 * j) * self.T) / (2 * w + 2 * j) for j in range(1, J))
        up = mp.exp(2 * w * self.tmax) / (2 * w)   # tail of the term -1 on [tmax, inf), Re w < 0 (L^-2 is negligible there)
        return 2 * (s + tail + up)
    def Ddot(self, w):
        hv = self.hvals()
        s = mp.fsum(wt * mp.exp(2 * w * t) * hh for t, wt, hh in zip(self.t, self.wt, hv))
        _, hs = series_pieces(0)
        tail = mp.fsum(hs[j] * mp.exp((2 * w + 2 * j) * self.T) / (2 * w + 2 * j) for j in range(1, J))
        return -4 * (s + tail)

def D0cf(w):
    return -2 * mp.pi * (2 * w + 1) * mp.zeta(-2 * w) / mp.sin(mp.pi * w)
def D0p(w):
    return 4 * mp.pi * (2 * w + 1) * mp.zeta(-2 * w, 1, 1) / mp.sin(mp.pi * w)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default='tlz_deform.json'); ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    mp.mp.dps = 40 if a.quick else 60
    tmax = mp.log(16)
    Q1 = Quad(-30, mp.mpf(1) / 2, 30 if not a.quick else 16, tmax)
    out = {'dps': mp.mp.dps}
    wc = mp.mpc(-0.3, 4)
    out['control_closed'] = mp.nstr(D0cf(wc), 12); out['control_quad'] = mp.nstr(Q1.D(wc, 0), 12)
    print('control', out['control_closed'], out['control_quad'], flush=True)
    ks = [1, 2, 3]
    gam = {k: mp.zetazero(k).imag for k in ks}
    der = {}
    Q2 = Quad(-30, mp.mpf(1) / 4, 40 if not a.quick else 20, tmax) if not a.quick else Q1
    Q3 = Quad(-40, mp.mpf(1) / 2, 30 if not a.quick else 16, tmax) if not a.quick else Q1
    for k in ks:
        w = mp.mpc(-0.25, gam[k] / 2)
        d = {'gamma': mp.nstr(gam[k], 20), 'D0_at_w': mp.nstr(abs(D0cf(w)), 3), 'D0prime': mp.nstr(D0p(w), 20),
             'D0prime_numeric_check': mp.nstr(mp.diff(D0cf, w), 15)}
        vals = []
        for name, Q in (('rule1', Q1), ('rule2', Q2), ('rule3_T-40', Q3)):
            dd = Q.Ddot(w); dwde = -dd / D0p(w)
            vals.append(dwde); d['Ddot_' + name] = mp.nstr(dd, 25); d['dw_' + name] = mp.nstr(dwde, 25)
        d['rule_diff'] = mp.nstr(max(abs(vals[0] - vals[1]), abs(vals[0] - vals[2])), 3)
        d['print'] = '%s %s %si' % (mp.nstr(vals[0].real, 10), '+' if vals[0].imag >= 0 else '-', mp.nstr(abs(vals[0].imag), 10))
        # D_eps at the zero for the control: quadrature D_0(w) small
        d['quadD0_at_w'] = mp.nstr(abs(Q1.D(w, 0)), 3)
        der[k] = vals[0]
        out['k%d' % k] = d
        print('k', k, json.dumps(d), flush=True)
    # central difference of the zeros
    for k in (ks if not a.quick else [1]):
        w0 = mp.mpc(-0.25, gam[k] / 2)
        zs = {}
        for e in (mp.mpf(10) ** -3, -mp.mpf(10) ** -3, mp.mpf(10) ** -4, -mp.mpf(10) ** -4, mp.mpf(2) / 1000):
            guess = w0 + e * der[k]
            z = mp.findroot(lambda w: Q1.D(w, e), (guess, guess + mp.mpf(10) ** -6), solver='secant', tol=mp.mpf(10) ** (-8 if a.quick else -35), maxsteps=40)
            zs[str(e)] = z
        cd = (zs[str(mp.mpf(10) ** -3)] - zs[str(-mp.mpf(10) ** -3)]) / (2 * mp.mpf(10) ** -3)
        cd4 = (zs[str(mp.mpf(10) ** -4)] - zs[str(-mp.mpf(10) ** -4)]) / (2 * mp.mpf(10) ** -4)
        out['k%d' % k]['central_difference_h1e-4'] = mp.nstr(cd4, 12); out['k%d' % k]['central_h1e-4_minus_formula'] = mp.nstr(abs(cd4 - der[k]), 3)
        out['k%d' % k]['central_difference'] = mp.nstr(cd, 12)
        out['k%d' % k]['central_minus_formula'] = mp.nstr(abs(cd - der[k]), 3)
        out['k%d' % k]['zero_eps_0.002'] = mp.nstr(zs[str(mp.mpf(2) / 1000)], 12)
        print('central', k, out['k%d' % k]['central_difference'], out['k%d' % k]['central_minus_formula'], out['k%d' % k]['zero_eps_0.002'], flush=True)
    json.dump(out, open(a.out, 'w'), indent=1)

if __name__ == '__main__':
    main()
