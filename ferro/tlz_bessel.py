# tlz_bessel.py: TLZ.4 of reconcile-TLZ.md. The Mellin transform D of the clock of the pair of two copies of the Bessel law q with (a, ell) = (pi, 4), m = 1,
# phi_q(theta) = K_{i theta/4}(pi)/K_0(pi) (fm:res:besselmellin), of a single copy, and (control) of the uniform pair pi U + pi U'. Real form: L(x) = E e^{xY} = K_{x/4}(pi)/K_0(pi)
# (Bessel) or sinh(pi x)/(pi x) (uniform), x = sqrt u. D(w) = int_0^inf u^(w-1) (L(sqrt u)^(-p) - 1) du for -1 < Re w < 0 (p = 2 pair, p = 1 single), v = log u:
# D(w) = int e^(wv) (f(v) - 1) dv, f = L^(-p).
# Quadrature: composite Gauss Legendre in v on [-14, 10] (u_1 = e^10 >= 10 as the plan requires; f is below 1e-80 at the upper end for the Bessel laws, checked and logged),
# 24 nodes per panel of width 1/2 (rule 1) and 32 nodes per panel of width 1/4 (rule 2); f computed with mpmath (besselk of real order, 40 digits) and stored in double
# precision. Below v = -14: the Taylor series of f in u = e^v from the moments m_2j = 4^(-2j) int r^2j e^(-pi cosh r) dr / K_0(pi) (8 terms, closed form integrals);
# above v = 10: the tail of -1, e^(w v1)/w. Control: the closed form -2 pi (2w+1) zeta(-2w)/sin(pi w) at -0.3+4i and at the uniform zero -1/4 + i gamma_1/2.
# Argument principle: boxes -0.95 < Re w < -0.05 with 1/2 < Im w < 17/2 and 17/2 < Im w < 11, adaptive sampling of the boundary (each phase step below 0.3, at least 2000
# points per side, then refined): winding number, min |D| on the contour, largest phase step (computed, sampled: not a ball proof). Rule 2 gives the same.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import json, argparse
import numpy as np
import mpmath as mp

V0, V1 = -14.0, 10.0
J = 9

def smul(a, b):
    c = [mp.mpf(0)] * J
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < J: c[i + j] += x * y
    return c
def spow_inv(a, p):
    """series of a^(-p) for a[0] = 1"""
    b = [mp.mpf(0)] * J; b[0] = mp.mpf(1)
    # (1+z)^(-p) with z = a - 1
    z = [mp.mpf(0)] + a[1:]
    zk = [mp.mpf(1)] + [mp.mpf(0)] * (J - 1)
    for k in range(1, J):
        zk = smul(zk, z)
        c = mp.binomial(-p, k)
        b = [bi + c * zi for bi, zi in zip(b, zk)]
    return b

def Lseries(kind):
    if kind == 'uniform':
        return [mp.pi ** (2 * j) / mp.factorial(2 * j + 1) for j in range(J)]
    k0 = mp.besselk(0, mp.pi)
    out = []
    for j in range(J):
        I = mp.quad(lambda r: r ** (2 * j) * mp.exp(-mp.pi * mp.cosh(r)), [0, 1, 2, 4, 8, 16, 30])
        out.append(I / k0 / mp.mpf(4) ** (2 * j) / mp.factorial(2 * j))
    return out

def Lfun(kind, x):
    if kind == 'uniform': return mp.sinh(mp.pi * x) / (mp.pi * x)
    return mp.besselk(x / 4, mp.pi) / mp.besselk(0, mp.pi)

class Rule:
    def __init__(self, kind, p, width, n):
        x, wgt = np.polynomial.legendre.leggauss(n)
        npan = int(round((V1 - V0) / width))
        self.v = np.concatenate([V0 + width * (i + (x + 1) / 2) for i in range(npan)])
        self.wt = np.concatenate([wgt * width / 2 for i in range(npan)])
        fm1 = []
        fmin_top = None
        for vi in self.v:
            L = Lfun(kind, mp.exp(mp.mpf(vi) / 2))
            fm1.append(float(L ** (-p) - 1))
        self.fm1 = np.array(fm1)
        self.ftop = self.fm1[-1] + 1
        g = spow_inv(Lseries(kind), p)
        self.g = [float(c) for c in g]
        self.coef = self.wt * self.fm1
    def D(self, w):
        w = np.asarray(w, dtype=complex)
        E = np.exp(np.outer(w, self.v))
        s = E @ self.coef
        tail = sum(self.g[j] * np.exp((w + j) * V0) / (w + j) for j in range(1, J))
        return s + tail + np.exp(w * V1) / w

def wind(D, box, n0=2000):
    x0, x1, y0, y1 = box
    corners = [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0)]
    tot, mn, mxstep, npts = 0.0, np.inf, 0.0, 0
    for a, b in zip(corners[:-1], corners[1:]):
        t = np.linspace(0, 1, n0 + 1)
        for _ in range(25):
            z = a + (b - a) * t
            d = D(z)
            st = np.angle(d[1:] / d[:-1])
            bad = np.abs(st) > 0.3
            if not bad.any(): break
            mids = (t[:-1][bad] + t[1:][bad]) / 2
            t = np.sort(np.concatenate([t, mids]))
        tot += st.sum(); mn = min(mn, np.abs(d).min()); mxstep = max(mxstep, np.abs(st).max()); npts += len(t)
    return tot / (2 * np.pi), mn, mxstep, npts

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default='tlz_bessel.json'); ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    mp.mp.dps = 30 if a.quick else 40
    cases = {'uniform_pair': ('uniform', 2), 'bessel_pair': ('bessel', 2), 'bessel_single': ('bessel', 1)}
    boxes = {'low (1/2, 17/2)': (-0.95, -0.05, 0.5, 8.5), 'high (17/2, 11)': (-0.95, -0.05, 8.5, 11.0)}
    out = {'V0': V0, 'V1': V1}
    gam1 = float(mp.zetazero(1).imag); gam2 = float(mp.zetazero(2).imag)
    D0cf = lambda w: -2 * mp.pi * (2 * w + 1) * mp.zeta(-2 * w) / mp.sin(mp.pi * w)
    for name, (kind, p) in cases.items():
        R1 = Rule(kind, p, 0.5, 24 if not a.quick else 12)
        R2 = Rule(kind, p, 0.25, 32 if not a.quick else 16)
        res = {'f_at_top_rule1': R1.ftop, 'rule_diff_at_test_points': float(np.abs(R1.D([-0.3 + 4j, -0.5 + 3j, -0.7 + 10j]) - R2.D([-0.3 + 4j, -0.5 + 3j, -0.7 + 10j])).max())}
        if name == 'uniform_pair':
            res['control_closed_at_-0.3+4i'] = str(complex(D0cf(mp.mpc(-0.3, 4))))
            res['control_quad_at_-0.3+4i'] = str(complex(R1.D([-0.3 + 4j])[0]))
            res['|D|_at_-1/4+i gamma1/2'] = float(abs(R1.D([complex(-0.25, gam1 / 2)])[0]))
            res['|D|_at_-1/4+i gamma2/2'] = float(abs(R1.D([complex(-0.25, gam2 / 2)])[0]))
        for bn, box in boxes.items():
            r1 = wind(R1.D, box); r2 = wind(R2.D, box)
            res[bn] = {'winding_rule1': r1[0], 'minabs_rule1': r1[1], 'maxstep_rule1': r1[2], 'points_rule1': r1[3],
                       'winding_rule2': r2[0], 'minabs_rule2': r2[1], 'maxstep_rule2': r2[2], 'points_rule2': r2[3]}
        out[name] = res
        print(name, json.dumps(res), flush=True)
    json.dump(out, open(a.out, 'w'), indent=1, default=float)

if __name__ == '__main__':
    main()
