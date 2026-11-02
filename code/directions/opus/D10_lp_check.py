#!/usr/bin/env python3
"""D10: floating point check of the inequality behind the conjectured lower bound for one pair among real zeros.

Normalise eta = 1.  Suppose the rising curve exists at height w < kappa_1 (point p=(xp,w)) and at the height
r_c = 0.94169... (point q=(xq,r_c)).  Then T(p)+A(p) = T(q)+A(q) = L = (m+1/2) pi, with m >= 0 after a mirror,
and the slope bound s(p) <= A_x(p) holds, where T = phi y + sum over zeros k_a + Gaussian term.
A real zero is allowed only if its own slope s_a(p) = w/((xp-a)^2+w^2) is at most A_x(p).
Eliminating phi gives the linear programme: sum E_a = D, sum s_a <= A_x(p), so the configuration is
impossible as soon as  D > mu* A_x(p), mu* = max over allowed a of E_a/s_a (and x_q-x_p for a at infinity or the
Gaussian term).  We test that on a grid and report the worst margin.  No proof: floating point on a grid.

Usage: timeout 300 python3 -I D10_lp_check.py [nw nxp nxq]
"""
import sys
import numpy as np

t = 2.798386045783887  # tan t = -1/t
omega = np.pi - t
rc = np.cos(omega)
beta_c = np.sqrt(1 + t * t)


def kappa1():
    lo, hi = 0.1, 0.9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.arccos(mid) / mid > beta_c:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


K1 = kappa1()


def A(x, y):
    return np.angle((x + 1j * y) ** 2 + 1)


def Ax(x, y):
    z = x + 1j * y
    return (2 * z / (z * z + 1)).imag


def run(nw, nxp, nxq, qheights=(rc,)):
    worst = (np.inf, None)
    dist = np.geomspace(1e-3, 1e5, 500)
    for w in np.linspace(0.02, K1 * 0.9999, nw):
        ap = np.sqrt(1 - w * w)
        for xp in np.linspace(-ap, ap, nxp + 2)[1:-1]:
            Au = Ax(xp, w)
            if Au <= 0:
                continue
            for hq in qheights:
                aq = np.sqrt(1 - hq * hq)
                for xq in np.linspace(-aq, aq, nxq):
                    D = (np.pi / 2 - A(xp, w)) - (w / hq) * (np.pi / 2 - A(xq, hq))
                    # zeros on the left and on the right of xp, at the allowed distances
                    s = w / (dist ** 2 + w ** 2)
                    ok = s <= Au
                    d = dist[ok]
                    if len(d) == 0:
                        mu = xq - xp
                    else:
                        sa = w / (d ** 2 + w ** 2)
                        # left zero a = xp - d
                        kpL = np.angle(d + 1j * w)
                        kqL = np.angle(xq - xp + d + 1j * hq)
                        EL = kpL - (w / hq) * kqL
                        # right zero a = xp + d: branch arg - pi
                        kpR = np.angle(-d + 1j * w) - np.pi
                        kqR = np.angle(xq - xp - d + 1j * hq) - np.pi
                        ER = kpR - (w / hq) * kqR
                        mu = max(np.max(EL / sa), np.max(ER / sa), xq - xp)
                    margin = D - max(mu, 0.0) * Au
                    if margin < worst[0]:
                        worst = (margin, (w, xp, xq, hq, D, mu, Au))
    return worst


def main():
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nxp = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    nxq = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    print("t", t, "omega", omega, "r_c", rc, "beta_c", beta_c, "kappa_1", K1)
    print("worst margin with q at r_c:", run(nw, nxp, nxq))


if __name__ == "__main__":
    main()
