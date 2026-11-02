#!/usr/bin/env python3
"""Track C, script 8. Zero completed laws for the Davenport and Heilbronn function: the same construction as f_K of
the book, with the zeros of DH on the line below the first off-line zero. If the density is positive the law is
ferromagnetic by construction, and the DH law (which has off-line zeros) lies within its total variation distance.
Run: timeout 500 python3 -I C_dh2.py   (computed; positivity checked on a grid only, not certified)
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 20
KAPPA = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [mp.mpf(1), KAPPA, -KAPPA, mp.mpf(-1)]
Af = [float(a) for a in A] + [0.0]


def f_dh(s):
    return mp.mpf(5) ** (-s) * sum(a * mp.zeta(s, mp.mpf(r + 1) / 5) for r, a in enumerate(A))


def Lre(th):
    s = mp.mpc(0.5, th)
    return mp.re((5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * f_dh(s))


def kernel(t):
    t = np.abs(np.asarray(t, dtype=float))
    y = np.exp(2 * t)
    s = np.zeros_like(y)
    for n in range(1, 40):
        a = Af[(n - 1) % 5]
        if a != 0:
            s = s + a * n * np.exp(-np.pi * n**2 * y / 5)
    return np.exp(1.5 * t) * s


def main():
    ts = np.arange(1.0, 85.05, 0.1)
    vals = [Lre(float(t)) for t in ts]
    zs = []
    for i in range(len(ts) - 1):
        if vals[i] * vals[i + 1] < 0:
            zs.append(float(mp.findroot(Lre, (float(ts[i]), float(ts[i + 1])), solver="anderson", tol=1e-30, verify=False)))
    zs = np.array(zs)
    print(f"{len(zs)} sign changes on the line below 85; last {zs[-1]:.4f}", flush=True)
    xs = np.arange(0, 4.0001, 0.001)
    xq = np.arange(0, 5.0001, 0.0005)
    kq = kernel(xq)
    nq = (2 * np.sum(kq) - kq[0]) * 0.0005
    var = 2 * np.sum(xq**2 * kq / nq) * 0.0005
    p = kernel(xs) / nq
    print(f"Var={var:.10f}, half variance {var/2:.8f}, sum over found zeros of g^-2 = {np.sum(zs**-2.0):.8f}, sum g^-4 = {np.sum(zs**-4.0):.6e}", flush=True)
    # fourth derivative norm through the transform of p itself, on the range where double precision suffices
    th = np.arange(0, 400.0001, 0.01)
    w = np.ones_like(th) * 0.01
    w[0] = w[-1] = 0.005
    for K in (5, 10, 20, 30, 40, len(zs)):
        g = zs[:K]
        beta = var / 2 - np.sum(g**-2.0)
        logabs = -beta * th**2
        sign = np.ones_like(th)
        for z in g:
            ff = 1 - th**2 / z**2
            logabs += np.log(np.maximum(np.abs(ff), 1e-300))
            sign *= np.sign(ff)
        phik = sign * np.exp(logabs)
        fk = np.array([np.sum(w * phik * np.cos(x * th)) / np.pi for x in xs])
        d = np.abs(fk - p)
        tv = (np.sum(d) - 0.5 * d[0] - 0.5 * d[-1]) * 0.001
        print(f"K={K:3d} g_K={g[-1]:8.3f} beta_K={beta:.6f}: TV(f_K, DH)={tv:.4e}; min f_K on [0,4]={fk.min():.3e}; last x with f_K>1e-12: {xs[fk>1e-12][-1]:.3f}", flush=True)


if __name__ == "__main__":
    main()
