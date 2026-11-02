"""Shared helpers for the N4 review scripts. Reads only repository data files, nothing from misc/reviews/."""
import os, json
import numpy as np
import mpmath as mp

ROOT = '.'
KAPPA = '0.2840790438404122960282918'

def zero_strings(n=8000):
    out = []
    with open(os.path.join(ROOT, 'data', 'zeros_8000.txt')) as f:
        for l in f:
            l = l.strip()
            if l:
                out.append(l)
            if len(out) >= n:
                break
    return out

def zeros_float(n=8000):
    return np.array([float(s) for s in zero_strings(n)])

def dh_json():
    return json.load(open(os.path.join(ROOT, 'data', 'wl_dh_zeros.json')))

def xi(s):
    s = mp.mpmathify(s)
    if s == 1:
        return mp.mpf(1) / 2
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

def Xi(t):
    return mp.re(xi(mp.mpf(1) / 2 + 1j * mp.mpmathify(t)))

def Phi_raw(u, N=16):
    """Riemann kernel in the Rodgers and Tao normalisation; the n sum converges for real u and small complex u."""
    e4 = mp.exp(4 * u)
    e9 = mp.exp(9 * u)
    e5 = mp.exp(5 * u)
    return mp.fsum((2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n ** 2 * e4) for n in range(1, N + 1))

def Phi_np(u, N=8):
    u = np.abs(np.asarray(u, dtype=float))
    out = np.zeros_like(u)
    e4 = np.exp(4 * u)
    for n in range(1, N + 1):
        ex = np.exp(np.maximum(-np.pi * n ** 2 * e4 + 9 * u, -700.0))
        ex5 = np.exp(np.maximum(-np.pi * n ** 2 * e4 + 5 * u, -700.0))
        out += 2 * np.pi ** 2 * n ** 4 * ex - 3 * np.pi * n ** 2 * ex5
    return out

def dh_f(s):
    k = mp.mpf(KAPPA)
    a = [mp.mpf(1), k, -k, mp.mpf(-1)]
    s = mp.mpmathify(s)
    if mp.re(s) >= 30:
        cyc = [1, k, -k, -1, 0]
        return mp.fsum(cyc[(n - 1) % 5] * mp.mpf(n) ** (-s) for n in range(1, 400))
    return mp.mpf(5) ** (-s) * mp.fsum(a[r] * mp.zeta(s, mp.mpf(r + 1) / 5) for r in range(4))

def dh_Lam(s):
    s = mp.mpmathify(s)
    return (5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * dh_f(s)
