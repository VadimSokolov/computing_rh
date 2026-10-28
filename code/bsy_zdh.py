# bsy_zdh.py: shared evaluator of the scripts bsy_cm.py, bsy_count.py and bsy_levy.py (Chapter 35, Balazard, Saias and Yor meet Thorin, group BSY-A).
# Not run by itself. It evaluates the Hurwitz zeta function zeta(s, a) and its s derivative, vectorised in numpy (double precision), by the Euler and Maclaurin formula
#   zeta(s, a) = sum_{k<N} (k+a)^(-s) + (N+a)^(1-s)/(s-1) + (N+a)^(-s)/2 + sum_{j=1}^{M} B_2j/(2j)! (s)_(2j-1) (N+a)^(-s-2j+1),
# with N = ceil(nf |s|) + 10 and M terms (default nf = 0.3, M = 28: the ratio of consecutive remainder terms is (|s|/(2 pi (N+a)))^2 <= 0.28, so the remainder is below 1e-15
# of the sum; the alternative setting nf = 0.45, M = 36 is used by the scripts as a second run). From it: zeta(s), the Davenport and Heilbronn function
# f(s) = 5^(-s) [zeta(s,1/5) - zeta(s,4/5) + kappa (zeta(s,2/5) - zeta(s,3/5))], kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1) (coefficients 1, kappa, -kappa, -1, 0 mod 5), and
# the function h_theta(s) = zeta(s) (s-1)/(s-2 theta+1) of the chapter. The function arg_ray continues arg h along the horizontal ray from Re s = 2 (where |h/R - 1| < 1, so
# the principal arguments of the factors are the continuous ones, Y(2+i tau) = Arg zeta + Arg(1+i tau) - Arg(3-2 theta+i tau); for f, Y = Arg f) down to a given Re s, by principal
# increments of Arg(h(s_(k+1))/h(s_k)) on a mesh that is bisected wherever an increment exceeds 0.9 in absolute value.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import math
import numpy as np
import mpmath as mp

MMAX = 60
_C = [0.0] + [float(mp.bernoulli(2 * j) / mp.factorial(2 * j)) for j in range(1, MMAX + 1)]
SQ5 = math.sqrt(5.0)
KAPPA = (math.sqrt(10 - 2 * SQ5) - 2) / (SQ5 - 1)
DHC = {1: 1.0, 4: -1.0, 2: KAPPA, 3: -KAPPA}


def _block(s, a, N, M, deriv):
    k = np.arange(N, dtype=float) + a
    lk = np.log(k)
    E = np.exp(-np.multiply.outer(s, lk))
    S = E.sum(1)
    dS = -(E * lk).sum(1) if deriv else None
    X = N + a
    lX = math.log(X)
    Xs = np.exp(-s * lX)
    T = X * Xs / (s - 1) + Xs / 2
    dT = (-lX * X * Xs / (s - 1) - X * Xs / (s - 1) ** 2 - lX * Xs / 2) if deriv else None
    p = s.copy()
    q = 1 / s
    Xp = 1.0 / X
    for j in range(1, M + 1):
        T = T + _C[j] * p * Xp * Xs
        if deriv:
            dT = dT + _C[j] * Xp * Xs * p * (q - lX)
        p = p * (s + 2 * j - 1) * (s + 2 * j)
        q = q + 1 / (s + 2 * j - 1) + 1 / (s + 2 * j)
        Xp = Xp / (X * X)
    return (S + T, dS + dT) if deriv else (S + T, None)


def hurwitz(s, a, nf=0.3, M=28, deriv=False, budget=1500000):
    """zeta(s, a) (and d/ds) for a 1-D complex array s, Re s > 0, s not near 1."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    n = len(s)
    order = np.argsort(np.abs(s))
    val = np.empty(n, dtype=complex)
    der = np.empty(n, dtype=complex) if deriv else None
    i = 0
    while i < n:
        bs = 256
        while True:
            j = min(n, i + bs)
            N = int(math.ceil(nf * abs(s[order[j - 1]]))) + 10
            if bs * N <= budget or bs == 1:
                break
            bs = max(1, budget // N)
        idx = order[i:j]
        v, d = _block(s[idx], a, N, M, deriv)
        val[idx] = v
        if deriv:
            der[idx] = d
        i = j
    return (val, der) if deriv else val


def zeta(s, deriv=False, **kw):
    return hurwitz(s, 1.0, deriv=deriv, **kw)


def dh(s, deriv=False, **kw):
    """Davenport and Heilbronn f(s) (entire) and f'(s)."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    f = np.zeros(len(s), dtype=complex)
    g = np.zeros(len(s), dtype=complex) if deriv else None
    for r, c in DHC.items():
        if deriv:
            v, d = hurwitz(s, r / 5.0, deriv=True, **kw)
            f += c * v
            g += c * d
        else:
            f += c * hurwitz(s, r / 5.0, **kw)
    ph = np.exp(-s * math.log(5.0))
    if deriv:
        return ph * f, ph * (g - math.log(5.0) * f)
    return ph * f


def hfun(kind, theta, s, **kw):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    if kind == "zeta":
        return zeta(s, **kw) * (s - 1) / (s - 2 * theta + 1)
    return dh(s, **kw)


def arg_ray(kind, theta, tau, s_end, L=24, tol=0.9, **kw):
    """Y = arg h continued along the horizontal ray from Re s = 2 to Re s = s_end at heights tau; returns (Y, h at s_end)."""
    tau = np.asarray(tau, dtype=float)
    sig = np.linspace(2.0, s_end, L + 1)
    H = np.empty((len(tau), L + 1), dtype=complex)
    for k in range(L + 1):
        H[:, k] = hfun(kind, theta, sig[k] + 1j * tau, **kw)
    if kind == "zeta":
        z0 = zeta(2.0 + 1j * tau, **kw)
        Y = np.angle(z0) + np.angle(1 + 1j * tau) - np.angle(3 - 2 * theta + 1j * tau)
    else:
        Y = np.angle(H[:, 0]).copy()
    D = np.angle(H[:, 1:] / H[:, :-1])
    good = np.abs(D) <= tol
    Y = Y + np.where(good, D, 0.0).sum(1)
    ii, kk = np.nonzero(~good)
    qi, qa, qb = ii, sig[kk], sig[kk + 1]
    qha, qhb = H[ii, kk], H[ii, kk + 1]
    for _ in range(45):
        if len(qi) == 0:
            break
        qm = (qa + qb) / 2
        hm = hfun(kind, theta, qm + 1j * tau[qi], **kw)
        d1 = np.angle(hm / qha)
        d2 = np.angle(qhb / hm)
        ok1 = np.abs(d1) <= tol
        ok2 = np.abs(d2) <= tol
        np.add.at(Y, qi[ok1], d1[ok1])
        np.add.at(Y, qi[ok2], d2[ok2])
        n1, n2 = ~ok1, ~ok2
        qi, qa, qb, qha, qhb = (np.concatenate([qi[n1], qi[n2]]), np.concatenate([qa[n1], qm[n2]]),
                                np.concatenate([qm[n1], qb[n2]]), np.concatenate([qha[n1], hm[n2]]),
                                np.concatenate([hm[n1], qhb[n2]]))
    if len(qi):  # accept what remains (segment shorter than 1e-12): add the principal increment
        np.add.at(Y, qi, np.angle(qhb / qha))
    return Y, H[:, -1]


def mp_h(kind, theta, s, dps=30):
    mp.mp.dps = dps
    s = mp.mpmathify(s)
    if kind == "zeta":
        return complex(mp.zeta(s) * (s - 1) / (s - 2 * theta + 1))
    k = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    h = lambda a: mp.zeta(s, mp.mpf(a) / 5)
    return complex(mp.mpf(5) ** (-s) * (h(1) - h(4) + k * (h(2) - h(3))))
