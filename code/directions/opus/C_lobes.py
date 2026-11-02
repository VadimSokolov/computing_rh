#!/usr/bin/env python3
"""Track C, script 2. What a total variation bound can see.

(A) lobes of phi_X between consecutive zeros, against 2*eps
(B) the height T(eps) beyond which the bound is blind
(C) the zero free box implied by TV < eps, with the constant of Lemma Q in C.md
(D) blindness: explicit positive densities p(x)(1 + c cos Tx) with nonreal zeros near height T
(E) the floor for sums of uniforms and of fair signs through periodisation
(F) check of the Putnam identity zeta(1/2+ikd) = -(1/2+ikd) * (Fourier coefficient of G-hat)

Run: timeout 550 python3 -I C_lobes.py
Status: computed in mpmath and double precision, not certified ball arithmetic.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 40
XI0 = None


def xi(s):
    s = mp.mpmathify(s)
    if abs(s - 1) < mp.mpf(10) ** (-30) or abs(s) < mp.mpf(10) ** (-30):
        return mp.mpf(1) / 2  # the removable points of (s-1) zeta(s) and of s Gamma(s/2)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def phi(theta):
    """phi_X(theta) = xi(1/2 + i theta)/xi(1/2), theta complex allowed."""
    return xi(mp.mpf(1) / 2 + 1j * theta) / XI0


def L(z):
    return xi(mp.mpf(1) / 2 + z) / XI0


def p_density(x):
    """Polya density in double precision, theta series, even in x."""
    x = np.abs(np.asarray(x, dtype=float))
    s = np.zeros_like(x)
    for n in range(1, 8):
        y = np.pi * n**2 * np.exp(2 * x)
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * x) - 3 * np.pi * n**2 * np.exp(2.5 * x)) * np.exp(-y)
    return 2 * s / float(XI0)


def lobes(zeros):
    out = []
    for a, b in zip(zeros[:-1], zeros[1:]):
        f = lambda t: abs(mp.re(phi(mp.mpf(t))))
        ts = [a + (b - a) * j / 40 for j in range(1, 40)]
        vals = [f(t) for t in ts]
        i = int(np.argmax([float(v) for v in vals]))
        lo, hi = ts[max(i - 1, 0)], ts[min(i + 1, 38)]
        for _ in range(30):  # golden style refinement
            m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
            if f(m1) < f(m2):
                lo = m1
            else:
                hi = m2
        T = (lo + hi) / 2
        out.append((a, b, T, f(T), mp.sign(mp.re(phi(T)))))
    return out


def part_ab(zeros):
    print("(A) lobes: k, gamma_k, gamma_k+1, argmax, max|phi_X|, sign")
    lb = lobes(zeros)
    for k, (a, b, T, m, s) in enumerate(lb, 1):
        if k <= 14 or k % 8 == 0:
            print(f"  {k:3d} {mp.nstr(a,8):>10s} {mp.nstr(b,8):>10s} {mp.nstr(T,7):>9s} {mp.nstr(m,4):>11s} {int(s):+d}")
    print("check of the amplitude law |phi_X(t)| = 2.2519 t^(7/4) e^(-pi t/4) |zeta(1/2+it)|:")
    for t in (20, 50, 100):
        t = mp.mpf(t) + mp.mpf("0.3")
        lhs = abs(phi(t))
        rhs = mp.mpf("2.2519") * t ** (mp.mpf(7) / 4) * mp.exp(-mp.pi * t / 4) * abs(mp.zeta(0.5 + 1j * t))
        print(f"  t={mp.nstr(t,5)} ratio {mp.nstr(lhs/rhs,8)}")
    print("(B) real line blindness: last lobe with max|phi_X| >= 2 eps, then the smooth law for tiny eps")
    for eps in ("3e-3", "1.48e-2", "1.266e-4", "6.024e-6", "1.932e-6", "1e-8", "1e-10", "1e-15", "1e-20"):
        e = mp.mpf(eps)
        last = [k for k, r in enumerate(lb, 1) if r[3] >= 2 * e]
        if last:
            k = last[-1]
            print(f"  eps={eps:>9s}: last visible lobe k={k:3d} ({mp.nstr(lb[k-1][0],6)},{mp.nstr(lb[k-1][1],6)}), so zeros above gamma_{k+1} are not pinned at all")
        else:
            print(f"  eps={eps:>9s}: no lobe visible")
    print("  smooth law: pi t/4 - (7/4) log t - log 2.2519 = log(1/(2 eps)) with |zeta| = 1:")
    for dec in (6, 10, 20, 50, 100, 1000, 10**6):
        L_ = dec * mp.log(10) - mp.log(2)
        t = mp.findroot(lambda t: mp.pi * t / 4 - 1.75 * mp.log(t) - mp.log(2.2519) - L_, 4 * L_ / mp.pi + 10)
        print(f"  eps=1e-{dec}: T = {mp.nstr(t,7)}")
    H = mp.mpf("3e12")
    dec = (mp.pi * H / 4 - 1.75 * mp.log(H)) / mp.log(10)
    print(f"  to reach the verified height 3e12: eps = 10^(-{mp.nstr(dec,6)}), {mp.nstr(dec*mp.log(10)/mp.log(2),5)} bits")
    return lb


def sigma_bootstrap(eps):
    """Upper bound for the standard deviation of a ferromagnetic law within TV eps of X (Lemma Q, step 1)."""
    xs = np.linspace(0, 3, 300001)
    pv = p_density(xs)
    tail = 2 * (np.cumsum(pv[::-1])[::-1]) * (xs[1] - xs[0])  # P(|X| > x)
    A0 = xs[np.argmax(tail < 3 / 16 - eps)]
    s2 = (2 * A0) ** 2
    hist = [float(s2)]
    L_ = np.log(1 / eps)
    var = 0.0462099862308379
    for _ in range(6):
        A2 = 2 * s2 * L_
        s2 = var + 4 * eps * A2 + 4 * s2 * eps + 2 * eps  # the last term covers the tail of X, far below eps
        hist.append(float(s2))
    return float(A0), np.sqrt(hist)


def part_c():
    print("(C) zero free boxes implied by TV < eps through Lemma Q: Delta = 5 eps exp(R sigma sqrt(2 log(1/eps))), R = 1/2")
    zeros = [mp.zetazero(k).imag for k in range(1, 18)]
    ts = sorted([mp.mpf(j) / 40 for j in range(0, 1801)] + zeros)
    ys = [mp.mpf(10) ** (-mp.mpf(j) / 8) for j in range(32, 0, -1)]
    ys = [y for y in ys if y < 0.5] + [mp.mpf(1) / 2]
    table = {y: np.array([float(abs(L(y + 1j * t))) for t in ts]) for y in ys}
    tsf = np.array([float(t) for t in ts])
    for eps in (1.266e-4, 6.024e-6, 1.932e-6, 1e-8, 1e-10):
        A0, sig = sigma_bootstrap(eps)
        sigma = sig[-1]
        Delta = 5 * eps * np.exp(0.5 * sigma * np.sqrt(2 * np.log(1 / eps)))
        best = None
        for T in np.arange(14.5, 45.01, 0.25):
            mask = tsf <= T
            if table[ys[-1]][mask].min() <= Delta:
                continue
            for yy in ys[:-1]:
                if table[yy][mask].min() > Delta:
                    itop = np.argmin(np.abs(tsf - T))
                    if all(table[y][itop] > Delta for y in ys if y >= yy):
                        best = (float(T), float(yy))
                        break
        print(f"  eps={eps:.3e}: A0={A0:.4f} sigma bound {sigma:.4f}, Delta={Delta:.3e}; largest box height T={best[0] if best else None}, least distance y0 at that T={best[1] if best else None}")
        # the least y0 for T just above each zero
        for g in zeros[:6]:
            T = float(g) + 0.05
            mask = tsf <= T
            found = None
            if table[ys[-1]][mask].min() > Delta:
                for yy in ys[:-1]:
                    itop = np.argmin(np.abs(tsf - T))
                    if table[yy][mask].min() > Delta and all(table[y][itop] > Delta for y in ys if y >= yy):
                        found = float(yy)
                        break
            print(f"     box up to gamma={float(g):8.3f}: y0={found}")


def part_d(lb):
    print("(D) blindness: p_c(x) = p(x)(1 + c cos Tx)/(1 + c phi(T)), c>0 on negative lobes, c<0 on positive lobes")
    print("    k  T  m  c  TV(p_c,p)  a nonreal zero of the transform of p_c   lobe real zero free?  zeros in box")
    xs = np.linspace(-3, 3, 240001)
    dx = xs[1] - xs[0]
    pv = p_density(xs)
    for k in (1, 2, 3, 4, 5, 6, 7, 10, 14):
        a, b, T, m, s = lb[k - 1]
        h = (b - a) / 2
        # level of the shifted copy at the ends of the lobe
        c = -s * mp.mpf("2.6") * m / mp.re(phi(h))
        Tf = float(T)
        phiT = mp.re(phi(T))
        norm = 1 + c * phiT
        fc = lambda th: (phi(th) + c / 2 * (phi(th - T) + phi(th + T))) / norm
        tv = 0.5 * float(abs(c) / norm) * np.sum(pv * np.abs(np.cos(Tf * xs) - float(phiT))) * dx
        z = None
        for start in (T + 0.5j * h, T + 0.9j * h, T + 0.25j * h):
            try:
                z = mp.findroot(fc, start, tol=mp.mpf(10) ** (-20) * abs(m), verify=False, maxsteps=60)
            except Exception:
                continue
            if abs(mp.im(z)) > 1e-3 and abs(fc(z)) < mp.mpf(10) ** (-6) * abs(m):
                break
            z = None
        grid = [a - h / 3 + (b - a + 2 * h / 3) * j / 80 for j in range(81)]
        vals = [mp.re(fc(t)) for t in grid]
        sc = sum(1 for u, v in zip(vals[:-1], vals[1:]) if u * v < 0)
        # argument principle on the box [a-h/3, b+h/3] x [-1.2h, 1.2h]
        x0, x1, Y = a - h / 3, b + h / 3, 1.2 * h
        pts = []
        n = 60
        pts += [x0 + (x1 - x0) * j / n - 1j * Y for j in range(n)]
        pts += [x1 + 1j * (-Y + 2 * Y * j / n) for j in range(n)]
        pts += [x1 - (x1 - x0) * j / n + 1j * Y for j in range(n)]
        pts += [x0 + 1j * (Y - 2 * Y * j / n) for j in range(n)]
        fv = [fc(q) for q in pts] + [fc(pts[0])]
        w = sum(mp.arg(fv[i + 1] / fv[i]) for i in range(len(pts))) / (2 * mp.pi)
        zs = mp.nstr(z, 8) if z is not None else "not found"
        print(f"  {k:3d} {mp.nstr(T,6):>8s} {mp.nstr(m,3):>9s} {mp.nstr(c,3):>10s} {tv:9.3e}  {zs:>34s}  real sign changes near lobe: {sc}  winding: {mp.nstr(w,3)}")


def periodisation_tv(theta0, kind):
    """TV(X mod 2c, uniform) for uniform (theta0 = pi/c), or half of TV(X mod 4a, shift 2a) for signs (theta0 = pi/(2a))."""
    ns = list(range(1, 14)) if kind == "u" else list(range(1, 14, 2))
    coeff = [float(mp.re(phi(n * theta0))) for n in ns]
    v = np.linspace(0, 2 * np.pi, 20001)
    s = sum(cf * np.cos(n * v) for cf, n in zip(coeff, ns))
    return float(np.sum(np.abs(s[:-1])) * (v[1] - v[0]) / (2 * np.pi))


def part_e():
    print("(E) floor for sums of uniforms and of fair signs (coefficients close to pi/gamma_1 are forced)")
    g1 = mp.zetazero(1).imag
    for kind, mult in (("u", 2), ("s", 3)):
        A = abs(phi(mult * g1))
        print(f"  kind={kind}: |phi_X({mult} gamma_1)| = {mp.nstr(A,6)}, (2/pi)|.| = {mp.nstr(2*A/mp.pi,6)}, (1/2)|.| = {mp.nstr(A/2,6)}")
        # scan theta_S near gamma_1; eps feasible only if M(theta_S) <= eps and theta_S in (theta_1(eps), theta_2(eps)]
        d1 = abs(mp.diff(lambda t: mp.re(phi(t)), g1))
        best = None
        for eps in [float(A) * f for f in np.linspace(0.70, 0.40, 61)]:
            w = 2 * eps / float(d1) * 1.05  # width of the window where |phi_X| <= 2 eps, slope estimate
            scan = np.linspace(float(g1) - w, float(g1) + w, 41)
            Mmin = min(periodisation_tv(mp.mpf(float(t)), kind) for t in scan)
            if Mmin > eps:
                best = (eps, Mmin)
                break
        print(f"  largest eps on the grid with min M > eps: {best}")
    print("  allowed scales for sums of uniforms: M(theta_0) at theta_0 = pi/c")
    for t0 in (10, 12, 14, 14.13472514, 15, 18, 21.02203964, 23, 25.01085758, 27, 30.42487613, 32.93506159):
        print(f"    theta_0={t0:>12}: M={periodisation_tv(mp.mpf(t0),'u'):.3e}")


def part_f():
    print("(F) Putnam identity: zeta(1/2+ikd)/(1/2+ikd) = -Fourier coefficient of G-hat(e^w) on [0,log a]")
    a = 2.0
    d = 2 * np.pi / np.log(a)
    C = 1 / (1 - a ** -0.5)
    N = 400000
    w = (np.arange(N) + 0.5) * np.log(a) / N
    u = np.exp(w)
    E = np.zeros_like(u)
    j = 0
    while a**j <= u.max() + 1e-12:
        E += a ** (j / 2) * np.floor(u / a**j)
        j += 1
    F = C * u - E
    G = F / np.sqrt(u)
    # tail series of {a^i u}/(a^i u)^(1/2); for a=2, keep i up to 40 with care for double precision
    R = np.zeros_like(u)
    for i in range(1, 41):
        au = a**i * u
        R += (au - np.floor(au)) / np.sqrt(au)
    Ghat = G + R
    for k in (0, 1, 2, 3):
        coef = np.sum(Ghat * np.exp(-1j * k * d * w)) * (np.log(a) / N)
        s = mp.mpc(0.5, k * d)
        print(f"  k={k}: Fourier coefficient {complex(coef):.6f}, -zeta(s)/s = {complex(-mp.zeta(s)/s):.6f}")


def main():
    global XI0
    XI0 = xi(mp.mpf(1) / 2)
    zeros = [mp.zetazero(k).imag for k in range(1, 42)]
    lb = part_ab(zeros)
    part_c()
    part_d(lb)
    part_e()
    part_f()


if __name__ == "__main__":
    main()
