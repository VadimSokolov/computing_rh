#!/usr/bin/env python3
"""Track C, script 7. The Hankel tests from cumulants alone, at two resolutions, with magnitudes.
Moments by the trapezoid rule on the half line: the integrand is analytic in |Im t| < pi/4 and decays
double exponentially, so the error is about exp(-2 pi y N) at step 1/N, y < pi/4. No zero of zeta enters
the three real sequences (u) uniform test, (s) sign test, (p) the sequence of inverse powers.
Synthetic part: a spectrum of 60 fixed real numbers (the first 60 ordinates to 15 digits, taken as exact)
with one pair moved off the line: at which n does the Hankel test first go negative, and how small is H_n?
Run: timeout 570 python3 -I C_hankel_check.py   (computed in mpmath, not certified)
"""
import time
import mpmath as mp


def phi_values(N, tmax):
    h = mp.mpf(1) / N
    pi = mp.pi
    tiny = mp.mpf(10) ** (-mp.mp.dps - 20)
    vals, ts = [], []
    for i in range(int(tmax * N) + 1):
        t = i * h
        q = mp.exp(-pi * mp.exp(2 * t))
        a, b = mp.exp(mp.mpf(9) / 2 * t), mp.exp(mp.mpf(5) / 2 * t)
        q2 = q * q
        f, cur = q, q
        s, n = mp.mpf(0), 1
        while cur > tiny:
            s += (2 * pi * pi * n**4 * a - 3 * pi * n**2 * b) * cur
            n += 1
            f = f * q2
            cur = cur * f
        vals.append(2 * s)
        ts.append(t)
    return h, ts, vals


def moments(N, jmax, tmax=3.2):
    h, ts, vals = phi_values(N, tmax)
    x2 = [t * t for t in ts]
    pw = list(vals)
    mom = []
    for j in range(jmax + 1):
        m = 2 * h * mp.fsum(pw[1:]) + (h * vals[0] if j == 0 else 0)
        mom.append(m)
        pw = [u * y for u, y in zip(pw, x2)]
    return mom


def cumulants(mom):
    n = 2 * (len(mom) - 1)
    m = [mp.mpf(0)] * (n + 1)
    for j, v in enumerate(mom):
        m[2 * j] = v / mom[0]
    k = [mp.mpf(0)] * (n + 1)
    for r in range(2, n + 1, 2):
        k[r] = m[r] - mp.fsum(mp.binomial(r - 1, i - 1) * k[i] * m[r - i] for i in range(2, r, 2))
    return k


def leading_minors(seq, nmax):
    """H_1,...,H_nmax of the Hankel matrix of seq, by elimination without pivoting (real or complex input)."""
    A = [[seq[i + j] for j in range(nmax)] for i in range(nmax)]
    out = []
    d = mp.mpf(1)
    for c in range(nmax):
        piv = A[c][c]
        d = d * piv
        out.append(d)
        for i in range(c + 1, nmax):
            f = A[i][c] / piv
            if f != 0:
                row, top = A[i], A[c]
                for j in range(c, nmax):
                    row[j] = row[j] - f * top[j]
    return out


def norm(seq):
    x = seq[-1] / seq[-2]
    return [seq[j] / x ** (j + 2) for j in range(len(seq))]


def sign_string(H):
    return "".join("+" if mp.re(d) > 0 else "-" for d in H)


def lg(d):
    return float(mp.log10(abs(d)))


def sequences(N, jmax):
    mom = moments(N, jmax)
    k = cumulants(mom)
    p = {m: (-1) ** (m + 1) * m * k[2 * m] / mp.factorial(2 * m) for m in range(1, jmax + 1)}
    z = {m: mp.zeta(2 * m) for m in range(2, jmax + 1)}
    q = [mp.pi ** (2 * m) * p[m] / z[m] for m in range(2, jmax + 1)]
    r = [(mp.pi / 2) ** (2 * m) * p[m] / ((1 - mp.mpf(2) ** (-2 * m)) * z[m]) for m in range(2, jmax + 1)]
    pl = [p[m] for m in range(2, jmax + 1)]
    return mom, k, q, r, pl


def main():
    jmax = 64
    nmax = 30
    t0 = time.time()
    mp.mp.dps = 1500
    mom, k, q, r, pl = sequences(700, jmax)
    print(f"resolution 1 done in {time.time()-t0:.0f}s; xi(1/2)={mp.nstr(mom[0],18)} k2={mp.nstr(k[2],15)} k4={mp.nstr(k[4],15)} top atom ratio pi^2/(q_last/q_prev)^(1)...", flush=True)
    Hu = leading_minors(norm(q), 16)
    Hs = leading_minors(norm(r), nmax)
    Hp = leading_minors(norm(pl), nmax)
    mp.mp.dps = 1900
    mom2, k2, q2, r2, pl2 = sequences(900, jmax)
    print(f"resolution 2 done at {time.time()-t0:.0f}s", flush=True)
    Hu2 = leading_minors(norm(q2), 16)
    Hs2 = leading_minors(norm(r2), nmax)
    Hp2 = leading_minors(norm(pl2), nmax)
    rel = lambda a, b: max(float(mp.log10(abs(x - y) / abs(y) + mp.mpf(10) ** (-2000))) for x, y in zip(a, b))
    print("log10 max relative difference between resolutions: uniform", round(rel(Hu, Hu2)), "sign", round(rel(Hs, Hs2)), "inverse powers", round(rel(Hp, Hp2)))
    print("uniform test signs n=1..16:", sign_string(Hu2))
    print("  H_n:", " ".join(f"{n+1}:{mp.nstr(Hu2[n],4)}" for n in range(9)))
    print("sign test signs n=1..30   :", sign_string(Hs2))
    print("  log10|H_n|:", " ".join(f"{n+1}:{lg(Hs2[n]):.1f}" for n in range(0, nmax)))
    print("inverse powers (p_m), must be + for all n exactly under RH, signs n=1..30:", sign_string(Hp2))
    print("  log10|H_n|:", " ".join(f"{n+1}:{lg(Hp2[n]):.1f}" for n in range(0, nmax)))
    # control: a true positive measure (harmonic coefficients, 2000 atoms)
    mp.mp.dps = 1500
    atoms = [(mp.mpf(1) / (2 * (j + mp.mpf("1.25839697")))) ** 2 for j in range(1, 2001)]
    pwr = [a * a for a in atoms]
    ctrl = []
    for m in range(2, 40):
        ctrl.append(mp.fsum(pwr))
        pwr = [u * a for u, a in zip(pwr, atoms)]
    print("control, harmonic uniform coefficients, n=1..18:", sign_string(leading_minors(norm(ctrl), 18)), flush=True)
    # synthetic spectra
    mp.mp.dps = 20
    zs = [mp.zetazero(j).imag for j in range(1, 61)]
    mp.mp.dps = 1200
    zs = [mp.mpf(mp.nstr(z, 15)) for z in zs]
    base = [z ** (-2) for z in zs]
    print("synthetic: first 60 ordinates as exact reals; pair (K,K+1) moved to g +- i beta, g their mean")
    for K, beta in ((3, "0.1"), (6, "0.1"), (10, "0.1"), (15, "0.1"), (3, "0.001"), (10, "0.001"), (20, "0.1")):
        g = (zs[K - 1] + zs[K]) / 2
        c = (g + 1j * mp.mpf(beta)) ** (-2)
        xs = [b for i, b in enumerate(base) if i not in (K - 1, K)] + [c, mp.conj(c)]
        seq = []
        for m in range(2, 62):
            seq.append(mp.re(mp.fsum(x**m for x in xs)))
        ref = [mp.fsum(b**m for b in base) for m in range(2, 62)]
        H = leading_minors(norm(seq), 28)
        Hr = leading_minors(norm(ref), 28)
        first = next((i + 1 for i, d in enumerate(H) if d < 0), None)
        rfirst = next((i + 1 for i, d in enumerate(Hr) if d < 0), None)
        at = lg(H[first - 1]) if first else float("nan")
        print(f"  K={K:2d} beta={beta}: first negative n={first}, log10|H_n| there={at:.1f}; unperturbed spectrum first negative: {rfirst}", flush=True)
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
