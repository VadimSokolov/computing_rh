#!/usr/bin/env python3
# Track A numerics. Run: timeout 600 python3 -I A_numerics.py SECTION [SECTION ...]
# Sections: surgery, threshold, spectrum, riemann, dh. Floating point (not certified) except where stated.
# Reads only the zero tables data/zeros_8000.txt and data/wl_dh_zeros.json of the repository.
import json
import sys

import numpy as np

ROOT = './data/'
Z = np.array([float(x) for x in open(ROOT + 'zeros_8000.txt').read().split()])


def lin(x):
    return np.asarray(x, dtype=float)


def sigma_online(a, t, zeros, reach=15.0):
    """Green potential Sigma(1/2+a+it) of on line zeros in the strip of width one."""
    z = zeros[np.abs(zeros - t) < reach]
    w = np.pi * (a + 1j * (t - z)) / 2
    return np.sum(np.log(np.abs(np.cos(w) / np.sin(w))))


def arg_sum(x, y, tj):
    return np.sum(np.arctan2(y, x - tj))


def curve_points(y, tj):
    """tau_k(y), k=1..N-1, where arg P(tau+iy)=k pi; arg_sum is decreasing in x."""
    N = len(tj)
    out = []
    for k in range(1, N):
        lo, hi = tj[0] - 1e5, tj[-1] + 1e5
        for _ in range(90):
            mid = 0.5 * (lo + hi)
            if arg_sum(mid, y, tj) > k * np.pi:
                lo = mid
            else:
                hi = mid
        out.append(0.5 * (lo + hi))
    return np.array(out)


def logP(t, tj):
    return np.sum(np.log(t[..., None] - tj), axis=-1)


class Surgery:
    def __init__(self, tj, y, gap_rule='min', line=0.5):
        self.tj, self.y, self.line = tj, y, line
        tau = curve_points(y, tj)
        lp = np.array([np.sum(np.log(np.abs(x + 1j * y - tj))) for x in tau])
        i = int(np.argmin(lp)) if gap_rule == 'min' else len(tau) // 2
        self.k = i + 1
        self.tau = tau[i]
        self.logc = lp[i]
        self.sign = (-1.0) ** (self.k + 1)  # c = -P(tau+iy), P has argument k pi
        self.lp_all = lp

    def R(self, t):
        return self.sign * np.exp(self.logc - logP(t, self.tj))

    def w(self, t):
        r1, r2 = self.R(t), self.R(-t)
        return (1 + r1) * (1 + r2) - 1

    def E(self, t):
        r = self.R(t)
        return -np.sum(1 / (t[..., None] - self.tj), axis=-1) * r / (1 + r)

    def data(self, pad=30.0, step=0.01):
        v = np.arange(self.tj[0] - pad, self.tj[-1] + pad, step)
        t = v - 1j * self.line
        w = self.w(t)
        D = -1j * (self.E(t) - self.E(-t))
        i = int(np.argmax(np.abs(w)))
        # hiding ratio r = |c| / inf |P(v-i line)|
        r = np.exp(self.logc - np.min(np.real(logP(t, self.tj))))
        return dict(eps=np.abs(w).max(), v=v[i], du=np.abs(D.real).max() / np.pi, dfull=np.abs(D).max(), r=r)


def sec_surgery():
    print('# surgery: eps_F = sup |F/xi-1| on Re s=1, Sigma at the created zero, Theorem 1 says eps_F*exp(Sigma)>=1')
    print('center  N      y     tau        r          eps_F      sup|du|    exp(-Sigma)  eps_F*exp(Sigma)  gridmax')
    for name, c0 in [('H100', 29), ('H1000', 649), ('H8000', 7900), ('Lehmer', None)]:
        if c0 is None:
            c0 = int(np.argmin(np.abs(Z - 7005.08))) + 1
        for N in [2, 4, 8, 16, 32, 64]:
            lo = max(0, c0 - N // 2)
            tj = Z[lo:lo + N]
            for y in [0.05, 0.1, 0.2, 0.3, 0.4]:
                S = Surgery(tj, y)
                d = S.data()
                sig = sigma_online(y, S.tau, Z)
                # full inequality of Theorem 1 on a grid of the strip
                A = np.array([-0.49, -0.4, -0.25, -0.1, -0.02, 0.02, 0.1, 0.25, 0.4, 0.49, y])
                T = np.arange(tj[0] - 3, tj[-1] + 3, 0.05)
                gm = 0.0
                for a in A:
                    wv = np.abs(S.w(T - 1j * a))
                    sv = np.array([sigma_online(abs(a), x, Z) for x in T])
                    gm = max(gm, np.max(wv * np.exp(-sv)) / d['eps'])
                print(f"{name:7s}{N:3d}  {y:5.2f} {S.tau:10.3f} {d['r']:.3e}  {d['eps']:.3e}  {d['du']:.3e}  {np.exp(-sig):.3e}   {d['eps']*np.exp(sig):9.3f}  {gm:.4f}")


def sec_threshold():
    print('# threshold: exp(-Sigma(1/2+y+it)) over t in a window of length 20, against (H/2pi)^(-(1/2-y)/2)')
    print('H       y     predicted   geom.mean   min         max         y_min(eps=1e-2)  y_min(1e-3)')
    for H in [100.0, 1000.0, 8000.0]:
        ts = np.arange(H - 10, H + 10, 0.02)
        for y in [0.02, 0.05, 0.1, 0.2, 0.3, 0.4]:
            v = np.array([sigma_online(y, t, Z) for t in ts])
            pred = (H / 2 / np.pi) ** (-(0.5 - y) / 2)
            e = np.exp(-v)
            print(f"{H:7.0f} {y:5.2f} {pred:.3e}   {np.exp(np.mean(np.log(e))):.3e}   {e.min():.3e}   {e.max():.3e}   {0.5-2*np.log(1e2)/np.log(H/2/np.pi):+.3f}  {0.5-2*np.log(1e3)/np.log(H/2/np.pi):+.3f}")
    print('# fraction of the level y line (window 20) excluded, i.e. Sigma<log(1/eps), at H=1000')
    ts = np.arange(990, 1010, 0.02)
    for eps in [1e-1, 1e-2, 1e-3]:
        row = []
        for y in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]:
            v = np.array([sigma_online(y, t, Z) for t in ts])
            row.append(f"{np.mean(v < np.log(1/eps)):.3f}")
        print(f"eps={eps:.0e} y=0.01,0.05,0.1,0.2,0.3,0.4 :", ' '.join(row))


def sec_spectrum():
    print('# spectrum on Re s=1 of the perturbation c/P: |c sum_j exp(-i t_j w)/P\'(t_j)| exp(-w/2), against w')
    for c0, name in [(649, 'H1000'), (7900, 'H8000')]:
        H = Z[c0]
        for N in [16, 32, 64]:
            tj = Z[c0 - N // 2:c0 + N // 2]
            S = Surgery(tj, 0.1)
            lpp = np.array([np.sum(np.log(np.abs(tj[i] - np.delete(tj, i)))) for i in range(N)])
            sgn = (-1.0) ** (N - 1 - np.arange(N))
            om = np.linspace(0.02, 2 * np.log(H), 2000)
            amp = np.array([abs(np.sum(sgn * np.exp(S.logc - lpp - 1j * tj * w))) for w in om]) * np.exp(-om / 2)
            i = int(np.argmax(amp))
            print(f"{name} N={N}: logH={np.log(H):.3f} peak at w={om[i]:.3f} (0.5*log(H/2pi)={0.5*np.log(H/2/np.pi):.3f}, 0.5*logH={0.5*np.log(H):.3f}), peak amp={amp[i]:.3e}")
            for w0 in [1.0, 2.0, np.log(H) / np.pi, 0.5 * np.log(H / 2 / np.pi), 0.5 * np.log(H), np.log(H)]:
                j = int(np.argmin(np.abs(om - w0)))
                print(f"    w={om[j]:.3f}  n=exp(w)={np.exp(om[j]):9.2f}  amp={amp[j]:.3e}")


def sec_riemann():
    import mpmath as mp
    mp.mp.dps = 80
    print("# Riemann's formula xi(s)=1/2+1/2 s(s-1) sum T_n(s); R_X keeps n<=X; relative error at 1/2+a+i tau, and the bound |s(s-1)| e^{-pi X^2}/(pi X^2)")

    def T(n, s):
        A = mp.pi * n * n
        return A ** (-s / 2) * mp.gammainc(s / 2, A) + A ** (-(1 - s) / 2) * mp.gammainc((1 - s) / 2, A)

    def xi(s):
        return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
    print('tau    a     X   tau/4  X^2   |xi|        |xi-R_X|    bound       rel.err')
    for tau in [20, 40, 80, 160]:
        for a in [0.1, 0.4]:
            s = mp.mpf('0.5') + a + 1j * tau
            x = xi(s)
            tot = mp.mpf(0)
            for X in range(1, 9):
                tot += T(X, s)
                RX = mp.mpf('0.5') + s * (s - 1) / 2 * tot
                bound = abs(s * (s - 1)) * mp.exp(-mp.pi * (X + 1) ** 2) / (mp.pi * (X + 1) ** 2) * 1.01
                if X in (1, 2, 3, 4, 5, 6, 7, 8):
                    print(f"{tau:4d} {a:4.1f} {X:3d} {tau/4:6.1f} {X*X:4d}  {mp.nstr(abs(x),3):10s}  {mp.nstr(abs(x-RX),3):10s}  {mp.nstr(bound,3):10s}  {mp.nstr(abs(x-RX)/abs(x),3)}")


def sec_dh():
    d = json.load(open(ROOT + 'wl_dh_zeros.json'))
    on = lin([float(x) for x in d['online']])
    off = [(float(a), float(b)) for a, b in d['offline']]
    rho = [complex(0.5, g) for g in on]
    for a, b in off:
        rho += [complex(a, b), complex(1 - a, b)]
    rho = np.array(rho)
    print('# Davenport and Heilbronn: reference Lambda, strip -1<sigma<2, data line sigma=2 (a=3/2); Theorem 1 with the Green function of width 3')

    def sig3(s):
        z = rho[np.abs(rho.imag - s.imag) < 30]
        return np.sum(np.log(np.abs(np.sin(np.pi * (s + np.conj(z) + 2) / 6) / np.sin(np.pi * (s - z) / 6))))
    print('level y (s=1/2+y+i tau), at tau=90.2 (no zero near): exp(-Sigma_3) vs (5 tau/2pi)^(-(3/2-y)/2)')
    for y in [0.05, 0.2, 0.31, 0.4, 0.8]:
        s = complex(0.5 + y, 90.2)
        print(f"  y={y:4.2f} exp(-Sigma3)={np.exp(-sig3(s)):.3e} predicted={(5*90.2/2/np.pi)**(-(1.5-y)/2):.3e}")
    print('exp(-Sigma3) on the real part sweep through the off line zero 0.80852+85.69935i, at tau=85.69935:')
    for sg in [0.55, 0.65, 0.75, 0.80, 0.8085, 0.82, 0.9, 1.1]:
        s = complex(sg, 85.69935)
        print(f"  sigma={sg:6.4f} exp(-Sigma3)={np.exp(-sig3(s)):.3e}")
    print('surgery on on line zeros near 90 with c from y, data line a=3/2, Theorem 1 check')
    idx = int(np.argmin(np.abs(on - 90.0)))
    for N in [4, 8, 16]:
        tj = on[max(0, idx - N // 2):idx + N // 2]
        for y in [0.1, 0.3, 0.6]:
            S = Surgery(tj, y, line=1.5)
            dd = S.data(pad=20.0, step=0.02)
            sstar = complex(0.5 + y, S.tau)
            sg = sig3(sstar)
            T = np.arange(tj[0] - 2, tj[-1] + 2, 0.1)
            gm = 0.0
            for a in [-1.45, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 1.45, y]:
                wv = np.abs(S.w(T - 1j * a))
                sv = np.array([sig3(complex(0.5 + a, x)) for x in T])
                gm = max(gm, np.max(wv * np.exp(-sv)) / dd['eps'])
            print(f"  N={N:2d} y={y:4.2f} tau={S.tau:8.3f} r={dd['r']:.3e} eps_F={dd['eps']:.3e} exp(-Sigma3)={np.exp(-sg):.3e} eps_F*exp(Sigma3)={dd['eps']*np.exp(sg):8.3f} gridmax={gm:.4f}")



def sec_recon():
    import mpmath as mp
    mp.mp.dps = 40
    print('# sech reconstruction: xi(s0) = int Re xi(1+iv) / sin(pi(s0-iv)) dv (Proposition 2), trapezoid h=0.05')

    def xi(s):
        if abs(s - 1) < mp.mpf(10) ** (-25):
            return mp.mpf(1) / 2  # removable: xi(1) = 1/2
        return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

    def recon(s0, center, half=16.0, h=0.05, eps=0.0, seed=1, cache={}):
        key = (center, half, h)
        if key not in cache:
            vs = [center - half + j * h for j in range(int(2 * half / h) + 1)]
            cache[key] = (vs, [xi(mp.mpc(1, v)) for v in vs])
        vs, xs = cache[key]
        rng = np.random.default_rng(seed)
        ph = rng.uniform(0, 2 * np.pi, 6)
        fr = rng.uniform(0.5, 4.0, 6)
        tot = mp.mpc(0)
        for v, x in zip(vs, xs):
            g = sum(np.exp(1j * (p + f * v)) for p, f in zip(ph, fr)) / 6.0
            x2 = x * (1 + eps * mp.mpc(g))
            tot += mp.re(x2) / mp.sin(mp.pi * (s0 - 1j * v))
        return tot * h
    s0 = mp.mpf('0.5')
    print('xi(1/2) direct', mp.nstr(xi(s0), 15), ' reconstruction', mp.nstr(mp.re(recon(s0, 0.0)), 15))
    for s0, cen in [(mp.mpc(0.7, 14), 14.0), (mp.mpc(0.5, 14.134725141734693), 14.0), (mp.mpc(0.55, 21.0), 21.0)]:
        r, d = recon(s0, cen, half=30.0), xi(s0)
        print('s0', mp.nstr(s0, 8), ' direct', mp.nstr(d, 8), ' recon', mp.nstr(r, 8), ' |diff|/|xi(1+it)|', mp.nstr(abs(r - d) / abs(xi(mp.mpc(1, s0.imag))), 3))
    print('height 1000: relative error of the reconstruction under relative data noise eps; last column eps*exp(Sigma)')
    for a in [0.05, 0.1, 0.3]:
        s0 = mp.mpc(0.5 + a, 1000.5)
        d = xi(s0)
        sg = sigma_online(a, 1000.5, Z)
        base = recon(s0, 1000.5)
        print(f'  a={a}: noiseless |recon-direct|/|xi|={mp.nstr(abs(base-d)/abs(d),3)}')
        for eps in [1e-6, 1e-3]:
            r = recon(s0, 1000.5, eps=eps)
            print(f'    eps={eps:.0e}  rel.err={mp.nstr(abs(r-d)/abs(d),3)}  eps*exp(Sigma)={eps*np.exp(sg):.3e}')


def sec_hidden():
    print('# zero hidden above the truncation height: target tau, data on v<=tau-D; compare -log eps_T with -log r_c + n D log(3+2 sqrt 2) (uniform model)')
    i0 = 7700
    tau0 = Z[i0]
    n = np.log(tau0 / 2 / np.pi) / 2 / np.pi
    print(f'tau0={tau0:.3f}, n={n:.4f}, n*log(3+2sqrt2)={n*np.log(3+2*np.sqrt(2)):.4f}')
    y = 0.1
    for D in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
        best = None
        for fac in [0.7, 1.0, 1.2, 1.414, 1.7, 2.0, 2.5, 3.0]:
            W = fac * D
            tj = Z[np.abs(Z - tau0) < W / 2]
            if len(tj) % 2:
                tj = tj[:-1]
            if len(tj) < 2:
                continue
            S = Surgery(tj, y, gap_rule='central')
            T = S.tau - D
            v = np.arange(T - 80, T, 0.01)
            m = np.min(np.real(logP(v - 0.5j, tj)))
            allv = np.arange(tj[0] - 40, tj[-1] + 40, 0.01)
            rc = S.logc - np.min(np.real(logP(allv - 0.5j, tj)))
            eT = S.logc - m
            if best is None or eT < best[0]:
                best = (eT, fac, len(tj), rc)
        if best is None:
            print(f'D={D:5.1f} window contains too few zeros')
            continue
        eT, fac, N, rc = best
        print(f'D={D:5.1f} best W/D={fac:5.3f} N={N:3d} log r_c={rc:8.3f} log eps_T={eT:9.3f} (log r_c-log eps_T)/D={(rc-eT)/D:7.3f}  uniform model {n*np.log(3+2*np.sqrt(2)):7.3f}')



def sec_iso():
    print('# isolation threshold: for each gap the saddle of Sigma on the critical line; below eps_iso every competitor keeps those zeros simple and on the line (Corollary 1(d))')
    print('window            gaps  (H/2pi)^(-1/4)  median e^{-saddle}  min e^{-saddle}  at gap (length)')
    for lo, hi in [(14, 100), (900, 1100), (7000, 7010), (7900, 8100)]:
        idx = np.where((Z > lo) & (Z < hi))[0]
        vals = []
        for i in idx[:-1]:
            ts = np.linspace(Z[i], Z[i + 1], 203)[1:-1]
            near = Z[np.abs(Z - 0.5 * (Z[i] + Z[i + 1])) < 15.0]
            tt = np.prod(np.tanh(np.pi * np.abs(ts[:, None] - near) / 2), axis=1)
            vals.append((tt.max(), Z[i + 1] - Z[i], Z[i]))
        vals = np.array(vals)
        j = int(np.argmin(vals[:, 0]))
        H = 0.5 * (lo + hi)
        print(f'[{lo:5d},{hi:5d}] {len(vals):5d}  {(H/2/np.pi)**-0.25:.3e}       {np.median(vals[:,0]):.3e}          {vals[j,0]:.3e}        {vals[j,2]:.3f} ({vals[j,1]:.4f})')


if __name__ == '__main__':
    table = dict(surgery=sec_surgery, threshold=sec_threshold, spectrum=sec_spectrum, riemann=sec_riemann, dh=sec_dh, recon=sec_recon, hidden=sec_hidden, iso=sec_iso)
    for name in sys.argv[1:]:
        table[name]()
