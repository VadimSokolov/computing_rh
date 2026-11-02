#!/usr/bin/env python3
# Track E. (g) the PF3 condition near 0 for L function laws, (h) smooth two sided truncations of the theta series,
# (i) the first zero of the range (Dyck height) approximants, (j) the kink of the folded truncation and the height T_M.
# Run: timeout 500 python3 -I E_more.py [g|h|i|j]
import sys, math
sys.path.insert(0, '/tmp/e_lib')
from flint import arb, acb, arb_series, ctx
import numpy as np
which = sys.argv[1] if len(sys.argv) > 1 else 'ghij'
ctx.prec = 300
ctx.cap = 8
pi = arb.pi()


def kron(D, n):
    # Kronecker symbol (D/n) for fundamental D, n>=1, by multiplicativity over the prime factors of n
    def leg(D, p):
        if p == 2:
            if D % 2 == 0:
                return 0
            return 1 if D % 8 in (1, 7) else -1
        v = pow(D % p, (p - 1) // 2, p)
        return -1 if v == p - 1 else v
    r = 1
    m = n
    p = 2
    while m > 1:
        if p * p > m:
            p = m
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            r *= leg(D, p)**e
            if r == 0:
                return 0
        p += 1
    return r


def V2V4(coef, q, odd):
    """V''(0), V''''(0) for the even density of Lambda(1/2+z)/Lambda(1/2), coef(n) real, conductor q."""
    T = arb_series([arb(0), arb(1)], prec=6)
    w = (2 * T).exp()
    tot = None
    N = int(14 * math.sqrt(q)) + 8
    for n in range(1, N):
        a = coef(n)
        if a == 0:
            continue
        if odd:
            term = a * n * ((-pi * n**2 * w * (arb(1) / q)).exp()) * ((T * (arb(3) / 2)).exp())
        else:
            # even primitive, q>1: theta(x)=sum a_n exp(-pi n^2 x/q), Lambda(s)=(q/pi)^{s/2}Gamma(s/2)L(s), weight x^{s/2-1}dx
            term = a * ((-pi * n**2 * w * (arb(1) / q)).exp()) * ((T * (arb(1) / 2)).exp())
        tot = term if tot is None else tot + term
    L = tot.log()
    return -2 * L[2], -24 * L[4]


if 'g' in which:
    print('--- (g) PF3 near 0: V4(0) <= 2 V2(0)^2, ratio V4/(2 V2^2) (must be <= 1); laws of L(1/2+z,chi)/L(1/2,chi) completed')
    # zeta
    T = arb_series([arb(0), arb(1)], prec=6)
    w = (2 * T).exp()
    tot = None
    for n in range(1, 9):
        term = (2 * pi**2 * n**4 * (T * (arb(9) / 2)).exp() - 3 * pi * n**2 * (T * (arb(5) / 2)).exp()) * ((-pi * n**2 * w).exp())
        tot = term if tot is None else tot + term
    L = tot.log()
    v2, v4 = -2 * L[2], -24 * L[4]
    print(' zeta      V2=%.5f V4=%.5f ratio=%.5f' % (float(v2), float(v4), float(v4 / (2 * v2**2))))
    sq5 = arb(5).sqrt()
    kappa = ((10 - 2 * sq5).sqrt() - 2) / (sq5 - 1)
    dh = [arb(0), arb(1), kappa, -kappa, arb(-1)]
    v2, v4 = V2V4(lambda n: dh[n % 5], 5, True)
    print(' DH (q=5)  V2=%.5f V4=%.5f ratio=%.5f' % (float(v2), float(v4), float(v4 / (2 * v2**2))))
    for D in (-3, -4, -7, -8, -11, -15, -19, -23):
        q = -D
        v2, v4 = V2V4(lambda n: kron(D, n) if (math.gcd(n, q) == 1) else 0, q, True)
        print(' D=%-4d odd  V2=%.5f V4=%.5f ratio=%.5f' % (D, float(v2), float(v4), float(v4 / (2 * v2**2))))
    for D in (5, 8, 12, 13, 17, 21):
        q = D
        v2, v4 = V2V4(lambda n: kron(D, n) if (math.gcd(n, q) == 1) else 0, q, False)
        print(' D=%-4d even V2=%.5f V4=%.5f ratio=%.5f' % (D, float(v2), float(v4), float(v4 / (2 * v2**2))))
    sys.stdout.flush()

if 'h' in which:
    print('--- (h) smooth two sided truncations Q_M = sum_{n<=M} (4 pi^2 n^4 cosh(9t/2) - 6 pi n^2 cosh(5t/2)) exp(-2 pi n^2 cosh 2t)')
    ctx.prec = 200
    h = arb(1) / 64
    ts = [h * k for k in range(0, 3 * 64 + 1)]
    pi = arb.pi()

    def kernel(M):
        out = []
        for t in ts:
            s = arb(0)
            for n in range(1, M + 1):
                s += (4 * pi**2 * n**4 * (arb(9) / 2 * t).cosh() - 6 * pi * n**2 * (arb(5) / 2 * t).cosh()) * (-2 * pi * n**2 * (2 * t).cosh()).exp()
            out.append(s)
        return out

    def F(ker, th):
        # th acb: trapezoid for int_R Q(t) cos(th t) dt
        s = ker[0] * acb(1)
        for k, v in zip(ts[1:], ker[1:]):
            s += 2 * v * (th * k).cos()
        return s * h

    TH = 80.0
    d = 0.05
    for M in (1, 2, 3, 5):
        ker = kernel(M)
        # real zeros by sign change
        zs = []
        prev = F(ker, acb(0.2)).real
        th = 0.2
        while th < TH:
            th2 = th + d
            cur = F(ker, acb(th2)).real
            if (prev > 0 and cur < 0) or (prev < 0 and cur > 0):
                a, b = th, th2
                fa = prev
                for _ in range(12):
                    m = (a + b) / 2
                    fm = F(ker, acb(m)).real
                    if (fa > 0) == (fm > 0):
                        a = m
                        fa = fm
                    else:
                        b = m
                zs.append((a + b) / 2)
            prev = cur
            th = th2
        # zeros off the axis: winding on the boundary of (0.2,TH)x(0.05,4)
        def phase_box(x0, x1, y0, y1, step=0.05):
            pts = []
            n = int((x1 - x0) / step)
            pts += [complex(x0 + (x1 - x0) * j / n, y0) for j in range(n)]
            m = int((y1 - y0) / step)
            pts += [complex(x1, y0 + (y1 - y0) * j / m) for j in range(m)]
            pts += [complex(x1 - (x1 - x0) * j / n, y1) for j in range(n)]
            pts += [complex(x0, y1 - (y1 - y0) * j / m) for j in range(m)]
            total = 0.0
            prevv = None
            for p in pts + pts[:1]:
                v = F(ker, acb(p.real, p.imag))
                ang = math.atan2(float(v.imag), float(v.real))
                if prevv is not None:
                    dd = ang - prevv
                    while dd > math.pi:
                        dd -= 2 * math.pi
                    while dd < -math.pi:
                        dd += 2 * math.pi
                    total += dd
                prevv = ang
            return total / (2 * math.pi)
        w = phase_box(0.2, 60.0, 0.05, 3.0)
        print(' M=%d: %d real zeros in (0.2,%g); first: %s; winding in (0.2,60)x(0.05,3): %.3f' % (M, len(zs), TH, [round(z, 3) for z in zs[:8]], w), flush=True)

if 'i' in which:
    print('--- (i) first zero of L_n(s)=sum_h P(H_n=h) (sqrt(2/pi) h/sqrt(2n))^s, H_n the height of a uniform Dyck path (exact law)')
    from math import lgamma
    import mpmath as mp
    mp.mp.dps = 40
    for n in (50, 200, 800, 3200, 12800):
        hmax = int(12 * math.sqrt(n)) + 10
        # P(H<=h)=4^{n+1}/(C_n (h+2)) sum_{j=1}^{h+1} sin^2(j pi/(h+2)) cos^{2n}(j pi/(h+2)), C_n Catalan
        logC = lgamma(2 * n + 1) - lgamma(n + 1) - lgamma(n + 2)
        F = [mp.mpf(0)] * (hmax + 1)
        for hh in range(1, hmax + 1):
            x = np.pi * np.arange(1, hh + 2) / (hh + 2)
            c = np.cos(x)
            s = np.sin(x)
            with np.errstate(divide='ignore'):
                lg = 2 * n * np.log(np.abs(c) + 1e-300) + (n + 1) * math.log(4.0) - logC - math.log(hh + 2)
            F[hh] = mp.mpf(float(np.sum(s * s * np.exp(lg))))
        P = [F[i] - F[i - 1] for i in range(1, hmax + 1)]
        scale = mp.sqrt(2 / mp.pi) / mp.sqrt(2 * n)
        hs = [mp.mpf(i) for i in range(1, hmax + 1)]
        def G(s):
            return sum(p * mp.power(scale * hk, s) for p, hk in zip(P, hs) if p > 1e-30)
        s0 = mp.mpc(0.5, 14.134725)
        try:
            r = mp.findroot(G, s0, tol=1e-12, maxsteps=30, verify=False)
            print(' n=%6d: first zero %.8f + %.8f i; Re-1/2=%.3e; Im-gamma_1=%.3e' % (n, r.real, r.imag, r.real - 0.5, r.imag - 14.134725141734694), flush=True)
        except Exception as e:
            print(' n=%d: failed %s' % (n, e))

if 'j' in which:
    print('--- (j) kink of the folded truncation: J_M = Phi_M\'(0+), and the height T_M where 2|J_M|/(xi(1/2) th^2) = amplitude of phi_X')
    ctx.prec = 200
    pi = arb.pi()
    T = arb_series([arb(0), arb(1)], prec=3)
    w = (2 * T).exp()
    xi12 = arb('0.4971207781883141')
    for M in (1, 2, 3, 4, 5, 6):
        tot = None
        for n in range(1, M + 1):
            term = 2 * (2 * pi**2 * n**4 * (T * (arb(9) / 2)).exp() - 3 * pi * n**2 * (T * (arb(5) / 2)).exp()) * ((-pi * n**2 * w).exp())
            tot = term if tot is None else tot + term
        J = float(tot[1]) if abs(float(tot[1])) > 1e-300 else 0.0
        c = 2 * abs(tot[1]) / xi12
        # solve 2.2519 th^{7/4} exp(-pi th/4) = c/th^2
        a, b = 1.0, 2000.0
        logc = float(c.log())
        for _ in range(100):
            m = (a + b) / 2
            if math.log(2.2519) + 3.75 * math.log(m) - math.pi * m / 4 > logc:
                a = m
            else:
                b = m
        print(' M=%d J_M=%s predicted T_M=%.2f  (4(M+1)^2=%d)' % (M, tot[1].str(5), a, 4 * (M + 1)**2))
