# pg3_fejer.py: the Polya law as a triangle (Fejer) mixture (Proposition fm:prop:triangle and the paragraph after it, Chapter 13, merged
# from note X3; placeholders N4a, N4b). Added for the book in October 2026, written from the mathematics (the collaborator's scripts were
# not read, run or copied).
# Conventions (the book's, Titchmarsh's variable, eq:ch2:Phi). Phi(u) = 2 sum_n (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u))
# = 2 e^(u/2) sum_n (2 y_n^2 - 3 y_n) e^(-y_n), y_n = pi n^2 e^(2u); Xi(0) = xi(1/2) = int_R Phi; p = Phi/Xi(0) the Polya density.
# The mixing measure of the Polya law over the triangle f_2 is r^2 p''(r) dr on (0, inf), of mass one.
# (a) CERTIFIED (python-flint, Hopper; skipped without flint). Phi''(u) = (1/2) e^(u/2) S(u), S(u) = sum_n R(y_n) e^(-y_n),
#     R(y) = 32 y^4 - 224 y^3 + 330 y^2 - 75 y. Terms n <= 6 in ball arithmetic (200 bits), the terms n >= 7 bounded by [0, T] with
#     T = sum_{n=7}^{30} 32 pi^4 n^8 e^(-pi n^2) + 1e-300 (0 < R(y) <= 32 y^4 for y >= 7 and y^4 e^(-y) decreases for y > 4, so the term of
#     index n is largest at u = 0; the terms n > 30 are below e^(-2800) each and below 1e-300 together); the derivative S' = sum
#     2 y (R'(y) - R(y)) e^(-y) has its tail bounded by 70 pi^5 n^10 e^(-pi n^2) in the same way (|2 y (R' - R)| <= 70 y^5 for y >= 150).
#     On [0, u_b], u_b = (1/2) log(7/pi) = 0.4007: with the derivative S' > 0 on [a, b] = [u* - eps, u* + eps] (eps = 1e-3, centred
#     on the root found by bisection), S < 0 on [0, a] and S > 0 on [b, u_b] by adaptive interval subdivision. Then S has exactly one
#     zero on [0, u_b]; it is bisected to width 1e-25 using the monotonicity of S on [a, b]; the ball is N4a. (Phi'' > 0 beyond u_b is
#     proved in the text.)
# (d) The constants of the proof of (ii) in ball arithmetic: R(pi) < -800, e^(-pi) > 0.043, 32 pi^4 sum_{n>=2} n^8 e^(-pi n^2) < 3 (the
#     terms n > 30 below 1e-300).
# (b) mpmath, 30 digits: Xi(0) = int_R Phi against xi(1/2); u* by findroot of the formula for Phi'' (also checked against mp.diff of Phi);
#     the mass M = int_0^u* r^2 p'' dr = u*^2 p'(u*) - 2 u* p(u*) + 2 int_0^u* p (N4b) against direct quadrature, and int_0^inf r^2 p'' = 1.
# (c) mpmath: the two identities of the text, p(x) = int_x^inf (r - x) p''(r) dr and E e^(i theta X) = (4/theta^2) int_0^inf p''(r) sin^2(theta r/2) dr
#     = 2 int_0^inf p(x) cos(theta x) dx at x in {0.04, 0.2, 0.6} and theta in {2.5, 10, 20}, and at theta = gamma_1 (unrounded) as an absolute residual.
# Options: --quick (laptop, (b) and (c) at 20 digits), --out FILE.
import sys, json, argparse, time
import mpmath as mp

def Phi(u):
    u = mp.mpf(u); s = mp.mpf(0); n = 1
    while True:
        y = mp.pi * n * n * mp.exp(2 * u)
        t = (2 * y * y - 3 * y) * mp.exp(-y)
        s += t
        if y > mp.mp.dps * 2.4 + 40 and abs(t) < mp.mpf(10) ** (-mp.mp.dps - 20): break
        n += 1
    return 2 * mp.exp(u / 2) * s

def R(y): return 32 * y ** 4 - 224 * y ** 3 + 330 * y ** 2 - 75 * y

def Phi2(u):
    u = mp.mpf(u); s = mp.mpf(0); n = 1
    while True:
        y = mp.pi * n * n * mp.exp(2 * u)
        t = R(y) * mp.exp(-y)
        s += t
        if y > mp.mp.dps * 2.4 + 40: break
        n += 1
    return mp.exp(u / 2) * s / 2

def cert(log, res):
    from flint import arb, ctx
    ctx.prec = 200
    pi = arb.pi()
    def Rb(y): return 32 * y ** 4 - 224 * y ** 3 + 330 * y ** 2 - 75 * y
    def Rpb(y): return 128 * y ** 3 - 672 * y ** 2 + 660 * y - 75
    T = sum((32 * pi ** 4 * n ** 8 * (-pi * n * n).exp() for n in range(7, 31)), arb(0)) + arb('1e-300')
    Tp = sum((70 * pi ** 5 * n ** 10 * (-pi * n * n).exp() for n in range(7, 31)), arb(0)) + arb('1e-300')
    def S(u):
        s = arb(0)
        for n in range(1, 7):
            y = pi * n * n * (2 * u).exp()
            s += Rb(y) * (-y).exp()
        return s + arb(0).union(T)
    def Sp(u):
        s = arb(0)
        for n in range(1, 7):
            y = pi * n * n * (2 * u).exp()
            s += 2 * y * (Rpb(y) - Rb(y)) * (-y).exp()
        return s + (-Tp).union(Tp)
    ub = ((arb(7) / pi).log()) / 2
    def all_sign(f, lo, hi, sign):
        stack = [(lo, hi, 0)]; cnt = 0
        while stack:
            a, b, d = stack.pop()
            v = f(a.union(b)); cnt += 1
            if (sign < 0 and v.upper() < 0) or (sign > 0 and v.lower() > 0): continue
            if d > 60: return False, cnt
            m = (a + b) / 2
            stack.append((a, m, d + 1)); stack.append((m, b, d + 1))
        return True, cnt
    # root by bisection on points (S at points has the tail interval; use the midpoint value without the tail for the search)
    def S0(u):
        return S(u).mid()
    lo, hi = arb(0), ub
    assert S(lo).upper() < 0
    for _ in range(80):
        m = (lo + hi) / 2
        if S(m).mid() > 0: hi = m
        else: lo = m
    root = (lo + hi) / 2
    eps = arb('1e-3')
    a, b = root - eps, root + eps
    ok1, c1 = all_sign(S, arb(0), a, -1)
    ok2, c2 = all_sign(Sp, a, b, +1)
    ok3, c3 = all_sign(S, b, ub, +1)
    log('(a) S < 0 on [0, u*-1e-3]: %s (%d boxes); S\' > 0 on [u*-1e-3, u*+1e-3]: %s (%d boxes); S > 0 on [u*+1e-3, (1/2)log(7/pi)]: %s (%d boxes)' % (ok1, c1, ok2, c2, ok3, c3))
    sa, sb = S(a), S(b)
    log('(a) S(u*-eps) = %s, S(u*+eps) = %s' % (sa.str(8), sb.str(8)))
    assert sa.upper() < 0 and sb.lower() > 0
    lo, hi = a, b
    for _ in range(90):
        m = (lo + hi) / 2
        v = S(m)
        if v.lower() > 0: hi = m
        elif v.upper() < 0: lo = m
        else: break
    ball = lo.union(hi)
    log('(a) exactly one zero of Phi\'\' on [0, (1/2)log(7/pi) = %s]: %s; u* in %s' % (ub.str(10, radius=False), ok1 and ok2 and ok3, ball.str(20)))
    res['a'] = dict(single_zero_certified=bool(ok1 and ok2 and ok3), N4a_ball=ball.str(25), N4a_mid=ball.mid().str(15, radius=False), radius=float(ball.rad()),
                    boxes=[c1, c2, c3], u_b=ub.str(15, radius=False))
    # (d)
    Rpi = Rb(pi); e = (-pi).exp()
    ssum = 32 * pi ** 4 * sum((n ** 8 * (-pi * n * n).exp() for n in range(2, 31)), arb(0)) + arb('1e-300')
    log('(d) R(pi) = %s (< -800: %s); exp(-pi) = %s (> 0.043: %s); 32 pi^4 sum_(n>=2) n^8 exp(-pi n^2) = %s (< 3: %s)'
        % (Rpi.str(10), Rpi.upper() < -800, e.str(10), e.lower() > arb('0.043'), ssum.str(10), ssum.upper() < 3))
    res['d'] = dict(R_pi=Rpi.str(10), exp_minus_pi=e.str(10), tail=ssum.str(10), ok=bool(Rpi.upper() < -800 and e.lower() > arb('0.043') and ssum.upper() < 3))

def laptop(dps, log, res):
    mp.mp.dps = dps
    X0 = 2 * mp.quad(Phi, [0, 0.5, 1, 2, 4])
    s = mp.mpf(1) / 2
    xi_half = s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2
    log('Xi(0) = int_R Phi = %s, xi(1/2) = %s, difference %s' % (mp.nstr(X0, 20), mp.nstr(xi_half, 20), mp.nstr(abs(X0 - xi_half), 3)))
    d = mp.diff(Phi, mp.mpf('0.1'), 2)
    log('formula for Phi\'\' against mp.diff at u = 0.1: %s vs %s' % (mp.nstr(Phi2(mp.mpf('0.1')), 15), mp.nstr(d, 15)))
    us = mp.findroot(Phi2, mp.mpf('0.2372'))
    log('u* (findroot) = %s ; Phi\'\'(0) = %s < 0; u*_dBN = u*/2 = %s' % (mp.nstr(us, 20), mp.nstr(Phi2(0), 8), mp.nstr(us / 2, 10)))
    p = lambda r: Phi(r) / X0
    p1 = mp.diff(Phi, us, 1) / X0
    M = us ** 2 * p1 - 2 * us * p(us) + 2 * mp.quad(p, [0, us])
    M2 = mp.quad(lambda r: r * r * Phi2(r) / X0, [0, us / 2, us])
    tot = mp.quad(lambda r: r * r * Phi2(r) / X0, [0, us, 0.5, 1, 2, 4])
    log('mass of r^2 p\'\'(r) dr on (0, u*): %s (direct quadrature %s), on (u*, inf): %s, total %s' % (mp.nstr(M, 20), mp.nstr(M2, 20), mp.nstr(1 - M, 20), mp.nstr(tot, 20)))
    Mr = mp.quad(lambda r: r * Phi2(r) / X0, [0, us])
    log('measure r p\'\'(r) dr (the note\'s omega p\'\'(omega) d omega is twice it): mass %s on (0, u*) and %s on (u*, inf), total p(0) = %s; the note\'s masses are %s and %s'
        % (mp.nstr(Mr, 10), mp.nstr(Phi(0) / X0 - Mr, 10), mp.nstr(Phi(0) / X0, 10), mp.nstr(2 * Mr, 10), mp.nstr(2 * (Phi(0) / X0 - Mr), 10)))
    res['b'] = dict(Xi0=mp.nstr(X0, 25), xi_half=mp.nstr(xi_half, 25), u_star=mp.nstr(us, 25), N4b=mp.nstr(M, 20), N4b_direct=mp.nstr(M2, 20),
                    outer_mass=mp.nstr(1 - M, 20), total=mp.nstr(tot, 20), p0=mp.nstr(Phi(0) / X0, 15), p_note_variable_p0=mp.nstr(2 * Phi(0) / X0, 15))
    # (c)
    pp = lambda r: Phi2(r) / X0
    cres = []
    for x in ['0.04', '0.2', '0.6']:
        x = mp.mpf(x)
        rhs = mp.quad(lambda r: (r - x) * pp(r), [x, 0.5, 1, 2, 4])
        cres.append(('p(x)', float(x), mp.nstr(p(x), 18), mp.nstr(rhs, 18)))
        log('(c) x = %s: p(x) = %s, int_x^inf (r-x) p\'\' dr = %s, difference %s' % (mp.nstr(x, 3), mp.nstr(p(x), 18), mp.nstr(rhs, 18), mp.nstr(abs(p(x) - rhs), 3)))
    g1 = mp.zetazero(1).imag
    for th in [mp.mpf('2.5'), mp.mpf(10), mp.mpf(20), g1]:
        lhs = 4 / th ** 2 * mp.quad(lambda r: pp(r) * mp.sin(th * r / 2) ** 2, mp.linspace(0, 4, 41))
        rr = 2 * mp.quad(lambda x: p(x) * mp.cos(th * x), mp.linspace(0, 4, 41))
        s = mp.mpf(1) / 2 + 1j * th
        xiv = mp.re(s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2) / xi_half
        cres.append(('cf', mp.nstr(th, 8), mp.nstr(lhs, 18), mp.nstr(rr, 18), mp.nstr(xiv, 18)))
        log('(c) theta = %s: (4/theta^2) int p\'\' sin^2 = %s, 2 int p cos = %s, xi(1/2+i theta)/xi(1/2) = %s; residuals %s, %s' % (mp.nstr(th, 12), mp.nstr(lhs, 18), mp.nstr(rr, 18), mp.nstr(xiv, 18), mp.nstr(abs(lhs - rr), 3), mp.nstr(abs(lhs - xiv), 3)))
    res['c'] = cres

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='pg3_fejer.json')
    a = ap.parse_args()
    L = []
    def log(s):
        print(s); sys.stdout.flush(); L.append(s)
    res = {}
    laptop(20 if a.quick else 30, log, res)
    try:
        import flint
        cert(log, res)
    except ImportError:
        log('python-flint not installed: parts (a) and (d) not run (PENDING, Hopper)')
    json.dump(res, open(a.out, 'w'), indent=1, default=str)
    open(a.out.replace('.json', '.log'), 'w').write('\n'.join(L) + '\n')
