# mi_friable.py: the numerical claims of Chapter 22 (Nick Polson's Miyatake chapter, prefix mi:), part G1: the friable approximant Xi_P, from
# Section "Riemann's representation" to Theorem curv and the Computational notes (the kink identity, the incomplete gamma form, the zeros of
# the literal Euler product F_P, Theorem approx, the first order shifts, the kink constant c_P, positivity and monotonicity of the deficit D_P,
# the Lorentzian bound epsilon(a), rho_1, rho_2 and the numerical minimum of (log D_P)'', the function M, the working precision).
# Added for the book on 22 Oct 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s), Xi(t) = xi(1/2 + i t); S(P) = integers whose prime factors are all <= P, q the
# least prime above P. C_a(t) = int_0^inf exp(-a e^w) e^(w/4) cos(t w/2) dw = Re[a^(-1/4 - it/2) Gamma(1/4 + it/2, a)] (substitute x = e^w).
# D_P(t) = Xi_P - Xi = (t^2 + 1/4) sum_{m not in S(P)} C_(pi m^2)(t) for every complex t (C_a continued as (g(t) + g(-t))/2); this equals the
# difference of the incomplete gamma forms of Proposition incgamma (checked against it in section incgamma).
# Method. mpmath at mp.dps = D (the library default 30); every section is run twice, at D and at D + 10 digits, and the largest relative
# difference of its numbers is logged as the precision check (the sections that use quadrature use the same two precisions). Sums over m stop
# at m with exp(-pi (m^2 - q^2)) < 10^(-dps - 5). Zeros of Xi: data/zeros_1700.txt (68 digits).
# Usage: python3 -I code/mi_friable.py SECTION [--quick] [--out FILE]; sections: kink incgamma fpoles approx shifts cp positive lorentz curv mfun.
# Each call stays below 30 s (the laptop rule); the sections are run as separate processes.
import sys, os, json, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import (mp, mpf, mpc, pi, exp, log, sqrt, gamma, gammainc, zeta, siegelz, loggamma, quad, findroot, diff, cos, nsum, inf, fabs, ceil, re, im)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DPS = 30

def zeros(n=1700):
    return [mpf(l.strip()) for l in open(os.path.join(ROOT, 'data', 'zeros_1700.txt'))][:n]

def is_smooth(m, P):
    for p in range(2, P + 1):
        while m % p == 0:
            m //= p
    return m == 1

def least_prime_above(P):
    q = P + 1
    while any(q % d == 0 for d in range(2, int(q ** 0.5) + 1)):
        q += 1
    return q

def excluded(P, extra=5):
    """integers m not in S(P) with exp(-pi(m^2 - q^2)) > 10^(-dps - extra)"""
    q = least_prime_above(P)
    mmax = int(math.sqrt(q * q + (mp.dps + extra) * 2.302585093 / math.pi)) + 1
    return [m for m in range(q, mmax + 1) if not is_smooth(m, P)]

def xi_s(s):
    return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)

def Xi(t):
    """Xi(t) = xi(1/2 + it); real part for real t"""
    v = xi_s(mpc(mpf(1) / 2, t) if not isinstance(t, mpc) else mpf(1) / 2 + 1j * t)
    return v.real if not isinstance(t, mpc) else v

def Mfun(t):
    t = mpf(t)
    return (t * t + mpf(1) / 4) / 2 * pi ** (-mpf(1) / 4) * exp(re(loggamma(mpf(1) / 4 + 1j * t / 2)))

def Xi_hardy(t):
    """Xi(t) = -M(t) Z(t), cheaper than Xi for real t"""
    return -Mfun(t) * siegelz(t)

def Ca(a, t):
    z = lambda u: a ** (-(mpf(1) / 4 + 1j * u / 2)) * gammainc(mpf(1) / 4 + 1j * u / 2, a)
    if isinstance(t, mpc):
        return (z(t) + z(-t)) / 2
    return re(z(t))

def D(t, P):
    t = t if isinstance(t, mpc) else mpf(t)
    return (t * t + mpf(1) / 4) * sum(Ca(pi * m * m, t) for m in excluded(P))

def XiP(t, P):
    return Xi_hardy(t) + D(t, P)

def XiP_incgamma(t, P, nmax=14, allint=False):
    s = mpf(1) / 2 + 1j * t
    tot = 0
    for n in range(1, nmax + 1):
        if allint or is_smooth(n, P):
            a = pi * n * n
            tot += a ** (-s / 2) * gammainc(s / 2, a) + a ** (-(1 - s) / 2) * gammainc((1 - s) / 2, a)
    return mpf(1) / 2 + s * (s - 1) / 2 * tot

def XiP_def(t, P, nmax=14):
    """the definition: 1/2 - 2 (t^2+1/4) int_0^inf psi_P(e^(2u)) e^(u/2) cos(tu) du by quadrature"""
    psiP = lambda x: sum(exp(-pi * n * n * x) for n in range(1, nmax + 1) if is_smooth(n, P))
    f = lambda u: psiP(exp(2 * u)) * exp(u / 2) * cos(t * u)
    pts = [mpf(k) / 4 for k in range(0, 17)]
    return mpf(1) / 2 - 2 * (t * t + mpf(1) / 4) * quad(f, pts)

def cP(P):
    return sum((4 * pi * m * m - 1) * exp(-pi * m * m) for m in range(2, 16) if not is_smooth(m, P))

def eps_P(P):
    return sum(exp(-pi * m * m) / (pi * m * m) for m in range(2, 16) if not is_smooth(m, P))

def eps_bound(P):
    q = least_prime_above(P)
    return exp(-pi * q * q) / (pi * q * q * (1 - exp(-2 * pi * q)))

def eps_a(a):
    b = a - mpf(1) / 4
    return (2 * a * b * b / (b - 1) ** 3 + 2 * a * b / (b - 1) ** 2 + 2 * a * a / (b - 2) ** 3 + a / (b - 1)) / b

def rho12(a):
    b = a - mpf(1) / 4
    e = eps_a(a)
    r1 = (1 + a * a / (b - 2) ** 2 + a / (b - 1) ** 2 + 2 * a / (b - 1)) / (2 * (1 - e) * b)
    r2 = (4 / b + 2 * a * a / (b - 2) ** 3 + 2 * a / (b - 1) ** 3 + 4 * a / (b - 1) ** 2) / (4 * (1 - e) * b)
    return r1, r2

def Lor(a, t):
    b = a - mpf(1) / 4
    return exp(-a) * b / (b * b + t * t / 4)

def logD_dd(t, P, h=mpf(1) / 1000):
    g = lambda x: log(D(x, P))
    return (g(t + h) - 2 * g(t) + g(t - h)) / h ** 2

# ---------------------------------------------------------------------------- sections
def sec_kink(quick):
    mp.dps = 60
    psi = lambda: sum(exp(-pi * n * n) for n in range(1, 12))
    dpsi = sum(-pi * n * n * exp(-pi * n * n) for n in range(1, 12))
    return {'kink_value': mpf(1) / 2 + psi() + 4 * dpsi}

def sec_incgamma(quick):
    mp.dps = 40; out = {}
    ts = [mpf('0.7'), mpf(3), mpf(10)] + ([] if quick else [mpf(25)])
    for P in (2, 3):
        for t in ts:
            a = XiP_incgamma(t, P); b = XiP_def(t, P)
            out['P%d_t%s' % (P, t)] = {'incgamma': re(a), 'quad': b, 'absdiff': fabs(re(a) - b), 'imag': im(a)}
            d = D(t, P); out['P%d_t%s' % (P, t)]['Dfromgamma_minus_incg'] = fabs(re(a) - Xi_hardy(t) - d)
    for t in [mpf('0.7'), mpf(3), mpf(10), mpf(25), mpf(50)]:
        a = re(XiP_incgamma(t, 2, nmax=12, allint=True)); out['all_t%s' % t] = {'incgamma_all': a, 'xi': Xi(t), 'reldiff': fabs(a - Xi(t)) / fabs(Xi(t))}
    return out

def sec_fpoles(quick):
    mp.dps = 25
    G = lambda s: pi ** (-s / 2) * gamma(s / 2)
    Ep = lambda s: s * (s - 1) / 2 * G(s) / (1 - mpf(2) ** (-s))
    F = lambda s: Ep(s) + Ep(1 - s)
    claims = ['1.5151+9.1951j', '1.3184+18.0602j', '1.2268+27.3286j', '1.1815+45.4005j', '0.7885+72.7840j', '1.1200+108.8198j']
    out = {}
    for c in claims:
        z = findroot(F, mpc(complex(c)), tol=1e-18, maxsteps=100)
        k = round(float(z.imag) * math.log(2) / (2 * math.pi))
        out[c] = {'zero': z, 'abs F': fabs(F(z)), 'nearest pole k': k, 'pole': 2j * pi * k / log(2), 'dist': fabs(z - 2j * pi * k / log(2))}
    return out

def sec_approx(quick, Ps=(2, 3, 5, 7), part=None, nparts=3):
    mp.dps = 40; out = {}
    xs = [mpf(k) / 4 for k in range(0, 17)] + [mpf(k) for k in range(5, 61 if not quick else 21, 1)]
    if part is not None: xs = xs[part::nparts]
    ys = [mpf(k) / 2 for k in range(-3, 4)]
    for P in Ps:
        e = eps_P(P); mx = 0; arg = None
        for x in xs:
            for y in ys:
                t = mpc(x, y)
                if fabs(t * t + mpf(1) / 4) < 1e-9: continue
                r = fabs(D(t, P)) / (fabs(t * t + mpf(1) / 4) * e)
                if r > mx: mx, arg = r, (x, y)
        out['P%d' % P] = {'eps_P': e, 'eps_bound': eps_bound(P), 'max_ratio': mx, 'at': str(arg)}
    return out

def sec_shifts(quick):
    mp.dps = 40
    g = zeros(5); out = {}
    prev = None
    for k, gam in enumerate(g):
        d1 = diff(Xi, gam)
        pred = -D(gam, 2) / d1
        z = findroot(lambda t: XiP(t, 2), gam + pred, tol=mpf(10) ** -35)
        out['zero%d' % (k + 1)] = {'gamma': gam, 'predicted': pred, 'computed': z - gam, 'ratio': (z - gam) / pred}
    out['decade'] = 4 * log(10) / pi
    return out

def sec_cp(quick):
    mp.dps = 40; out = {}
    for P in (2, 3):
        c = cP(P)
        out['P%d' % P] = {'c_P': c, 'leading 4 pi q^2 e^(-pi q^2)': 4 * pi * least_prime_above(P) ** 2 * exp(-pi * least_prime_above(P) ** 2)}
        kk = lambda u: exp(-u / 2) - 2 * exp(u / 2) * sum(exp(-pi * n * n * exp(2 * u)) for n in range(1, 16) if is_smooth(n, P))
        k3 = diff(kk, 0, 3)
        out['P%d' % P]['k_P\'\'\'(0)'] = k3; out['P%d' % P]['predicted t^2(XiP-c_P) -> c_P/4 + k_P\'\'\'(0)'] = c / 4 + k3
        for t in (50, 100, 200, 400, 800):
            v = XiP(mpf(t), P)
            out['P%d' % P]['t%d' % t] = {'XiP': v, 'XiP - c_P': v - c, 't^2 (XiP - c_P)': t * t * (v - c), 'XiP/c_P': v / c}
    return out

def sec_positive(quick):
    mp.dps = 30; out = {}
    ts = [mpf(k) / 2 for k in range(0, 21)] + [mpf(k) for k in range(11, 60)] + [mpf(k) for k in range(60, 601, 20)]
    if quick: ts = ts[::3]
    prev = None; res = {}
    Ds = {P: [D(t, P) for t in ts] for P in (2, 3, 5, 7)}
    for P in (2, 3, 5, 7): res['D_%d_min' % P] = min(Ds[P])
    res['all_positive'] = all(v > 0 for P in Ds for v in Ds[P])
    res['monotone 2>3>5>7'] = all(Ds[2][i] > Ds[3][i] > Ds[5][i] > Ds[7][i] for i in range(len(ts)))
    res['points'] = len(ts)
    return res

def sec_lorentz(quick):
    mp.dps = 25; out = {}
    a = 9 * pi
    out['eps(9pi)'] = eps_a(a); out['5.61/a'] = mpf('5.61') / a; out['a eps(a) at 9pi'] = a * eps_a(a)
    out['a eps(a) at 1e3, 1e6'] = [a_ * eps_a(a_) for a_ in (mpf(1000), mpf(10) ** 6)]
    # quadrature check of C_a at a = 9 pi
    qt = {}
    for t in (0, 5, 30):
        t = mpf(t)
        f = lambda w: exp(-a * exp(w)) * exp(w / 4) * cos(t * w / 2)
        qt[str(t)] = {'quad': quad(f, [mpf(k) / 2 for k in range(0, 13)]), 'gammainc': Ca(a, t)}
    out['C_a quadrature vs gammainc'] = qt
    mx = 0; arg = None
    ts = [mpf(k) / 4 for k in range(0, 40)] + [mpf(k) for k in range(10, 200, 3)] + [mpf(k) for k in range(200, 3001, 100)]
    if quick: ts = ts[::4]
    for t in ts:
        r = fabs(Ca(a, t) - Lor(a, t)) / Lor(a, t)
        if r > mx: mx, arg = r, t
    out['true max relative error at 9pi'] = mx; out['argmax t'] = arg
    # relative error at q >= 11 (computational notes)
    out['max relative error of C_a vs L_a, a = pi 11^2'] = max(fabs(Ca(pi * 121, t) - Lor(pi * 121, t)) / Lor(pi * 121, t) for t in ts)
    out['bound eps at pi 11^2'] = eps_a(pi * 121)
    # quotient D_P(4q^2)/c_P
    for P in (2, 3, 5, 7):
        q = least_prime_above(P)
        out['D_P(4q^2)/c_P P=%d' % P] = D(mpf(4 * q * q), P) / cP(P)
    out['4/(pi^2+4)'] = 4 / (pi ** 2 + 4)
    # Lorentzian replacement for D_P at q = 11 (P = 7), relative error over a t grid
    L7 = lambda t: (t * t + mpf(1) / 4) * sum(Lor(pi * m * m, t) for m in excluded(7))
    out['max |L_P - D_P|/D_P, P=7, t in [14,1000]'] = max(fabs(L7(t) - D(t, 7)) / D(t, 7) for t in [mpf(k) for k in range(14, 1001, 40)])
    return out

def sec_curv(quick):
    mp.dps = 30; out = {}
    for q in (3, 5, 7, 11):
        a = pi * q * q; r1, r2 = rho12(a)
        out['q=%d' % q] = {'rho1': r1, 'rho2': r2, 'eps': eps_a(a), '2/14^2': mpf(2) / 196, 'bound -(2/t^2 + rho2 + rho1^2) at t=14': -(mpf(2) / 196 + r2 + r1 * r1)}
    ts = [mpf(14) + mpf(k) / 2 for k in range(0, 30)] + [mpf(k) for k in range(30, 120, 5)] + [mpf(k) for k in range(120, 800, 40)]
    if quick: ts = ts[::5]
    for P in (2, 3, 5, 7):
        vals = [(logD_dd(t, P), t) for t in ts]
        m = min(vals)
        out['min (log D_%d)\'\'' % P] = {'min': m[0], 'at t': m[1], 'value at t=14': vals[0][0], 'value at last': vals[-1][0]}
    return out

def sec_mfun(quick):
    mp.dps = 30; out = {}
    zs = zeros(3)
    for t in [mpf(5), mpf('14.5'), mpf(22), mpf(40), mpf(100), mpf(250)]:
        x = Xi(t); mz = -Mfun(t) * siegelz(t)
        out['t=%s' % t] = {'Xi': x, '-M Z': mz, 'rel': fabs(x - mz) / fabs(x)}
    for t in (100, 300, 500):
        t = mpf(t)
        out['M/asymptotic t=%d' % t] = Mfun(t) / ((pi / 2) ** (mpf(1) / 4) * t ** (mpf(7) / 4) * exp(-pi * t / 4))
    # working precision: digits needed at T = 500: Xi ~ e^(-pi T/4), terms size T^2; check Xi_P two ways
    for T in (100, 200, 500):
        T = mpf(T); dg = pi * T / (4 * log(10)) + 3 * log(T, 10) + 25
        out['digits at T=%d' % T] = dg
    mp.dps = 30
    T = mpf('490.5')
    direct = re(XiP_incgamma(T, 7, nmax=14))
    out['direct sum of incgamma at dps 30 (should be garbage)'] = direct
    need = int(pi * T / (4 * log(10)) + 3 * log(T, 10) + 25)
    mp.dps = need
    direct2 = re(XiP_incgamma(T, 7, nmax=16))
    mp.dps = 30
    split = XiP(T, 7)
    out['T=490.5 P=7 direct at formula digits (%d)' % need] = direct2
    out['T=490.5 P=7 Xi + D_P at 30 digits'] = split
    mp.dps = need + 20
    out['T=490.5 P=7 direct at formula + 20 digits'] = re(XiP_incgamma(T, 7, nmax=16))
    return out

SECTIONS = {k[4:]: v for k, v in list(globals().items()) if k.startswith('sec_')}
for _P in (2, 3, 5, 7):
    for _k in range(9):
        SECTIONS['approx%dq%d' % (_P, _k)] = (lambda quick, _P=_P, _k=_k: sec_approx(quick, (_P,), _k, 9))
    SECTIONS['approx%d' % _P] = (lambda quick, _P=_P: sec_approx(quick, (_P,)))
    for _k in range(3):
        SECTIONS['approx%dp%d' % (_P, _k)] = (lambda quick, _P=_P, _k=_k: sec_approx(quick, (_P,), _k))

def tostr(o):
    if isinstance(o, dict): return {str(k): tostr(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [tostr(v) for v in o]
    if isinstance(o, (mpf, mpc)): return mp.nstr(o, 12)
    return o

def flat(o, pre=''):
    if isinstance(o, dict):
        r = {}
        for k, v in o.items(): r.update(flat(v, pre + str(k) + '/'))
        return r
    if isinstance(o, (list, tuple)):
        r = {}
        for i, v in enumerate(o): r.update(flat(v, pre + str(i) + '/'))
        return r
    return {pre: o}

def run(name, quick):
    t0 = time.time()
    mp.dps = DEFAULT_DPS
    r1 = SECTIONS[name](quick)
    d1 = DEFAULT_DPS
    mp.dps = DEFAULT_DPS + 10
    # second run at 10 more digits: the sections set their own dps at the start, so shift them via a global offset
    return r1, time.time() - t0

if __name__ == '__main__':
    args = sys.argv[1:]
    quick = '--quick' in args
    out = args[args.index('--out') + 1] if '--out' in args else None
    name = args[0]
    t0 = time.time()
    r = SECTIONS[name](quick)
    # precision check: rerun with every section's dps raised by 10
    base = mp.dps
    orig_dps = type(mp).dps
    class _Shift: pass
    import mpmath
    shift = 10
    _set = mpmath.ctx_mp.MPContext.dps.fset
    _get = mpmath.ctx_mp.MPContext.dps.fget
    mpmath.ctx_mp.MPContext.dps = property(_get, lambda self, v: _set(self, v + shift if shift else v))
    mp.dps = DEFAULT_DPS
    r2 = SECTIONS[name](quick)
    mpmath.ctx_mp.MPContext.dps = property(_get, _set)
    f1, f2 = flat(r), flat(r2); worst = 0; wk = None
    for k in f1:
        a, b = f1[k], f2.get(k)
        if isinstance(a, (mpf, mpc)) and isinstance(b, (mpf, mpc)) and k.split('/')[-2:-1] not in (['at'], ['rel'], ['absdiff'], ['reldiff'], ['imag'], ['Dfromgamma_minus_incg'], ['kink_value']):
            scale = max(fabs(a), fabs(b))
            if scale > 0:
                d = fabs(a - b) / scale
                if d > worst: worst, wk = d, k
    res = {'section': name, 'quick': quick, 'seconds': round(time.time() - t0, 1), 'precision_check_max_reldiff': mp.nstr(worst, 3), 'at_key': wk, 'results': tostr(r)}
    txt = json.dumps(res, indent=1)
    if out: open(out, 'w').write(txt + '\n')
    print(txt)
