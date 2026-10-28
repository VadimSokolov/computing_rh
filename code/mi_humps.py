# mi_humps.py: the numerical claims of Chapter 22 (Nick Polson's Miyatake chapter, prefix mi:), part G1: the hump criterion and the real zeros of
# Xi_P (Theorem humps, Theorem count, Table count with the breakdown height T_P^*, T_P^+, N(T), R_P, the 215 digit check of the humps of P = 7 in
# [481, 501], the breakdown formula, and the claim that the primes 2, 3, 5, 7 reproduce the first R_7 zeros of Xi within 0.4).
# Added for the book on 22 Oct 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. Xi, M, D_P, S(P), q as in code/mi_friable.py (imported). Negative humps I_j = (gamma_(2j-1), gamma_(2j)) of Xi (Xi(0) > 0, simple
# zeros); T_P^* solves M(T) = D_P(T) (the root nearest 4 q^2); T_P^+ solves D_P(T) = 0.63 M(T) T^(1/6) log T; N(T) counts the zeros of
# data/zeros_1700.txt up to T; R_P = number of positive real zeros of Xi_P = Xi + D_P.
# Method. Hump criterion: h(t) = log(-Xi(t)) - log D_P(t) is maximised on each hump by golden section (60 steps, 30 digits; h is concave there
# under RH, Theorem curv), the hump holds two zeros of Xi_P exactly when max h > 0; the zeros are then located by bisection on
# (gamma_(2j-1), t*) and (t*, gamma_(2j)). Full sign scan (P <= 5): Xi_P on a grid of step 0.05 and the sign changes counted, all humps up to
# T_P^+ + 25 are in the scan; the scan cannot see a pair closer than the step, so the hump criterion is run for every P as the second count.
# Error control. The hump criterion is run at 30 digits and every Table count number at 30 and 40 digits (the difference is logged); the
# P = 7 humps in [481, 501] are recomputed at 215 digits (task check215: maximum of |Xi| - D_P by golden section at 215 digits, both zeros of Xi_P by
# the secant method at 215 digits).
# Usage: python3 -I code/mi_humps.py TASK [--P P] [--part i --nparts n] [--lo a --hi b] [--out FILE] with TASK in table humps scan check215 merge
# Each call stays below 30 s; the parts are run as separate processes and merged by `merge`.
import sys, os, json, time, math, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mi_friable as F
from mpmath import mp, mpf, log, exp, pi, findroot, fabs, floor, ceil

ZS = None
def zs():
    global ZS
    if ZS is None:
        d = mp.dps; mp.dps = 80; ZS = F.zeros(1700); mp.dps = d
    return ZS

def Ncount(T):
    return sum(1 for z in zs() if z <= T)

def Tstar(P):
    q = F.least_prime_above(P)
    f = lambda T: log(F.Mfun(T)) - log(F.D(T, P))
    return findroot(f, (4 * q * q, 4 * q * q + 4), solver='secant', tol=mpf(10) ** -(mp.dps - 5), maxsteps=100)

def Tplus(P):
    q = F.least_prime_above(P)
    f = lambda T: log(F.D(T, P)) - log(mpf('0.63') * F.Mfun(T) * T ** (mpf(1) / 6) * log(T))
    return findroot(f, (4 * q * q, 4 * q * q + 4), solver='secant', tol=mpf(10) ** -(mp.dps - 5), maxsteps=100)

def hump_max(a, b, P, iters=60):
    h = lambda t: log(-F.Xi_hardy(t)) - log(F.D(t, P))
    g = (math.sqrt(5) - 1) / 2
    c = b - (b - a) * g; d = a + (b - a) * g; fc = h(c); fd = h(d)
    for _ in range(iters):
        if fc > fd:
            b, d, fd = d, c, fc; c = b - (b - a) * g; fc = h(c)
        else:
            a, c, fc = c, d, fd; d = a + (b - a) * g; fd = h(d)
    t = (a + b) / 2
    return t, h(t)

def bisect(f, a, b, n=40):
    fa = f(a)
    for _ in range(n):
        m = (a + b) / 2; fm = f(m)
        if (fm > 0) == (fa > 0): a, fa = m, fm
        else: b = m
    return (a + b) / 2

def humps_range(P):
    T = Tstar(P) if mp.dps else None
    n = Ncount(T + 30)
    return (n + 1) // 2 + 1

def task_table(args):
    out = {}
    for P in (2, 3, 5, 7):
        q = F.least_prime_above(P)
        mp.dps = 30; ts = Tstar(P); tp = Tplus(P)
        mp.dps = 40; ts2 = Tstar(P); tp2 = Tplus(P)
        out['P%d' % P] = {'q': q, 'T*': ts, 'T*(40 digits)': ts2, 'N(T*)': Ncount(ts), 'T+': tp, 'T+(40 digits)': tp2, 'N(T+)': Ncount(tp),
                          'formula 4q^2+(6/pi)log q+1.61': 4 * q * q + 6 / pi * log(q) + mpf('1.61'),
                          'formula with exact constant': 4 * q * q + 6 / pi * log(q) + 4 / pi * log((pi ** 2 + 4) / 4) + mpf('0.01'),
                          'T* - 4q^2': ts - 4 * q * q, 'reldiff 30/40 digits': fabs(ts - ts2) / ts}
    return out

def task_humps(args):
    P = int(args['P']); part = int(args.get('part', 0)); nparts = int(args.get('nparts', 1))
    mp.dps = 30
    J = humps_range(P)
    js = [j for j in range(1, J + 1) if (j - 1) % nparts == part]
    if 'js' in args: js = [int(x) for x in args['js'].split(',')]
    G = zs(); res = []
    for j in js:
        a, b = G[2 * j - 2], G[2 * j - 1]
        t, hv = hump_max(a, b, P, iters=40 if args.get('quick') else 60)
        r = {'j': j, 'gamma_a': a, 'gamma_b': b, 'tstar': t, 'maxh': hv, 'survives': bool(hv > 0)}
        if hv > 0:
            f = lambda x: F.Xi_hardy(x) + F.D(x, P)
            # near a zero of Xi the sum Xi + D_P is resolved only if the working digits exceed log10(M/D_P) (D_P is far below Xi's rounding
            # noise at low heights), so the zeros are located at 40 + log10(M(a)/D_P(a)) digits (the criterion itself is at 30 digits)
            lr = math.log10(F.Mfun(a) / F.D(a, P))
            if lr > 50:
                # the shift of the zero is below 1e-50 (first order: D_P/Xi' ~ D_P/M); the ends of the hump are the zeros to that accuracy
                r['z1'], r['z2'], r['dps_zeros'] = a, b, 'ends (shift < 1e-50)'
            else:
                dh = 70 + int(max(0, lr))
                mp.dps = dh
                a2, b2 = F.zeros(2 * J + 2)[2 * j - 2], F.zeros(2 * J + 2)[2 * j - 1]
                r['dps_zeros'] = dh
                r['z1'] = bisect(f, a2, t, 30); r['z2'] = bisect(f, t, b2, 30)
            mp.dps = 30
        res.append(r)
    return {'P': P, 'part': part, 'humps': res, 'J': J}

def task_scan(args):
    P = int(args['P']); lo = mpf(args['lo']); hi = mpf(args['hi']); step = mpf('0.05')
    mp.dps = 30
    n = int(round((hi - lo) / step)); prev = F.XiP(lo, P); ch = []
    for i in range(1, n + 1):
        t = lo + i * step; v = F.XiP(t, P)
        if (v > 0) != (prev > 0): ch.append(t - step / 2)
        prev = v
    return {'P': P, 'lo': lo, 'hi': hi, 'step': step, 'sign_changes': len(ch), 'at': ch}

def task_check215(args):
    P = 7; mp.dps = 215
    G = zs(); out = []
    lo, hi = mpf(481), mpf(501)
    js = [j for j in range(1, 200) if G[2 * j - 2] >= lo and G[2 * j - 1] <= hi]
    part = int(args.get('part', 0)); nparts = int(args.get('nparts', 1)); js = js[part::nparts]
    for j in js:
        a, b = G[2 * j - 2], G[2 * j - 1]
        mp.dps = 30
        t30, h30 = hump_max(a, b, P, iters=40)
        mp.dps = 215
        # refine: golden section at 215 digits on [t30 - 0.02, t30 + 0.02] of Delta = |Xi| - D (same sign criterion)
        lo_, hi_ = t30 - mpf('0.02'), t30 + mpf('0.02')
        d = lambda t: -F.Xi_hardy(t) - F.D(t, P)
        g = (mp.sqrt(5) - 1) / 2
        c = hi_ - (hi_ - lo_) * g; e = lo_ + (hi_ - lo_) * g; fc, fe = d(c), d(e)
        for _ in range(14):
            if fc > fe: hi_, e, fe = e, c, fc; c = hi_ - (hi_ - lo_) * g; fc = d(c)
            else: lo_, c, fc = c, e, fe; e = lo_ + (hi_ - lo_) * g; fe = d(e)
        tm = (lo_ + hi_) / 2; dm = d(tm)
        r = {'j': j, 'gamma_a': G[2 * j - 2], 'gamma_b': G[2 * j - 1], 't*(30 digits)': t30, 'maxh (30 digits)': h30, 't*(215 digits)': tm,
             'max(|Xi|-D) at 215 digits': dm, 'D at t*': F.D(tm, P), 'ratio |Xi|/D - 1': dm / F.D(tm, P), 'survives215': bool(dm > 0), 'survives30': bool(h30 > 0)}
        if dm > 0:
            f = lambda x: F.Xi_hardy(x) + F.D(x, P)
            z1 = findroot(f, (a + (tm - a) * mpf('0.1'), a + (tm - a) * mpf('0.9')), solver='anderson', tol=mpf(10) ** -200, maxsteps=200) if False else None
            mp.dps = 60
            zz1 = bisect(f, a, tm, 25); zz2 = bisect(f, tm, b, 25)
            mp.dps = 215
            z1 = findroot(f, (zz1, zz1 + mpf('1e-6')), solver='secant', tol=mpf(10) ** -400, maxsteps=40)
            z2 = findroot(f, (zz2, zz2 - mpf('1e-6')), solver='secant', tol=mpf(10) ** -400, maxsteps=40)
            r['z1'] = z1; r['z2'] = z2; r['|Xi_P(z)| residual'] = max(fabs(f(z1)), fabs(f(z2)))
        out.append(r)
    mp.dps = 30
    return {'P': 7, 'humps': out}

def task_merge(args):
    out = {}
    for P in (2, 3, 5, 7):
        hs = []
        for f in sorted(glob.glob(os.path.join(F.ROOT, 'data', 'mi_humps_P%d_p*.json' % P))):
            hs += json.load(open(f))['humps']
        hs.sort(key=lambda r: r['j'])
        mp.dps = 30
        surv = [r for r in hs if r['survives']]
        zl = sorted([mpf(r['z1']) for r in surv] + [mpf(r['z2']) for r in surv])
        G = zs(); R = len(zl)
        disp = [(fabs(zl[k] - G[k]), k + 1) for k in range(R)]
        dead = [r['j'] for r in hs if not r['survives']]
        # nearest-zero pairing
        pair = []
        for k, z in enumerate(zl):
            n = min(range(len(G)), key=lambda i: fabs(G[i] - z))
            pair.append(n + 1)
        out['P%d' % P] = {'humps examined': len(hs), 'js': (hs[0]['j'], hs[-1]['j']), 'surviving humps': len(surv), 'R_P (2 x surviving)': 2 * len(surv),
                          'dead humps j': dead, 'max |z_k - gamma_k| (index pairing)': max(disp), 'k of max': max(disp)[1],
                          'zeros displaced > 0.4 (index pairing)': [d[1] for d in disp if d[0] > mpf('0.4')],
                          'nearest-gamma index differs from k': [(k + 1, pair[k]) for k in range(R) if pair[k] != k + 1][:20],
                          'last surviving hump': (surv[-1]['z1'], surv[-1]['z2']) if surv else None,
                          'maxh margin of last survivor / first dead': (surv[-1]['maxh'] if surv else None, [r['maxh'] for r in hs if not r['survives']][:3]),
                          'max over humps of tstar-less: largest displacement of a zero from its own hump end': max([max(fabs(mpf(r['z1']) - mpf(r['gamma_a'])), fabs(mpf(r['z2']) - mpf(r['gamma_b']))) for r in surv]) if surv else None}
    return out

TASKS = {'merge': task_merge, 'table': task_table, 'humps': task_humps, 'scan': task_scan, 'check215': task_check215}

def tostr(o):
    if isinstance(o, dict): return {str(k): tostr(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [tostr(v) for v in o]
    if isinstance(o, (mpf, F.mpc)): return mp.nstr(o, 30 if mp.dps > 100 else 14)
    return o

if __name__ == '__main__':
    a = sys.argv[1:]; task = a[0]
    args = {a[i][2:]: a[i + 1] for i in range(1, len(a) - 1) if a[i].startswith('--') and a[i] not in ('--quick',)}
    args['quick'] = '--quick' in a
    t0 = time.time()
    r = TASKS[task](args)
    r = tostr(r); r['seconds'] = round(time.time() - t0, 1)
    txt = json.dumps(r, indent=1)
    if 'out' in args: open(args['out'], 'w').write(txt + '\n')
    print(txt[:3000])
