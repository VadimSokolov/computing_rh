# mi_dh.py: the Davenport Heilbronn tests of Chapter 22 (Nick Polson's Miyatake chapter, Section "The ladder and the quasi Riemann abscissa"):
# the definition and functional equation of Lambda, its zeros off the line, the cosh ladder of Lambda (nonreal pair, Im^2 against c,
# collision, box count), the order of the zeros of Re and Im Lambda(c + i t) (paired form), the real zero counts of Re Lambda(c + i t) on
# (0.5, 200) against the zeros to the right of c (Theorem density), and Table thr (critical rung b_c and critical heat time lambda_c).
# Claims supported: chapter lines C:591 (first off line zero, Im^2 against c, collision near 0.798, box count 4), C:605 to C:609
# (ABAABBAB, ABABABAB), C:620 to C:627 (130, 130, 126, 124), C:635 to C:652 (Table thr, local model). Added for the book in October 2026.
# Conventions. chi = the character mod 5 with chi(1) = 1, chi(2) = i, chi(3) = -i, chi(4) = -1 (odd). L(s, chi) = 5^(-s) sum_a chi(a) zeta(s, a/5)
# (Hurwitz). f(s) = ((1 - i kappa)/2) L(s, chi) + ((1 + i kappa)/2) L(s, chibar) = 5^(-s) [ (h1 - h4) + kappa (h2 - h3) ], h_a = zeta(s, a/5),
# kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1). Lambda(s) = (5/pi)^(s/2) Gamma((s + 1)/2) f(s). The task defn checks Lambda(s) = Lambda(1 - s) and
# that Lambda is real on Re s = 1/2, and solves the (linear in kappa) functional equation residual for kappa, which pins the combination.
# The ladder: G(z) = Lambda(1/2 + i z) is real on R, F_b(z) = (G(z + i b) + G(z - i b))/2 = (Lambda(c + i z) + Lambda(c - i z))/2 with c = 1/2 + b
# (Lambda(1 - s) = Lambda(s)); on R, F_b(x) = Re Lambda(c + i x). Heat flow: Lambda(1/2 + i z) = int phi(u) e^(i z u) du with
# phi(u) = 2 (5/pi)^(-1/2) e^(3u/2) psi(e^(2u)), psi(y) = sum_n n (a chi(n) + abar chibar(n)) exp(-pi n^2 y/5) (the theta kernels of the two L functions,
# Gamma((s+1)/2) (pi/5)^(-(s+1)/2) L(s, chi) = int_0^inf psi_chi(y) y^((s+1)/2) dy/y, y = e^(2u)); the flow is phi(u) e^(lambda u^2), which acts on
# the entire function by e^(-lambda d^2/dz^2). Evenness of phi and the identity with Lambda are checked in task kernel.
# Method. Zeros: Newton / secant (mp.findroot) and the argument principle on a rectangle (phase increments of the ratio of consecutive values
# on a polygon of step h, integer count; the count is run at h and h/2 and both must agree). Real zero counts and orders: sign changes on a
# grid of step h refined by regula falsi, run at h and a second step. Critical rung and critical heat time: the double root (x, c) of
# G(x) = 0 = G'(x) for G(x) = Re Lambda(c + i x), resp. of the flow integral (Newton in two variables, seeded by a tracked nonreal pair
# whose Im^2 is extrapolated to 0). Heat integral: Gauss Legendre (48 points per panel of width 0.05 on [0, 3], phi even) at 50 digits.
# Error control: every task repeats its key numbers at a different precision (dps 20 and 30 for the Hurwitz evaluations, panel width 0.05
# and 0.04 and 50 and 60 digits for the heat integral) and logs the largest difference; step refinement for counts.
# Options: --task {defn,zeros,ladder,order,counts,thr,kernel,heat,all} --quick --out FILE (base name of the json; the log is next to it).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, time, argparse, os
from multiprocessing import Pool
import mpmath as mp

SQ5 = None
OFF = [(0.8085171825, 85.6993484854), (0.6508, 114.163), (0.5744, None), (0.7243, 176.702)]  # chapter values (seeds only)


def kappa():
    return (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)


def Lam(s, kap=None):
    if kap is None:
        kap = kappa()
    s = mp.mpmathify(s)
    h = [mp.zeta(s, mp.mpf(a) / 5) for a in (1, 2, 3, 4)]
    f = mp.mpf(5) ** (-s) * ((h[0] - h[3]) + kap * (h[1] - h[2]))
    return (5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * f


def Fb(c, z):
    z = mp.mpmathify(z)
    return (Lam(c + 1j * z) + Lam(c - 1j * z)) / 2


def ReLam(c, t):
    return mp.re(Lam(mp.mpc(c, t)))


def worker(job):
    kind, dps, args = job
    mp.mp.dps = dps
    if kind == 'Lam':
        return complex(Lam(mp.mpc(*args)))
    if kind == 'LamH':  # high precision value as string pair, args (re, im)
        v = Lam(mp.mpc(mp.mpf(args[0]), mp.mpf(args[1])))
        return (mp.nstr(v.real, 30), mp.nstr(v.imag, 30))
    if kind == 'Fb':
        c, x, y = args
        v = Fb(mp.mpf(c), mp.mpc(x, y))
        return (float(v.real), float(v.imag), mp.nstr(v.real, 25), mp.nstr(v.imag, 25))
    raise ValueError(kind)


def ncores():
    return min(8, os.cpu_count() or 1)


def pmap(jobs, chunk=4):
    with Pool(ncores()) as p:
        return p.map(worker, jobs, chunksize=chunk)


class Log:
    def __init__(self, path):
        self.f = open(path, 'a')

    def __call__(self, *a):
        s = ' '.join(str(x) for x in a)
        print(s)
        self.f.write(s + '\n')
        self.f.flush()


def polygon(x0, x1, y0, y1, h):
    pts = []
    def seg(a, b):
        n = max(1, int(abs(b - a) / h + 0.999))
        return [a + (b - a) * k / n for k in range(n)]
    for x in seg(complex(x0, y0), complex(x1, y0)): pts.append(x)
    for x in seg(complex(x1, y0), complex(x1, y1)): pts.append(x)
    for x in seg(complex(x1, y1), complex(x0, y1)): pts.append(x)
    for x in seg(complex(x0, y1), complex(x0, y0)): pts.append(x)
    pts.append(pts[0])
    return pts


def count_box(kind, c, x0, x1, y0, y1, h, dps=20):
    """argument principle count of zeros of Lambda (kind 'L', c ignored) or F_b (kind 'F') in the rectangle."""
    pts = polygon(x0, x1, y0, y1, h)
    if kind == 'L':
        jobs = [('LamH', dps, (str(p.real), str(p.imag))) for p in pts]
    else:
        jobs = [('Fb', dps, (c, p.real, p.imag)) for p in pts]
    res = pmap(jobs)
    if kind == 'L':
        vals = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in res]
    else:
        vals = [mp.mpc(mp.mpf(r[2]), mp.mpf(r[3])) for r in res]
    tot = mp.mpf(0)
    worst = 0
    for a, b in zip(vals[:-1], vals[1:]):
        d = mp.arg(b / a)
        tot += d
        worst = max(worst, abs(float(d)))
    return float(tot / (2 * mp.pi)), worst


def sign_scan(c, a, b, h, part='re', dps=15):
    ts = [a + k * h for k in range(int((b - a) / h) + 1)]
    res = pmap([('Lam', dps, (c, t)) for t in ts], chunk=8)
    v = [(z.real if part == 're' else z.imag) for z in res]
    return ts, v


def refine(c, t0, t1, part, dps=20):
    mp.mp.dps = dps
    g = (lambda t: mp.re(Lam(mp.mpc(c, t)))) if part == 're' else (lambda t: mp.im(Lam(mp.mpc(c, t))))
    return float(mp.findroot(g, (t0, t1), solver='illinois', tol=1e-16, maxsteps=40))


def zeros_from_scan(c, ts, v, part):
    out = []
    for i in range(len(ts) - 1):
        if v[i] == 0 or v[i] * v[i + 1] < 0:
            out.append((ts[i], ts[i + 1]))
    return out


def task_defn(log, R, quick):
    mp.mp.dps = 30
    kap = kappa()
    log('kappa =', mp.nstr(kap, 20))
    out = {'kappa': mp.nstr(kap, 25)}
    pts = [mp.mpc('0.3', 2), mp.mpc('0.7', 30), mp.mpc('0.9', 85.7), mp.mpc('0.2', 150)]
    err = max(abs(Lam(s) - Lam(1 - s)) / abs(Lam(s)) for s in pts)
    log('max rel |Lambda(s) - Lambda(1-s)| at 4 points, 30 digits:', mp.nstr(err, 3))
    out['fe_rel_err'] = mp.nstr(err, 3)
    real = max(abs(mp.im(Lam(mp.mpc('0.5', t)))) / (abs(mp.re(Lam(mp.mpc('0.5', t)))) + mp.mpf(10) ** -40) for t in [1, 5, 20, 60, 100, 150])
    log('max |Im Lambda|/|Re Lambda| on the critical line (6 points):', mp.nstr(real, 3))
    out['line_imag_ratio'] = mp.nstr(real, 3)
    # pin kappa: R(kappa) = Lambda_kappa(s) - Lambda_kappa(1-s) is affine in kappa
    sol = []
    for s in pts[:3]:
        r0 = Lam(s, mp.mpf(0)) - Lam(1 - s, mp.mpf(0))
        r1 = Lam(s, mp.mpf(1)) - Lam(1 - s, mp.mpf(1))
        sol.append(-r0 / (r1 - r0))
    log('kappa solved from the functional equation at 3 points:', [mp.nstr(z, 18) for z in sol])
    out['kappa_solved'] = [mp.nstr(z, 18) for z in sol]
    wrong = max(abs(Lam(s, kap + mp.mpf('0.1')) - Lam(1 - s, kap + mp.mpf('0.1'))) / abs(Lam(s)) for s in pts)
    log('with kappa + 0.1 the relative functional equation error is', mp.nstr(wrong, 3))
    out['wrong_kappa_err'] = mp.nstr(wrong, 3)
    # second precision
    mp.mp.dps = 20
    d = max(abs(Lam(s) - Lam(1 - s)) / abs(Lam(s)) for s in pts)
    log('same at 20 digits:', mp.nstr(d, 3))
    R['defn'] = out


def task_zeros(log, R, quick):
    out = []
    for beta0, g0 in OFF:
        if g0 is None:
            continue
        zs = []
        for dps in (20, 30):
            mp.mp.dps = dps
            z = mp.findroot(lambda s: Lam(s), mp.mpc(beta0, g0), tol=mp.mpf(10) ** (-dps + 6))
            zs.append(z)
        log('zero of Lambda near', beta0, g0, '->', mp.nstr(zs[1], 14), ' dps20-dps30 diff', mp.nstr(abs(zs[0] - zs[1]), 3))
        out.append({'seed': [beta0, g0], 'zero': [mp.nstr(zs[1].real, 14), mp.nstr(zs[1].imag, 14)], 'diff': mp.nstr(abs(zs[0] - zs[1]), 3)})
    R['zeros_newton'] = out
    # argument principle in a box around the first zero and (full) in strips for the complete list
    n1, w1 = count_box('L', 0, 0.55, 1.2, 84.0, 87.5, 0.05)
    n2, w2 = count_box('L', 0, 0.55, 1.2, 84.0, 87.5, 0.025)
    log('zeros of Lambda in [0.55,1.2] x [84,87.5] (h=0.05, 0.025):', n1, n2, ' max phase step', w1, w2)
    R['zeros_box_first'] = [n1, n2]
    if not quick:
        pass


def task_ladder(log, R, quick):
    out = {}
    cs = [0.6, 0.7, 0.78, 0.79, 0.795] if not os.environ.get('MI_SKIPPAIR') else []
    z0 = mp.mpc(85.69, 0.5)
    rows = []
    for c in cs:
        for dps in (20, 30):
            mp.mp.dps = dps
            # start near gamma + i(beta - c)
            z = mp.findroot(lambda z: Fb(mp.mpf(c), z), mp.mpc(85.69, 0.8085171825 - c), tol=mp.mpf(10) ** (-dps + 8))
            if dps == 20:
                z20 = z
        log('c =', c, 'pair root', mp.nstr(z, 14), 'Im^2', mp.nstr(z.imag ** 2, 8), 'dps diff', mp.nstr(abs(z - z20), 3))
        rows.append({'c': c, 'z': [mp.nstr(z.real, 14), mp.nstr(z.imag, 14)], 'im2': mp.nstr(z.imag ** 2, 8)})
    out['pair'] = rows
    if not cs:
        cs = []
    # linear fit of Im^2 in c through the first and last, residuals at the middle
    if rows:
      c_ = [r['c'] for r in rows]
    else:
      c_ = []
    if rows:
        y_ = [float(r['im2']) for r in rows]
        slope = (y_[-1] - y_[0]) / (c_[-1] - c_[0])
        res = [y_[i] - (y_[0] + slope * (c_[i] - c_[0])) for i in range(len(c_))]
        log('Im^2 slope between c=0.6 and 0.795:', slope, ' residuals from the chord:', res, ' chord zero at c =', c_[0] - y_[0] / slope)
        out['im2_chord'] = {'slope': slope, 'resid': res, 'c_zero': c_[0] - y_[0] / slope}
        # collision: double root
        mp.mp.dps = 20
        x0 = float(mp.re(mp.findroot(lambda z: Fb(mp.mpf(0.795), z), mp.mpc(85.69, 0.2))))
        def G(x, c): return mp.re(Lam(mp.mpc(c, x)))
        def H(x, c): return mp.diff(lambda t: G(t, c), x, h=mp.mpf('1e-6'))
        S = 1 / abs(Lam(mp.mpc('0.798', x0)))   # scale: Lambda is about 1e-38 here, so unscaled residuals look converged at once
        sol = mp.findroot(lambda x, c: (S * G(x, c), S * H(x, c)), (mp.mpf(x0), mp.mpf('0.798')), tol=1e-24, maxsteps=40)
        log('double root (x, c) of Re Lambda(c + i x):', mp.nstr(sol[0], 12), mp.nstr(sol[1], 12))
        out['collision'] = [mp.nstr(sol[0], 12), mp.nstr(sol[1], 12)]
    # box counts
    box = []
    cl = [0.6, 0.7, 0.78, 0.79, 0.795, 0.80, 0.85] if not quick else [0.7, 0.8]
    if os.environ.get('MI_CS'):
        cl = [float(x) for x in os.environ['MI_CS'].split(',')]
    for c in cl:
        n1, _ = count_box('F', c, 83, 88.5, -0.9, 0.9, 0.05)
        n2, _ = count_box('F', c, 83, 88.5, -0.9, 0.9, 0.03)
        ts, v = sign_scan(c, 83, 88.5, 0.01)
        nreal = len(zeros_from_scan(c, ts, v, 're'))
        log('c =', c, 'zeros of F_b in box (h=0.05, 0.03):', n1, n2, ' real zeros of Re Lambda on (83,88.5):', nreal)
        box.append({'c': c, 'count_h05': n1, 'count_h03': n2, 'real_zeros': nreal})
    out['box'] = box
    R['ladder'] = out


def task_order(log, R, quick):
    out = []
    for c in ([0.800, 0.805, 0.82, 0.85] if not quick else [0.8, 0.85]):
        seqs = []
        for h in (0.01, 0.02):
            ts, va = sign_scan(c, 83, 88.5, h, 're')
            _, vb = sign_scan(c, 83, 88.5, h, 'im')
            za = [refine(c, a, b, 're') for a, b in zeros_from_scan(c, ts, va, 're')]
            zb = [refine(c, a, b, 'im') for a, b in zeros_from_scan(c, ts, vb, 'im')]
            ev = sorted([(z, 'A') for z in za] + [(z, 'B') for z in zb])
            seqs.append((''.join(e[1] for e in ev), za, zb))
        log('c =', c, 'order (h=0.01, 0.02):', seqs[0][0], seqs[1][0], ' A zeros', [round(z, 4) for z in seqs[0][1]], ' B zeros', [round(z, 4) for z in seqs[0][2]])
        out.append({'c': c, 'order': seqs[0][0], 'order_h02': seqs[1][0], 'A': seqs[0][1], 'B': seqs[0][2]})
    R['order'] = out


def task_counts(log, R, quick):
    T = 40 if quick else 200
    out = []
    cl = [0.95, 0.85, 0.7, 0.6]
    if os.environ.get('MI_CS'):
        cl = [float(x) for x in os.environ['MI_CS'].split(',')]
    for c in cl:
        n = []
        for h in (0.1, 0.07):
            ts, v = sign_scan(c, 0.5, T, h, 're')
            n.append(len(zeros_from_scan(c, ts, v, 're')))
        log('c =', c, 'real zeros of Re Lambda(c + it) on (0.5,', T, ') h=0.1, 0.07:', n)
        out.append({'c': c, 'T': T, 'counts': n})
    R['counts'] = out


# ---- heat flow
def kernel_nodes(dps, h=0.05, U=3.0, deg=5):
    mp.mp.dps = dps
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = gl.calc_nodes(deg, mp.mp.prec)
    kap = kappa()
    a = (1 - 1j * kap) / 2
    chi = {1: 1, 2: 1j, 3: -1j, 4: -1}
    cn = {r: 2 * mp.re(a * chi[r]) for r in (1, 2, 3, 4)}  # a chi + abar chibar; conj of residue class: chi(n) depends on n mod 5
    cn[0] = mp.mpf(0)
    def phi(u):
        y = mp.exp(2 * u)
        nmax = int(mp.sqrt(5 * (dps * 2.4 + 20) / (mp.pi * y))) + 2
        s = mp.mpf(0)
        for n in range(1, nmax + 1):
            s += n * cn[n % 5] * mp.exp(-mp.pi * n * n * y / 5)
        return 2 * mp.sqrt(mp.pi / 5) * mp.exp(1.5 * u) * s
    us, ws, ph = [], [], []
    npan = int(round(U / h))
    h = mp.mpf(int(round(h * 1000))) / 1000
    for k in range(npan):
        lo = k * h
        for (x, w) in nodes:
            u = lo + (x + 1) * h / 2
            us.append(u)
            ws.append(w * h / 2)
            ph.append(phi(u))
    return us, ws, ph, phi


def heat_val(K, lam, z, d=0):
    us, ws, ph = K[:3]
    lam = mp.mpf(lam)
    s = mp.mpf(0)
    for u, w, p in zip(us, ws, ph):
        e = w * p * mp.exp(lam * u * u)
        if d == 0:
            s += e * mp.cos(z * u)
        elif d == 1:
            s += -e * u * mp.sin(z * u)
        elif d == 2:
            s += -e * u * u * mp.cos(z * u)
    return 2 * s


def heat_dlam(K, lam, z):
    us, ws, ph = K[:3]
    lam = mp.mpf(lam)
    s = mp.mpf(0)
    for u, w, p in zip(us, ws, ph):
        s += w * p * u * u * mp.exp(lam * u * u) * mp.cos(z * u)
    return 2 * s


def heat_dxlam(K, lam, z):
    us, ws, ph = K[:3]
    lam = mp.mpf(lam)
    s = mp.mpf(0)
    for u, w, p in zip(us, ws, ph):
        s += -w * p * u ** 3 * mp.exp(lam * u * u) * mp.sin(z * u)
    return 2 * s


def task_kernel(log, R, quick):
    K = kernel_nodes(90)
    phi = K[3]
    out = {}
    mp.mp.dps = 90
    errs = []
    for x in (3, 20, 50, 85, 177):
        v = heat_val(K, 0, x)
        L = mp.re(Lam(mp.mpc('0.5', x)))
        errs.append(mp.nstr(abs(v - L) / abs(L), 3))
        log('x =', x, 'int phi cos vs Lambda(1/2 + i x):', mp.nstr(v, 12), mp.nstr(L, 12), 'rel', errs[-1])
    out['rel_err'] = errs
    out['even'] = [mp.nstr(abs(phi(u) - phi(-u)) / phi(u), 3) for u in (mp.mpf('0.1'), mp.mpf('0.4'), mp.mpf('0.8'))]
    log('phi(u) vs phi(-u):', out['even'])
    R['kernel'] = out


def seed_pair(K, zeta, lam_hi):
    pass


def heat_critical(job):
    """Double root (x, lam) of the heat flow for one zero of Lambda. job = (beta, gamma, dps, h)."""
    beta, gam, dps, h = job
    K = kernel_nodes(dps, h=h)
    mp.mp.dps = dps
    eta = beta - 0.5
    # locate the pair above the zero: root z = x + i y of the flow integral, tracked in lam
    F = lambda lam, z: heat_val(K, lam, z)
    lam_t = [0.55 * eta * eta / 2, 0.70 * eta * eta / 2, 0.80 * eta * eta / 2]
    zs = []
    zprev = mp.mpc(gam, eta * 0.8)
    for lam in lam_t:
        zprev = mp.findroot(lambda z: F(lam, z), zprev, tol=mp.mpf(10) ** (-dps // 2 + 4), maxsteps=40)
        zs.append(zprev)
    y2 = [z.imag ** 2 for z in zs]
    slope = (y2[2] - y2[0]) / (lam_t[2] - lam_t[0])
    lam0 = lam_t[2] - y2[2] / slope
    x0 = zs[2].real
    xs = mp.mpf(x0); ls = mp.mpf(lam0)
    for it in range(12):
        g1 = heat_val(K, ls, xs); g2 = heat_val(K, ls, xs, 1)
        a = heat_dlam(K, ls, xs)           # dg1/dlam
        b = heat_val(K, ls, xs, 2)         # dg2/dx
        c_ = heat_dxlam(K, ls, xs)         # dg2/dlam
        # dg1/dx = g2
        det = g2 * c_ - a * b
        dx = (-g1 * c_ + a * g2) / det
        dl = (-g2 * g2 + g1 * b) / det
        xs += dx; ls += dl
        if abs(dx) + abs(dl) < mp.mpf(10) ** (-dps // 2 + 3):
            break
    return (mp.nstr(xs, 14), mp.nstr(ls, 14), mp.nstr(heat_val(K, ls, xs), 3), mp.nstr(heat_val(K, ls, xs, 1), 3), [mp.nstr(z, 10) for z in zs], mp.nstr(lam0, 8))


def task_heat(log, R, quick):
    # env MI_HR = 1 (dps 80, panel 0.05) or 2 (dps 100, panel 0.04, the second run); MI_HZ = comma list of zero indices (default all)
    zs = R.get('thr_zeros')
    hr = int(os.environ.get('MI_HR', '1'))
    idx = [int(x) for x in os.environ.get('MI_HZ', '0,1,2').split(',')]
    par = (80, 0.05) if hr == 1 else (100, 0.04)
    jobs = [(zs[i][0], zs[i][1], par[0], par[1]) for i in idx]
    with Pool(len(jobs)) as p:
        r1 = p.map(heat_critical, jobs, chunksize=1)
    out = []
    for i, a in zip(idx, r1):
        z = zs[i]
        eta = z[0] - 0.5
        lc = mp.mpf(a[1])
        log('run', hr, 'dps', par[0], 'panel', par[1], 'zero', z[:2], 'heat double root x =', a[0], 'lambda_c =', a[1], ' residuals', a[2], a[3], ' lambda_c/(eta^2/2) =', mp.nstr(lc / (eta * eta / 2), 6), ' seed track', a[4], a[5])
        out.append({'zero': z[:2], 'run': hr, 'x': a[0], 'lam_c': a[1], 'ratio': mp.nstr(lc / (eta * eta / 2), 6), 'seed_track': a[4]})
    R['heat'] = out


def crit_rung(job):
    beta, gam, dps = job
    mp.mp.dps = dps
    eta = beta - 0.5
    # track the pair in c
    cs = [0.5 + 0.55 * eta, 0.5 + 0.7 * eta, 0.5 + 0.8 * eta]
    zprev = mp.mpc(gam, eta * 0.8)
    zs = []
    for c in cs:
        zprev = mp.findroot(lambda z: Fb(mp.mpf(c), z), zprev, tol=mp.mpf(10) ** (-dps + 8))
        zs.append(zprev)
    y2 = [z.imag ** 2 for z in zs]
    slope = (y2[2] - y2[0]) / (cs[2] - cs[0])
    c0 = cs[2] - y2[2] / slope
    def G(x, c): return mp.re(Lam(mp.mpc(c, x)))
    def H(x, c): return mp.diff(lambda t: G(t, c), x, h=mp.mpf('1e-6'))
    S = 1 / abs(Lam(mp.mpc(c0, zs[2].real)))
    sol = mp.findroot(lambda x, c: (S * G(x, c), S * H(x, c)), (zs[2].real, c0), tol=mp.mpf(10) ** (-20), maxsteps=40)
    return (mp.nstr(sol[0], 12), mp.nstr(sol[1] - mp.mpf('0.5'), 12), [mp.nstr(z, 10) for z in zs], mp.nstr(S * G(sol[0], sol[1]), 3), mp.nstr(S * H(sol[0], sol[1]), 3))


def task_thr(log, R, quick):
    zs = []
    for beta0, g0 in OFF:
        if g0 is None:
            continue
        mp.mp.dps = 25
        z = mp.findroot(lambda s: Lam(s), mp.mpc(beta0, g0), tol=1e-18)
        zs.append((float(z.real), float(z.imag), mp.nstr(z.real, 12), mp.nstr(z.imag, 12)))
    R['thr_zeros'] = zs
    jobs = [(z[0], z[1], 20) for z in zs]
    with Pool(3) as p:
        r1 = p.map(crit_rung, jobs, chunksize=1)
    jobs = [(z[0], z[1], 28) for z in zs]
    with Pool(3) as p:
        r2 = p.map(crit_rung, jobs, chunksize=1)
    rows = []
    for z, a, b in zip(zs, r1, r2):
        eta = z[0] - 0.5
        L = mp.log(z[1] / (2 * mp.pi))
        bc = mp.mpf(a[1])
        log('zero', z[2], z[3], 'eta =', mp.nstr(eta, 6), 'L =', mp.nstr(L, 5), 'b_c =', a[1], '(dps28:', b[1], ') residuals', a[3], a[4], ' b_c/eta =', mp.nstr(bc / eta, 6), 'x_c =', a[0])
        rows.append({'beta': z[2], 'gamma': z[3], 'eta': eta, 'L': float(L), 'b_c': a[1], 'b_c_dps28': b[1], 'bc_over_eta': float(bc / eta), 'x_c': a[0]})
    R['thr'] = rows
    # local model: cos(b d/dz) z^2 = z^2 - b^2, exp(-lam d^2/dz^2) z^2 = z^2 - 2 lam
    import sympy as sp
    z, b, lam = sp.symbols('z b lam')
    p = z ** 2 + sp.Rational(1, 10)
    cosb = sum((-1) ** k * b ** (2 * k) / sp.factorial(2 * k) * sp.diff(p, z, 2 * k) for k in range(3))
    heat = sum((-lam) ** k / sp.factorial(k) * sp.diff(p, z, 2 * k) for k in range(3))
    log('local model: cos(b d/dz)(z^2 + eta^2) =', sp.expand(cosb), '; exp(-lam d^2/dz^2)(z^2 + eta^2) =', sp.expand(heat))
    R['local_model'] = [str(sp.expand(cosb)), str(sp.expand(heat))]


def task_zerosall(log, R, quick):
    # all zeros of Lambda with 0.52 <= Re s <= 1.5, 0.1 <= Im s <= T (T = 200, quick 100): strips of height 10 by the argument
    # principle at step H (env MI_H, default 0.1; run it also at 0.05), strips with a zero are halved down to height 0.625, then Newton.
    T = 100.0 if quick else 200.0
    H = float(os.environ.get('MI_H', '0.1'))
    x0, x1 = 0.52, 1.5
    edges = [0.1] + [10.0 * k for k in range(1, int(T / 10) + 1)]
    if os.environ.get('MI_RANGE'):
        lo_, hi_ = [float(x) for x in os.environ['MI_RANGE'].split(',')]
        edges = [lo_ if lo_ > 0 else 0.1] + [10.0 * k for k in range(int(lo_ / 10) + 1, int(hi_ / 10) + 1)]
        T = hi_
    found, tot = [], 0
    def rec(y0, y1):
        n, w = count_box('L', 0, x0, x1, y0, y1, H, dps=15)
        ni = int(round(n))
        if ni == 0:
            return
        if y1 - y0 > 0.7:
            m = (y0 + y1) / 2
            rec(y0, m); rec(m, y1)
            return
        found.append((y0, y1, ni))
    for a, b in zip(edges[:-1], edges[1:]):
        n, w = count_box('L', 0, x0, x1, a, b, H, dps=15)
        log('strip', a, b, 'count', n, 'max phase step', round(w, 3))
        if int(round(n)) > 0:
            tot += int(round(n))
            rec(a, b)
    zs = []
    for y0, y1, ni in found:
        mp.mp.dps = 25
        for sx in (0.6, 0.8, 1.0):
            try:
                z = mp.findroot(lambda s: Lam(s), mp.mpc(sx, (y0 + y1) / 2), tol=1e-18, maxsteps=30)
            except Exception:
                continue
            if x0 <= z.real <= x1 and y0 - 1 <= z.imag <= y1 + 1 and not any(abs(z - w) < 1e-8 for w in zs):
                zs.append(z)
    log('total zeros in [0.52,1.5] x [0.1,%g] by the argument principle at step %g: %d' % (T, H, tot))
    for z in sorted(zs, key=lambda z: z.imag):
        log('  zero', mp.nstr(z, 14))
    R['zerosall'] = {'T': T, 'h': H, 'total': tot, 'zeros': [mp.nstr(z, 14) for z in zs]}


TASKS = ['defn', 'zeros', 'ladder', 'order', 'counts', 'thr', 'kernel', 'heat']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', default='all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    base = a.out or ('data/mi_dh' if not a.quick else 'misc/reviews/nick-miyatake/quick/G4/mi_dh')
    os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
    tasks = TASKS if a.task == 'all' else a.task.split(',')
    for t in tasks:
        log = Log(base + '.log')
        R = {}
        if t == 'heat':
            try:
                R['thr_zeros'] = json.load(open('data/mi_dh_thr.json'))['thr_zeros']
            except Exception:
                R['thr_zeros'] = [[0.8085171825, 85.6993484854], [0.6508, 114.163], [0.7243, 176.702]]
        t0 = time.time()
        log('== task', t, 'quick' if a.quick else 'full', time.strftime('%Y-%m-%d %H:%M'))
        globals()['task_' + t](log, R, a.quick)
        log('== task', t, 'done in %.1f s' % (time.time() - t0))
        json.dump(R, open(base + '_' + t + '.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
