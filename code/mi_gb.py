# mi_gb.py: the numerical claims of Chapter 22 (Nick Polson's Miyatake chapter) about g_b(z) = int Phi(u) sech(b u) e^(izu) du and the families
# built from the kernel: the shift identity, the Levy area and logistic formulas, the real zero counts of g_b and H_(-b^2/2), the thresholds in
# b, the transition at b = 2 (zeros of Miyatake's g, the cone, the escape law), the clocks PG(1), PG(4), PG(16) and the flow H_(-1/2), the strips,
# and the zeros of g_2 = 2 sum_{n in S(2)} n^4 K_(iw)(2 pi n^2). Group G3 of the computation briefs of 22 Oct 2026. Written from the mathematics
# (the collaborator's scripts were not read, run or copied). Uses the library of code/mi_miyatake.py (Phi, sweep, gl_panels, C_b).
# Conventions. Phi(u) = sum (4 pi^2 n^4 e^(9u/2) - 6 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)), Xi(z) = int Phi(u) e^(izu) du. A clock with weight s(u) means
# F(z) = int Phi(u) s(u) e^(izu) du: g_b has s = sech(bu); H_lambda has s = e^(lambda u^2); the clock omega_h = PG(h,0)/h (mean 1/4, E e^(-s omega_h) =
# cosh(sqrt(s/(2h)))^(-h)) at radius b gives s(u) = sech(b u/sqrt(h))^h (h = 1: g_b, h -> infinity: H_(-b^2/2)). This definition of PG(4), PG(16) is
# ours: the chapter does not define it.
# Method. F is evaluated by the trapezoid rule of step 0.01 in u (|u| <= 2.4, where Phi < 1e-165): even integrand, so the error is the Fourier transform at
# 2 pi/h - z; on a grid of constant step by the fixed point recurrence (400 bits) of mi_miyatake.sweep. Real zeros: sign changes of F(x) on the grid of step
# 0.05 (a second run at 330 bits and step 0.04). Nonreal zeros: argument principle on [0,X] x [-Y,Y] using F(conj z) = conj F(z) (top row, verticals), minus
# the real sign changes; their location by the winding number of the cells of a grid (rows y = 0.125 + 0.25 m) and Newton's method in mpmath.
# Output: data/mi_gb_<task>.json and .log (quick runs: --out FILE).
import sys, os, json, time, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mpm
from mpmath import mp, mpf, mpc
import mi_miyatake as M

LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def weight_fun(kind, par):
    if kind == 'sech':  # g_b
        b = mpf(par)
        return lambda u: mpm.sech(b * u)
    if kind == 'pg':   # clock PG(h)/h at radius b: par = (b, h)
        b, h = mpf(par[0]), mpf(par[1])
        return lambda u: mpm.sech(b * u / mpm.sqrt(h)) ** h
    if kind == 'gauss':  # H_lambda, lambda = -b^2/2
        b = mpf(par)
        return lambda u: mpm.exp(-b * b * u * u / 2)
    if kind == 'one':
        return lambda u: mpf(1)
    raise ValueError(kind)


def kernel(kind, par, h=mpf(1) / 100, umax=2.4, dps=130):
    mp.dps = dps
    f = weight_fun(kind, par)
    N = int(math.ceil(umax / float(h)))
    nodes = [h * j for j in range(N + 1)]
    w = [(1 if j == 0 else 2) * h * M.Phi(nodes[j]) * f(nodes[j]) for j in range(N + 1)]
    return nodes, w


def mirror(nodes, w):
    if nodes[0] != 0:
        return nodes + [-u for u in nodes], [x / 2 for x in w] * 2
    mn = nodes + [-u for u in nodes[1:]]
    mw = [w[0]] + [x / 2 for x in w[1:]] * 2
    return mn, mw


def direct(nodes, w, z):
    """F(z) for complex z by a direct sum (even weights)."""
    return mpm.fsum(wj * mpm.cos(z * u) for u, wj in zip(nodes, w))


def real_count(nodes, w, zmax, dz, P):
    K = int(zmax / dz) + 1
    v = M.sweep(nodes, [w], mpf(0), mpf(dz), K, P)[0]
    r = [x.real for x in v]
    zr = [dz * i + dz * r[i] / (r[i] - r[i + 1]) for i in range(len(r) - 1) if (r[i] > 0) != (r[i + 1] > 0)]
    return zr


def nonreal_pairs(nodes, w, X, Y, P, dxt=0.1, dyv=0.1):
    """(total zeros in [0,X] x [-Y,Y]) by the argument principle; uses the conjugation symmetry."""
    mn, mw = mirror(nodes, w)
    K = int(round(X / dxt)) + 1
    top = M.sweep(mn, [mw], mpc(0, Y), mpf(dxt), K, P)[0]
    dtop = M.unwrap_change(top)
    nv = int(round(Y / dyv)) + 1
    vX = M.sweep(mn, [mw], mpf(X), mpc(0, dyv), nv, P)[0]
    v0 = M.sweep(mn, [mw], mpf(0), mpc(0, dyv), nv, P)[0]
    W = (-2 * dtop + 2 * M.unwrap_change(vX) - 2 * M.unwrap_change(v0)) / (2 * math.pi)
    return W


def cell_zeros(nodes, w, x0, x1, dx, ys, P):
    mn, mw = mirror(nodes, w)
    K = int(round((x1 - x0) / dx)) + 1
    rows = []
    for y in ys:
        v = M.sweep(mn, [mw], mpc(x0, y), mpf(dx), K, P)[0]
        rows.append(np.array(v))
    ang = lambda a, b: np.angle(b / a)
    zs = []
    for m in range(len(ys) - 1):
        A, B = rows[m], rows[m + 1]
        wnd = (ang(A[:-1], A[1:]) + ang(A[1:], B[1:]) + ang(B[1:], B[:-1]) + ang(B[:-1], A[:-1])) / (2 * np.pi)
        for i in np.nonzero(np.abs(wnd) > 0.5)[0]:
            zs.append((x0 + dx * (i + 0.5), 0.5 * (ys[m] + ys[m + 1]), int(round(wnd[i]))))
    return zs


def refine(nodes, w, z0, dps=120, its=14):
    mp.dps = dps
    mn, mw = mirror(nodes, w)
    z = mpc(z0)
    for _ in range(its):
        f = mpc(0); d = mpc(0)
        for u, wj in zip(mn, mw):
            e = wj * mpm.exp(1j * z * u)
            f += e; d += 1j * u * e
        step = f / d
        z -= step
        if abs(step) < mpf(10) ** (-40):
            break
    return z


def task_ident(R, part='abc'):
    out = {}
    M.setdps(40)
    h = mpf(1) / 100
    # shift identity
    for b in ((mpf('0.7'), mpf('1.3'), mpf(9) / 2) if 'a' in part else ()):
        nodes, w = kernel('sech', b, h, dps=60)
        for z in (mpf(3), mpf('14.134725141734693790457251983562')):
            lhs = M.Xi(z)
            rhs = (direct(nodes, w, z + 1j * b) + direct(nodes, w, z - 1j * b)) / 2
            out['shift identity Xi(z) vs (g_b(z+ib)+g_b(z-ib))/2, b=%s, z=%s: ratio' % (mpm.nstr(b, 3), mpm.nstr(z, 5))] = mpm.nstr(mpm.re(rhs) / lhs if lhs != 0 else rhs, 14) + ' (imag %s)' % mpm.nstr(mpm.im(rhs), 3)
    # Levy area convolution
    for b in ((mpf('0.5'), mpf('1.5')) if 'b' in part else ()):
        nodes, w = kernel('sech', b, h, dps=60)
        for x in (mpf(3), mpf(10)):
            g = direct(nodes, w, x)
            vals = []
            for step in (mpf(1) / 10, mpf(1) / 14):
                n = int(45 / step)
                I = mpf(0)
                for j in range(-n, n + 1):
                    y = j * step
                    I += M.Xi(x - y) * mpm.sech(mp.pi * y / (2 * b)) / (2 * b)
                vals.append(I * step)
            out['Levy area: g_b(x) vs int Xi(x-y) sech(pi y/2b)/(2b) dy, b=%s x=%s' % (mpm.nstr(b, 3), mpm.nstr(x, 3))] = [mpm.nstr(mpm.re(g), 14), mpm.nstr(vals[0], 14), mpm.nstr(vals[1], 14)]
    # Thorin product and PG Laplace transform
    b = mpf('0.9'); u = mpf('0.8')
    if 'a' not in part and 'c' not in part:
        R['ident_' + part] = out
        for k, v in out.items():
            log('ident', k, v)
        return
    prod = mpm.nprod(lambda k: 1 / (1 + 4 * b * b * u * u / ((2 * k - 1) * mp.pi) ** 2), [1, mpm.inf])
    out['Thorin product vs sech(bu), b=0.9 u=0.8'] = [mpm.nstr(prod, 14), mpm.nstr(mpm.sech(b * u), 14)]
    s = 2 * b * b * u * u
    out['sech(sqrt(s/2)) with s = 2 b^2 u^2 vs sech(bu)'] = [mpm.nstr(mpm.sech(mpm.sqrt(s / 2)), 14), mpm.nstr(mpm.sech(b * u), 14)]
    if 'c' not in part:
        R['ident_' + part] = out
        for k, v in out.items():
            log('ident', k, v)
        return
    # logistic partition: E[Y^s/(1+Y^c)] + E[Y^(s+c)/(1+Y^c)] with log Y of density 2 Phi(u) e^(-u/2), s = 1/2 + iz, c = 2b
    b = mpf('0.9'); c = 2 * b; z = mpf(7)
    st = mpf(1) / 10
    N = int(2.4 / st)
    t1 = mpc(0); t2 = mpc(0)
    for j in range(-N, N + 1):
        uu = j * st
        dens = 2 * M.Phi(uu) * mpm.exp(-uu / 2)
        Ys = mpm.exp(uu * (mpf(1) / 2 + 1j * z))
        Yc = mpm.exp(c * uu)
        t1 += dens * Ys / (1 + Yc) * st
        t2 += dens * Ys * Yc / (1 + Yc) * st
    nodes, w = kernel('sech', b, h, dps=60)
    out['logistic: E[Y^s/(1+Y^c)] vs g_b(z+ib); E[Y^(s+c)/(1+Y^c)] vs g_b(z-ib); sum vs 2 xi(s) (step 0.1)'] = [
        mpm.nstr(t1, 12), mpm.nstr(direct(nodes, w, z + 1j * b), 12), mpm.nstr(t2, 12), mpm.nstr(direct(nodes, w, z - 1j * b), 12), mpm.nstr(t1 + t2, 12), mpm.nstr(2 * M.Xi(z), 12)]
    R['ident'] = out
    for k, v in out.items():
        log('ident', k, v)


def task_counts(R, quick, blist=('0.5', '0.75', '1', '1.5')):
    zmax = 100 if quick else 250
    out = {}
    M.setdps(60)
    nodes, w = kernel('one', 0)
    zr = real_count(nodes, w, zmax, 0.05, 400)
    out['Xi real zeros (0,%d) (kernel with s = 1)' % zmax] = len(zr)
    zs = [mpf(l.strip()) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_1700.txt')) if l.strip()]
    out['zeta zeros with ordinate < %d (file)' % zmax] = sum(1 for z in zs if z < zmax)
    out['max |first 3 zeros of Xi from kernel - file|'] = max(abs(a - float(b)) for a, b in zip(zr[:3], zs[:3]))
    for b in blist:
        for kind in ('sech', 'gauss'):
            nodes, w = kernel(kind, mpf(b))
            c1 = len(real_count(nodes, w, zmax, 0.05, 400))
            c2 = len(real_count(nodes, w, zmax, 0.04, 330))
            out['%s b=%s real zeros in (0,%d), P400 dz0.05 | P330 dz0.04' % ('g_b' if kind == 'sech' else 'H_(-b^2/2)', b, zmax)] = [c1, c2]
            log('counts', kind, b, c1, c2)
    R['counts'] = out
    for k, v in out.items():
        log('counts', k, v)


def task_thresh(R, blist, quick):
    """for each b the number of nonreal zeros of g_b with real part in (0,T) for T = 40, 100, 160 (Y = 3)."""
    out = {}
    Y = 3.0
    for b in blist:
        M.setdps(60)
        nodes, w = kernel('sech', mpf(b))
        zr = real_count(nodes, w, 160, 0.05, 400)
        row = {}
        for T in (40, 100, 160):
            c = sum(1 for z in zr if z < T)
            W = nonreal_pairs(nodes, w, T, Y, 400)
            row[T] = [c, round(W, 3), int(round(W)) - c]
        out[str(b)] = row
        log('thresh b=%s' % b, row, ' (real count, winding, nonreal zeros counted with their conjugates)')
    R['thresh_%s_%s' % (blist[0], blist[-1])] = out


def task_nonreal(R, case, zmax, quick, nrows=13):
    kinds = {'b0.5': ('sech', 0.5), 'b1': ('sech', 1), 'b1.5': ('sech', 1.5), 'pg4': ('pg', (1, 4)), 'pg16': ('pg', (1, 16)), 'h': ('gauss', 1), 'pg1': ('pg', (1, 1))}
    kind, par = kinds[case]
    M.setdps(60)
    nodes, w = kernel(kind, par)
    Y = 3.0
    t0 = time.time()
    zr = real_count(nodes, w, zmax, 0.05, 400)
    W = nonreal_pairs(nodes, w, zmax, Y, 400)
    W2 = nonreal_pairs(nodes, w, zmax, Y, 330, 0.08, 0.1)
    c = len(zr)
    out = {'case': case, 'zmax': zmax, 'real zeros in (0,zmax)': c, 'winding P400': round(W, 3), 'winding P330 step 0.08': round(W2, 3),
           'nonreal pairs (z and conj z counted once)': (int(round(W)) - c) / 2}
    log('nonreal', out)
    ys = [0.125 + 0.25 * m for m in range(nrows)]
    cells = cell_zeros(nodes, w, 0.0, zmax, 0.2, ys, int(os.environ.get('MI_CELLP', '400')))
    out['cells with a zero (Im > 0.125)'] = len(cells)
    if cells:
        top = max(cells, key=lambda t: t[1])
        zz = refine(nodes, w, mpc(top[0], top[1]))
        out['highest cell'] = top
        out['refined highest zero'] = mpm.nstr(zz, 10)
        out['max Im over cell centres'] = max(t[1] for t in cells)
        # refine all zeros with Im above 0.9 * highest cell
        hz = [t for t in cells if t[1] >= float(os.environ.get('MI_YREF', top[1] - 0.2))][:int(os.environ.get('MI_NREF', '3'))]
        ref = [mpm.nstr(refine(nodes, w, mpc(t[0], t[1])), 8) for t in hz]
        out['refined zeros near the top'] = ref
        out['max Im of the refined zeros'] = max(float(mpm.im(refine(nodes, w, mpc(t[0], t[1])))) for t in hz) if False else None
    out['seconds'] = round(time.time() - t0, 1)
    R['nonreal_' + case] = out
    for k, v in out.items():
        log('nonreal', case, k, v)


def task_g2(R, quick):
    """zeros of g_2(w) = 2 sum_{n in S(2)} n^4 K_(iw)(2 pi n^2), K_(iw)(a) = int_0^inf e^(-a cosh t) cos(wt) dt, in 0 < Re w < 60."""
    M.setdps(145)
    ns = [1, 2, 4, 8, 16]
    sn, sw = M.gl_panels(4.7, 0.2, 24)
    nodes = list(sn)
    wts = [2 * 2 * sum(n ** 4 * mpm.exp(-2 * mp.pi * n * n * mpm.cosh(t)) for n in ns) * wj for t, wj in zip(sn, sw)]  # g_2 = 2 sum n^4 int_0^inf = (int_R ... ) so weights for sum w_j cos(w t_j): 2 * n^4 * wj
    wts = [x / 2 for x in wts]
    X = 60
    zr = real_count(nodes, wts, X, 0.05, 400)
    Wt = nonreal_pairs(nodes, wts, X, 4.0, 400)
    out = {'real zeros in (0,60)': len(zr), 'winding on [0,60]x[-4,4]': round(Wt, 3), 'nonreal in right half plane (pairs)': (round(Wt) - len(zr)) / 2}
    ys = [0.1 + 0.2 * m for m in range(13)]
    cells = cell_zeros(nodes, wts, 0.0, X, 0.2, ys, 340)
    out['cells'] = len(cells)
    zs = []
    for (x, y, k) in cells:
        z = refine(nodes, wts, mpc(x, y), dps=120)
        zs.append((float(mpm.re(z)), float(mpm.im(z))))
    zs = sorted(set((round(a, 6), round(b, 6)) for a, b in zs))
    out['refined zeros (Re, Im), Im > 0'] = [z for z in zs if z[1] > 0]
    out['max Im'] = max([z[1] for z in zs] or [0])
    R['g2'] = out
    for k, v in out.items():
        log('g2', k, v)


def escape(b, x, Cb):
    L = mpm.log(x / (2 * mp.pi))
    return b - mpf(1) / 2 + 2 / L * (mp.pi * (b - 2) * x / (4 * b) - mpf(3) / 2 * mpm.log(x) + mpm.log(Cb) - mpm.log(4 * mp.pi) / 2)


def task_g92(R, ylo, yhi, quick, xmax=160):
    """zeros of g_{9/2} (variable z = 2x of g) with real part below xmax: cell winding in y in [ylo,yhi]; the escape law."""
    M.setdps(60)
    nodes, w = kernel('sech', mpf(9) / 2, dps=130)
    out = {}
    ys = [ylo + 0.5 * m for m in range(int(round((yhi - ylo) / 0.5)) + 1)]
    t0 = time.time()
    cells = cell_zeros(nodes, w, 0.0, xmax, 0.4, ys, 400)
    out['cells with a zero, y in [%s,%s]' % (ylo, yhi)] = len(cells)
    zs = []
    for (x, y, k) in cells:
        z = refine(nodes, w, mpc(x, y), dps=120)
        zs.append((float(mpm.re(z)), float(mpm.im(z))))
    zs = sorted(set((round(a, 4), round(b, 4)) for a, b in zs))
    out['zeros (Re,Im)'] = zs
    M.setdps(40)
    Cb = M.C_b(mpf(9) / 2).real
    dev = [(zr, zi, float(escape(mpf(9) / 2, mpf(zr), Cb)) - zi) for zr, zi in zs if zr > 8]
    out['C_9/2'] = float(Cb)
    out['escape law minus zero Im (Re, Im, diff)'] = [(a, b, round(c, 3)) for a, b, c in dev]
    out['max |diff|'] = max([abs(c) for a, b, c in dev] or [0])
    out['seconds'] = round(time.time() - t0, 1)
    R['g92_%s_%s' % (ylo, yhi)] = out
    for k, v in out.items():
        log('g92', k, v)


def task_g92rect(R):
    """Miyatake's g(x) = g_{9/2}(2x)/(4 pi^2): no zeros in (0,26) x (-4,4) i.e. z = 2x in (0,52) x (-8,8); none on (0,60] i.e. z in (0,120]."""
    M.setdps(60)
    nodes, w = kernel('sech', mpf(9) / 2, dps=130)
    out = {}
    zr = real_count(nodes, w, 120, 0.05, 400)
    zr2 = real_count(nodes, w, 120, 0.04, 330)
    out['real zeros of g_{9/2} in (0,120] (x in (0,60]), step 0.05 P400 | step 0.04 P330'] = [len(zr), len(zr2)]
    W = nonreal_pairs(nodes, w, 52, 8.0, 400)
    W2 = nonreal_pairs(nodes, w, 52, 8.0, 330, 0.08, 0.1)
    out['winding on [0,52]x[-8,8] (zeros of g in (0,26)x(-4,4)), P400 | P330'] = [round(W, 4), round(W2, 4)]
    nz = mpc(17.9168, 10.3988)
    z1 = refine(nodes, w, nz)
    out['zero near 17.9168+10.3988i'] = mpm.nstr(z1, 12)
    out['in the scale of g (z/2)'] = mpm.nstr(z1 / 2, 8)
    out['g(0) = g_{9/2}(0)/(4 pi^2)'] = mpm.nstr(direct(nodes, w, mpf(0)) / (4 * mp.pi ** 2), 10)
    R['g92rect'] = out
    for k, v in out.items():
        log('g92rect', k, v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', default='ident')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--case', default='b1')
    ap.add_argument('--nrows', type=int, default=13)
    ap.add_argument('--zmax', type=int, default=250)
    ap.add_argument('--blist', default='0.25,0.5')
    ap.add_argument('--ylo', type=float, default=2.0)
    ap.add_argument('--yhi', type=float, default=20.0)
    a = ap.parse_args()
    R = {}
    here = os.path.dirname(os.path.abspath(__file__))
    t = a.task
    tag = t + ('_' + a.case if t == 'nonreal' else '')
    out = a.out or os.path.join(here, '..', 'data', 'mi_gb_%s.json' % tag)
    if t == 'ident': task_ident(R, a.case if a.case in ('a', 'b', 'c') else 'abc')
    if t == 'counts': task_counts(R, a.quick, a.blist.split(',') if a.blist != '0.25,0.5' else ('0.5', '0.75', '1', '1.5'))
    if t == 'thresh': task_thresh(R, [float(x) for x in a.blist.split(',')], a.quick)
    if t == 'nonreal': task_nonreal(R, a.case, a.zmax, a.quick, a.nrows)
    if t == 'g2': task_g2(R, a.quick)
    if t == 'g92': task_g92(R, a.ylo, a.yhi, a.quick)
    if t == 'g92rect': task_g92rect(R)
    json.dump(R, open(out, 'w'), indent=1, default=str)
    open(out.replace('.json', '.log'), 'w').write('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
