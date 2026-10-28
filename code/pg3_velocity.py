# pg3_velocity.py: the velocity of the zeros of H_t at t = 0 (Proposition prop:ch21:velocity and the paragraph after it, Chapter 21,
# merged from note X3; placeholders N2a, N2b, N2c, N2d, N2e, N2f, N2g, N2h). Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. Xi(theta) = xi(1/2 + i theta), xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, gamma_n the n-th positive zero, and the
# zeros of H_0(z) = Xi(z/2)/8 are z_n = 2 gamma_n. The flow dz_k/dt = H''/H' gives zdot_n(0) = Xi''(gamma_n)/(2 Xi'(gamma_n)), the book's
# units (z = 2 theta), so the values need no conversion from the note's.
# Method A (python-flint, Hopper, certified). gamma_n are balls of acb.zeta_zeros(1, N) at prec bits (isolated and counted by Arb, so
# every zero is simple); zdot_n is c2/c1 for the Taylor jet c0 + c1 h + c2 h^2 of Xi(gamma_n + h), from acb_series of s(s-1)/2,
# pi^(-s/2), Gamma(s/2), zeta(s) at s = 1/2 + i gamma_n + i h over the ball of gamma_n (an enclosure of the value at the true zero).
# The signs are certified when the ball excludes 0. The census (N = 1500) uses the same balls; its statistics (count of upward zeros,
# mean, extremes and their indices, the distance of the extreme pair) are computed from the midpoints and the extremes are certified
# separated from the next values. The run is repeated at a second precision and the maximum width of the balls and the maximum
# difference of the midpoints are logged. The ordinates are compared with data/zeros_8000.txt.
# Method B (mpmath, 45 digits, computed). For selected n the formula of the proposition,
#   zdot = (1/2)[Z''/Z' + 4 g/(g^2+1/4) - Im psi(1/4 + i g/2)],  Z''/Z' = i (2 theta'(g) + zeta''/zeta'(1/2 + i g)),
#   2 theta'(g) = Re psi(1/4 + i g/2) - log pi  (equivalently Im psi enters through the derivative of |Gamma|),
# is compared with Xi''/(2 Xi') by numerical differentiation of Xi, and gamma_1 zdot_1 with adot_1(0) = -8.193 (ch21.tex:620).
# Options: --quick (n <= 35 by method B at 30 digits, no flint, seconds), --cert (method A; needs python-flint), --flow (task N2c at t = 0.2, own output file pg3_velocity_flow.json), --N N (census size),
# --out FILE. Method B alone runs on a laptop in the 30 second limit with --quick.
import sys, json, argparse, os, time
import mpmath as mp

def zfile(n):
    for d in ('.', '../data', 'data'):
        f = os.path.join(d, 'zeros_8000.txt')
        if os.path.exists(f):
            return [l.strip() for l in open(f).read().split('\n') if l.strip()][:n]
    raise SystemExit('zeros_8000.txt not found')

def zdot_formula(g, dps):
    mp.mp.dps = dps
    g = mp.mpf(g)
    s = mp.mpf(1) / 2 + 1j * g
    z1 = mp.zeta(s, derivative=1); z2 = mp.zeta(s, derivative=2)
    psi = mp.digamma(mp.mpf(1) / 4 + 1j * g / 2)
    two_th = mp.re(psi) - mp.log(mp.pi)
    zz = (1j * (two_th + z2 / z1))          # Z''/Z' at a zero, real
    return mp.re(mp.mpf(1) / 2 * (zz + 4 * g / (g * g + mp.mpf(1) / 4) - mp.im(psi))), mp.im(zz)

def zdot_xi(g, dps):
    mp.mp.dps = dps + 15
    def Xi(th):
        s = mp.mpf(1) / 2 + 1j * th
        return mp.re(s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2)
    g = mp.mpf(g)
    d1 = mp.diff(Xi, g, 1); d2 = mp.diff(Xi, g, 2)
    r = d2 / (2 * d1)
    mp.mp.dps = dps
    return +r

def cert(N, prec, zs):
    from flint import arb, acb, acb_series, ctx
    ctx.prec = prec
    Z = acb.zeta_zeros(1, N)
    half = arb(1) / 2
    assert all(z.real == half for z in Z)
    G = [z.imag for z in Z]
    lp = arb.pi().log()
    out = []
    for g in G:
        s0 = acb(half, g)
        s = acb_series([s0, acb(0, 1)], prec=3)
        f = s * (s - 1) / 2 * (-(s / 2) * lp).exp() * (s / 2).gamma() * s.zeta()
        c = f.coeffs()
        out.append(c[2] / c[1])
    return G, out

def run_cert(N, log):
    from flint import arb, ctx
    zs = zfile(N)
    res = {}
    t0 = time.time()
    G, Z1 = cert(N, 256, zs)
    t1 = time.time()
    G2, Z2 = cert(N, 384, zs)
    log('Arb zeros and jets for n <= %d: %.0f s (256 bits), %.0f s (384 bits)' % (N, t1 - t0, time.time() - t1))
    mp.mp.dps = 40
    for n in range(N):
        assert abs(mp.mpf(G[n].mid().str(45, radius=False)) - mp.mpf(zs[n])) < mp.mpf(10) ** -40, n
    log('ordinates agree with data/zeros_8000.txt to 1e-40 for n <= %d' % N)
    imag_ok = all(z.imag.contains(0) for z in Z1)
    log('every jet quotient is real to within its ball: %s' % imag_ok)
    zr = [z.real for z in Z1]
    zr2 = [z.real for z in Z2]
    maxrad = max(float(z.rad()) for z in zr)
    mpv = [mp.mpf(z.mid().str(30, radius=False)) for z in zr]
    mpv2 = [mp.mpf(z.mid().str(30, radius=False)) for z in zr2]
    maxdiff = max(abs(a - b) for a, b in zip(mpv, mpv2))
    log('largest ball radius of zdot at 256 bits %.3g; largest difference of midpoints between 256 and 384 bits %s' % (maxrad, mp.nstr(maxdiff, 3)))
    down = [n + 1 for n in range(min(N, 35)) if zr[n].upper() < 0]
    up = [n + 1 for n in range(N) if zr[n].lower() > 0]
    undecided = [n + 1 for n in range(N) if not (zr[n].upper() < 0 or zr[n].lower() > 0)]
    log('certified downward for n = 1..34: %s; certified upward n = 35: %s; signs excluding 0 as balls: %d of %d (the certificate is claimed for n <= 35 only); undecided: %s'
        % (all(n in down for n in range(1, 35)), 35 in up, len(up) + len([n for n in range(N) if zr[n].upper() < 0]), N, undecided))
    log('upward zeros among n <= 35 (certified): %s' % [n for n in up if n <= 35])
    g35 = G[34]; z35 = zr[34]
    log('gamma_35 = %s ; zdot_35(0) = %s' % (g35.str(15, radius=False), z35.str(15)))
    res['gamma_35'] = g35.mid().str(20, radius=False)
    res['zdot_35'] = z35.mid().str(20, radius=False)
    res['zdot_35_ball'] = z35.str(25)
    res['signs_n_le_35_certified'] = all(n in down for n in range(1, 35)) and 35 in up
    res['max_ball_radius_256'] = maxrad
    res['max_midpoint_diff_256_384'] = mp.nstr(maxdiff, 3)
    # census
    cnt = len(up)
    mean = mp.fsum(mpv) / N
    imax = max(range(N), key=lambda i: mpv[i]); imin = min(range(N), key=lambda i: mpv[i])
    second_max = max(mpv[i] for i in range(N) if i != imax); second_min = min(mpv[i] for i in range(N) if i != imin)
    sep_max = zr[imax].lower() > max(zr[i].upper() for i in range(N) if i != imax)
    sep_min = zr[imin].upper() < min(zr[i].lower() for i in range(N) if i != imin)
    gi, gj = mp.mpf(G[imax].mid().str(30, radius=False)), mp.mpf(G[imin].mid().str(30, radius=False))
    log('census n <= %d: %d upward (N2c); mean %s (N2d); maximum %s at n = %d (N2e, N2h), minimum %s at n = %d (N2f, N2h); '
        'extremes certified separated from all other values: %s, %s; distance of the two zeros in gamma %s (N2g)'
        % (N, cnt, mp.nstr(mean, 8), mp.nstr(mpv[imax], 8), imax + 1, mp.nstr(mpv[imin], 8), imin + 1, sep_max, sep_min, mp.nstr(abs(gi - gj), 6)))
    log('second largest %s, second smallest %s; drift -pi/4 = %s' % (mp.nstr(second_max, 6), mp.nstr(second_min, 6), mp.nstr(-mp.pi / 4, 6)))
    ups = [n + 1 for n in range(N) if zr[n].lower() > 0]
    log('first upward indices: %s' % ups[:12])
    log('min |zdot_n| for n <= 34: %s at n = %d' % (mp.nstr(min(abs(v) for v in mpv[:34]), 6), 1 + min(range(34), key=lambda i: abs(mpv[i]))))
    res.update(N2c=cnt, N2d=mp.nstr(mean, 8), N2e=mp.nstr(mpv[imax], 8), N2f=mp.nstr(mpv[imin], 8), N2h=[imax + 1, imin + 1],
               N2g=mp.nstr(abs(gi - gj), 6), extremes_separated=bool(sep_max and sep_min), first_upward=ups[:12],
               undecided=undecided, N=N, certified=False, certified_n_le_35=bool(res["signs_n_le_35_certified"]))
    res['second_max'] = mp.nstr(second_max, 6); res['second_min'] = mp.nstr(second_min, 6)
    res['zdot_first12'] = [mp.nstr(v, 8) for v in mpv[:12]]
    return res

def run_mp(ns, dps, log):
    zs = zfile(max(ns))
    rows = {}
    for n in ns:
        a, im = zdot_formula(zs[n - 1], dps)
        b = zdot_xi(zs[n - 1], dps)
        log('n=%d gamma=%s  formula %s  Xi\'\'/(2Xi\') %s  |diff| %s (imag part of Z\'\'/Z\' %s)'
            % (n, mp.nstr(mp.mpf(zs[n - 1]), 12), mp.nstr(a, 15), mp.nstr(b, 15), mp.nstr(abs(a - b), 3), mp.nstr(im, 3)))
        rows[n] = dict(gamma=zs[n - 1][:30], formula=mp.nstr(a, 25), xi=mp.nstr(b, 25), diff=mp.nstr(abs(a - b), 3))
    if 1 in ns:
        mp.mp.dps = dps
        a1 = mp.mpf(zs[0]) * mp.mpf(rows[1]['formula'])
        log('gamma_1 zdot_1 = %s (ch21.tex:620: adot_1(0) = -8.193)' % mp.nstr(a1, 10))
        rows['a1dot'] = mp.nstr(a1, 12)
    return rows

def main1():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--cert', action='store_true')
    ap.add_argument('--N', type=int, default=1500); ap.add_argument('--out', default='pg3_velocity.json')
    a = ap.parse_args()
    L = []
    def log(s):
        print(s); sys.stdout.flush(); L.append(s)
    res = {}
    if a.quick:
        res['mp'] = run_mp([1, 10, 34, 35], 30, log)
    else:
        res['mp'] = run_mp([1, 10, 34, 35, 1496, 1497], 45, log)
    if a.cert and not a.quick:
        res['cert'] = run_cert(a.N, log)
    json.dump(res, open(a.out, 'w'), indent=1)
    open(a.out.replace('.json', '.log'), 'w').write('\n'.join(L) + '\n')

# ---- task N2c: the velocities at t = 0.2 (computed, not certified) ------------------------------------------------------------------------
# Xi_t(theta) = E xi(1/2 + i theta + sqrt(t/2) Z) (Z standard normal; the same flow, e^(-(t/4) d^2/dtheta^2) Xi, as dz/dt = H''/H' in z = 2 theta),
# by Gauss Hermite quadrature with K nodes (double precision nodes of numpy, weights summing to sqrt(pi)): the sum of the Taylor jets of xi
# at the shifted points (acb_series, python-flint) gives Xi_t, Xi_t', Xi_t''/2 at once. The zeros of H_t are tracked from t = 0 to 0.2 in 4 steps
# (predictor theta += dt zdot/2, then Newton), and zdot_n(0.2) = Xi_t''/(2 Xi_t') there. K = 60 against K = 90 gives the quadrature error;
# the order of the tracked zeros is checked to stay increasing. Run with --flow (python-flint needed; uses the first NF = 100 zeros).
def run_flow(NF, log, K1=60, K2=90, steps=4, tend=0.2):
    import numpy as np
    from flint import arb, acb, acb_series, ctx
    ctx.prec = 160
    lp = arb.pi().log()
    zs = zfile(NF)
    mp.mp.dps = 30
    cache = {}
    def nodes(K):
        if K not in cache:
            cache[K] = np.polynomial.hermite.hermgauss(K)
        return cache[K]
    def jets(theta, t, K):
        y, w = nodes(K)
        c0 = c1 = c2 = arb(0)
        th = arb(mp.nstr(theta, 28))
        for yj, wj in zip(y, w):
            s0 = acb(arb(1) / 2 + arb(float(np.sqrt(t) * yj)), th)
            s = acb_series([s0, acb(0, 1)], prec=3)
            f = s * (s - 1) / 2 * (-(s / 2) * lp).exp() * (s / 2).gamma() * s.zeta()
            c = f.coeffs()
            ww = arb(float(wj))
            c0 += ww * c[0].real; c1 += ww * c[1].real; c2 += ww * c[2].real
        sq = arb.pi().sqrt()
        conv = lambda x: mp.mpf((x / sq).mid().str(25, radius=False))
        return conv(c0), conv(c1), conv(c2)
    def track(K):
        out = []
        for n in range(NF):
            th = mp.mpf(zs[n]); t = 0.0
            for i in range(1, steps + 1):
                tn = tend * i / steps
                c0, c1, c2 = jets(th, t, K)
                th = th + (tn - t) * (c2 / c1) / 2
                for _ in range(8):
                    c0, c1, c2 = jets(th, tn, K)
                    d = c0 / c1; th = th - d
                    if abs(d) < mp.mpf(10) ** -13: break
                t = tn
            c0, c1, c2 = jets(th, tend, K)
            out.append((th, c2 / c1))
        return out
    t0 = time.time()
    A = track(K1)
    log('flow: K = %d nodes, %d zeros tracked to t = %s in %.0f s' % (K1, NF, tend, time.time() - t0))
    B = track(K2)
    dz = max(abs(a[1] - b[1]) for a, b in zip(A, B)); dth = max(abs(a[0] - b[0]) for a, b in zip(A, B))
    increasing = all(A[i][0] < A[i + 1][0] for i in range(NF - 1))
    ups = [(i + 1, A[i][1]) for i in range(NF) if A[i][1] > 0]
    zd0 = [mp.mpf(0)] * 0
    log('max difference of zdot_n(0.2) between K = %d and K = %d: %s; of the zero positions: %s; order of the zeros preserved: %s' % (K1, K2, mp.nstr(dz, 3), mp.nstr(dth, 3), increasing))
    mx = max(range(NF), key=lambda i: A[i][1])
    log('zdot_n(0.2), n <= %d: %d upward, %s; largest %s at n = %d (theta = %s); first ten: %s' %
        (NF, len(ups), [(n, mp.nstr(v, 6)) for n, v in ups][:12], mp.nstr(A[mx][1], 8), mx + 1, mp.nstr(A[mx][0], 12), [mp.nstr(a[1], 5) for a in A[:10]]))
    return dict(K=[K1, K2], diff_zdot=mp.nstr(dz, 3), diff_theta=mp.nstr(dth, 3), order_preserved=bool(increasing),
                upward=[(n, mp.nstr(v, 10)) for n, v in ups], N2i=[(n, mp.nstr(v, 8)) for n, v in ups][:5],
                zdot_02=[mp.nstr(a[1], 10) for a in A], theta_02=[mp.nstr(a[0], 14) for a in A], certified=False)

def main2():
    ap = argparse.ArgumentParser()
    ap.add_argument('--flow', action='store_true'); ap.add_argument('--NF', type=int, default=100)
    ap.add_argument('--out', default='pg3_velocity_flow.json')
    a, _ = ap.parse_known_args()
    L = []
    def log(s):
        print(s); sys.stdout.flush(); L.append(s)
    res = run_flow(a.NF, log)
    json.dump(res, open(a.out, 'w'), indent=1)
    open(a.out.replace('.json', '.log'), 'w').write('\n'.join(L) + '\n')

if __name__ == '__main__':
    if '--flow' in sys.argv:
        main2()
    else:
        main1()
