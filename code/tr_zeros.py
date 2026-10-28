# tr_zeros.py: the zero and Gram claims of Chapter 28 (Truncations and Hermite Coefficients), Section tr: "The side of approach is Gram's law":
# the 379 zeros of L(s,chi_4) below height 500, the five Gram anomalies (138.750, 290.252, 364.420, 400.826, 422.510) and their close pair
# gaps, the drift coefficients delta_m (thm:drift formula) for m = 0..5 at every zero, the signs of Re delta_m, and the phase law thm:phase
# (sign Re delta_m = sign(-sin Theta_m Z')) tested zero by zero. Also delta_0, delta_1 at the first zero (the table of Section tr).
# Conventions. L(s) = L(s,chi_4) = 4^(-s) (zeta(s,1/4) - zeta(s,3/4)); Lambda(s) = (4/pi)^(s/2) Gamma((s+1)/2) L(s);
# vartheta(x) = (x/2) log(4/pi) + Im loggamma((3/2 + i x)/2); Z(x) = Re(exp(i vartheta) L(1/2 + i x)) (the imaginary part is logged);
# Gram points: vartheta(g_j) = j pi; zero index k = 1, 2, ...; Gram interval index of a zero j(k) = floor(vartheta(gamma_k)/pi); the offset is
# j(k) - k, the expected offset is the most frequent one, the anomalies are the zeros with another offset.
# delta_m = (-1)^(m+1) (z0)_(m+1) Lambda(rho+2m+2) / (2 (m+1)(2m+1) pi^(m+1) Lambda'(rho)), z0 = rho/2, Lambda'(rho) = (4/pi)^(rho/2)
# Gamma((rho+1)/2) L'(rho) at the zero. Z'(gamma) = i Lambda'(rho)/|G(rho)| (from Lambda' = -i |G| Z').
# Method. Z on two independent grids (steps 0.02 and 0.0137), every sign change refined by the Illinois method to 35 digits and recomputed at
# 45 digits (the difference is logged); the two grids must give the same zeros. Not certified: no Turing bound, no ball arithmetic; the
# count is compared with vartheta(T)/pi and the minimum gap (0.134) with the step.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Output: data/tr_zeros.json, log data/tr_zeros.log. Usage: python3 tr_zeros.py [--quick] [--out FILE] [--cpus N] [--height H]
# --quick: height 60, one grid, seconds. Full: height 500 on Hopper, 16 cpus, about 10 minutes:
# bash misc/tools/hopper_run.sh -c 16 -t 60 -g tr_zeros.json -g tr_zeros.log code/tr_zeros.py -- --out tr_zeros.json
import argparse, json, multiprocessing as mpr, sys, time
import mpmath as mp

def L(s):
    return mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf(1) / 4) - mp.zeta(s, mp.mpf(3) / 4))
def theta(x):
    return x / 2 * mp.log(4 / mp.pi) + mp.im(mp.loggamma((mp.mpf(3) / 2 + mp.mpc(0, x)) / 2))
def Zc(x):
    return mp.exp(mp.mpc(0, 1) * theta(x)) * L(mp.mpf(1) / 2 + mp.mpc(0, x))
def Z(x):
    return mp.re(Zc(x))
def Lam(s):
    return (4 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * L(s)
def Lprime(s):
    return -mp.log(4) * L(s) + mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf(1) / 4, 1) - mp.zeta(s, mp.mpf(3) / 4, 1))
def Lamprime_at(rho):
    return (4 / mp.pi) ** (rho / 2) * mp.gamma((rho + 1) / 2) * Lprime(rho)

def grid_chunk(args):
    xs, dps = args
    mp.mp.dps = dps
    return [float(Z(mp.mpf(x))) for x in xs]
def refine(args):
    a, b, dps = args
    mp.mp.dps = dps
    a, b = mp.mpf(a), mp.mpf(b)
    r = mp.findroot(Z, (a, b), solver='illinois', tol=mp.mpf(10) ** (-(dps - 8)), maxsteps=200)
    return str(r)

def zeros_on_grid(pool, T, step, dps, log):
    n = int(T / step) + 1
    xs = [0.01 + i * step for i in range(n)]
    chunks = [xs[i:i + 50] for i in range(0, n, 50)]
    vals = [v for c in pool.map(grid_chunk, [(c, 20) for c in chunks]) for v in c]
    br = [(xs[i], xs[i + 1]) for i in range(n - 1) if vals[i] * vals[i + 1] < 0]
    log(f"grid step {step}: {n} points, {len(br)} sign changes up to {xs[-1]:.2f}")
    return br

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None)
    ap.add_argument('--cpus', type=int, default=8); ap.add_argument('--height', type=float, default=None)
    a = ap.parse_args()
    T = a.height or (40.0 if a.quick else 500.0)
    out = a.out or ('/tmp/tr_zeros.json' if a.quick else 'data/tr_zeros.json')
    logf = open(out.replace('.json', '.log'), 'w')
    def log(s):
        print(s); logf.write(s + '\n'); logf.flush()
    t0 = time.time()
    pool = mpr.Pool(a.cpus)
    steps = [0.05] if a.quick else [0.02, 0.0137]
    brs = [zeros_on_grid(pool, T, st, 20, log) for st in steps]
    mp.mp.dps = 35
    r1 = [mp.mpf(s) for s in pool.map(refine, [(x, y, 35) for x, y in brs[0]])]
    r2 = [mp.mpf(s) for s in pool.map(refine, [(x, y, 45) for x, y in brs[0]])]
    dprec = max(abs(x - y) for x, y in zip(r1, r2))
    log(f"refined {len(r1)} zeros; max difference 35 vs 45 digits = {mp.nstr(dprec, 3)}")
    if len(brs) > 1:
        rb = [mp.mpf(s) for s in pool.map(refine, [(x, y, 35) for x, y in brs[1]])]
        same = len(rb) == len(r1) and max(abs(x - y) for x, y in zip(rb, r1)) < mp.mpf(10) ** -25
        log(f"second grid: {len(rb)} zeros; identical to the first: {same}")
    gam = [x for x in r1 if x < T]
    N = len(gam)
    log(f"N({T}) = {N}; vartheta(T)/pi = {mp.nstr(theta(T) / mp.pi, 8)}; min gap = {mp.nstr(min(gam[i + 1] - gam[i] for i in range(N - 1)), 4)}")
    log(f"Z imaginary part at gamma_1, gamma_50: {mp.nstr(mp.im(Zc(gam[0])), 3)}, {mp.nstr(mp.im(Zc(gam[min(49, N - 1)])), 3)}")
    # Gram
    off = [int(mp.floor(theta(g) / mp.pi)) - (k + 1) for k, g in enumerate(gam)]
    from collections import Counter
    cnt = Counter(off); mode = cnt.most_common(1)[0][0]
    anom = [k for k in range(N) if off[k] != mode]
    log(f"offsets {dict(cnt)}; expected offset {mode}; {N - len(anom)} with it, {len(anom)} anomalous")
    gaps = {}
    for k in anom:
        gl = gam[k] - gam[k - 1] if k > 0 else None; gr = gam[k + 1] - gam[k] if k + 1 < N else None
        gaps[k] = (gl, gr)
        log(f"anomaly zero #{k + 1} gamma = {mp.nstr(gam[k], 10)} offset {off[k]}; gap left {mp.nstr(gl, 4)} gap right {mp.nstr(gr, 4)}")
    # delta_m
    M = 6
    neg = {m: [] for m in range(M)}; phase_ok = {m: 0 for m in range(M)}; rows = []
    for k, g in enumerate(gam):
        rho = mp.mpf(1) / 2 + mp.mpc(0, g); z0 = rho / 2
        lp = Lamprime_at(rho); G = abs((4 / mp.pi) ** (rho / 2) * mp.gamma((rho + 1) / 2))
        Zp = mp.re(mp.mpc(0, 1) * lp / G)
        row = {'gamma': str(g)}
        for m in range(M):
            d = (-1) ** (m + 1) * mp.rf(z0, m + 1) * Lam(rho + 2 * m + 2) / (2 * (m + 1) * (2 * m + 1) * mp.pi ** (m + 1) * lp)
            Th = theta(g) + mp.arg(L(rho + 2 * m + 2)) - sum(mp.atan((mp.mpf(1) / 2 + 2 * i) / g) + mp.atan((mp.mpf(3) / 2 + 2 * i) / g) for i in range(m + 1))
            pred = -mp.sin(Th) * Zp
            if d.real < 0: neg[m].append(k)
            if (d.real < 0) == (pred < 0): phase_ok[m] += 1
            row[f'd{m}'] = [str(d.real), str(d.imag)]
        rows.append(row)
        if k == 0:
            log(f"first zero gamma_1 = {mp.nstr(g, 12)}: delta_0 = {mp.nstr(mp.mpc(*map(mp.mpf, rows[0]['d0'])), 8)}, delta_1 = {mp.nstr(mp.mpc(*map(mp.mpf, rows[0]['d1'])), 8)}")
    for m in range(M):
        gl = [float(gam[k]) for k in neg[m]]
        log(f"m = {m}: Re delta_m < 0 at {len(neg[m])} zeros; phase law sign agrees at {phase_ok[m]} of {N}; "
            + (f"ordinates {[round(x, 3) for x in gl]}" if len(gl) <= 12 else f"{sum(1 for x in gl if x < 105)} below 105, last {[round(x, 3) for x in gl[-5:]]}"))
    json.dump({'N': N, 'gamma': [str(x) for x in gam], 'offsets': off, 'mode': mode, 'anomalies': [k + 1 for k in anom],
               'neg': {m: [k + 1 for k in neg[m]] for m in neg}, 'phase_ok': phase_ok, 'rows': rows, 'prec_diff': str(dprec)}, open(out, 'w'))
    log(f"time {time.time() - t0:.1f} s")
if __name__ == '__main__':
    main()
