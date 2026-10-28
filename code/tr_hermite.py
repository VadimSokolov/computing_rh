# tr_hermite.py: the numerical claims of the Hermite coefficient part of Chapter 28 (open_problems note): the moment formula of thm:heat, the
# saddle point ratio table tab:ratio, the first negative index m*(a), the zero counting table tab:count with its edge lag, the sign agreement
# (82 percent), the amplitude ratio at a = 4, the run lengths table, the grid spacing arithmetic and the Lehmer pair.
# Conventions. Phi(u) = 2 sum_{n>=1} (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)) (the book's kernel, ch/ch02.tex:37), Xi(t) = int_R Phi(u) cos(t u) du = xi(1/2 + i t) with
# xi(s) = s (s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2; h_n the orthonormal Hermite functions (h_0 = pi^(-1/4) e^(-x^2/2), h_{n+1} = sqrt(2/(n+1)) x h_n
# - sqrt(n/(n+1)) h_{n-1}); c_n(a) = sqrt(a) int Phi(u) h_n(a u) du; b_m(a) = sqrt(2 pi) (-1)^m c_{2m}(a); t_m = a sqrt(4m+1); nu_k(a) = int Phi(u)
# exp(-a^2 u^2/2) u^(2k) du. A "sign change" is an index m with b_(m-1) b_m < 0, with t_m <= T; the count is compared with N(T) (mpmath nzeros).
# Method. All integrals by the trapezoid rule on the whole line in u (Phi even, so nodes u >= 0 with doubled weights), step h = 2 pi/(3T + 100)
# for the counting rows (the text's rule) and h/2 with 20 more digits as the second run; for the moment formula check h = 1/64 and 1/128 at 60
# digits. Digits 0.34 T + 40 for the counting rows. The nodes are spread over processes; every b_m for m < M comes from the recurrence in n.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Output: data/tr_hermite.json, log data/tr_hermite.log. Usage: python3 tr_hermite.py [--quick] [--out FILE] [--cpus N]; --quick: checks 1 to 3 and
# the T = 100 rows at a = 4 only, seconds. Full on Hopper (16 cpus, about 30 minutes):
# bash misc/tools/hopper_run.sh -c 16 -t 240 -g tr_hermite.json -g tr_hermite.log code/tr_hermite.py -- --out tr_hermite.json --cpus 16
import argparse, json, math, multiprocessing as mpr, time
import mpmath as mp

def Phi(u):
    e4 = mp.exp(2 * u); e5 = mp.exp(5 * u / 2); e9 = mp.exp(9 * u / 2)
    s = mp.mpf(0); n = 1
    while True:
        t = 2 * (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        s += t
        if abs(t) < mp.mpf(10) ** (-mp.mp.dps - 5) * (abs(s) + mp.mpf(10) ** -mp.mp.dps) and n > 2: break
        n += 1
    return s
def Xi(t):
    s = mp.mpf(1) / 2 + mp.mpc(0, t)
    return mp.re(s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2)

def nodes(h, dps, umax_extra=0):
    mp.mp.dps = dps
    out = []; j = 0
    while True:
        u = j * h; f = Phi(u)
        out.append((u, f * h * (1 if j else mp.mpf(1) / 2) * 2))
        if j > 5 and abs(f) < mp.mpf(10) ** (-dps - 5): break
        j += 1
    return out

def chunk_work(args):
    nd, a, M, dps = args
    mp.mp.dps = dps
    a = mp.mpf(a); N = 2 * M
    c1 = [mp.sqrt(mp.mpf(2) / (n + 1)) for n in range(N)]; c2 = [mp.sqrt(mp.mpf(n) / (n + 1)) for n in range(N)]
    S = [mp.mpf(0)] * M; p = mp.pi ** (-mp.mpf(1) / 4)
    for u, w in nd:
        x = a * u; hm = p * mp.exp(-x * x / 2); h = mp.sqrt(2) * x * hm
        S[0] += w * hm
        for n in range(1, N):
            if n % 2 == 0: S[n // 2] += w * h
            hm, h = h, c1[n] * x * h - c2[n] * hm
        # note: h now is h_{n+1}; the loop adds h_n before advancing
    return S

def coefficients(pool, a, M, h, dps):
    mp.mp.dps = dps
    nd = nodes(h, dps)
    k = pool._processes if hasattr(pool, '_processes') else 4
    chunks = [nd[i::k] for i in range(k)]
    parts = pool.map(chunk_work, [(c, a, M, dps) for c in chunks])
    mp.mp.dps = dps
    c = [mp.fsum(p[m] for p in parts) * mp.sqrt(a) for m in range(M)]      # c_{2m}
    return [mp.sqrt(2 * mp.pi) * (-1) ** m * c[m] for m in range(M)]

def moment_b(a, M, h, dps):
    mp.mp.dps = dps
    nd = nodes(h, dps)
    a = mp.mpf(a)
    nu = [mp.fsum(w * mp.exp(-a * a * u * u / 2) * u ** (2 * k) for u, w in nd) for k in range(M)]
    out = []
    for m in range(M):
        s = mp.fsum((-1) ** k * 2 ** k * a ** (2 * k) * nu[k] / (mp.factorial(2 * k) * 2 ** (m - k) * mp.factorial(m - k)) for k in range(m + 1))
        out.append(mp.sqrt(2 * mp.pi) * mp.pi ** (-mp.mpf(1) / 4) * mp.sqrt(a) * mp.sqrt(mp.factorial(2 * m)) * s)
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None); ap.add_argument('--cpus', type=int, default=4)
    a = ap.parse_args()
    out = a.out or ('/tmp/tr_hermite.json' if a.quick else 'data/tr_hermite.json')
    logf = open(out.replace('.json', '.log'), 'w'); R = {}
    def log(s): print(s); logf.write(s + '\n'); logf.flush()
    t0 = time.time(); pool = mpr.Pool(a.cpus)
    # 1. moment formula vs quadrature
    for aa in ('3', '1.5'):
        M = 16
        b1 = coefficients(pool, mp.mpf(aa), M, mp.mpf(1) / 64, 60); b2 = coefficients(pool, mp.mpf(aa), M, mp.mpf(1) / 128, 70)
        bm = moment_b(mp.mpf(aa), M, mp.mpf(1) / 64, 60)
        mp.mp.dps = 30
        d12 = max(abs(x - y) / abs(y) for x, y in zip(b1, b2)); dmo = max(abs(x - y) / abs(y) for x, y in zip(b1, bm))
        log(f"moment formula a={aa}, m<=15: max rel diff quadrature vs moment formula {mp.nstr(dmo, 3)}; quadrature h=1/64 vs 1/128 {mp.nstr(d12, 3)}; b_0 = {mp.nstr(b1[0], 12)}, b_15 = {mp.nstr(b1[15], 12)}")
        R[f'moment_{aa}'] = (str(dmo), str(d12))
    # Xi relation: b_m = int Xi(t) h_{2m}(t/a)/sqrt(a) dt for a = 1.5, m = 3 (30 digits)
    mp.mp.dps = 30
    a15 = mp.mpf('1.5'); b = coefficients(pool, a15, 6, mp.mpf(1) / 64, 30)
    def hfun(n, x):
        hm = mp.pi ** (-mp.mpf(1) / 4) * mp.exp(-x * x / 2); h = mp.sqrt(2) * x * hm
        if n == 0: return hm
        for k in range(1, n): hm, h = h, mp.sqrt(mp.mpf(2) / (k + 1)) * x * h - mp.sqrt(mp.mpf(k) / (k + 1)) * hm
        return h
    q = 2 * mp.quad(lambda t: Xi(t) * hfun(6, t / a15) / mp.sqrt(a15), mp.linspace(0, 60, 31))
    log(f"b_3(1.5) from Phi {mp.nstr(b[3], 15)}; from int Xi(t) h_6(t/a)/sqrt(a) dt {mp.nstr(q, 15)}")
    # 2. ratio table
    for aa in ('0.5', '0.75', '1'):
        af = mp.mpf(aa); M = int(((14.5 / af) ** 2 - 1) / 4) + 2
        mp.mp.dps = 50
        bb = coefficients(pool, af, M, mp.mpf(1) / 64, 50); mp.mp.dps = 50
        row = []
        for t in (3, 6, 9, 12, 13.5):
            m0 = int(round(float(((t / af) ** 2 - 1) / 4)))
            vals = []
            for m in (m0 - 1, m0, m0 + 1):
                if m < 1: continue
                tm = af * mp.sqrt(4 * m + 1)
                vals.append((m, float(tm), float(bb[m] / (mp.sqrt(2 * af) * mp.mpf(m) ** (-mp.mpf(1) / 4) * Xi(tm)))))
            row.append(vals)
        for t, v in zip((3, 6, 9, 12, 13.5), row):
            log(f"ratio a={aa} t~{t}: " + "; ".join(f"m={m} t_m={tm:.3f} ratio={r:.4f}" for m, tm, r in v))
        r9 = row[2][1][2]; log(f"  error at t~9 over a^2: (1 - ratio)/a^2 = {(1 - r9) / float(af) ** 2:.4f}")
        R[f'ratio_{aa}'] = row
    # 3. first negative index
    g1 = mp.findroot(Xi, 14.1347)
    for aa in ('0.5', '0.75', '1', '1.5', '2'):
        af = mp.mpf(aa); M = int((g1 / af) ** 2 / 4) + 60 if aa != '2' else 40
        mp.mp.dps = 50
        bb = coefficients(pool, af, M, mp.mpf(1) / 64, 50); mp.mp.dps = 50
        neg = next(m for m in range(M) if bb[m] < 0)
        bb2 = coefficients(pool, af, M, mp.mpf(1) / 128, 60); mp.mp.dps = 50
        neg2 = next(m for m in range(M) if bb2[m] < 0)
        pos_before = all(bb[m] > 0 for m in range(neg))
        log(f"a={aa}: (g1^2/a^2-1)/4 = {mp.nstr((g1 ** 2 / af ** 2 - 1) / 4, 6)}; first negative m = {neg} (second run {neg2}); b_m>0 for all m below it: {pos_before}; ceil formula {int(mp.ceil((g1 ** 2 / af ** 2 - 1) / 4))}; b_neg = {mp.nstr(bb[neg], 5)}")
        R[f'firstneg_{aa}'] = neg
    # 4. counting
    rows = [(100, [2, 3, 4, 5]), (200, [4, 5, 6, 7, 8]), (500, [6, 7, 8, 9, 10, 14, 18]), (1000, [5, 6]), (1004, [6])]
    if a.quick: rows = [(100, [4])]
    mp.mp.dps = 30
    keep = {}
    for T, As in rows:
        N = int(mp.nzeros(T)); dps = int(0.34 * T + 40)
        for aa in As:
            M = int(((T / aa) ** 2 - 1) / 4) + 1
            h = 2 * mp.pi / (3 * T + 100)
            b = coefficients(pool, mp.mpf(aa), M, h, dps); mp.mp.dps = 30
            tm = [aa * math.sqrt(4 * m + 1) for m in range(M)]
            sc = [m for m in range(1, M) if b[m - 1] * b[m] < 0]
            line = f"T={T} N(T)={N} a={aa} digits {dps} coefficients {M} (formula T^2/(4a^2) = {T * T / (4 * aa * aa):.1f}): sign changes {len(sc)}"
            if T <= 200 or (T == 500 and aa == 6):
                b2 = coefficients(pool, mp.mpf(aa), M, h / 2, dps + 20); mp.mp.dps = 30
                sc2 = [m for m in range(1, M) if b2[m - 1] * b2[m] < 0]
                line += f"; second run (h/2, +20 digits) {len(sc2)}, same positions {sc == sc2}"
            log(line + f"; grid spacing at the top 2a^2/T = {2 * aa * aa / T:.3f}, mean gap 2pi/log(T/2pi) = {2 * math.pi / math.log(T / (2 * math.pi)):.3f}, ratio {2 * aa * aa / T / (2 * math.pi / math.log(T / (2 * math.pi))):.3f}")
            R[f'count_{T}_{aa}'] = (N, len(sc))
            if (T, aa) in ((200, 4), (1004, 6)): keep[(T, aa)] = (b, tm, sc)
    # 5. sign agreement and runs at T=200, a=4; amplitude ratio; edge lag at 1004
    if (200, 4) in keep:
        b, tm, sc = keep[(200, 4)]
        mp.mp.dps = 30
        for c in (0, 1, 2, 3):
            ag = sum(1 for m in range(len(b)) if (b[m] > 0) == (Xi(4 * mp.sqrt(4 * m + c)) > 0)) / len(b)
            log(f"sign agreement a=4, T=200, sampling a sqrt(4m+{c}): {100 * ag:.1f} percent of {len(b)} indices")
        rat = [float(b[m] / (mp.sqrt(8) * mp.mpf(m) ** (-mp.mpf(1) / 4) * Xi(tm[m]))) for m in range(1, len(b)) if 100 < tm[m] < 200]
        log(f"amplitude ratio b_m/(sqrt(2a) m^(-1/4) Xi(t_m)) at a=4, 100<t<200: min {min(rat):.3f}, max {max(rat):.3f}, abs min {min(abs(r) for r in rat):.3f}, abs max {max(abs(r) for r in rat):.3f}")
        for lo, hi in ((40, 80), (80, 120), (120, 160), (160, 200)):
            idx = [m for m in range(len(b)) if lo <= tm[m] < hi]
            blocks = 1 + sum(1 for m in idx[1:] if (b[m] > 0) != (b[m - 1] > 0))
            tmid = (lo + hi) / 2
            log(f"runs a=4 window {lo} to {hi}: {len(idx)} indices, {blocks} sign blocks, mean run {len(idx) / blocks:.2f}; predicted pi t/(a^2 log(t/2pi)) at t={tmid}: {math.pi * tmid / (16 * math.log(tmid / (2 * math.pi))):.2f}")
    if (1004, 6) in keep:
        b, tm, sc = keep[(1004, 6)]
        last = sc[-5:]
        log(f"T=1004 a=6 last sign changes (t_(m-1), t_m): {[(round(tm[m - 1], 2), round(tm[m], 2)) for m in last]}")
        mp.mp.dps = 30
        log(f"zeros 649..652: {[mp.nstr(mp.zetazero(k).imag, 8) for k in range(649, 653)]}")
    mp.mp.dps = 30
    if not a.quick:
        z = [mp.zetazero(k).imag for k in (6709, 6710)]
        log(f"Lehmer pair zeros 6709, 6710: {mp.nstr(z[0], 10)}, {mp.nstr(z[1], 10)}, gap {mp.nstr(z[1] - z[0], 4)}")
    json.dump(R, open(out, 'w'), default=str); log(f"time {time.time() - t0:.1f} s")
if __name__ == '__main__':
    main()
