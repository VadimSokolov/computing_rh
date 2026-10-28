# ev_g1_fold.py: the numerical claims of Nick Polson's extreme value hitting times chapter (E; Chapter 25 after the merge) about the folded
# Gumbel family H^F_N and the Polya kernel: tab:approx (rows folded N = 2, 3, 4 and two sided), E:103 (kink law kappa_N and the breakdown
# heights), E:294 (B_N), prop:rouche and tab:rouche (E:296 to E:322), the margins at E:322, the statement of prop:reduction that the
# hypothesis holds for N <= 4 (E:271), and the question of A1 M3 whether the "Re" form of E:74 is right off the real axis.
# Added for the book in October 2026 by computation group G1, written from the mathematics (the collaborator's scripts were not read, run
# or copied).
# Conventions. Phi(u) = sum_n phi_n(u), phi_n(u) = (2 a^2 e^{9u} - 3 a e^{5u}) exp(-a e^{4u}), a = pi n^2, u >= 0 (E:49). Phi^F_N = sum_{n<=N} phi_n(|u|),
# R_N = sum_{n>N} phi_n on u >= 0, H^F_N(z) = int_0^inf Phi^F_N(u) cos(zu) du, E_N = H - H^F_N, H(z) = Xi(z/2)/8, Xi(t) = xi(1/2 + i t),
# xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2. F_n(z) = int_0^inf phi_n(u) e^{izu} du = a^(2-w9) Gamma(w9,a)/2 - 3 a^(1-w5) Gamma(w5,a)/4 with
# w_j = (j + i z)/4 (substitution y = a e^{4u}); the cosine transform is (F_n(z) + F_n(-z))/2 for every complex z, which equals Re F_n(z) only for
# real z (this script also evaluates the "Re" form of E:74 at complex z to see what it would give). Ordinates gamma_k from data/zeros_1700.txt.
# t = z/2 is the variable of the zeros gamma_k; "below 70" in tab:approx is read in t (the count in z would be 5 for N = 2, not 7).
# Method. Everything is in ball arithmetic (python-flint acb, incomplete gamma with complex parameter). E_N(z) = sum_{n=N+1}^{N+5} cosine transform +
# a rigorous tail: for n > N+5, a_n >= 2(w-1) with w = (9 + |Im z|)/4, so Gamma(w,a) <= 2 a^(w-1) e^(-a) and |tail_n| <= (a_n + 3/2) e^(-a_n) (checked in
# the code), summed with a factor 2. Sign changes of H^F_N(2t) on a grid in t: a sign is used only if the ball excludes 0, so the count is a certified
# lower bound for the number of real zeros (gaps between grid points are not excluded: 'computed'). The Rouche inequality |E_N| < |H| on the quarter circle
# z = 2T e^{i theta}, theta in [0, pi/2], is certified by interval cells in theta (z in the enclosing box), bisected until |E_N|_upper < |H|_lower on every
# cell: this is a CERTIFIED check of the inequality (given the cosine transform formulas above), reported next to the sampled 721 point minimum. The
# count of N(T) sign changes is certified as a lower bound by the grid; that there are no more than N(T) real zeros follows from the Rouche count
# 2 N(T) zeros of H in the disc once the inequality holds (prop:rouche). B_N = |R_N'(0)| + int |R_N''| is certified by adaptive cells on which R_N'' has a
# certified sign (Arb integration on each) plus width * sup on the cells at which the sign cannot be decided, plus the rigorous majorant of the tail.
# Error control. The 721 point values are computed at two precisions (bits and 2 bits) and with N+5 and N+8 omitted terms; the logged difference is
# the numerical error; the cell certificates carry their own radii.
# Output: data/ev_g1_fold.json and data/ev_g1_fold.log (--out FILE for the json, the log goes next to it).
# Usage: python3 ev_g1_fold.py [--quick] [--out FILE] [--part approx|rouche|bn|polya|all]; --quick uses the coarse grid and fewer rows and runs in a few
# seconds; full run: bash misc/tools/hopper_run.sh -c 1 -t 60 -g ev_g1_fold.json -g ev_g1_fold.log code/ev_g1_fold.py data/zeros_1700.txt -- --out ev_g1_fold.json
import argparse, json, math, os, sys, time
from flint import arb, acb, ctx

START = time.time()
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


def load_zeros(path):
    return [float(l.split()[0]) for l in open(path) if l.strip()]


def Nof(T, gam):
    return sum(1 for g in gam if g < T)


def pi_():
    return arb.pi()


def F_n(n, z):
    """int_0^inf phi_n(u) e^{izu} du for an acb z."""
    a = pi_() * n * n
    ac = acb(a)
    w9 = (acb(9) + acb(0, 1) * z) / 4
    w5 = (acb(5) + acb(0, 1) * z) / 4
    return ac ** (2 - w9) * ac.gamma_upper(w9) / 2 - ac ** (1 - w5) * ac.gamma_upper(w5) * 3 / 4


def cos_n(n, z):
    return (F_n(n, z) + F_n(n, -z)) / 2


def cos_n_real(n, x):
    """real z = x: cosine transform = Re F_n(x)."""
    return F_n(n, acb(x)).real


def HF_real(N, x):
    return sum((cos_n_real(n, x) for n in range(1, N + 1)), arb(0))


def xi(s):
    return s * (s - 1) / 2 * ((-s / 2) * pi_().log()).exp() * (s / 2).gamma() * s.zeta()


def Hfun(z):
    return xi(acb(0.5) + acb(0, 1) * z / 2) / 8


def tail_bound(N, M, y):
    """rigorous bound for |sum_{n>M} cosine transform of phi_n| at |Im z| <= y (arb y)."""
    w1 = (5 + y) / 4          # w - 1 for the 9 part: (9+y)/4 - 1
    a1 = pi_() * (M + 1) ** 2
    assert float(a1) >= 2 * float(w1.upper()), "tail bound hypothesis a >= 2(w-1) fails"
    return 2 * (a1 + arb(3) / 2) * (-a1).exp()


def ball_zero(rad_arb):
    r = float(rad_arb.upper())
    return arb(0, max(r, 1e-300))


def E_N(N, z, M=None, y=None):
    """omitted terms n = N+1..M plus a rigorous tail ball; z acb (possibly a box)."""
    M = M or N + 5
    s = acb(0)
    for n in range(N + 1, M + 1):
        s += cos_n(n, z)
    if y is None:
        y = abs(z.imag).upper()
    t = tail_bound(N, M, y)
    return s + acb(ball_zero(t), ball_zero(t))



def crude_E(N, y, M=None):
    """rigorous bound for |E_N(z)| for every z with |Im z| <= y (arb y): int R_N(u) e^{y u} du with |phi_n| majorised by (2a^2 e^{9u} + 3a e^{5u}) e^{-a e^{4u}}
    (closed forms with the incomplete gamma function at real parameters) plus the rigorous tail for n > M."""
    M = M or N + 5
    w9 = (9 + y) / 4
    w5 = (5 + y) / 4
    tot = arb(0)
    for n in range(N + 1, M + 1):
        a = pi_() * n * n
        ac = acb(a)
        tot += a ** (2 - w9) * ac.gamma_upper(acb(w9)).real / 2 + a ** (1 - w5) * ac.gamma_upper(acb(w5)).real * 3 / 4
    return tot + tail_bound(N, M, y)


def E_N_Re(N, z, M=None):
    """what the 'Re' form of E:74 gives at complex z: Re of sum_{n>N} F_n(z)."""
    M = M or N + 5
    s = acb(0)
    for n in range(N + 1, M + 1):
        s += F_n(n, z)
    return s.real


# ---------- B_N ----------
def phi_derivs(a, u):
    """phi, phi', phi'' of phi(u) = (2a^2 e^{9u} - 3a e^{5u}) exp(-a e^{4u}) for an arb or acb u."""
    e4 = (4 * u).exp()
    e5 = (5 * u).exp()
    e9 = (9 * u).exp()
    q = (-a * e4).exp()
    p = 2 * a * a * e9 - 3 * a * e5
    p1 = 18 * a * a * e9 - 15 * a * e5
    p2 = 162 * a * a * e9 - 75 * a * e5
    q1 = -4 * a * e4 * q
    q2 = (16 * a * a * e4 * e4 - 16 * a * e4) * q
    return p * q, p1 * q + p * q1, p2 * q + 2 * p1 * q1 + p * q2


def Rpp(N, M, u):
    s = u * 0
    for n in range(N + 1, M + 1):
        s += phi_derivs(pi_() * n * n, u)[2]
    return s


def majorant_int(a, U):
    """int_U^inf of the positive majorant of |phi_n''| (rigorous): 32a^4 e^{17u} + 224 a^3 e^{13u} + 330 a^2 e^{9u} + 75 a e^{5u}, times e^{-a e^{4u}}."""
    tot = arb(0)
    A = acb(a * (4 * U).exp())
    for c, k, pw in ((32, 17, 4), (224, 13, 3), (330, 9, 2), (75, 5, 1)):
        # int_U^inf e^{ku} e^{-a e^{4u}} du = (1/4) a^{-k/4} Gamma(k/4, a e^{4U})
        g = A.gamma_upper(acb(arb(k) / 4)).real
        tot += c * a ** pw * g * a ** (-arb(k) / 4) / 4
    return tot


def L1_Rpp(N, M, Umax, ncells0, mindepth_w):
    """certified-sign integral of |R_N''| on [0, Umax] plus rigorous majorant beyond Umax and for n > M."""
    total = arb(0)
    undecided = arb(0)
    cells = [0]
    stack = [(arb(Umax) * i / ncells0, arb(Umax) * (i + 1) / ncells0) for i in range(ncells0)]
    ncell = 0
    while stack:
        a, b = stack.pop()
        ncell += 1
        ub = arb(0.5) * (a + b)
        ub = arb(ub.mid(), (b - a).mid() / 2 + 0)  # ball covering [a, b]
        f = Rpp(N, M, ub)
        if not (f.contains(0)):
            I = acb.integral(lambda x, _: Rpp(N, M, x), acb(a), acb(b)).real
            total += abs(I)
        elif float((b - a).mid()) < mindepth_w:
            undecided += (b - a) * abs(f).upper()
        else:
            m = (a + b) / 2
            stack.append((a, m))
            stack.append((m, b))
    tailU = arb(0)
    for n in range(N + 1, M + 1):
        tailU += majorant_int(pi_() * n * n, arb(Umax))
    a1 = pi_() * (M + 1) ** 2
    tailn = 2 * 350 * a1 ** 3 * (-a1).exp()
    return total, undecided, tailU + tailn, ncell


def Rp0(N, M):
    s = arb(0)
    for n in range(N + 1, M + 1):
        s += phi_derivs(pi_() * n * n, arb(0))[1]
    a1 = pi_() * (M + 1) ** 2
    tail = 2 * 100 * a1 ** 3 * (-a1).exp()
    return s, tail


# ---------- sign changes and zeros ----------
def sign_scan(fun, t0, t1, step, gam=None):
    """certified-sign grid scan of fun(x) (x real arb ball value) on t in [t0, t1]; returns list of (t_left, t_right) brackets and #uncertain."""
    brackets, unc = [], 0
    prev = None
    j = 0
    n = int(round((t1 - t0) / step))
    for j in range(n + 1):
        t = t0 + j * step
        v = fun(t)
        if v.contains(0):
            unc += 1
            continue
        sg = 1 if v > 0 else -1
        if prev is not None and sg != prev[1]:
            brackets.append((prev[0], t))
        prev = (t, sg)
    return brackets, unc


def refine(fun, lo, hi, iters=40):
    slo = fun(lo) > 0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        v = fun(mid)
        if v.contains(0):
            break
        if (v > 0) == slo:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run_approx(args, gam, out):
    res = {}
    step = 0.2 if args.quick else 0.05
    for N, tmax in ((2, 70), (3, 100), (4, 140)):
        if args.quick:
            tmax = min(tmax, 45)
        fun = lambda t, N=N: HF_real(N, 2 * arb(t))
        t_start = time.time()
        br, unc = sign_scan(fun, 0.1, tmax + 2, step)
        zs = [refine(fun, a, b, 30) for a, b in br]
        cnt = sum(1 for z in zs if z < tmax)
        ngam = Nof(tmax, gam)
        errs = [abs(zs[k] - gam[k]) if k < len(gam) else None for k in range(len(zs))]
        log(f"[approx] N={N} grid step {step}, t in (0.1, {tmax + 2}], uncertain {unc}, certified sign changes below t={tmax}: {cnt} (true zeta zeros below: {ngam}), time {time.time() - t_start:.1f}s")
        track = []
        for k, z in enumerate(zs):
            track.append((k + 1, z, errs[k]))
        log(f"[approx] N={N} first zeros (k, t_k, |t_k - gamma_k|): " + "; ".join(f"{k}:{z:.9f}:{e:.2e}" for k, z, e in track[:36]))
        for k in (3, 11, 25):
            if k <= len(zs):
                log(f"[approx] N={N} error at gamma_{k}: {errs[k - 1]:.3e}")
        # breakdown: first k with error > 1e-2, 1e-1, first missing sign change
        bd = {}
        for thr in (1e-3, 1e-2, 1e-1):
            kk = next((k for k, z, e in track if e is not None and e > thr), None)
            bd[str(thr)] = (kk, gam[kk - 1] if kk else None)
        first_deficit = None
        for k in range(len(zs)):
            if k < len(gam) and zs[k] - gam[k] > 0.5 or (k >= len(zs)):
                first_deficit = k + 1
                break
        # missing sign change: the first k with |t_k - gamma_k| > 0.5 (a skipped zero shifts every later index)
        log(f"[approx] N={N} first k with error above threshold (k, gamma_k): {bd}; first k with a shift above 0.5 (zero lost): {first_deficit}")
        res[f"N{N}"] = dict(step=step, tmax=tmax, certified_count_below_tmax=cnt, true_count=ngam, uncertain=unc, zeros=zs, errs=errs, thresholds=bd,
                            err_gamma3=errs[2] if len(errs) > 2 else None, err_gamma11=errs[10] if len(errs) > 10 else None,
                            err_gamma25=errs[24] if len(errs) > 24 else None)
    out["approx"] = res


def run_polya(args, gam, out):
    ctx.prec = 200
    pi = pi_()

    def HP(t):
        z = 2 * arb(t)
        zc = acb(z)
        k9 = acb(2 * pi).bessel_k((acb(9) + acb(0, 1) * zc) / 4)
        k5 = acb(2 * pi).bessel_k((acb(5) + acb(0, 1) * zc) / 4)
        return (pi * pi * k9 - 3 * pi / 2 * k5).real

    br, unc = sign_scan(HP, 0.1, 40 if not args.quick else 25, 0.1)
    zs = [refine(HP, a, b, 50) for a, b in br]
    log(f"[polya] sign changes of H^P(2t) in t (grid 0.1, certified signs, uncertain {unc}): first zeros in t: " + ", ".join(f"{z:.9f}" for z in zs[:8]) + "; z = 2t: " + ", ".join(f"{2 * z:.6f}" for z in zs[:3]))
    log(f"[polya] number of sign changes below t=40: {len(zs)}; zeta zeros below 40: {Nof(40, gam)} (the Polya kernel is not meant to track them; first zero t_1 vs gamma_1 = {gam[0]:.6f})")
    out["polya"] = dict(zeros_t=zs, first_zero_z=2 * zs[0])


def run_bn(args, gam, out, prec):
    ctx.prec = prec
    res = {}
    for N in (1, 2, 3, 4):
        M = N + 6
        p0, p0t = Rp0(N, M)
        tot, und, tl, nc = L1_Rpp(N, M, 1.0, 40, 1e-10)
        Bup = abs(p0).upper() + p0t.upper() + tot.upper() + und.upper() + tl.upper()
        kform = 8 * pi_() ** 3 * arb(N + 1) ** 6 * (-pi_() * (N + 1) ** 2).exp()
        log(f"[bn] N={N}: |R_N'(0)| = {abs(p0).str(12, radius=False)} (kink formula 8 pi^3 (N+1)^6 e^(-pi (N+1)^2) = {kform.str(6, radius=False)}), int|R''| = {tot.str(12, radius=False)} (undecided cells {und.str(3)}, tail bound {tl.str(3)}), cells {nc}")
        log(f"[bn] N={N}: B_N = {Bup.str(8, radius=False)}  B_N/|R'(0)| = {(Bup / abs(p0)).str(5, radius=False)}")
        res[f"N{N}"] = dict(Rp0=abs(p0).str(15, radius=False), L1=tot.str(15, radius=False), B=Bup.str(10, radius=False), kappa_formula=kform.str(8, radius=False))
    out["bn"] = res
    return res


def run_rouche(args, gam, out, bn):
    """modes: ident (identity check), signs (sign changes of H^F_N(2t) on (0,T)), sampled (points jlo..jhi-1 of the 721 point grid of the quarter circle,
    min |H|/|E_N| for the cosine form and for the 'Re' form, at two precisions), cert (interval cells), margin (real axis margin from B_N)."""
    ctx.prec = 300
    rows = [(1, 14.4), (1, 17.5), (2, 31.7), (2, 35.3), (3, 51.4), (3, 57.9), (4, 86.1), (4, 97.3)]
    if args.rows:
        rows = [(int(a.split(":")[0]), float(a.split(":")[1])) for a in args.rows.split(",")]
    mode = args.mode
    res = {}
    if mode == "ident":
        for z in (acb(0), acb(10.5), acb(30, 20), acb(60, -40), acb(150, 150)):
            d = Hfun(z) - sum((cos_n(n, z) for n in range(1, 16)), acb(0))
            log(f"[rouche] identity H(z) = Xi(z/2)/8 vs direct sum of 15 cosine transforms at z = {z.real.str(4, radius=False)} + {z.imag.str(4, radius=False)}i: |difference| = {abs(d).str(3)}")
        out["rouche"] = res
        return
    for N, T in rows:
        t0 = time.time()
        key = f"N{N}_T{T}"
        if mode == "signs":
            ngam = Nof(T, gam)
            lo, hi = gam[ngam - 1], gam[ngam]
            fun = lambda t, N=N: HF_real(N, 2 * arb(t))
            br, unc = sign_scan(fun, 0.05, T, 0.05)
            log(f"[rouche-signs] N={N} T={T}: zeta zeros below T: {ngam} (T between gamma_{ngam}={lo:.4f} and gamma_{ngam + 1}={hi:.4f}); certified sign changes of H^F_N(2t) on (0,T) at grid 0.05: {len(br)} (uncertain grid points {unc}); c = T/(N+1)^2 = {T / (N + 1) ** 2:.4f}")
            res[key] = dict(NT=ngam, sign_changes=len(br), uncertain=unc, c=T / (N + 1) ** 2)
        elif mode == "sampled":
            npts = 721
            def pt(pr, M, j):
                ctx.prec = pr
                th = arb.pi() / 2 * j / (npts - 1)
                z = acb(2 * T * th.cos(), 2 * T * th.sin())
                h = abs(Hfun(z))
                e = abs(E_N(N, z, M=M))
                rr = h / abs(E_N_Re(N, z, M=M))
                return h / e, rr
            best = None
            bestRe = None
            for j in range(args.jlo, args.jhi):
                r, rr = pt(250, N + 5, j)
                if best is None or r < best[0]:
                    best = (r, j)
                if bestRe is None or rr < bestRe[0]:
                    bestRe = (rr, j)
            # second precision / more terms on every 10th point of the chunk
            maxrel = arb(0)
            for j in range(args.jlo, args.jhi, 10):
                r1, rr1 = pt(250, N + 5, j)
                r2, rr2 = pt(450, N + 8, j)
                maxrel = max(maxrel, abs(r1 - r2) / r1, key=lambda v: float(v.mid())) if False else (maxrel if float(maxrel.mid()) > float((abs(r1 - r2) / r1).mid()) else abs(r1 - r2) / r1)
            ctx.prec = 300
            log(f"[rouche-sampled] N={N} T={T} points {args.jlo}..{args.jhi - 1} of 721: min |H|/|E_N| (cosine form) = {best[0].str(10, radius=False)} at index {best[1]}; 'Re' form min = {bestRe[0].str(10, radius=False)} at index {bestRe[1]}; max relative change of the ratio between (250 bits, N+5 terms) and (450 bits, N+8 terms) over every 10th point: {maxrel.str(3)}; time {time.time() - t0:.1f}s")
            res[key] = dict(jlo=args.jlo, jhi=args.jhi, min=best[0].str(12, radius=False), at=best[1], min_re=bestRe[0].str(12, radius=False), at_re=bestRe[1])
        elif mode == "cert":
            cells, worst, fail = 0, None, 0
            stack = [(arb(0), arb.pi() / 2)]
            while stack:
                a, b = stack.pop()
                cells += 1
                m = (a + b) / 2
                rad = (b - a) / 2
                th = arb(m.mid(), rad.mid() * 1.0000001 + 1e-60)
                z = acb(2 * T * th.cos(), 2 * T * th.sin())
                h = abs(Hfun(z)).lower()
                e = crude_E(N, abs(z.imag).upper()).upper()
                if h > e:
                    q = h / e
                    if worst is None or q < worst:
                        worst = q
                elif float(rad.mid()) > 1e-6:
                    stack.append((a, m))
                    stack.append((m, b))
                else:
                    fail += 1
            log(f"[rouche-cert] N={N} T={T}: CERTIFIED (|E_N| by the absolute integral majorant on each cell of theta, |H| by its box value) cells {cells}, uncertified cells {fail}, certified lower bound of |H|/|E_N| on the quarter circle: {worst.str(6, radius=False) if worst is not None else None}; time {time.time() - t0:.1f}s")
            res[key] = dict(cells=cells, fail=fail, lower=worst.str(8, radius=False) if worst is not None else None)
        elif mode == "margin":
            Bup = arb(bn[f"N{N}"]["B"])
            marg = abs(Hfun(acb(2 * T))).real.lower() * (2 * T) ** 2 / Bup
            hv = Hfun(acb(2 * T)).real
            log(f"[rouche-margin] N={N} T={T}: H(2T) = {hv.str(8)}; margin |H(2T)|(2T)^2/B_N = {marg.str(6, radius=False)}")
            res[key] = dict(margin=marg.str(8, radius=False))
    out["rouche_" + mode] = res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--part", default="all")
    ap.add_argument("--zeros", default=None)
    ap.add_argument("--rows", default=None)
    ap.add_argument("--mode", default="ident")
    ap.add_argument("--jlo", type=int, default=0)
    ap.add_argument("--jhi", type=int, default=721)
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    zp = args.zeros or next(p for p in (os.path.join(here, "..", "data", "zeros_1700.txt"), "zeros_1700.txt") if os.path.exists(p))
    gam = load_zeros(zp)
    ctx.prec = 300
    out = {"quick": args.quick}
    part = args.part
    bn = None
    if part in ("bn", "all") or (part == "rouche" and args.mode == "margin"):
        bn = run_bn(args, gam, out, 300)
        bn2 = run_bn(args, gam, {}, 450)
        for k in bn:
            d = abs(arb(bn[k]["B"]) - arb(bn2[k]["B"])) / arb(bn[k]["B"])
            log(f"[bn] precision check {k}: relative difference of B_N between 300 and 450 bits (and with the same cell rule): {d.str(3)}")
        ctx.prec = 300
    if part in ("approx", "all"):
        run_approx(args, gam, out)
    if part in ("polya", "all"):
        run_polya(args, gam, out)
    if part in ("rouche", "all"):
        run_rouche(args, gam, out, bn)
    out["seconds"] = time.time() - START
    outp = args.out or (("/tmp/ev_g1_fold_quick.json") if args.quick else os.path.join(here, "..", "data", "ev_g1_fold.json"))
    json.dump(out, open(outp, "w"), indent=1, default=str)
    open(os.path.splitext(outp)[0] + ".log", "w").write("\n".join(LOG) + "\n")
    log(f"done in {time.time() - START:.1f}s; wrote {outp}")


if __name__ == "__main__":
    main()
