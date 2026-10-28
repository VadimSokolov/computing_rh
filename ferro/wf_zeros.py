# Note W (Chapter 13, Section fm:sec:bessel, paragraph "A small distance does not pin the zeros"), computations N3 and N4.
# N3: the first twenty positive zeros of the characteristic function phi_q of the laws of rows 2 and 4 of Table fm:tab:bessel
# (exact decimals of wf_fit.json, book units). By Theorem fm:thm:bessel(ii), phi_q(theta) is proportional to
#   F(theta) = sum_{eta in {+-1}^N} K_{(i theta + sum_k eta_k lam_k)/l}(a)     (mu = 0; real for real theta),
# with K the Bessel function of the second kind, evaluated in ball arithmetic (python-flint acb.bessel_k, complex order,
# working precision 256 bits, doubled until the sign or the ball is decided, up to 2048).
# Method. (1) Scan [0.25, 100] with step 1/20; a sign change is certified when the balls of F at both ends exclude 0
# with opposite signs; each sign change is bisected 34 times (balls at the midpoints exclude 0), width below 1e-10, and the
# midpoint is reported. (2) Theta = (twenty first zero) + 1/4 (checked: F(Theta) excludes 0, and the twenty second zero is
# above Theta + 1/4). (3) Completeness and multiplicity: the number of zeros of the entire function F in the rectangle
# (1/4, Theta) x (-2, 2), counted by the argument principle on its boundary.
# NOT RIGOROUS: a ball enclosure of F on a segment is far too wide (Arb loses a factor e^{pi nu} on a ball order, F being
# about 1e-25 near theta = 80), so the half plane criterion is not used; the count is the unwrapped phase of F at samples
# (values certified to 1e-12 relative, step 1/40 and again 1/80, largest increment of the phase below 1 radian),
# as in fm_polya.py. N must equal the number of certified
# sign changes in (1/4, Theta): then every zero in the strip is real and simple. (4) gamma_1..gamma_20 from Arb's certified
# zeros (acb.zeta_zeros, 128 bits); mean and largest |theta_n - gamma_n|, n = 1..20, from the unrounded midpoints.
# N4: theta* = the least theta_0 such that |phi_X(theta)| = |xi(1/2 + i theta)|/xi(1/2) < 2 TV_2, TV_2 the certified bound
# of row 2 (wf_tvcert.json, upper end), for every theta in [theta_0, 100]: mpmath at 30 digits, grid step 1e-3 downward from
# 100 to the first point where |phi_X| >= 2 TV_2, then bisection of the crossing to 1e-12. Computed, not certified.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 wf_zeros.py [--quick] [rowname ...]  reads wf_fit.json, wf_tvcert.json; writes wf_zeros.json
import json, math, os, sys, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wf_tvcert as W
from flint import arb, acb, ctx

QUICK = '--quick' in sys.argv
CPUS = int(os.environ.get('SLURM_CPUS_PER_TASK', '4'))
PRECS = (256, 512, 1024, 2048)
SPEC = {k: v['exact'] for k, v in json.load(open('wf_fit.json')).items() if isinstance(v, dict) and 'exact' in v}
_cache = {}


def member(row, prec):
    key = (row, prec)
    if key not in _cache:
        ctx.prec = prec
        _cache[key] = W.Member.from_row(row, SPEC[row])
    return _cache[key]


def ysum(m):
    out = []
    for y in m.ys():
        if y.contains(0):
            assert y.rad() < 1e-30
            y = arb(0)
        out.append(y)
    return out


def Fz(row, z, prec):
    """F at an acb z (a point or a box) at the working precision prec."""
    m = member(row, prec)
    ctx.prec = prec
    z = acb(z)
    tot = acb(0)
    for y in ysum(m):
        tot += acb(m.a).bessel_k((acb(0, 1)*z + y)/m.l)
    return tot


def sign(row, th):
    """+1, -1 (certified sign of the real function F at the real ball th), or 0 if undecided up to 2048 bits."""
    for p in PRECS:
        ctx.prec = p
        v = Fz(row, acb(th), p).real
        if v > 0:
            return 1
        if v < 0:
            return -1
    return 0


def scan_chunk(args):
    row, i0, i1 = args
    ctx.prec = 256
    return [sign(row, arb(i)/20) for i in range(i0, i1)]


def bisect(args):
    row, i = args
    ctx.prec = 256
    lo, hi = arb(i)/20, arb(i + 1)/20
    slo = sign(row, lo)
    for _ in range(34):
        mid = (lo + hi)/2
        s = sign(row, mid)
        assert s != 0
        if s == slo:
            lo = mid
        else:
            hi = mid
    return lo.mid().str(20, radius=False), hi.mid().str(20, radius=False)


def ap_chunk(args):
    """Sum of the principal increments arg(F(z_{j+1})/F(z_j)) over the samples z_j = (pa + (pb - pa) j/n)/20 (units of 1/20)."""
    row, pa, pb, n = args
    p = 512
    ctx.prec = p
    def pt(j):
        return acb(arb(pa[0]) + (arb(pb[0]) - arb(pa[0]))*j/n, arb(pa[1]) + (arb(pb[1]) - arb(pa[1]))*j/n)/20
    vals = []
    for j in range(n + 1):
        v = Fz(row, pt(j), p)
        ctx.prec = p
        assert float(max(v.real.rad(), v.imag.rad())) < 1e-12*float(abs(v).mid()), 'inaccurate value'
        vals.append(v)
    tot = arb(0); big = 0.0
    for j in range(n):
        d = (vals[j + 1]/vals[j]).arg()
        big = max(big, abs(float(d.mid())))
        tot += d
    return tot.str(30), big


def run_row(row, pool):
    t0 = time.time()
    top = 100*20 if not QUICK else 60*20
    chunks = [(row, i, min(i + 40, top)) for i in range(5, top, 40)]       # i/20 from 0.25
    signs = [s for c in pool.map(scan_chunk, chunks) for s in c]
    assert 0 not in signs, 'undecided sign on the scan'
    idx = [5 + k for k in range(len(signs) - 1) if signs[k] != signs[k + 1]]
    zer = pool.map(bisect, [(row, i) for i in idx])
    zmid = [(arb(a) + arb(b))/2 for a, b in zer]
    nz = 21 if not QUICK else 5
    Theta20 = arb(zmid[nz - 1]) + arb(1)/4
    assert arb(zmid[nz]) > Theta20 + arb(1)/4, 'zeros %d and %d too close' % (nz, nz + 1)
    s_th = sign(row, Theta20)
    assert s_th != 0
    ThetaN = int((Theta20*20).floor().mid().str(5, radius=False).split('.')[0]) + 1       # Theta on the grid 1/20 above that
    Theta = arb(ThetaN)/20
    assert arb(zmid[nz - 1]) < Theta < arb(zmid[nz])
    assert sign(row, Theta) != 0
    inside = sum(1 for z in zmid if z < Theta)
    n_sign = inside
    # argument principle on (1/4, Theta) x (-2, 2) by the unwrapped phase of F at samples (corners in units of 1/20), with
    # the step 1/40 and again 1/80 (the boxes of a ball enclosure are too wide: see the header)
    L, R, Bt, Tp = 5, ThetaN, -40, 40
    counts = {}
    for per in (2, 4):                                                          # samples per 1/20
        jobs = []
        for (pa, pb) in (((L, Bt), (R, Bt)), ((R, Bt), (R, Tp)), ((R, Tp), (L, Tp)), ((L, Tp), (L, Bt))):
            length = max(abs(pb[0] - pa[0]), abs(pb[1] - pa[1]))
            nchunk = max(1, length//20)
            for c in range(nchunk):
                a = (pa[0] + (pb[0] - pa[0])*c/nchunk, pa[1] + (pb[1] - pa[1])*c/nchunk)
                b = (pa[0] + (pb[0] - pa[0])*(c + 1)/nchunk, pa[1] + (pb[1] - pa[1])*(c + 1)/nchunk)
                jobs.append((row, a, b, int(length/nchunk*per) + 1))
        res = pool.map(ap_chunk, jobs)
        ctx.prec = 256
        total = sum((arb(r[0]) for r in res), arb(0))
        N = total/(2*arb.pi())
        counts[per] = (N, max(r[1] for r in res))
    Ns = {per: (int(round(float(N.mid()))), float(abs(N.mid() - round(float(N.mid())))), mx) for per, (N, mx) in counts.items()}
    ok = len(set(v[0] for v in Ns.values())) == 1 and all(v[1] < 0.01 and v[2] < 1.0 for v in Ns.values())
    Nint = Ns[2][0]
    N = counts[4][0]
    print('%s: %d sign changes in (0.25, %s), argument principle count %s (steps 1/40 and 1/80: %s; %s)'
          % (row, n_sign, Theta.str(8, radius=False), N.str(8), Ns, 'consistent' if ok and Nint == n_sign else 'MISMATCH'), flush=True)
    NM = 20 if not QUICK else 4
    ctx.prec = 128
    from flint import acb as A
    G = [z.imag for z in A.zeta_zeros(1, NM)]
    th = zmid[:NM]
    diffs = [abs(t - g) for t, g in zip(th, G)]
    mean = sum(diffs, arb(0))/NM
    mx = max(diffs, key=lambda x: float(x.mid()))
    out = dict(row=row, spec=SPEC[row], Theta=Theta.str(10, radius=False), n_sign_changes=n_sign, ap_count=N.str(10),
               ap_ok=bool(ok and Nint == n_sign), ap_detail={str(k): list(v) for k, v in Ns.items()}, zeros=[z.str(15, radius=False) for z in zmid[:nz]],
               zero_width=[float((arb(b) - arb(a)).mid()) for a, b in zer[:nz]],
               gamma=[g.str(15, radius=False) for g in G], first_zero=th[0].str(10, radius=False),
               mean_abs_diff=mean.str(10, radius=False), max_abs_diff=mx.str(10, radius=False),
               argmax=int(max(range(NM), key=lambda k: float(diffs[k].mid()))) + 1, secs=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    return out


def n4(tv2):
    import mpmath as mp
    mp.mp.dps = 30
    half = mp.mpf(1)/2
    xi = lambda s: s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    X0 = xi(half).real
    phi = lambda th: abs(xi(half + 1j*th))/X0
    thr = 2*mp.mpf(tv2)
    step = mp.mpf(1)/1000
    th = mp.mpf(100) if not QUICK else mp.mpf(40)
    while phi(th) < thr:
        th -= step
        assert th > 0
    hi_, lo_ = th + step, th                       # |phi| < thr at hi_ (checked below), >= thr at lo_
    assert phi(hi_) < thr
    for _ in range(45):
        mid = (hi_ + lo_)/2
        if phi(mid) >= thr:
            lo_ = mid
        else:
            hi_ = mid
    ths = (hi_ + lo_)/2
    return dict(threshold=mp.nstr(thr, 10), theta_star=mp.nstr(ths, 12), phi_at=mp.nstr(phi(ths), 10),
                phi_above_max_check='grid 1e-3 from 100 down; first point with |phi| >= threshold at %s' % mp.nstr(lo_, 8))


if __name__ == '__main__':
    t0 = time.time()
    rows = [a for a in sys.argv[1:] if not a.startswith('--')] or ['row2', 'row4']
    ctx.prec = 256
    R = {}
    with Pool(CPUS) as pool:
        for r in rows:
            R[r] = run_row(r, pool)
    tvj = json.load(open('wf_tvcert.json')) if os.path.exists('wf_tvcert.json') else None
    if tvj:
        tv2 = [x for x in tvj['results'] if x['name'] == 'row2'][0]['TV'][1]
        R['N4'] = n4(tv2)
        R['N4']['tv2_certified_upper'] = tv2
        print('N4:', json.dumps(R['N4']), flush=True)
    R['secs'] = round(time.time() - t0, 1)
    json.dump(R, open('/tmp/wf_zeros_quick.json' if QUICK else 'wf_zeros.json', 'w'), indent=1)
    print('done in %.0f s' % R['secs'])
