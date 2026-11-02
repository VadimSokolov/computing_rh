"""Independent D2 calculations. Reads book data only; uses no book scripts.

Run: timeout 30 python3 -I check_margin.py
All newly calculated quantities are floating point, not certificates.
"""
import json
import math
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
BOOK = Path('.')
ZERO_PATH = BOOK / 'data/zeros_1700.txt'
FILES = ['wl_pcert.json', 'wl_grid.json', 'wl_grid2.json', 'wl_grid3.json']
zeros = [float(s) for s in ZERO_PATH.read_text().split()]
records = {}
for filename in FILES:
    for row in json.loads((BOOK / 'data' / filename).read_text()):
        if row.get('certified_positive') and 'eps_up_log10' in row:
            key = (row['X'], row['N'])
            records[key] = dict(row, source=filename)


def logsumexp(values):
    values = list(values)
    m = max(values)
    return m + math.log(math.fsum(math.exp(v-m) for v in values))


def tail_bound(T, count):
    """HSW Corollary 1.2 integrated by parts, proved in report.md."""
    u = math.log(T)
    return (math.log(T/(2*math.pi))/(math.pi*T)
            + (0.1038*(u+0.5)+0.2573*(math.log(u)+0.5/u)
               + 9.3675-count)/(T*T))


def interpolation(X, n):
    L = math.log(X)/2
    h = math.pi/L
    g = zeros[:n]
    q = max(h*j/z for j, z in enumerate(g, 1))
    if q >= 1:
        return dict(q=q, valid=False)
    log_r0 = 2*math.fsum(math.log(z) for z in g)-2*n*math.log(h)-2*math.lgamma(n+1)
    logr = [math.log(2)+2*log_r0]
    for j in range(1, n+1):
        w2 = (h*j)**2
        r = (math.fsum(math.log(abs(w2-z*z)) for z in g)
             - 2*n*math.log(h)-math.lgamma(n-j+1)
             - math.lgamma(n+j+1)+math.log(2))
        logr.append(2*r)
    log_norm = math.log(L)+logsumexp(logr)
    T, K = zeros[-1], len(zeros)
    remaining = tail_bound(T, K)
    total_tail = math.fsum(z**-2 for z in zeros[n:])+remaining
    log_product = math.log(4/L)-2*log_r0+math.log(total_tail)
    log_residue = math.log(8)-log_norm+math.log(total_tail)
    sample_logs = []
    for t in zeros[n:]:
        log_ratios = math.fsum(math.log((t*t-z*z)/(t*t-(h*j)**2))
                              for j, z in enumerate(g, 1))
        sample_logs.append(math.log(8)-log_norm-2*math.log(t)
                           +2*math.log(abs(math.sin(L*t)))+2*log_ratios)
    sample_log = logsumexp(sample_logs)
    upper_log = logsumexp([sample_log, math.log(8)-log_norm+math.log(remaining)])
    return dict(valid=True, q=q, c_q=4*math.log(1/q),
                product_bound_log10=log_product/math.log(10),
                residue_bound_log10=log_residue/math.log(10),
                sampled_trial_log10=sample_log/math.log(10),
                sampled_plus_tail_bound_log10=upper_log/math.log(10))


tasks = [(2, 8), (2, 100), (2, 200), (2, 400), (2, 600),
         (13, 8), (13, 12), (13, 20), (13, 40), (13, 100),
         (13, 200), (13, 400), (13, 600),
         (50, 20), (50, 40), (50, 100), (50, 200), (50, 400),
         (1000, 20), (1000, 30), (1000, 40), (1000, 100),
         (1000, 200), (1000, 400), (1000, 600),
         (1000000, 8), (1000000, 100)]
rows = []
for X, n in tasks:
    if (X, n) not in records:
        continue
    r = records[(X, n)]
    row = dict(X=X, N=n, book_log10=r['eps_up_log10'],
               book_upper=r['eps_up'], book_lower=r['eps_lo'],
               source=r['source'], precision_digits=r['D'],
               interpolation=interpolation(X, n))
    rows.append(row)

# Floating-point checks of exact algebra, independent of the grid matrices.
X, n = 13, 8
L, h = math.log(X)/2, 2*math.pi/math.log(X)
g = zeros[:n]
r0 = math.prod((z/(h*j))**2 for j,z in enumerate(g, 1))
rs = [r0]
for j in range(1,n+1):
    w2 = (h*j)**2
    numerator = math.prod(w2-z*z for z in g)
    denominator = w2*math.prod(w2-(h*k)**2 for k in range(1,n+1) if k != j)
    rs.append(numerator/denominator)
algebra_errors = []
for t in [3.7, 11.2, 27.8, 80.1]:
    rational = math.prod((t*t-z*z)/(t*t-(h*j)**2) for j,z in enumerate(g,1))/t**2
    partial = math.fsum(r/(t*t-(h*j)**2) for j,r in enumerate(rs))
    algebra_errors.append(abs(rational-partial)/max(1,abs(rational)))

mp.mp.dps = 90
lm = mp.log(13)/2
hm = mp.pi/lm
gm = [mp.mpf(s) for s in ZERO_PATH.read_text().split()[:8]]
rm = [mp.fprod((z/(hm*j))**2 for j,z in enumerate(gm, 1))]
for j in range(1, 9):
    w2 = (hm*j)**2
    rm.append(mp.fprod(w2-z*z for z in gm)
              / (w2*mp.fprod(w2-(hm*k)**2 for k in range(1,9) if k != j)))
mp_errors = []
for ts in ['3.7','11.2','27.8','80.1']:
    t = mp.mpf(ts)
    rational = mp.fprod((t*t-z*z)/(t*t-(hm*j)**2)
                        for j,z in enumerate(gm,1))/t**2
    partial = mp.fsum(r/(t*t-(hm*j)**2) for j,r in enumerate(rm))
    mp_errors.append(abs(rational-partial)/max(1,abs(rational)))
assert max(mp_errors) < mp.mpf('1e-80')

eligible = [r for (X,n),r in records.items() if n <= X]
weakest = min(eligible, key=lambda r: -math.log(10)*r['eps_up_log10']/r['N'])
finite_c = dict(c=-math.log(10)*weakest['eps_up_log10']/weakest['N'],
                X=weakest['X'], N=weakest['N'], cells=len(eligible))

scaled_boundary = []
for exponent in [10, 100, 1000, 10000]:
    # Asymptotic diagnostic only. gamma_N solved from the RVM MAIN TERM.
    logn = float(exponent)
    # W(N/e) is v with logn-1 on the right.
    v = logn-math.log(logn)-1
    for _ in range(20):
        v -= (v+math.log(v)-(logn-1))/(1+1/v)
    qmodel = v/logn
    scaled_boundary.append(dict(logN=exponent, q_model=qmodel,
                                c_model=4*math.log(1/qmodel)))

result = dict(status='new calculations are floating point, not certified',
              grid_cells=len(records), zeros_used=len(zeros),
              zero_height=zeros[-1], tail_bound=tail_bound(zeros[-1],len(zeros)),
              partial_fraction_max_relative_error=max(algebra_errors),
              partial_fraction_mp_error=mp.nstr(max(mp_errors),8),
              finite_grid_exponential=finite_c,
              crossover_scales={str(X): X*math.log(X) for X in (13, 1000)},
              rows=rows, boundary_main_term_model=scaled_boundary)
(HERE/'comparison.json').write_text(json.dumps(result, indent=2)+'\n')
lines = ['X N book_log10 product_bound_log10 q sampled_plus_tail_bound_log10']
for r in rows:
    t = r['interpolation']
    lines.append(f"{r['X']} {r['N']} {r['book_log10']:.8f} "
                 + (f"{t['product_bound_log10']:.8f} {t['q']:.6f} "
                    f"{t['sampled_plus_tail_bound_log10']:.8f}" if t['valid']
                    else f"not_applicable {t['q']:.6f} not_applicable"))
lines += ['','Boundary main-term diagnostic (not exact zeros):']
lines += [str(r) for r in scaled_boundary]
lines += [f"Tail beyond {len(zeros)} ordinates: {result['tail_bound']:.12g}",
          f"Partial fraction relative residual: {max(algebra_errors):.3g}",
          f"90-digit partial fraction residual: {mp.nstr(max(mp_errors),8)}",
          f"Finite grid exponential diagnostic: {finite_c}"]
(HERE/'numbers.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
