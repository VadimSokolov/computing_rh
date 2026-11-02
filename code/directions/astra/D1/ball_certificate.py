"""Proposed ball verification of the first prime rung.

Not run in this review: python-flint is unavailable locally.
Requires python-flint 0.9 API. Run with a bounded timeout.
All quadrature and operator errors are proved in report.md.
"""
import json
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / '.deps'))
try:
    from flint import arb, acb, ctx
except ImportError:
    raise SystemExit('NOT RUN: install python-flint before running this verifier.')

ctx.prec = 256
L = arb(3).log()/2
T = 60
N = 64
QORDER = 40
PANELS = 240
pi = arb.pi()
c = arb(2).sqrt()*arb(2).log()
beta = (arb(T)/(2*pi)).log()-arb(1)/(3*T)-arb(1)/(2*T*T)-c
target = arb('2e-8')
quadrature_entry_error = arb(8)*T*20000 / (arb(2)**(2*QORDER))
quadrature_norm_error = N*quadrature_entry_error
nodes = [arb.legendre_p_root(QORDER, j, weight=True) for j in range(QORDER)]
scales = [(2*L*(2*n+1)).sqrt() for n in range(N)]
x = L/2
pole = [scales[n]*(pi/(2*x)).sqrt()*x.bessel_i(arb(n)+arb(1)/2) for n in range(N)]
blocks = {}
for parity in (0,1):
    ids = list(range(parity,N,2))
    blocks[parity] = [[(beta if i==j else arb(0))
                       +(2 if parity==0 else -2)*pole[n]*pole[m]
                       for j,m in enumerate(ids)] for i,n in enumerate(ids)]

for panel in range(PANELS):
    mid = (arb(panel)+arb(1)/2)/4
    half = arb(1)/8
    for node, weight in nodes:
        t = mid+half*node
        x = L*t
        vals = [scales[n]*(pi/(2*x)).sqrt()*x.bessel_j(arb(n)+arb(1)/2)
                *(-1 if (n//2)%2 else 1) for n in range(N)]
        symbol = acb(arb(1)/4,t/2).digamma().real-pi.log()-c*(t*arb(2).log()).cos()-beta
        weight = half*weight*symbol/pi
        for parity in (0,1):
            ids = list(range(parity,N,2))
            M = blocks[parity]
            for i,n in enumerate(ids):
                v = weight*vals[n]
                for j in range(i+1):
                    M[i][j] += v*vals[ids[j]]
    if panel % 40 == 39:
        print(f'Assembled {panel+1}/{PANELS} panels', flush=True)


def verify_ldl(M):
    n = len(M)
    lower = [[arb(0) for _ in range(n)] for _ in range(n)]
    d = []
    for k in range(n):
        pivot = M[k][k]-target-quadrature_norm_error
        for j in range(k):
            pivot -= lower[k][j]**2*d[j]
        if not pivot > 0:
            raise ArithmeticError(f'Unproved positive pivot {k}: {pivot}')
        d.append(pivot)
        for i in range(k+1,n):
            v = M[i][k]
            for j in range(k):
                v -= lower[i][j]*lower[k][j]*d[j]
            lower[i][k] = v/pivot
    return [str(x) for x in d]


double_factorial = prod(range(1,2*N+2,2))
r = (L*T)**2/((2*N+1)*(2*N+3))
S = 2*L*(2*N+1)*(L*T)**(2*N)/(arb(double_factorial)**2*(1-r))
rp = (L/2)**2/((2*N+1)*(2*N+3))
Sp = 2*L*(2*N+1)*L.exp()*(L/2)**(2*N)/(arb(double_factorial)**2*(1-rp))
kappa = arb(8)*T/pi*(2*(2*L*S).sqrt()+S)
kappa += 2*(2*(2*L.sinh()*Sp).sqrt()+Sp)
assert 0 < beta and r < 1 and rp < 1
assert kappa < arb('3e-9')
pivots = {str(p):verify_ldl(blocks[p]) for p in (0,1)}
assert target-kappa > arb('1.7e-8')
result = {'status':'ball certificate passed', 'precision_bits':ctx.prec,
          'L':str(L), 'T':T, 'orders':N, 'beta':str(beta),
          'quadrature_entry_error':str(quadrature_entry_error),
          'quadrature_norm_error':str(quadrature_norm_error),
          'operator_tail':str(kappa), 'proved_lower_bound':'1.7e-8',
          'ldl_pivots':pivots}
Path('ball_certificate_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='ldl_pivots'},indent=2))
