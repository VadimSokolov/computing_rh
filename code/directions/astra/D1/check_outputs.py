"""Check local deliverables and record exact arithmetic tail checks."""
import ast
from fractions import Fraction as F
from math import prod
from pathlib import Path
import json

root = Path(__file__).resolve().parent
scripts = ['pilot.py', 'zero_tail_pilot.py', 'ball_certificate.py',
           'probe_env.py', 'install_flint.py', 'check_outputs.py']
for name in scripts:
    ast.parse((root/name).read_text(), filename=name)
for name in scripts+['report.md', 'empirical.md']:
    txt = (root/name).read_text()
    assert chr(8212) not in txt and chr(45)*3 not in txt, name


def exp_bounds(x, order=40):
    term = F(1)
    total = term
    for k in range(1,order+1):
        term *= x/k
        total += term
    nxt = term*x/(order+1)
    return total, total+nxt/(1-x/(order+2))


lo = F(5493061443340,10**13)
hi = F(5493061443341,10**13)
assert exp_bounds(2*lo)[1] < 3 < exp_bounds(2*hi)[0]
N,T = 64,60
df = prod(range(1,2*N+2,2))
r = (hi*T)**2/((2*N+1)*(2*N+3))
S = 2*hi*(2*N+1)*(hi*T)**(2*N)/(df**2*(1-r))
rp = (hi/2)**2/((2*N+1)*(2*N+3))
Sp = 4*hi*(2*N+1)*(hi/2)**(2*N)/(df**2*(1-rp))
assert r < 1 and rp < 1
assert S < F(67,10**24)
sqrt_bound = F(86,10**13)
assert 2*hi*S < sqrt_bound**2
assert Sp < F(1,10**250)
# pi>3, 2*sinh(L)<2, and sqrt(2*Sp)<2e-125 suffice here.
kappa_upper = F(8*T,3)*(2*sqrt_bound+F(67,10**24))
kappa_upper += 2*(4*F(1,10**125)+F(1,10**250))
assert kappa_upper < F(2754,10**12) < F(3,10**9)
quadrature_norm = F(N*8*T*20000,2**80)
assert quadrature_norm < F(6,10**16)
records = {
    'geometric_pilot': json.loads((root/'pilot_results.json').read_text()),
    'zero_tail_pilot': json.loads((root/'zero_tail_results.json').read_text()),
    'exact_rational_checks': {
        'L_lower':str(lo), 'L_upper':str(hi), 'tail_S_upper':'6.7e-23',
        'pole_tail_upper':'1e-250', 'operator_tail_upper':'2.754e-9',
        'quadrature_norm_upper':'6e-16', 'all_assertions_passed':True,
        'note':'Exact arithmetic checks of analytic error allowances, not a ball matrix certificate.'
    }
}
(root/'results').mkdir(exist_ok=True)
(root/'results'/'numbers.txt').write_text(json.dumps(records,indent=2)+'\n')
print('PASS: script syntax, output punctuation, exact rational tail allowances, numerical records.')
print('NOT RUN: ball matrix factorization, because python-flint is unavailable.')
