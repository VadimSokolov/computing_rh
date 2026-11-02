"""Ball values accompanying the exact counterexamples proved in report.md."""
from pathlib import Path
import sys
import json

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'pylib'))
from flint import arb, acb, ctx
ctx.prec = 200
h = arb(1)/10
z0 = acb(2*h, arb.pi()/2)

def M(z):
    return z.cosh()*(z-2*h).cosh()*(z+2*h).cosh()/(2*h).cosh()**2

def A(z):
    return z.tanh()+(z-2*h).tanh()+(z+2*h).tanh()

out = {
    'example_h': '1/10',
    'tested_level': '3/5',
    'right_zero_ball': (z0+arb(1)/2).str(45),
    'M_at_explicit_zero_ball': M(z0).str(45),
    'global_boundary_density_lower_bound_ball': ((3*h).tanh()/arb.pi()).str(45),
    'boundary_values': []
}
for w in [arb(0), arb(1), arb.pi()/2, arb(10)]:
    actual = A(acb(h,w)).real/arb.pi()
    expected = (6*h).sinh()/((6*h).cosh()+(2*w).cos())/arb.pi()
    assert (actual-expected).contains(0)
    assert actual > 0
    out['boundary_values'].append({'w':w.str(30),'actual_ball':actual.str(45),'formula_ball':expected.str(45)})
n = json.loads((ROOT/'results/numbers.json').read_text())
c = n['dh_root']
rho = acb(arb(c['center_real'], c['radius']), arb(c['center_imag'],c['radius']))
pole = (rho-arb(3)/4)**2
assert pole.imag > 0
out['DH_level_3over4_pole_ball'] = pole.str(35)
out['DH_principal_sqrt_real_ball'] = pole.sqrt().real.str(35)
(ROOT/'results/counterexamples.json').write_text(json.dumps(out,indent=2)+'\n')
ledger = ROOT/'results/numbers.txt'
previous = ledger.read_text().split('\nexact_counterexample_ball_values\n', 1)[0]
ledger.write_text(previous+'\nexact_counterexample_ball_values\n'+json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
