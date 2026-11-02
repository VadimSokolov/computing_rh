"""Targeted checks of numerical discrepancies, using only independently saved data."""
import bisect
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=60
OUT=Path(__file__).resolve().parent
zeros=[mp.mpf(x) for x in json.loads((OUT/'independent_zeta_zeros.json').read_text())]
bessel=json.loads((OUT/'independent_bessel_zeros.json').read_text())
data={'precision_dps':60,'status':'computed, not certified'}


def interval_counts(ref):
    ref=[r for r in ref if r<zeros[-1]]
    count=[bisect.bisect_left(zeros,y)-bisect.bisect_right(zeros,x) for x,y in zip(ref[:-1],ref[1:])]
    bad=[j for j,n in enumerate(count) if n!=1]
    return {'intervals':len(count),'failures':len(bad),'first_failure':str(ref[bad[0]])}


data['ideal_phase_references']={str(c):interval_counts([mp.grampoint(n+c) for n in range(405)])
                                for c in [mp.mpf('.5'),mp.mpf('.625')]}
data['actual_first_roots_above_height10']={k:[str(2*mp.mpf(x)) for x in v[:4]] for k,v in bessel.items()}
two_term=lambda nu:(2*mp.pi**2*mp.re(mp.besselk(mp.mpf(9)/4+1j*nu,2*mp.pi))
                    -3*mp.pi*mp.re(mp.besselk(mp.mpf(5)/4+1j*nu,2*mp.pi)))
data['two_term_lowest_height']=str(2*mp.findroot(two_term,(mp.mpf('3.7'),mp.mpf('3.9'))))
data['base_refined_roots']={}
for index in [0,1,2,489]:
    start=mp.mpf(bessel['base'][index])
    fun=lambda x:mp.re(mp.besselk(1j*x,2*mp.pi))*mp.exp(mp.pi*x/2)
    root=mp.findroot(fun,(start-mp.mpf('.001'),start+mp.mpf('.001')))
    data['base_refined_roots'][str(index+1)]={'root':str(root),'change':str(root-start),'scaled_residual':str(abs(fun(root)))}
t=mp.mpf('8.039316023094285')
data['screw_counterexample']={'t':str(t),'value_60dps':str(sum(2*(1-mp.cos(g*t))/g**2 for g in zeros))}
(OUT/'verification_numbers.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
