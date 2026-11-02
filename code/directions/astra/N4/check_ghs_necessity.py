"""Independent countercheck of a necessity inference, no external code read."""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps=60
sig=mp.mpf('.2')
def density(x):
    return sum(p*mp.exp(-(x-a)**2/(2*sig**2)) for a,p in [(-2,mp.mpf('.25')),(0,mp.mpf('.5')),(2,mp.mpf('.25'))])
def U(x):
    return -mp.log(density(x))
values={str(x):str(mp.diff(U,mp.mpf(x),3)) for x in ['.5','.8','.9','1','1.1','1.5']}
res={'law':'two independent fair signs plus Gaussian standard deviation .2',
     'Uthird':values,
     'mgf':'exp(.02*z^2)*cosh(z)^2',
     'Kthird_at_1':str(-4/mp.cosh(1)**2*mp.tanh(1))}

# Convolution of two quartic Ising-reachable spins, with all Gaussian moments.
quartic=[]
for a in [20,50]:
    s=mp.mpf('.5')
    c=4-3*s*s
    peak=mp.sqrt(c/4)
    def integral(k):
        return mp.quad(lambda t:t**(2*k)*mp.exp(-2*a*(t*t-c/4)**2),[0,peak,1,2,mp.inf])
    norm=integral(0)
    m=[integral(k)/norm for k in range(1,4)]
    var=m[1]-m[0]**2
    k3=m[2]-3*m[0]*m[1]+2*m[0]**3
    third=3*a*s+216*s**3*a**3*k3-108*s*a*a*var
    quartic.append({'A':a,'s':str(s),'Uthird':str(third),'leading_asymptotic':str(-24*a*s)})
res['sum_of_two_quartic_spins']=quartic

# The positive Davenport and Heilbronn theta kernel, in the note's u scale.
kap=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
aa=[0,1,kap,-kap,-1]
def dhpotential(u):
    theta=mp.fsum(aa[n%5]*n*mp.exp(-mp.pi*n*n*mp.exp(4*u)/5) for n in range(1,31))
    return -3*u-mp.log(theta)
res['DH_potential']={'U2_at0':str(mp.diff(dhpotential,0,2)),
                     'U4_at0':str(mp.diff(dhpotential,0,4)),
                     'U3_samples':{str(x):str(mp.diff(dhpotential,mp.mpf(x),3)) for x in ['.01','.1','.2','.5','1','2']}}
Path(__file__).with_name('check_ghs_necessity_results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
