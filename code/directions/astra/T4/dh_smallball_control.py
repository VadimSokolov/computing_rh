"""Davenport-Heilbronn control for the small-ball rate, using Arb balls."""
from pathlib import Path
import sys
import json
sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import mpmath as mp
from flint import arb, arb_series, ctx

mp.mp.dps = 75
ctx.dps = 85
BASE = Path(__file__).resolve().parent
kb = ((10 - 2*arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
km = (mp.sqrt(10 - 2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)


def dm(z):
    return 5**(-z)*sum(c*mp.zeta(z, mp.mpf(j)/5)
                       for j,c in ((1,1),(2,km),(3,-km),(4,-1)))


def lm(z):
    return (z+1)/2*mp.log(5/mp.pi)+mp.loggamma((z+1)/2)+mp.log(dm(z))


def lb(z):
    v=arb_series([z,1],prec=3)
    d=(-v*arb(5).log()).exp()*sum(
        c*v.zeta(arb(j)/5) for j,c in ((1,arb(1)),(2,kb),(3,-kb),(4,arb(-1))))
    return (v+1)/2*(5/arb.pi()).log()+((v+1)/2).lgamma()+d.log()


rows=[]
for es in ("0.004","0.003","0.001","0.0001","0.000001"):
    e=mp.mpf(es)
    x0=-mp.lambertw(-8*mp.pi*e/5,-1).real/(4*e)
    x=mp.findroot(lambda x: mp.diff(lm,2+x)-2*e*x,(x0,x0+1))
    r=x*mp.mpf("1e-45")
    lo,hi=[arb(mp.nstr(t,70)) for t in (x-r,x+r)]
    E=arb(es)
    assert lb(2+lo)[1]-2*E*lo>0
    assert lb(2+hi)[1]-2*E*hi<0
    X=arb(mp.nstr(x,70),mp.nstr(2*r,70))
    I=lb(2+X)[0]-lb(arb(2))[0]-E*X*X
    leading=(1/E).log()**2/(16*E)
    rows.append({"epsilon":es,"alpha":2,"I_DH_ball":str(I),
                 "I_DH_over_universal_leading_ball":str(I/leading)})
data={"status":"computed","arithmetic":"Arb balls", "rows":rows,
      "kappa_ball":str(kb),
      "zero_free_halfplane_Re_ge2_bound":str(arb(2).zeta()-1)}
(BASE/"dh_smallball_results.json").write_text(json.dumps(data,indent=2)+"\n")
print(json.dumps(data,indent=2))
