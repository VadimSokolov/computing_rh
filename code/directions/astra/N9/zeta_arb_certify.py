"""Well-conditioned Arb residual certificates for zeta grid sections."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['VECLIB_MAXIMUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys
from pathlib import Path
sys.path.insert(0,'/private/tmp/d1lib')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import json
from flint import arb,ctx
from zeta_arb import psi_arb,covariance
from dh_certify import certify_matrix


def main():
    ctx.prec=256
    result={'arithmetic':'Arb 256 bits and rigorous floating congruence residual bounds','grid':[],'omissions':[],'table':[]}
    vals=psi_arb(40,482)
    c=covariance(vals)
    for t in [.3,.7,1,2,3,5,8,12,12.05]:
        result['table'].append({'t':t,'psi_ball':str(vals[round(40*t)])})
    for a in [1,2,2.2,3,4,6,8,10,12]:
        item=certify_matrix(c[:round(40*a)])
        item.update({'a':a,'delta':.025})
        result['grid'].append(item)
        print('grid',a,item['positive_definite'],item['lower_bound_display'],item['upper_bound_display'],flush=True)
    for denom in [80,200]:
        cc=covariance(psi_arb(denom,12*denom))
        item=certify_matrix(cc)
        item.update({'a':12,'delta':1/denom})
        result['grid'].append(item)
        print('finer',denom,item['positive_definite'],item['lower_bound_display'],item['upper_bound_display'],flush=True)
    rows=json.loads(Path('zeta_results.json').read_text())['single_omission']
    for row in rows:
        q=row['n']
        m=round(row['first_failed_a']*40)
        cc=covariance(psi_arb(40,m,omit=q))
        before=certify_matrix(cc[:-1])
        after=certify_matrix(cc)
        result['omissions'].append({'n':q,'m':m,'a':m/40,'delta':.025,'before':before,'after':after})
        print('omission',q,before['positive_definite'],after['negative_rayleigh'],m/40,flush=True)
    Path('zeta_arb_certificates.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
