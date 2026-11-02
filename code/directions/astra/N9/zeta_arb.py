"""Independent Arb enclosure of the prime formula and Toeplitz signs."""
import sys
from pathlib import Path
sys.path.insert(0, '/private/tmp/d1lib')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import json
import math
import flint
from flint import arb, ctx
from zeta_review import prime_powers


def psi_arb(denom, m, cutoff=None, omit=None):
    d=arb(1)/denom
    aa=arb(1)/4
    pi=arb.pi()
    linear=(aa.digamma()-pi.log())/2
    tri=(pi*pi+8*arb.const_catalan())/4
    pp=prime_powers(math.floor(math.exp(m/denom))+1)
    data=[(arb(n).log(), arb(p).log()/arb(n).sqrt()) for n,p in pp if (cutoff is None or n<=cutoff) and n!=omit]
    out=[arb(0)]
    sw=arb(0)
    sl=arb(0)
    pos=0
    for j in range(1,m+1):
        t=j*d
        while pos<len(data) and data[pos][0]<t:
            x,w=data[pos]
            sw+=w
            sl+=x*w
            pos+=1
        if pos<len(data) and not data[pos][0]>t:
            raise ArithmeticError('ambiguous boundary')
        count=max(10,math.ceil((ctx.prec*math.log(2)+40)/(2*j/denom)))
        q=(-2*t).exp()
        e=(-t/2).exp()
        series=arb(0)
        for k in range(count):
            u=k+aa
            series+=e/(4*u*u)
            e*=q
        u=count+aa
        tail=e/(4*u*u*(1-q))
        series+=arb(0,tail.upper())
        out.append(4*((t/2).exp()+(-t/2).exp()-2)+linear*t+tri-series-t*sw+sl)
    return out


def covariance(psi):
    return [psi[1]]+[(psi[j+1]+psi[j-1]-2*psi[j])/2 for j in range(1,len(psi)-1)]


def levinson(c):
    coeff=[]
    vs=[c[0]]
    ks=[]
    if not vs[0]>0:
        return vs,ks,1
    for k in range(1,len(c)):
        refl=(c[k]-sum((coeff[j]*c[k-j-1] for j in range(k-1)),arb(0)))/vs[-1]
        ks.append(refl)
        coeff=[coeff[j]-refl*coeff[-j-1] for j in range(k-1)]+[refl]
        vs.append(vs[-1]*(1-refl*refl))
        if not vs[-1]>0:
            return vs,ks,k+1
    return vs,ks,None


def record(x):
    return {'ball':str(x),'mid':float(x.mid()),'radius_upper':math.nextafter(float(x.rad()),math.inf)}


def main():
    ctx.prec=2048
    vals=psi_arb(40,482)
    c=covariance(vals)
    vs,ks,fail=levinson(c[:480])
    result={'library':'python-flint '+flint.__version__,'precision_bits':ctx.prec,'max_a':12,'delta':'1/40','m':480,'first_nonpositive':fail,'all_pivots_positive':fail is None,'v_ratio':record(vs[-1]/vs[0]),'psi_12':record(vals[480]),'psi_12_05':record(vals[482]),'eigen_brackets':[],'truncation':[]}
    print('full',fail,'ratio',vs[-1]/vs[0],flush=True)
    rows=json.loads(Path('zeta_results.json').read_text())['grid']
    for row in rows:
        if row['delta']!=.025:
            continue
        m=row['m']
        lo=format(row['lambda_min']*.995,'.9g')
        hi=format(row['lambda_min']*1.005,'.9g')
        lv,lk,lf=levinson([c[0]-arb(lo)]+c[1:m])
        hv,hk,hf=levinson([c[0]-arb(hi)]+c[1:m])
        res={'a':row['a'],'m':m,'lo':lo,'hi':hi,'lower_proved':lf is None,'upper_proved':hf is not None and hv[-1]<0}
        result['eigen_brackets'].append(res)
        print(res,flush=True)
    for cut,denom,m in [(0,200,149),(7,200,422),(31,40,142),(1000,40,278)]:
        vv,kk,ff=levinson(covariance(psi_arb(denom,m,cutoff=cut)))
        res={'cutoff':cut,'delta':1/denom,'m':m,'first_nonpositive':ff,'negative_proved':vv[-1]<0,'last_positive':record(vv[-2]),'negative':record(vv[-1])}
        result['truncation'].append(res)
        print({k:v for k,v in res.items() if k not in ['last_positive','negative']},flush=True)
    Path('zeta_arb_results.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
