"""Ball certificates for the independent finite grid calculation."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import json
import math
import mpmath as mp
import ball_review
from ball_review import B, blog, bexp, bsqrt, set_precision, ball_levinson
from zeta_review import prime_powers


def psi_balls(delta, m, cutoff=None, omit=None):
    pi=B.from_iv(mp.iv.pi)
    psi0=-B.from_iv(mp.iv.euler)-pi/2-3*blog(2)
    trigamma=pi*pi+8*B.from_iv(mp.iv.catalan)
    logs=[]
    weights=[]
    for n,p in prime_powers(math.floor(math.exp(float(delta)*m))+1):
        if (cutoff is None or n<=cutoff) and n!=omit:
            logs.append(blog(n))
            weights.append(blog(p)/bsqrt(n))
    out=[B(0)]
    sw=B(0)
    slw=B(0)
    pos=0
    d=B(delta)
    for j in range(1,m+1):
        t=d*j
        while pos<len(logs) and (t-logs[pos]).positive():
            sw+=weights[pos]
            slw+=weights[pos]*logs[pos]
            pos+=1
        if pos<len(logs) and not (logs[pos]-t).positive():
            raise ArithmeticError('ambiguous prime boundary')
        count=max(12,math.ceil((ball_review.BITS*math.log(2)+40)/(2*j*float(delta))))
        q=bexp(-2*t)
        expterm=bexp(-t/2)
        series=B(0)
        for k in range(count):
            u=B(k)+B('0.25')
            series+=expterm/(4*u*u)
            expterm*=q
        u=B(count)+B('0.25')
        tail=expterm/(4*u*u*(1-q))
        series+=B.raw(tail.m//2,(tail.m+1)//2+tail.r)
        value=4*(bexp(t/2)+bexp(-t/2)-2)+t/2*(psi0-blog(pi))+trigamma/4-series-t*sw+slw
        out.append(value)
    return out


def cov(psi):
    return [psi[1]]+[(psi[j+1]+psi[j-1]-2*psi[j])/2 for j in range(1,len(psi)-1)]


def main():
    set_precision(2048)
    vals=psi_balls('0.025',482)
    c=cov(vals)
    vs,ks,fail=ball_levinson(c[:480])
    print('full grid',fail,'v ratio',vs[-1]/vs[0],'max kappa radius',max(k.radius() for k in ks),flush=True)
    result={'precision_bits':2048,'full_grid_first_nonpositive':fail,'v_ratio':(vs[-1]/vs[0]).bounds(),'max_psi_radius':max(v.radius() for v in vals),'psi_12':vals[480].bounds(),'psi_12_05':vals[482].bounds(),'eigen_brackets':[],'truncation':[]}
    floats=json.loads(Path('zeta_results.json').read_text())
    for row in floats['grid']:
        if row['delta']!=.025:
            continue
        m=row['m']
        lo=B(format(row['lambda_min']*.995,'.9g'))
        hi=B(format(row['lambda_min']*1.005,'.9g'))
        lowc=[c[0]-lo]+c[1:m]
        highc=[c[0]-hi]+c[1:m]
        lv,lk,lf=ball_levinson(lowc)
        hv,hk,hf=ball_levinson(highc)
        entry={'a':row['a'],'m':m,'lo':lo.mid(),'hi':hi.mid(),'lower_proved':lf is None,'upper_proved':hf is not None and hv[-1].negative()}
        result['eigen_brackets'].append(entry)
        print(entry,flush=True)
    for cut,d,m in [(0,'0.005',149),(7,'0.005',422),(31,'0.025',142),(1000,'0.025',278)]:
        cv=cov(psi_balls(d,m,cutoff=cut))
        vv,kk,ff=ball_levinson(cv)
        entry={'cutoff':cut,'delta':d,'m':m,'first_nonpositive':ff,'negative_proved':vv[-1].negative(),'last_positive_pivot':vv[-2].bounds(),'negative_pivot':vv[-1].bounds()}
        result['truncation'].append(entry)
        print(entry,flush=True)
    Path('zeta_ball_results.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
