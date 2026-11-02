"""Independent evaluation of the screw function and finite increment tests.

No inputs or code are taken from the manuscript directory.
"""
import json
import math
from pathlib import Path
import numpy as np
from scipy.special import digamma, polygamma
from scipy.linalg import toeplitz, eigh


def prime_powers(limit):
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
    terms = []
    for p in np.flatnonzero(sieve):
        n = int(p)
        while n <= limit:
            terms.append((n, int(p)))
            n *= int(p)
    return sorted(terms)


def arch(t):
    if t == 0:
        return 0.0
    k = np.arange(max(20, math.ceil(38 / (2*t))), dtype=float) + .25
    series = np.sum(np.exp(-2*k*t)/(4*k*k))
    return 8*(np.cosh(t/2)-1)+t/2*(digamma(.25)-math.log(math.pi))+polygamma(1,.25)/4-series


def psi_grid(delta, length, cutoff=None, omit=None):
    nmax = math.floor(math.exp(length*delta))
    pp = prime_powers(nmax)
    if cutoff is not None:
        pp = [(n,p) for n,p in pp if n<=cutoff]
    if omit is not None:
        pp = [(n,p) for n,p in pp if n!=omit]
    logs = np.array([math.log(n) for n,p in pp])
    weights = np.array([math.log(p)/math.sqrt(n) for n,p in pp])
    sw = np.r_[0.,np.cumsum(weights)]
    slw = np.r_[0.,np.cumsum(weights*logs)]
    out=[]
    for j in range(length+1):
        t=j*delta
        idx=np.searchsorted(logs,t,side='right')
        out.append(arch(t)-t*sw[idx]+slw[idx])
    return np.array(out)


def covariance(psi):
    return np.r_[psi[1], (psi[2:]+psi[:-2]-2*psi[1:-1])/2]


def levinson(c):
    coeff=np.empty(0)
    v=c[0]
    vs=[v]
    ks=[]
    fail=None
    for k in range(1,len(c)):
        refl=(c[k]-np.dot(coeff,c[k-1:0:-1]))/v
        ks.append(float(refl))
        coeff=np.r_[coeff-refl*coeff[::-1],refl]
        v=v*(1-refl*refl)
        vs.append(float(v))
        if v<=0:
            fail=k+1
            break
    return np.array(vs),np.array(ks),fail


def main():
    out={'method':'float64, independent exponential gamma series; not certified','table':[], 'grid':[], 'truncation':[], 'single_omission':[]}
    pp=prime_powers(math.floor(math.exp(12.1)))
    for t in [.3,.7,1,2,3,5,8,12,12.05]:
        value=arch(t)-math.fsum(math.log(p)/math.sqrt(n)*(t-math.log(n)) for n,p in pp if math.log(n)<=t)
        out['table'].append({'t':t,'psi':value,'arch':arch(t)})
    for d in [.025,.0125,.005]:
        c=covariance(psi_grid(d,round(12/d)))
        vs,ks,fail=levinson(c)
        aa=[1,2,3,4,6,8,10,12] if d==.025 else [12]
        for a in aa:
            m=round(a/d)
            matrix=toeplitz(c[:m])
            eig=eigh(matrix,subset_by_index=[0,0],eigvals_only=True)[0]
            out['grid'].append({'a':a,'delta':d,'m':m,'lambda_min':float(eig),'weil_min':float(2*eig/d),'v_ratio':float(vs[m-1]/vs[0]),'max_kappa':float(max(abs(ks[:m-1]))),'T_Nyquist':math.pi/d,'X':math.exp(a),'T_star':2*math.pi*math.exp(a),'heuristic_digits':4*math.pi*math.exp(a)/math.log(10),'needed_m':math.ceil(2*a*math.exp(a)),'first_failure':fail})
    for cut,maxa in [(0,1),(7,2.5),(31,4),(1000,7.5)]:
        for d in [.025,.005]:
            c=covariance(psi_grid(d,round(maxa/d),cutoff=cut))
            vs,ks,fail=levinson(c)
            out['truncation'].append({'cutoff':cut,'delta':d,'first_failed_m':fail,'first_failed_a':None if fail is None else fail*d,'last_positive_pivot':float(vs[-2]),'negative_pivot':float(vs[-1])})
    for q in [2,3,4,5,7,8,9,11,13,32,1009]:
        d=.025
        c=covariance(psi_grid(d,480,omit=q))
        vs,ks,fail=levinson(c)
        out['single_omission'].append({'n':q,'first_failed_a':None if fail is None else fail*d,'pivot':float(vs[-1])})
    Path('zeta_results.json').write_text(json.dumps(out,indent=2))
    for key,rows in out.items():
        if key != 'method':
            print(key,rows,flush=True)


if __name__=='__main__':
    main()
