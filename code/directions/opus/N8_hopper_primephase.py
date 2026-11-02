#!/usr/bin/env python3
"""Proposal H1 (Hopper, not run there). For a table of zeta ordinates, the first height at which a reference phase with the primes up to P stops admitting rank one.
Input: text file, one ordinate per line, or a JSON list. For increasing phase, levels alternate with the zeros for some constant c iff e_n = n - Theta_P(gamma_n)/pi has width < 1, where Theta_P = theta + Im log zeta_P(1/2+it).
Usage: timeout 3000 python3 -I N8_hopper_primephase.py ZEROFILE [n0]. Memory O(N); time O(N * number of prime powers used)."""
import sys, json, math
import numpy as np
Bn=[1/12,-1/360,1/1260,-1/1680,1/1188]
def theta(t,shift=12):
    z=0.25+0.5j*np.asarray(t,dtype=float); s=np.zeros_like(z); w=z.copy()
    for k in range(shift):
        s=s-np.log(w); w=w+1
    iw=1/w; iw2=iw*iw; ser=0; p=iw
    for b in Bn:
        ser=ser+b*p; p=p*iw2
    return (s+(w-0.5)*np.log(w)-w+ser).imag-0.5*np.asarray(t)*math.log(math.pi)
path=sys.argv[1]; n0=int(sys.argv[2]) if len(sys.argv)>2 else 8
Z=np.array(json.load(open(path)),dtype=float) if path.endswith('.json') else np.loadtxt(path)
idx=np.arange(1,len(Z)+1); e=idx-theta(Z)/math.pi
PR=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]
sel=idx>=n0
def report(tag,e):
    x=e[sel]; hi=np.maximum.accumulate(x); lo=np.minimum.accumulate(x); w=hi-lo; bad=np.where(w>=1)[0]
    first=('n=%d t=%.3f'%(idx[sel][bad[0]],Z[sel][bad[0]])) if len(bad) else 'none'
    print('%-8s width %.4f rms %.4f first height with width>=1: %s'%(tag,w[-1],x.std(),first),flush=True)
report('P=1',e)
for p in PR:
    lp=math.log(p); kmax=max(1,int(2*math.log(1e4)/lp))
    for k in range(1,kmax+1):
        e=e+np.sin(k*Z*lp)/(k*p**(k/2))/math.pi
    report('P=%d'%p,e)
