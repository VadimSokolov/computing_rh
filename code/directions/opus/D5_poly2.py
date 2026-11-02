# D5 experiment 6: conjecture check for identical spins, local model, with every apparent excess rechecked at 100 digits. Not certified.
import numpy as np, math, mpmath as mp, time
P=np.polynomial.polynomial; mp.mp.dps=100
f=lambda J:(1-2*J+math.sqrt(1+2*J*J))/2 if J<2 else 0.0
def Zpoly(q,J):
    out=np.zeros(2*len(q)-1,dtype=complex); d=q.copy(); fac=1.0; k=0
    while True:
        m=P.polymul(d,d); out[:len(m)]+=fac*m
        if len(d)==1: break
        d=P.polyder(d); k+=1; fac*=J/k
    return out
def Qfrom(ts,ds):
    q=np.array([1.0+0j])
    for t,d in zip(ts,ds): q=P.polymul(q,np.array([(1j*t)**2-d**2,2j*t,1.0]))
    return q
def ymax_mp(ts,ds,J):
    # exact coefficient arithmetic in mpmath
    q=[mp.mpc(1)]
    def mul(a,b):
        c=[mp.mpc(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b): c[i+j]+=x*y
        return c
    for t,d in zip(ts,ds): q=mul(q,[(mp.mpc(0,float(t)))**2-mp.mpf(float(d))**2,mp.mpc(0,2*float(t)),mp.mpc(1)])
    out=[mp.mpc(0)]*(2*len(q)-1); dq=q; fac=mp.mpf(1); k=0; J=mp.mpf(float(J))
    while True:
        m=mul(dq,dq)
        for i,x in enumerate(m): out[i]+=fac*x
        if len(dq)==1: break
        dq=[(i+1)*dq[i+1] for i in range(len(dq)-1)]; k+=1; fac*=J/k
    r=mp.polyroots(out[::-1],maxsteps=4000,extraprec=1200)
    return max(mp.re(x) for x in r)
rng=np.random.default_rng(11); t0=time.time(); tot=0; rechecked=0; realexcess=0; worst=-9; closest=[]
Js=[0.02,0.05,0.1,0.2,0.3,0.5,0.75,1.0,1.25,1.5,1.75,1.9,2.0,2.2,3.0]
while time.time()-t0<150:
    J=float(rng.choice(Js)); n=int(rng.integers(2,6))
    s=float(rng.choice([0.03,0.1,0.3,1,3]))*math.sqrt(max(J,0.05))
    ts=rng.normal(0,s,n); ts[0]=0
    ds=np.ones(n) if rng.random()<0.6 else np.concatenate([[1.0],rng.uniform(0,1,n-1)])
    r=P.polyroots(Zpoly(Qfrom(ts,ds),J)); y=max(r.real.max(),0); tot+=1
    ex=y*y-f(J)
    if ex>-2e-5:
        rechecked+=1; y2=float(ymax_mp(ts,ds,J)); y2=max(y2,0); ex2=y2*y2-f(J); worst=max(worst,ex2)
        if ex2>1e-30: realexcess+=1; print("EXCESS",J,ts,ds,ex2)
    else: worst=max(worst,ex)
print(f"tests={tot}, rechecked at 100 digits={rechecked}, true excesses={realexcess}, largest y^2-f(J) seen={worst:.3e}")
print("multiplicity: Q=(w^2-1)^m, identical: y^2 at J; closure J* (first J on a grid with y=0)")
for m in (1,2,3,4,6):
    row=[]; 
    for J in (0.02,0.1,0.25,0.5,1.0):
        y=ymax_mp([0.0]*m,[1.0]*m,J); row.append(max(float(y),0)**2)
    lo,hi=0.0,2.0
    for _ in range(14):
        mid=(lo+hi)/2
        if float(ymax_mp([0.0]*m,[1.0]*m,mid))>1e-20: lo=mid
        else: hi=mid
    print(f" m={m}: y^2(J=.02,.1,.25,.5,1)="+" ".join(f"{x:.4f}" for x in row)+f"  J*~{hi:.4f}")
