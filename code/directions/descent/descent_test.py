# Route A descent test: H_a(t)=Re xi'/xi(a+it) exact vs the descended prime side
# Descent from level 2 by the Poisson multiplier e^{(2-a)w} at frequency w=log n is term by term
# Lambda(n) n^-2 -> Lambda(n) n^-a, so the descended series IS the level a prime series cut at log n <= wmax.
import numpy as np, mpmath as mp, json, sys
mp.mp.dps=20
def H_exact(a,t):
    s=mp.mpc(a,t); return float(mp.re(1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)))
def arch(a,t):
    s=mp.mpc(a,t); return float(mp.re(1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2))
M=int(np.exp(12))+1
sieve=np.ones(M,bool); sieve[:2]=False
for i in range(2,int(M**.5)+1):
    if sieve[i]: sieve[i*i::i]=False
ps=np.nonzero(sieve)[0]; ns=[];L=[]
for p in ps:
    q=p
    while q<M: ns.append(q);L.append(np.log(p)); q*=p
ns=np.array(ns,float);L=np.array(L)
ts=np.linspace(10,100,361)
out={}
for a in [1.0,0.875,0.75]:
    He=np.array([H_exact(a,t) for t in ts]); Ar=np.array([arch(a,t) for t in ts])
    res={'min_exact':float(He.min())}
    for wmax in [6,9,12]:
        for taper in ['sharp','fejer']:
            m=np.log(ns)<=wmax; ln=np.log(ns[m]); c=L[m]*ns[m]**(-a)
            w=np.ones_like(ln) if taper=='sharp' else (1-ln/wmax)
            S=np.cos(np.outer(ts,ln))@(c*w)
            Hd=Ar-S; err=Hd-He
            res[f'{taper}{wmax}']={'max_abs_err':float(abs(err).max()),'min_desc':float(Hd.min()),'gain_e^(2-a)w':float(np.exp((2-a)*wmax))}
    out[str(a)]=res
json.dump(out,open(sys.argv[1],'w'),indent=1); print(json.dumps(out,indent=1))
