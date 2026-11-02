# D5 experiment 10: adversarial hill climb against the conjecture Y^2 <= delta^2 f(J/delta^2) for identical local polynomials.
# Parameters: heights t_j and c_j <= 1 with c_0=1 (delta=1). Objective Y^2-f(J). Floating point search, best points rechecked at 100 digits. Not certified.
import numpy as np, mpmath as mp, time
def f(J): return 0.5*(1-2*J+np.sqrt(1+2*J*J)) if J<2 else 0.0
def Qpoly(ts,cs):
    Q=np.array([1.0+0j])
    for t,c in zip(ts,cs):
        Q=np.convolve(Q,np.array([1,2j*t,-t*t-c],dtype=complex))
    return Q
def Zpoly(Q,J):
    d=Q.copy(); s=np.convolve(d,d); term=1.0; k=0
    while len(d)>1:
        k+=1; d=np.polyder(d); term=term*J/k
        s=np.polyadd(s,term*np.convolve(d,d))
    return s
def Y2(ts,cs,J):
    r=np.roots(Zpoly(Qpoly(ts,cs),J)); return float(np.max(r.real)**2) if np.max(r.real)>0 else 0.0
def clip(c):
    c=np.minimum(c,1.0); c[0]=1.0; return c
rng=np.random.default_rng(55); t0=time.time(); best_all=[]
for n in (2,3,4):
    for J in (0.1,0.5,1.0,1.5,1.9):
        best=-9; arg=None
        for rs in range(30):
            ts=rng.normal(0,rng.choice([0.05,0.3,1,3]),n); cs=clip(rng.uniform(-1,1,n)); cur=Y2(ts,cs,J)-f(J); sig=0.5
            for it in range(250):
                t2=ts+sig*rng.normal(0,1,n)*rng.choice([.03,.3,1]); c2=clip(cs+sig*rng.normal(0,1,n)*rng.choice([.03,.3,1])); v=Y2(t2,c2,J)-f(J)
                if v>cur: ts,cs,cur=t2,c2,v; sig=min(sig*1.2,1)
                else: sig=max(sig*.985,1e-3)
            if cur>best: best,arg=cur,(ts.copy(),cs.copy())
        mp.mp.dps=100
        ts,cs=arg; Q=[mp.mpc(1)]
        for t,c in zip(ts,cs):
            q=[mp.mpc(1),mp.mpc(0,2*mp.mpf(float(t))),-(mp.mpf(float(t))**2)-mp.mpf(float(c))]
            R=[mp.mpc(0)]*(len(Q)+2)
            for i,x in enumerate(Q):
                for j,y in enumerate(q): R[i+j]+=x*y
            Q=R
        def dp(p):
            d=len(p)-1; return [p[i]*(d-i) for i in range(d)]
        d=Q; term=mp.mpf(1); k=0; S=[mp.mpc(0)]*(2*len(Q)-1)
        while True:
            D=[mp.mpc(0)]*(2*len(d)-1)
            for i,x in enumerate(d):
                for j,y in enumerate(d): D[i+j]+=x*y
            off=len(S)-len(D)
            for i,x in enumerate(D): S[off+i]+=term*x
            if len(d)==1: break
            d=dp(d); k+=1; term=term*mp.mpf(J)/k
        r=mp.polyroots(S,maxsteps=5000,extraprec=1500); Ym=max(mp.re(z) for z in r); ex=Ym**2-(1-2*mp.mpf(J)+mp.sqrt(1+2*mp.mpf(J)**2))/2 if Ym>0 else -1
        print(f"n={n} J={J}: float best Y^2-f={best:.3e}, 100 digit recheck {mp.nstr(ex,5)}, min height gap {min(abs(np.diff(np.sort(ts)))) :.2f}, c={np.round(cs,3)}",flush=True)
    if time.time()-t0>420: print("time cap"); break
print("elapsed",round(time.time()-t0))
