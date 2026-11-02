# D5 experiment 13, independent checks of hand results. (A) Theorem 4b: predicted zero of Z for Q=exp(v h^2/2)(h^2-c). (B) Theorem 1 with diagonal, J not positive definite.
# (C) Proposition 13: zero of (w^2-1)^mu against r-Jk+i sqrt(J z). mpmath and numpy, not certified.
import mpmath as mp, numpy as np, random
mp.mp.dps=60
def Zgauss(h,v,c,J,K=160):
    P=[mp.mpf(-c),mp.mpf(0),mp.mpf(1)]; S=mp.mpc(0); t=mp.mpf(1)
    for k in range(K):
        val=mp.polyval(P[::-1],h); S+=t*val*val
        dP=[(i+1)*P[i+1] for i in range(len(P)-1)]+[mp.mpf(0),mp.mpf(0)]
        vP=[mp.mpf(0)]+[v*x for x in P]
        P=[(dP[i] if i<len(dP) else 0)+(vP[i] if i<len(vP) else 0) for i in range(len(P)+1)]
        t=t*J/(k+1)
    return S*mp.exp(v*h*h)
print("(A) Theorem 4b")
for (v,c,J) in ((mp.mpf('.5'),mp.mpf(1),mp.mpf('.3')),(mp.mpf('.2'),mp.mpf(1),mp.mpf('.5')),(mp.mpf('1'),mp.mpf(2),mp.mpf('.25'))):
    C=1/(1-J*J*v*v); a=v*J*J*C; Jp=J*C; cp=c-a
    u=cp-2*Jp+1j*mp.sqrt(2*Jp*(2*cp-Jp)); x=mp.sqrt(u); h=x*(1-J*v)
    near=Zgauss(h,v,c,J); off=Zgauss(h+mp.mpf('.05'),v,c,J)
    print(f"v={v} c={c} J={J} predicted h={mp.nstr(h,8)} |Z(h)|={mp.nstr(abs(near),3)} |Z(h+.05)|={mp.nstr(abs(off),3)} bound (1-Jv)sqrt(c-a)={mp.nstr((1-J*v)*mp.sqrt(cp),6)}")
print("(B) Theorem 1, two spins with diagonal, J12 free, worst Re zero minus bound")
rng=np.random.default_rng(3); worst=-9
def heat(Q,a):
    # exp((a/2) d^2) on a polynomial with high-to-low coefficients
    out=np.array(Q,dtype=complex); d=np.array(Q,dtype=complex); t=1.0; k=0
    while len(d)>2:
        d=np.polyder(np.polyder(d)); k+=1; t*=a/(2*k); out=np.polyadd(out,t*d)
    return out
for it in range(4000):
    n1,n2=rng.integers(1,4),rng.integers(1,4); Q1=np.array([1+0j]); Q2=np.array([1+0j]); c1=rng.uniform(-1,1,n1); c2=rng.uniform(-1,1,n2); t1=rng.normal(0,1,n1); t2=rng.normal(0,1,n2)
    for t,c in zip(t1,c1): Q1=np.convolve(Q1,[1,2j*t,-t*t-c])
    for t,c in zip(t2,c2): Q2=np.convolve(Q2,[1,2j*t,-t*t-c])
    J11,J22,J12=rng.uniform(0,.8),rng.uniform(0,.8),rng.uniform(0,2)
    a,b=heat(Q1,J11),heat(Q2,J22); s=np.convolve(a,b); t=1.0;k=0
    while len(a)>1 and len(b)>1:
        k+=1;a,b=np.polyder(a),np.polyder(b);t*=J12/k;s=np.polyadd(s,t*np.convolve(a,b))
    y=np.roots(s).real.max(); bound=max(np.sqrt(max(c1.max()-J11,0)),np.sqrt(max(c2.max()-J22,0)))
    worst=max(worst,y-bound)
print(f"4000 tests, worst excess {worst:.3e}")
print("(C) Proposition 13, Q=(w^2-1)^mu, rightmost zero at small J against 1-J mu/2")
for mu in (1,2,3):
    for J in (mp.mpf('1e-4'),mp.mpf('1e-3')):
        Q=[mp.mpf(1)]
        for _ in range(mu):
            R=[mp.mpf(0)]*(len(Q)+2)
            for i,x in enumerate(Q):
                R[i]+=x; R[i+2]-=x
            Q=R
        d=Q;S=[mp.mpf(0)]*(2*len(Q)-1);t=mp.mpf(1);k=0
        while True:
            D=[mp.mpf(0)]*(2*len(d)-1)
            for i,x in enumerate(d):
                for j,u in enumerate(d): D[i+j]+=x*u
            off=len(S)-len(D)
            for i,x in enumerate(D): S[off+i]+=t*x
            if len(d)==1: break
            m=len(d)-1;d=[d[i]*(m-i) for i in range(m)];k+=1;t=t*J/k
        y=max(mp.re(z) for z in mp.polyroots(S,maxsteps=4000,extraprec=1200))
        print(f"mu={mu} J={mp.nstr(J,2)} (1-y)/(J mu/2)={mp.nstr((1-y)/(J*mu/2),6)}")
