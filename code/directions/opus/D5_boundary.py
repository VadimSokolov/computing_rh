# D5 experiment 11, delta=1. (i) exact boundary criterion: Z(h0,h0)=0 iff a>deg Q2 or b>deg Q1, and Y<1 otherwise. (ii) max property Y(Q1 Q2)^2<=max(Y(Q1)^2,Y(Q2)^2) for identical spins.
# (iii) three spin chain: Y^2<=u(J) for identical Q? (iv) large multiplicity. mpmath and numpy, not certified.
import numpy as np, mpmath as mp, time, random
mp.mp.dps=60
random.seed(5); rng=np.random.default_rng(5)
def pm(a,b):
    c=[mp.mpc(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def fac(t,c): return [mp.mpf(1),mp.mpc(0,2*mp.mpf(t)),-(mp.mpf(t)**2)-mp.mpf(c)]  # high to low
def build(L):
    Q=[mp.mpc(1)]
    for t,c in L: Q=pm(Q,fac(t,c))
    return Q
def dp(p):
    d=len(p)-1
    return [p[i]*(d-i) for i in range(d)] if d>0 else [mp.mpc(0)]
def Zp(Q1,Q2,J):
    S=[mp.mpc(0)]*(len(Q1)+len(Q2)-1); a,b=Q1,Q2; term=mp.mpf(1); k=0
    while True:
        D=pm(a,b); off=len(S)-len(D)
        for i,x in enumerate(D): S[off+i]+=term*x
        k+=1
        if len(a)==1 or len(b)==1: break
        a,b=dp(a),dp(b); term=term*J/k
    return S
def Ymax(S):
    r=mp.polyroots(S,maxsteps=4000,extraprec=1500)
    return max(mp.re(z) for z in r)
# (i)
bad=0;tot=0;worstY=-9;t0=time.time()
for it in range(400):
    T1,T2=random.randint(1,4),random.randint(1,3); t0h=random.uniform(-2,2); J=mp.mpf(random.choice([.05,.2,.5,1,2,4]))
    a=random.randint(0,T1); b=random.randint(0,T2)
    L1=[(-t0h,1)]*a+[(random.uniform(-3,3),random.uniform(-2,.99)) for _ in range(T1-a)]
    L2=[(-t0h,1)]*b+[(random.uniform(-3,3),random.uniform(-2,.99)) for _ in range(T2-b)]
    Q1,Q2=build(L1),build(L2); S=Zp(Q1,Q2,J); val=mp.polyval(S,mp.mpc(1,t0h))
    pred=(a>2*T2) or (b>2*T1); isz=abs(val)<mp.mpf(10)**(-40)
    tot+=1; bad+=(pred!=isz)
    if not pred and tot<=150:
        worstY=max(worstY,Ymax(S))
print(f"(i) {tot} cases, disagreements with the criterion {bad}; largest Re zero among the first 150 nonzero cases {mp.nstr(worstY,8)}; {time.time()-t0:.0f}s",flush=True)
# (ii)
mp.mp.dps=40
def rnd(n): return [(random.gauss(0,random.choice([.1,.5,2])),min(1,random.uniform(-1,1.3))) for _ in range(n)]
def nYsq(L,J):
    Q=build(L); y=Ymax(Zp(Q,Q,J)); return mp.mpf(0) if y<1e-20 else y*y
worst=-9;cnt=0;t0=time.time();wc=None
while time.time()-t0<150:
    J=mp.mpf(random.choice([.05,.2,.5,1,1.5,1.9])); A,B=rnd(random.randint(1,2)),rnd(random.randint(1,2))
    e=nYsq(A+B,J)-max(nYsq(A,J),nYsq(B,J)); cnt+=1
    if e>worst: worst,wc=e,(A,B,J)
print(f"(ii) {cnt} random pairs, largest Y(AB)^2-max(Y(A)^2,Y(B)^2)={mp.nstr(worst,6)} at J={wc[2]} A={[(round(t,3),round(c,3)) for t,c in wc[0]]} B={[(round(t,3),round(c,3)) for t,c in wc[1]]}",flush=True)
# (iii)
def uJ(J): 
    return mp.findroot(lambda u:(u-1)**3+8*J*u*(u-1)+J*J*(12*u-4),(mp.mpf(1)/3,mp.mpf(1)),solver='anderson')
worst=-9;cnt=0;t0=time.time()
while time.time()-t0<60:
    J=mp.mpf(random.choice([.1,.5,1,2,4,16])); L=[(0,1)]+rnd(random.randint(0,2)); Q=build(L); S=Zp(pm(Q,Q),Q,J); y=Ymax(S); e=y*y-uJ(J); cnt+=1
    if e>worst: worst=e;wl=(J,[(round(t,3),round(c,3)) for t,c in L])
print(f"(iii) three spin chain, {cnt} tests, largest Y^2-u(J)={mp.nstr(worst,6)} at {wl}",flush=True)
# (iv)
mp.mp.dps=120
for m in (4,8,16,24):
    row=[]
    for J in (mp.mpf('.1'),mp.mpf('.25'),mp.mpf('.5')):
        Q=[mp.mpc(1)]
        for _ in range(m): Q=pm(Q,fac(0,1))
        S=Zp(Q,Q,J)[::2]  # even polynomial in w: coefficients in u=w^2
        r=mp.polyroots(S,maxsteps=6000,extraprec=4000)
        y2=max((abs(z)+mp.re(z))/2 for z in r); row.append(float(y2))
    print(f"(iv) m={m} y^2 at J=.1,.25,.5: {row[0]:.4f} {row[1]:.4f} {row[2]:.4f}  (1-sqrt J)^2: .4675 .2500 .0858",flush=True)
