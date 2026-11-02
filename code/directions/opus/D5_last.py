# D5 experiment 12. (a) three spin chain with Q=(w^2-1)^m: pair (2m,m), y^2 against the square of the largest Legendre zero of P_{2m}.
# (b) hill climb against the max property Y(AB)^2<=max(Y(A)^2,Y(B)^2), floating point, excess rechecked at 60 digits. (c) two copies Curie Weiss check. Not certified.
import numpy as np, mpmath as mp, time
mp.mp.dps=150
def pm(a,b):
    c=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def pw(m):
    p=[mp.mpf(1)]
    for _ in range(m): p=pm(p,[mp.mpf(1),mp.mpf(0),mp.mpf(-1)])
    return p
def dp(p):
    d=len(p)-1
    return [p[i]*(d-i) for i in range(d)] if d>0 else [mp.mpf(0)]
def Zp(Q1,Q2,J):
    S=[mp.mpf(0)]*(len(Q1)+len(Q2)-1); a,b=Q1,Q2; term=mp.mpf(1); k=0
    while True:
        D=pm(a,b); off=len(S)-len(D)
        for i,x in enumerate(D): S[off+i]+=term*x
        k+=1
        if len(a)==1 or len(b)==1: break
        a,b=dp(a),dp(b); term=term*J/k
    return S
print("(a) m, largest Legendre zero squared of P_{2m}, y^2 at J=.01,.1,1,10,100")
for m in (1,2,3,4,6,8):
    p=mp.legendre; x=max(mp.re(z) for z in mp.polyroots(mp.taylor(lambda t:mp.legendre(2*m,t),0,2*m)[::-1],maxsteps=4000,extraprec=2000))
    row=[]
    for J in ('.01','.1','1','10','100'):
        S=Zp(pw(2*m),pw(m),mp.mpf(J))[::2]
        r=mp.polyroots(S,maxsteps=8000,extraprec=6000)
        row.append(float(max((abs(z)+mp.re(z))/2 for z in r)))
    print(f"m={m} x^2={float(x*x):.5f}  "+" ".join(f"{v:.5f}" for v in row),flush=True)
def Qpoly(L):
    Q=np.array([1.0+0j])
    for t,c in L: Q=np.convolve(Q,np.array([1,2j*t,-t*t-c],dtype=complex))
    return Q
def Zpoly(Q,J):
    d=Q.copy(); s=np.convolve(d,d); term=1.0; k=0
    while len(d)>1:
        k+=1; d=np.polyder(d); term=term*J/k; s=np.polyadd(s,term*np.convolve(d,d))
    return s
def Y2(L,J):
    r=np.roots(Zpoly(Qpoly(L),J)); y=np.max(r.real); return float(y*y) if y>1e-7 else 0.0
rng=np.random.default_rng(77); t0=time.time(); res=[]
print("(b) hill climb, J, n_A, n_B, best excess Y(AB)^2-max(Y(A)^2,Y(B)^2)")
for J in (0.1,0.5,1.0,1.5):
    for (na,nb) in ((1,1),(1,2),(2,2)):
        best=-9;ba=None
        for rs in range(25):
            ta,ca=rng.normal(0,1,na),np.minimum(1,rng.uniform(-1,1.2,na)); tb,cb=rng.normal(0,1,nb),np.minimum(1,rng.uniform(-1,1.2,nb))
            def obj(ta,ca,tb,cb):
                A=list(zip(ta,ca));B=list(zip(tb,cb)); return Y2(A+B,J)-max(Y2(A,J),Y2(B,J))
            cur=obj(ta,ca,tb,cb);sig=.5
            for it in range(150):
                a,b,c,d=ta+sig*rng.normal(0,1,na),np.minimum(1,ca+sig*rng.normal(0,1,na)),tb+sig*rng.normal(0,1,nb),np.minimum(1,cb+sig*rng.normal(0,1,nb))
                v=obj(a,b,c,d)
                if v>cur: ta,ca,tb,cb,cur=a,b,c,d,v; sig=min(sig*1.2,1)
                else: sig=max(sig*.98,1e-3)
            if cur>best: best=cur;ba=(ta.copy(),ca.copy(),tb.copy(),cb.copy())
        print(f"J={J} ({na},{nb}) best {best:.3e} ta={np.round(ba[0],3)} ca={np.round(ba[1],3)} tb={np.round(ba[2],3)} cb={np.round(ba[3],3)}",flush=True)
print("elapsed",round(time.time()-t0))
try:
    import flint; print("python-flint available",flint.__version__)
except Exception as e: print("python-flint not on the laptop:",type(e).__name__)
