# Own numerics for review N10. Nothing imported from misc/reviews/nick-*.
import numpy as np, math
from scipy.special import digamma, polygamma, exp1
PSI14=digamma(.25); PSI114=polygamma(1,.25); PSI34=digamma(.75); PSI134=polygamma(1,.75)
XI=(math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
def E(t,b0=.5,J=400):
    t=np.asarray(t,float); out=np.zeros_like(t)
    for j in range(J):
        b=2*j+b0; out+=np.exp(-b*t)/b**2
    B=2*J+b0
    tp=np.maximum(t,1e-300)
    out+=.5*(np.exp(-B*t)/B-t*exp1(B*tp))      # integral tail
    out+=.5*np.exp(-B*t)/B**2                  # +g(J)/2
    out+=(2*t/B**2+4/B**3)*np.exp(-B*t)/12     # -g'(J)/12
    return out
def Parch(t):
    t=np.asarray(t,float)
    return 8*(np.cosh(t/2)-1)+(t/2)*(PSI14-math.log(math.pi))+PSI114/4-E(t)
def pp(N):
    s=np.ones(N+1,bool); s[:2]=False
    for p in range(2,int(N**.5)+1):
        if s[p]: s[p*p::p]=False
    ns=[];Ls=[]
    for p in np.nonzero(s)[0]:
        p=int(p);q=p;lp=math.log(p)
        while q<=N: ns.append(q);Ls.append(lp);q*=p
    o=np.argsort(ns); return np.array(ns,float)[o],np.array(Ls)[o]
def comb(t,ln,u):
    "sum u_n (t-l_n)_+ , ln increasing"
    W0=np.concatenate([[0],np.cumsum(u)]);W1=np.concatenate([[0],np.cumsum(u*ln)])
    j=np.searchsorted(ln,t,side='right')
    return t*W0[j]-W1[j]
def c_of(v):
    v=np.asarray(v); return .5*(v[2:]+v[:-2]-2*v[1:-1]) if False else None
def incr(P):
    # P[j]=Psi(j d), j=0..M ; c(k)=.5(Psi(k+1)+Psi(|k-1|)-2Psi(k))
    c=np.empty(len(P)-1); c[0]=P[1]; c[1:]=.5*(P[2:]+P[:-2]-2*P[1:-1]); return c
def lev(c,stop=True):
    M=len(c); a=np.zeros(M); v=c[0]; K=np.zeros(M); V=np.zeros(M); V[0]=v
    for k in range(1,M):
        acc=c[k]-(a[:k-1] @ c[k-1:0:-1] if k>1 else 0.)
        kap=acc/v; K[k]=kap
        if k>1:
            a[:k-1]=a[:k-1]-kap*a[:k-1][::-1]
        a[k-1]=kap; v=v*(1-kap*kap); V[k]=v
        if v<=0:
            if stop: return K[:k+1],V[:k+1],k
    return K,V,None
def dh():
    "Lambda_f(n) for n<=N via a(n)log n = sum_{d|n} Lf(d) a(n/d)"
    def a(n): return {0:0.,1:1.,2:XI,3:-XI,4:-1.}[n%5]
    N=3000; L=np.zeros(N+1)
    for n in range(2,N+1):
        s=a(n)*math.log(n)
        for d in range(2,int(n**.5)+1) if False else []: pass
        tot=0.
        for d in range(2,n):
            if n%d==0: tot+=L[d]*a(n//d)
        L[n]=s-tot
    ns=np.nonzero(np.abs(L)>1e-12)[0].astype(float)
    return ns,L[ns.astype(int)]
def Pdh_arch(t):
    t=np.asarray(t,float)
    return (t/2)*(math.log(5/math.pi)+PSI34)+PSI134/4-E(t,1.5)
