# N9_common.py  Own numerics for review N9 (nothing imported from misc/reviews/nick-*).
# Psi(t) = sum_rho (1-cos(gamma t))/gamma^2, closed form of N9.md eq (1.2):
#   Psi(t) = 8(cosh(t/2)-1) + (t/2)(psi(1/4)-log pi) + psi'(1/4)/4 - E(t) - sum_{n<=e^t} Lambda(n) n^{-1/2} (t-log n),
#   E(t) = sum_{j>=0} exp(-b_j t)/b_j^2, b_j = 2j+1/2.
# T_m(delta) = Toeplitz of c(m) = (Psi((m+1)d)+Psi((m-1)d)-2Psi(md))/2.
import sys, math, json, time
sys.path.insert(0,'/tmp/d1lib')
import numpy as np
from flint import arb, acb, ctx
from scipy.special import digamma, polygamma
from scipy.linalg import eigvalsh, eigh, toeplitz, solve_toeplitz

ROOT='.'
OPUS=ROOT+'/misc/research/directions/opus'
KAPPA=(math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)

def load_zeros(K=8000):
    with open(ROOT+'/data/zeros_8000.txt') as f:
        strs=[l.strip() for l in f if l.strip()]
    return np.array([float(s) for s in strs[:K]])

def prime_powers(N):
    s=np.ones(N+1,bool); s[:2]=False
    for p in range(2,int(N**.5)+1):
        if s[p]: s[p*p::p]=False
    ns=[];Ls=[]
    for p in np.nonzero(s)[0]:
        p=int(p); q=p; lp=math.log(p)
        while q<=N:
            ns.append(q); Ls.append(lp); q*=p
    idx=np.argsort(ns)
    return np.array(ns,dtype=np.int64)[idx], np.array(Ls)[idx]

def series_E(t,b0):
    J=int(23.0/t)+5
    b=2.0*np.arange(J)+b0
    return float(np.sum(np.exp(-b*t)/b**2))

def arch_cov(delta,M,b0):
    c=np.zeros(M)
    s0=0.25*float(polygamma(1,b0/2))
    E1=series_E(delta,b0); E2=series_E(2*delta,b0)
    if M>1: c[1]=-0.5*(E2+s0-2*E1)
    for m in range(2,M):
        J=int(23.0/((m-1)*delta))+5
        b=2.0*np.arange(J)+b0
        c[m]=-0.5*np.sum(np.exp(-b*(m-1)*delta)*(1-np.exp(-b*delta))**2/b**2)
    return c,s0,E1

def add_primes(c,delta,ln,u):
    M=len(c)
    x=ln/delta; m0=np.floor(x).astype(np.int64); eta=x-m0
    ok=m0<M
    np.add.at(c,m0[ok],-0.5*u[ok]*(1-eta[ok])*delta)
    ok2=(m0+1)<M
    np.add.at(c,m0[ok2]+1,-0.5*u[ok2]*eta[ok2]*delta)

def cov_arch_zeta(delta,M):
    c,s0,E1=arch_cov(delta,M,0.5)
    mm=np.arange(M)
    c[1:]+=16*np.cosh(mm[1:]*delta/2)*math.sinh(delta/4)**2
    c[0]=8*(math.cosh(delta/2)-1)+0.5*delta*(float(digamma(0.25))-math.log(math.pi))+s0-E1
    return c

def cov_zeta(delta,M,ln,u):
    c=cov_arch_zeta(delta,M)
    add_primes(c,delta,ln,u)
    return c

def dh_coeffs(N):
    a5=np.array([0.,1.,KAPPA,-KAPPA,-1.])
    an=a5[np.arange(N+1)%5]
    logn=np.zeros(N+1); logn[1:]=np.log(np.arange(1,N+1))
    S=np.zeros(N+1); Lam=np.zeros(N+1)
    for d in range(1,N+1):
        Lam[d]=an[d]*logn[d]-S[d]
        k=N//d
        if k>=2 and Lam[d]!=0.0:
            S[2*d:d*k+1:d]+=Lam[d]*an[2:k+1]
    return Lam

def cov_dh(delta,M,Lam):
    c,s0,E1=arch_cov(delta,M,1.5)
    c[0]=0.5*delta*(math.log(5/math.pi)+float(digamma(0.75)))+s0-E1
    n=np.arange(2,len(Lam)); mask=Lam[2:]!=0
    add_primes(c,delta,np.log(n[mask]),Lam[2:][mask]/np.sqrt(n[mask]))
    return c

def cov_from_zeros(delta,M,gam):
    w=4*np.sin(gam*delta/2)**2/gam**2
    return np.array([np.sum(w*np.cos(m*gam*delta)) for m in range(M)])

def levinson(c,stop=True):
    n=len(c); a=np.zeros(n); v=c[0]; ks=np.zeros(n); vs=np.zeros(n); vs[0]=v
    if v<=0: return ks[:1],vs[:1],0
    for k in range(1,n):
        acc=c[k]-(np.dot(a[:k-1],c[k-1:0:-1]) if k>1 else 0.0)
        kap=acc/v
        prev=a[:k-1].copy()
        a[:k-1]=prev-kap*prev[::-1]
        a[k-1]=kap
        v=v*(1-kap*kap)
        ks[k]=kap; vs[k]=v
        if v<=0 and stop: return ks[:k+1],vs[:k+1],k
    return ks,vs,None

def lam_min(c,m,want_vec=False,parity=None,iters=40):
    """least eigenvalue of T_m; dense below 1600, else inverse iteration (parity keeps the sector)."""
    r=c[:m].copy()
    if m<=1600 and parity is None:
        T=toeplitz(r)
        w,V=eigh(T,subset_by_index=[0,0])
        return (w[0],V[:,0]) if want_vec else w[0]
    rng=np.random.default_rng(1)
    x=rng.standard_normal(m)
    if parity=='e': x=x+x[::-1]
    if parity=='o': x=x-x[::-1]
    x/=np.linalg.norm(x); lam=None
    for it in range(iters):
        y=solve_toeplitz(r,x)
        if parity=='e': y=0.5*(y+y[::-1])
        if parity=='o': y=0.5*(y-y[::-1])
        lamn=1.0/np.dot(x,y)
        y/=np.linalg.norm(y); x=y
        if lam is not None and abs(lamn-lam)<1e-7*abs(lamn): lam=lamn; break
        lam=lamn
    return (lam,x) if want_vec else lam

def peak_freq(u,delta):
    """frequency (in gamma units) where |sum u_k e^{ik theta}| peaks, theta in (0,pi)."""
    n=1<<int(math.ceil(math.log2(8*len(u)+1)))
    F=np.abs(np.fft.rfft(u,n))
    k=int(np.argmax(F))
    return 2*math.pi*k/n/delta

# ------- ball arithmetic
def consts():
    pi=arb.pi(); g=arb.const_euler(); l2=arb.const_log2(); G=arb.const_catalan()
    return dict(pi=pi, psi14=-g-pi/2-3*l2, psi34=-g+pi/2-3*l2, tri14=pi*pi+8*G, tri34=pi*pi-8*G)

def fl(x):
    return float(x.mid())

def E_ball(t,b0):
    """sum_{j>=0} exp(-(2j+b0)t)/(2j+b0)^2 with the tail in the radius. t arb > 0, b0 arb."""
    x=(-t).exp(); x2=x*x
    term=(-(t*b0)).exp()
    J=int((ctx.prec*0.6932+30)/(2*fl(t)))+2
    s=arb(0); b=b0
    for j in range(J):
        s+=term/(b*b)
        term*=x2; b+=2
    return s+(term/(b*b)/(1-x2))*arb(0.5,0.5)   # tail in [0, first omitted/(1-x2)]

def arch_ball(t,C=None):
    """8(cosh(t/2)-1)+(t/2)(psi(1/4)-log pi)+psi'(1/4)/4-E(t), the archimedean part of N9.md (1.2)."""
    C=C or consts()
    return 8*((t/2).cosh()-1)+(t/2)*(C['psi14']-C['pi'].log())+C['tri14']/4-E_ball(t,arb(1)/2)

def pp_ball(N):
    ns,Ls=prime_powers(N)
    out=[]
    for n,lp in zip(ns,Ls):
        n=int(n)
        # p from n: log p = Lambda(n)
        p=int(round(math.exp(lp)))
        out.append((math.log(n),arb(n).log(),arb(p).log()/arb(n).sqrt(),n))
    return out

def prime_sum_ball(t,pps):
    tf=fl(t); s=arb(0)
    for lf,ln,u,n in pps:
        if lf<tf-1e-9: s+=u*(t-ln)
        elif lf>tf+1e-9: break
        else:
            d=t-ln
            if d>0: s+=u*d
            elif d<0: pass
            else: raise RuntimeError('tie at n=%d'%n)
    return s

def Psi_grid_ball(delta,K,pps,C=None,kind='zeta',weights=None):
    """Psi(k delta) for k=0..K, cumulative prime sum. kind='zeta'. delta arb."""
    C=C or consts(); P=[arb(0)]
    A=arb(0); B=arb(0); ip=0; np_=len(pps)
    for k in range(1,K+1):
        t=delta*k; tf=fl(t)
        while ip<np_ and pps[ip][0]<tf:
            lf,ln,u,n=pps[ip]
            if abs(lf-tf)<1e-9:
                d=t-ln
                if not (d>0 or d<0): raise RuntimeError('tie')
                if d<0: break
            A+=u; B+=u*ln; ip+=1
        P.append(arch_ball(t,C)-(t*A-B))
    return P

def cov_from_Psi_ball(P):
    return [(P[m+1]+P[abs(m-1)]-2*P[m])/2 for m in range(len(P)-1)]

def levinson_ball(c,shift=None,record=True):
    n=len(c)
    v=c[0] if shift is None else c[0]-shift
    a=[]; ks=[arb(0)]; vs=[v]
    if not (v>0): return ks,vs,0
    for k in range(1,n):
        acc=c[k]
        for j in range(k-1): acc=acc-a[j]*c[k-1-j]
        kap=acc/v
        a=[a[j]-kap*a[k-2-j] for j in range(k-1)]+[kap]
        v=v*(1-kap*kap)
        if record: ks.append(kap); vs.append(v)
        else: ks=[kap]; vs=[v]
        if not (v>0): return ks,vs,k
    return ks,vs,None

def E_law(N,L):
    """D2.md: E(N,L), saturating at 4 pi X. Returns E,b,phi."""
    X=math.exp(2*L); T=math.pi*N/L
    if T>=4*math.pi*X: return 4*math.pi*X,4*math.pi*X,math.pi/2
    lo,hi=1e-9,math.pi/2-1e-12
    f=lambda ph:(4*L/math.pi)*ph-math.log(T/(4*math.pi*math.sin(ph)))
    for _ in range(200):
        mid=0.5*(lo+hi)
        if f(mid)>0: hi=mid
        else: lo=mid
    ph=0.5*(lo+hi); b=T/math.sin(ph)
    return b+4*N*math.log(1/math.tan(ph/2)),b,ph

# ------- verified positivity of a Toeplitz matrix of balls: dyadic Cholesky with residual (opus/D1.md Section 3.9)
from flint import arb_mat
def ub(x):
    return (abs(float(x.mid()))+float(x.rad()))*(1+1e-12)
def chol_cert(c,mu=0.0):
    """True if T(c)-mu I is certified positive definite. c list of arb. mu float (exact dyadic)."""
    m=len(c)
    mid=np.array([fl(x) for x in c]); mid[0]-=mu
    L=np.linalg.cholesky(toeplitz(mid))
    X=np.linalg.inv(L)
    Lm=arb_mat(m,m,[arb(float(v)) for v in L.ravel()])
    Xm=arb_mat(m,m,[arb(float(v)) for v in X.ravel()])
    Am=arb_mat(m,m,[(c[abs(i-j)]-arb(mu)) if i==j else c[abs(i-j)] for i in range(m) for j in range(m)])
    E=Am-Lm*Lm.transpose()
    S=Xm*Lm
    fE=0.0; fX=0.0; fS=0.0
    for i in range(m):
        for j in range(m):
            fE+=ub(E[i,j])**2
            fX+=float(X[i,j])**2
            fS+=ub(S[i,j]-(1 if i==j else 0))**2
    fE=math.sqrt(fE)*(1+1e-12); fX=math.sqrt(fX)*(1+1e-12); fS=math.sqrt(fS)*(1+1e-12)
    if fS>=1: return False,(fE,fX,fS)
    Lin=fX/(1-fS)
    return (fE*Lin*Lin<1),(fE,fX,fS,fE*Lin*Lin)
def rayleigh_ball(c,u):
    """ball value of u^T T u / u^T u, u exact dyadic (double)."""
    m=len(c); s=arb(0); uu=arb(0)
    ua=[arb(float(x)) for x in u]
    for i in range(m):
        uu+=ua[i]*ua[i]
        for j in range(m):
            s+=ua[i]*ua[j]*c[abs(i-j)]
    return s/uu
def S_tail(H):
    """upper bound for sum_{|gamma|>H} gamma^-2, both signs, from N(T) (HSW22 error term)."""
    return (math.log(H/(2*math.pi))+1)/(math.pi*H)+60/H**2
