# Local band model: f=e^{-iTx}h(x), h on [0,a], modes nu_k=2 pi k/a, |nu_k|<=B.
# W(f) ~ L|h|^2 - sum_{n<=e^a} u_n (phase terms) + planted pair T +- i beta. No Nick code.
import sys; sys.path.insert(0,'/tmp/n10rev')
from core import *
from scipy.optimize import brentq
ROOT='.'
ns,Ls=pp(300000); ln=np.log(ns); u=Ls/np.sqrt(ns)
lnl=np.log(ns.astype(np.longdouble))
TP=2*np.longdouble(np.pi)
z=np.array([float(x) for x in open(ROOT+'/data/zeros_8000.txt').read().split()])
def base(T,a,B):
    K=int(B*a/(2*math.pi)); k=np.arange(-K,K+1); nu=2*math.pi*k/a; D=len(k)
    dn=nu[None,:]-nu[:,None]   # nu_k-nu_j
    ph=np.array(np.fmod(np.longdouble(T)*lnl[ns<=math.exp(a)],TP),dtype=float)
    M=math.log(T/(2*math.pi))*np.eye(D,dtype=complex)
    for i,(l,up,p) in enumerate(zip(ln[ns<=math.exp(a)],u,ph)):
        s=a-l
        with np.errstate(divide='ignore',invalid='ignore'):
            I=np.where(np.abs(dn)<1e-12,s+0j,(np.exp(1j*dn*s)-1)/(1j*dn))/a
        A=np.exp(-1j*p)*np.exp(1j*nu*l)[None,:]*I
        M-=up*(A+A.conj().T)
    return M,nu
def plant(nu,a,b):
    p=(np.exp((-b+1j*nu)*a)-1)/(-b+1j*nu)/math.sqrt(a); q=(np.exp((b+1j*nu)*a)-1)/(b+1j*nu)/math.sqrt(a)
    return np.outer(q.conj(),p)+np.outer(p.conj(),q)
def mean_root(L,b):
    f=lambda a:math.sinh(b*a)/b-a-L
    return brentq(f,1e-3,600/b)
# validation of the base against the zero sum, T=2000
T=1999.99;a=4.0;M,nu=base(T,a,30.)
Z=z[(z>T-400)&(z<T+400)]; G=np.array([[ (np.exp(1j*(n+g-T)*a)-1)/(1j*(n+g-T))/math.sqrt(a) for n in nu] for g in Z])
Mz=G.conj().T @ G
print('validate T=2000 a=4: base lmin',np.linalg.eigvalsh(M)[0],'zero sum lmin',np.linalg.eigvalsh(Mz)[0],' rel diff (fro)',np.linalg.norm(M-Mz)/np.linalg.norm(Mz),'diag mean',M.diagonal().real.mean(),Mz.diagonal().real.mean())
print('book table zh:tab:det windows (2a): g=199.64: .01 3.1 .05 3.0 .2 3.0 .45 3.0; g=1999.99 <=: 6.1 6.0 5.9 4.9')
def first(T,b,amin,amax,da,B,mode='model'):
    out=None
    for a in np.arange(amin,amax+1e-9,da):
        M,nu=base(T,a,B); lb=np.linalg.eigvalsh(M)[0]; lp=np.linalg.eigvalsh(M+plant(nu,a,b))[0]
        if lp<-0.02: return round(float(a),3),round(float(lb),3),round(float(lp),3)
    return None,round(float(lb),3),round(float(lp),3)
for T,rng in [(199.64,(1.5,4.5,.1)),(1999.99,(3.0,7.0,.1))]:
    for b in (.2,.45):
        print('T',T,'beta',b,'model first neg (a, base lmin, plant lmin):',first(T,b,*rng,30.),' L=',round(math.log(T/2/math.pi),3),' mean root',round(mean_root(math.log(T/2/math.pi),b),3))
import time
for T0 in (2.99e12,2.9937e12,2.87e12):
    for b in (.5,.3):
        t0=time.time(); r=first(T0,b,4.0,12.5,.25,20.); print('T',T0,'beta',b,'model:',r,'mean root',round(mean_root(math.log(T0/2/math.pi),b),3),'L',round(math.log(T0/2/math.pi),2),'time',round(time.time()-t0),flush=True)
print('mean roots at T0, L=26.89:',[(b,round(mean_root(26.89,b),2)) for b in (.5,.45,.4,.3,.2,.15,.1,.05,.03)])
print('DH: L for gamma=85.70, q=5: log(5 g/2pi)=',math.log(5*85.70/2/math.pi),' break window <=3.6281; mean root beta .3085:',round(mean_root(math.log(5*85.70/2/math.pi),.3085),3))
