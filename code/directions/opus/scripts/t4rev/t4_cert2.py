import sys,time,json; sys.path.insert(0,'/tmp/t4rev')
from lib import *
from scipy.interpolate import CubicSpline
g=zeros8000()
_gx,_gw=np.polynomial.legendre.leggauss(16)
def Phi(x,L,b):
    G=np.where(x<L,(1-x/L),0.0)/np.cosh(x/2)
    return 2*(1-G)*np.exp(-b*x)/(-np.expm1(-2*x))
def corr(w,L,b):
    out=np.zeros_like(w)
    for i,ww in enumerate(w):
        Xm=L+40/b
        ed=np.unique(np.concatenate([np.linspace(0,L,max(int(L*max(ww,1)/2)+30,60)),np.linspace(L,Xm,max(int((Xm-L)*max(ww,1)/2)+30,60))]))
        x=(ed[:-1,None]+ed[1:,None])/2+(ed[1:,None]-ed[:-1,None])/2*_gx; wt=(ed[1:,None]-ed[:-1,None])/2*_gw
        out[i]=np.sum(wt*Phi(x,L,b)*np.cos(ww*x))
    return out
def archS(wg,L,b):
    co=np.arange(0,45.05,.1); cs=CubicSpline(co,corr(co,L,b))
    c=np.where(wg<=40,cs(np.minimum(wg,45)),-(1/8+(1-b)/L)/np.maximum(wg,1)**2)
    return np.real(digamma(b/2+.5j*wg))+c
def H0S(L):
    x=np.linspace(0,L,200001); return 2*np.trapezoid((1-x/L)/np.cosh(x/2),x)
for N in [1e4,2.2e4,2e6]:
    L=np.log(N); ns,lam=sieve_pp(int(N)); ln=np.log(ns)
    wts=2*lam*np.exp(-.5*ln)/np.cosh(ln/2)*(1-ln/L)
    dw=0.01 if N<1e6 else 0.02
    wg=np.arange(0.01,1500,dw)
    pole=4*(1-np.cos(wg*L))/(wg*wg*L)
    O=pole-np.log(np.pi)+archS(wg,L,.5)-primesum(wg,ln,wts); H0=H0S(L)
    R=O/H0
    first=wg[R>=2][0] if (R>=2).any() else None
    print(f'sech zeta N={N:.3g} L={L:.3f} #pp={len(ns)} H0={H0:.4f} first w with O>=2H0: {first} sup[0,1000]={R[wg<1000].max():.4f} at {wg[wg<1000][R[wg<1000].argmax()]:.2f} inf={R.min():.3f} sup on [0,first)={(R[wg<first].max() if first else R.max()):.4f}',flush=True)
    # zero side check at few w, Lipschitz constant
    w0=np.array([14.1347,100.,1000.])
    # numeric H for sech
    x=np.linspace(0,L,200001); Gx=(1-x/L)/np.cosh(x/2)
    zs=[2*np.trapezoid(Gx*np.cos(np.outer(w,x))).sum() for w in []]
    print('  sum of |wts*log n| (Lipschitz prime part)',float((wts*ln).sum()),' max slope on grid',float(np.abs(np.diff(O)).max()/dw))
# Direct check at the exp kernel L=log 2e6 on [0,1500] with primes to 2e4 kernel exp
for N in [2.2e4]:
    L=np.log(N); ns,lam=sieve_pp(int(N)); wg=np.arange(0.01,700,0.01)
    O=O_zeta_prime(wg,.5,L,ns,lam); H0=float(np.real(H(0.,.5,L))); R=O/H0
    print('exp zeta N',N,'H0',H0,'first',wg[R>=2][0] if (R>=2).any() else None)
# DH with the sech kernel (positivity needs no zero with Re>=1 in range)
c=np.load('/tmp/t4rev/dh_c.npy'); nz=np.nonzero(c[1:])[0]+1; ln=np.log(nz.astype(float)); cc=c[nz]
for N in [1e3,1e4,1e5,3e5]:
    L=np.log(N);k=ln<=L; wts=2*cc[k]*np.exp(-.5*ln[k])/np.cosh(ln[k]/2)*(1-ln[k]/L)
    wg=np.arange(1,300,0.02); O=np.log(5/np.pi)+archS(wg,L,1.5)-primesum(wg,ln[k],wts);R=O/H0S(L)
    idx=np.nonzero((R[1:-1]>R[:-2])&(R[1:-1]>=R[2:]))[0]+1
    print('sech DH N',N,'H0',round(H0S(L),3),'peaks>=2:',[(round(float(wg[i]),2),round(float(R[i]),2)) for i in idx if R[i]>=2][:14],flush=True)
