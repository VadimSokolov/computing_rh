import numpy as np
from scipy.special import digamma
ROOT='./'
def zeros8000():
    return np.array([float(l) for l in open(ROOT+'data/zeros_8000.txt') if l.strip()])
def I_(z,L):
    # int_0^L e^{-zx}(1-x/L)dx, analytic in z
    z=np.asarray(z,complex)
    small=np.abs(z*L)<1e-4
    zz=np.where(small,1.0,z)
    v=(1-np.exp(-zz*L))/zz-(1-np.exp(-zz*L)*(1+zz*L))/(zz*zz*L)
    return np.where(small,L/2-z*L*L/6,v)
def H(r,d,L):
    r=np.asarray(r,complex)
    return I_(d-1j*r,L)+I_(d+1j*r,L)
def sieve_pp(N):
    s=np.ones(N+1,bool); s[:2]=False
    for p in range(2,int(N**.5)+1):
        if s[p]: s[p*p::p]=False
    ps=np.nonzero(s)[0]
    ns=[ps.astype(float)]; lam=[np.log(ps.astype(float))]
    for p in ps[ps<=int(N**.5)]:
        q=int(p)*int(p)
        while q<=N:
            ns.append(np.array([float(q)])); lam.append(np.array([np.log(float(p))])); q*=int(p)
    return np.concatenate(ns),np.concatenate(lam)
_gx,_gw=np.polynomial.legendre.leggauss(16)
def arch_corr(w,d,L,b):
    """(1/L)int_0^L 2x e^{-(d+b)x}cos(wx)/(1-e^{-2x}) dx + int_L^inf 2e^{-(d+b)x}cos(wx)/(1-e^{-2x}) dx"""
    w=np.atleast_1d(np.asarray(w,float)); out=np.zeros_like(w)
    for i,ww in enumerate(w):
        if ww>40:
            out[i]=(1/L)*(-(1-(d+b)))/ww**2*(1.0)  # -phi'(0)/w^2/L, phi'(0)=1-(d+b)
            continue
        Xmax=max(L,L+60/(d+b))
        np_=max(int(Xmax*max(ww,1)/2)+40,60)
        ed=np.linspace(0,Xmax,np_+1)
        x=(ed[:-1,None]+ed[1:,None])/2+(ed[1:,None]-ed[:-1,None])/2*_gx
        wt=(ed[1:,None]-ed[:-1,None])/2*_gw
        den=-np.expm1(-2*x)
        wgt=np.where(x<=L,2*x/L,2.0)*np.exp(-(d+b)*x)/den
        # panel edge at L: add L as a node boundary approximately by fine panels
        out[i]=np.sum(wt*wgt*np.cos(ww*x))
    return out
def arch(w,d,L,b):
    w=np.atleast_1d(np.asarray(w,float))
    return np.real(digamma(b/2+d/2+0.5j*w))+arch_corr(w,d,L,b)
def arch_spline(wg,d,L,b):
    from scipy.interpolate import CubicSpline
    lo,hi=wg.min(),wg.max()
    if hi<=40 or lo>40:
        pass
    co=np.arange(max(lo-1,0),min(hi,45)+1.01,0.1) if lo<45 else np.array([])
    out=np.real(digamma(b/2+d/2+0.5j*wg))
    corr=np.zeros_like(wg)
    m=wg>40
    corr[m]=(1/L)*(-(1-(d+b)))/wg[m]**2
    if len(co):
        cs=CubicSpline(co,arch_corr(co,d,L,b))
        corr[~m]=cs(wg[~m])
    return out+corr
def primesum(wg,ns_log,wts,chunk=400):
    out=np.empty_like(wg)
    for i in range(0,len(wg),chunk):
        out[i:i+chunk]=np.cos(np.outer(wg[i:i+chunk],ns_log)) @ wts
    return out
def O_zeta_prime(wg,d,L,ns,lam):
    keep=np.log(ns)<=L
    ln=np.log(ns[keep]); 
    wts=2*lam[keep]*np.exp(-(0.5+d)*ln)*(1-ln/L)
    pole=2*np.real(I_(d-0.5-1j*wg,L)+I_(d+0.5-1j*wg,L))
    return pole-np.log(np.pi)+arch_spline(wg,d,L,0.5)-primesum(wg,ln,wts)
def O_zeros(wg,d,L,g,tail_density=lambda x:np.log(x/(2*np.pi))/(2*np.pi),chunk=200):
    out=np.empty_like(wg)
    for i in range(0,len(wg),chunk):
        ww=wg[i:i+chunk][:,None]
        out[i:i+chunk]=np.real(H(ww-g,d,L)+H(ww+g,d,L)).sum(axis=1)
    T=g[-1]
    # tail: H ~ 2(d+1/L)/r^2, integrate density
    xs=T*np.exp(np.linspace(0,12,4001)); fx=tail_density(xs)*2*(d+1/L)*(1/(xs-wg[:,None])**2+1/(xs+wg[:,None])**2)
    out+=np.trapezoid(fx*xs,np.linspace(0,12,4001),axis=1)
    return out
