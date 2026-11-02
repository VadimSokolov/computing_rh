import sys,time; sys.path.insert(0,'/tmp/t4rev')
from lib import *
g=zeros8000()
# Kernels positive on the whole strip: exp: G=e^{-|x|/2}(1-|x|/L); sech: G=sech(x/2)(1-|x|/L)
# Zero side for all zeros on line; H kernel numerically via fine FFT-free quadrature in x.
def kern(name,L):
    x=np.linspace(0,L,60001)
    if name=='exp': Gx=np.exp(-x/2)
    elif name=='sech': Gx=1/np.cosh(x/2)
    return x,Gx*(1-x/L)
def Hgrid(name,L,rmax=60,dr=0.01):
    x,G=kern(name,L); dx=x[1]-x[0]
    r=np.arange(0,rmax,dr)
    v=np.empty_like(r); w=np.full_like(x,dx); w[0]=w[-1]=dx/2
    for i in range(0,len(r),500):
        v[i:i+500]=2*(np.cos(np.outer(r[i:i+500],x))*(w*G)).sum(axis=1)
    return r,v
for name in ['exp','sech']:
  for L in [10,12,14.5,18,25,40]:
    r,v=Hgrid(name,L); H0=v[0]; c=(v[-1]*r[-1]**2)
    dw=0.02; wg=np.arange(0.5,g[-1]-5,dw)
    O=np.zeros_like(wg)
    # convolve zeros by binned add: for each zero add kernel on grid
    idx=np.round((g)/dw).astype(int)
    tip=np.arange(-int(60/dw),int(60/dw)+1)
    ker=np.interp(np.abs(tip*dw),r,v)
    # exact placement by nearest grid node: error dw/2 in position
    dep=np.zeros(len(wg)+2*len(tip));
    pos=np.round((g-wg[0])/dw).astype(int)+len(tip)
    np.add.at(dep,pos[pos<len(dep)],1.0)
    from scipy.signal import fftconvolve
    O=fftconvolve(dep,ker,mode='same')[len(tip):len(tip)+len(wg)]
    # far tail by mean: two sided beyond 60: c/r^2*dens
    dens=np.log(np.maximum(wg,10)/(2*np.pi))/(2*np.pi)
    O+=2*c*dens/60
    N=O/H0; thr=2.0
    # first crossing
    m=N>=thr
    first=wg[m][0] if m.any() else None
    sups=[(T,round(float(N[wg<=T].max()),3)) for T in (160,500,1000,3000,9000)]
    print(name,'L',L,'H0',round(H0,4),'2H0',round(2*H0,4),'first w with O>=2H0:',first,'sup/H0 up to T:',sups)
