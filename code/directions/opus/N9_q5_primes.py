# Q5: every prime power is indispensable: truncation rows by delta, single omission, tolerance, forecast, spectrum.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
ns,Ls=prime_powers(1200000); ln=np.log(ns); u=Ls/np.sqrt(ns); w=2*u
G=load_zeros(8000)
ds=[0.025,0.0125,0.00625,0.003125,0.0015625,0.00078125]
amax=7.6; arch={d:cov_arch_zeta(d,int(amax/d)+2) for d in ds}; print('arch ready %.1f'%(time.time()-t0))
def sigma(T0): return 2*(math.log(T0/(2*math.pi))+1)/(math.pi*T0)
def model(n0,d,idx):
    T0=min(math.pi/d,4*math.pi*n0); S=sigma(T0); l0=ln[idx]
    lo,hi=0.0,2.0
    for _ in range(60):
        s=0.5*(lo+hi); q=(ln[idx:]<=l0+s)
        if np.sum(w[idx:][q]*(l0+s-ln[idx:][q]))>S: hi=s
        else: lo=s
    return hi
print('A. truncation rows: keep prime powers n<=N, omit all others. delay=first non PD window - log n0 [model]; last PD window in the note\'s convention')
for N in [1,2,3,4,5,7,8,9,11,13,16,17,31,100,1000]:
    i0=int(np.searchsorted(ns,N,side='right')); n0=int(ns[i0]); l0=ln[i0]
    out=[]
    for d in ds:
        M=min(len(arch[d]),int((l0+0.7)/d)+2); c=arch[d][:M].copy()
        add_primes(c,d,ln[:i0],u[:i0])
        ks,vs,f=levinson(c)
        if f is None: out.append('d=%.5f none'%d); continue
        out.append('%.5f:%.4f[%.4f]'%(d,(f+1)*d-l0,model(n0,d,i0)))
    print(' N=%d n0=%d log n0=%.4f w=%.4f | '%(N,n0,l0,w[i0])+' '.join(out))
print('time %.1f'%(time.time()-t0))
print('\nB. omit one prime power, keep all others, windows to 12.05: first non PD window - log n0')
base={}
for d in [0.025,0.00625]:
    M=int(12.05/d)+1; c=cov_arch_zeta(d,M); add_primes(c,d,ln,u); base[d]=c
def spike(d,M,l,uu):
    e=np.zeros(M); x=l/d; m0=int(math.floor(x)); eta=x-m0
    if m0<M: e[m0]=0.5*uu*(1-eta)*d
    if m0+1<M: e[m0+1]=0.5*uu*eta*d
    return e
lst=[2,3,4,5,7,8,9,11,13,16,17,25,27,32,64,101,256,1009,2048,10007,30011,50021,70001,100003,120011,130003,140009,150001]
for n0 in lst:
    i=int(np.searchsorted(ns,n0)); assert ns[i]==n0
    row=[]
    for d in [0.025,0.00625]:
        cc=base[d]+spike(d,len(base[d]),ln[i],u[i]); ks,vs,f=levinson(cc)
        row.append('%g: %s'%(d,('%.4f'%((f+1)*d-ln[i])) if f is not None else 'none to 12.05'))
    print(' n0=%d log=%.3f w=%.4f | %s | 4v/w-type model at T0=pi/d: %.3f'%(n0,ln[i],w[i],' ; '.join(row),sigma(math.pi/0.025)/w[i]))
print('time %.1f'%(time.time()-t0))
print('\nC. tolerance on the grid delta=0.025: T + eps E stays PD for eps in (lo,hi), Lambda(n0) -> (1+eps) Lambda(n0)')
from scipy.linalg import eigh as geigh
for a in [6.0,12.0]:
    m=int(round(a/0.025)); T=toeplitz(base[0.025][:m])
    print(' a=%g:'%a)
    for n0 in [2,3,4,5,7,8,9,11,13,16,101,1009,10007,100003]:
        i=int(np.searchsorted(ns,n0))
        if ln[i]>=a: continue
        E=toeplitz(-spike(0.025,m,ln[i],u[i]))
        mu=geigh(E,T,eigvals_only=True,subset_by_index=[0,m-1])
        print('   n0=%d: eps in (%.3e, %.3e)'%(n0,-1/mu.max(),-1/mu.min()))
print('time %.1f'%(time.time()-t0))
print('\nD. forecast: with primes below n0, what does the predictor expect at the lag of the spike of n0?')
print(' n0, delta, k*, spike size, c at k*, predicted c at k*, kappa at k*, naive spike/v, largest |kappa| within 3 lags and its lag, v(+)/v(-) two lags after, first non PD window without n0..')
for d in [0.025,0.00625]:
    for n0 in [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]:
        i=int(np.searchsorted(ns,n0)); l=ln[i]
        M=int((l+0.6)/d)+2
        A=arch[d][:M].copy(); add_primes(A,d,ln[:i],u[:i]); cm=A.copy(); cp=A.copy(); add_primes(cp,d,ln[i:i+1],u[i:i+1])
        kp,vp,fp=levinson(cp[:int((l+0.2)/d)+2]); km,vm,fm=levinson(cm)
        x=l/d; m0=int(math.floor(x)); eta=x-m0; ks=m0+1 if eta>=0.5 else m0
        sp=0.5*u[i]*(eta if eta>=0.5 else 1-eta)*d
        pred=cp[ks]-kp[ks]*vp[ks-1]
        win=np.abs(kp[ks-3:ks+4]); jm=ks-3+int(np.argmax(win))
        j2=min(ks+2,len(vm)-1,len(vp)-1)
        print(' %d %g k*=%d eta=%.2f spike=%.5f c=%.5f pred=%.5f kappa=%.3f naive=%.3f max|k| %.3f at lag %+d; v+/v- %.3f; no-n0 break +%.4f'%(n0,d,ks,eta,sp,cp[ks],pred,kp[ks],sp/vp[ks-1],win.max(),jm-ks,vp[j2]/vm[j2],((fm+1)*d-l) if fm is not None else float('nan')))
print('time %.1f'%(time.time()-t0))
print('\nE. approach to the break, archimedean only, delta=0.001 (odd sector): a, max|kappa| last 10 lags, min(1-kappa^2) so far, v/v0, 2 lam_min/delta')
d=0.001; c=cov_arch_zeta(d,760); ks,vs,f=levinson(c); print(' break at k=',f,' window',(f+1)*d)
for a in [0.3,0.5,0.6,0.65,0.69,0.7,0.72,0.73,0.735,0.74,0.741,0.742,0.743]:
    m=int(round(a/d)); print('  a=%.3f max|k|=%.3f min(1-k^2)=%.4f v/v0=%.4f eps_odd=%.4e'%(a,np.abs(ks[max(1,m-10):m]).max(),(1-ks[1:m]**2).min(),vs[m-1]/vs[0],2*lam_min(c,m,parity='o')/d))
print('\nF. spectrum: Hann lag window on [0,12], peaks against the 13 ordinates below 60')
d=0.025; Mw=480; c=base[d]; hann=0.5*(1+np.cos(np.pi*np.arange(Mw+1)/Mw))
r=np.arange(5,70,0.002); lags=np.arange(1,Mw+1)*d
S=c[0]+2*np.sum((hann[1:]*c[1:Mw+1])[None,:]*np.cos(r[:,None]*lags[None,:]),axis=1)
q=(ln<=12)
Sp=-np.array([np.sum(u[q]*0.5*(1+np.cos(np.pi*ln[q]/12))*np.cos(x*ln[q])) for x in r])
def report(S,name):
    idx=np.nonzero((S[1:-1]>S[:-2])&(S[1:-1]>S[2:]))[0]+1; pk=r[idx]; gg=G[G<60]
    dist=[pk[np.argmin(np.abs(pk-g))]-g for g in gg]
    print(' %s: peaks in (5,60): %d; max|distance|=%.4f; distances %s'%(name,int(np.sum(pk<60)),np.max(np.abs(dist)),np.round(dist,3)))
report(S,'covariance of the process'); report(Sp,'primes only, no process, same window')
print('total %.1f'%(time.time()-t0))
