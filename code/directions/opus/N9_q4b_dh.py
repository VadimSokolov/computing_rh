# Q4, second part: negative inertia on longer windows, aliasing, and a search for the zeros it points at.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
N=int(math.exp(14.2))+2; Lam=dh_coeffs(N)
def top_peaks(u,d,k=3):
    n=1<<int(math.ceil(math.log2(16*len(u)+1)))
    F=np.abs(np.fft.rfft(u,n))
    idx=np.nonzero((F[1:-1]>F[:-2])&(F[1:-1]>F[2:]))[0]+1
    idx=idx[np.argsort(F[idx])[::-1]][:k]
    return [round(float(2*math.pi*i/n/d),2) for i in idx]
for d in [0.025,0.02,0.0175]:
    M=int(14.0/d); c=cov_dh(d,M,Lam)
    print('delta=%g, sampling frequency %.3f'%(d,2*math.pi/d))
    for a in [4,5,5.5,6,7,8,9,10,12,14]:
        m=int(round(a/d))
        if m>M: continue
        w,V=np.linalg.eigh(toeplitz(c[:m]))
        neg=np.nonzero(w<0)[0]
        print('  a=%g m=%d negative eigenvalues %d; most negative %.3e; top: %s'%(a,m,len(neg),w[0],[ (round(float(w[i]),5),top_peaks(V[:,i],d,2)) for i in neg[:8]]))
print('time %.1f'%(time.time()-t0))
# zeros of f off the line for 200<gamma<300 by Newton from seeds
ctx.prec=128
ka=arb(10-2*arb(5).sqrt()).sqrt(); ka=(ka-2)/(arb(5).sqrt()-1)
ar=[None,arb(1),ka,-ka,arb(-1)]; hur=[None]+[acb(arb(r)/5) for r in range(1,5)]; l5=arb(5).log()
def f(s): return sum((acb(ar[r])*s.zeta(hur[r]) for r in range(1,5)),acb(0))*(-s*acb(l5)).exp()
def cf(s):
    v=f(acb(s.real,s.imag)); return complex(float(v.real.mid()),float(v.imag.mid()))
found=[]
def newton(s0):
    s=s0
    for it in range(25):
        h=1e-7; fs=cf(s); df=(cf(s+h)-cf(s-h))/(2*h)
        if df==0: return None
        ds=fs/df
        if abs(ds)>1.0: ds=ds/abs(ds)
        s=s-ds
        if abs(ds)<1e-11: return s if abs(cf(s))<1e-8 else None
    return None
for lo,hi in [(200,300),(960,1040)]:
    for g in np.arange(lo,hi,1.0):
        for b in [0.55,0.65,0.75,0.85,0.95,1.05,1.15]:
            r=newton(complex(b,g))
            if r is None or r.real<0.5005 or not (lo-1<=r.imag<=hi+1): continue
            if not any(abs(r-q)<1e-6 for q in found): found.append(r)
found.sort(key=lambda z:z.imag)
print('off line zeros found by Newton in 200<gamma<300 and 960<gamma<1040 (seeded, not proved complete):')
for z in found: print('  beta=%.6f gamma=%.6f (Re s=%.6f)'%(z.real-0.5,z.imag,z.real), ' alias at delta=.025: %.2f, at .02: %.2f, at .0175: %.2f'%tuple(abs(((z.imag+math.pi/d)%(2*math.pi/d))-math.pi/d) for d in (0.025,0.02,0.0175)))
print('total %.1f'%(time.time()-t0))
