import sys,time; sys.path.insert(0,'/tmp/t4rev')
from lib import *
g=zeros8000(); print('gamma_8000',g[-1])
print('planted pair')
for T0,W in [(100,50),(1000,100),(6000,100)]:
    z=g[(g>T0)&(g<T0+W)]; k=np.argmax(np.diff(z)); g0=(z[k]+z[k+1])/2
    sel=g[(g>T0-90)&(g<T0+W+90)]
    for L in [14.5,25,50]:
        d=1/L; H0=float(np.real(H(0.,d,L)))
        wg=np.arange(g0-1,g0+1,0.002); base=np.real(H(wg[:,None]-sel,d,L)).sum(axis=1)
        out=[]
        for eta in [0,0.003,0.01,0.03,0.05,0.1,0.2,0.3]:
            own=np.real(H(wg-g0-1j*eta,d,L)+H(wg-g0+1j*eta,d,L)); out.append((eta,round(float(((base+own)/H0).max()),3)))
        hz=np.real(H(z[:,None]-sel,d,L)).sum(axis=1)/H0
        # height of an on-line close pair at gap .05
        print(T0,'g0',round(float(g0),2),'L',L,'single heights mean %.3f sd %.3f max %.3f'%(hz.mean(),hz.std(),hz.max()),'pair',out,'2+eta^2L^2/6 at .03:',round(2+.0009*L*L/6,3))
print('prime side at height, primes<2e6')
NM=2_000_000; ns,lam=sieve_pp(NM); L=np.log(NM)
for T0,W in [(1000,60),(3000,60),(6000,60),(8000,60)]:
    d=1/L; wg=np.arange(T0,T0+W,0.02)
    O=O_zeta_prime(wg,d,L,ns,lam); H0=float(np.real(H(0.,d,L))); R=O/H0
    idx=np.nonzero((R[1:-1]>R[:-2])&(R[1:-1]>=R[2:]))[0]+1
    pk=wg[idx][R[idx]>.45]; z=g[(g>T0+1)&(g<T0+W-1)]
    err=np.array([np.min(np.abs(pk-x)) for x in z]); sp=[round(float(p),2) for p in pk if T0+1<p<T0+W-1 and np.min(np.abs(z-p))>.25]
    print(T0,'zeros',len(z),'peaks>.45',len(pk),'found within .1:',(err<.1).sum(),'med err',np.median(err),'mean',R.mean(),'log/L',np.log(T0/2/np.pi)/L,'sup',R.max(),'spurious',sp[:6],'zero-side diff max',float(np.abs(O-O_zeros(wg,d,L,g)).max()) if T0+W<g[-1]-500 else None,flush=True)
# gaps table for law
gaps=np.diff(g); s=gaps*np.log(g[:-1]/2/np.pi)/(2*np.pi)
for a,b in [(14,150),(150,1000),(1000,3000),(3000,g[-1])]:
    m=(g[:-1]>=a)&(g[:-1]<b); print('range',a,b,'n',m.sum(),'min gap',gaps[m].min(),'min normalised',s[m].min(),'10pct gap norm',np.quantile(s[m],.1),'2pi/mingap',2*np.pi/gaps[m].min())
