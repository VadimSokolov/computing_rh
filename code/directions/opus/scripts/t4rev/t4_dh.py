import sys,time,json; sys.path.insert(0,'/tmp/t4rev')
from lib import *
t0=time.time()
kap=0.2840790438404122960282918
NM=int(float(sys.argv[1])) if len(sys.argv)>1 else 20_000_000
per=np.array([0.,1.,kap,-kap,-1.])
a=per[np.arange(NM+1)%5]; a[0]=0
c=a*np.log(np.maximum(np.arange(NM+1),1).astype(float))
for n in range(1,NM//2+1):
    cn=c[n]
    if cn!=0.0:
        q=NM//n
        c[2*n::n]-=cn*a[2:q+1]
print('c done',time.time()-t0,flush=True)
print('c(1..12)',np.round(c[1:13],5).tolist())
nz=np.nonzero(c[1:])[0]+1; ln=np.log(nz.astype(float)); cc=c[nz]
print('nonzero',len(nz),'max|c|n^-.5',np.abs(cc*np.exp(-.5*ln)).max())
# fit growth
for lo,hi in [(1e3,1e4),(1e5,1e6),(2e6,2e7)]:
    m=(nz>=lo)&(nz<hi)
    if m.any(): print('rms c(n)/sqrt n in',lo,hi,float(np.sqrt(np.mean((cc[m]*np.exp(-.5*ln[m]))**2))))
d=json.load(open(ROOT+'data/wl_dh_zeros.json'))
on=np.array([float(x) for x in d['online']]); print('n online <200',len(on),'offline',d['offline'][:5] if isinstance(d['offline'],list) else d['offline'])
np.save('/tmp/t4rev/dh_c.npy',c[:2000001])
off=[(0.30852,85.69935),(0.15083,114.16334),(0.07436,166.47931),(0.22426,176.70246),(0.36953,240.40467),(0.26250,983.75219),(0.19140,1012.01998),(0.19810,1018.79576),(0.08614,1029.00469)]
def Oprime(wg,dl,L):
    k=ln<=L
    wts=2*cc[k]*np.exp(-(0.5+dl)*ln[k])*(1-ln[k]/L)
    return np.log(5/np.pi)+arch_spline(wg,dl,L,1.5)-primesum(wg,ln[k],wts)
def Ozero(wg,dl,L,T=200.):
    O=np.empty_like(wg)
    for i in range(0,len(wg),200):
        ww=wg[i:i+200][:,None]
        O[i:i+200]=np.real(H(ww-on,dl,L)+H(ww+on,dl,L)).sum(axis=1)
    return O
# check: explicit formula at low height with offline set; compare near w=20 etc, far tail via mean density
L=np.log(2e6)
for dl in (0.5,0.05):
    wt=np.array([20.,50.,85.7,114.16,150.])
    p=Oprime(wt,dl,L)
    z=Ozero(wt,dl,L)
    for e,g0 in off:
        z=z+np.real(H(wt-g0-1j*e,dl,L)+H(wt-g0+1j*e,dl,L)+H(wt+g0-1j*e,dl,L)+H(wt+g0+1j*e,dl,L))
    print('EF check dl',dl,'prime',np.round(p,4).tolist(),'zero(<200 listed, offline listed, no other tail)',np.round(z,4).tolist())
# detection sweep
for dl in (0.5,0.2,0.05,'1/L'):
  print('delta',dl)
  for N in [10,20,30,40,60,100,200,500,1000,3000,1e4,3e4,1e5,3e5,1e6,3e6,1e7,NM]:
    L=np.log(N); dd=1/L if dl=='1/L' else dl
    wg=np.arange(5,300,0.02)
    O=Oprime(wg,dd,L); H0=float(np.real(H(0.,dd,L))); R=O/H0
    idx=np.nonzero((R[1:-1]>R[:-2])&(R[1:-1]>=R[2:]))[0]+1
    pk=wg[idx]; hp=R[idx]
    s=[]
    for e,g0 in off[:5]:
        j=np.argmin(np.abs(pk-g0)); s.append((round(float(pk[j]),2),round(float(hp[j]),2)))
    # all with height >=2
    fl=pk[hp>=2]; 
    onl=np.concatenate([on,[x[1] for x in off[:5]]])
    tall_spur=[round(float(x),2) for x in fl if min(abs(x-y[1]) for y in off[:5])>0.5]
    neg=float(R.min())
    print(f' N={N:.3g} L={L:.2f} H0={H0:.3f} mean={R.mean():.3f} min={neg:.2f} sup={R.max():.2f} near(85.70,114.16,166.48,176.70,240.40)={s} n>=2:{len(fl)} tall not at offline:{tall_spur[:12]}',flush=True)
print(time.time()-t0)
