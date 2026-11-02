import sys,time,json; sys.path.insert(0,'/tmp/t4rev')
from lib import *
c=np.load('/tmp/t4rev/dh_c.npy'); 
nz=np.nonzero(c[1:])[0]+1; ln=np.log(nz.astype(float)); cc=c[nz]
d=json.load(open(ROOT+'data/wl_dh_zeros.json')); on=np.array([float(x) for x in d['online']])
off=[(0.30852,85.69935),(0.15083,114.16334),(0.07436,166.47931),(0.22426,176.70246),(0.36953,240.40467)]
def Oprime(wg,dl,L):
    k=ln<=L; wts=2*cc[k]*np.exp(-(0.5+dl)*ln[k])*(1-ln[k]/L)
    return np.log(5/np.pi)+arch_spline(wg,dl,L,1.5)-primesum(wg,ln[k],wts)
wg=np.arange(5,300,0.03)
for dl in (0.2,0.05,'1/L'):
  print('delta',dl,flush=True)
  for N in [20,30,40,60,100,200,500,1000,3000,1e4,3e4,1e5,3e5,1e6,2e6]:
    if N>3e5 and dl==0.2: continue
    L=np.log(N); dd=1/L if dl=='1/L' else dl
    O=Oprime(wg,dd,L); H0=float(np.real(H(0.,dd,L))); R=O/H0
    idx=np.nonzero((R[1:-1]>R[:-2])&(R[1:-1]>=R[2:]))[0]+1
    pk=wg[idx]; hp=R[idx]
    s=[]
    for e,g0 in off:
        j=np.argmin(np.abs(pk-g0)); s.append((round(float(pk[j]),2),round(float(hp[j]),2)))
    fl=pk[hp>=2]
    tall=[(round(float(x),2),round(float(h),1)) for x,h in zip(pk,hp) if h>=2 and min(abs(x-y[1]) for y in off)>0.5]
    print(f' N={N:.3g} L={L:.2f} H0={H0:.3f} mean={R.mean():.3f} min={R.min():.2f} sup={R.max():.2f} off={s} nge2={len(fl)} spur_n={len(tall)} spur={tall[:7]}',flush=True)
# shape at 85.7: exact pair kernel, delta=.05, L=log 2e6, side lobes
L=np.log(2e6);dd=.05;e,g0=off[0]
u=np.arange(-6,6,0.002); P=np.real(H(u-1j*e,dd,L)+H(u+1j*e,dd,L))/float(np.real(H(0.,dd,L)))
idx=np.nonzero((P[1:-1]>P[:-2])&(P[1:-1]>=P[2:]))[0]+1
print('pair kernel own eta=.3085 d=.05 L=14.5: local max (u,height)',[(round(float(u[i]),2),round(float(P[i]),2)) for i in idx if abs(u[i])<3],'min',round(float(P.min()),2),'R theory',P.max()/2)
print('interference at distance d from pair 240.40 eta=.3695, d=.05: L=10,14.5,20,30', [(LL,[round(float(np.real(H(dd_-1j*.36953,.05,LL)+H(dd_+1j*.36953,.05,LL))/np.real(H(0.,.05,LL))),3) for dd_ in (5.03,10.3,50.1)]) for LL in (10,14.5,20,30)])
