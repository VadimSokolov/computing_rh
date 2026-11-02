import sys,time; sys.path.insert(0,'/tmp/t4rev')
from lib import *
from scipy.optimize import minimize_scalar
g=zeros8000(); t0=time.time()
NM=2_000_000; ns,lam=sieve_pp(NM); L=np.log(NM)
print('L',L,'n prime powers',len(ns),'zeros<150',(g<150).sum(),'zeros in [5,160]',(g<160).sum())
for d in [0.5,0.05]:
    wt=np.array([14.1347,20.,50.,100.])
    p=O_zeta_prime(wt,d,L,ns,lam); z=O_zeros(wt,d,L,g)
    print('check d',d,[ (float(a),round(float(b),6),round(float(c),6),float(b-c)) for a,b,c in zip(wt,p,z)])
for d in [0.5,0.2,0.1,0.05,1/L]:
    H0=float(np.real(H(0.,d,L))); wg=np.arange(5,160,0.01)
    O=O_zeta_prime(wg,d,L,ns,lam)
    idx=np.nonzero((O[1:-1]>O[:-2])&(O[1:-1]>=O[2:]))[0]+1
    allpk=len(idx); pk=idx[O[idx]>0.25*H0]
    # parabolic refine
    a,b,c=O[pk-1],O[pk],O[pk+1]; est=wg[pk]+0.01*0.5*(a-c)/(a-2*b+c)
    hts=(b-0.125*(a-c)*(a-c)/(a-2*b+c)*0)/H0
    tr=g[g<150]; err=np.array([np.min(np.abs(est-x)) for x in tr])
    tr58=g[g<160]; err58=np.array([np.min(np.abs(est-x)) for x in tr58])
    sp=[round(float(e),3) for e in est if np.min(np.abs(tr58-e))>0.5]
    print(f'd={d:.4f} H0={H0:.4f} all_max={allpk} peaks>H0/4={len(pk)} max_err={err.max():.5f} med_err={np.median(err):.5f} heights {hts.min():.4f} to {hts.max():.4f} sup_O/H0_on[5,160]={O.max()/H0:.4f} spurious={sp} min_next_above_threshold_sidelobe={np.sort(O[idx])[:0]}')
    # highest non zero local max
    nz=[O[i]/H0 for i in idx if np.min(np.abs(tr58-wg[i]))>0.5]
    print('   highest non-zero local max /H0',max(nz) if nz else None,' lowest zero peak',hts.min())
print(time.time()-t0)
