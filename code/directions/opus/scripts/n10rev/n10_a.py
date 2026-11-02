import sys; sys.path.insert(0,'/tmp/n10rev')
from core import *
import time, json
ns,Ls=pp(2000000); ln=np.log(ns); u=Ls/np.sqrt(ns)
print('check Psi(8)=',Parch(8.)-comb(np.array([8.]),ln,u),' arch(8)=',Parch(8.))
# --- arch break
print('A. arch first non PD by Levinson')
for inv in [200,1000,2000,5000,10000,20000]:
    d=1/inv; t=np.arange(int(1.0/d)+2)*d; P=Parch(t); P[0]=0
    K,V,k=lev(incr(P)); print(inv,'last PD',k*d,'first nonPD',(k+1)*d)
# --- delays
PP=[n for n in ns if n<=37]; 
def fit(x,y):
    A=np.vstack([np.ones_like(x),x]).T; b,_,_,_=np.linalg.lstsq(A,y,rcond=None); r=y-A @ b
    R2=1-(r @ r)/((y-y.mean()) @ (y-y.mean())); s2=(r @ r)/(len(y)-2); se=np.sqrt(np.diag(s2*np.linalg.inv(A.T @ A)))
    return b,se,R2
print('B. delays, omit n and larger')
res={}
for inv in [200,800,3200]:
    d=1/inv; rows=[]
    for idx,n in enumerate(PP):
        l=math.log(n); t=np.arange(int((l+.35)/d)+2)*d
        P=Parch(t)-comb(t,ln[:idx],u[:idx]); P[0]=0
        K,V,k=lev(incr(P))
        rows.append((n,u[idx],k*d-l,(k+1)*d-l))
    rows=np.array(rows); res[inv]=rows
    wn=rows[:,1]; nn=rows[:,0]
    for name,dd in [('lastPD',rows[:,2]),('firstNonPD',rows[:,3])]:
        pos=dd>0
        b,se,R2=fit(np.log(nn[pos]),np.log(wn[pos]*dd[pos]))
        b2,se2,R22=fit(np.log(nn[pos]),np.log(dd[pos]))
        A=np.vstack([np.ones(pos.sum()),np.log(wn[pos]),np.log(nn[pos])]).T
        cc,_,_,_=np.linalg.lstsq(A,np.log(dd[pos]),rcond=None); r=np.log(dd[pos])-A @ cc; R23=1-(r @ r)/np.sum((np.log(dd[pos])-np.log(dd[pos]).mean())**2)
        print(f'd=1/{inv} {name} N={pos.sum()} w d = {np.exp(b[0]):.4f} n^{b[1]:.3f} (se {se[1]:.3f}) R2={R2:.3f}; log d = {cc[0]:.2f}{cc[1]:+.2f}log w{cc[2]:+.2f}log n R2={R23:.3f}; d only: slope {b2[1]:.2f} R2={R22:.3f}')
    print(' delays lastPD',np.round(rows[:,2],4).tolist()); print(' delays first',np.round(rows[:,3],4).tolist()); print(' n',rows[:,0].astype(int).tolist())
    print(' frac of lastPD delay below step:',np.mean(rows[:,2]<d),' spread of d/step',np.round(rows[:,3]/d,1).tolist())
# 1/w ordering
r=res[200]; print('rank corr d vs 1/w, d=1/200:',np.corrcoef(np.argsort(np.argsort(r[:,3])),np.argsort(np.argsort(1/r[:,1])))[0,1])
r=res[3200]; print('rank corr d vs 1/w, d=1/3200:',np.corrcoef(np.argsort(np.argsort(r[:,3])),np.argsort(np.argsort(1/r[:,1])))[0,1])
# DH control
print('C. DH control')
dn,dL=dh(); dl=np.log(dn); du=dL/np.sqrt(dn)
for inv in [200,800,3200]:
    d=1/inv; rows=[]
    for idx in range(len(dn)):
        n=dn[idx]
        if n>32: break
        l=math.log(n); t=np.arange(int((l+.5)/d)+2)*d
        P=Pdh_arch(t)-comb(t,dl[:idx],du[:idx]); P[0]=0
        K,V,k=lev(incr(P))
        if k is None: rows.append((n,du[idx],np.nan,np.nan)); continue
        rows.append((n,du[idx],k*d-l,(k+1)*d-l))
    rows=np.array(rows); ok=~np.isnan(rows[:,3])&(rows[:,2]>0)
    b,se,R2=fit(np.log(rows[ok,0]),np.log(np.abs(rows[ok,1])*rows[ok,2]))
    b1,se1,R21=fit(np.log(rows[ok,0]),np.log(np.abs(rows[ok,1])*rows[ok,3]))
    print(f'DH d=1/{inv} N={ok.sum()} lastPD |w|d={np.exp(b[0]):.4f} n^{b[1]:.3f} se {se[1]:.3f} R2={R2:.3f}; first: {np.exp(b1[0]):.4f} n^{b1[1]:.3f} R2={R21:.3f}')
    print(' n',rows[:,0].astype(int).tolist()); print(' u',np.round(rows[:,1],3).tolist()); print(' d',np.round(rows[:,3],4).tolist())
# DH break at .025 and zeta long
print('D. long grid .025')
d=.025; t=np.arange(int(16.2/d)+2)*d
P=Pdh_arch(t)-comb(t,dl,du); P[0]=0
K,V,k=lev(incr(P)); print('DH (n<=3000) d=.025 last PD',k*d,'first nonPD',(k+1)*d, 'kappa last 4',np.round(K[-4:],3),'min(1-k^2) before last 3',np.min(1-K[1:-3]**2),'at',1+np.argmin(1-K[1:-3]**2),'v/v0 at 4,5',V[int(4/d)]/V[0],V[int(5/d)]/V[0])
