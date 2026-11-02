import sys; sys.path.insert(0,'/tmp/n10rev')
from core import *
from scipy.linalg import eigh, toeplitz
from numpy.polynomial import legendre as Lg
ns,Ls=pp(11000000); ln=np.log(ns); u=Ls/np.sqrt(ns)
ROOT='.'
# 1 mechanism test
print('1. mechanism: margin at log n, w d, overlap at break')
for inv in [200,800]:
    d=1/inv
    for n in [3,5,7,13,16,32,37]:
        idx=int(np.searchsorted(ns,n)); l=ln[idx]; w=u[idx]
        t=np.arange(int((l+.35)/d)+2)*d
        Pt=Parch(t)-comb(t,ln,u); Pt[0]=0
        Pn=Parch(t)-comb(t,ln[:idx],u[:idx]); Pn[0]=0
        Kn,Vn,k=lev(incr(Pn)); m=k+1   # first non PD size
        ct=incr(Pt); cn=incr(Pn)
        lt=eigh(toeplitz(ct[:m]),eigvals_only=True,subset_by_index=[0,0])[0]*2/d
        ev,vv=eigh(toeplitz(cn[:m]),subset_by_index=[0,0]); f=vv[:,0]/np.sqrt(d)  # unit L2 step fn
        s=int(round(l/d)); Rl=sum(d*f[i]*f[i+s] for i in range(m-s)); dbr=m*d-l
        endv=max(abs(f[0]),abs(f[-1]))
        print(f'd=1/{inv} n={n} delay={dbr:.4f} w*d={w*dbr:.2e} true margin (2/d)lmin at break window={lt:.3e} law .034n^-.7={.034*n**-.7:.2e}; eigvec: R(l)={Rl:.2e} R/d={Rl/dbr:.3f} |f(end)|^2={endv**2:.2f}')
# 2 kappa dips, log v steps
print('2. Krein statistics')
pc=[n for n in ns if n<44]
for inv in [200,800,3200]:
    d=1/inv; t=np.arange(int(4.0/d)+2)*d; P=Parch(t)-comb(t,ln,u); P[0]=0
    K,V,_=lev(incr(P),stop=False); lv=np.log(V)
    dips=[];steps=[]
    for n in pc:
        idx=int(np.searchsorted(ns,n)); l=ln[idx]; k=int(math.ceil(l/d))
        if k>=len(K)-1: continue
        base=np.median(K[max(k-int(.3/d),1):k+int(.3/d)])
        dips.append((n,K[k]-base,K[k+1]-base,(K[k]-base)/u[idx]))
        dv=np.diff(lv); h=int(round(.01/d)) if inv>200 else 2
        st=lv[k+h]-lv[k-h]; tr=np.median(dv[max(k-int(.25/d),0):k+int(.25/d)])*2*h
        steps.append((n,u[idx],st-tr,lv[k]-lv[k-1]))
    dips=np.array(dips);steps=np.array(steps)
    print(f'd=1/{inv} first dips k[entry]-base n=2,3,4:',np.round(dips[:3,1],3),'next',np.round(dips[:3,2],3),'ratio/w',np.round(dips[:3,3],3),' min kappa near after log2,3,4',[round(float(K[int(math.ceil(math.log(n)/d)):int(math.ceil(math.log(n)/d))+6].min()),3) for n in (2,3,4)])
    r=np.corrcoef(steps[:,1],steps[:,2])[0,1]; pf=np.polyfit(steps[:,1],steps[:,2],1)
    r0=steps[:,0]>2; r2=np.corrcoef(steps[r0,1],steps[r0,2])[0,1]
    print(f'  logv step vs w: corr={r:.3f} (without n=2 {r2:.3f}) fit {pf[0]:.3f} w {pf[1]:+.3f}; single lag step corr {np.corrcoef(steps[:,1],steps[:,3])[0,1]:.3f}; size of steps max {steps[:,2].min():.3f}')
    # control: composites not prime powers
    comp=[n for n in range(6,44) if n not in set(pc)]
    cs=[]
    for n in comp:
        l=math.log(n);k=int(math.ceil(l/d)); h=int(round(.01/d)) if inv>200 else 2
        dv=np.diff(lv); cs.append(lv[k+h]-lv[k-h]-np.median(dv[max(k-int(.25/d),0):k+int(.25/d)])*2*h)
    print(f'  control composites: mean {np.mean(cs):+.4f} sd {np.std(cs):.4f}; prime powers: mean {steps[:,2].mean():+.4f} sd {steps[:,2].std():.4f}')
# DH steps correlate?
dn,dL=dh(); dl=np.log(dn); du=dL/np.sqrt(dn)
d=1/200; t=np.arange(int(3.4/d)+2)*d; P=Pdh_arch(t)-comb(t,dl,du); P[0]=0
K,V,kk=lev(incr(P),stop=True); lv=np.log(V); dv=np.diff(lv); s=[]
for i,n in enumerate(dn):
    k=int(math.ceil(dl[i]/d))
    if k+2>=len(lv) or n<3: continue
    s.append((n,du[i],lv[k+2]-lv[k-2]-np.median(dv[max(k-50,0):k+50])*4))
s=np.array(s); print('DH d=1/200 n<=30: corr step vs |u|',np.corrcoef(np.abs(s[:,1]),s[:,2])[0,1],'vs signed u',np.corrcoef(s[:,1],s[:,2])[0,1],'N',len(s))
# 3 long
print('3. long grid')
d=.025; t=np.arange(int(16/d)+2)*d; P=Parch(t)-comb(t,ln,u); P[0]=0; c=incr(P)
K,V,k=lev(c,stop=False); print(' breaks',np.any(V<=0),'min 1-k^2',np.min(1-K[1:]**2),'at k',1+np.argmin(1-K[1:]**2),'excl k=1..5',np.min(1-K[6:]**2),'v/v0',[round(float(V[int(a/d)]/V[0]),4) for a in (4,8,12,16)])
for a in (4,8,12.05,16):
    m=int(round(a/d)); lam=eigh(toeplitz(c[:m]),eigvals_only=True,subset_by_index=[0,0])[0]
    tau=np.sinh(a*.5)/np.tanh(.5*d/2)*2.96e-12
    print(f' a={a} lam_min={lam:.4e} (2/d)lam={2*lam/d:.3e} tau(beta0=1/2,H=3e12)={tau:.3e} ratio={lam/tau:.2f}')
# 4 N(T) bound
T=3e12; main=(math.log(T/(2*math.pi))+1)/(2*math.pi*T)
R=0.1038*math.log(T)+0.2573*math.log(math.log(T))+9.3675
print('4. sum_{gamma>T}gamma^-2 main',main,'both signs',2*main,'error both signs bound',2*(R/T**2+R/T**2*1.1),' t=16 bound',2*main*(1+math.cosh(8)),'relative to mean .0462',2*main*(1+math.cosh(8))/.0462,'Psi(16)=',P[int(16/d)])
print(' log(T/2pi)',math.log(T/(2*math.pi)),'log(T/4pi)',math.log(T/(4*math.pi)),'Nyq',math.pi/d,'log(Nyq/2pi)',math.log(math.pi/d/(2*math.pi)))
print(' e^27',math.exp(27))
# 5 A_a
A=np.cumsum(2*u); 
for L0 in (2*.8,2*1.6,4.0,26.89):
    pass
j=np.searchsorted(A,math.log(T/(2*math.pi))); print('5. window where 2 sum u_n = log(T0/2pi):',ln[j],'n=',ns[j],' A at window 1.6,3.2,4:',A[np.searchsorted(ln,1.6)-1],A[np.searchsorted(ln,3.2)-1],A[np.searchsorted(ln,4.)-1])
# 6 explicit formula at t=0 term
print('6. delta at zero test')
z=np.array([float(x) for x in open(ROOT+'/data/zeros_8000.txt').read().split()])
from scipy.special import digamma as dg
for s in (.3,.15,.6):
    lhs=2*np.sum(np.exp(-z**2*s**2/2))
    r=np.linspace(0,400/s/10,800001); h=np.exp(-s**2*r**2/2)
    gam=np.trapezoid(h*np.real(dg(.25+.5j*r)),r)/np.pi
    g=lambda x:np.exp(-x**2/(2*s**2))/(s*math.sqrt(2*math.pi))
    pole=2*math.exp(s**2/8); pr=2*np.sum(u[:5000]*g(ln[:5000]))
    print(f' s={s}: zeros {lhs:.8f}; note display {pole+gam-pr:.8f}; with -log pi delta {pole+gam-pr-math.log(math.pi)*g(0.):.8f}')
# 7 Single omission at 1.6
print('7. single omission, window 1.6')
for inv in [320,1280,5120]:
    d=1/inv; t=np.arange(int(1.62/d)+2)*d
    out=[]
    for om in (None,2,3,4):
        uu=u.copy()
        if om: uu[int(np.searchsorted(ns,om))]=0.
        P=Parch(t)-comb(t,ln,uu); P[0]=0; K,V,k=lev(incr(P)); out.append((om,None if k is None else round((k+1)*d,5)))
    print(' d=1/%d'%inv,out)
# 8 Galerkin arch
print('8. arch Galerkin, Legendre, sectors')
def gal(a,N):
    L=a/2; X,WX=np.polynomial.legendre.leggauss(N+8)
    # orthonormal Legendre on [-L,L]
    def basis(x):
        y=x/L; V=Lg.legvander(y,N); return V*np.sqrt((2*np.arange(N+1)+1)/(2*L))
    def dbasis(x):
        y=x/L; out=np.zeros((len(y),N+1))
        for i in range(N+1):
            cc=np.zeros(i+1);cc[i]=1; out[:,i]=Lg.legval(y,Lg.legder(cc)) if i>0 else 0
        return out*np.sqrt((2*np.arange(N+1)+1)/(2*L))/L
    u_,wu=np.polynomial.legendre.leggauss(240); u_=(u_+1)/2; wu=wu/2
    tt=a*u_**4; wt=a*4*u_**3*wu; Pt=Parch(tt)
    W=np.zeros((N+1,N+1))
    for q in range(len(tt)):
        lo=-L;hi=L-tt[q]; y=(hi+lo)/2+(hi-lo)/2*X; wy=(hi-lo)/2*WX
        C=(dbasis(y+tt[q])*wy[:,None]).T @ dbasis(y)
        W-=wt[q]*Pt[q]*(C+C.T)
    Em=basis(np.array([-L]))[0];Ep=basis(np.array([L]))[0]
    # boundary: Psi(x+L) with x in (-L,L)= tt-L ; Psi(L-x) x=L-tt
    Dm=dbasis(tt-L); Dp=dbasis(L-tt)
    I1=(wt*Pt) @ Dm; I2=(wt*Pt) @ Dp   # int Psi(x+L) g'(x), int Psi(L-x) g'(x)
    W-=np.outer(Em,I1)+np.outer(I1,Em)-np.outer(Ep,I2)-np.outer(I2,Ep)
    Pa=Parch(np.array([a]))[0]
    W+=Pa*(np.outer(Em,Ep)+np.outer(Ep,Em))
    return W
for N in (20,40,80):
    def root(par):
        lo,hi=(.5,1.0)
        for _ in range(26):
            a=(lo+hi)/2; W=gal(a,N); idx=np.arange(par,N+1,2); lam=np.linalg.eigvalsh(W[np.ix_(idx,idx)])[0]
            if lam>0: lo=a
            else: hi=a
        return lo
    print(f' N={N}: even root {root(0):.6f} odd root {root(1):.6f}; lam even,odd at log2:',[float(np.linalg.eigvalsh(gal(math.log(2),N)[np.ix_(np.arange(p,N+1,2),np.arange(p,N+1,2))])[0]) for p in (0,1)])
