# Q3: how the grid margin depends on delta and a; the aliasing floor; the zero side implication.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
ns,Ls=prime_powers(250000); ln=np.log(ns); u=Ls/np.sqrt(ns)
G=load_zeros(8000)
# ---- A. convergence to the ladder values: exact window a, delta=a/m
print('A. 2 lam_min/delta on the window a=log X, by parity, against opus/D1 (upper bounds)')
D1={2:(1.329433e-3,7.372781e-2),3:(5.535508e-8,1.528437e-5)}
for X,mlist in [(2,[28,56,112,224,448,896,1792,3584]),(3,[44,88,176,352,704,1408,2816,5632,11264])]:
    a=math.log(X); rec={'e':[],'o':[]}
    for m in mlist:
        d=a/m
        c=cov_zeta(d,m+1,ln[ns<X],u[ns<X])
        le=lam_min(c,m,parity='e'); lo=lam_min(c,m,parity='o')
        rec['e'].append((d,2*le/d)); rec['o'].append((d,2*lo/d))
        print(' X=%d m=%d delta=%.3e even %.6e odd %.6e'%(X,m,d,2*le/d,2*lo/d))
    for s,idx in [('e',0),('o',1)]:
        r=rec[s][-3:]
        Am=np.array([[1,x[0]**2*math.log(1/x[0]),x[0]**2] for x in r]); y=np.array([x[1] for x in r])
        sol=np.linalg.solve(Am,y)
        print('  X=%d sector %s: Richardson eps0=%.6e, coefficient of delta^2 log(1/delta)=%.4f ; D1 upper %.6e'%(X,s,sol[0],sol[1],D1[X][idx]))
print('time %.1f'%(time.time()-t0))
# ---- B. thresholds of the archimedean form
print('\nB. archimedean form, first non PD window by Levinson; D1: 2L = 0.74386 (odd), 0.82022 (even)')
for d in [0.025,0.005,0.001,0.0002]:
    M=int(1.0/d)+2
    c=cov_arch_zeta(d,M); ks,vs,f=levinson(c)
    print(' delta=%g: last PD window %.5f first non PD %.5f; kappa before break %s'%(d,f*d,(f+1)*d,np.round(ks[max(1,f-4):f+1],3)))
print('time %.1f'%(time.time()-t0))
# ---- C. floor table
A=lambda d: math.log(1/d)/12+0.9375482543/(2*math.pi**2)
def dir_model(a,d,J=60):
    nu=np.arange(1,J+1)*np.pi/a; sg=(-1.0)**np.arange(1,J+1)
    Phi=np.sqrt(2/a)*nu[:,None]*(1-sg[:,None]*np.exp(1j*G[None,:]*a))/(nu[:,None]**2-G[None,:]**2)
    W=2*np.real(Phi @ Phi.conj().T)
    H=W+d*d*A(d)*np.diag(nu**2)
    return eigvalsh(H,subset_by_index=[0,0])[0]
print('\nC. eps_delta(a)=2 lam_min(T_m)/delta [model: least eigenvalue of W + delta^2 A_delta |phi\'|^2 on Dirichlet modes], peak frequency, ratio to previous delta')
for a in [2,3,4,6,8,10,12]:
    prev=None; line=[]
    for d in [0.025,0.0125,0.00625,0.003125]:
        m=int(round(a/d)); c=cov_zeta(d,m+1,ln,u)
        lam=lam_min(c,m); eps=2*lam/d; mod=dir_model(a,d)
        line.append('d=%.5f: %.4e [%.4e] r=%s'%(d,eps,mod,('%.2f'%(prev/eps)) if prev else '-'))
        prev=eps
    print(' a=%g: '%a+' | '.join(line))
print('time %.1f'%(time.time()-t0))
# ---- D. kappa and v depend on delta
print('\nD. Levinson statistics at a=12 by delta: kappa_1, max|kappa|, max over lags >= 0.25, min v/v0, c(0), v_last')
for d in [0.025,0.0125,0.00625,0.003125,0.0015625]:
    m=int(round(12/d)); c=cov_zeta(d,m+1,ln,u); ks,vs,f=levinson(c)
    k0=int(0.25/d); print(' delta=%.7f m=%d fail=%s kappa1=%.4f max|k|=%.4f max|k|(lag>=.25)=%.4f at lag %.3f min v/v0=%.4f v0=%.5e v_last=%.5e 2v_last/delta=%.4f'%(d,m,f,ks[1],np.abs(ks).max(),np.abs(ks[k0:]).max(),(k0+int(np.argmax(np.abs(ks[k0:]))))*d,vs.min()/vs[0],vs[0],vs[-1],2*vs[-1]/d))
print('time %.1f'%(time.time()-t0))
# ---- E. zero side: T_m >= T^(K)_m - tau I, tau = sinh(a beta0) coth(beta0 delta/2) S(H) from the zeros above H
print('\nE. zero side implication at delta=0.025, K=8000 zeros; H=3e12 (Platt, Trudgian); beta0=1/2 unless stated')
d=0.025; cK=cov_from_zeros(d,482,G); cP=cov_zeta(d,483,ln,u)
def tau(a,b0,H): return math.sinh(a*b0)/math.tanh(b0*d/2)*S_tail(H)
def Hstar(val,a,b0):
    lo,hi=10.0,1e30
    for _ in range(300):
        mid=math.sqrt(lo*hi)
        if tau(a,b0,mid)>val: lo=mid
        else: hi=mid
    return hi
print(' a | m | lam_min(T) prime side | lam_min(T^(8000)) | tau(H=3e12,1/2) | implied? | H* for beta0=1/2, 0.1, 0.01, on the line limit')
for a in [1,2,3,4,6,8,10,11,12,12.05]:
    m=int(round(a/d)); lK=lam_min(cK,m); lP=lam_min(cP,m); t=tau(a,0.5,3e12)
    print(' %g | %d | %.4e | %.4e | %.3e | %s | %.2e %.2e %.2e'%(a,m,lP,lK,t,'yes' if lK>t else 'NO',Hstar(lK,a,0.5),Hstar(lK,a,0.1),Hstar(lK,a,0.01)))
print('total %.1f'%(time.time()-t0))
