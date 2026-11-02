# Q5: (A) forecast fraction and reflection coefficients near a spike against delta, with the spike exactly on a lag;
# (B) the truncated break across long prime power gaps; (C) the same truncation test for Davenport and Heilbronn; (D) DH break by delta.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
ns,Ls=prime_powers(40000); ln=np.log(ns); u=Ls/np.sqrt(ns)
print('A. true process, grid tuned so log n0=k d. f=forecast fraction, kappa at lags k-2..k+1, s/v, truncated break delay and its size in lags')
for d0 in [0.025,0.0125,0.00625,0.003125,0.0015625,0.00078125]:
    print(' d0=%g'%d0)
    for n0 in [2,3,5,7,11,13,17,23,31,47]:
        k=int(round(math.log(n0)/d0)); d=math.log(n0)/k; i=int(np.searchsorted(ns,n0))
        A=cov_arch_zeta(d,k+int(0.12/d)+6); T=A.copy(); add_primes(T,d,ln[:i],u[:i]); add_primes(A,d,ln,u)
        ks,vs,fl_=levinson(A[:k+3]); kt,vt,ft=levinson(T)
        s=0.5*u[i]*d; v=vs[k-1]; pred=A[k]-ks[k]*v; cno=A[k]+s
        dl=((ft+1)*d-math.log(n0)) if ft is not None else float('nan')
        print('  n=%3d k=%5d f=%6.3f  kappa[k-2..k+1]= %+.3f %+.3f %+.3f %+.3f  s/v=%.3f  delay=%.4f (%.1f lags)'%(n0,k,(cno-pred)/s,ks[k-2],ks[k-1],ks[k],ks[k+1],s/v,dl,dl/d))
print('time %.1f'%(time.time()-t0))
print('B. keep prime powers up to N; next prime power n0. First non PD window a, e^a, against the gap')
for d in [0.003125,0.00078125]:
    for N,n0 in [(13,16),(19,23),(89,97),(113,121),(139,149),(199,211),(293,307),(523,529),(887,907)]:
        i0=int(np.searchsorted(ns,N,side='right')); assert ns[i0]==n0
        c=cov_arch_zeta(d,int((ln[i0]+0.12)/d)+2); add_primes(c,d,ln[:i0],u[:i0]); ks,vs,f=levinson(c)
        a=(f+1)*d; print('  d=%.6f N=%d n0=%d log gap=%.4f first non PD %.4f e^a=%.2f delay=%.4f'%(d,N,n0,ln[i0]-math.log(N),a,math.exp(a),a-ln[i0]))
print('time %.1f'%(time.time()-t0))
print('C. DH: keep Lambda_f(n) for n<=N, omit the rest. Same test. u0=Lambda_f(n0)/sqrt(n0)')
Lam=dh_coeffs(200); nz=[n for n in range(2,60) if abs(Lam[n])>1e-12]
print('  nonzero n below 60:',nz)
for N in [1]+nz[:26]:
    n0=[n for n in nz if n>N][0]; lo=math.log(n0); row=[]
    for d in [0.00625,0.0015625,0.00078125]:
        L2=Lam.copy(); L2[N+1:]=0.0; c=cov_dh(d,int((lo+0.6)/d)+2,L2); ks,vs,f=levinson(c)
        row.append('%.5f: %s'%(d,('%+.4f'%((f+1)*d-lo)) if f is not None else 'none'))
    print('  N=%d n0=%d log n0=%.4f u0=%+.4f | delay %s'%(N,n0,lo,Lam[n0]/math.sqrt(n0),' ; '.join(row)))
print('time %.1f'%(time.time()-t0))
print('D. DH first non PD window against delta')
Lam=dh_coeffs(int(math.exp(6.2))+2)
for d in [0.025,0.0125,0.00625,0.003125,0.0015625,0.00078125,0.000390625]:
    c=cov_dh(d,int(6.0/d),Lam); ks,vs,f=levinson(c); print('  d=%.7f first non PD %.5f last PD %.5f X=%.2f'%(d,(f+1)*d,f*d,math.exp((f+1)*d)))
print('time %.1f'%(time.time()-t0))
print('E. ball certificate of the DH negative direction at delta=1/320')
ctx.prec=256
d=1/320; Lamd=dh_coeffs(100); c=cov_dh(d,int(4.0/d),Lamd); ks,vs,f=levinson(c); m1=f+1
w,V=eigh(toeplitz(c[:m1]),subset_by_index=[0,0])
k_b=arb(10-2*arb(5).sqrt()).sqrt(); k_b=(k_b-2)/(arb(5).sqrt()-1)
A5=[arb(0),arb(1),k_b,-k_b,arb(-1)]
Nb=int(math.exp((m1+1)*d))+2
S=[arb(0)]*(Nb+1); Lb=[arb(0)]*(Nb+1)
for dd in range(1,Nb+1):
    Lb[dd]=A5[dd%5]*arb(dd).log()-S[dd]
    for q in range(2,Nb//dd+1):
        if q%5: S[dd*q]=S[dd*q]+Lb[dd]*A5[q%5]
delta=arb(1)/320; pi=arb.pi(); g=arb.const_euler(); l2=arb.const_log2(); G_=arb.const_catalan()
psi34=-g+pi/2-3*l2; tri34=pi*pi-8*G_
P=[arb(0)]
for kk in range(1,m1+2):
    t=delta*kk; ps=arb(0)
    for n in range(2,Nb+1):
        dn=t-arb(n).log()
        if dn>0: ps+=Lb[n]/arb(n).sqrt()*dn
        elif dn<0: break
    P.append((t/2)*((arb(5)/pi).log()+psi34)+tri34/4-E_ball(t,arb(3)/2)-ps)
cb=cov_from_Psi_ball(P)[:m1]
R=rayleigh_ball(cb,V[:,0]); print('  m=%d window %.5f X=%.3f: Rayleigh %.6e +- %.1e, certified negative: %s; (double lam_min %.6e)'%(m1,m1*d,math.exp(m1*d),fl(R),fl(R.rad()),R<0,w[0]))
print('total %.1f'%(time.time()-t0))
