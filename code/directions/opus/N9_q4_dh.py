# Q4: the Davenport and Heilbronn control: same process, coefficients of -f'/f, no pole, conductor 5.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
ctx.prec=128
ka=arb(10-2*arb(5).sqrt()).sqrt(); ka=(ka-2)/(arb(5).sqrt()-1)
ar=[None,arb(1),ka,-ka,arb(-1)]
hur=[None]+[acb(arb(r)/5) for r in range(1,5)]
l5=arb(5).log(); lq=(arb(5)/arb.pi()).log()
def LamDH(s):
    fs=sum((acb(ar[r])*s.zeta(hur[r]) for r in range(1,5)),acb(0))*(-s*acb(l5)).exp()
    return (s/2*acb(lq)).exp()*((s+1)/2).gamma()*fs
ts=np.arange(0.02,200.0001,0.02); vals=[]
for t in ts:
    v=LamDH(acb(0.5,float(t))); vals.append(float(v.real.mid()) if True else 0)
vals=np.array(vals); sc=np.nonzero(vals[:-1]*vals[1:]<0)[0]
print('sign changes on the line in (0,200]:',len(sc),' scan %.1fs'%(time.time()-t0))
Gon=[]
for i in sc:
    lo,hi=float(ts[i]),float(ts[i+1]); flo=vals[i]
    for _ in range(40):
        mid=0.5*(lo+hi); fm=float(LamDH(acb(0.5,mid)).real.mid())
        if flo*fm<=0: hi=mid
        else: lo=mid; flo=fm
    Gon.append(0.5*(lo+hi))
Gon=np.array(Gon); print('first on line:',np.round(Gon[:6],5),'...',np.round(Gon[-2:],4))
off=[(0.8085171824566373855534-0.5,85.69934848537759217193),(0.6508300806097370824038-0.5,114.1633427307569809042),(0.5743560504508059907215-0.5,166.4793059131681558765),(0.7242576946268097802112-0.5,176.7024612428558250545)]
print('off line (beta,gamma):',off)
# mean sum_rho gamma_rho^-2 = (log xi_f)''(1/2) by a central difference in balls
ctx.prec=512
ka=arb(10-2*arb(5).sqrt()).sqrt(); ka=(ka-2)/(arb(5).sqrt()-1)
ar=[None,arb(1),ka,-ka,arb(-1)]; hur=[None]+[acb(arb(r)/5) for r in range(1,5)]; l5=arb(5).log(); lq=(arb(5)/arb.pi()).log()
h=arb(10)**(-30)
fv=lambda x: LamDH(acb(x)).real
lg=lambda x: abs(fv(x)).log()
mean=fl((lg(arb(1)/2+h)+lg(arb(1)/2-h)-2*lg(arb(1)/2))/(h*h)); print('sum_rho gamma_rho^-2 =',mean,' Lambda(1/2)=',fl(fv(arb(1)/2)))
ctx.prec=128
# prime side
N=int(math.exp(14.2))+2; Lam=dh_coeffs(N); print('coefficients to %d, %.1fs; Lambda_f(2,3,4,6)='%(N,time.time()-t0),Lam[2],Lam[3],Lam[4],Lam[6],' max|Lam(n)|/sqrt n up to N: %.1f'%np.max(np.abs(Lam[2:])/np.sqrt(np.arange(2,N+1))))
nn=np.arange(2,N+1); mk=Lam[2:]!=0; lnn=np.log(nn[mk]); uu=Lam[2:][mk]/np.sqrt(nn[mk])
def Psi_prime(t):
    q=lnn<=t
    return 0.5*t*(math.log(5/math.pi)+float(digamma(0.75)))+0.25*float(polygamma(1,0.75))-series_E(t,1.5)-np.sum(uu[q]*(t-lnn[q]))
def Psi_zero(t,raw=False):
    s=2*np.sum(np.cos(Gon*t)/Gon**2)
    for b,g in off:
        z=g+1j*b; s+=4*(np.cos(z*t)/z**2).real
    return mean-s
print('t | prime formula | zero sum to T=200 (completed) | difference')
for t in [0.1,0.5,1,2,3,4,5,6,8,10,12]:
    print(' %5.2f %.8f %.8f %+.2e'%(t,Psi_prime(t),Psi_zero(t),Psi_prime(t)-Psi_zero(t)))
tt=np.arange(0.1,14.0,0.01); vv=np.array([Psi_prime(t) for t in tt]); neg=np.nonzero(vv<0)[0]
print('first t with Psi_DH<0 on [0.1,14]:',tt[neg[0]] if len(neg) else None,' min %.5f at %.2f, max %.5f at %.2f'%(vv.min(),tt[int(np.argmin(vv))],vv.max(),tt[int(np.argmax(vv))]))
def plant(d,m,gam,beta):
    z=gam+1j*beta
    w=(1-np.cos(z*d))/z**2
    return 4*np.real(w*np.cos(np.arange(m)*z*d))
def top_peaks(u,d,k=4):
    n=1<<int(math.ceil(math.log2(16*len(u)+1)))
    F=np.abs(np.fft.rfft(u,n)); pk=[]
    idx=np.nonzero((F[1:-1]>F[:-2])&(F[1:-1]>F[2:]))[0]+1
    idx=idx[np.argsort(F[idx])[::-1]][:k]
    return [(round(2*math.pi*i/n/d,2),round(float(F[i]/F.max()),2)) for i in idx]
print('\nFirst non PD window on the step grid; eigenvector of the lowest eigenvalue at that window and 0.25 later')
for d in [0.025,0.02,0.0125,0.00625,0.003125]:
    M=int(14.0/d)
    c=cov_dh(d,M,Lam)
    ks,vs,f=levinson(c)
    if f is None: print(' delta=%g no break to 14; min v/v0=%.4f max|k|=%.3f'%(d,vs.min()/vs[0],np.abs(ks).max())); continue
    m1=f+1; aa=m1*d
    line=' delta=%g Nyquist %.1f: last PD %.4f, first non PD %.4f; kappa last 6: %s; v/v0 before break %.4f; max|k| before %.3f'%(d,math.pi/d,f*d,aa,np.round(ks[max(1,f-5):f+1],3),vs[f-1]/vs[0],np.abs(ks[1:f]).max())
    print(line)
    if m1<=3200:
        w,V=eigh(toeplitz(c[:m1]),subset_by_index=[0,1]); print('   at first non PD: lam=%.3e (next %.3e), peaks (gamma, rel):'%(w[0],w[1]),top_peaks(V[:,0],d))
        m2=min(M,m1+int(0.25/d))
        if m2<=3200:
            w,V=eigh(toeplitz(c[:m2]),subset_by_index=[0,1]); print('   at window %.3f: lam=%.3e, %.3e; peaks:'%(m2*d,w[0],w[1]),top_peaks(V[:,0],d))
print('time %.1f'%(time.time()-t0))
# peeling: remove the known quadruples one at a time on the grid
print('\nPeeling on delta=0.025 and delta=0.02: remove quadruples from c, break window and peak')
for d in [0.025,0.02]:
    M=int(14.0/d); c=cov_dh(d,M,Lam)
    for rem in [[],[0],[0,1],[0,3],[0,1,3],[0,1,2,3]]:
        cc=c.copy()
        for j in rem: cc-=plant(d,M,off[j][1],off[j][0])
        ks,vs,f=levinson(cc)
        if f is None:
            print(' delta=%g removed %s: no break to a=14; lam_min %.3e'%(d,[round(off[j][1],1) for j in rem],lam_min(cc,M))); continue
        w,V=eigh(toeplitz(cc[:f+1]),subset_by_index=[0,0])
        print(' delta=%g removed %s: first non PD %.3f, peaks %s'%(d,[round(off[j][1],1) for j in rem],(f+1)*d,top_peaks(V[:,0],d,3)))
# lam_min approaching the break, delta=0.025
d=0.025; M=int(14.0/d); c=cov_dh(d,M,Lam); ks,vs,f=levinson(c)
print('\nApproach on delta=0.025: a, 2 lam_min/delta, |kappa| local max over last 0.25, v/v0')
for a in np.arange(1,(f+1)*d+0.001,0.5):
    m=int(round(a/d)); print('  a=%.2f eps=%.4e max|k| window=%.3f v/v0=%.4f'%(a,2*lam_min(c,m)/d,np.abs(ks[max(1,m-10):m]).max(),vs[m-1]/vs[0]))
# ball certificate of the negative direction at delta=1/40
ctx.prec=512
k_b=arb(10-2*arb(5).sqrt()).sqrt(); k_b=(k_b-2)/(arb(5).sqrt()-1)
A5=[arb(0),arb(1),k_b,-k_b,arb(-1)]
m1=f+1; Nb=int(math.exp((m1+1)*d))+2
S=[arb(0)]*(Nb+1); Lb=[arb(0)]*(Nb+1)
for dd in range(1,Nb+1):
    Lb[dd]=A5[dd%5]*arb(dd).log()-S[dd]
    for q in range(2,Nb//dd+1):
        if q%5: S[dd*q]=S[dd*q]+Lb[dd]*A5[q%5]
delta=arb(1)/40; pi=arb.pi(); g=arb.const_euler(); l2=arb.const_log2(); G_=arb.const_catalan()
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
w,V=eigh(toeplitz(np.array([fl(x) for x in cb])),subset_by_index=[0,0])
R=rayleigh_ball(cb,V[:,0]); print('\nball Rayleigh quotient at first non PD window %.3f (m=%d, delta=1/40): %.6e +- %.1e  certified negative: %s'%(m1*d,m1,fl(R),fl(R.rad()),R<0))
ok,info=chol_cert(cb[:f],0.0); print('ball Cholesky of T_%d (window %.3f) positive: %s'%(f,f*d,ok))
print('total %.1f'%(time.time()-t0))
