#!/usr/bin/env python3
"""N8 check 3: (A) what the rows of the note's table really are; (B) levels per Gram interval; (C) reference phase with a finite Euler product,
interlacing and an explicit rank one matrix; (D) onset of the DH kernel failure, decomposed on the zeros; parts of Psi; (E) traces at small b.
Run: timeout 600 python3 -I N8_check3.py"""
import json, math, time
import numpy as np, mpmath as mp
ROOT='./'; HERE=ROOT+'misc/research/directions/opus/'
Z=np.array([float(x) for x in json.load(open(HERE+'N8_zeros.json'))])
PI=math.pi; T0=time.time()
Bn=[1/12,-1/360,1/1260,-1/1680,1/1188,-691/360360,1/156]
def lg(z,shift=12):
    z=np.asarray(z,dtype=complex); s=np.zeros_like(z); w=z.copy()
    for k in range(shift):
        s=s-np.log(w); w=w+1
    iw=1/w; iw2=iw*iw; ser=0; p=iw
    for b in Bn:
        ser=ser+b*p; p=p*iw2
    return s+(w-0.5)*np.log(w)-w+0.5*math.log(2*PI)+ser
def digam(z,shift=12):
    z=np.asarray(z,dtype=complex); s=0; w=z.copy()
    for k in range(shift):
        s=s-1/w; w=w+1
    x=1/(w*w)
    return s+np.log(w)-0.5/w-x*(1/12-x*(1/120-x*(1/252-x*(1/240-x*(1/132)))))
def theta(t):
    t=np.asarray(t,dtype=float); return lg(0.25+0.5j*t).imag-0.5*t*math.log(PI)
def F01(b,z,K=70):
    term=np.ones_like(b,dtype=complex); s=term.copy()
    for k in range(K):
        term=term*z/((k+1)*(b+k)); s=s+term
    return s
def Ksc(a,t,r=2*PI):
    t=np.asarray(t,dtype=float); mu=a+0.5j*t; z=(r/2)**2
    A=np.exp(lg(mu)+0.25*PI*t-mu*math.log(r/2))*F01(1-mu,z)
    Bb=np.exp(lg(-mu)+0.25*PI*t+mu*math.log(r/2))*F01(1+mu,z)
    return 0.5*(A+Bb)
def levels(fun,lo=1.0,hi=705.0,h=0.01):
    t=np.arange(lo,hi,h); v=fun(t); idx=np.where(v[:-1]*v[1:]<0)[0]
    a=t[idx].copy(); b=t[idx+1].copy(); fa=v[idx].copy()
    for _ in range(45):
        m=0.5*(a+b); fm=fun(m); left=fa*fm<=0
        b=np.where(left,m,b); a=np.where(left,a,m); fa=np.where(left,fa,fm)
    return 0.5*(a+b)
base=levels(lambda t:Ksc(0.,t).real); two=levels(lambda t:(2*PI**2*Ksc(2.25,t)-3*PI*Ksc(1.25,t)).real)
tg=np.arange(6.5,725,0.001); thg=theta(tg)
def lev_c(c):
    k=np.arange(math.ceil(thg[0]/PI-c),math.floor(thg[-1]/PI-c)+1); return np.interp((k+c)*PI,thg,tg)
print('(A) constant phase surrogates')
for c,nm in ((0.5,'c=1/2, note row two term (397,207,20.65)'),(0.625,'c=5/8, note row Liouville (397,188,26.05)')):
    L0=lev_c(c)
    for lo,hi in ((14.,680.),(17.,680.),(20.,680.),(14.,679.74)):
        L=L0[(L0>=lo)&(L0<=hi)]; cnt=np.searchsorted(Z,L[1:])-np.searchsorted(Z,L[:-1]); bad=np.where(cnt!=1)[0]
        print(' ',nm,'levels in [%g,%g]: intervals %d bad %d first failure %.2f first level %.3f'%(lo,hi,len(L)-1,len(bad),L[bad[0]],L[0]))
print(' true first levels: two term',two[:5],' base',base[:5],' surrogates c=1/2',lev_c(.5)[:6],' c=5/8',lev_c(.625)[:6])
print(' true levels minus surrogate at t about 30, 100, 300, 680:')
for nm,L,c in (('two',two,.5),('base',base,.625)):
    fr=(theta(L)/PI-c+0.5)%1-0.5
    for t0 in (30,100,300,680):
        i=np.argmin(abs(L-t0)); print('   %s t=%.2f theta/pi-c mod 1 = %.4f (in t units %.4f)'%(nm,L[i],fr[i],fr[i]*PI/(0.5*math.log(L[i]/(2*PI)))))
print('(B) levels per Gram interval')
g=lev_c(0.0); g=g[g<=700]
for nm,L in (('base',base),('two',two)):
    cnt=np.searchsorted(L,g[1:])-np.searchsorted(L,g[:-1])
    print(' ',nm,'from g_{-1}=%.3f: counts per interval'%g[0],np.bincount(cnt),'first bad',g[:-1][cnt!=1][:3],'levels below g_-1',(L<g[0]).sum())
print('(C) reference phase Theta_P=theta+Im log zeta_P(1/2+it)')
t=np.arange(30,680.00005,0.0005); th=theta(t); Nz=np.searchsorted(Z,t,side='right').astype(float)
def corr(primes):
    s=np.zeros_like(t)
    for p in primes:
        lp=math.log(p); kmax=max(1,int(2*math.log(1e4)/lp))
        for k in range(1,kmax+1): s-=np.sin(k*t*lp)/(k*p**(k/2))
    return s
for pr in ([],[2],[2,3],[2,3,5],[2,3,5,7,11,13]):
    Th=th+corr(pr); sp=[]
    for c in np.arange(0,1,1/200):
        env=np.maximum.accumulate(np.floor(Th/PI-c)); D=Nz-env; sp.append(int(D.max()-D.min()))
    sp=np.array(sp); ok=np.arange(0,1,1/200)[sp==1]
    env=np.maximum.accumulate(np.floor(Th/PI)); jmp=np.where(np.diff(env)>0)[0]+1; L=t[jmp]
    cnt=np.searchsorted(Z,L[1:])-np.searchsorted(Z,L[:-1])
    print('  primes',pr,'span counts over 200 phases',dict(zip(*[x.tolist() for x in np.unique(sp,return_counts=True)])),'phases with span 1: %d (%s..%s)'%(len(ok),ok.min() if len(ok) else '-',ok.max() if len(ok) else '-'),'| c=0 levels: intervals %d, without exactly one zero %d'%(len(L)-1,(cnt!=1).sum()))
    if pr==[2,3]:
        idx=np.searchsorted(Z,L[:-1]); ok1=(cnt==1).all()
        mu=L[:-1]**2; nu=Z[idx]**2; m=len(mu)
        lv=np.zeros(m); sg=np.ones(m)
        for j in range(m):
            num=nu-mu[j]; den=mu-mu[j]; den[j]=1.0
            lv[j]=np.log(np.abs(num)).sum()-np.log(np.abs(den)).sum(); sg[j]=np.sign(num).prod()*np.sign(den).prod()
        v=np.sqrt(np.exp(lv))
        ev=np.linalg.eigvalsh(np.diag(mu)+np.outer(v,v))
        print('   explicit rank one matrix diag(levels^2)+vv^T, size',m,': all signs positive',(sg>0).all(),'; max |sqrt(eig)-gamma| = %.2e'%np.abs(np.sqrt(ev)-np.sqrt(nu)).max(),'; first gamma used %.3f last %.3f'%(math.sqrt(nu[0]),math.sqrt(nu[-1])))
print('t',round(time.time()-T0,1))
print('(D) kernels')
kap=(math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
N=int(math.exp(14.0))+2; an=np.arange(N+1,dtype=float)
a5=np.zeros(N+1); a5[1::5]=1.0; a5[2::5]=kap; a5[3::5]=-kap; a5[4::5]=-1.0
r=a5*np.log(np.maximum(an,1.0)); c=np.zeros(N+1)
for dd in range(2,N//2+1):
    cd=r[dd]; c[dd]=cd
    if cd!=0.0:
        m=N//dd; r[2*dd::dd]-=cd*a5[2:m+1]
c[N//2+1:]=r[N//2+1:]
nz=np.nonzero(c)[0]; L_dh=np.log(nz.astype(float)); W_dh=c[nz]/np.sqrt(nz.astype(float)); o=np.argsort(L_dh); L_dh=L_dh[o]; W_dh=W_dh[o]; A_dh=np.cumsum(W_dh); B_dh=np.cumsum(W_dh*L_dh)
def sumexp(tt,a0,Kmax=400000):
    tt=np.asarray(tt,dtype=float); res=np.zeros_like(tt)
    for k in range(Kmax):
        a=2*k+a0; j=np.searchsorted(tt,45/a)
        if j==0: break
        res[:j]+=np.exp(-a*tt[:j])/a**2
    return res
psi34=float(mp.digamma(0.75)); tri34=float(mp.psi(1,0.75)); psi14=float(mp.digamma(0.25)); tri14=float(mp.psi(1,0.25))
def Psi_dh(tt):
    tt=np.asarray(tt,dtype=float); idx=np.searchsorted(L_dh,tt,side='right')
    pr=np.where(idx>0,tt*A_dh[np.maximum(idx-1,0)]-B_dh[np.maximum(idx-1,0)],0.0)
    return 0.5*tt*(psi34+math.log(5/PI))+0.25*tri34-sumexp(tt,1.5)-pr
def sieve(Nn):
    s=np.ones(Nn+1,bool); s[:2]=False
    for i in range(2,int(Nn**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
pr=sieve(N); L1=np.log(pr.astype(float)); W1=L1/np.sqrt(pr.astype(float)); eL=[];eW=[]
for p in pr[pr<=int(N**.5)]:
    p=int(p); q=p*p
    while q<=N:
        eL.append(math.log(q)); eW.append(math.log(p)/math.sqrt(q)); q*=p
LL=np.concatenate([L1,np.array(eL)]); WW=np.concatenate([W1,np.array(eW)]); o=np.argsort(LL); LL=LL[o]; WW=WW[o]; Az=np.cumsum(WW); Bz=np.cumsum(WW*LL)
def part_prime(tt):
    tt=np.asarray(tt,dtype=float); idx=np.searchsorted(LL,tt,side='right')
    return -np.where(idx>0,tt*Az[np.maximum(idx-1,0)]-Bz[np.maximum(idx-1,0)],0.0)
def part_pole(tt): return 8*(np.cosh(np.asarray(tt)/2)-1)
def part_arch(tt): return 0.5*np.asarray(tt)*(psi14-math.log(PI))+0.25*tri14-sumexp(tt,0.5)
def Psi_z(tt): return part_prime(tt)+part_pole(tt)+part_arch(tt)
def kern(Pfun,ell,n):
    u=ell*np.arange(1,n+1)/n; P0=np.concatenate([[0.0],Pfun(u)]); i=np.arange(1,n+1); Ii,Jj=np.meshgrid(i,i,indexing='ij')
    return P0[Ii]+P0[Jj]-P0[np.abs(Ii-Jj)],u
def rayleigh(G,ell,n):
    """matrix of Q on the cell values F_j of a step function, normalised by the L2 norm: F=sum c_i 1_[0,u_i], c_i=F_i-F_{i+1}"""
    h=ell/n; D=np.eye(n)-np.eye(n,k=1)
    return D.T @ G@D/h
print(' ell | DH least eig of G, least Q/||F||^2 on step functions | zeta same')
for ell in (3.0,3.2,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,5.0,6.0):
    n=1200; G,_=kern(Psi_dh,ell,n); Gz,_=kern(Psi_z,ell,n)
    print(' %.2f | %.3e %.3e | %.3e %.3e'%(ell,np.linalg.eigvalsh(G)[0],np.linalg.eigvalsh(rayleigh(G,ell,n))[0],np.linalg.eigvalsh(Gz)[0],np.linalg.eigvalsh(rayleigh(Gz,ell,n))[0]))
print(' parts of Psi at ell=3, n=600: least eigenvalue of the kernel of prime part, pole part, archimedean part, sum')
for ell in (1.0,3.0,6.0):
    print('  ell=%g'%ell,['%.3e'%np.linalg.eigvalsh(kern(f,ell,600)[0])[0] for f in (part_prime,part_pole,part_arch,Psi_z)],' largest',['%.3e'%np.linalg.eigvalsh(kern(f,ell,600)[0])[-1] for f in (part_prime,part_pole,part_arch)])
# decomposition of the DH negative direction on zeros
tgr=np.arange(0.5,200,0.02)
BEM=[1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510]
def hurw(s,a,M=500):
    s=np.asarray(s,dtype=complex); w=M+a; tot=np.zeros_like(s); ln=np.log(np.arange(M)+a)
    for i in range(0,len(s),400):
        ss=s[i:i+400]; tot[i:i+400]=np.exp(-np.outer(ss,ln)).sum(axis=1)
    tot=tot+w**(1-s)/(s-1)+0.5*w**(-s)+BEM[0]/2*s*w**(-s-1); prod=s.copy(); pw=w**(-s-1); f=2.0
    for j in range(2,len(BEM)+1):
        prod=prod*(s+2*j-3)*(s+2*j-2); pw=pw/w**2; f=f*(2*j-1)*(2*j); tot=tot+BEM[j-1]/f*prod*pw
    return tot
def f_dh(s):
    s=np.asarray(s,dtype=complex); return 5.0**(-s)*(hurw(s,.2)+kap*hurw(s,.4)-kap*hurw(s,.6)-hurw(s,.8))
def thdh(tt):
    tt=np.asarray(tt,dtype=float); return lg(0.75+0.5j*tt).imag+0.5*tt*math.log(5/PI)
v=(np.exp(1j*thdh(tgr))*f_dh(0.5+1j*tgr)).real; idx=np.where(v[:-1]*v[1:]<0)[0]; a=tgr[idx].copy(); b=tgr[idx+1].copy(); fa=v[idx].copy()
for _ in range(40):
    m=0.5*(a+b); fm=(np.exp(1j*thdh(m))*f_dh(0.5+1j*m)).real; left=fa*fm<=0
    b=np.where(left,m,b); a=np.where(left,a,m); fa=np.where(left,fa,fm)
on=0.5*(a+b)
d=json.load(open(ROOT+'data/cf_dhzeros.json'))['T2']['zeros']; off=np.array([[float(z['center_re']),float(z['center_im'])] for z in d if 'center_re' in z]); off=off[off[:,1]<200]
for ell in (4.0,6.0):
    n=1200; G,u=kern(Psi_dh,ell,n); ev,V=np.linalg.eigh(G); cc=V[:,0]
    def E(gm): return (cc*(np.exp(1j*gm*u)-1)).sum()/(1j*gm)
    Qon=sum(2*abs(E(g))**2 for g in on)
    Qoff=[4*(E(g-1j*(x-0.5))*np.conj(E(g+1j*(x-0.5)))).real for x,g in off]
    print(' DH ell=%g: c^T G c = %.5e; on line zeros to 200 give +%.5e; off line pairs at'%(ell,ev[0],Qon),[round(float(g),2) for g in off[:,1]],'give',['%.4e'%q for q in Qoff],'; remainder (zeros above 200) %.3e'%(ev[0]-Qon-sum(Qoff)))
print('(E) traces at small b: b, W, A/2-W_wall, A/2-W_twoterm')
def Arch(b):
    R=math.sqrt(45/b); hh=0.005; rr=np.arange(0,R,hh); f=np.exp(-b*rr*rr)*digam(0.25+0.5j*rr).real
    return 2*math.exp(b/4)-math.log(PI)/math.sqrt(4*PI*b)+(hh*(f.sum()-0.5*f[0]))/PI
for b in (8e-5,1.5e-4,3e-4,6e-4,1e-3,2e-3):
    Aa=Arch(b); print('  %.5f %.6f %.6f %.6f'%(b,np.exp(-b*Z*Z).sum(),0.5*Aa-np.exp(-b*base*base).sum(),0.5*Aa-np.exp(-b*two*two).sum()))
print('t total',round(time.time()-T0,1))
