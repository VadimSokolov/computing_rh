#!/usr/bin/env python3
"""N8 control: the Davenport and Heilbronn function against every claim of the note, and the variogram kernel on a window
for zeta and for DH, from the coefficient side, no zero as input.  Off line zeros are read from data/cf_dhzeros.json.
Run: timeout 600 python3 -I N8_dh.py"""
import json, math, time
import numpy as np, mpmath as mp
ROOT='./'
PI=math.pi; T0=time.time()
kap=(math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
Bn=[1/12,-1/360,1/1260,-1/1680,1/1188,-691/360360,1/156]
def lg(z,shift=12):
    z=np.asarray(z,dtype=complex); s=np.zeros_like(z); w=z.copy()
    for k in range(shift):
        s=s-np.log(w); w=w+1
    iw=1/w; iw2=iw*iw; ser=0; p=iw
    for b in Bn:
        ser=ser+b*p; p=p*iw2
    return s+(w-0.5)*np.log(w)-w+0.5*math.log(2*PI)+ser
# Bernoulli coefficients for Euler Maclaurin
BEM=[1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510]
def hurw(s,a,M=500):
    """Hurwitz zeta for an array s by Euler Maclaurin."""
    s=np.asarray(s,dtype=complex); w=M+a
    tot=np.zeros_like(s)
    n=np.arange(M)+a; ln=np.log(n)
    # block to limit memory
    for i in range(0,len(s),400):
        ss=s[i:i+400]; tot[i:i+400]=np.exp(-np.outer(ss,ln)).sum(axis=1)
    tot=tot+w**(1-s)/(s-1)+0.5*w**(-s)
    term=s*w**(-s-1); fct=2.0
    # B_2/2! s w^{-s-1}, B_4/4! s(s+1)(s+2) w^{-s-3},...
    tot=tot+BEM[0]/2*term
    prod=s.copy(); pw=w**(-s-1); f=2.0
    for j in range(2,len(BEM)+1):
        prod=prod*(s+2*j-3)*(s+2*j-2); pw=pw/w**2; f=f*(2*j-1)*(2*j)
        tot=tot+BEM[j-1]/f*prod*pw
    return tot
def f_dh(s):
    s=np.asarray(s,dtype=complex)
    return 5.0**(-s)*(hurw(s,.2)+kap*hurw(s,.4)-kap*hurw(s,.6)-hurw(s,.8))
def thdh(t):
    t=np.asarray(t,dtype=float)
    return lg(0.75+0.5j*t).imag+0.5*t*math.log(5/PI)
mp.mp.dps=25
def fmp(s):
    return 5**(-s)*(mp.zeta(s,mp.mpf(1)/5)+kap*mp.zeta(s,mp.mpf(2)/5)-kap*mp.zeta(s,mp.mpf(3)/5)-mp.zeta(s,mp.mpf(4)/5))
def Lamp(s):
    return (5/mp.pi)**(s/2)*mp.gamma((1+s)/2)*fmp(s)
s0=mp.mpc(0.3,20); print('functional equation check',mp.nstr(Lamp(s0),10),mp.nstr(Lamp(1-s0),10),' f numpy vs mpmath',complex(f_dh(np.array([0.5+150j]))[0]),mp.nstr(fmp(mp.mpc(0.5,150)),12),' kappa',kap)
tg=np.arange(0.5,200,0.02)
Zt=(np.exp(1j*thdh(tg))*f_dh(0.5+1j*tg))
print('max |Im Z| on the line',np.abs(Zt.imag).max())
v=Zt.real; idx=np.where(v[:-1]*v[1:]<0)[0]; a=tg[idx].copy(); b=tg[idx+1].copy(); fa=v[idx].copy()
for _ in range(40):
    m=0.5*(a+b); fm=(np.exp(1j*thdh(m))*f_dh(0.5+1j*m)).real; left=fa*fm<=0
    b=np.where(left,m,b); a=np.where(left,a,m); fa=np.where(left,fa,fm)
on=0.5*(a+b); print('on line zeros below 200:',len(on),'(certified log: 122), first',on[:3])
d=json.load(open(ROOT+'data/cf_dhzeros.json'))['T2']['zeros']
off=np.array([[float(z['center_re']),float(z['center_im'])] for z in d if 'center_re' in z])
off=off[off[:,1]<200]; print('off line zeros read',off.tolist(),' |f| there',[mp.nstr(abs(fmp(mp.mpc(x,y))),3) for x,y in off])
# on the note's Psi for DH: gamma_rho=gamma-i delta
del_=off[:,0]-0.5; gam=off[:,1]
def Psi_zero(t,offline=True,T=200.):
    t=np.asarray(t,dtype=float); s=np.zeros_like(t)
    for g in on: s+=2*(1-np.cos(g*t))/g**2
    if offline:
        for g,dl in zip(gam,del_):
            gr=g-1j*dl; s+=4*((1-np.cos(gr*t))/gr**2).real
    return s+(math.log(5*T/(2*PI))+1)/(PI*T)
# coefficient side
N=int(math.exp(14.0))+2
an=np.arange(N+1,dtype=float)
a5=np.zeros(N+1); a5[1::5]=1.0; a5[2::5]=kap; a5[3::5]=-kap; a5[4::5]=-1.0
r=a5*np.log(np.maximum(an,1.0)); c=np.zeros(N+1)
for dd in range(2,N//2+1):
    cd=r[dd]; c[dd]=cd
    if cd!=0.0:
        m=N//dd
        r[2*dd::dd]-=cd*a5[2:m+1]
c[N//2+1:]=r[N//2+1:]
nz=np.nonzero(c)[0]; L_dh=np.log(nz.astype(float)); W_dh=c[nz]/np.sqrt(nz.astype(float)); o=np.argsort(L_dh); L_dh=L_dh[o]; W_dh=W_dh[o]
A_dh=np.cumsum(W_dh); B_dh=np.cumsum(W_dh*L_dh)
print('coefficient c(2),c(3),c(4),c(6):',c[2],c[3],c[4],c[6],' max|c(n)|/log n near top',np.abs(c[N-5000:]).max()/math.log(N),' t',round(time.time()-T0,1))
def sumexp(t,a0,Kmax=400000):
    t=np.asarray(t,dtype=float); out=np.zeros_like(t); o=np.argsort(t); ts=t[o]; res=np.zeros_like(ts)
    for k in range(Kmax):
        a=2*k+a0; j=np.searchsorted(ts,45/a)
        if j==0: break
        res[:j]+=np.exp(-a*ts[:j])/a**2
    out[o]=res
    return out
psi34=float(mp.digamma(0.75)); tri34=float(mp.psi(1,0.75)); psi14=float(mp.digamma(0.25)); tri14=float(mp.psi(1,0.25))
def Psi_dh(t):
    t=np.asarray(t,dtype=float); idx=np.searchsorted(L_dh,t,side='right')
    pr=np.where(idx>0,t*A_dh[np.maximum(idx-1,0)]-B_dh[np.maximum(idx-1,0)],0.0)
    return 0.5*t*(psi34+math.log(5/PI))+0.25*tri34-sumexp(t,1.5)-pr
ts=np.array([0.5,1.,2.,3.,5.,7.,10.,12.])
print('Psi_DH coefficient side',np.round(Psi_dh(ts),6)); print('Psi_DH zero side (130 zeros + mean tail)',np.round(Psi_zero(ts),6))
tt=np.arange(5e-4,14.0,5e-4); P=Psi_dh(tt)
def llogd(s):
    return mp.log(5/mp.pi)/2+mp.digamma((s+1)/2)/2+(-mp.log(5)*fmp(s)+5**(-s)*(mp.zeta(s,mp.mpf(1)/5,1)+kap*mp.zeta(s,mp.mpf(2)/5,1)-kap*mp.zeta(s,mp.mpf(3)/5,1)-mp.zeta(s,mp.mpf(4)/5,1)))/fmp(s)
for z in (2.,3.):
    I=5e-4*((np.exp(-z*tt)*P).sum()-0.5*math.exp(-z*tt[-1])*P[-1])
    print('second identity for DH, z=%g: Laplace of coefficient side %.8f ; (Lambda\'/Lambda)(1/2+z)/z^2 %s'%(z,I,mp.nstr(llogd(mp.mpf(1)/2+z)/z**2,10)))
print('first identity for DH at lambda=4: (Lambda\'/Lambda)(5/2)/(2*2)=',mp.nstr(llogd(mp.mpf(5)/2)/4,10),' sum over zeros to 200 + tail:',
  (1/(4+on**2)).sum()+sum(2*(1/(4+(g-1j*dl)**2)).real for g,dl in zip(gam,del_))+float(mp.quad(lambda g:mp.log(5*g/(2*mp.pi))/(2*mp.pi)/(4+g*g),[200,mp.inf])))
print('DH mean of Psi: sum gamma^-2 over zeros to 200 + tail',2*(1/on**2).sum()+4*sum((1/(g-1j*dl)**2).real for g,dl in zip(gam,del_))+(math.log(5*200/(2*PI))+1)/(PI*200))
tz=np.arange(0.05,80,0.002); Pz=Psi_zero(tz); Pz_on=Psi_zero(tz,offline=False)
neg=np.where(Pz<0)[0]
print('Psi_DH zero side: min on [0.1,12] %.5f; first negative at t=%s; book predicted scale (2 log gamma+log p1)/delta with p1=mean/2:'%(Pz[(tz>=.1)&(tz<=12)].min(),tz[neg[0]] if len(neg) else None))
print('Psi_DH exact coefficient side min on [0.1,12] %.5f at %.3f'%(P[(tt>=.1)&(tt<=12)].min(),tt[(tt>=.1)&(tt<=12)][P[(tt>=.1)&(tt<=12)].argmin()]))
# phase scan for DH
tgr=np.arange(1.0,215,0.001); thg=thdh(tgr)
def lev_c(c):
    k=np.arange(math.ceil(thg[0]/PI-c),math.floor(thg[-1]/PI-c)+1); return np.interp((k+c)*PI,thg,tgr)
all_z=np.sort(np.concatenate([on,gam,gam]))
def span(levs,zeros,a=10.,b=199.):
    ez=zeros[(zeros>=a)&(zeros<=b)]; el=levs[(levs>=a)&(levs<=b)]
    t=np.concatenate([ez,el]); s=np.concatenate([np.ones(len(ez)),-np.ones(len(el))])
    o=np.argsort(t,kind='mergesort'); D=np.concatenate([[0],np.cumsum(s[o])])
    return int(D.max()-D.min())
for nm,zz in (('all zeros',all_z),('on line only',on)):
    sp=np.array([span(lev_c(j/200),zz) for j in range(200)])
    print('DH phase scan on [10,199],',nm,dict(zip(*np.unique(sp,return_counts=True))))
    nn=np.arange(1,len(zz)+1); th=thdh(zz)/PI; Sp=nn-th; Sm=nn-1-th; m=(zz>=10)&(zz<=199)
    print('  N-theta/pi range: [%.3f,%.3f], width %.3f'%(Sm[m].min(),Sp[m].max(),Sp[m].max()-Sm[m].min()))
# zeta coefficient side and the kernels
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
def Psi_z(t):
    t=np.asarray(t,dtype=float); idx=np.searchsorted(LL,t,side='right')
    pp=np.where(idx>0,t*Az[np.maximum(idx-1,0)]-Bz[np.maximum(idx-1,0)],0.0)
    return 8*(np.cosh(t/2)-1)+0.5*t*(psi14-math.log(PI))+0.25*tri14-sumexp(t,0.5)-pp
def kernel_eigs(Pfun,ell,n):
    t=ell*np.arange(1,n+1)/n; P0=np.concatenate([[0.0],Pfun(t)])
    i=np.arange(1,n+1); Ii,Jj=np.meshgrid(i,i,indexing='ij'); G=P0[Ii]+P0[Jj]-P0[np.abs(Ii-Jj)]
    ev=np.linalg.eigvalsh(G); return ev
print('window length ell | zeta: least eigenvalue, largest | DH: least eigenvalue, number below -1e-9, largest (grid n=600; the kernel Psi(u)+Psi(v)-Psi(u-v) on [0,ell])')
for ell in (0.5,1.0,1.6,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,14.0):
    e1=kernel_eigs(Psi_z,ell,600); e2=kernel_eigs(Psi_dh,ell,600)
    print('%5.2f | %.3e %.3f | %.3e %d %.3f'%(ell,e1[0],e1[-1],e2[0],(e2<-1e-9).sum(),e2[-1]))
for ell in (4.0,6.0,8.0):
    e2=kernel_eigs(Psi_dh,ell,1500); print('DH check n=1500 ell=%g: least %.3e, number below -1e-9: %d'%(ell,e2[0],(e2<-1e-9).sum()))
print('t total',round(time.time()-T0,1))
