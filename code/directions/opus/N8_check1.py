#!/usr/bin/env python3
"""N8 check 1: independent recheck of the numbers of the hitting time note.
No code from misc/reviews is used.  Run: timeout 600 python3 -I N8_check1.py"""
import json, math, time
import numpy as np, mpmath as mp
HERE='./misc/research/directions/opus/'
Z=np.array([float(x) for x in json.load(open(HERE+'N8_zeros.json'))])
PI=math.pi
T0=time.time()
B=[1/12,-1/360,1/1260,-1/1680,1/1188,-691/360360,1/156]
def lg(z,shift=12):
    z=np.asarray(z,dtype=complex); s=np.zeros_like(z); w=z.copy()
    for k in range(shift):
        s=s-np.log(w); w=w+1
    iw=1/w; iw2=iw*iw; ser=0; p=iw
    for b in B:
        ser=ser+b*p; p=p*iw2
    return s+(w-0.5)*np.log(w)-w+0.5*math.log(2*PI)+ser
def digam(z,shift=12):
    z=np.asarray(z,dtype=complex); s=0; w=z.copy()
    for k in range(shift):
        s=s-1/w; w=w+1
    x=1/(w*w)
    return s+np.log(w)-0.5/w-x*(1/12-x*(1/120-x*(1/252-x*(1/240-x*(1/132)))))
def theta(t):
    t=np.asarray(t,dtype=float)
    return lg(0.25+0.5j*t).imag-0.5*t*math.log(PI)
print('theta check',float(theta(100.0)),float(mp.siegeltheta(100)),float(theta(500.0)),float(mp.siegeltheta(500)))
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
for a,t in ((2.25,10.),(0.,40.),(1.25,300.)):
    print('K check',a,t,complex(Ksc(a,np.array([t]))[0]),complex(mp.besselk(mp.mpc(a,t/2),2*mp.pi)*mp.exp(mp.pi*t/4)))
def levels(fun,lo=1.0,hi=705.0,h=0.01):
    t=np.arange(lo,hi,h); v=fun(t); idx=np.where(v[:-1]*v[1:]<0)[0]
    a=t[idx].copy(); b=t[idx+1].copy(); fa=v[idx].copy()
    for _ in range(45):
        m=0.5*(a+b); fm=fun(m); left=fa*fm<=0
        b=np.where(left,m,b); a=np.where(left,a,m); fa=np.where(left,fa,fm)
    return 0.5*(a+b)
base=levels(lambda t:Ksc(0.,t).real)
two=levels(lambda t:(2*PI**2*Ksc(2.25,t)-3*PI*Ksc(1.25,t)).real)
xis=levels(lambda t:Ksc(2.25,t).real)
print('base first',base[:4],'nu1',base[0]/2,'count',len(base))
print('two term first',two[:4],'count',len(two))
print('Xi* first',xis[:4],'count',len(xis))
for nm,L in (('base',base),('two',two),('Xi*',xis)):
    fr=(theta(L[-40:])/PI)%1
    print(nm,'theta/pi mod 1 at the last 40 levels: mean %.5f min %.5f max %.5f'%(fr.mean(),fr.min(),fr.max()))
tg=np.arange(6.5,725,0.001); thg=theta(tg)
def lev_c(c):
    k=np.arange(math.ceil(thg[0]/PI-c),math.floor(thg[-1]/PI-c)+1)
    return np.interp((k+c)*PI,thg,tg)
gram=lev_c(0.0)
print('gram points',gram[:4],'g126 index check: first >280:',gram[gram>280][:2])
def table(levs,name):
    for start in (0.,14.,17.):
        L=levs[(levs>=start)&(levs<=680.0)]
        cnt=np.searchsorted(Z,L[1:])-np.searchsorted(Z,L[:-1])
        bad=np.where(cnt!=1)[0]
        print(name,'start>=%g: intervals %d, without exactly one zero %d, first failure at left end %.2f, histogram %s, first level %.3f, last %.3f'%(start,len(L)-1,len(bad),L[bad[0]] if len(bad) else float('nan'),np.bincount(cnt),L[0],L[-1]))
print('== table of the note: (Gram 398,22,280.80), (two term 397,207,20.65), (Liouville 397,188,26.05)')
table(gram,'Gram'); table(two,'two term'); table(base,'Liouville'); table(xis,'Xi*')
def span(levs,zeros,a=30.,b=680.):
    ez=zeros[(zeros>=a)&(zeros<=b)]; el=levs[(levs>=a)&(levs<=b)]
    t=np.concatenate([ez,el]); s=np.concatenate([np.ones(len(ez)),-np.ones(len(el))])
    o=np.argsort(t,kind='mergesort'); D=np.concatenate([[0],np.cumsum(s[o])])
    return int(D.max()-D.min())
print('actual level spans on [30,680]: Gram',span(gram,Z),'two',span(two,Z),'base',span(base,Z),'Xi*',span(xis,Z))
n=np.arange(1,len(Z)+1); th=theta(Z)/PI; Sp=n-th-1; Sm=n-2-th
for a,b in ((30.,680.),(14.,680.),(14.,706.)):
    m=(Z>=a)&(Z<=b)
    print('S range on [%g,%g]: sup %.4f at gamma=%.3f (n=%d), inf %.4f at gamma=%.3f (n=%d), width %.4f'%(a,b,Sp[m].max(),Z[m][Sp[m].argmax()],n[m][Sp[m].argmax()],Sm[m].min(),Z[m][Sm[m].argmin()],n[m][Sm[m].argmin()],Sp[m].max()-Sm[m].min()))
for grid in (200,20000):
    sp=np.array([span(lev_c(j/grid),Z) for j in range(grid)])
    print('phase scan grid',grid,'span counts',dict(zip(*np.unique(sp,return_counts=True))),'fractions',{int(k):float((sp==k).mean()) for k in np.unique(sp)})
    if grid==20000:
        ch=np.where(np.diff(sp)!=0)[0]
        print('  breakpoints c:',[round((j+1)/grid,4) for j in ch][:40],' span at c=0:',sp[0])
print('bad intervals against the phase c (levels theta=(k+c)pi, window [30,680]):')
for c in (0,.05,.1,.2,.3,.4,.5,.6,.625,.7,.8,.9,.95):
    L=lev_c(c); L=L[(L>=30)&(L<=680)]
    cnt=np.searchsorted(Z,L[1:])-np.searchsorted(Z,L[:-1])
    print('  c=%.3f intervals %d bad %d span %d'%(c,len(L)-1,(cnt!=1).sum(),span(lev_c(c),Z)))
print('t so far',round(time.time()-T0,1))
# == prime side of Psi
def sieve(N):
    s=np.ones(N+1,bool); s[:2]=False
    for i in range(2,int(N**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
N=9000000
pr=sieve(N); L1=np.log(pr.astype(float)); W1=L1/np.sqrt(pr.astype(float))
eL=[];eW=[]
for p in pr[pr<=int(N**.5)]:
    p=int(p); q=p*p
    while q<=N:
        eL.append(math.log(q)); eW.append(math.log(p)/math.sqrt(q)); q*=p
LL=np.concatenate([L1,np.array(eL)]); WW=np.concatenate([W1,np.array(eW)]); o=np.argsort(LL); LL=LL[o]; WW=WW[o]
A=np.cumsum(WW); Bc=np.cumsum(WW*LL)
psi14=float(mp.digamma(0.25)); tri14=float(mp.psi(1,0.25))
def sumexp(t,a0=0.5,Kmax=200000):
    t=np.asarray(t,dtype=float); out=np.zeros_like(t)
    for k in range(Kmax):
        a=2*k+a0; j=np.searchsorted(t,45/a)
        if j==0: break
        out[:j]+=np.exp(-a*t[:j])/a**2
    return out
def Psi_exact(t):
    """t ascending, positive.  Psi from the prime side, no zero is used (valid for t<=log N)."""
    t=np.asarray(t,dtype=float); idx=np.searchsorted(LL,t,side='right')
    pr_=np.where(idx>0,t*A[np.maximum(idx-1,0)]-Bc[np.maximum(idx-1,0)],0.0)
    return 8*(np.cosh(t/2)-1)+0.5*t*(psi14-math.log(PI))+0.25*tri14-sumexp(t)-pr_
xs=np.array(sorted([0.2,0.5,math.log(2),1.0,1.5,math.log(5),2.0,3.0,4.0,6.0]))
print('D=2Psi against the book Table yo:tab:D (prime side):',[(round(float(x),4),round(float(2*v),7)) for x,v in zip(xs,Psi_exact(xs))])
tgr=np.arange(0.1,12+1e-9,0.0005)
g400=Z[:400]
P400=np.zeros_like(tgr)
for g in g400: P400+=2*(1-np.cos(g*tgr))/g**2
Tc=g400[-1]; tail=(math.log(Tc/(2*PI))+1)/(PI*Tc)
Pex=Psi_exact(tgr)
print('Psi on [0.1,12]: from first 400 zeros min %.5f at t=%.4f; plus mean tail %.5f min %.5f; exact prime side min %.5f at t=%.4f; endpoint Psi(0.1)=%.5f'%(P400.min(),tgr[P400.argmin()],tail,(P400+tail).min(),Pex.min(),tgr[Pex.argmin()],Pex[0]))
tgr2=np.arange(0.05,16+1e-9,0.0005); Pex2=Psi_exact(tgr2)
print('exact on [0.05,16]: min D=2Psi %.5f at x=%.3f (book 0.0550 at 8.04); mean of D on [2,16] %.6f (book 0.092409)'%(2*Pex2.min(),tgr2[Pex2.argmin()],2*Pex2[tgr2>=2].mean()))
mp.mp.dps=30
def xld(s):
    s=mp.mpmathify(s)
    return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,1,1)/mp.zeta(s)
ss=mp.diff(lambda z:xld(mp.mpf(1)/2+z),0)
print('sum over all gamma of gamma^-2 = d/dz (xi\'/xi)(1/2+z) at 0 =',mp.nstr(ss,12),' 2+gammaE-log4pi=',mp.nstr(2+mp.euler-mp.log(4*mp.pi),12),' partial 400 zeros',2*(1/g400**2).sum())
for lam in (1.,10.):
    r=mp.sqrt(lam); lhs=xld(mp.mpf(1)/2+r)/(2*r)
    Tz=float(Z[-1]); rhs=(1/(lam+Z**2)).sum()+float(mp.quad(lambda g:mp.log(g/(2*mp.pi))/(2*mp.pi)/(lam+g*g),[Tz,mp.inf]))
    print('first identity lambda=%g: LHS %s RHS(420 zeros + mean tail) %.8f'%(lam,mp.nstr(lhs,10),rhs))
h=5e-4; tL=h*np.arange(1,32001); PL=Psi_exact(tL)
for z in (1.,2.):
    I=h*((np.exp(-z*tL)*PL).sum()-0.5*np.exp(-z*tL[-1])*PL[-1])+0.5*h*h*0  # trapezoid, the first panel contributes ~0.5 h exp(-z h) Psi(h), included as half of first node
    I=I+0.0
    bar=float(ss); I2=I+bar*math.exp(-z*16)/z
    print('second identity z=%g: int e^{-zt}Psi (prime side, t<=16, mean tail beyond) %.8f; (xi\'/xi)(1/2+z)/z^2 = %s; 2/z sum 1/(z^2+gamma^2)= %.8f'%(z,I2,mp.nstr(xld(mp.mpf(1)/2+z)/z**2,10),2/z*((1/(z*z+Z**2)).sum()+float(mp.quad(lambda g:mp.log(g/(2*mp.pi))/(2*mp.pi)/(z*z+g*g),[float(Z[-1]),mp.inf])))))
def Arch(b):
    R=math.sqrt(45/b); hh=0.005; r=np.arange(0,R,hh); f=np.exp(-b*r*r)*digam(0.25+0.5j*r).real
    return 2*math.exp(b/4)-math.log(PI)/math.sqrt(4*PI*b)+(hh*(f.sum()-0.5*f[0]))/PI
def Pr(b):
    return -2/math.sqrt(4*PI*b)*(WW*np.exp(-LL*LL/(4*b))).sum()
print('b, W(zeros), (A+P)/2, P/2, A/2-W_wall(base), A/2-W_twoterm, A/2-W_Xi*')
for b in (0.0015,0.002,0.003,0.005,0.01,0.02,0.05,0.1,0.2):
    W=np.exp(-b*Z*Z).sum(); Aa=Arch(b); Pp=Pr(b)
    print('%.4f %.10f %.10f %.3e %.8f %.8f %.8f'%(b,W,0.5*(Aa+Pp),0.5*Pp,0.5*Aa-np.exp(-b*base*base).sum(),0.5*Aa-np.exp(-b*two*two).sum(),0.5*Aa-np.exp(-b*xis*xis).sum()))
nn=np.arange(1,2000001.)
for lam in (1.0,):
    prod=np.exp(-np.log1p(2*lam/(nn*PI)**2).sum()); cl=math.sqrt(2*lam)/math.sinh(math.sqrt(2*lam))
    print('Williams: product of the note field at lambda=1 %.8f; sqrt(2l)/sinh sqrt(2l)= %.8f; its square %.8f; mean of H_W sum 2/(n pi)^2= %.6f'%(prod,cl,cl*cl,(2/(nn*PI)**2).sum()))
rng=np.random.default_rng(8)
g=g400; dm=g[:,None]-g[None,:]; sm=g[:,None]+g[None,:]
def snc(x,T):
    y=x*T; out=np.ones_like(y); m=y!=0; out[m]=np.sin(y[m])/y[m]; return out
def csn(x,T):
    y=x*T; out=np.zeros_like(y); m=y!=0; out[m]=(1-np.cos(y[m]))/y[m]; return out
for rep in range(3):
    a,b=rng.normal(size=(2,400)); cc=a/g; dd=b/g; target=((a*a+b*b)/(2*g*g)).sum()
    res=[]
    for T in (1e2,1e3,1e4,1e5,1e6):
        Mcc=0.5*(snc(dm,T)+snc(sm,T)); Mss=0.5*(snc(dm,T)-snc(sm,T)); Mcs=0.5*(csn(sm,T)+csn(-dm,T))
        val=cc @ Mcc @ cc+dd @ Mss @ dd+2*cc @ Mcs @ dd
        res.append('%.1e'%abs(val-target))
    print('Besicovitch draw',rep,'target %.6f, |mean over [0,T]-target| at T=1e2..1e6:'%target,res)
print('t total',round(time.time()-T0,1))
