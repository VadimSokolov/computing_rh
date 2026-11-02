#!/usr/bin/env python3
"""N8 check 2: (a) counting function of the wall plus a real potential, by Pruefer phase and Sturm oscillation;
(b) reference levels corrected by a finite Euler product: is the rank one test stable?
Run: timeout 600 python3 -I N8_check2.py"""
import json, math, time
import numpy as np
HERE='./misc/research/directions/opus/'
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
def theta(t):
    t=np.asarray(t,dtype=float)
    return lg(0.25+0.5j*t).imag-0.5*t*math.log(PI)
def F01(b,z,K=70):
    term=np.ones_like(b,dtype=complex); s=term.copy()
    for k in range(K):
        term=term*z/((k+1)*(b+k)); s=s+term
    return s
# (a) Pruefer: u=R sin phi, u'=nu R cos phi, phi'=nu cos^2 phi+(nu^2-Q-V)/nu sin^2 phi, phi(x0)=0.
Tl=np.linspace(40.137,800.137,153); nu=(Tl/2)[:,None]
names=['V=0','k=+1 (2e^x)','k=-1 (-2e^x)','k=+3 (6e^x)','k=-3 (-6e^x)','6e^x sin3x','50 sin3x','200 exp(-(x-3)^2)','e^{1.5x} sin3x','x e^{x} sin3x']
def Vs(x):
    e=math.exp(x); s=math.sin(3*x)
    return np.array([0.,2*e,-2*e,6*e,-6*e,6*e*s,50*s,200*math.exp(-(x-3)**2),e**1.5*s,x*e*s])[None,:]
x0=math.log(2*PI); numax=float(nu.max()); xend=math.log(5*numax); h=0.04/numax
nst=int((xend-x0)/h)+1; h=(xend-x0)/nst
phi=np.zeros((len(Tl),len(names)))
def rhs(x,phi):
    Q=math.exp(2*x); act=(Q<25*nu*nu)
    c=np.cos(phi); s=np.sin(phi)
    return act*(nu*c*c+((nu*nu-Q-Vs(x))/nu)*s*s)
x=x0
for i in range(nst):
    k1=rhs(x,phi); k2=rhs(x+0.5*h,phi+0.5*h*k1); k3=rhs(x+0.5*h,phi+0.5*h*k2); k4=rhs(x+h,phi+h*k3)
    phi=phi+h/6*(k1+2*k2+2*k3+k4); x+=h
ps=math.atan(1/math.sqrt(24))
NV=np.floor((phi+ps)/PI)
nuu=Tl/2
exact=np.floor((lg(1+1j*nuu).imag-nuu*math.log(PI)-np.angle(F01(1+1j*nuu,PI**2)))/PI)
print('Pruefer against the exact rule N_P=floor(vartheta_P/pi) of the book (V=0): max abs diff',np.abs(NV[:,0]-exact).max(),' steps',nst)
th=theta(Tl)/PI
hi=Tl>420
print('potential | min and max of N_V(T)-theta(T)/pi on T in [40,800] | mean on T<420 | mean on T>420 | predicted shift -k/2 for Morse')
for j,nm in enumerate(names):
    d=NV[:,j]-th
    print('%-20s min %8.3f max %8.3f mean_low %8.3f mean_high %8.3f'%(nm,d.min(),d.max(),d[~hi].mean(),d[hi].mean()))
print('t',round(time.time()-T0,1))
# (b) finite Euler product phase
def sieve(N):
    s=np.ones(N+1,bool); s[:2]=False
    for i in range(2,int(N**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
t=np.arange(30,680.00005,0.001); th=theta(t); Nz=np.searchsorted(Z,t,side='right').astype(float)
S=Nz-th/PI-1
corr=np.zeros_like(t); pr=[int(p) for p in sieve(1100)]
Plist=[1,2,3,5,7,13,31,101,1009]
done=0
print('P | min span over 40 phases (best c) | max span | range and rms of S-S_P | rms of S_P | share of the grid where the phase decreases')
for P in Plist:
    while done<len(pr) and pr[done]<=P:
        p=pr[done]; lp=math.log(p); kmax=max(1,int(2*math.log(1e4)/lp))
        for k in range(1,kmax+1):
            corr-=np.sin(k*t*lp)/(k*p**(k/2))
        done+=1
    Th=th+corr
    best=(99,0); worst=0
    for c in np.arange(0,1,1/40):
        env=np.maximum.accumulate(np.floor(Th/PI-c)); D=Nz-env; sp=int(D.max()-D.min())
        if sp<best[0]: best=(sp,c)
        worst=max(worst,sp)
    E=S-corr/PI
    print('%5d | %d (c=%.3f) | %d | [%.3f,%.3f] rms %.4f | %.4f | %.4f'%(P,best[0],best[1],worst,E.min(),E.max(),E.std(),(corr/PI).std(),(np.diff(Th)<0).mean()))
print('t total',round(time.time()-T0,1))
