# Own checks for polson_above_one.tex. Run: timeout 600 python3 -I /tmp/t3rev/t3a.py
import json, numpy as np, mpmath as mp
from scipy.special import zeta, gammaln
from scipy.interpolate import CubicSpline
mp.mp.dps=30
def xi(s):
    if s==1: return mp.mpf(1)/2
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def logxi_np(t):
    t=np.asarray(t,float); out=np.full_like(t,np.log(.5)); m=np.abs(t-1)>1e-13; v=t[m]
    out[m]=np.log(v*(v-1)/2)-v/2*np.log(np.pi)+gammaln(v/2)+np.log(zeta(v,1)); return out
print("== A. table, true values against note's rule")
xs=np.concatenate([np.linspace(0,5,201),np.linspace(5.25,200,400)])
for a in [1.0,1.5,2.0,3.0]:
    c=mp.diff(lambda s: mp.log(xi(s)),mp.mpf(a)+mp.mpf('0.0'),1) if a!=1 else mp.mpf(0)
    if a==1: c=mp.diff(lambda s: mp.log(xi(s)),mp.mpf(1),1,h=mp.mpf('1e-6'),direction=0) if False else (mp.diff(lambda s: mp.log(xi(s)),mp.mpf(1)+mp.mpf('1e-8')))
    # phi''' at 0 by g derivatives
    aa=mp.mpf(a)+(mp.mpf('1e-9') if a==1 else 0)
    g=lambda u: mp.log(xi(aa+u))-mp.log(xi(aa))
    g1,g2,g3=[mp.diff(g,mp.mpf('0.0'),k) for k in (1,2,3)]
    ph3=-(g3-3*g1*g2+g1**3)
    # Gauss-Legendre panels on [0,60]
    gx,gw=np.polynomial.legendre.leggauss(12)
    def f(x):
        p=min(0.05,0.5/max(x,1)); e=np.arange(0,60+p,p); mid=(e[1:]+e[:-1])/2; hw=(e[1:]-e[:-1])/2
        u=(mid[:,None]+hw[:,None]*gx); w=hw[:,None]*gw
        ph=np.exp(logxi_np(np.array([a]))[0]-logxi_np(a+u))
        return float(np.sum(w*ph*np.cos(u*x))/np.pi)
    row=[]
    for x in (50,100,200):
        xp=xs[np.argmin(abs(xs-x))]
        h=1e-3; alias=(x*0+xp*h/2)**2/np.sin(xp*h/2)**2
        row.append((x,round(float(xp),4),np.pi*x*x*f(x),float(g1+ph3/x**2), float(g1)*(x/xp)**2*alias))
    print("alpha",a,"c",mp.nstr(g1,10),"phi'''(0)",mp.nstr(ph3,6))
    for r in row: print("  x=%d gridx=%.4f true pi x^2 f=%.8f asym c+phi3/x^2=%.8f predicted-note-entry=%.6f"%r)
print("== B. Prop V: Thorin integral vs closed form, tail")
for th,ga,s in [(0.5,14.13,3.0),(0.9,21.0,50.0),(2,0,1.0)]:
    u=lambda w: (th/(th**2+(w-ga)**2)+th/(th**2+(w+ga)**2))/mp.pi
    pts=sorted(set([0,max(ga-1,0.0),ga,ga+1,ga+50,1e4]))
    I=mp.quad(lambda w: mp.log(1+s/w**2)*u(w),pts)+mp.quad(lambda w: mp.log(1+s/w**2)*u(w),[1e4,mp.inf])
    print(th,ga,s,"int",mp.nstr(I,12),"closed",mp.nstr(mp.log(((th+mp.sqrt(s))**2+ga**2)/(th**2+ga**2)),12))
print("== C. thm polson product vs direct")
Z=[mp.mpf(z) for z in json.load(open('./misc/research/directions/opus/N8_zeros.json'))]
print("zeros",len(Z),"last",mp.nstr(Z[-1],8))
T=(Z[-1]+mp.zetazero(len(Z)+1).imag)/2
dens=lambda x: mp.log(x/(2*mp.pi))/(2*mp.pi)
for a,s in [(1,1),(1,25),(1.5,4),(3,9)]:
    th=mp.mpf(a)-mp.mpf(1)/2; rs=mp.sqrt(s)
    term=lambda g: mp.log(((th+rs)**2+g**2)/(th**2+g**2))
    L=sum(term(g) for g in Z)+mp.quad(lambda x: dens(x)*term(x),[T,10*T,mp.inf])
    print(a,s,"product",mp.nstr(mp.exp(-L),10),"direct",mp.nstr(xi(mp.mpf(a))/xi(mp.mpf(a)+rs),10))
S1=sum(1/(mp.mpf(1)/4+g**2) for g in Z)+mp.quad(lambda x: dens(x)/(mp.mpf(1)/4+x**2),[T,10*T,mp.inf])
print("sum 2theta/(theta^2+g^2) at alpha=1",mp.nstr(S1,10))
print("== D. noindep")
g1=Z[0]
for m in (2,3):
    print(m,"m g1",mp.nstr(m*g1,10),"|zeta(1/2+i m g1)|",mp.nstr(abs(mp.zeta(mp.mpc(.5,m*g1))),8),"nearest ordinate gap",mp.nstr(min(abs(m*g1-g) for g in Z),6))
print("== E. harmonic")
k2=mp.diff(lambda s: mp.log(xi(s)),mp.mpf(.5),2); k4=mp.diff(lambda s: mp.log(xi(s)),mp.mpf(.5),4)
print("k2X",mp.nstr(k2,12),"k4X",mp.nstr(k4,12),"z2k2",mp.nstr(mp.zeta(2)*k2,12),"z4k4",mp.nstr(mp.zeta(4)*k4,12))
print("exc kurt X",mp.nstr(k4/k2**2,8),"S",mp.nstr(mp.zeta(4)*k4/(mp.zeta(2)*k2)**2,8),"ratio",mp.nstr(mp.zeta(4)/mp.zeta(2)**2,6))
print("sd X",mp.nstr(mp.sqrt(k2),8),"sd S",mp.nstr(mp.sqrt(mp.zeta(2)*k2),8))
print("1-1/zeta2",mp.nstr(1-1/mp.zeta(2),6),"gross k>=2 squarefree over k=1: zeta2/zeta4-1",mp.nstr(mp.zeta(2)/mp.zeta(4)-1,6))
print("standardised kappa ratios zeta(2m)/zeta(2)^m",[mp.nstr(mp.zeta(2*m)/mp.zeta(2)**m,5) for m in (2,3,4,5)])
mp.mp.dps=18
h=0.01; ts=np.arange(0,70+h,h); xi0=xi(mp.mpf(.5))
phi=np.array([float(mp.re(xi(mp.mpc(.5,t))/xi0)) if t>0 else 1.0 for t in ts])
sp=CubicSpline(ts,phi)
def phiS(t,N=400):
    t=np.asarray(t); p=np.ones_like(t)
    for n in range(1,N+1): p=p*np.where(t/n<=70,sp(np.minimum(t/n,70)),0.0)
    tail=float(mp.zeta(2,N+1)); return p*np.exp(-float(k2)*t**2*tail/2)
tg=np.arange(0,70,0.005); pX=sp(tg); pS=phiS(tg)
sX=float(mp.sqrt(k2)); sS=float(mp.sqrt(mp.zeta(2)*k2))
xg=np.linspace(-3,3,3001); dx=xg[1]-xg[0]
def dens_from(p,t,x,scale=1.0):
    dt=t[1]-t[0]; w=np.full_like(t,dt); w[0]=dt/2
    return np.array([np.sum(w*p*np.cos(t*y*scale)) for y in x])/np.pi*scale
fX=dens_from(pX,tg,xg); fS=dens_from(pS,tg,xg)
print("mass",fX.sum()*dx,fS.sum()*dx,"min",fX.min(),fS.min())
print("TV(X,S) raw",0.5*np.sum(abs(fX-fS))*dx,"W1 raw",np.sum(abs(np.cumsum(fX-fS)*dx))*dx)
zg=np.linspace(-6,6,3001); dz=zg[1]-zg[0]
fXs=dens_from(pX,tg,zg,sX); fSs=dens_from(pS,tg,zg,sS)
print("TV standardised",0.5*np.sum(abs(fXs-fSs))*dz,"TV(S std, N01)",0.5*np.sum(abs(fSs-np.exp(-zg**2/2)/np.sqrt(2*np.pi)))*dz,"TV(X std,N01)",0.5*np.sum(abs(fXs-np.exp(-zg**2/2)/np.sqrt(2*np.pi)))*dz)
print("density at 0 std X,S,N",fXs[1500],fSs[1500],1/np.sqrt(2*np.pi))
# balanced Mobius identity, no zeros: phi_X * prod_{mu=-1} phi_S(t/k) = prod_{mu=+1} phi_S(t/k)
def mu(n):
    r=1;p=2
    while p*p<=n:
        if n%p==0:
            n//=p
            if n%p==0: return 0
            r=-r
        p+=1
    return -r if n>1 else r
K=2000
for t0 in (2.0,5.0,10.0):
    lhs=float(sp(t0)); 
    L=0.0
    for k in range(1,K+1):
        m=mu(k)
        if m: L+=m*np.log(float(phiS(np.array([t0/k]))[0]))
    print('t',t0,'phi_X',lhs,'exp(sum mu(k) log phi_S(t/k))',np.exp(L))
