import numpy as np, mpmath as mp, json, os
from scipy.special import loggamma, digamma, polygamma
from scipy.optimize import brentq
mp.mp.dps=15
fn='/tmp/n8rev/zeros.npy'
if os.path.exists(fn): g=np.load(fn)
else:
    g=np.array([float(mp.zetazero(n).imag) for n in range(1,421)]); np.save(fn,g)
print("gamma_400,420,289,379",g[399],g[419],g[288],g[378])
th=lambda t: np.imag(loggamma(0.25+0.5j*np.asarray(t,float)))-np.asarray(t,float)/2*np.log(np.pi)
print("S band: S(gamma_m^+)=m-1-th/pi; S^-=S^+-1")
m=np.arange(1,421);Sp=m-1-th(g)/np.pi;Sm=Sp-1
mask=(g>30)&(g<680)
print(" sup S",Sp[mask].max(),"at",g[mask][Sp[mask].argmax()],"inf S",Sm[mask].min(),"at",g[mask][Sm[mask].argmin()])
sup,inf=Sp[mask].max(),Sm[mask].min();print(" breakpoints -inf-1",-inf-1," 2-sup",2-sup)
def inv(v,lo,hi):
    return brentq(lambda t:float(th(t))-v,lo,hi,xtol=1e-11)
def phase_levels(c,tmax=720,tmin=8):
    k0=int(np.floor(float(th(tmin))/np.pi-c))+1;L=[]
    k=k0
    while True:
        v=(k+c)*np.pi
        if float(th(tmax))<v:break
        L.append(inv(v,tmin-1,tmax+1));k+=1
    return np.array(L)
# wall levels via book's exact rule, theta_P(nu)=n pi
def f01(b,z,N=200):
    s=1+0j;t=1+0j
    for j in range(N):
        t=t*z/((b+j)*(j+1));s+=t
    return s
def thP(nu):
    return np.imag(loggamma(1+1j*nu))-nu*np.log(np.pi)-np.angle(f01(1+1j*nu,np.pi**2))
wall=[];n=1
while True:
    v=n*np.pi
    hi=2*np.pi
    while thP(hi)<v: hi+=1.0
    if hi>400:break
    wall.append(2*brentq(lambda x:thP(x)-v,hi-1.0 if n>1 else 0.1,hi,xtol=1e-12));n+=1
wall=np.array(wall);print("wall first",wall[:5].round(3),"K check:",float(abs(mp.besselk(1j*wall[9]/2,2*mp.pi))*mp.exp(mp.pi*wall[9]/4)))
# gram
gr=[]
k=0
while True:
    v=k*np.pi
    if float(th(720))<v:break
    gr.append(inv(v,9,721));k+=1
gr=np.array(gr);print("g0",gr[0],"g126",gr[126],gr[127])
def table(L,name,gamma=g[:400],top=None):
    top=gamma[-1]
    L=L[(L>=17)&(L<=top)]
    cnt=np.histogram(gamma,bins=L)[0]
    bad=np.where(cnt!=1)[0]
    # span on [30,680]
    return name,L[0],len(L)-1,len(bad),(L[bad[0]] if len(bad) else None),np.bincount(np.minimum(cnt,2)).tolist()
def span(L,a=30.,b=680.,G=g):
    ev=np.sort(np.concatenate([L[(L>a)&(L<b)],G[(G>a)&(G<b)]]))
    pts=np.concatenate([[a],(ev[:-1]+ev[1:])/2,[b]])
    d=np.searchsorted(G,pts,side='right')-np.searchsorted(L,pts,side='right')
    return d.max()-d.min(),d.min(),d.max()
for L,nm in [(gr,'Gram'),(wall,'wall'),(phase_levels(.5),'c=1/2'),(phase_levels(.625),'c=5/8')]:
    print(" ",table(L,nm),"span",span(L))
print(" gram intervals g0 to 700, wall levels per interval:",end=' ')
Gg=gr[gr<=700];print(len(Gg)-1,np.unique(np.histogram(wall,bins=Gg)[0]))
for tt in [100,680]:
    i=np.argmin(abs(wall-tt));x=th(wall[i])/np.pi;print(f" wall level {wall[i]:.2f} theta/pi mod 1 - 5/8 = {(x%1)-.625:.4f}")
i=np.argmin(abs(wall-29.37));print(" near 30: level",wall[i],"ideal 5/8 level",phase_levels(.625,tmin=10)[:3],"ideal 1/2",phase_levels(.5,tmin=10)[:3])
lastm=np.mean(th(wall[-40:])/np.pi%1);print(" mean last 40 wall levels theta/pi mod 1:",lastm)
print("span by c scan");cs=np.linspace(0,1,2001)[:-1];sp=np.array([span(phase_levels(c,tmin=20))[0] for c in cs[::20]])
print(" c grid 100: spans",np.unique(sp,return_counts=True));c3=cs[::20][sp==3];print(" span3 c range:",c3[c3<.5].max(),c3[c3>.5].min())
for c in [0.145,0.146,0.835,0.836,0.837]: print("  c",c,span(phase_levels(c,tmin=20))[0])
print("Theta_P")
primes={2:[2],3:[2,3]}
def ThP(t,P):
    t=np.asarray(t,float);s=th(t)
    for p in primes.get(P,[]):
        for k in range(1,60):
            s=s-np.sin(k*t*np.log(p))/(k*p**(k/2))
    return s
for P in [0,2,3]:
    e=np.arange(8,421)-ThP(g[7:],P)/np.pi
    r=np.array([e[:i+1].max()-e[:i+1].min() for i in range(len(e))])
    h=g[7:][np.argmax(r>=1)] if (r>=1).any() else None
    tg=np.arange(30,680,0.001);dd=np.diff(ThP(tg,P));
    print(f" P={P}: width {r[-1]:.4f}, first height {h}, share decreasing {np.mean(dd<0):.4f}")
    # span by phases with running max
    tg2=np.arange(30,680,0.002);T=np.maximum.accumulate(ThP(tg2,P));Ng=np.searchsorted(g,tg2,side='right')
    cnt=0
    for c in np.arange(200)/200:
        d=Ng-np.floor(T/np.pi-c);cnt+=(d.max()-d.min()==1)
    print(f"   phases of 200 with span 1: {cnt}")
print("P=3,c=0 levels and rank one matrix; rank=span matrices")
def levels_theta(P,c,a,b):
    tg=np.arange(a,b,0.001);v=ThP(tg,P)/np.pi-c;R=np.maximum.accumulate(v)
    idx=np.where(np.floor(R[1:])>np.floor(R[:-1]))[0];out=[]
    for i in idx:
        out.append(brentq(lambda t:float(ThP(t,P))/np.pi-c-np.floor(R[i+1]),tg[i],tg[i+1]))
    return np.array(out)
def interlace_step(x,y):
    # x_j<=y_j<=x_{j+1}: return symmetric F with spec(diag x+vv^T)=y
    keep=np.abs(x-y)>1e-9*np.maximum(x,1)
    # delete values common
    xs=x.copy();ys=y.copy()
    common=np.intersect1d(np.round(x,6),np.round(y,6))
    kx=~np.isin(np.round(x,6),common);ky=~np.isin(np.round(y,6),common)
    a=x[kx];b=y[ky];v=np.zeros(len(x))
    w=np.zeros(len(a))
    for i in range(len(a)):
        w[i]=np.sum(np.log(np.abs(b-a[i])))-np.sum(np.log(np.abs(np.delete(a,i)-a[i])))
    v[kx]=np.exp(w/2)
    return v
def count(L,pts): return np.searchsorted(L,pts,side='right')
def build(al,be):
    """al,be sorted, equal length; returns matrix A=diag al+F, in the Remark's scheme"""
    n=len(al);pts=np.sort(np.concatenate([al,be]))
    Na=count(al,pts);Nb=count(be,pts);M=np.minimum(Na,Nb);d=Nb-Na;p=(-d).max();q=d.max()
    def lst(C):
        return np.array([pts[np.searchsorted(C,j,side='left')] for j in range(1,n+1)])
    Cs=[np.minimum(Na,M+s) for s in range(p,-1,-1)]  # from N_alpha down to M
    Es=[np.minimum(Nb,M+s) for s in range(0,q+1)]
    seq=Cs+Es[1:];lists=[lst(C) for C in seq]
    A=np.diag(al.astype(float));Q=np.eye(n);ok=True
    for i in range(len(lists)-1):
        x,y=lists[i],lists[i+1];sgn=1
        if np.sum(y)<np.sum(x): x,y=y,x;sgn=-1
        ok&=bool(np.all(y>=x-1e-9) and np.all(y[:-1]<=x[1:]+1e-9))
        v=interlace_step(x,y)
        # eigenbasis of current list
        lam,U=np.linalg.eigh(A)
        # lam ascending = lists[i]
        if sgn==1: A=A+U @ np.outer(v,v) @ U.T
        else:
            # from y (current) down to x: go via x then add: A=U diag(x)U^T with x+vv^T has spectrum y; use inverse construction
            B=np.diag(y)-np.diag(y)  # placeholder
            # spec(diag y - uu^T)=x: weights from secular equation with roles swapped
            u=interlace_step_neg(y,x);A=A-U @ np.outer(u,u) @ U.T
    return A,p,q,ok
def interlace_step_neg(y,x):
    common=np.intersect1d(np.round(x,6),np.round(y,6))
    ky=~np.isin(np.round(y,6),common);kx=~np.isin(np.round(x,6),common)
    a=y[ky];b=x[kx];u=np.zeros(len(y));w=np.zeros(len(a))
    for i in range(len(a)):
        w[i]=np.sum(np.log(np.abs(b-a[i])))-np.sum(np.log(np.abs(np.delete(a,i)-a[i])))
    u[ky]=np.exp(w/2);return u
def rec(L,name,a=30,b=680):
    z=g[(g>a)&(g<b)];l=L[(L>a)&(L<b)];n=min(len(z),len(l));z=z[:n];l=l[:n]
    try:
        A,p,q,ok=build(l**2,z**2);ev=np.linalg.eigvalsh(A);err=np.max(np.abs(np.sqrt(ev)-z))
        F=A-np.diag(l**2);sv=np.linalg.svd(F,compute_uv=False)
        print(f" {name}: n={n} p={p} q={q} interlace_ok={ok} rank(1e-6 rel)={(sv>1e-8*sv[0]).sum()} recovery err in t {err:.2e}, z range {z[0]:.3f} {z[-1]:.3f}")
    except Exception as ex: print(" ",name,"failed",ex)
rec(gr,'Gram');rec(wall,'wall');rec(phase_levels(.5),'c=1/2');rec(phase_levels(.625),'c=5/8')
L3=levels_theta(3,0,20,720);print(" P=3 c=0 first levels",L3[:4].round(3)," span",span(L3),"bad on [30,680]",(np.histogram(g,bins=L3[(L3>30)&(L3<680)])[0]!=1).sum());rec(L3,'P=3 c=0')
print("heat trace at b")
for b in [0.005,0.002,0.0003]:
    W=np.sum(np.exp(-g**2*b));Ww=np.sum(np.exp(-wall**2*b))
    P=-(1/(2*np.sqrt(np.pi*b)))*sum(np.log(p)/np.sqrt(n)*np.exp(-np.log(n)**2/(4*b)) for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3)])
    print(f" b={b}: W-W_wall={W-Ww:.4f} (420 zeros, tail matters if small b) P/2={P:.3g}")
print("Laplace (a) at z=1,2")
for z in [1,2]:
    lhs=sum(1/(z*(z**2+x**2)) for x in g)*2
    rhs=(2*z/(z**2-.25)-.5*np.log(np.pi)+.5*digamma(.25+z/2)+float(mp.diff(mp.zeta,.5+z)/mp.zeta(.5+z)))/z**2
    print(f" z={z}: rhs {rhs:.8f}; zeros partial {lhs:.8f}")
print("DH")
# a(n) real DH: f=((1-i k)/2)L(chi)+((1+i k)/2)L(chibar), chi(2)=i, mod 5
kap=0.28407904384041
chi={0:0,1:1,2:1j,4:-1,3:-1j}
N=200000
A=np.zeros(N+1,float)
for n in range(1,N+1):
    c=chi[n%5];A[n]=(((1-1j*kap)/2)*c+((1+1j*kap)/2)*np.conj(c)).real
print(" a(1..6)",A[1:7].round(4))
# c(n) of -f'/f : n a(n) log... use recursion: -f'/f = sum c(n) n^-s, with f'=-sum a log n: a(n) log n = sum_{d|n} c(d) a(n/d)
M=3000;c=np.zeros(M+1)
for n in range(2,M+1):
    s=A[n]*np.log(n)
    for d in range(2,n):
        if n%d==0: s-=c[d]*A[n//d]
    c[n]=s
print(" c(6)",c[6])
b=2*np.arange(0,200000)+1.5
C0=.5*digamma(.75)+.5*np.log(5/np.pi);C1=.25*polygamma(1,.75)
lnn=np.log(np.arange(2,M+1));cn=c[2:]/np.sqrt(np.arange(2,M+1))
def Psi(t):
    t=np.atleast_1d(t);out=C0*t+C1
    b0=b[:400]
    out=out-np.exp(-np.outer(t,b0)).dot(1/b0**2)
    for i,tt in enumerate(t):
        m=lnn<tt
        out[i]-=np.sum(cn[m]*(tt-lnn[m]))
    return out
tg=np.arange(.1,12,.001);Pv=Psi(tg);print(f" Psi_DH min on [.1,12] {Pv.min():.4f} at {tg[Pv.argmin()]:.3f}")
tg=np.arange(12,25,.001);Pv=Psi(tg);print(" first negative",tg[np.argmax(Pv<0)] if (Pv<0).any() else None)
for l in [3.0,3.6,3.7,3.8,4.0]:
    n=1200;u=np.arange(1,n+1)*l/n;Pu=Psi(np.arange(0,n+1)*l/n);Pu[0]=0
    i=np.arange(1,n+1);Gm=Pu[i][:,None]+Pu[i][None,:]-Pu[np.abs(i[:,None]-i[None,:])]
    print(f" l={l}: least eig G {np.linalg.eigvalsh(Gm)[0]:.3e}")
