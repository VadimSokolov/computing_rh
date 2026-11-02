import numpy as np, mpmath as mp
from scipy.special import digamma, polygamma, loggamma
g=np.load('/tmp/n8rev/zeros.npy')
th=lambda t: np.imag(loggamma(0.25+0.5j*np.asarray(t,float)))-np.asarray(t,float)/2*np.log(np.pi)
# width restricted to height 680
for name,ps in [('none',[]),('2',[2]),('3',[2,3])]:
    def T(t):
        s=th(t)
        for p in ps:
            for k in range(1,60): s=s-np.sin(k*t*np.log(p))/(k*p**(k/2))
        return s
    idx=np.arange(8,421);e=idx-T(g[7:])/np.pi
    mask=g[7:]<680
    print(name,"width m=8..420",round(np.ptp(e),4),"m=8..400 (to 680)",round(np.ptp(e[mask]),4))
    i30=np.arange(5,421);e2=i30-T(g[4:])/np.pi;m2=g[4:]<680
    print("   from gamma_5 (first >30) to 680:",round(np.ptp(e2[m2]),4))
# DH prime side
kap=0.28407904384041;chi={0:0,1:1,2:1j,4:-1,3:-1j}
M=170000
A=np.array([0.]+[(((1-1j*kap)/2)*chi[n%5]+((1+1j*kap)/2)*np.conj(chi[n%5])).real for n in range(1,M+1)])
rhs=A*np.log(np.maximum(np.arange(M+1),1));c=np.zeros(M+1)
for d in range(2,M+1):
    c[d]=rhs[d]
    if 2*d<=M: rhs[2*d::d]-=c[d]*A[2:M//d+1]
lnn=np.log(np.arange(2,M+1));cn=c[2:]/np.sqrt(np.arange(2,M+1))
b0=2*np.arange(400)+1.5
C0=.5*digamma(.75)+.5*np.log(5/np.pi);C1=.25*polygamma(1,.75)
def Psi(t):
    t=np.atleast_1d(t);o=C0*t+C1-np.exp(-np.outer(t,b0)).dot(1/b0**2)
    o=np.array(o)
    for i,tt in enumerate(t):
        j=np.searchsorted(lnn,tt);o[i]-=np.sum(cn[:j]*(tt-lnn[:j]))
    return o
tg=np.arange(.1,12,.002);P=Psi(tg);print("Psi_DH min [.1,12]",round(P.min(),4),"at",round(tg[P.argmin()],3))
# zero side beyond 12
mp.mp.dps=15
chiL=[0,1,1j,-1j,-1]
chiB=[0,1,-1j,1j,-1]
def f(s): return (1-1j*kap)/2*mp.dirichlet(s,chiL)+(1+1j*kap)/2*mp.dirichlet(s,chiB)
def Z(t):
    s=mp.mpc(.5,t); gf=(5/mp.pi)**(s/2)*mp.gamma((s+1)/2)
    return f(s)*gf/abs(gf)
ts=np.arange(.5,200,.03);v=[]
z0=Z(1.0);
vals=np.array([complex(Z(t)) for t in ts])
# Lambda(1/2+it) is real up to a constant phase: rotate
ph=np.angle(vals[np.argmax(np.abs(vals))]);r=(vals*np.exp(-1j*ph)).real
print(" imag check",np.max(np.abs((vals*np.exp(-1j*ph)).imag)))
on=[]
for i in range(len(ts)-1):
    if r[i]*r[i+1]<0: on.append(ts[i]-r[i]*(ts[i+1]-ts[i])/(r[i+1]-r[i]))
on=np.array(on);print(" on line below 200:",len(on),on[:3].round(4))
off=[0.80851718+85.69934849j,0.65083008+114.16334273j,0.57435605+166.47930591j,0.72425769+176.70246124j]
for s0 in off[:2]:
    print(" |f| at listed off line zero",float(abs(f(mp.mpc(s0)))))
gam=[-1j*(s-.5) for s in off]
def Pz(t):
    s=np.sum(2*(1-np.cos(on*t))/on**2)
    for gg in gam:
        # four zeros: g, -g, conj g, -conj g -> 2 Re each pair
        s+=4*((1-np.cos(gg*t))/gg**2).real
    # tail above 200: mean
    Td=200.;tail=2*(np.log(5*Td/(2*np.pi))+1)/(2*np.pi*Td)*2/2*2/2
    return s+tail,s
for t in [.5,6.,12.]:
    print(f" t={t}: prime side {Psi(t)[0]:.5f} zero side (with mean tail) {Pz(t)[0]:.5f}")
tt=np.arange(12,30,.002);vz=np.array([Pz(x)[0] for x in tt]);print(" first negative (zero side) ",tt[np.argmax(vz<0)] if (vz<0).any() else None)
# zeta kernels from book formula (c), sign of screw kernel
LAM={}
def vm(n):
    for p in range(2,int(n**.5)+1):
        if n%p==0:
            m=n
            while m%p==0:m//=p
            return np.log(p) if m==1 else 0.
    return np.log(n)
ns=np.arange(2,3000);LN=np.array([vm(int(n)) for n in ns])/np.sqrt(ns);lg=np.log(ns)
bz=2*np.arange(400)+.5
def Pzeta(t):
    t=np.atleast_1d(t);o=8*(np.cosh(t/2)-1)+t/2*(digamma(.25)-np.log(np.pi))+.25*polygamma(1,.25)-np.exp(-np.outer(t,bz)).dot(1/bz**2)
    o=np.array(o)
    for i,x in enumerate(t):o[i]-=np.sum((LN*(x-lg))[lg<x])
    return o
print(" Pzeta(log2)*2",2*Pzeta(np.log(2))[0],"book D(log2)=0.18338*log2=",.18338*np.log(2))
n=400;l=3.;u=np.arange(0,n+1)*l/n;PP=Pzeta(u);PP[0]=0
i=np.arange(1,n+1);Gm=PP[i][:,None]+PP[i][None,:]-PP[np.abs(i[:,None]-i[None,:])]
# Krein Langer kernel for g=-Psi equals G, for g=+Psi equals -G
print(" zeta window 3, n=400: least eig G (g=-Psi kernel)",np.linalg.eigvalsh(Gm)[0]," least eig -G (g=+Psi kernel)",np.linalg.eigvalsh(-Gm)[0])
# Suzuki form as recalled: g(t)=sum (e^{i gamma t}-1)/gamma^2
t0=2.3;sz=2*np.sum((np.cos(g*t0)-1)/g**2);pz=2*np.sum((1-np.cos(g*t0))/g**2);print(" partial sums, sum(e^{i g t}-1)/g^2",sz," Psi partial",pz,"sum",sz+pz)
# rank one nonselfadjoint statement
al=np.array([0.,1.]);be=np.array([.5,2.]);u=np.array([-np.prod(a-be)/np.prod(np.delete(al,i)-a if False else a-np.delete(al,i)) for i,a in enumerate(al)])
Mx=np.diag(al)+np.outer(u,np.ones(2));print(" 2x2 interlacing example u",u,"eig",np.linalg.eigvals(Mx),"symmetric?",np.allclose(Mx,Mx.T))
al=np.array([0.,1.]);be=np.array([3.,5.]);u=np.array([-np.prod(a-be)/np.prod(a-np.delete(al,i)) for i,a in enumerate(al)]);Mx=np.diag(al)+np.outer(u,np.ones(2));print(" non interlacing u",u,"eig",np.linalg.eigvals(Mx))
