# BN flow: first zeros of int e^{t u^2}Phi cos(Eu) move; those of e^{tE^2}xi do not; backward e^{tH^2}Phi undefined.
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
def Phi(u):
    u=abs(u);s=0.
    for n in range(1,8):
        a=np.pi*n*n*np.exp(2*u)
        if a>700:break
        s+=(2*np.pi**2*n**4*np.exp(4.5*u)-3*np.pi*n*n*np.exp(2.5*u))*np.exp(-a)
    return s
U=np.linspace(0,3,60001)
P=np.array([Phi(u) for u in U]);dU=U[1]-U[0]
w=np.ones_like(U)*2;w[0]=1;w[-1]=1;w[1::2]=4;w=w*dU/3
def H(E,t): return 4*np.sum(w*np.exp(t*U**2)*P*np.cos(E*U))
Es=np.arange(1,40,.01)
for t in [-2,-1,-.5,0,.1,.5,1,2,4]:
    v=np.array([H(E,t) for E in Es]);i=np.where(v[:-1]*v[1:]<0)[0][:3]
    print('t',t,'first zeros',[round(float(brentq(H,Es[j],Es[j+1],args=(t,))),6) for j in i])
print('e^{tE^2}xi zeros at every t: identical to t=0 (nonvanishing factor)')
# Langevin
d={}
for N,LL in [(16000,2.),(8000,2.),(16000,2.4)]:
    from scipy.linalg import eigh_tridiagonal
    x=np.linspace(-LL,LL,N+2)[1:-1];h=x[1]-x[0]
    def lp(u):
        return np.array([0.])
    # log Phi derivatives from series with mp-free ratios: use high order finite differences on logPhi with analytic n-sum
    def logPhi(u):
        u=np.abs(u);e=np.exp(2*u);out=np.full_like(u,-np.inf,dtype=float)
        n1=(2*np.pi**2*np.exp(4.5*u)-3*np.pi*np.exp(2.5*u))
        r=np.ones_like(u)
        for n in range(2,8):
            r+=(2*np.pi**2*n**4*np.exp(4.5*u)-3*np.pi*n*n*np.exp(2.5*u))/n1*np.exp(-np.pi*(n*n-1)*e)
        return np.log(n1)-np.pi*e+np.log(r)
    e=1e-3
    l1=(logPhi(x+e)-logPhi(x-e))/(2*e); l2=(logPhi(x+e)-2*logPhi(x)+logPhi(x-e))/e**2
    V=.25*l1**2+.5*l2
    ev=eigh_tridiagonal(2/h**2+V,-np.ones(N-1)/h**2,select='i',select_range=(0,8))[0]
    print('N,L',N,LL,'gs',round(float(ev[0]),5),[round(float(a),3) for a in ev[1:]])
xx=np.linspace(0,1.5,31)
e=1e-3
print('-(logPhi)\'\' on [0,1.5]',[ (round(float(a),2),round(float(-(logPhi(np.array([a+e]))-2*logPhi(np.array([a]))+logPhi(np.array([a-e])))[0]/e**2),2)) for a in xx[:8]])
