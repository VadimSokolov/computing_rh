import mpmath as mp, numpy as np, json
mp.mp.dps=25
def dl(s):
    s=mp.mpc(s); return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)
def xi(s): return (s-1)*mp.pi**(-s/2)*mp.gamma(s/2+1)*mp.zeta(s)
def u(a,w): return mp.re(dl(mp.mpf(a)+1j*mp.mpf(w)))/mp.pi
def Phi(u0):
    u0=abs(u0); return 2*sum((2*mp.pi**2*n**4*mp.e**(4.5*u0)-3*mp.pi*n**2*mp.e**(2.5*u0))*mp.e**(-mp.pi*n**2*mp.e**(2*u0)) for n in range(1,12))
print("int Phi",mp.quad(Phi,[-3,-1,0,1,3]),"xi(1/2)",xi(mp.mpf(1)/2))
print("mgf z=.5",mp.quad(lambda x:Phi(x)*mp.e**(.5*x),[-3,-1,0,1,3]),"xi(1)",mp.mpf(1)/2)
print("mgf z=1.3",mp.quad(lambda x:Phi(x)*mp.e**(1.3*x),[-3,-1,0,1,3]),"xi(1.8)",xi(mp.mpf('1.8')))
print("E X1 quad",mp.quad(lambda x:2*x*Phi(x)*mp.e**(.5*x),[-3,-1,0,1,3]),"closed",1+mp.euler/2-mp.log(4*mp.pi)/2)
print("u1(0) limit",(1+mp.euler/2-mp.log(4*mp.pi)/2)/mp.pi)
# var at 50, mean at 50
a=mp.mpf(50)
print("mean50",dl(a),"var50",mp.diff(lambda x:mp.re(dl(x)),a), "approx",mp.log(54/(2*mp.pi))/2,1/mp.mpf(108))
# table rows
gam=[mp.im(mp.zetazero(n)) for n in range(1,302)]
g=np.array([float(x) for x in gam[:300]]); T=(float(gam[299])+float(gam[300]))/2
def zs(al,w):
    A=al-.5
    return np.sum(A/(A*A+(w-g)**2)+A/(A*A+(w+g)**2))/np.pi
def tail(al,w):
    A=al-.5
    f=lambda x: mp.log(x/(2*mp.pi))/(2*mp.pi)*(A/(A*A+(w-x)**2)+A/(A*A+(w+x)**2))/mp.pi
    return float(mp.quad(f,[T,10*T,mp.inf]))
for w in [1e-9,14.1347,30.4249,45.0]:
    print("zero sum",w,mp.nstr(u(1,w),9),round(zs(1.0,w)+tail(1.0,w),8))
LIM=200000
s=np.ones(LIM+1,bool);s[:2]=False
for p in range(2,448):
    if s[p]: s[p*p::p]=False
ns=[];lam=[]
for p in np.nonzero(s)[0]:
    q=int(p)
    while q<=LIM: ns.append(q);lam.append(np.log(p));q*=int(p)
ns=np.array(ns,float);lam=np.array(lam)
for w in [0.,14.1347,25.]:
    S=mp.mpf(2)+1j*w
    sm=mp.re(1/S+1/(S-1)-mp.log(mp.pi)/2+mp.digamma(S/2)/2)
    print("prime",w,mp.nstr(u(2,w),9),round(float((sm-np.sum(lam*ns**-2.0*np.cos(w*np.log(ns))))/mp.pi),8),"tail est 1/(pi N)",1/(np.pi*LIM))
# semigroup independently with Gauss-style exact tail: split near part numerically, far part with exact zero-free formula
def conv(w,h=0.25,alpha0=0.75,L=250):
    pts=[-L,-100,-50]+[ -x for x in [40,30,20,10]]+[0,10,20,30,40,50,100,L]
    Wg=np.arange(-L,L+1e-9,0.025)
    return Wg
# use zero-only exact test: convolution of closed form on grid vs closed; tail with mean density vs. exact digamma-mean
Wg=np.arange(-250,250.0001,0.05)
ua=np.array([float(u(.75,abs(x)+1e-12)) for x in Wg])
for w in [1e-9,14.1347,40.0]:
    h=.25
    ker=(h/np.pi)/(h*h+(w-Wg)**2)
    c=np.trapezoid(ua*ker,Wg)
    # note's tail: folds g>250 twice, with the kernel centred at +w only
    tl=2*float(mp.quad(lambda x:mp.log(x/(2*mp.pi))/(2*mp.pi)*(h/mp.pi)/(h*h+(x-w)**2),[250,mp.inf]))
    tl2=float(mp.quad(lambda x:mp.log(x/(2*mp.pi))/(2*mp.pi)*(h/mp.pi)*(1/(h*h+(x-w)**2)+1/(h*h+(x+w)**2)),[250,mp.inf]))
    # the mean of u at level 3/4 above 250 includes the digamma mean: same density. window edge term
    print("semi",w,mp.nstr(u(1,w),9),"trap no tail",round(c,7),"+note tail",round(c+tl,7),"+even tail",round(c+tl2,7),"exact gap",round(float(u(1,w))-c,7),"tl",round(tl,7),"tl2",round(tl2,7))
# min on grid, step .02 on [0,60]
for a in [.55,.75,1,1.5,2]:
    W=np.arange(0,60.0001,0.05)
    v=[float(u(a,x+1e-12)) for x in W]
    print("min",a,min(v),W[int(np.argmin(v))])
# Decomposition: mgf at alpha=2, z=.7
al=mp.mpf(2);z=mp.mpf('.7')
lhs=xi(al+z)/xi(al)*(al-1)/(al-1+z)
rhs=mp.pi**(-z/2)*mp.gamma(al/2+1+z/2)/mp.gamma(al/2+1)*mp.zeta(al+z)/mp.zeta(al)
print("decomp",lhs,rhs)
# symmetric identity: alpha=.9, v=.8, k=2 sum_{g>0} e^{-g x}
al=.9;v=.8;th=al-.5
lhs=mp.log(xi(mp.mpf(al)+v)*xi(mp.mpf(al)-v)/xi(mp.mpf(al))**2)
gs=[mp.im(mp.zetazero(n)) for n in range(1,301)]
# use termwise closed form of the integral with complex c'=g+i th: pair gives log(1+v^2/c^2)+conj
rhs=sum(2*mp.re(mp.log(1+v**2/(x+1j*th)**2)) for x in gs)
Tt=gs[-1]
# tail: integrate density of zeros of 2 Re log(1+v^2/(x+i th)^2) ~ 2 v^2/x^2
tailv=mp.quad(lambda x:mp.log(x/(2*mp.pi))/(2*mp.pi)*2*mp.re(mp.log(1+v**2/(x+1j*th)**2)),[Tt,10*Tt,mp.inf])
print("sym termwise",lhs,rhs+tailv)
# direct quadrature with k=2 sum e^{-g x}, x>0 doubled, for x>=x0 and compare factor
def k(x): return 2*sum(mp.e**(-y*x) for y in gs[:300])
I=2*mp.quad(lambda x:2*mp.sin(v*x/2)**2*mp.cos(th*x)*k(x)/x,[0.02,.1,.5,2,6])
print("sym quad (x>=.02) with k=2Sum",I,"lhs",lhs)
# level 1/2: Re dl = 0
print("Re xi'/xi at 1/2+iw",[mp.nstr(mp.re(dl(mp.mpf(1)/2+1j*w)),3) for w in [3,10,17.3]])
