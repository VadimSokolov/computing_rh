# Independent checks for x1_uniform_spins.tex. Does not import the note's scripts.
import mpmath as mp, json, math
mp.mp.dps=30
R='./data/zeros_8000.txt'
Gs=[mp.mpf(l.strip()) for l in open(R) if l.strip()]
g1=Gs[0]; print('zeros',len(Gs),mp.nstr(g1,12),mp.nstr(Gs[-1],12))
K0=200000
mu=[1]*(K0+1); mu[0]=0
isp=[True]*(K0+1)
for p in range(2,K0+1):
    if isp[p]:
        for m in range(p,K0+1,p):
            if m>p: isp[m]=False
            mu[m]=-mu[m]
        for m in range(p*p,K0+1,p*p): mu[m]=0
print('mu check',mu[1:13])
# sum mu/k^2 tails
tot=6/mp.pi**2; s=mp.mpf(0); tail={}
Ks=[100,1000,2000,10000,20000,100000,200000]
for k in range(1,K0+1):
    if mu[k]: s+=mp.mpf(mu[k])/k**2
    if k in Ks: tail[k]=tot-s
def lsinc(x):
    # log sinh(x)/x, series for small x
    if abs(x)<mp.mpf('1e-3'):
        return x**2/6-x**4/180+x**6/2835
    return mp.log(mp.sinh(x)/x)
def perzero(g,z,K,th=mp.mpf(1)/2):
    l=mp.mpc(0)
    for k in range(1,K+1):
        if mu[k]:
            a=mp.pi/(k*g)
            l+=mu[k]*(lsinc(a*(z+th))-lsinc(a*th))
    return mp.exp(l)
def target(g,z,th=mp.mpf(1)/2): return ((th+z)**2+g**2)/(th**2+g**2)
print('== truncation at z=0.7, gamma_1')
z=mp.mpf('0.7'); prev=None; c=(mp.pi**2/6)*((z+.5)**2-.25)/g1**2
for K in Ks:
    e=mp.re(perzero(g1,z,K)-target(g1,z))
    print(K,'err',mp.nstr(e,5),'tail_mu',mp.nstr(tail[K],5),'pred',mp.nstr(c*tail[K]*target(g1,z),5),'err*K',mp.nstr(e*K,4),'err*K^1.5',mp.nstr(e*K**1.5,4))
print('== arbitrary and DH ordinates, K=20000')
dh=mp.mpc('85.699348485378','-0.308517182457')
for name,g,zz in [('gamma1 z=.2+13i',g1,mp.mpc('.2','13')),('g=pi z=.7',mp.pi,mp.mpf('.7')),('g=1 z=.3+.2i',mp.mpf(1),mp.mpc('.3','.2')),('DH offline z=.7',dh,mp.mpf('.7')),('DH offline z=.2+13i',dh,mp.mpc('.2','13'))]:
    v=perzero(g,zz,20000); t=target(g,zz)
    print(name,mp.nstr(v,12),mp.nstr(t,12),'abs err',mp.nstr(abs(v-t),3))
# theta=0 for DH
v=mp.mpf(0)
print('== hypothesis not needed: cosh family alpha_j=(2j-1)pi/2 violates no-multiple')
zz=mp.mpf('0.7'); J=1500; K=4000
lp=mp.mpf(0); lt=mp.mpf(0)
for j in range(1,J+1):
    al=(2*j-1)*mp.pi/2
    lt+=mp.log(1+zz**2/al**2)
    for k in range(1,K+1):
        if mu[k]: lp+=mu[k]*lsinc(mp.pi*zz/(k*al))
print('double product',mp.nstr(mp.exp(lp),12),'truncated cosh prod',mp.nstr(mp.exp(lt),12),'cosh',mp.nstr(mp.cosh(zz),12),'reduced sinh2z/2z / sinhz/z',mp.nstr((mp.sinh(2*zz)/(2*zz))/(mp.sinh(zz)/zz),12))
print('== c(x) double count: ordinates {1,2,6}, z=1.3')
zz=mp.mpf('1.3'); O=[1,2,6]; M=60000
def cf(m): return (mu[m] if m<=K0 else 0)+(mu[m//2] if m%2==0 else 0)+(mu[m//6] if m%6==0 else 0)
tgt=mp.mpf(1)
for o in O: tgt*=1+zz**2/o**2
lred=mp.mpf(0); ldbl=mp.mpf(0); lmu=mp.mpf(0)
for m in range(1,M+1):
    c=cf(m)
    if c: lred+=c*lsinc(mp.pi*zz/m)
for o in O:
    for k in range(1,M//o+1):
        lmu+=mu[k]*lsinc(mp.pi*zz/(k*o)) if mu[k] else 0
        c=cf(k*o)
        if c: ldbl+=c*lsinc(mp.pi*zz/(k*o))
print('target',mp.nstr(tgt,10),'mu(k) labelled',mp.nstr(mp.exp(lmu),10),'one factor per x, c(x)',mp.nstr(mp.exp(lred),10),'note: pairs with c(k gamma)',mp.nstr(mp.exp(ldbl),10))
print('c(6)=',cf(6),'visited by (1,6),(2,3),(6,1)')
print('== rational hypothesis: 2g=3g-prime family {2,3}: c at 6')
print('c(6) for N at 2,3 =',(mu[3]+mu[2]), 'while mu(k)=', mu[3],'(gamma=2,k=3) and',mu[2],'(gamma=3,k=2)')
print('== moments, own implementation')
th=mp.mpf(1)/2
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
# direct: from logs of xi via Taylor at 1 using Stieltjes-free contour derivative around s=1
def lx(s): return mp.log(xi(s))
d=[mp.diff(lx,mp.mpf('1.0'),n,h=mp.mpf('1e-3'),direction=1) if False else None for n in (1,2)]
# Cauchy integral on circle radius .25 around 1 for accuracy
def cauchy(c,n,r=mp.mpf('.25'),N=400):
    tot=mp.mpc(0)
    for j in range(N):
        w=mp.expjpi(2*mp.mpf(j)/N)
        tot+=lx(c+r*w)/w**n
    return tot/N/r**n*mp.factorial(n)
m1=mp.re(cauchy(mp.mpf(1),1)); v1=mp.re(cauchy(mp.mpf(1),2))
print('direct mean',mp.nstr(m1,12),'var',mp.nstr(v1,12))
dens=lambda x: mp.log(x/(2*mp.pi))/(2*mp.pi)
for N in (1000,8000):
    gs=Gs[:N]; T=(Gs[N-1]+Gs[N])/2 if N<len(Gs) else Gs[-1]+mp.mpf('.2')
    msum=sum(1/(th**2+g**2) for g in gs); vsum=sum(2*(g**2-th**2)/(th**2+g**2)**2 for g in gs)
    mt=mp.quad(lambda x:dens(x)/(th**2+x**2),[T,mp.inf]); vt=mp.quad(lambda x:dens(x)*2*(x**2-th**2)/(th**2+x**2)**2,[T,mp.inf])
    print(N,'mean partial',mp.nstr(msum,10),'plus tail',mp.nstr(msum+mt,10),'var partial',mp.nstr(vsum,10),'plus tail',mp.nstr(vsum+vt,10))
print('== spin side per zero moments at gamma_1, K=100000')
g=g1; mm=vv=vabs=mp.mpf(0)
for k in range(1,100001):
    if mu[k]:
        a=mp.pi/(k*g); c=a*th
        # mean and variance of a T_c, series for small c
        if c<mp.mpf('1e-3'):
            me=a*(c/3-c**3/45); va=a**2*(mp.mpf(1)/3-c**2/15)
        else:
            me=a*(mp.coth(c)-1/c); va=a**2*(1/c**2-1/mp.sinh(c)**2)
        mm+=mu[k]*me; vv+=mu[k]*va; vabs+=va
print('mean',mp.nstr(mm,10),mp.nstr(2*th/(th**2+g**2),10),'var',mp.nstr(vv,10),mp.nstr(2*(g**2-th**2)/(th**2+g**2)**2,10),'literal sign var',mp.nstr(vabs,10),'5/g^2 approx',mp.nstr(5/g**2,10))
print('== literal reading bookwide')
S2=sum(1/g**2 for g in Gs)+mp.quad(lambda x:dens(x)/x**2,[Gs[-1]+mp.mpf('.2'),mp.inf])
print('sum 1/g^2',mp.nstr(S2,10),'VarX=2S2',mp.nstr(2*S2,10),'S+ var',mp.nstr(3.5*S2,8),'S- var',mp.nstr(1.5*S2,8),'literal var 5S2',mp.nstr(5*S2,8),'k=1 sum var',mp.nstr(mp.pi**2/3*S2,8),'W2 lower',mp.nstr(mp.sqrt(mp.pi**2/3*S2)-mp.sqrt(2*S2),6))
print('sums: mu=+1',mp.nstr(21/(2*mp.pi**2),8),'mu=-1',mp.nstr(9/(2*mp.pi**2),8),'k=2 share',mp.nstr(mp.mpf(1)/4/(9/(2*mp.pi**2)),5))
print('== kappa4')
k4z=-12*sum(1/g**4 for g in Gs)
k4d=mp.re(cauchy(mp.mpf(1)/2,4))
B4=mp.bernoulli(4); ks=2**4*B4/4
k4s=mp.mpf(0)
for k in range(1,20001):
    if mu[k]: k4s+=mp.mpf(mu[k])/k**4
k4spin=ks*mp.pi**4*k4s*sum(1/g**4 for g in Gs)
print('direct',mp.nstr(k4d,10),'-12 sum g^-4',mp.nstr(k4z,10),'spin kappa4 sum',mp.nstr(k4spin,10))
print('== multiples')
print('2g1',mp.nstr(2*g1,10),'3g1',mp.nstr(3*g1,10))
print('neighbours 2g1',[mp.nstr(x,8) for x in Gs if abs(x-2*g1)<6],'3g1',[mp.nstr(x,8) for x in Gs if abs(x-3*g1)<3])
def xi1(t): return xi(mp.mpc('.5',t))
print('xi(1/2+2i g1)',mp.nstr(xi1(2*g1),6),'zeta',mp.nstr(mp.zeta(mp.mpc('.5',2*g1)),6),'xi(1/2+3i g1)',mp.nstr(xi1(3*g1),6),'zeta',mp.nstr(mp.zeta(mp.mpc('.5',3*g1)),6))
