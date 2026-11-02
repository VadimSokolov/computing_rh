# Identities: xi-mellin factor, Parseval, Table 1, Table 2 (Cauchy), length, Turan, class (strip), BN flow.
import mpmath as mp, sys, json
mp.mp.dps=30
xi=lambda s: mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Phi(u):
    u=mp.mpf(u); s=mp.mpf(0)
    for n in range(1,12):
        a=mp.pi*n*n*mp.e**(2*abs(u))
        if a>400: break
        s+=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)/2*abs(u))-3*mp.pi*n*n*mp.e**(mp.mpf(5)/2*abs(u)))*mp.e**(-a)
    return s
def I(E): return mp.quad(lambda u:Phi(u)*mp.cos(E*u),[0,.25,.5,1,1.5,2.5])
print('E, xi(1/2+iE), 4 int_0^inf Phi cos, 4sqrtpi*M+g, 2sqrtpi*M+g')
for E in [0,5,10,15]:
    v=4*I(E); Mg=(1/(2*mp.sqrt(mp.pi)))*2*(2*I(E))  # (1/2sqrtpi)*2*int_R e^{iEu}Phi du
    print(E,mp.nstr(mp.re(xi(mp.mpf(1)/2+1j*E)),11),mp.nstr(v,11),mp.nstr(4*mp.sqrt(mp.pi)*Mg,11),mp.nstr(2*mp.sqrt(mp.pi)*Mg,11))
# direct M+g by x integral at E=5 (no shortcut)
E=5
Md=(1/(2*mp.sqrt(mp.pi)))*2*mp.quad(lambda x:x**(-1)*mp.cos(E*mp.log(x))*Phi(mp.log(x)),[mp.e**-2.5,mp.e**-1,1,mp.e,mp.e**2.5])
print('direct M+g(5)',mp.nstr(Md,11),'xi/(2sqrtpi)',mp.nstr(mp.re(xi(.5+5j))/(2*mp.sqrt(mp.pi)),11))
# Parseval spectral measure |M g|^2
n2=2*2*mp.quad(lambda u:Phi(u)**2,[0,.25,.5,1,1.5,2.5])
sp=mp.quad(lambda E:2*(mp.re(xi(.5+1j*E))/(2*mp.sqrt(mp.pi)))**2,[0,5,10,20,40,80])
print('||g||^2=',mp.nstr(n2,10),' int_R|M+g|^2 dE=',mp.nstr(sp,10))
# Length, b_k
b=mp.taylor(lambda z:xi(mp.mpf(1)/2+z),0,22,radius=mp.mpf(1)/2) if 0 else None
mp.mp.dps=60
cs=[]
R=mp.mpf(6)
N=256
vals=[xi(mp.mpf(1)/2+R*mp.expjpi(2*mp.mpf(j)/N)) for j in range(N)]
for k in range(0,24):
    cs.append(mp.re(sum(vals[j]*mp.expjpi(-2*mp.mpf(j)*k/N) for j in range(N))/N)/R**k)
bk=[cs[2*k] for k in range(0,12)]
print('L=b1/b0',mp.nstr(bk[1]/bk[0],16),' odd coeff',mp.nstr(cs[1],3))
for k in range(0,11):
    t=None if k==0 else bk[k]**2-bk[k-1]*bk[k+1]
    r=None if k==0 else bk[k]**2/(bk[k-1]*bk[k+1])
    print(k,mp.nstr(bk[k],12),'' if t is None else mp.nstr(t,3),'' if r is None else mp.nstr(r,6), '' if k==0 else mp.nstr(mp.mpf(k+1)/k,5))
mp.mp.dps=25
print('xi1/xi(1) =1+g/2-log(4pi)/2',mp.nstr(1+mp.euler/2-mp.log(4*mp.pi)/2,10))
# Table 2: exact
zz=[mp.mpf(str(mp.im(mp.zetazero(k)))) for k in range(1,301)]
def exact(t,s):
    z=mp.mpc(s,t); return mp.re(mp.diff(lambda w:mp.log(xi(w)),z))
def sum_(t,s,N):
    a=s-.5; return sum(a/(a*a+(t-g)**2)+a/(a*a+(t+g)**2) for g in zz[:N])
def tail(t,s,G):
    a=s-.5;f=lambda g:(a/(a*a+(t-g)**2)+a/(a*a+(t+g)**2))*mp.log(g/(2*mp.pi))/(2*mp.pi)
    return mp.quad(f,[G,2*G,10*G,100*G,mp.inf])
for t,s in [(0,1),(0,.75),(10,1),(10,.75),(14.1347,1),(14.1347,.75),(30,1),(30,.75)]:
    t=mp.mpf(t);s=mp.mpf(s);e=exact(t,s);a=sum_(t,s,100);c=sum_(t,s,300);d=tail(t,s,zz[299])
    print(t,s,mp.nstr(e,8),mp.nstr(a,8),mp.nstr(c,8),mp.nstr(c+d,8),'relerr',mp.nstr(abs(c+d-e)/e,2))
print('S300',mp.nstr(sum(1/g**2 for g in zz),8))
# Exact value check: Sigma_{all rho}|rho|^-2 vs Sigma_{gamma>0}
print('Sigma_gamma>0 (1/4+g^2)^-1 300 terms',mp.nstr(sum(1/(mp.mpf(1)/4+g**2) for g in zz),8))
# Phi not entire: partial sums at u=i*0.5 (beyond pi/4=.785? Im u must exceed pi/4)
mp.mp.dps=20
def term(n,u): return (2*mp.pi**2*n**4*mp.exp(4.5*u)-3*mp.pi*n*n*mp.exp(2.5*u))*mp.exp(-mp.pi*n*n*mp.exp(2*u))
for u in [mp.mpc(0,.7),mp.mpc(0,.9)]:
    print('u',u,[mp.nstr(abs(term(n,u)),3) for n in (1,2,3,5,8)])
# decay of xi on the line vs e^{-pi E/4}: an entire even Phi of Schwartz class would force faster than exponential
for E in [50,100,200,400]:
    x=abs(mp.pi**(-mp.mpf(1)/4)*mp.gamma(mp.mpc(.25,E/2.))*(.25+E*E)/2)  # |xi|/|zeta|
    print('E',E,'log(|xi/zeta|)/E',mp.nstr(mp.log(x)/E,5),' -pi/4=',mp.nstr(-mp.pi/4,5))
