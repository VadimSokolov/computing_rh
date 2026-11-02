import mpmath as mp
mp.mp.dps=110
def xi(s):
    if s==1: return mp.mpf(1)/2
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def dlx(z): return 1/z+1/(z-1)-mp.log(mp.pi)/2+mp.digamma(z/2)/2+mp.zeta(z,derivative=1)/mp.zeta(z)
def talbot(F,t,M=220):
    # fixed Talbot (Abate-Valko)
    s=mp.mpf(0); r=2*M/(5*t)
    for k in range(M):
        if k==0:
            d=r; 
            s+=mp.mpf(1)/2*mp.exp(r*t)*F(r)
        else:
            th=k*mp.pi/M; c=mp.cot(th); d=r*th*(c+1j)
            sg=1+1j*th*(1+c*c)-1j*c
            s+=mp.re(mp.exp(t*d)*F(d)*sg)
    return r/M*s
def c_alpha(a): return -mp.log(2)+mp.log(2*mp.pi)/2-(a/2)*mp.log(mp.pi)-((a-1)/2)*mp.log(2)-mp.log(xi(a))
for a in [mp.mpf(1),mp.mpf(2)]:
    xa=xi(a); ca=c_alpha(a); print('alpha',a,'c_alpha',mp.nstr(ca,8))
    for e in ['0.004','0.003','0.001','1e-4','1e-6','1e-8','1e-12']:
        e=mp.mpf(e)
        X=mp.findroot(lambda x:dlx(a+x)/(2*x)-e,-mp.lambertw(-8*mp.pi*e,-1).real/(4*e))
        s=X*X; I=mp.log(xi(a+X)/xa)-s*e
        xe=-mp.lambertw(-8*mp.pi*e,-1).real/(4*e)
        W=xe*xe*e-xe/2+(a+3)/2*mp.log(xe)
        lead=mp.log(1/e)**2/(16*e)
        out=[mp.nstr(e,3),mp.nstr(I,10),mp.nstr(W,10),'I-W',mp.nstr(I-W,6),'I/lead',mp.nstr(I/lead,4)]
        if e>=mp.mpf('0.0029'):
            dLam=lambda q: dlx(a+mp.sqrt(q))/(2*mp.sqrt(q))
            d2=mp.diff(dLam,s)
            sp=mp.exp(-I)/(s*mp.sqrt(2*mp.pi*abs(d2)))
            mp.mp.dps=110
            P=talbot(lambda q:xa/(q*xi(a+mp.sqrt(q))),e)
            out+=['Talbot',mp.nstr(P,8),'sp/exact',mp.nstr(sp/P,6)]
        print(*out)
# median scale of Y_1: s with Lambda(s)=log 2, and mean finite?
a=mp.mpf(1); print('xi prime/xi at 1',mp.nstr(dlx(mp.mpf(1)+mp.mpf(10)**-30),8))
x=mp.findroot(lambda x:mp.log(xi(1+x)/xi(a))-mp.log(2),5); print('Lambda=log2 at x',mp.nstr(x,6),'1/s',mp.nstr(1/x**2,5))
for e in ['0.01','0.02','0.05','0.1','0.2','0.5']:
    print('P(Y1<=',e,')',mp.nstr(talbot(lambda q:xi(a)/(q*xi(1+mp.sqrt(q))),mp.mpf(e),M=120),6))
