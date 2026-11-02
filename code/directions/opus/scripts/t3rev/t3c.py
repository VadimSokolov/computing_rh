# Run: timeout 590 python3 -I /tmp/t3rev/t3c.py
import mpmath as mp
mp.mp.dps=20
KAPPA=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
A=[mp.mpf(1),KAPPA,-KAPPA,mp.mpf(-1)]
def f_dh(s): return mp.mpf(5)**(-s)*sum(a*mp.zeta(s,mp.mpf(r+1)/5) for r,a in enumerate(A))
def xif(s): return (mp.mpf(5)/mp.pi)**(s/2)*mp.gamma((s+1)/2)*f_dh(s)
dl=lambda s: mp.diff(lambda z: mp.log(xif(z)),s)
r=mp.findroot(f_dh,mp.mpc('0.8085','85.6993'),tol=1e-12,verify=False); print('DH zero',mp.nstr(r,10),abs(f_dh(r)))
be,ga=r.real,r.imag
for al in (0.6,0.75,0.8,0.82,0.9,1.0000001,1.2):
    u=mp.re(dl(mp.mpc(al,ga)))/mp.pi
    print('alpha',al,'u_alpha(gamma)',mp.nstr(u,6),' b_alpha',mp.nstr(mp.re(dl(mp.mpf(al))),6))
th=mp.mpf('0.75')-be
print('block at alpha=.75 small w:',mp.nstr((th**2+ga**2)/((th+mp.mpf('1e-3'))**2+ga**2),15),'> 1 means not a Laplace transform')
print('ordinates test: 2g1=10.1883, 3g1=15.2825 vs 8.94,12.13,14.40,17.13')
mp.mp.dps=25
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
N=2000
g=[mp.im(mp.zetazero(n)) for n in range(1,N+2)]
T=(g[N-1]+g[N])/2; dens=lambda x: mp.log(x/(2*mp.pi))/(2*mp.pi)
Nt=lambda T: T/(2*mp.pi)*mp.log(T/(2*mp.pi*mp.e))+mp.mpf(7)/8
for name,fn in (('k2',lambda x: mp.pi**2/(3*x**2)),('k4',lambda x:-2*mp.pi**4/(15*x**4))):
    s=sum(fn(x) for x in g[:N]); tail=mp.quad(lambda x: dens(x)*fn(x),[T,mp.inf]); bd=fn(T)*(Nt(T)-N)
    print(name,'note style',mp.nstr(s+tail,12),'with boundary term',mp.nstr(s+tail+bd,12))
print('z2k2',mp.nstr(mp.zeta(2)*mp.diff(lambda s: mp.log(xi(s)),mp.mpf(.5),2),12))
