# DH control. Lambda(s)=(5/pi)^{s/2} Gamma((s+1)/2) f(s), f=5^{-s}[zeta(s,1/5)+k zeta(s,2/5)-k zeta(s,3/5)-zeta(s,4/5)].
import mpmath as mp
mp.mp.dps=40
k=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
def f(s): return 5**(-s)*(mp.zeta(s,mp.mpf(1)/5)+k*mp.zeta(s,mp.mpf(2)/5)-k*mp.zeta(s,mp.mpf(3)/5)-mp.zeta(s,mp.mpf(4)/5))
La=lambda s:(5/mp.pi)**(s/2)*mp.gamma((s+1)/2)*f(s)
print('FE check',mp.nstr(abs(La(mp.mpc(.3,7))/La(mp.mpc(.7,-7))),12),mp.nstr(abs(La(mp.mpc(.3,7))-La(mp.mpc(.7,7))),3))
# Taylor at 1/2
N=128;R=mp.mpf(5)
vals=[La(mp.mpf(1)/2+R*mp.expjpi(2*mp.mpf(j)/N)) for j in range(N)]
cs=[sum(vals[j]*mp.expjpi(-2*mp.mpf(j)*n/N) for j in range(N))/N/R**n for n in range(0,24)]
print('La(1/2)',mp.nstr(cs[0],10),'odd',mp.nstr(abs(cs[1]),3))
b=[mp.re(cs[2*n]) for n in range(11)]
print('b_k',[mp.nstr(a,6) for a in b])
print('Turan ratio b_k^2/(b_{k-1}b_{k+1}), compare (k+1)/k',[mp.nstr(b[n]**2/(b[n-1]*b[n+1]),5) for n in range(1,10)])
print('b1/b0 = "L"',mp.nstr(b[1]/b[0],10))
# Off line zero near 0.8085+85.699i
r=mp.findroot(La,mp.mpc(0.8085,85.699))
print('zero',mp.nstr(r,10))
for sig,t in [(.75,mp.im(r)),(.7,mp.im(r)),(mp.re(r)+.05,mp.im(r)),(1,mp.im(r)),(.75,80)]:
    z=mp.mpc(sig,t);print('u_sigma',mp.nstr(sig,4),mp.nstr(t,6),mp.nstr(mp.re(mp.diff(La,z)/La(z))/mp.pi,6))
# Cauchy identity with on line zeros only fails for DH at sigma=1: compare u_1 with P_{1/2}*mu_on near the zero using a window
