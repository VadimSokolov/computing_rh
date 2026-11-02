# Q1: the prime formula in balls against the note's table and against the zeros.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
import mpmath as mp
ctx.prec=256
C=consts()
pps=pp_ball(250000)
# the mean sum_rho gamma^-2 = (log xi)''(1/2) = -8 + psi'(1/4)/4 + (log zeta)''(1/2)
mp.mp.dps=40
z0=mp.zeta(0.5); z1=mp.zeta(0.5,derivative=1); z2=mp.zeta(0.5,derivative=2)
mean=-8+mp.polygamma(1,0.25)/4+z2/z0-(z1/z0)**2
print('mean sum_rho gamma^-2 =',mp.nstr(mean,20),' (log zeta)\'\'(1/2)=',mp.nstr(z2/z0-(z1/z0)**2,12))
sumrho=2+mp.euler-mp.log(4*mp.pi)
print('sum_rho 1/(rho(1-rho)) =',mp.nstr(sumrho,20))
G=load_zeros(8000); mean=float(mean)
note={0.3:(0.047366,0.047367,0.047366),0.7:(0.061767,0.061734,0.065126),1.0:(0.044008,0.044004,0.194406),2.0:(0.053344,0.053340,1.799195),3.0:(0.039715,0.039711,6.167731),5.0:(0.045388,0.045406,31.598866),8.0:(0.030992,0.030960,193.203242)}
print('\nt | Psi ball (radius) | arch part ball | note: prime, zero400, arch | raw400 | note tail add | completed400 | completed8000')
def tailnote(T): return 2*(1/(2*math.pi*T))*(math.log(T/(2*math.pi))+1)
for t,(a,b,cc) in note.items():
    tb=arb(str(t))
    ar=arch_ball(tb,C); P=ar-prime_sum_ball(tb,pps)
    raw400=2*np.sum((1-np.cos(G[:400]*t))/G[:400]**2)
    comp=lambda K: mean-2*np.sum(np.cos(G[:K]*t)/G[:K]**2)
    print('%.2f %.9f (r=%.1e) %.6f | %.6f %.6f %.6f | %.6f %.6f %.7f %.7f | d(prime)=%+.1e d(arch)=%+.1e d(zero)=%+.1e'%(t,fl(P),fl(P.rad())*1.0,fl(ar),a,b,cc,raw400,raw400+tailnote(G[399]),comp(400),comp(8000),a-fl(P),cc-fl(ar),b-fl(P)))
print('tail of the raw sum at gamma_400=%.2f: %.2e ; note says about 1e-4'%(G[399],tailnote(G[399])))
for t in [12.05]:
    tb=arb(str(t)); P=arch_ball(tb,C)-prime_sum_ball(tb,pps)
    print('Psi(12.05) =',fl(P),'radius',fl(P.rad()),' note 0.03780; completed 8000:',mean-2*np.sum(np.cos(G*t)/G**2))
# min of Psi on [log 2, 12] and on [0.1,12]
ts=np.arange(0.1,12.0001,0.005)
ns,Ls=prime_powers(250000); ln=np.log(ns); u=Ls/np.sqrt(ns)
def Psi_f(t):
    return 8*(math.cosh(t/2)-1)+0.5*t*(float(digamma(0.25))-math.log(math.pi))+0.25*float(polygamma(1,0.25))-series_E(t,0.5)-np.sum(u[ln<=t]*(t-ln[ln<=t]))
vals=np.array([Psi_f(t) for t in ts])
i=int(np.argmin(vals)); print('min Psi on [0.1,12] = %.6f at t=%.3f; max=%.6f at t=%.3f; time mean=%.5f (mean sum=%.5f)'%(vals[i],ts[i],vals.max(),ts[int(np.argmax(vals))],vals[ts>1].mean(),mean))
# check of double formula against zero sum on the grid and closed form c(m)
d=0.025; M=482; cz=cov_zeta(d,M,ln,u); c0=cov_from_zeros(d,M,G)
print('max|c_prime - c_8000zeros| over m>=1: %.2e ; at m=0: %.2e (tail of the aliased zeros)'%(np.abs(cz[1:]-c0[1:]).max(),abs(cz[0]-c0[0])))
# series check: A(t) against the integral of the note, Gauss panels of length 1 to r=6000 plus the tail
xg,wg=np.polynomial.legendre.leggauss(24)
def A_int(t,R=6000):
    tot=0.0
    for a in range(0,R):
        r=a+0.5+0.5*xg
        tot+=0.5*np.sum(wg*(1-np.cos(r*t))/r**2*np.real(digamma(0.25+0.5j*r)))
    tail=(math.log(R/2)+1)/R+math.sin(R*t)*math.log(R/2)/(t*R*R)
    return (tot+tail)/math.pi
for t in [0.3,1.0,2.0,8.0]:
    ser=fl(arb(str(t))/2*C['psi14']+C['tri14']/4-E_ball(arb(str(t)),arb(1)/2))
    print('A(%.1f): Gauss integral %.10f series %.10f diff %.1e'%(t,A_int(t),ser,A_int(t)-ser))
