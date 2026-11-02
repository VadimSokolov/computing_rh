"""Independent calculations for the four-note audit. No review-folder imports.

Run: timeout 600 python3 -I independent_checks.py base|zeros|synth|prime|dh
All output is ordinary multiprecision or floating point, not ball arithmetic.
"""
import json
import math
import sys
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent
mp.mp.dps = 55

def save(name, data):
    (ROOT / (name + '.json')).write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(data, indent=2), flush=True)

def st(x):
    return mp.nstr(x, 32)

def xi(s):
    if s == 0 or s == 1:
        return mp.mpf('0.5')
    return s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)/2

def phi(u):
    u = abs(u)
    t = mp.pi*mp.exp(4*u)
    return mp.exp(u)*mp.fsum((2*(t*n*n)**2-3*t*n*n)*mp.exp(-t*n*n)
                             for n in range(1, 12))

def integ(f):
    return 2*mp.quad(f, [0, mp.mpf('.125'), mp.mpf('.25'), mp.mpf('.5'), 1, 2])

def base():
    mass = integ(phi)
    moments = [mp.mpf(1)]
    for k in range(1, 11):
        moments.append(integ(lambda u: u**k*phi(u))/mass if k % 2 == 0 else mp.mpf(0))
    cumulants = [mp.mpf(0)]
    for n in range(1, 11):
        cumulants.append(moments[n]-mp.fsum(mp.binomial(n-1,j-1)*cumulants[j]*moments[n-j]
                                          for j in range(1,n)))
    v = moments[2]
    tv_slope = 2*mp.quad(lambda u: (v-u*u)*phi(u)/mass, [0,mp.sqrt(v)])
    rows = []
    for t in map(mp.mpf, ['.01','.05','.1','.2','.22','.5']):
        norm = integ(lambda u: mp.exp(t*u*u)*phi(u))
        cross = mp.sqrt(mp.log(norm/mass)/t)
        tv = 2*mp.quad(lambda u: phi(u)/mass-mp.exp(t*u*u)*phi(u)/norm,[0,cross])
        var = integ(lambda u: u*u*mp.exp(t*u*u)*phi(u))/norm
        rows.append(dict(t=st(t),tv=st(tv),kolmogorov=st(tv/2),variance=st(var)))
    kfun = lambda z: mp.log(xi(mp.mpf('.5')+z/2)/xi(mp.mpf('.5')))
    save('base_numbers',dict(mpmath=mp.__version__,dps=mp.mp.dps,H0=st(mass),xi_half=st(xi(mp.mpf('.5'))),
        moments={str(k):st(moments[k]) for k in range(2,11,2)},
        cumulants={str(k):st(cumulants[k]) for k in range(2,11,2)},tv_slope=st(tv_slope),tilts=rows,
        K={str(z):st(kfun(mp.mpf(z))) for z in ['.5','3','4','6','10','20']},
        K3={str(z):st(mp.diff(kfun,mp.mpf(z),3)) for z in ['.05','.25','1','12','21','80']}))

def zeros():
    nmax = int(sys.argv[2]) if len(sys.argv)>2 else 400
    vals=[]
    for j in range(1,nmax+1):
        vals.append(st(mp.im(mp.zetazero(j))))
        if j%50==0:
            print('independently generated zeros:',j,flush=True)
    (ROOT/'generated_zeros.json').write_text(json.dumps(vals)+'\n')
    save('zero_checks',dict(dps=mp.mp.dps,n=nmax,last=vals[-1],
        Xi_at_multiples={str(k):st(xi(mp.mpf('.5')+1j*k*mp.mpf(vals[0]))) for k in [2,3]},
        finite_cumulants={str(2*k):st((-1)**(k+1)*mp.factorial(2*k)/k*
                         mp.fsum((2*mp.mpf(g))**(-2*k) for g in vals)) for k in range(1,6)}))

def readzeros():
    return list(map(mp.mpf,json.loads((ROOT/'generated_zeros.json').read_text())))

def synth():
    zs=readzeros()
    cases=[(2,'.000619'),(8,'.000688'),(30,'.000586'),(30,'.000588'),(60,'.000491'),
           (120,'.000376'),(160,'.000209'),(200,'.000117'),(300,'.000045')]
    rows=[]
    for n,b in cases:
        if n>len(zs):
            continue
        mp.mp.dps=80+3*n
        beta=mp.mpf(b)
        e=[mp.mpf(1)]
        for g in zs[:n]:
            a=1/(4*g*g)
            e.append(mp.mpf(0))
            for k in range(len(e)-1,0,-1):
                e[k]+=a*e[k-1]
        coeff=[]
        for j in range(n+1):
            coeff.append(mp.fsum(e[k]*(2*beta)**(-k)*mp.factorial(2*k)*(-1)**(k-j)/
                          (2**(k-j)*mp.factorial(k-j)*mp.factorial(2*j)) for k in range(j,n+1)))
        def density(u):
            y2=u*u/(2*beta)
            val=mp.polyval(list(reversed(coeff)),y2)
            return mp.exp(-y2/2)*val/mp.sqrt(4*mp.pi*beta)
        grid=[mp.mpf(j)/1000 for j in range(0,1001)]
        samples=[density(u) for u in grid]
        idx=min(range(len(samples)),key=lambda j:samples[j])
        rows.append(dict(N=n,beta=b,dps=mp.mp.dps,grid_min=st(samples[idx]),u=st(grid[idx]),f0=st(samples[0])))
        print('synthesis',rows[-1],flush=True)
    save('synthesis_checks',rows)

def necessary():
    import numpy as np
    zs=readzeros()
    rows=[]
    for n,b in [(2,.000619),(8,.000688),(30,.000586),(60,.000491),
                (120,.000376),(160,.000209),(200,.000117),(300,.000045)]:
        gg=np.array([float(g) for g in zs[:n]])
        xx=np.linspace(2*float(gg[-1])+1,12*float(gg[-1]),20001)
        vv=np.log(np.abs(1-(xx[:,None]/(2*gg))**2)).sum(axis=1)/xx**2
        j=int(np.argmax(vv)); x=mp.mpf(float(xx[j]))
        logprod=mp.fsum(mp.log(abs(1-x*x/(4*g*g))) for g in zs[:n])
        rows.append(dict(N=n,claimed_beta=b,x=st(x),necessary_beta=st(logprod/(x*x)),
                         log_abs_characteristic=st(logprod-mp.mpf(str(b))*x*x)))
    save('necessary_synthesis_bounds',rows)

def extra():
    import numpy as np
    from scipy.optimize import minimize_scalar
    from scipy.signal import fftconvolve
    zs=readzeros()
    rows=[]
    for n,guess in [(2,.000619),(4,.000667),(6,.000685),(8,.000688),(10,.000676),
                    (12,.000670),(16,.000650),(20,.000629),(24,.000610),(30,.000588),
                    (40,.000550),(60,.000491),(80,.000447),(120,.000376)]:
        mp.mp.dps=100
        e=[mp.mpf(1)]
        for g in zs[:n]:
            e.append(mp.mpf(0))
            for k in range(len(e)-1,0,-1):e[k]+=e[k-1]/(4*g*g)
        def p(u,beta):
            y=u/mp.sqrt(2*beta); h0=mp.mpf(1);h1=y
            ans=e[0]
            for k in range(2,2*n+1):
                h2=y*h1-(k-1)*h0
                if k%2==0:ans+=e[k//2]*(2*beta)**(-k//2)*h2
                h0,h1=h1,h2
            return ans
        ugrid=[mp.mpf(j)/2000 for j in range(1,401)]
        b=mp.mpf(guess)
        minj=min(range(len(ugrid)),key=lambda j:p(ugrid[j],b))
        try:
            ustar,bstar=mp.findroot(lambda u,bb:(p(u,bb),mp.diff(lambda v:p(v,bb),u)),
                                   (ugrid[minj],b),tol=mp.mpf('1e-60'),maxsteps=80)
            rows.append(dict(N=n,beta=st(bstar),touch_u=st(ustar)))
        except Exception as exc:rows.append(dict(N=n,error=type(exc).__name__))
    # Independent convolution on a fine real-space grid. Signed Hermite sums
    # are evaluated at high precision before conversion to float.
    mp.mp.dps=100
    dx=.0005; xx=np.arange(-1,1+dx/2,dx);beta=mp.mpf('.001')
    pp=np.array([float(phi(mp.mpf(abs(u)))) for u in xx]); pp/=pp.sum()*dx
    gauss=np.exp(-xx*xx/(4*float(beta)))/math.sqrt(4*math.pi*float(beta))
    target=fftconvolve(pp,gauss,mode='same')*dx
    tvs=[]
    for n in [2,10,40]:
        e=[mp.mpf(1)]
        for g in zs[:n]:
            e.append(mp.mpf(0))
            for k in range(len(e)-1,0,-1):e[k]+=e[k-1]/(4*g*g)
        cs=[mp.fsum(e[k]*(2*beta)**(-k)*mp.factorial(2*k)*(-1)**(k-j)/
                  (2**(k-j)*mp.factorial(k-j)*mp.factorial(2*j)) for k in range(j,n+1))
                  for j in range(n+1)]
        def den(u):
            y2=mp.mpf(u)**2/(2*beta)
            return mp.polyval(list(reversed(cs)),y2)*mp.exp(-y2/2)/mp.sqrt(4*mp.pi*beta)
        vals=np.array([float(den(u)) for u in xx])
        tvs.append(dict(N=n,tv=float(np.sum(abs(vals-target))*dx/2),mass=float(vals.sum()*dx),
                        grid_min=float(vals.min())))
    # The finite arithmetic-progression check is numerical, not an enclosure.
    g1=zs[0];eligible=[(j,g) for j,g in enumerate(zs,1) if 2*g1<g<=41*g1]
    near=[]
    for j,g in eligible:
        t=2*g-g1
        # Hardy Z at every transformed point is more costly than nearest-zero
        # comparison; compute it only for the nearest candidate in our list.
        nearest=min(zs,key=lambda a:abs(a-t))
        near.append((abs(nearest-t),j,t,nearest))
    best=min(near)
    dist=abs(mp.siegelz(best[2])/mp.diff(mp.siegelz,best[2]))
    save('extra_numbers',dict(threshold_contacts=rows,tv_synthesis=tvs,progression=dict(
        candidates=len(eligible),index=best[1],target=st(best[2]),nearest=st(best[3]),
        distance=st(best[0]),newton_distance=st(dist))))

def highcontacts():
    zs=readzeros();rows=[]
    for n,b0,x0 in [(120,'.000386','678.574827'),(160,'.000344','833.13727'),
                    (200,'.000312','980.01618'),(300,'.000258','1324.14724')]:
        mp.mp.dps=80+3*n
        e=[mp.mpf(1)]
        for g in zs[:n]:
            e.append(mp.mpf(0))
            for k in range(len(e)-1,0,-1):e[k]+=e[k-1]/(4*g*g)
        def equations(u,b):
            sigma=mp.sqrt(2*b);y=u/sigma;hs=[mp.mpf(1),y]
            for k in range(2,2*n+1):hs.append(y*hs[-1]-(k-1)*hs[-2])
            p=mp.fsum(e[k]*(2*b)**(-k)*hs[2*k] for k in range(n+1))
            dp=mp.fsum(e[k]*(2*b)**(-k)*2*k*hs[2*k-1]/sigma for k in range(1,n+1))
            return p,dp
        try:
            u,b=mp.findroot(equations,(mp.pi/mp.mpf(x0),mp.mpf(b0)),maxsteps=50,tol=mp.mpf('1e-30'))
            rows.append(dict(N=n,beta=st(b),touch_u=st(u),global_positivity='not established'))
        except Exception as exc:rows.append(dict(N=n,error=type(exc).__name__))
        print(rows[-1],flush=True)
    save('high_contacts',rows)

def dh():
    kap=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
    def d(s):
        return mp.fsum(a*mp.zeta(s,mp.mpf(r)/5) for r,a in [(1,1),(2,kap),(3,-kap),(4,-1)])/5**s
    def completed(s):
        return (5/mp.pi)**(s/2)*mp.gamma((s+1)/2)*d(s)
    root=mp.findroot(d,mp.mpc('.8085','85.6993'))
    def density(u):
        u=abs(u); t=mp.exp(4*u); coeff=[0,1,kap,-kap,-1]
        return mp.exp(3*u)*mp.fsum(coeff[n%5]*n*mp.exp(-mp.pi*n*n*t/5) for n in range(1,35))
    norm=integ(density)
    ms=[mp.mpf(1)]+[integ(lambda u:u**k*density(u))/norm if k%2==0 else mp.mpf(0) for k in range(1,11)]
    ks=[mp.mpf(0)]
    for n in range(1,11):
        ks.append(ms[n]-mp.fsum(mp.binomial(n-1,j-1)*ks[j]*ms[n-j] for j in range(1,n)))
    # Avoid abs at zero for differentiation; moments suffice for cumulants.
    U=lambda u:-mp.log(density(u))
    ugrid=[mp.mpf(j)/50 for j in range(1,101)]
    u3=[mp.diff(U,u,3) for u in ugrid]
    kfun=lambda z:mp.log(completed(mp.mpf('.5')+z/2)/completed(mp.mpf('.5')))
    checks=[(z,st(integ(lambda u:mp.cosh(mp.mpf(z)*u)*density(u))/norm),
              st(mp.exp(kfun(mp.mpf(z))))) for z in ['.5','2','5']]
    save('dh_control',dict(kappa=st(kap),root=st(root),residual=st(abs(d(root))),
         functional_equation_residual=st(abs(completed(mp.mpc('.7','3'))-completed(mp.mpc('.3','-3')))),
         moments={str(k):st(ms[k]) for k in range(2,11,2)},
         cumulants={str(k):st(ks[k]) for k in range(2,11,2)},
         U3_grid_min=st(min(u3)),U3_grid_min_u=st(ugrid[u3.index(min(u3))]),
         K3={z:st(mp.diff(kfun,mp.mpf(z),3)) for z in ['1','12','21']},mgf_checks=checks))

def prime():
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import digamma
    zs=readzeros()
    X=5_000_000
    sieve=np.ones(X+1,dtype=bool); sieve[:2]=False
    for p in range(2,math.isqrt(X)+1):
        if sieve[p]: sieve[p*p::p]=False
    primes=np.flatnonzero(sieve)
    ns=[]; weights=[]
    for p in primes:
        q=int(p)
        while q<=X:
            ns.append(q);weights.append(math.log(int(p))/math.sqrt(q));q*=int(p)
    ns=np.asarray(ns);weights=np.asarray(weights);logs=np.log(ns)
    heat=[]
    for b in [.001,.002,.005,.01,.02,.05]:
        A=2*math.exp(b/4)-math.log(math.pi)/math.sqrt(4*math.pi*b)+quad(
            lambda r:math.exp(-b*r*r)*digamma(.25+.5j*r).real,0,400,
            epsabs=2e-13,epsrel=2e-13,limit=400)[0]/math.pi
        P=-2*np.sum(weights*np.exp(-logs**2/(4*b)))/math.sqrt(4*math.pi*b)
        Z=2*mp.fsum(mp.exp(-mp.mpf(str(b))*g*g) for g in zs)
        heat.append(dict(b=b,all_sign_zero_sum=st(Z),arch=A,prime=float(P),explicit=A+P))
    # Shifted Gaussian test, half the all-zero formula is the positive spectrum.
    spectral=[]
    for t in [5,14.15,21,30.4,43.35,59.35,60]:
        pole=2*math.exp(-2*t*t+.5)*math.cos(2*t)
        g0=1/math.sqrt(2*math.pi)
        A=pole-g0*math.log(math.pi)+quad(lambda r:(math.exp(-2*(r-t)**2)+
            math.exp(-2*(r+t)**2))*digamma(.25+.5j*r).real,0,90,
            epsabs=1e-12,epsrel=1e-12,points=[t],limit=500)[0]/math.pi
        P=-2*np.sum(weights*np.exp(-logs**2/8)*np.cos(t*logs))/math.sqrt(2*math.pi)
        target=mp.fsum(mp.exp(-2*(mp.mpf(t)-g)**2)+mp.exp(-2*(mp.mpf(t)+g)**2) for g in zs)
        spectral.append(dict(t=t,positive_spectrum=st(target),half_explicit=(A+P)/2,
                             residual=(A+P)/2-float(target)))
    # Cut error at z=3 is obtained directly via E1 for each known zero.
    errors=[]; gamma1=zs[0]; c=mp.mpf(9)/4
    for cut in [2,3,5,10,20,50,100]:
        B=mp.log(cut)/(2*gamma1)
        tail=mp.fsum(mp.e1(B*g*g)-mp.e1(B*(g*g+c)) for g in zs)
        # Stable positive formula for the b integral of each missing prime.
        def one(L):
            aa=mp.sqrt(c); root=mp.sqrt(B)
            J0=mp.erfc(L/(2*root))/L
            Jc=(mp.exp(-aa*L)*mp.erfc(L/(2*root)-aa*root)+
                mp.exp(aa*L)*mp.erfc(L/(2*root)+aa*root))/(2*L)
            return J0-Jc
        missing=mp.mpf(0)
        # Beyond exp(16*sqrt(B)) Gaussian suppression is far below displayed errors.
        for n,w in zip(ns,weights):
            if cut<n<max(1000,math.exp(16*math.sqrt(float(B)))):
                missing+=mp.mpf(float(w))*one(mp.log(int(n)))
        errors.append(dict(X=cut,B=st(B),cut_tail=st(tail),missing_primes=st(missing),
                           net_error=st(missing-tail),reference_rate=st(mp.mpf(cut)**(-gamma1/2))))
    save('prime_checks',dict(heat=heat,spectral=spectral,truncation=errors))

if __name__=='__main__':
    globals()[sys.argv[1]]()
