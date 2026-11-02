import mpmath as mp, numpy as np
from scipy.integrate import solve_ivp
mp.mp.dps=40
print("A. DLMF 13.14.33 against inward ODE integration from z0 (asymptotic series 13.19.3)")
def asym(k,nu,z0,terms=25):
    kap=-k; m=1j*nu; s=mp.mpf(0); term=mp.mpf(1)
    for j in range(terms):
        s+=term; term=term*(mp.mpf(1)/2+m-kap+j)*(mp.mpf(1)/2-m-kap+j)/(-(z0))/(j+1)
    return s
def ode(k,nu,z=4*np.pi,z0=400.0):
    # W''=(1/4+k/z-(1/4+nu^2)/z^2)W ; scale W=e^{-z/2} z^{-k} f to avoid underflow
    h=mp.mpf('1e-6')
    F=lambda zz: mp.exp(-zz/2)*zz**(-k)*asym(k,nu,zz)
    f0=float(mp.re(asym(k,nu,mp.mpf(z0))))
    dW=mp.diff(F,mp.mpf(z0))
    w0=F(mp.mpf(z0))
    # work with g=W*e^{z/2}: g''-g'+(1/4 - q)g... q=1/4+k/z-(1/4+nu^2)/z^2 ; g''=g'+(q-1/4)g
    g0=float(mp.re(w0*mp.exp(mp.mpf(z0)/2))); gp0=float(mp.re((dW+w0/2)*mp.exp(mp.mpf(z0)/2)))
    q=lambda z:k/z-(0.25+nu**2)/z**2
    # integrate inward with u=-z: stable enough
    sol=solve_ivp(lambda z,Y:[Y[1],Y[1]+q(z)*Y[0]],[z0,z],[g0,gp0],method='DOP853',rtol=1e-12,atol=1e-30)
    return sol.y[0,-1]*np.exp(-z/2)
def conn(k,nu,z):
    z=mp.mpf(z);n=mp.mpf(nu)
    a=mp.mpf(1)/2+k-1j*n;b=1-2j*n
    br=mp.gamma(2j*n)/mp.gamma(mp.mpf(1)/2+k+1j*n)*z**(mp.mpf(1)/2-1j*n)*mp.exp(-z/2)*mp.hyp1f1(a,b,z)
    return 2*mp.re(br)
for k,nu in [(0,3.3),(0,12.),(1,5.),(1,12.),(-1,7.),(3,9.),(-3,9.),(2.5,6.)]:
    c=conn(k,nu,4*mp.pi);o=ode(k,nu);w=mp.whitw(-k,1j*nu,4*mp.pi)
    print(f" k={k} nu={nu}: formula {mp.nstr(c,10)} ODE {o:.10g} mp.whitw {mp.nstr(mp.re(w),10)}")
print("B. levels: zeros in nu of the formula, residual of the rule, times nu")
mp.mp.dps=30
theta=lambda t: mp.im(mp.loggamma(mp.mpf(1)/4+0.5j*t))-t/2*mp.log(mp.pi)
for k in [-3,-1,0,1,3]:
    out=[]
    for nu0 in [50,100,200,400]:
        Phi=lambda n: conn(k,n,4*mp.pi)
        # locate sign change starting at nu0
        a=mp.mpf(nu0);d=mp.pi/mp.log(a/mp.pi)/6
        fa=Phi(a)
        while True:
            b=a+d;fb=Phi(b)
            if mp.sign(fa)!=mp.sign(fb):break
            a,fa=b,fb
        r=mp.findroot(Phi,(a,b),solver='anderson',tol=1e-15,verify=False)
        v=(theta(2*r)+3*mp.pi/8-k*mp.pi/2)/mp.pi
        res=(v-mp.nint(v))*mp.pi
        out.append(f"nu={mp.nstr(r,8)} res*nu={mp.nstr(res*r,5)}")
    print(f" k={k}: "+"; ".join(out))
print("C. M limit z=4pi k=1:")
for n in [50,200,800]:
    a=mp.mpf(3)/2-1j*n;b=1-2j*n
    v=mp.hyp1f1(a,b,4*mp.pi)/mp.exp(2*mp.pi)
    print(f" nu={n}: |M e^-z/2 -1|*nu={mp.nstr(abs(v-1)*n,5)} arg*nu={mp.nstr(mp.arg(v)*n,5)}")
print("D. independent eigenvalue check by ODE in x, k=1 and 0, count against theta")
def levels_shoot(k,numax,x0=np.log(2*np.pi)):
    # Pruefer-free: integrate from large x to x0 at each nu grid; sign of y(x0)
    def y0(nu):
        xm=np.log(max(nu,6)*3+20)
        s=solve_ivp(lambda x,Y:[Y[1],(np.exp(2*x)+2*k*np.exp(x)-nu**2)*Y[0]],[xm,x0],[1e-200,-1e-200],method='DOP853',rtol=1e-10,atol=1e-300)
        return s.y[0,-1]
    g=np.arange(0.5,numax,0.02);v=[y0(n) for n in g]
    z=[]
    for i in range(len(g)-1):
        if v[i]*v[i+1]<0: z.append(g[i]-v[i]*(g[i+1]-g[i])/(v[i+1]-v[i]))
    return np.array(z)
th=lambda t: float(theta(mp.mpf(t)))
for k in [0,1,-1]:
    z=levels_shoot(k,45)
    print(f" k={k}: first levels {np.round(z[:5],3)} count {len(z)}; residual ({'theta(2nu)+3pi/8-k pi/2'}) mod pi at last 3:",[round(((th(2*r)+3*np.pi/8-k*np.pi/2+np.pi/2)%np.pi-np.pi/2),4) for r in z[-3:]], "N-theta(2numax)/pi=",round(len(z)-th(90)/np.pi,3))
