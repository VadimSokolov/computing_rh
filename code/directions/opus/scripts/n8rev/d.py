import mpmath as mp, numpy as np
from scipy.integrate import solve_ivp
mp.mp.dps=25
theta=lambda t: mp.im(mp.loggamma(mp.mpf(1)/4+0.5j*t))-t/2*mp.log(mp.pi)
def Phi(k,nu,z=4*mp.pi):
    n=mp.mpf(nu);a=mp.mpf(1)/2+k-1j*n;b=1-2j*n
    br=mp.gamma(2j*n)/mp.gamma(mp.mpf(1)/2+k+1j*n)*z**(mp.mpf(1)/2-1j*n)*mp.exp(-z/2)*mp.hyp1f1(a,b,z)
    return mp.arg(br)   # phase of the bracket, principal
def F(k,nu):
    return mp.cos(Phi(k,nu))
print("rule residual with secant-refined roots of cos(phase of bracket)")
for k in [-3,-1,0,1,3]:
    o=[]
    for nu0 in [100,200,400,800]:
        a=mp.mpf(nu0);d=mp.pi/mp.log(a/mp.pi)/5;fa=F(k,a)
        while True:
            b=a+d;fb=F(k,b)
            if fa*fb<0:break
            a,fa=b,fb
        for _ in range(36):
            c=(a+b)/2;fc=F(k,c)
            if fa*fc<=0:b=c
            else:a,fa=c,fc
        r=(a+b)/2
        v=(theta(2*r)+3*mp.pi/8-k*mp.pi/2)/mp.pi
        res=(v-mp.nint(v))*mp.pi
        o.append(mp.nstr(res*r,6))
    print(f" k={k}: res*nu at nu near 100,200,400,800:",o)
print("independent shooting in x")
th=lambda t: float(theta(mp.mpf(t)))
def levels(k,numax,x0=np.log(2*np.pi),step=0.05):
    def y0(nu):
        xm=np.log(max(nu,6)*3+20)
        s=solve_ivp(lambda x,Y:[Y[1],(np.exp(2*x)+2*k*np.exp(x)-nu**2)*Y[0]],[xm,x0],[1e-30,-1e-30],method='DOP853',rtol=1e-10,atol=1e-40)
        return np.sign(s.y[0,-1]),s.y[0,-1]
    g=np.arange(0.5,numax,step);v=[y0(n) for n in g];z=[]
    for i in range(len(g)-1):
        if v[i][0]!=v[i+1][0]:
            a,b=g[i],g[i+1]
            for _ in range(14):
                c=(a+b)/2
                if y0(c)[0]==y0(a)[0]:a=c
                else:b=c
            z.append((a+b)/2)
    return np.array(z)
for k in [0,1,-3]:
    z=levels(k,30)
    rr=[((th(2*r)+3*np.pi/8-k*np.pi/2+np.pi/2)%np.pi-np.pi/2) for r in z]
    cnt=[len(z[z<=n])-th(2*n)/np.pi for n in [15,22,29.9]]
    print(f" k={k}: {len(z)} levels below 30, first {np.round(z[:4],4)}; rule residual at last 3 {np.round(rr[-3:],4)}; N-theta(2nu)/pi at nu=15,22,29.9: {np.round(cnt,3)}")
