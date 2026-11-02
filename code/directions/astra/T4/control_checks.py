"""Additional independent checks and derived numbers for report.md."""
from pathlib import Path
import sys
import json
sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import mpmath as mp
from flint import arb, acb, acb_series, ctx

BASE = Path(__file__).resolve().parent
ctx.dps = 75
mp.mp.dps = 65


def mobius(nmax):
    mu=[1]*(nmax+1)
    prime=[True]*(nmax+1)
    prime[0]=prime[1]=False
    for p in range(2,nmax+1):
        if prime[p]:
            for n in range(p,nmax+1,p):
                mu[n] *= -1
                if n>p:
                    prime[n]=False
            for n in range(p*p,nmax+1,p*p):
                mu[n]=0
    return mu


def F(a,L):
    return 1/a-(1-(-a*L).exp())/(L*a*a)


def theta_density_kernel(u):
    return sum((4*mp.pi**2*n**4*mp.exp(mp.mpf("4.5")*u)
                -6*mp.pi*n*n*mp.exp(mp.mpf("2.5")*u))
               *mp.exp(-mp.pi*n*n*mp.exp(2*u)) for n in range(1,11))


flow=json.loads((BASE/"flow_results.json").read_text())
z=flow["davenport_heilbronn"]
q=acb(arb(z["center_imag"],"1e-20"),
      -(arb(z["center_real"],"1e-20")-arb(1)/2))
mu=mobius(500)
checks=[]
for K in (20,100,500):
    t=arb(10)
    prod=acb(1)
    for k in range(1,K+1):
        if mu[k]:
            v=arb.pi()*t/(k*q)
            prod *= (v.sin()/v)**mu[k]
    exact=1-(t/q)**2
    mk=6/(arb.pi()**2)-sum(arb(mu[k])/(k*k) for k in range(1,K+1))
    bound=arb.pi()**2*t*t/6*abs(mk)*abs(1/q**2)
    bound+=2*arb(4).zeta()/(9*K**3)*t**4*abs(1/q**4)
    absolute_bound=abs(exact)*(bound.exp()-1)
    error=abs(prod-exact)
    assert error < absolute_bound
    checks.append({"K":K,"t":10,"product_ball":str(prod),
                   "algebraic_target_ball":str(exact),
                   "absolute_error_ball":str(error),
                   "proved_truncation_bound_ball":str(absolute_bound)})

L=arb(20)
I0=F(arb(1),L)
I2=2-6/L+(-L).exp()*(6/L+4+L)
quad=I2/(2*I0)
small=json.loads((BASE/"smallball_results.json").read_text())["rows"]
balls=json.loads((BASE/"ball_inversion_results.json").read_text())
derived=[]
for i in range(5):
    derived.append({"epsilon":small[i]["epsilon"],
                    "I2_minus_I1_ball":str(arb(small[i+5]["I_ball"])-arb(small[i]["I_ball"]))})
ratios=[str(arb(small[i]["saddle_probability_ball"])/arb(balls[i]["probability_ball"])) for i in range(2)]

fourier=[]
for t in (0,3,10):
    integral=2*mp.quad(lambda u:theta_density_kernel(u)*mp.cos(t*u),[0,.25,.5,1,2,3])
    s=mp.mpc("0.5",t)
    exact=mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
    fourier.append({"t":t,"two_sided_kernel_transform":mp.nstr(integral,50),
                    "xi_direct":mp.nstr(exact,50),"difference":mp.nstr(integral-exact,12)})

def dh_level_two():
    """Whole finite control flow at a level with all positive widths."""
    nmax=10000
    k=((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)
    aa=[arb(0),arb(1),k,-k,arb(-1)]
    coeff=[arb(0)]*(nmax+1)
    for n in range(2,nmax+1):
        coeff[n]=aa[n%5]*arb(n).log()
    for d in range(2,nmax//2+1):
        for m in range(2,nmax//d+1):
            coeff[d*m]-=coeff[d]*aa[m%5]
    length=arb(nmax).log()
    delta=arb("1.5")
    w=arb(z["center_imag"],"1e-20")
    a=acb(delta+arb("1.5"),-w)
    gam=acb_series([a/2,1],prec=3).lgamma()
    arch=(5/arb.pi()).log()+gam[1].real+gam[2].real/length
    stop=4
    for n in range(stop):
        b=a+2*n
        arch-=2/length*((-b*length).exp()/(b*b)).real
    first=delta+arb("1.5")+2*stop
    tail=2/length*(-first*length).exp()/(first*first*(1-(-2*length).exp()))
    arch+=arb(0,tail.upper())
    prime=arb(0)
    for n in range(2,nmax+1):
        ln=arb(n).log()
        prime+=2*coeff[n]*(-2*ln).exp()*(1-ln/length)*(w*ln).cos()
    value=arch-prime
    h0=2*F(delta,length)
    assert value/h0>2
    return {"alpha":2,"delta":"1.5","L":str(length),
            "w_ball":str(w),"whole_flow_ball":str(value),
            "H0_ball":str(h0),"normalized_flow_ball":str(value/h0),
            "criterion_rejects_known_offline_ordinate":True}


data={"status":"computed","DH_complex_Mobius_balls":checks,
      "positive_width_DH_control":dh_level_two(),
      "finite_L_quadratic_coefficient_delta1_L20_ball":str(quad),
      "coefficient_is_less_than_one":bool(quad<1),
      "smallball_rate_differences":derived,"saddle_probability_ratios_ball":ratios,
      "kernel_normalization_checks_mpmath":fourier}
(BASE/"control_results.json").write_text(json.dumps(data,indent=2)+"\n")
print("Wrote control_results.json")
