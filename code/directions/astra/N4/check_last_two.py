"""Independent bounded audit calculations. No manuscript-directory imports or reads."""
import json
import math
import time
from pathlib import Path
import mpmath as mp
import numpy as np

started = time.time()
mp.mp.dps = 45
pi = mp.pi

def phi(u):
    return mp.fsum((2*pi*pi*n**4*mp.exp(9*u)-3*pi*n*n*mp.exp(5*u))*
                   mp.exp(-pi*n*n*mp.exp(4*u)) for n in range(1,13))

def potential(u):
    return -mp.log(phi(u))

def star_potential(u):
    return 2*pi*mp.cosh(4*u)-mp.log(2*pi*pi*mp.cosh(9*u)-3*pi*mp.cosh(5*u))

def xi(s):
    if s == 0 or s == 1:
        return mp.mpf('0.5')
    return s*(s-1)*mp.power(pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)/2

xic = xi(mp.mpf('0.5'))
def q(w):
    return xic/xi(mp.mpf('0.5')+w/2)

out = {'method':'mpmath 45 digits; 12 theta terms; floating point, not interval arithmetic'}
out['U2_at_0']=str(mp.diff(potential,0,2))
out['U4_at_0']=str(mp.diff(potential,0,4))
out['Ustar2_at_0']=str(mp.diff(star_potential,0,2))
out['Ustar4_at_0']=str(mp.diff(star_potential,0,4))
grid=[mp.mpf(j)/50 for j in range(1,101)]
u3=[mp.diff(potential,u,3) for u in grid]
s3=[mp.diff(star_potential,u,3) for u in grid]
out['U3_grid']={'min':str(min(u3)),'argmin':str(grid[u3.index(min(u3))]),'at2':str(u3[-1])}
out['Ustar3_grid']={'min':str(min(s3)),'argmin':str(grid[s3.index(min(s3))]),'at2':str(s3[-1])}
out['Ustar3_at_0_01']=str(mp.diff(star_potential,mp.mpf('0.01'),3))
out['variance']=str(-mp.diff(q,0,2))
out['q400']=str(q(mp.mpf(400)))
zeros=[mp.im(mp.zetazero(j)) for j in range(1,4)]
coeff=[2*xic/abs(mp.diff(xi,mp.mpf('0.5')+1j*g)) for g in zeros]
out['first_three_zeros']=[str(g) for g in zeros]
out['first_three_residue_coefficients']=[str(c) for c in coeff]
out['three_term_spline']={str(x):str(mp.fsum((-1)**j*coeff[j]*mp.exp(-2*zeros[j]*x) for j in range(3))) for x in [mp.mpf('.5'),mp.mpf('.55'),mp.mpf('.6')]}

# Fixed Gaussian quadrature on the reciprocal transform, followed by cheap cosine sums.
nodes,weights=np.polynomial.legendre.leggauss(24)
ww=np.concatenate([5*k+2.5*(nodes+1) for k in range(80)])
gw=np.tile(2.5*weights,80)
qw=np.array([float(q(mp.mpf(float(w)))) for w in ww])
base=gw*qw/math.pi
def density(x,order=0):
    x=np.atleast_1d(x)
    if order==0:
        return np.cos(np.outer(x,ww)) @ base
    if order==1:
        return -np.sin(np.outer(x,ww)) @ (base*ww)
    return -np.cos(np.outer(x,ww)) @ (base*ww*ww)
out['fourier_density']={str(x):float(density(x)[0]) for x in [.5,.55,.6]}
xgrid=np.linspace(0,.6,301)
p,p1,p2=density(xgrid),density(xgrid,1),density(xgrid,2)
lc=p2/p-(p1/p)**2
out['log_curvature_grid']={'max':float(lc.max()),'argmax':float(xgrid[lc.argmax()]),'at0':float(lc[0])}

# Deterministic separated configurations avoid unreproducible random minima.
out['fixed_minors']={}
for n in range(2,6):
    xx=np.linspace(-.20,.20,n)
    yy=np.linspace(-.19,.21,n)
    mat=density((xx[:,None]-yy[None,:]).ravel()).reshape(n,n)
    out['fixed_minors'][str(n)]={'normalized':float(np.linalg.det(mat)/np.diag(mat).prod())}

# Independent density quadrature for the two total variation distances.
gx,wt=np.polynomial.legendre.leggauss(400)
uu=.6*(gx+1)
uw=.6*wt
pf=np.array([float(phi(mp.mpf(float(u)))) for u in uu])
sf=np.array([float(mp.exp(-star_potential(mp.mpf(float(u))))) for u in uu])
zp=2*np.dot(uw,pf)
zs=2*np.dot(uw,sf)
vp=2*np.dot(uw,pf*uu**2)/zp
vs=2*np.dot(uw,sf*uu**2)/zs
scale=math.sqrt(vp/vs)
sm=np.array([float(mp.exp(-star_potential(mp.mpf(float(u/scale)))))/scale for u in uu])
out['star_comparison']={'var_polya':float(vp),'var_star':float(vs),'scale':scale,
                        'TV':float(np.dot(uw,np.abs(pf/zp-sf/zs))),
                        'TV_matched':float(np.dot(uw,np.abs(pf/zp-sm/zs)))}

# Davenport and Heilbronn: independent Hurwitz-zeta expression and off-axis zero.
kappa=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
def dh(s):
    return mp.power(5,-s)*mp.fsum(a*mp.zeta(s,mp.mpf(j)/5) for j,a in [(1,1),(2,kappa),(3,-kappa),(4,-1)])
root=mp.findroot(dh,(mp.mpc('.80','85.69'),mp.mpc('.82','85.71')))
out['DH']={'root':str(root),'residual':str(abs(dh(root))),
           'reciprocal_fourier_pole_w':str(2*(root-mp.mpf('.5'))),
           'reciprocal_laplace_pole_s':str(-2j*(root-mp.mpf('.5')))}

# Independent high precision checks on the contested tail and curvature numbers.
mp.mp.dps=65
xic=xi(mp.mpf('.5'))
xx=mp.mpf('.6')
breaks=[0,1,5,10,20,40,80,120,200,400]
p0=mp.quad(lambda w:mp.cos(w*xx)*q(w),breaks)/pi
p1=mp.quad(lambda w:-w*mp.sin(w*xx)*q(w),breaks)/pi
p2=mp.quad(lambda w:-w*w*mp.cos(w*xx)*q(w),breaks)/pi
out['multiprecision_at_0_6']={'density':str(p0),'log_curvature':str(p2/p0-(p1/p0)**2),
                             'three_term_relative_error':str((mp.mpf(out['three_term_spline']['0.6'])-p0)/p0)}
out['variance_discrepancy_vs_note']=str(mp.mpf(out['variance'])-mp.mpf('.0115516'))

# A Lee Yang probability density that violates GHS at a positive field.
# Its transform is exp(b*z*z)*(1+z*z/a^2), with a=b=1.
# Its density is the N(0,2) density multiplied by 1/2+x^2/4.
out['lee_yang_not_GHS']={'a':1,'b':1,'K3_at_2':str(mp.diff(lambda z:z*z+mp.log(1+z*z),mp.mpf(2),3))}
out['elapsed_seconds']=time.time()-started
target=Path(__file__).with_name('check_last_two_results.json')
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
