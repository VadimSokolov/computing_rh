import numpy as np, matplotlib, textwrap
matplotlib.use('Agg'); matplotlib.rcParams['mathtext.fontset']='cm'
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
# Usage (from cover/): python3 coverwrap.py [pages]; the page count defaults to that of ../main-kdp.pdf, the KDP interior.
# plotcurves.npz is written by code/analyse.py; the copy here comes from the September 2026 rerun on Hopper.
import sys, subprocess, re
C=np.load('plotcurves.npz'); Z=C['zeros']
if len(sys.argv)>1: PAGES=int(sys.argv[1])
else: PAGES=int(re.search(r'Pages:\s+(\d+)',subprocess.run(['pdfinfo','../main-kdp.pdf'],capture_output=True,text=True).stdout).group(1))
PPI=0.002252                        # KDP white paper, per page
TRIM_W,TRIM_H,BL=7.0,10.0,0.125
K=TRIM_W/8.5                        # the design was drawn for an 8.5 in wide panel; type and spacing scale with the width
FADE=lambda c,a: tuple(a*np.array(matplotlib.colors.to_rgb(c))+(1-a)*np.array([0.52,0.07,0.08]))   # alpha blended into the background, so the PDF has no transparency
S=round(PAGES*PPI,4)                # spine width in inches
W=2*BL+2*TRIM_W+S; H=TRIM_H+2*BL
print('pages',PAGES,'spine',S,'full size',round(W,4),H)
fig=plt.figure(figsize=(W,H),dpi=300)
bg=fig.add_axes([0,0,1,1]); bg.set_axis_off(); bg.set_xlim(0,W); bg.set_ylim(0,H)
# matte red background with soft vignettes centred on each panel
y,x=np.mgrid[0:H:700j,0:W:1100j]
backc=BL+TRIM_W/2; frontc=BL+TRIM_W+S+TRIM_W/2
d=np.minimum(np.hypot((x-backc)/TRIM_W,(y-0.55*H)/TRIM_H),np.hypot((x-frontc)/TRIM_W,(y-0.55*H)/TRIM_H))
t=np.clip(d/0.75,0,1)[...,None]
base=np.array([0.52,0.07,0.08]); edge=np.array([0.36,0.04,0.05])
bg.imshow(base*(1-t**2)+edge*t**2,extent=[0,W,0,H],aspect='auto',origin='lower',zorder=0)
cream='#F3E6C8'; gold='#E0B04F'
def lmfont(fname,family):          # the Latin Modern OpenType files of the TeX installation, else the installed family
    p=subprocess.run(['kpsewhich',fname],capture_output=True,text=True).stdout.strip()
    return fm.FontProperties(fname=p) if p else fm.FontProperties(family=family)
caps=lmfont('lmromancaps10-regular.otf','Latin Modern Roman Caps'); serif=lmfont('lmroman10-regular.otf','Latin Modern Roman'); ital=lmfont('lmromanslant10-regular.otf','Latin Modern Roman Slanted')
# ---------------- front panel ----------------
fx0=BL+TRIM_W+S; fxc=fx0+TRIM_W/2
def fy(v): return BL+v*TRIM_H     # fraction of trim height
bg.text(fxc,fy(0.905),'Computing the',ha='center',va='center',color=cream,fontsize=40*K,fontproperties=caps)
bg.text(fxc,fy(0.845),'Riemann Hypothesis',ha='center',va='center',color=cream,fontsize=50*K,fontproperties=caps)
bg.plot([fxc-0.22*TRIM_W,fxc+0.22*TRIM_W],[fy(0.805)]*2,color=gold,lw=1.2)
bg.text(fxc,fy(0.772),'Thorin Measures, P\u00f3lya Densities and Hitting Times',ha='center',va='center',color=gold,fontsize=19*K,fontproperties=ital)
ax=fig.add_axes([(fx0+0.1*TRIM_W)/W,fy(0.23)/H,0.8*TRIM_W/W,0.49*TRIM_H/H]); ax.set_axis_off()
eps=[1.0,0.5,0.25,0.1,0.05,0.02,0.01,0.005]; cmap=LinearSegmentedColormap.from_list('g',['#C98A3A','#F7EBD0'])
tmin,tmax=10,60
for i,e in list(enumerate(eps))[::-1]:
    th=C['th_%s'%e]; f=C['F_%s'%e]/np.pi; m=(th>=tmin)&(th<=tmax)
    lv=np.log10(f[m]); yv=(lv-lv.min())/(lv.max()-lv.min()); y0=i*0.55
    ax.fill_between(th[m],y0,y0+yv*1.6,color=(0.47,0.06,0.07),zorder=2*(len(eps)-i))
    ax.plot(th[m],y0+yv*1.6,color=cmap(i/(len(eps)-1)),lw=0.8+0.14*i,zorder=2*(len(eps)-i)+1)
for z in Z[(Z>tmin)&(Z<tmax)]: ax.plot([z,z],[-0.25,-0.05],color=FADE(gold,0.9),lw=1.0)
ax.set_xlim(tmin,tmax); ax.set_ylim(-0.4,(len(eps)-1)*0.55+1.75)
bg.text(fxc,fy(0.212),r'the Thorin density of $\xi(\frac{1}{2})/\xi(\frac{1}{2}+\sqrt{u})$ sharpening onto the zeros',ha='center',va='center',color=FADE(cream,0.85),fontsize=11.5*K,fontproperties=ital)
bg.plot([fxc-0.14*TRIM_W,fxc+0.14*TRIM_W],[fy(0.165)]*2,color=gold,lw=0.8)
bg.text(fxc,fy(0.125),'Nicholas G. Polson',ha='center',va='center',color=cream,fontsize=24*K,fontproperties=serif)
bg.text(fxc,fy(0.082),'Vadim Sokolov',ha='center',va='center',color=cream,fontsize=24*K,fontproperties=serif)
# ---------------- spine ----------------
sx=BL+TRIM_W+S/2
bg.text(sx,BL+TRIM_H*0.60,'Computing the Riemann Hypothesis',rotation=-90,ha='center',va='center',color=cream,fontsize=20,fontproperties=caps)
bg.text(sx,BL+TRIM_H*0.15,'Polson  \u00b7  Sokolov',rotation=-90,ha='center',va='center',color=gold,fontsize=14,fontproperties=serif)
# ---------------- back panel ----------------
bx0=BL; bxc=BL+TRIM_W/2; lm=BL+0.75; tw=TRIM_W-1.5
bg.text(bxc,fy(0.915),'Where computation meets the Riemann Hypothesis',ha='center',va='center',color=gold,fontsize=20*K,fontproperties=ital)
blurb=("The Riemann Hypothesis can be stated as a positivity condition. The reciprocal Xi function is the Laplace transform "
"of a generalised gamma convolution exactly when every zero lies on the critical line, and then its Thorin measure is the "
"discrete measure on the squared ordinates. This book asks what computation can say about that condition.",
"It develops the P\u00f3lya route, which evaluates the Xi function in certified ball arithmetic from P\u00f3lya's kernel alone, "
"and uses it to compute the Thorin measure at every resolution, to measure its distance to the measure demanded by RH, "
"to read the zero density from a Thorin form of Turing's method, and to calibrate exactly what an off line zero would leave behind. "
"The sine squared lemma emerges as the bridge from the primes to the zeros, and as the reason the primes alone cannot compute near the line.",
"Hitting times and Brownian boundaries, Keiper Li coefficients, the de Bruijn Newman flow, the Dedekind eta kernel of $L(s,\\chi_{12})$ "
"and a study of the reciprocal zeta function complete the picture, and Chapter 12 sets out a path to RH with every step marked "
"proved, computed or open.",
"The last chapter reads RH as a statement of statistical mechanics: the P\u00f3lya law, whose characteristic function is "
"the normalised Xi function, is ferromagnetic in the sense of Lee and Yang. With no hypothesis on the zeros, it builds "
"a law within total variation 0.000002 of the P\u00f3lya law and proves by computer, in ball arithmetic, that this law is "
"ferromagnetic. The code and data behind the book's computations are at github.com/VadimSokolov/computing_rh.")
yy=fy(0.865); ASC=0.105                # lines sit on a common baseline, ASC below the top of the first, so inline math does not shift them
for para in blurb:
    for line in textwrap.wrap(para,76):
        bg.text(lm,yy-ASC,line,ha='left',va='baseline',color=cream,fontsize=12.5*K,fontproperties=serif); yy-=0.235*K
    yy-=0.14*K
# small Thorin flow heat map strip
yy-=0.05
axh=fig.add_axes([lm/W,(yy-1.55*K)/H,tw/W,1.45*K/H])
tg=np.linspace(10,60,2000)
img=np.array([np.interp(tg,C['th_%s'%e],np.log10(C['F_%s'%e]/np.pi)) for e in eps])
hm=LinearSegmentedColormap.from_list('h',['#3a0507','#8c1c1f','#d08a3a','#f7ebd0'])
axh.pcolormesh(tg,np.arange(len(eps)),img,shading='auto',cmap=hm); axh.set_axis_off()
bg.text(bxc,yy-1.72*K,r'the same measure at eight resolutions, from $\varepsilon=1$ (bottom) to $\varepsilon=0.005$ (top)',ha='center',va='top',color=FADE(cream,0.85),fontsize=10*K,fontproperties=ital)
yy-=2.15*K
bios=["Nicholas G. Polson is the Robert Law Jr. Professor of Econometrics and Statistics at the Booth School of Business, University of Chicago.",
"Vadim Sokolov is on the faculty of the Department of Systems Engineering and Operations Research at George Mason University."]
bg.plot([lm,lm+tw*0.35],[yy+0.08]*2,color=gold,lw=0.8)
for b in bios:
    for line in textwrap.wrap(b,62):
        bg.text(lm,yy,line,ha='left',va='top',color=cream,fontsize=11*K,fontproperties=serif); yy-=0.21*K
    yy-=0.1*K
# barcode area (kept clear; KDP places the barcode here)
bx1=BL+TRIM_W-0.25-2.0; by1=BL+0.25
print('lowest text y',yy,'barcode top',by1+1.2)
fig.savefig('cover_full.pdf',dpi=300); fig.savefig('cover_full.png',dpi=300)
# guide overlay for checking only (not for upload)
g=[]
for xv in [BL, BL+TRIM_W, BL+TRIM_W+S, W-BL]: g.append(bg.axvline(xv,color='cyan',lw=0.8))
for yv in [BL, H-BL]: g.append(bg.axhline(yv,color='cyan',lw=0.8))
from matplotlib.patches import Rectangle
r=Rectangle((BL+TRIM_W-0.25-2.0,BL+0.25),2.0,1.2,fill=False,ec='white',lw=1,ls='--'); bg.add_patch(r)
bg.text(BL+TRIM_W-1.25,BL+0.85,'barcode area',color='white',ha='center',va='center',fontsize=9)
fig.savefig('cover_full_guides.png',dpi=100)
print('ok')
