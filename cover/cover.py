import numpy as np, matplotlib
matplotlib.rcParams['mathtext.fontset']='cm'
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
C=np.load('plotcurves.npz'); Z=C['zeros']
import sys
MODE=sys.argv[1] if len(sys.argv)>1 else 'print'
W,H=(8.625,11.25) if MODE=='print' else (6.25,10.0)
K=W/8.625
fig=plt.figure(figsize=(W,H),dpi=300)
bg=fig.add_axes([0,0,1,1]); bg.set_axis_off()
# matte red with a very soft vignette
y,x=np.mgrid[0:1:600j,0:1:460j]
r=np.sqrt((x-0.5)**2+(y-0.55)**2)
base=np.array([0.52,0.07,0.08]); edge=np.array([0.36,0.04,0.05])
t=np.clip(r/0.75,0,1)[...,None]
img=base*(1-t**2)+edge*t**2
bg.imshow(img,extent=[0,1,0,1],aspect='auto',origin='lower')
bg.set_xlim(0,1); bg.set_ylim(0,1)
cream='#F3E6C8'; gold='#E0B04F'
titlefont=fm.FontProperties(family='Latin Modern Roman Caps' if any('Latin Modern Roman Caps' in f.name for f in fm.fontManager.ttflist) else 'DejaVu Serif')
serif=fm.FontProperties(family='Latin Modern Roman')
ital=fm.FontProperties(family='Latin Modern Roman Slanted')
# title block
bg.text(0.5,0.905,'Computing the',ha='center',va='center',color=cream,fontsize=40*K,fontproperties=titlefont)
bg.text(0.5,0.845,'Riemann Hypothesis',ha='center',va='center',color=cream,fontsize=50*K,fontproperties=titlefont)
bg.plot([0.28,0.72],[0.805,0.805],color=gold,lw=1.2)
bg.text(0.5,0.772,'Thorin Measures, P\u00f3lya Densities and Hitting Times',ha='center',va='center',color=gold,fontsize=19*K,fontproperties=ital)
# the Thorin plot: ridgelines of the Pólya computed Thorin density
ax=fig.add_axes([0.1,0.23,0.8,0.49]); ax.set_axis_off()
eps=[1.0,0.5,0.25,0.1,0.05,0.02,0.01,0.005]
cmap=LinearSegmentedColormap.from_list('g',['#C98A3A','#F7EBD0'])
tmin,tmax=10,60; off=1.05
order=list(enumerate(eps))[::-1]   # sharpest at the back (top), coarsest in front (bottom)
for i,e in order:
    th=C['th_%s'%e]; f=C['F_%s'%e]/np.pi; m=(th>=tmin)&(th<=tmax)
    lv=np.log10(f[m]); yv=(lv-lv.min())/(lv.max()-lv.min())
    y0=i*0.55
    col=cmap(i/(len(eps)-1))
    ax.fill_between(th[m],y0,y0+yv*1.6,color=(0.47,0.06,0.07),zorder=2*(len(eps)-i))
    ax.plot(th[m],y0+yv*1.6,color=col,lw=0.8+0.14*i,zorder=2*(len(eps)-i)+1)
for z in Z[(Z>tmin)&(Z<tmax)]:
    ax.plot([z,z],[-0.25,-0.05],color=gold,lw=1.0,alpha=0.9)
ax.set_xlim(tmin,tmax); ax.set_ylim(-0.4,(len(eps)-1)*0.55+1.75)
bg.text(0.5,0.212,r'the Thorin density of $\xi(\frac{1}{2})/\xi(\frac{1}{2}+\sqrt{u})$ sharpening onto the zeros',ha='center',va='center',color=cream,fontsize=11.5*K,alpha=0.85,fontproperties=ital)
# authors
bg.plot([0.36,0.64],[0.165,0.165],color=gold,lw=0.8)
bg.text(0.5,0.125,'Nicholas G. Polson',ha='center',va='center',color=cream,fontsize=24*K,fontproperties=serif)
bg.text(0.5,0.082,'Vadim Sokolov',ha='center',va='center',color=cream,fontsize=24*K,fontproperties=serif)
if MODE=='print':
    fig.savefig('cover_front.png',dpi=300); fig.savefig('cover_front.pdf',dpi=300)
else:
    fig.savefig('cover_ebook.jpg',dpi=256,pil_kwargs={'quality':95})
print('ok')
