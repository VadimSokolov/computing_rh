import numpy as np, matplotlib, textwrap, sys, subprocess, re
matplotlib.use('Agg'); matplotlib.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
# Wrapped KDP covers (back, spine, front) for the two volume edition, 7 by 10 inches, white paper.
# Usage (from cover/): python3 cover_volumes.py 1|2 [pages]. The page count defaults to that of ../main_vol{1,2}-kdp.pdf.
# Writes cover_vol{1,2}_full.pdf (upload), .png (preview) and _guides.png (trim, spine and barcode guides; not for upload).
VOL = int(sys.argv[1])
C = np.load('plotcurves.npz'); Z = C['zeros']
if len(sys.argv) > 2: PAGES = int(sys.argv[2])
else: PAGES = int(re.search(r'Pages:\s+(\d+)', subprocess.run(['pdfinfo', '../main_vol%d-kdp.pdf' % VOL], capture_output=True, text=True).stdout).group(1))
PPI = 0.002252                      # KDP white paper, inches per page
TRIM_W, TRIM_H, BL = 7.0, 10.0, 0.125
K = TRIM_W/8.5
BASE = {1: np.array([0.52, 0.07, 0.08]), 2: np.array([0.07, 0.15, 0.32])}[VOL]
EDGE = {1: np.array([0.36, 0.04, 0.05]), 2: np.array([0.03, 0.08, 0.20])}[VOL]
FADE = lambda c, a: tuple(a*np.array(matplotlib.colors.to_rgb(c)) + (1-a)*BASE)
S = round(PAGES*PPI, 4)
W = 2*BL + 2*TRIM_W + S; H = TRIM_H + 2*BL
print('volume', VOL, 'pages', PAGES, 'spine', S, 'full size', round(W, 4), H)
fig = plt.figure(figsize=(W, H), dpi=300)
bg = fig.add_axes([0, 0, 1, 1]); bg.set_axis_off(); bg.set_xlim(0, W); bg.set_ylim(0, H)
y, x = np.mgrid[0:H:700j, 0:W:1100j]
backc = BL + TRIM_W/2; frontc = BL + TRIM_W + S + TRIM_W/2
d = np.minimum(np.hypot((x-backc)/TRIM_W, (y-0.55*H)/TRIM_H), np.hypot((x-frontc)/TRIM_W, (y-0.55*H)/TRIM_H))
t = np.clip(d/0.75, 0, 1)[..., None]
bg.imshow(BASE*(1-t**2) + EDGE*t**2, extent=[0, W, 0, H], aspect='auto', origin='lower', zorder=0)
cream = '#F3E6C8'; gold = '#E0B04F'
def lmfont(fname, family):
    p = subprocess.run(['kpsewhich', fname], capture_output=True, text=True).stdout.strip()
    return fm.FontProperties(fname=p) if p else fm.FontProperties(family=family)
caps = lmfont('lmromancaps10-regular.otf', 'Latin Modern Roman Caps'); serif = lmfont('lmroman10-regular.otf', 'Latin Modern Roman')
ital = lmfont('lmromanslant10-regular.otf', 'Latin Modern Roman Slanted')
eps = [1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005]
tmin, tmax = 10, 60
SUB = {1: 'Representations of Xi, Thorin Measures\nand the Reciprocal Xi Function',
       2: 'Cauchy Deconvolution\nand Certified Computation'}[VOL]
VOLNAME = {1: 'Volume I', 2: 'Volume II'}[VOL]
# ---------------- front panel ----------------
fx0 = BL + TRIM_W + S; fxc = fx0 + TRIM_W/2
def fy(v): return BL + v*TRIM_H
bg.text(fxc, fy(0.915), 'Computing the', ha='center', va='center', color=cream, fontsize=40*K, fontproperties=caps)
bg.text(fxc, fy(0.858), 'Riemann Hypothesis', ha='center', va='center', color=cream, fontsize=50*K, fontproperties=caps)
bg.text(fxc, fy(0.805), VOLNAME, ha='center', va='center', color=gold, fontsize=24*K, fontproperties=caps)
bg.plot([fxc-0.22*TRIM_W, fxc+0.22*TRIM_W], [fy(0.776)]*2, color=gold, lw=1.2)
bg.text(fxc, fy(0.735), SUB, ha='center', va='center', color=gold, fontsize=18*K, fontproperties=ital, linespacing=1.35)
ax = fig.add_axes([(fx0+0.1*TRIM_W)/W, fy(0.235)/H, 0.8*TRIM_W/W, 0.445*TRIM_H/H]); ax.set_axis_off()
if VOL == 1:
    cmap = LinearSegmentedColormap.from_list('g', ['#C98A3A', '#F7EBD0'])
    for i, e in list(enumerate(eps))[::-1]:
        th = C['th_%s' % e]; f = C['F_%s' % e]/np.pi; m = (th >= tmin) & (th <= tmax)
        lv = np.log10(f[m]); yv = (lv-lv.min())/(lv.max()-lv.min()); y0 = i*0.55
        ax.fill_between(th[m], y0, y0+yv*1.6, color=FADE('#000000', 0.15), zorder=2*(len(eps)-i))
        ax.plot(th[m], y0+yv*1.6, color=cmap(i/(len(eps)-1)), lw=0.8+0.14*i, zorder=2*(len(eps)-i)+1)
    for z in Z[(Z > tmin) & (Z < tmax)]: ax.plot([z, z], [-0.25, -0.05], color=FADE(gold, 0.9), lw=1.0)
    ax.set_xlim(tmin, tmax); ax.set_ylim(-0.4, (len(eps)-1)*0.55+1.75)
    cap = r'the Thorin density of $\xi(\frac{1}{2})/\xi(\frac{1}{2}+\sqrt{u})$ sharpening onto the zeros'
else:
    # the backward Cauchy flow: the broad level one density (top) deconvolved towards the zeros (bottom)
    cmap = LinearSegmentedColormap.from_list('b', ['#7FB3E0', '#F3E6C8'])
    for i, e in enumerate(eps):
        th = C['th_%s' % e]; f = C['F_%s' % e]/np.pi; m = (th >= tmin) & (th <= tmax)
        lv = np.log10(f[m]); yv = (lv-lv.min())/(lv.max()-lv.min()); y0 = (len(eps)-1-i)*0.55
        ax.fill_between(th[m], y0, y0+yv*1.6, color=FADE('#000000', 0.15), zorder=2*i)
        ax.plot(th[m], y0+yv*1.6, color=cmap(i/(len(eps)-1)), lw=0.8+0.14*i, zorder=2*i+1)
    for z in Z[(Z > tmin) & (Z < tmax)]: ax.plot([z, z], [-0.25, -0.05], color=FADE(gold, 0.9), lw=1.0)
    ax.annotate('', xy=(61.2, 0.3), xytext=(61.2, (len(eps)-1)*0.55+1.2), arrowprops=dict(arrowstyle='->', color=gold, lw=1.0))
    ax.set_xlim(tmin, 62); ax.set_ylim(-0.4, (len(eps)-1)*0.55+1.75)
    cap = r'the backward Cauchy flow: from the level one density down to the zeros'
bg.text(fxc, fy(0.212), cap, ha='center', va='center', color=FADE(cream, 0.85), fontsize=11.5*K, fontproperties=ital)
bg.plot([fxc-0.14*TRIM_W, fxc+0.14*TRIM_W], [fy(0.165)]*2, color=gold, lw=0.8)
bg.text(fxc, fy(0.125), 'Nicholas G. Polson', ha='center', va='center', color=cream, fontsize=24*K, fontproperties=serif)
bg.text(fxc, fy(0.082), 'Vadim Sokolov', ha='center', va='center', color=cream, fontsize=24*K, fontproperties=serif)
# ---------------- spine ----------------
sx = BL + TRIM_W + S/2
bg.text(sx, BL+TRIM_H*0.585, 'Computing the Riemann Hypothesis', rotation=-90, ha='center', va='center', color=cream, fontsize=19, fontproperties=caps)
bg.text(sx, BL+TRIM_H*0.27, VOLNAME, rotation=-90, ha='center', va='center', color=gold, fontsize=16, fontproperties=caps)
bg.text(sx, BL+TRIM_H*0.11, 'Polson  ·  Sokolov', rotation=-90, ha='center', va='center', color=gold, fontsize=14, fontproperties=serif)
# ---------------- back panel ----------------
bxc = BL + TRIM_W/2; lm = BL + 0.75; tw = TRIM_W - 1.5
TAG = {1: 'The Riemann Hypothesis as the positivity of a measure',
       2: 'Running the Cauchy flow backwards'}[VOL]
BLURB = {1: (
"The Riemann Hypothesis can be stated as a positivity condition. The reciprocal Xi function is the Laplace transform of a "
"generalised gamma convolution exactly when every zero lies on the critical line, and then its Thorin measure is the "
"discrete measure on the squared ordinates. This volume is the theory of that condition.",
"It gives every representation of the Xi function the two volumes use, the classes of infinitely divisible laws in which "
"they live, from generalised gamma convolutions and their extended form to hyperbolic complete monotonicity, and the "
"reciprocal Xi function at every basepoint. It follows the Thorin condition through the reciprocal zeta function, the "
"critical strip, Dirichlet L functions and the Dedekind eta function, and through hitting times, extreme values and "
"ferromagnetic random variables.",
"The last chapters read the Thorin measure as a mixture of Cauchy laws, prove a strip version of the Lee and Yang theorem "
"for spins of Xi type, and assemble a route to the Riemann Hypothesis in six equivalent forms, each proved at the line "
"Re s = 1, with curves over finite fields as the benchmark where the same route is a proof. Every statement is marked "
"proved, computed, certified, conditional, heuristic or open."),
2: (
"The Thorin measure of the reciprocal Xi function at any level is the Cauchy smoothing of the zeros, and at level one it "
"is computable from the Euler product. The Riemann Hypothesis is the statement that this smoothing can be undone, down to "
"the critical line, through positive measures. This volume computes that inverse problem as well as it can be done.",
"It evaluates the Xi function in ball arithmetic from P\u00f3lya's kernel alone, computes the de Bruijn Newman flow and its "
"barrier method in Thorin form, and runs the backward Cauchy flow from the level one density, locating zeros and "
"certifying window masses, by particle filters, nested sampling and the Monte Carlo criterion of Balazard, Saias and Yor.",
"It certifies Weil positivity in windows, follows truncations, friable, ferromagnetic and extreme value approximants and "
"their zeros, and gathers the certificates that hold under the quasi Riemann hypothesis announced in October 2026. Each "
"chapter says what a finite computation can decide and why no finite computation can decide the hypothesis itself. The "
"code and data are at github.com/VadimSokolov/computing_rh.")}[VOL]
bg.text(bxc, fy(0.915), TAG, ha='center', va='center', color=gold, fontsize=20*K, fontproperties=ital)
yy = fy(0.865); ASC = 0.105
for para in BLURB:
    for line in textwrap.wrap(para, 76):
        bg.text(lm, yy-ASC, line, ha='left', va='baseline', color=cream, fontsize=12.5*K, fontproperties=serif); yy -= 0.235*K
    yy -= 0.14*K
yy -= 0.05
axh = fig.add_axes([lm/W, (yy-1.55*K)/H, tw/W, 1.45*K/H])
tg = np.linspace(10, 60, 2000)
img = np.array([np.interp(tg, C['th_%s' % e], np.log10(C['F_%s' % e]/np.pi)) for e in eps])
if VOL == 1: hm = LinearSegmentedColormap.from_list('h', ['#3a0507', '#8c1c1f', '#d08a3a', '#f7ebd0'])
else: hm = LinearSegmentedColormap.from_list('h', ['#061226', '#173a6b', '#7FB3E0', '#f3e6c8'])
axh.pcolormesh(tg, np.arange(len(eps)), img, shading='auto', cmap=hm); axh.set_axis_off()
bg.text(bxc, yy-1.72*K, r'the same measure at eight resolutions, from $\varepsilon=1$ (bottom) to $\varepsilon=0.005$ (top)', ha='center', va='top', color=FADE(cream, 0.85), fontsize=10*K, fontproperties=ital)
yy -= 2.15*K
bios = ["Nicholas G. Polson is the Robert Law Jr. Professor of Econometrics and Statistics at the Booth School of Business, University of Chicago.",
        "Vadim Sokolov is on the faculty of the Department of Systems Engineering and Operations Research at George Mason University."]
bg.plot([lm, lm+tw*0.35], [yy+0.08]*2, color=gold, lw=0.8)
for b in bios:
    for line in textwrap.wrap(b, 62):
        bg.text(lm, yy, line, ha='left', va='top', color=cream, fontsize=11*K, fontproperties=serif); yy -= 0.21*K
    yy -= 0.1*K
by1 = BL + 0.25
print('lowest text y', round(yy, 3), 'barcode top', by1 + 1.2)
out = 'cover_vol%d_full' % VOL
fig.savefig(out + '.pdf', dpi=300); fig.savefig(out + '.png', dpi=150)
for xv in [BL, BL+TRIM_W, BL+TRIM_W+S, W-BL]: bg.axvline(xv, color='cyan', lw=0.8)
for yv in [BL, H-BL]: bg.axhline(yv, color='cyan', lw=0.8)
bg.add_patch(Rectangle((BL+TRIM_W-0.25-2.0, BL+0.25), 2.0, 1.2, fill=False, ec='white', lw=1, ls='--'))
bg.text(BL+TRIM_W-1.25, BL+0.85, 'barcode area', color='white', ha='center', va='center', fontsize=9)
fig.savefig(out + '_guides.png', dpi=80)
print('ok', out)
