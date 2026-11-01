import mpmath as mp, json
mp.mp.dps = 25
def xi(s):
    if abs(s-1) < mp.mpf("1e-20") or abs(s) < mp.mpf("1e-20"): return mp.mpf(1)/2
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def dlogxi(s): return 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)
def F(w, t): return mp.re(xi(mp.mpf(1)/2 + w + 1j*t))
T = 60
# zero trajectories
traj = {}
for i in range(0, 21):
    w = mp.mpf(i)/40; zs = []
    h = mp.mpf('0.02'); t = mp.mpf('0.01'); f0 = F(w, t) if w > 0 else mp.re(xi(mp.mpf(1)/2+1j*t))
    G = (lambda x, w=w: F(w, x))
    while t < T:
        t1 = t + h; f1 = G(t1)
        if f0*f1 < 0: zs.append(float(mp.findroot(G, (t, t1), solver='anderson')))
        t, f0 = t1, f1
    traj[float(w)] = zs
    print(float(w), len(zs))
# characteristic functions phi_w(t) normalised
cf = {}
for w in ['0', '0.25', '0.5']:
    wm = mp.mpf(w); n0 = mp.re(xi(mp.mpf(1)/2 + wm))
    cf[w] = [(float(t), float(F(wm, t)/n0)) for t in mp.linspace(0, T, 601)]
# Thorin density at level y of p_w: u = (1/pi) Re L_w'/L_w(y+iv),  L_w(s) ∝ xi(1/2+s+w)+xi(1/2+s-w)
def u_w(w, y, v):
    s = y + 1j*v
    a = mp.mpf(1)/2 + s + w; b = mp.mpf(1)/2 + s - w
    xa, xb = xi(a), xi(b)
    return mp.re((xa*dlogxi(a) + xb*dlogxi(b))/(xa + xb))/mp.pi
th = {}
for w in ['0', '0.5']:
    wm = mp.mpf(w); rows = []
    for v in mp.linspace(mp.mpf('0.001'), T, 600):
        rows.append((float(v), float(u_w(wm, mp.mpf('0.3'), v))))
    th[w] = rows
json.dump(dict(traj=traj, cf=cf, th=th), open('figs_data.json', 'w'))
