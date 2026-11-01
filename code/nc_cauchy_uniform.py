# Problem 2 on exact laws: X ~ q U[-1,1] + (1-q) U[-3,3]. phi(z) = (sin z/z)(1 - (4/3)(1-q) sin^2 z).
# Cosine tilt: psi_w(z) = phi(z+iw)+phi(z-iw). Count zeros in a rectangle by the argument principle,
# real zeros by sign changes; nonreal = difference. Good set G(q) = {w : no nonreal zeros}.
import numpy as np, json
def phi(z, q): 
    s = np.sin(z); return np.where(np.abs(z) < 1e-12, 1.0, s/np.where(np.abs(z)<1e-12,1,z))*(1 - (4/3)*(1-q)*s*s)
def psi(z, q, w): return phi(z+1j*w, q) + phi(z-1j*w, q)
def count(q, w, X=40*np.pi+0.37, Y=4.0, n=200000):
    # rectangle [-X,X]x[-Y,Y]
    t = np.linspace(0, 1, n, endpoint=False)
    pts = np.concatenate([-X+2*X*t - 1j*Y, X+1j*(-Y+2*Y*t), X-2*X*t + 1j*Y, -X+1j*(Y-2*Y*t)])
    v = psi(pts, q, w); ang = np.unwrap(np.angle(v))
    total = int(round((ang[-1]-ang[0] + (np.angle(v[0]) - np.angle(v[-1]) + np.pi) % (2*np.pi) - np.pi)/(2*np.pi)))
    x = np.linspace(-X, X, 400001); f = psi(x.astype(complex), q, w).real
    real = int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
    return total, real
out = {}
for q in [0.3, 0.5, 0.8]:
    A = 3/(4*(1-q)); delta = float(np.arccosh(np.sqrt(A)))
    rows = []
    for w in np.round(np.arange(0.0, 1.01, 0.02), 3):
        tot, re = count(q, w); rows.append((float(w), tot, re, tot-re))
    good = [r[0] for r in rows if r[3] == 0]
    out[q] = dict(delta=delta, rows=rows)
    print(f"q={q} delta={delta:.4f}  good w: from {min(good) if good else None}; nonreal counts:", [(r[0], r[3]) for r in rows if r[0] in (0.0,0.2,0.4,0.5,0.6,0.7,0.8,1.0)])
    # find threshold and interval check
    bad_after = [r[0] for r in rows if r[3] > 0 and good and r[0] > min(good)]
    print("   bad w above first good w:", bad_after)
json.dump(out, open("toy_uniform.json", "w"), indent=0)
