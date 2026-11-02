"""Track E, splines, second script. Knots -q, -(q-k), q-k, q (integers). With C = cosh z, s = sinh z:
k = 1: f_q(z) = sinh(qz)/q - sinh((q-1)z)/(q-1) = s (C-1) Q_q(C), deg Q_q = q-2, so the B spline is
       U_1 + U_{1/2} + U_{1/2} + S_q with E e^{i theta S_q} = Q_q(cos theta)/Q_q(1) (if the Fourier coefficients are >= 0);
k = 2, q odd: f = s^3 Q(C), deg Q = q-3, B spline = 3 U_1 + S;
k = 3: f = s (C-1) Q(C), deg Q = q-2.
For each: zeros of f on the imaginary axis in [0.5,120] against all zeros in [0.5,120]x[-8,8] (argument principle),
and the Fourier coefficients of Q(cos theta) by FFT (minimum coefficient), so that S is a probability law or not."""
import numpy as np

def winding(fun, x0, x1, y0, y1, m=80000):
    xs = np.linspace(x0, x1, m); ys = np.linspace(y0, y1, m//4)
    path = np.concatenate([xs + 1j*y0, x1 + 1j*ys, xs[::-1] + 1j*y1, x0 + 1j*ys[::-1]])
    v = fun(path); d = np.angle(v[1:]/v[:-1])
    return -1 if np.max(np.abs(d)) > 2.5 else int(round(d.sum()/(2*np.pi)))

def report(q, k):
    a, b = q, q - k
    f = lambda th: np.sin(a*th)/a - np.sin(b*th)/b        # theta is the variable of the integer knots
    nfull = winding(f, 0.5, 120.0, -8.0, 8.0); nthin = winding(f, 0.5, 120.0, -0.002, 0.002)
    # Q(cos theta): divide by the lattice factor on a grid avoiding its zeros, then FFT
    M = 4096
    th = (2*np.pi*(np.arange(M) + 0.5))/M
    if k == 2 and q % 2 == 1:
        den = np.sin(th)**3; deg = q - 3
    else:
        den = np.sin(th)*(np.cos(th) - 1); deg = q - 2
    Q = f(th)/den                                            # a trigonometric polynomial of degree deg in theta
    coef = np.fft.fft(Q*np.exp(1j*th*0))/M                   # coefficients on the half shifted grid
    # undo the half shift: coefficient of e^{i m theta} picks up e^{i pi m/M}
    m = np.fft.fftfreq(M, 1.0/M)
    c = (coef*np.exp(1j*np.pi*m/M)).real
    Qcos = c[:deg+1]; Qcos[1:] *= 2                           # cosine coefficients: Q = c0 + sum_{m>=1} 2 c_m cos(m theta)
    Qcos /= Qcos.sum()                                       # normalise Q(1) = 1
    beyond = np.max(np.abs(c[deg+1:M//2])) if deg + 1 < M//2 else 0.0
    print("  q=%2d k=%d  r=%.4f  zeros full=%3d thin=%3d %s   S law: min Fourier coef = %+.4e (max beyond degree %.1e)  %s"
          % (q, k, b/a, nfull, nthin, "ALL REAL" if nfull == nthin and nfull >= 0 else "nonreal ", Qcos.min(), beyond, "probability law" if Qcos.min() > -1e-12 else "signed"))
    return Qcos

print("family k = 1, knots -q, -(q-1), q-1, q")
for q in list(range(2, 13)) + [16, 20, 30, 40]:
    Qc = report(q, 1)
    if q <= 5: print("     coefficients of S (weights of 0, +-1, +-2, ...):", np.round(Qc, 6))
print("family k = 2, q odd, knots -q, -(q-2), q-2, q")
for q in [3, 5, 7, 9, 11, 13, 15, 21, 31]:
    Qc = report(q, 2)
    if q <= 7: print("     coefficients of S:", np.round(Qc, 6))
print("family k = 3, knots -q, -(q-3), q-3, q")
for q in [4, 5, 7, 8, 10, 11, 13, 20]:
    report(q, 3)
print("family k = 4")
for q in [5, 6, 7, 9, 11, 13]:
    report(q, 4)
