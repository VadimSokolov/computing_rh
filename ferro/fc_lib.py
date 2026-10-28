# fc_lib.py: shared routines of group FC (Knauf's arithmetic ferromagnet): counts c_N(m) from the births saved by fc_stern.py,
# the depth profile w_N = c_N / phi, Z_N(beta) = sum c_N(m) m^(-beta), the Mellin transform of the profile
# W_N(z) = int_0^inf w_N(x) x^(z-1) dx = (1/z) sum_{m >= 2} (w_N(m-1) - w_N(m)) m^z   (w_N(1) = 1, w_N(m) = 0 for m > F_(N+2); summation by
# parts of 1/z + sum_{m >= 1} w_N(m) ((m+1)^z - m^z)/z), the Perron remainder R_N and the profile corrected chain C_N = Z_N - W_N(2-beta)/zeta(2).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import numpy as np

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def phi_sieve(n):
    phi = np.arange(n + 1, dtype=np.int64)
    is_p = np.ones(n + 1, dtype=bool); is_p[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if is_p[p]:
            is_p[p * p::p] = False
    for p in np.nonzero(is_p)[0]:
        phi[p::p] -= phi[p::p] // p
    return phi

class Counts:
    """c_N(m) for N = 1..nmax from the births file."""
    def __init__(self, births_path, nmax=None):
        self.z = np.load(births_path)
        self.nmax = max(int(k[1:]) for k in self.z.files) if nmax is None else nmax
        self.phi = phi_sieve(fib(self.nmax + 2) + 2)
        self._c = np.zeros(fib(self.nmax + 2) + 2, dtype=np.int64)
        self._c[1] = 1
        self._N = 0
    def c(self, N):
        """array c_N(0..F_(N+2)) as int64 (N must be requested in increasing order or the first time)"""
        assert N >= self._N or N == 0
        while self._N < N:
            self._N += 1
            h = self.z['H%d' % self._N]
            self._c[:len(h)] += h.astype(np.int64)
        return self._c[:fib(N + 2) + 1].copy()
    def c_dense(self, N):
        """recompute from scratch (any order)"""
        c = np.zeros(fib(N + 2) + 1, dtype=np.int64); c[1] = 1
        for j in range(1, N + 1):
            h = self.z['H%d' % j].astype(np.int64); c[:len(h)] += h
        return c

def profile_terms(c, phi):
    """returns m (float), cm (float), d2 (float) with d2(m) = (w(m-1) - w(m)) m^2 for m = 2..F+1; w(F+1) = 0"""
    Mn = len(c) - 1
    w = np.zeros(Mn + 2)
    w[1:Mn + 1] = c[1:Mn + 1] / phi[1:Mn + 1]
    m = np.arange(2, Mn + 2, dtype=np.float64)
    d = w[1:Mn + 1] - w[2:Mn + 2]
    return d, w

def Z_and_W(c, phi, beta_list, with_deriv=False):
    """double precision values of A = sum c m^-beta, B = sum d m^(2-beta) (so W(2-beta) = B/(2-beta)) at each beta; optionally derivatives
    in beta (A', B'). Returns arrays."""
    Mn = len(c) - 1
    d, w = profile_terms(c, phi)
    m1 = np.arange(1, Mn + 1, dtype=np.float64)
    lm1 = np.log(m1)
    cm = c[1:Mn + 1].astype(np.float64)
    m2 = np.arange(2, Mn + 2, dtype=np.float64)
    lm2 = np.log(m2)
    d2 = d * m2 * m2
    out = []
    for b in beta_list:
        e1 = np.exp(-b * lm1)
        e2 = np.exp(-b * lm2)
        A = np.dot(cm, e1)
        B = np.dot(d2, e2)
        if with_deriv:
            Ap = -np.dot(cm * lm1, e1)
            Bp = -np.dot(d2 * lm2, e2)
            out.append((A, B, Ap, Bp))
        else:
            out.append((A, B))
    return out

def zeros_list(path, n):
    zs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                zs.append(line)
            if len(zs) >= n:
                break
    return zs
