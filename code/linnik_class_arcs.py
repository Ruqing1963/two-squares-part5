"""
The mechanism behind Var S_c(N) ~ (1/2) kappa_D pi(N): the exponential sum F(alpha) = sum_m c(m) w(m/X) e(m alpha) of the
weight-one dihedral form (smooth weight w(t) = t e^{-t}... here w(t) = e^{-t}) is small at a/q for q << sqrt(X):
the ratio rho(q) = mean_a |F(a/q)|^2 / (kappa_D sum w^2) as a function of q^2/X.  Then the predicted variance ratio
V/(kappa pi) = 1 - sum_q mu^2(q)/phi(q) (1 - rho(q)) / (log X + C), C = gamma + sum_p log p/(p(p-1)) = 1.3325,
from the mass mu^2(q)/phi(q) X/log^2 X of |S_p|^2 at the arcs a/q.
Usage: python linnik_class_arcs.py [D [X]]
"""
import os, sys, time
from math import pi, sqrt, log, gcd
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linnik_class
from linnik_class import r_form, FORMS

D = int(sys.argv[1]) if len(sys.argv) > 1 else 23
X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 6
linnik_class.D = D
M = 12 * X
(a0, b0, c0), (a1, b1, c1) = FORMS[D]
r0, r1 = r_form(a0, b0, c0, M), r_form(a1, b1, c1, M)
c = (r0 - r1) / 2
c[0] = 0.0
m = np.arange(M + 1, dtype=float)
w = np.exp(-m / X)
cw = c * w
theta = max(np.roots({23: [1, 0, -1, -1], 31: [1, -1, 0, -1]}[D]).real)
kap = 36 / (D + 1) * log(theta)
norm = kap * np.sum(w[1:] ** 2)
print(f"D = {D}, X = {X}; sum |c w|^2 / (kappa sum w^2) = {np.sum(cw ** 2) / norm:.4f}")


def phi(q):
    r, n, p = q, q, 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            r -= r // p
        p += 1
    if n > 1:
        r -= r // n
    return r


def sqfree(q):
    p = 2
    while p * p <= q:
        if q % (p * p) == 0:
            return False
        p += 1
    return True


c2 = c[1:] ** 2
nn = np.arange(1, M + 1, dtype=float)


def Phi(Y):
    """scaling profile Phi_D(Y) = (1/(kappa Y)) sum_n c(n)^2 e^{-n/Y} (Fricke expansion at the cusp a/q, averaged over a)"""
    return np.sum(c2 * np.exp(-nn / Y)) / (kap * Y)


rows = []
qs = sorted(set([int(round(v)) for v in np.geomspace(1, 6 * sqrt(X), 60)] + [D * k for k in range(1, int(6 * sqrt(X) / D) + 1, 3)]))
print("   Phi: Y = D q^2/(8 pi^2 X) for D not| q, Y = q^2/(8 pi^2 X) for D | q (cusp infinity, width 1)")
for q in qs:
    if not sqfree(q):
        continue
    # sum over m of c w e(a m/q) for all a at once: group m by residue mod q
    res = np.bincount((np.arange(M + 1) % q), weights=cw, minlength=q)
    F = np.fft.fft(res)            # F[a] = sum_r res[r] e(-a r/q)
    a_ok = [a for a in range(q) if gcd(a, q) == 1]
    rho = np.mean(np.abs(F[a_ok]) ** 2) / norm
    Y = (q * q if q % D == 0 else D * q * q) / (8 * pi ** 2 * X)
    ph = Phi(Y)
    if q % D != 0:
        rows.append((q, rho))
    print(f"   q = {q:6d}{' (D|q)' if q % D == 0 else '      '}, q^2/X = {q * q / X:8.4f}, Y = {Y:8.4f}: rho(q) = {rho:.4f},"
          f"  Phi(Y) = {ph:.4f}")
# the maximum of Phi and the Mellin identity  int_0^inf (1_{Y>1} - Phi(Y)) dY/Y = gamma - c_0/kappa,
# c_0 the constant term of sum c(n)^2 n^{-s} = kappa/(s-1) + c_0 + ..., i.e. c_0 = lim (sum_{n<=x} c(n)^2/n - kappa log x)
Ys = np.geomspace(0.05, 200, 400)
ph = np.array([Phi(Y) for Y in Ys])
i = np.argmax(ph)
print(f"max Phi = {ph[i]:.4f} at Y = {Ys[i]:.3f} (q^2/X = {Ys[i] * 8 * pi ** 2 / D:.2f} for D not| q)")
cs = np.cumsum(c2 / nn)
xs = np.array([M // 8, M // 4, M // 2, M])
c0s = cs[xs - 1] - kap * np.log(xs)
print("c_0 estimates  sum_{n<=x} c^2/n - kappa log x at x = " + ", ".join(f"{x}: {v:.5f}" for x, v in zip(xs, c0s)))
c0 = c0s[-1]
# J = int_0^1 (1 - Phi) dY/Y ... in t = 1/Y: J = int_0^1 (1 - t G/kappa) dt/t - int_1^inf G/kappa dt, G(t) = sum c^2 e^{-nt}
ts = np.geomspace(1.0 / (M / 40), 60, 3000)          # t >= 40/M: e^{-n t} negligible beyond n = M
G = np.array([np.sum(c2 * np.exp(-nn * t)) for t in ts])
lt = np.log(ts)
f1 = np.where(ts < 1, 1 - ts * G / kap, 0.0)
f2 = np.where(ts >= 1, ts * G / kap, 0.0)
J = np.trapz(f1, lt) - np.trapz(f2, lt)
gam = 0.5772156649015329
print(f"Mellin identity: int (1_(Y>1) - Phi) dY/Y = {J:.4f}   gamma - c_0/kappa = {gam - c0 / kap:.4f}   (kappa = {kap:.6f}, c_0 = {c0:.5f})")
# predicted variance ratio, interpolating rho in q^2/X
qq = np.array([r[0] for r in rows], float)
rr = np.array([r[1] for r in rows])
Ls = log(X)
tot = 0.0
QMAX = int(40 * sqrt(X))
for q in range(1, QMAX + 1):
    if not sqfree(q):
        continue
    rho = np.interp(q * q / X, qq ** 2 / X, rr, right=1.0)
    tot += (1 - rho) / phi(q)
pred = 1 - tot / (Ls + 1.3325)
print(f"predicted E S_c^2 / (kappa pi(N)) at N ~ X: 1 - {tot:.3f}/(log X + 1.3325) = {pred:.4f}  (log X = {Ls:.2f})")
