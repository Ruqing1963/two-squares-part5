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


rows = []
qs = sorted(set([int(round(v)) for v in np.geomspace(1, 6 * sqrt(X), 60)]))
for q in qs:
    if not sqfree(q):
        continue
    # sum over m of c w e(a m/q) for all a at once: group m by residue mod q
    res = np.bincount((np.arange(M + 1) % q), weights=cw, minlength=q)
    F = np.fft.fft(res)            # F[a] = sum_r res[r] e(-a r/q)
    a_ok = [a for a in range(q) if gcd(a, q) == 1]
    rho = np.mean(np.abs(F[a_ok]) ** 2) / norm
    rows.append((q, rho))
    print(f"   q = {q:6d}, q^2/X = {q * q / X:8.4f}: rho(q) = {rho:.4f}")
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
