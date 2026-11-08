"""
Class number five (D = 47 by default): the two dihedral modes of the single-class Linnik problem.
Classes C0 = (1,1,12), C1^{+-1} = (2,+-1,6), C1^{+-2} = (3,+-1,4); psi(C1) = e(1/5).
  c_1 = (1/2)(r0 + 2cos(2pi/5) r1 + 2cos(4pi/5) r2),  c_2 = (1/2)(r0 + 2cos(4pi/5) r1 + 2cos(8pi/5) r2)
(the real coefficient sequences of the two weight-one dihedral forms), S_j(N) = sum_{p<N} c_j(N-p), and
  S_C(N) - S_gen(N) = (2/5) sum_j 2Re psi_j(C) S_j(N).
Prints, on dyadic ranges of N:
  * Diag_j / (kappa_j pi) (kappa_j = mean of c_j^2) and the mean square E S_j^2 / (kappa_j pi)  (Conjecture 3.2: -> 1/2);
  * the cross diagonal sum_p c_1 c_2 (N-p) / pi(N) and the correlation of S_1, S_2 (orthogonality of the modes);
  * the mean square of the three class deviations against (4/25) sum_j (2Re psi_j(C))^2 E S_j^2 (no cross term).
Usage: python linnik_class_h5.py [D [NMAX]]   (D with h(-D) = 5 and the forms below: 47, 79)
"""
import os, sys, time
from math import cos, pi, sqrt
import numpy as np
from scipy.signal import fftconvolve
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from linnik_class import r_form, sieve   # noqa: E402

FORMS5 = {47: [(1, 1, 12), (2, 1, 6), (3, 1, 4)], 79: [(1, 1, 20), (2, 1, 10), (4, 1, 5)]}
D = int(sys.argv[1]) if len(sys.argv) > 1 else 47
NMAX = int(float(sys.argv[2])) if len(sys.argv) > 2 else 4 * 10 ** 6
t0 = time.time()
r = [r_form(a, b, c, NMAX).astype(float) for (a, b, c) in FORMS5[D]]
for x in r:
    x[0] = 0.0
# r[1] is the form of C1 (C1 and its inverse have the same r), r[2] that of C1^{+-2}; only cos(2 pi k/5) enters
w1, w2 = 2 * cos(2 * pi / 5), 2 * cos(4 * pi / 5)
c1 = 0.5 * (r[0] + w1 * r[1] + w2 * r[2])
c2 = 0.5 * (r[0] + w2 * r[1] + w1 * r[2])
isp = sieve(NMAX).astype(float)
piN = np.cumsum(isp)
S = [fftconvolve(isp, c)[:NMAX + 1] for c in (c1, c2)]
Dg = [fftconvolve(isp, c * c)[:NMAX + 1] for c in (c1, c2)]
D12 = fftconvolve(isp, c1 * c2)[:NMAX + 1]
kap = [np.mean(c[1:] ** 2) for c in (c1, c2)]
print(f"D = {D}, NMAX = {NMAX} [{time.time() - t0:.0f}s]; kappa_1, kappa_2 (mean c_j^2) = {kap[0]:.5f}, {kap[1]:.5f};"
      f" mean c_1 c_2 = {np.mean(c1[1:] * c2[1:]):+.5f}")
ranges = [(10 ** 5, 2 * 10 ** 5), (3 * 10 ** 5, 6 * 10 ** 5), (10 ** 6, 2 * 10 ** 6), (2 * 10 ** 6, 4 * 10 ** 6)]
print("range of N            Dg1/k1pi  ES1^2/k1pi   Dg2/k2pi  ES2^2/k2pi   Dg12/pi   corr(S1,S2)")
for lo, hi in ranges:
    if hi > NMAX:
        continue
    sl = slice(lo, hi)
    P = np.mean(piN[sl])
    v = [np.mean(s[sl] ** 2) for s in S]
    print(f"[{lo:.0e},{hi:.0e})   {np.mean(Dg[0][sl]) / (kap[0] * P):8.4f}  {v[0] / (kap[0] * P):9.4f}   {np.mean(Dg[1][sl]) / (kap[1] * P):8.4f}"
          f"  {v[1] / (kap[1] * P):9.4f}   {np.mean(D12[sl]) / P:+8.4f}   {np.mean(S[0][sl] * S[1][sl]) / sqrt(v[0] * v[1]):+8.4f}")
coef = {'C0': (2.0, 2.0), 'C1': (w1, w2), 'C2': (w2, w1)}
print("class deviations S_C - S_gen = (2/5)(a_1 S_1 + a_2 S_2):  E dev^2 against the uncorrelated value (4/25)(a_1^2 E S_1^2 + a_2^2 E S_2^2),")
print("   and E dev^2 / ((2/25)(a_1^2 kappa_1 + a_2^2 kappa_2) pi)   (= 1 under Conjecture 3.2 for both modes)")
for lo, hi in ranges:
    if hi > NMAX:
        continue
    sl = slice(lo, hi)
    P = np.mean(piN[sl])
    v = [np.mean(s[sl] ** 2) for s in S]
    row = []
    for name, (a1, a2) in coef.items():
        dev = 0.4 * (a1 * S[0][sl] + a2 * S[1][sl])
        m = np.mean(dev ** 2)
        row.append(f"{name}: {m / (0.16 * (a1 * a1 * v[0] + a2 * a2 * v[1])):.4f}, {m / (0.08 * (a1 * a1 * kap[0] + a2 * a2 * kap[1]) * P):.4f}")
    print(f"[{lo:.0e},{hi:.0e})   " + ";  ".join(row))
print(f"[{time.time() - t0:.0f}s]")
