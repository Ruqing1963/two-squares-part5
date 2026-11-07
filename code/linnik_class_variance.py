"""
Second moment of the class deviation S_c(N) = sum_{p<N} c(N-p) (D = 23 or 31): the diagonal Dg(N) = sum_{p<N} c(N-p)^2,
the full mean square, the Wilton-type sup of F(alpha) = sum_{m<=X} c(m) e(m alpha) over a grid, and the identity
sum_{N<=X} S_c(N)^2 <= int |S_p F|^2 <= sup|F|^2 pi(X) used for the almost-all theorem.
Usage: python linnik_class_variance.py [D [X]]
"""
import os, sys, time
from math import pi, sqrt, log
import numpy as np
from scipy.signal import fftconvolve
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from linnik_class import r_form, sieve, FORMS, kron   # noqa: E402

D = int(sys.argv[1]) if len(sys.argv) > 1 else 23
X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 4 * 10 ** 6
import linnik_class
linnik_class.D = D
t0 = time.time()
(a0, b0, c0), (a1, b1, c1) = FORMS[D]
r0, r1 = r_form(a0, b0, c0, X), r_form(a1, b1, c1, X)
r0[0] = r1[0] = 0.0
c = (r0 - r1) / 2
isp = sieve(X).astype(float)
Sc = np.rint(fftconvolve(isp, c)[:X + 1])
Dg = np.rint(fftconvolve(isp, c * c)[:X + 1])
piN = np.cumsum(isp)
theta = max(np.roots({23: [1, 0, -1, -1], 31: [1, -1, 0, -1]}[D]).real)
kap = 36 / (D + 1) * log(theta)
print(f"D = {D}, X = {X}, kappa_D = {kap:.6f} [{time.time() - t0:.0f}s]")
print(f"mean of c(m)^2 for m <= X: {np.mean(c[1:] ** 2):.5f} (-> kappa_D)")
for lo in (10 ** 5, 3 * 10 ** 5, 10 ** 6, 2 * 10 ** 6):
    if 2 * lo > X:
        continue
    sl = slice(lo, 2 * lo)
    V = np.mean(Sc[sl] ** 2)
    Dm = np.mean(Dg[sl])
    P = np.mean(piN[sl])
    print(f"N in [{lo},{2 * lo}): E S_c^2 = {V:.1f}, E Dg = {Dm:.1f}, E pi(N) = {P:.1f};  E S_c^2 / (kappa pi) = {V / (kap * P):.4f},"
          f" E Dg/(kappa pi) = {Dm / (kap * P):.4f}, E S_c^2 / E Dg = {V / Dm:.4f}")
# by residue of N mod D
sl = np.arange(X // 2, X)
for t in (1, -1, 0):
    m = np.array([kron(int(n)) == t for n in sl])
    print(f"   chi(N) = {t:+d}: E S_c^2/(kappa pi) = {np.mean(Sc[sl][m] ** 2) / (kap * np.mean(piN[sl][m])):.4f}, E Dg/(kappa pi) = {np.mean(Dg[sl][m]) / (kap * np.mean(piN[sl][m])):.4f}")
# Wilton sup over a fine grid of alpha (FFT of length 2^k >= X)
for Y in (10 ** 5, 10 ** 6, X):
    L = 1 << int(np.ceil(np.log2(8 * Y)))
    F = np.fft.fft(c[:Y + 1], n=L)
    print(f"sup_alpha |sum_(m<=Y) c(m) e(m alpha)| over {L} points, Y = {Y}: {np.max(np.abs(F)):.1f} = {np.max(np.abs(F)) / sqrt(Y):.3f} sqrt(Y) = {np.max(np.abs(F)) / (sqrt(Y) * log(Y)):.4f} sqrt(Y) log Y")
