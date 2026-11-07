"""
Split of the smoothed variance V = int |P|^2 |F|^2 (P = sum_p e^{-p/X} e(p a), F = sum c(m) e^{-m/X} e(m a)) into the
major arcs M(Q, K) = {|alpha - a/q| <= K/X, q <= Q} and the rest m:
  * mass fraction of |P|^2 on M, against the Hardy-Littlewood value sum_{q<=Q} mu^2/phi(q) int_{|b|<=K/X}|T|^2 / int|P|^2;
  * the share of V on M (predicted ~ 0 for Q << sqrt X);
  * the decorrelation ratio on m:  int_m |P|^2|F|^2 / (mean|F|^2 . int_m |P|^2),  with mean|F|^2 = int_0^1 |F|^2.
Conjecture 3.2 <=> this ratio -> 1 with M of order X^{1/2}.
Usage: python linnik_class_split.py [D [X]]
"""
import os, sys, time
from math import pi, sqrt, log, gcd
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linnik_class
from linnik_class import r_form, sieve, FORMS

D = int(sys.argv[1]) if len(sys.argv) > 1 else 23
X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 2 * 10 ** 5
linnik_class.D = D
M = 20 * X
L = 1 << int(np.ceil(np.log2(40 * X)))
t0 = time.time()
(a0, b0, c0), (a1, b1, c1) = FORMS[D]
r0, r1 = r_form(a0, b0, c0, M), r_form(a1, b1, c1, M)
c = (r0 - r1) / 2
c[0] = 0.0
n = np.arange(M + 1, dtype=float)
w = np.exp(-n / X)
isp = sieve(M).astype(float)
P = np.fft.fft(isp * w, n=L)
F = np.fft.fft(c * w, n=L)
P2 = np.abs(P) ** 2
F2 = np.abs(F) ** 2
V = np.mean(P2 * F2)
meanF2 = np.mean(F2)
totP = np.mean(P2)
print(f"D = {D}, X = {X}, grid L = {L} (step = 1/{L / X:.0f}X) [{time.time() - t0:.0f}s]")
print(f"V = {V:.4e}, int|P|^2 = {totP:.4e} (sum_p e^-2p/X = {np.sum(isp * w * w):.4e}), mean|F|^2 = {meanF2:.4e}; V/(mean|F|^2 int|P|^2) = {V / (meanF2 * totP):.4f}")
# T(beta) with weights e^{-n/X}/log n for the HL arc mass
tw = np.zeros(M + 1)
tw[2:] = w[2:] / np.log(n[2:])
T2 = np.abs(np.fft.fft(tw, n=L)) ** 2


def phi(q):
    r, m, p = q, q, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def mu2(q):
    p = 2
    while p * p <= q:
        if q % (p * p) == 0:
            return 0
        p += 1
    return 1


KS = [float(v) for v in os.environ.get('KS', '1,3').split(',')]
CQS = [float(v) for v in os.environ.get('CQS', '0.25,0.5,1,2').split(',')]
for K in KS:
    half = int(K / X * L)
    T2arc = np.sum(T2[:half + 1]) + np.sum(T2[L - half:])
    for cQ in CQS:
        Q = int(cQ * sqrt(X))
        mask = np.zeros(L, dtype=bool)
        hl = 0.0
        for q in range(1, Q + 1):
            a = np.array([a for a in range(q) if gcd(a, q) == 1], dtype=np.int64)
            centres = (a * L + q // 2) // q          # int64: a * L overflows int32 on Windows
            offs = np.arange(-half, half + 1)
            idx = (centres[:, None] + offs[None, :]) % L
            mask[idx.ravel()] = True
            if mu2(q):
                hl += 1.0 / phi(q)
        massM = np.sum(P2[mask]) / L
        VM = np.sum(P2[mask] * F2[mask]) / L
        massm = totP - massM
        Vm = V - VM
        print(f"K = {K}, Q = {cQ} sqrt X = {Q}: |P|^2 mass on M = {massM / totP:.4f} (HL {hl * (T2arc / L) / totP:.4f}); share of V on M = {VM / V:.4f};"
              f"  decorrelation on m: int_m|P F|^2/(mean|F|^2 int_m|P|^2) = {Vm / (meanF2 * massm):.4f};  measure of M = {mask.mean():.4f}", flush=True)
