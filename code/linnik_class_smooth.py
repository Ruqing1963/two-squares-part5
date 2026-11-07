"""
Smooth variance of the class deviation: V(X) = sum_N e^{-2N/X} S_c(N)^2 = int_0^1 |P(alpha)|^2 |F(alpha)|^2 d alpha with
P = sum_p e^{-p/X} e(p alpha), F = sum_m c(m) e^{-m/X} e(m alpha), against the diagonal
Diag(X) = sum_N e^{-2N/X} sum_{p<N} c(N-p)^2, and against the major-arc model
    V/Diag ~ 1 - sum_q mu^2(q)/phi(q) (1 - rho_beta(q)) / sum_q mu^2(q)/phi(q) [over q <= X],
where rho_beta(q) = mean over a of int |T(beta)|^2 |F(a/q+beta)|^2 d beta / (int |T|^2 . mean|F|^2), T = sum e^{-n/X} e(n beta)/log n,
computed exactly on a grid (the arcs of |P|^2 have width ~ 1/X).
Usage: python linnik_class_smooth.py [D [X ...]]
"""
import os, sys, time
from math import pi, sqrt, log, gcd
import numpy as np
from scipy.signal import fftconvolve
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linnik_class
from linnik_class import r_form, sieve, FORMS

D = int(sys.argv[1]) if len(sys.argv) > 1 else 23
XS = [int(float(a)) for a in sys.argv[2:]] or [10 ** 5, 2 * 10 ** 5, 4 * 10 ** 5]
linnik_class.D = D
M = 12 * max(XS)
(a0, b0, c0), (a1, b1, c1) = FORMS[D]
r0, r1 = r_form(a0, b0, c0, M), r_form(a1, b1, c1, M)
c = (r0 - r1) / 2
c[0] = 0.0
isp = sieve(M).astype(float)
Sc = np.rint(fftconvolve(isp, c)[:M + 1])
Dg = np.rint(fftconvolve(isp, c * c)[:M + 1])
n = np.arange(M + 1, dtype=float)
theta = max(np.roots({23: [1, 0, -1, -1], 31: [1, -1, 0, -1]}[D]).real)
kap = 36 / (D + 1) * log(theta)


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


for X in XS:
    wN = np.exp(-2 * n / X)
    V = np.sum(wN * Sc ** 2)
    Dm = np.sum(wN * Dg)
    # model: F at a/q + beta on a beta-grid of width 6/X; T(beta) for the prime weights e^{-n/X}/log n (n >= 2)
    cw = c * np.exp(-n / X)
    tw = np.zeros(M + 1)
    tw[2:] = np.exp(-n[2:] / X) / np.log(n[2:])
    betas = np.linspace(-3.0 / X, 3.0 / X, 13)
    Tb = np.array([np.sum(tw[:10 * X] * np.exp(2j * pi * n[:10 * X] * b)) for b in betas])
    wT = np.abs(Tb) ** 2
    meanF2 = np.sum(cw ** 2)                 # int_0^1 |F|^2
    num = den = 0.0
    qs = [q for q in range(1, int(4 * sqrt(X)) + 1) if mu2(q)]
    rhos = {}
    for q in qs:
        # F(a/q + beta) for all a and the beta grid: sum over residues of the twisted partial sums
        idx = np.arange(10 * X) % q
        vals = []
        for b in betas:
            res = np.bincount(idx, weights=cw[:10 * X] * np.cos(2 * pi * n[:10 * X] * b), minlength=q) \
                + 1j * np.bincount(idx, weights=cw[:10 * X] * np.sin(2 * pi * n[:10 * X] * b), minlength=q)
            Fa = np.fft.ifft(res) * q                   # F[a] = sum_r res[r] e(a r/q)
            a_ok = [a for a in range(q) if gcd(a, q) == 1]
            vals.append(np.mean(np.abs(Fa[a_ok]) ** 2))
        rho = np.sum(wT * np.array(vals)) / (np.sum(wT) * meanF2)
        rhos[q] = rho
    tot = sum(1 / phi(q) for q in qs) + 0.0
    # mass of |P|^2 on arcs q > 8 sqrt(X): rho = 1 there; total mass ~ sum_{q <= X} mu^2/phi
    full = sum(1.0 / phi(q) for q in range(1, X + 1) if mu2(q)) if X <= 2 * 10 ** 5 else log(X) + 1.3325
    defic = sum((1 - rhos[q]) / phi(q) for q in qs)
    pred = 1 - defic / full
    print(f"D = {D}, X = {X}: V/Diag = {V / Dm:.4f}  (Diag/(kappa sum e^-2N/X pi(N)) = {Dm / (kap * np.sum(wN * np.cumsum(isp))):.4f});"
          f"  model 1 - {defic:.3f}/{full:.3f} = {pred:.4f};  rho at q^2/X = 0.1, 1, 4: "
          f"{np.interp(0.1, [q * q / X for q in qs], [rhos[q] for q in qs]):.3f}, {np.interp(1, [q * q / X for q in qs], [rhos[q] for q in qs]):.3f}, {np.interp(4, [q * q / X for q in qs], [rhos[q] for q in qs]):.3f}", flush=True)
