"""
The Linnik problem for a single ideal class: N = p + Q_C(x, y), D = 23 (h = 3) and D = 31.

r_C = r_gen + (w/h) sum_psi psi-bar(C) c_psi, so for h = 3
    S_C(N) = sum_{p<N} r_C(N-p) = S_gen(N) + alpha_C S_c(N),   S_gen = sum_p r_gen(N-p),  S_c = sum_p c(N-p),
alpha_0 = 4/3, alpha_1 = -2/3, c = (r_0 - r_1)/2 the coefficients of the weight-one dihedral cusp form.
Prediction: S_gen(N) = (2 pi/sqrt D) SS_D(N) pi(N) (1 + o(1)) with the singular series of linnik_heegner.tex (the genus
form has the local densities of the principal form), and S_c(N) = o(pi(N)), so every class gets the same main term.
The script computes S_0, S_1, S_gen, S_c for N in ranges up to NMAX by FFT (prime indicator * r tables), the ratios
S_C/((2pi/sqrt D) SS pi(N)), and the size of S_c relative to sqrt(N), N/log N, and fits |S_c| ~ N^theta.
Usage: python linnik_class.py [D [NMAX]]
Output: results/linnik_class_{D}.txt
"""
import os
import sys
import time
from math import pi, sqrt, log
import numpy as np
from scipy.signal import fftconvolve

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
os.makedirs(RES, exist_ok=True)
D = int(sys.argv[1]) if len(sys.argv) > 1 else 23
NMAX = int(float(sys.argv[2])) if len(sys.argv) > 2 else 4 * 10 ** 6
FORMS = {23: [(1, 1, 6), (2, 1, 3)], 31: [(1, 1, 8), (2, 1, 4)]}


def r_form(a, b, c, X):
    Dq = 4 * a * c - b * b
    ymax = int(sqrt(4 * a * X / Dq)) + 1
    r = np.zeros(X + 1, dtype=np.int64)
    for y in range(-ymax, ymax + 1):
        disc = b * b * y * y - 4 * a * (c * y * y - X)
        if disc < 0:
            continue
        s = sqrt(disc)
        x = np.arange(int(np.floor((-b * y - s) / (2 * a))) - 1, int(np.ceil((-b * y + s) / (2 * a))) + 2, dtype=np.int64)
        k = a * x * x + b * x * y + c * y * y
        k = k[(k >= 0) & (k <= X)]
        r += np.bincount(k, minlength=X + 1)
    return r.astype(float)


def kron(n):
    """Kronecker symbol (-D/n) for prime D = 3 mod 4: equals the Legendre symbol (n/D)"""
    n = int(n) % D
    if n == 0:
        return 0
    return 1 if pow(n, (D - 1) // 2, D) == 1 else -1


def sieve(X):
    s = np.ones(X + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def singular_series(N, primes_small, chi_small):
    """SS_D(N) = (1 + c_D(N)) prod_{l | N} (1 - chi(l)/l) prod_{l not| N} (1 + chi(l)/(l(l-1))), c_D(N) = -chi(N)/(D-1);
    the product over l not| N is truncated at the table and corrected by the tail estimate (tiny)."""
    N = int(N)
    val = 1.0 - kron(N) / (D - 1)
    for l, c in zip(primes_small, chi_small):
        l = int(l)
        if N % l == 0:
            val *= 1 - c / l
        else:
            val *= 1 + c / (l * (l - 1))
    # primes dividing N beyond the table
    m = N
    for l in primes_small:
        l = int(l)
        while m % l == 0:
            m //= l
    if m > 1:   # m is a product of large primes (at most a few); factor by trial division
        q = int(primes_small[-1]) + 2
        while q * q <= m:
            if m % q == 0:
                cq = kron(q)
                val *= (1 - cq / q) / (1 + cq / (q * (q - 1)))
                while m % q == 0:
                    m //= q
            q += 2
        if m > 1:
            cq = kron(m)
            val *= (1 - cq / m) / (1 + cq / (m * (m - 1)))
    return val


def main():
    out = []

    def say(s=""):
        out.append(s)
        print(s, flush=True)
    t0 = time.time()
    (a0, b0, c0), (a1, b1, c1) = FORMS[D]
    r0 = r_form(a0, b0, c0, NMAX)
    r1 = r_form(a1, b1, c1, NMAX)
    r0[0] = r1[0] = 0.0                      # N - p >= 1
    rg = (r0 + 2 * r1) / 3
    cc = (r0 - r1) / 2
    isp = sieve(NMAX).astype(float)
    say(f"D = {D}, h = 3: forms {FORMS[D]}, tables to {NMAX} [{time.time() - t0:.0f}s]")
    # S(N) = sum_p r(N - p) = (isp * r)(N)
    S0 = np.rint(fftconvolve(isp, r0)[:NMAX + 1])
    S1 = np.rint(fftconvolve(isp, r1)[:NMAX + 1])
    Sg = (S0 + 2 * S1) / 3
    Sc = (S0 - S1) / 2
    pi_N = np.cumsum(isp)
    say(f"convolutions done [{time.time() - t0:.0f}s]")
    ps = np.nonzero(sieve(5000))[0]
    chis = [kron(int(l)) for l in ps]
    const = 2 * pi / sqrt(D)
    rows = []
    rng = np.random.default_rng(1)
    for lo in [10 ** 4, 3 * 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6, 2 * 10 ** 6, NMAX // 2]:
        if lo >= NMAX:
            continue
        Ns = rng.integers(lo, min(2 * lo, NMAX), size=400)
        rat0, rat1, ratg, scn = [], [], [], []
        for N in Ns:
            N = int(N)
            main = const * singular_series(N, ps, chis) * pi_N[N - 1]
            rat0.append(S0[N] / main)
            rat1.append(S1[N] / main)
            ratg.append(Sg[N] / main)
            scn.append(Sc[N])
        rat0, rat1, ratg, scn = map(np.array, (rat0, rat1, ratg, scn))
        rows.append((lo, np.sqrt(np.mean(scn ** 2))))
        say(f"N in [{lo}, {2 * lo}): S_0/main = {rat0.mean():.4f} (sd {rat0.std():.4f}), S_1/main = {rat1.mean():.4f} (sd {rat1.std():.4f}),"
            f" S_gen/main = {ratg.mean():.4f} (sd {ratg.std():.4f});  rms S_c = {np.sqrt(np.mean(scn ** 2)):.1f},"
            f" rms S_c/sqrt(N) = {np.sqrt(np.mean(scn ** 2)) / sqrt(1.5 * lo):.3f}, rms S_c/(N/log N) = {np.sqrt(np.mean(scn ** 2)) / (1.5 * lo / log(1.5 * lo)):.5f},"
            f" mean S_c/sqrt N = {np.mean(scn) / sqrt(1.5 * lo):+.3f}")
    L = np.log([1.5 * r[0] for r in rows])
    y = np.log([r[1] for r in rows])
    th = np.polyfit(L, y, 1)[0]
    say(f"\nfit rms |S_c(N)| ~ N^theta: theta = {th:.3f}  (main term ~ N/log N)")
    # dependence of S_c on N mod 23 and on the class of N
    sel = np.arange(NMAX // 2, NMAX)
    say(f"\nmean of S_c(N)/sqrt(N) over N in [{NMAX // 2}, {NMAX}) by the residue of N mod {D} (quadratic character):")
    for t in (1, -1, 0):
        m = np.array([kron(int(N)) == t for N in sel[::7]])
        vals = Sc[sel[::7]][m] / np.sqrt(sel[::7][m])
        say(f"   chi(N) = {t:+d}: mean {vals.mean():+.3f}, sd {vals.std():.3f}, n = {m.sum()}")
    path = os.path.join(RES, f'linnik_class_{D}.txt')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write("\n".join(out) + "\n")
    print("written", path)


if __name__ == "__main__":
    main()
