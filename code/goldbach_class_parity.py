"""
Sieve weights and Omega-strata of the class deviation (Section 3.3 of notes/linnik_class.tex), D = 23 (h = 3), D = 31.

For a weight u(m) on m = N - p and the two classes C_0 (principal form), C_1 of D (h = 3) we compute, for all N <= NMAX,
    S_C^u(N) = sum_{p<N} u(N-p) r_C(N-p),   S_c^u(N) = sum_{p<N} u(N-p) c(N-p),   Diag^u(N) = sum_{p<N} u(N-p) c(N-p)^2,
c = (r_0 - r_1)/2 the coefficients of the weight-one dihedral cusp form, by FFT (prime indicator * weighted tables).
Weights: u = 1 (the unsifted deviation of the note), u = 1_{spf(m) >= z} with z = lo^{1/beta}, beta = 3, 4, 6, lo the start
of the range of N (the sharp sieve by the primes below z; beta = 3 leaves N - p with Omega <= 2 up to a few Omega = 3,
the Chen set), u = 1_{Omega(m) = k} for k = 1 (Goldbach pairs N = p + q with q = Q_C(x, y) prime), 2, 3, and the sifted
strata 1_{spf >= z} 1_{Omega = k} for beta = 3.
Reported per range: the class shares S_{C_0}^u / (S_{C_0}^u + 2 S_{C_1}^u) and S_{C_1}^u / (...) (equidistribution: 1/3 each),
the bias mean S_c^u / rms S_c^u, the mean diagonal against pi(N), and the variance ratio E (S_c^u)^2 / E Diag^u, which is
0.42-0.43 (conjecturally -> 1/2) for u = 1: the question is whether the depletion survives the sieve.
Also the distribution of c(m) on the sifted strata: c(m) = (-1)^Omega(m) on the m whose prime factors are all split and
non-principal, but the mean of c over every stratum is 0 (Chebotarev and multiplicativity): no class-group grading of the
parity obstruction.
Usage: python goldbach_class_parity.py [D [NMAX]]      (D = 23 or 31; default NMAX = 4e6; ~3 min)
Output: results/goldbach_class_parity_{D}.txt
"""
import os
import sys
import time
from math import sqrt, log
import numpy as np
from scipy.fft import rfft, irfft, next_fast_len
from scipy.integrate import quad

from linnik_class import r_form, sieve, FORMS

C2 = 0.6601618158468696          # twin prime constant prod_{l>2} (1 - 1/(l-1)^2)


def kron_D(n):
    """(-D/n) for prime D = 3 mod 4 and n > 0: the Legendre symbol (n/D)"""
    n = int(n) % D
    if n == 0:
        return 0
    return 1 if pow(n, (D - 1) // 2, D) == 1 else -1


def singular_series_2(N, spf):
    """SS_2(N) = 2 C_2 prod_{l | N, l > 2} (l-1)/(l-2) for even N (ordered Goldbach pairs p + q = N)"""
    val = 2 * C2
    m = N
    while m > 1:
        l = int(spf[m])
        if l > 2:
            val *= (l - 1) / (l - 2)
        while m % l == 0:
            m //= l
    return val


def li2(N):
    """int_2^{N-2} dt / (log t log(N - t))"""
    f = lambda t: 1.0 / (log(t) * log(N - t))
    v1, _ = quad(f, 2, N / 2, limit=200)
    return 2 * v1

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
D = int(sys.argv[1]) if len(sys.argv) > 1 else 23
NMAX = int(float(sys.argv[2])) if len(sys.argv) > 2 else 4 * 10 ** 6
RANGES = [(10 ** 5, 2 * 10 ** 5), (3 * 10 ** 5, 6 * 10 ** 5), (10 ** 6, 2 * 10 ** 6), (2 * 10 ** 6, 4 * 10 ** 6)]
BETAS = (3, 4, 6)


def spf_omega(X):
    """smallest prime factor (spf[1] = X + 1, so that 1 is coprime to everything) and Omega(m), m <= X."""
    spf = np.zeros(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if spf[p] == 0:
            sl = spf[p * p::p]
            sl[sl == 0] = p
    idx = np.arange(X + 1, dtype=np.int64)
    spf[spf == 0] = idx[spf == 0]
    spf[0] = 0
    spf[1] = X + 1
    om = np.zeros(X + 1, dtype=np.int64)
    primes = np.nonzero(sieve(X))[0]
    for p in primes:
        p = int(p)
        pk = p
        while pk <= X:
            om[pk::pk] += 1
            pk *= p
    return spf, om


def main():
    out = []

    def say(s=""):
        out.append(s)
        print(s, flush=True)

    t0 = time.time()
    (a0, b0, c0), (a1, b1, c1) = FORMS[D]
    r0 = r_form(a0, b0, c0, NMAX)
    r1 = r_form(a1, b1, c1, NMAX)
    r0[0] = r1[0] = 0.0
    cc = (r0 - r1) / 2
    isp = sieve(NMAX).astype(float)
    pi_N = np.cumsum(isp)                      # pi(N) ; sum_{p<N} uses pi_N[N-1]
    spf, om = spf_omega(NMAX)
    say(f"D = {D}, h = 3, forms {FORMS[D]}, tables and spf/Omega to {NMAX} [{time.time() - t0:.0f}s]")

    n = next_fast_len(2 * NMAX + 1)
    ISP = rfft(isp, n)

    def conv(w):
        return np.rint(irfft(ISP * rfft(w, n), n)[:NMAX + 1])

    def stats(u, lo, hi, label):
        S0 = conv(r0 * u)
        S1 = conv(r1 * u)
        Dg = conv(cc * cc * u)
        Ns = np.arange(lo, hi)
        s0, s1, dg = S0[Ns], S1[Ns], Dg[Ns]
        sc = (s0 - s1) / 2
        tot = np.mean(s0 + 2 * s1)          # ratio of means: the mean of the ratios is dominated by the N with few pairs
        sh0 = np.mean(s0) / tot
        sh1 = np.mean(s1) / tot
        ms = np.mean(sc ** 2)
        bias = np.mean(sc) / sqrt(ms) if ms > 0 else 0.0
        pin = pi_N[Ns - 1]
        ratio = ms / np.mean(dg) if np.mean(dg) > 0 else float('nan')
        say(f"  {label:<28s} mean S_0/2 = {np.mean(s0) / 2:9.1f}  S_1/2 = {np.mean(s1) / 2:9.1f}   shares C_0, C_1: {sh0:.4f} {sh1:.4f}"
            f"   bias {bias:+.3f}   Diag/pi(N) = {np.mean(dg / pin):.4f}   E S_c^2 / E Diag = {ratio:.3f}")
        return ratio

    say("\nColumns: mean numbers of ideals of norm N - p in C_0 and in C_1 (r_C/2); class shares E S_{C_0}/E(S_{C_0}+2S_{C_1}),"
        " E S_{C_1}/E(...) (1/3 each if equidistributed);\n  bias = mean S_c / rms S_c; Diag^u / pi(N) (= kappa^u);"
        " variance ratio E (S_c^u)^2 / E Diag^u (0.42-0.43 for u = 1, Conjecture 3.2: 1/2).")
    ones = np.ones(NMAX + 1)
    ones[0] = 0.0
    strata = {k: (om == k).astype(float) for k in (1, 2, 3)}
    table = []
    for lo, hi in RANGES:
        if lo >= NMAX:
            continue
        hi = min(hi, NMAX + 1)
        say(f"\nN in [{lo}, {hi}):")
        row = {'lo': lo}
        row['all'] = stats(ones, lo, hi, "u = 1")
        for k in (1, 2, 3):
            row[f'om{k}'] = stats(strata[k], lo, hi, f"Omega(N-p) = {k}")
        for beta in BETAS:
            z = lo ** (1.0 / beta)
            u = (spf >= z).astype(float)
            u[0] = 0.0
            row[f'b{beta}'] = stats(u, lo, hi, f"spf >= z = {z:.1f} (beta = {beta})")
            if beta == 3:
                for k in (1, 2):
                    row[f'b3om{k}'] = stats(u * strata[k], lo, hi, f"  ... and Omega = {k}")
        table.append(row)

    say("\nSummary: E S_c^2 / E Diag by range and weight")
    hdr = f"{'N from':>8s} {'u=1':>6s} {'Om=1':>6s} {'Om=2':>6s} {'Om=3':>6s} {'b=6':>6s} {'b=4':>6s} {'b=3':>6s} {'b3,Om1':>7s} {'b3,Om2':>7s}"
    say(hdr)
    for row in table:
        say(f"{row['lo']:8d} {row['all']:6.3f} {row['om1']:6.3f} {row['om2']:6.3f} {row['om3']:6.3f} {row['b6']:6.3f} {row['b4']:6.3f}"
            f" {row['b3']:6.3f} {row['b3om1']:7.3f} {row['b3om2']:7.3f}")

    # Goldbach pairs with q = Q_C(x, y): G_C(N) = sum_{p+q=N} r_C(q) against (2/h)(1 - chi(N)/(D-2)) SS_2(N) li_2(N)
    say("\nGoldbach stratum (even N): G_C(N) = sum_{p+q=N} r_C(q) against the main term (2/3)(1 - chi(N)/(D-2)) SS_2(N) li_2(N),"
        "\n  G(N) = #{p+q=N} against SS_2(N) li_2(N); 4000 random even N per range; means (sd of the mean in brackets).")
    G0 = conv(r0 * strata[1])
    G1 = conv(r1 * strata[1])
    Gall = conv(isp)
    rng = np.random.default_rng(7)
    for lo, hi in RANGES:
        if lo >= NMAX:
            continue
        hi = min(hi, NMAX + 1)
        Ns = rng.integers(lo // 2, hi // 2, size=4000) * 2
        rows = {1: [], -1: [], 0: []}
        base = []
        for N in Ns:
            N = int(N)
            ss = singular_series_2(N, spf) * li2(N)
            chiN = kron_D(N)
            fac = (2 / 3) * (1 - chiN / (D - 2)) * ss
            base.append(Gall[N] / ss)
            rows[chiN].append((G0[N] / fac, G1[N] / fac, (G0[N] + 2 * G1[N]) / 4 / ss))
        base = np.array(base)
        say(f"  N in [{lo}, {hi}): G/(SS_2 li_2) = {base.mean():.4f} [{base.std() / sqrt(len(base)):.4f}]")
        R = {}
        for t in (1, -1, 0):
            a = np.array(rows[t])
            if len(a) == 0:
                continue
            R[t] = a[:, 2].mean()
            say(f"     chi(N) = {t:+d} (n = {len(a)}): G_C0/main = {a[:, 0].mean():.4f} [{a[:, 0].std() / sqrt(len(a)):.4f}],"
                f"  G_C1/main = {a[:, 1].mean():.4f} [{a[:, 1].std() / sqrt(len(a)):.4f}],  split share (G_C0+2G_C1)/(4G) vs (1-chi/(D-2))/2:"
                f" {a[:, 2].mean() / base.mean():.4f} vs {(1 - t / (D - 2)) / 2:.4f}")
        if 1 in R and -1 in R:
            say(f"     ratio of the split-pair counts chi(N)=+1 : chi(N)=-1 = {R[1] / R[-1]:.4f};  (D-3)/(D-1) = {(D - 3) / (D - 1):.4f}"
                f" [local factor 1 - chi(N)/(D-2)],  (D-2)/D = {(D - 2) / D:.4f} [1 - chi(N)/(D-1)]")

    # cross-class coincidence C(m) = r_0(m) r_1(m) on the Chen set (beta = 3), by Omega: C = 0 on primes, 8 on the semiprimes
    # with two non-principal split factors (density 4/9), 0 on the other semiprimes. 'HBT' weights w_0 = r_0 - C/2,
    # w_1 = r_1 - C/4 keep the primes and kill those semiprimes; purity = share of the weight on Omega = 1.
    say("\nCross-class coincidence C(m) = r_0(m) r_1(m) on the Chen set spf(N-p) >= N^{1/3}, means over N of sum_p ... by Omega(N-p):"
        "\n  predicted: C = 0 on Omega = 1; on Omega = 2, C-sum / r_0-sum = 4/3 and (C/4)-sum / r_1-sum = 1/3;"
        " w_0 = r_0 - C/2, w_1 = r_1 - C/4.")
    for lo, hi in RANGES:
        if lo >= NMAX:
            continue
        hi = min(hi, NMAX + 1)
        z = lo ** (1.0 / 3)
        u = (spf >= z).astype(float)
        u[0] = 0.0
        Ns = np.arange(lo, hi)
        R0, R1, CC = {}, {}, {}
        for k in (1, 2, 3):
            uk = u * strata[k]
            R0[k] = conv(r0 * uk)[Ns].mean()
            R1[k] = conv(r1 * uk)[Ns].mean()
            CC[k] = conv(r0 * r1 * uk)[Ns].mean()
        W0 = {k: R0[k] - CC[k] / 2 for k in (1, 2, 3)}
        W1 = {k: R1[k] - CC[k] / 4 for k in (1, 2, 3)}
        say(f"  N in [{lo}, {hi}): z = {z:.1f}")
        say(f"     Omega = 1, 2, 3:  r_0-sums {R0[1]:9.1f} {R0[2]:9.1f} {R0[3]:8.1f};  r_1-sums {R1[1]:9.1f} {R1[2]:9.1f} {R1[3]:8.1f};"
            f"  C-sums {CC[1]:6.1f} {CC[2]:9.1f} {CC[3]:8.1f}")
        say(f"     Omega = 2: C/r_0 = {CC[2] / R0[2]:.4f} (4/3 = 1.3333),  (C/4)/r_1 = {CC[2] / 4 / R1[2]:.4f} (1/3 = 0.3333)")
        say(f"     HBT weights w_0 by Omega: {W0[1]:9.1f} {W0[2]:9.1f} {W0[3]:8.1f};  w_1: {W1[1]:9.1f} {W1[2]:9.1f} {W1[3]:8.1f}")
        pur0 = R0[1] / (R0[1] + R0[2] + R0[3])
        pur0h = W0[1] / (W0[1] + W0[2] + W0[3])
        pur1 = R1[1] / (R1[1] + R1[2] + R1[3])
        pur1h = W1[1] / (W1[1] + W1[2] + W1[3])
        say(f"     prime purity (share of Omega = 1): C_0 {pur0:.3f} -> {pur0h:.3f} with w_0;  C_1 {pur1:.3f} -> {pur1h:.3f} with w_1")

    # the distribution of c(m) on the sifted strata (the 'parity grading')
    lo = RANGES[-1][0]
    z = lo ** (1.0 / 3)
    say(f"\nDistribution of c(m) on m < {NMAX} with spf(m) >= {z:.1f}, by Omega(m)  (c = (r_0 - r_1)/2):")
    for k in (1, 2, 3):
        sel = (spf >= z) & (om == k) & (np.arange(NMAX + 1) >= 2)
        vals = cc[sel]
        if len(vals) == 0:
            say(f"  Omega = {k}: none below {NMAX}")
            continue
        uniq, cnt = np.unique(vals, return_counts=True)
        dist = ", ".join(f"c={int(v):+d}: {c / len(vals):.4f}" for v, c in zip(uniq, cnt) if c / len(vals) >= 0.002)
        say(f"  Omega = {k}: n = {len(vals)},  mean c = {vals.mean():+.4f},  mean c^2 = {np.mean(vals ** 2):.4f},  {dist}")
    say("  (Omega = 1: c = +2 on the principal split primes, -1 on the non-principal ones, 0 on the inert ones, densities"
        " 1/6, 1/3, 1/2;\n   Omega = 2, both primes split: c = c(q_1) c(q_2) = +4, -2, +1 with densities 1/9, 4/9, 4/9 of"
        " the split pairs; the mean is 0 on every stratum.)")
    say(f"\n[{time.time() - t0:.0f}s]")
    os.makedirs(RES, exist_ok=True)
    path = os.path.join(RES, f'goldbach_class_parity_{D}.txt')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write("\n".join(out) + "\n")
    print("written", path)


if __name__ == "__main__":
    main()
