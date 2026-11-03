r"""
Part V, first stage: the Linnik-type sum S_D(N) = sum_{p <= N} r_D(N - p) for the nine discriminants of class number
one, against the main term (2 pi / sqrt D) S_D(N) pi(N) with the singular series

    S_D(N) = prod_{l | N} (1 - chi_{-D}(l)/l) * prod_{l not | N} (1 + chi_{-D}(l)/(l (l - 1)))

(manuscript, Proposition 1.2: w_D sum_{(d,N)=1} chi(d)/phi(d) = w_D L(1,chi) S_D(N), and w_D L(1,chi) = 2 pi/sqrt D).
Also counted: the primes p < N with N - p a norm (r_D(N - p) > 0), and those with N - p a norm with at most two
prime factors (what a 'prime + P_2-norm' theorem would bound from below); and, for D = 4, the proportion of p with
(N - p, P(z)) = 1 and with no inert prime factor below y, against the dimension-1 and dimension-1/2 products
(manuscript, Section 3).  Exact integer counts (r_D(n) = w_D sum_{d | n} chi_{-D}(d) for class number one);
floating point only in the Euler products.   Usage: python linnik_heegner.py [NMAX]   (default 200000; ~2 min)
"""
import math
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
os.makedirs(RES, exist_ok=True)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
DISCS = [3, 4, 7, 8, 11, 19, 43, 67, 163]
W = {3: 6, 4: 4}


def kronecker(a, n):
    """Kronecker symbol (a/n), n >= 1."""
    if n == 1:
        return 1
    res = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            res = -res
    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def chi_table(D, M):
    mod = D if D % 4 == 3 else 4 * D
    base = np.array([kronecker(-D, n) if n else 0 for n in range(mod)], dtype=np.int64)
    return base[np.arange(M + 1) % mod]


def r_table(D, M):
    """r_D(n) = w_D sum_{d | n} chi_{-D}(d), 0 <= n <= M (r_D(0) set to 0)."""
    chi = chi_table(D, M)
    r = np.zeros(M + 1, dtype=np.int64)
    for d in range(1, M + 1):
        c = int(chi[d])
        if c:
            r[d::d] += c
    r[0] = 0
    return W.get(D, 2) * r


def primes_upto(M):
    s = np.ones(M + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(M ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def bigomega(M, P):
    """Omega(n), 0 <= n <= M"""
    Om = np.zeros(M + 1, dtype=np.int64)
    for p in P:
        q = int(p)
        while q <= M:
            Om[q::q] += 1
            q *= int(p)
    return Om


def conductor(D):
    return D if D % 4 == 3 else 4 * D


def ramified_factor(D, N):
    """1 + c_D(N), c_D(N) = phi(f)^{-1} sum_{a mod f, (a,f)=1} chi_{-D}(N - a), f the conductor: the second half of the
    hyperbola carries chi(N - p), whose mean over the primes is c_D(N).  For odd D this is 1 - chi_{-D}(N)/(D - 1);
    for D = 4, 8 it is 1."""
    f = conductor(D)
    units = [a for a in range(1, f) if math.gcd(a, f) == 1]
    return 1 + sum(kronecker(-D, (N - a) % f) if (N - a) % f else 0 for a in units) / len(units)


def euler_part(D, N, P):
    """prod_{l | N} (1 - chi(l)/l) prod_{l not | N} (1 + chi(l)/(l(l-1))) over the primes in P (tail O(1/max P))."""
    val = 1.0
    for l in P:
        l = int(l)
        c = kronecker(-D, l)
        if c == 0:
            continue
        val *= (1 - c / l) if N % l == 0 else (1 + c / (l * (l - 1)))
    return val


def singular_series(D, N, P):
    """S_D(N) = (1 + c_D(N)) * Euler part (manuscript, Proposition 1.2)."""
    return ramified_factor(D, N) * euler_part(D, N, P)


def factor_small(m, P):
    out = {}
    for q in P:
        q = int(q)
        if q * q > m:
            break
        while m % q == 0:
            out[q] = out.get(q, 0) + 1
            m //= q
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def mixed_sieve_table(D, N, P, lines):
    """The sieve of Section 3: among the primes p < N, those with N - p free of inert primes below y (and, in the second
    block, also free of all primes below z), against the products prod (1 - 1/(q-1)); and the structure of the survivors:
    norms, norm x one inert prime, norm x two inert primes (inert primes to an odd power)."""
    pp = P[P < N]
    m = N - pp
    small = [int(q) for q in P if q < 10 ** 5]
    inert = [q for q in small if kronecker(-D, q) == -1]
    lines.append(f"   D = {D}, N = {N}, pi(N) = {len(pp)}: survivors of 'no inert prime < y' (and 'no prime < z'), with their structure")
    lines.append(f"   {'z':>6} {'y':>6} {'surv':>6} {'surv/pi':>8} {'product':>8} {'ratio':>6} {'norm':>6} {'q0*norm':>8} {'q0q1*norm':>9} {'norm,Om<=2':>10} {'norm,Om<=3':>10}")
    for z, y in ((1, 10), (1, 100), (1, 447), (1, 1000), (1, 10000), (10, 447), (30, 447), (100, 447), (100, 1000), (300, 1000)):
        mask = np.ones(len(m), dtype=bool)
        prod = 1.0
        for q in small:
            if q >= max(z, y):
                break
            if q < z or (q < y and kronecker(-D, q) == -1):
                mask &= (m % q != 0)
                if N % q:
                    prod *= 1 - 1 / (q - 1)
        surv = m[mask]
        cls = {0: 0, 1: 0, 2: 0}
        om2 = om3 = 0
        for mm in surv:
            f = factor_small(int(mm), P)
            ninert = sum(1 for q, e in f.items() if e % 2 == 1 and kronecker(-D, q) == -1)
            cls[min(ninert, 2)] += 1
            if ninert == 0:
                tot = sum(f.values())
                om2 += tot <= 2
                om3 += tot <= 3
        n = len(pp)
        ratio = len(surv) / n / prod if prod > 0 else float('nan')       # prod = 0: 2 inert and N odd (parity obstruction)
        lines.append(f"   {z:>6} {y:>6} {len(surv):>6} {len(surv) / n:>8.4f} {prod:>8.4f} {ratio:>6.3f} {cls[0]:>6} {cls[1]:>8} {cls[2]:>9} {om2:>10} {om3:>10}")


def main():
    P = primes_upto(NMAX)
    PE = primes_upto(2 * 10 ** 6)
    Om = bigomega(NMAX, P)
    Ns = [n for n in sorted({10 ** 4, 10 ** 4 + 1, 2 * 10 ** 4, 3 ** 9, 5 * 10 ** 4, 2 ** 16, 10 ** 5, 10 ** 5 + 3, 2 * 3 * 5 * 7 * 11 * 13 * 2,
                             150000, 2 * 10 ** 5, 2 * 10 ** 5 + 1}) if 1000 < n <= NMAX]
    lines = [f"S_D(N) = sum_(p<N) r_D(N-p) against the main term (2 pi/sqrt D) S_D(N) pi(N), N <= {NMAX}; exact counts.",
             "S_D(N) = (1 + c_D(N)) prod_(l|N)(1 - chi(l)/l) prod_(l not|N)(1 + chi(l)/(l(l-1))), c_D(N) the mean of chi_{-D}(N - a)",
             "over the units a modulo the conductor (= -chi_{-D}(N)/(D-1) for odd D, 0 for D = 4, 8).",
             "Columns: N, pi(N), S_D(N), 1 + c_D(N), singular series, main term, ratio S_D/main, #{p: N-p norm}, that over pi(N),",
             "#{p: N-p norm with Omega(N-p) <= 2}, that over pi(N).", ""]
    last = {}
    for D in DISCS:
        r = r_table(D, NMAX)
        lines.append(f"==== D = {D}, w_D = {W.get(D, 2)}, 2 pi/sqrt D = {2 * math.pi / math.sqrt(D):.6f} ====")
        lines.append(f"{'N':>8} {'pi(N)':>7} {'S_D(N)':>10} {'1+c_D':>6} {'sing.ser':>9} {'main':>10} {'ratio':>7} {'#norm':>6} {'/pi':>7} {'#norm,Om<=2':>11} {'/pi':>7}")
        for N in Ns:
            pp = P[P < N]
            m = N - pp
            rv = r[m]
            S = int(rv.sum())
            rf = ramified_factor(D, N)
            ss = singular_series(D, N, PE)
            main_t = 2 * math.pi / math.sqrt(D) * ss * len(pp)
            isnorm = rv > 0
            n1 = int(isnorm.sum())
            n2 = int((isnorm & (Om[m] <= 2)).sum())
            last[D] = S / main_t
            lines.append(f"{N:>8} {len(pp):>7} {S:>10} {rf:>6.3f} {ss:>9.5f} {main_t:>10.1f} {S / main_t:>7.4f} {n1:>6} {n1 / len(pp):>7.4f} {n2:>11} {n2 / len(pp):>7.4f}")
        lines.append("")
    lines.append("ratio at the largest N: " + ", ".join(f"D={D}: {last[D]:.4f}" for D in DISCS))
    lines.append("")
    # sieve densities for D = 4 at N = NMAX: dimension 1 (all primes < z) and 1/2 (inert primes < y)
    D, N = 4, NMAX
    pp = P[P < N]
    m = N - pp
    lines.append(f"Sieve densities of N - p, D = 4, N = {N}: (N-p, P(z)) = 1 against prod_(q<z, q not | N)(1 - 1/(q-1)) [dimension 1],")
    lines.append("and (N-p, product of inert primes q < y) = 1 against prod_(q<y inert, q not | N)(1 - 1/(q-1)) [dimension 1/2].")
    for z in (10, 30, 100, 300, 1000):
        small = [int(q) for q in P[P < z]]
        mask = np.ones(len(m), dtype=bool)
        prod = 1.0
        for q in small:
            mask &= (m % q != 0)
            if N % q:
                prod *= 1 - 1 / (q - 1)
        lines.append(f"   all primes < {z:>4}: observed {mask.mean():.4f}, product {prod:.4f}, ratio {mask.mean() / prod:.4f}")
    for y in (100, 1000, 10000, 100000):
        inert = [int(q) for q in P[P < y] if kronecker(-D, int(q)) == -1]
        mask = np.ones(len(m), dtype=bool)
        prod = 1.0
        for q in inert:
            mask &= (m % q != 0)
            if N % q:
                prod *= 1 - 1 / (q - 1)
        lines.append(f"   inert primes < {y:>6}: observed {mask.mean():.4f}, product {prod:.4f}, ratio {mask.mean() / prod:.4f}")
    lines.append("")
    lines.append("Section 3, the sieve on the unweighted sequence {N - p}: survivors of the inert-prime sieve and their structure.")
    lines.append("Columns: z (all primes < z sifted), y (inert primes < y sifted), survivors, their proportion of pi(N), the product")
    lines.append("prod_(q sifted, q not | N)(1 - 1/(q-1)), ratio, and among the survivors: norms, norm x one inert prime (odd power),")
    lines.append("norm x two or more inert primes, and norms with Omega <= 2, <= 3.  sqrt(N) = 447 for N = 2*10^5.")
    for D in (3, 4, 7, 8, 11, 163):
        for N in (NMAX, NMAX // 2 + 3 if (NMAX // 2 + 3) % 2 else NMAX // 2 + 4):
            mixed_sieve_table(D, N, P, lines)
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RES, 'linnik_heegner.txt'), 'w', encoding='utf-8') as fh:
        fh.write(text + "\n")


if __name__ == '__main__':
    main()
