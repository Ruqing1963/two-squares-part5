# Sums of two squares and Heegner norms with few prime factors, V: a prime plus a Heegner norm

First-stage draft (not archived). Part V of the series after Parts I–IV
([two-squares-theta-level](https://github.com/Ruqing1963/two-squares-theta-level),
[two-squares-mean-square](https://github.com/Ruqing1963/two-squares-mean-square),
[two-squares-short-intervals](https://github.com/Ruqing1963/two-squares-short-intervals),
[two-squares-part4](https://github.com/Ruqing1963/two-squares-part4)). The repository holds

* `paper/linnik_heegner.tex` (+ PDF): the first-stage manuscript;
* `code/linnik_heegner.py`, `results/linnik_heegner.txt`: the counting script and its output.

## The problem and what the draft contains

For the nine discriminants D ∈ {3, 4, 7, 8, 11, 19, 43, 67, 163} of class number one, write N = p + N_D(x, y) with p a
rational prime and N_D the principal form, through S_D(N) = Σ_{p<N} r_D(N − p), r_D(n) = w_D Σ_{d|n} χ_{−D}(d).
For D = 4 this is the Hardy–Littlewood problem N = p + x² + y², solved by Linnik with the dispersion method.

1. **Theorem 1.2.** S_D(N) = (2π/√D) 𝔖_D(N) π(N) + O(N (log log N)²/(log N)²) with the singular series in closed form,
   𝔖_D(N) = (1 + c_D(N)) ∏_{ℓ|N}(1 − χ(ℓ)/ℓ) ∏_{ℓ∤N}(1 + χ(ℓ)/(ℓ(ℓ−1))), c_D(N) = −χ_{−D}(N)/(D−1) for odd D and 0 for
   D = 4, 8 (the factor at the ramified prime, from the second half of the hyperbola). Proved for all nine D by the
   divisor-switching argument of Friedlander–Iwaniec (arXiv:1807.06648) **under Hypothesis H(N)**: the
   Bombieri–Friedlander–Iwaniec level √N (log N)^A for the single residue class N (proved in the literature for a fixed
   residue class, open as stated). Unconditional for D = 4 by Linnik. The data confirm the formula to 1% for all nine D
   at N ≤ 2·10⁵; without the factor 1 + c_D(N) the ratios for D = 3 are 0.5, 1, 1.5 according to N mod 3.
2. **Proposition 1.4 (the sieve).** {N − p} has dimension 1 and level 1/2 (Chen's setting), but the condition "N − p is
   a norm" is not a sieve condition. The norm-weighted sequence r_D(N − p) has dimension 1 on average; its level of
   distribution at N^θ is equivalent to the equidistribution of primes p ≡ N (mod qb) with q ≤ N^θ, b ≍ √N, in the single
   class N, which no known theorem gives for θ > 0. So no "N = p + P_k with N − p a norm" is proved, for any k, and the
   square-root barrier of Parts I–IV reappears in the residue class.
3. **The sieve on the unweighted sequence (Section 3).** Sifting {N − p} by the inert primes below N^{1/2−ε} with
   Iwaniec's half-dimensional sieve (sieving limit 1) and the Bombieri–Vinogradov level 1/2 gives, unconditionally,
   ≫ N/(log N)^{3/2} primes p with N − p a norm up to at most two inert primes above N^{1/2−ε} (Theorem 3.1; N even
   when 2 is inert). Removing those inert primes is a Chen-type comparison of constants: the exceptional set is of the
   same order (Proposition 3.2), the crude Selberg bound loses by an order of magnitude, and the data (Table 2: at
   N = 2·10⁵ the exceptional survivors are 21–51 % of all survivors) suggest the refined comparison may close; it is
   not done here, and no "N = p + P_k with N − p a norm" is proved. The reversed problem N = N_D(x,y) + P_3 is
   unconditional from the lattice-point level 1/2 and Richert's weighted sieve (Proposition 3.3), weaker than Linnik
   for D = 4.
4. **Numerics.** Table 1: S_D(N) against the main term for eleven N in [10⁴, 2·10⁵]; the proportion of primes p < N with
   N − p a norm (0.15–0.50 at N = 2·10⁵) and a norm with Ω ≤ 2 (0.11–0.34; zero at odd N for D ≡ 3 mod 8, where 2 is
   inert and an even norm is divisible by 4); sieve densities of {N − p} against the dimension-1 and dimension-1/2
   products. Table 2: the inert-prime sieve at N = 2·10⁵ (survivors against the product V, their structure n₁ / n₁q₀,
   and the P₂ norms left when the primes below 100 are sifted too).

## A single ideal class (`notes/linnik_class.tex`, working note, 13 pages)

For a prime D with h(−D) > 1 and a class C, S_C(N) = Σ_{p<N} r_C(N − p) splits as S_gen + (w/h) Σ_{ψ≠1} ψ̄(C) S_ψ with
S_ψ(N) = Σ_{p<N} c_ψ(N − p), c_ψ the coefficients of the weight-one dihedral cusp form of Part V.

* The genus part is the Linnik problem for (w/h)(1*χ): under H(N) its main term is (2π/√D) 𝔖_D(N) π(N) for every h
  (the class number cancels), the same for every class.
* **Theorem 1.3 (unconditional):** Σ_{N≤X} |S_ψ(N)|² ≪ X² log X (Wilton-type bound for weight one + Parseval). So all
  N ≤ X outside O(X^{1−2η} log X) have |S_C(N) − S_gen(N)| ≤ N^{1/2+η} for every class, and the classes are asymptotically
  equidistributed for all N ≤ X but O((log X)⁴).
* **Theorem 1.5:** under H(N) and the dual Titchmarsh bound S_ψ(N) = o(N/log N) (Assing–Blomer–Li, Theorem 1.3, if it
  extends to weight one with nebentypus), every class has the asymptotic for every large N.
* Numerics (D = 23, 31, N ≤ 4·10⁶): both classes match the genus formula to 10⁻⁴; S_c(N) ≈ 0.12 √N, unbiased.
* **The variance (Section 3, Conjecture 3.2):** the diagonal Σ_p c_ψ(N−p)² is κ_ψ π(N) (Rankin–Selberg residue of
  Part V), but the mean square of S_ψ(N) is 0.41–0.43 of it and rising (0.433 at N ≈ 1.2·10⁷). Lemma 3.1: by modularity
  the exponential sum of ϑ_ψ is exponentially small at a/q for q² ≪ X, and these major arcs carry half the L² mass
  of the prime exponential sum. Conjecture: mean square ~ ½ κ_ψ π(N).
* **Section 3.1 (what is proved):** with the smoothed V = Σ e^{−2N/X} S_ψ(N)² = ∫|P|²|F|² and the exact identity
  Diag = ∫|P|²·∫|F|², Theorem 3.4 proves that the major arcs of order X^{1/2}(log X)^{−B} (width (log X)^{(B−1)/2}/X)
  contribute O(X^{−A}) to V. The conjecture is reduced to (M) these arcs carry half the mass of |P|² (Hardy–Littlewood;
  needs a short-interval BDH estimate with weight r for moduli up to √X, open even under GRH, which reaches X^{1/3}) and
  (D) |P|² and |F|² decorrelate on the minor arcs (a pair-correlation statement, as in Goldston–Montgomery). Numerics
  (`code/linnik_class_split.py`, X = 5·10⁴, 2·10⁵, 6·10⁵): mass of the arcs within 0.006 of HL; (D)-ratio 0.71–0.88,
  increasing; V/Diag = 0.397, 0.408, 0.415, fitted slope 0.504, 0.505.
* **Section 3.2 (class number five):** two dihedral modes S_1, S_2 (D = 47, 79; `code/linnik_class_h5.py`, N ≤ 4·10⁶).
  Each shows the same deficit (mean square 0.42–0.45 of its diagonal), the cross diagonal Σ c_1 c_2(N − p) is below
  10⁻⁴ π(N) (no pole of the cross Rankin–Selberg series), corr(S_1, S_2) = 0.02–0.04 and slowly decreasing; the class
  deviations then have variances (2/25)(a_1²κ_1 + a_2²κ_2)π(N) with golden-ratio weights a_j (ratio 1.307 predicted,
  1.309 measured, for C₁^{±1} against C₁^{±2} at D = 47), up to the factor ½ of Conjecture 3.2.
* **Section 3.3 (sieve weights, towards Chen in a single class):** for sieve weights λ_d of level X^θ on N − p, the
  class deviation is a combination of dilates ϑ_ψ(sdz) (Hecke, Lemma 3.5), each with the Wilton bound d^{−1/2+ε}√Y log Y
  (Proposition 3.6), so Σ_{N≤X} |S^λ_ψ(N)|² ≪ X^{2+θ+ε} and the class deviation of every sieve functional of level θ < 1 is
  O(N^{1/2+θ/2+η}) for almost all N, unconditionally (Theorem 3.7): a sieve result for the genus holds for each class,
  for almost all N (the switching step of Chen's argument needs a Vinogradov-type bound for Σ Λ(m)c_ψ(m)e(mα), not
  covered). Major arcs of order X^{1/2−θ} keep the depletion (Proposition 3.13). Remark 3.12: there is no class-group
  grading of the parity obstruction (c(q₁q₂) = (−1)² only on integers with non-principal prime factors; the mean of c on
  every Ω-stratum is 0), only a Chebyshev bias of order (log N)^{−3/2} with alternating sign in Ω. Numerics
  (`code/goldbach_class_parity.py`, D = 23, 31, N ≤ 4·10⁶, Table 2): the variance ratio rises monotonically with the
  sifting, unsifted 0.41–0.43, sharp sieve by the primes below N^{1/6}, N^{1/4}, N^{1/3}: 0.41–0.45, 0.43–0.46,
  0.45–0.50, strata Ω = 1, 2, 3: 0.54–0.61, 0.62–0.68, 0.73–0.81; the classes stay equidistributed in every stratum.
* **Proposition 3.8 (Vinogradov–Siegel–Walfisz for ϑ_ψ):** Σ_{n≤X} Λ(n)c_ψ(n)e(nα) ≪ X(log X)^{−A} uniformly in α.
  Major arcs: the twists L(s, ϑ_ψ ⊗ ξ) are the Hecke L-functions L_K(s, ψ·(ξ∘N)), and ψ·(ξ∘N) is never real because h is
  odd (Ψ(𝔭)/Ψ(𝔭̄) = ψ(𝔭)² ≠ 1), so there is no Siegel zero and Siegel–Walfisz for K is effective. Minor arcs: Vaughan's
  identity, the Hecke inversion c(mn) = Σ_{g|(m,n)} μ(g)χ(g)c(m/g)c(n/g), Wilton for type I, the standard bilinear
  bound for type II. Corollary 3.9: the same for the strata Ω(m) = k and for Chen's switching set q₁q₂q₃ (peel off the
  largest prime; the square-divisible part is O(X^{1−1/k})).
* **Theorem 3.10:** Σ_{N≤X} |S^{(k)}_ψ(N)|² ≪ X³(log X)^{−A}, likewise for the switching sum with prime or sieve weights on
  N − m: the class deviation of every Ω-stratum and of Chen's switching step is ≤ N(log N)^{−A} for all N but
  O(X(log X)^{−A}). So, for almost all N, Chen's theorem in a single class is equivalent to Chen's theorem for the genus.
* **Theorem 3.11 (binary Goldbach in a single ideal class, almost all even N, unconditional):** G_C(N) = Σ_{p+q=N} r_C(q) =
  (2/h)(1 − χ(N)/(D−2)) 𝔖₂(N) li₂(N) + O(N(log N)^{−A}) for all even N ≤ X but O(X(log X)^{−A}), the same main term for
  every class (genus part: van der Corput–Chudakov–Estermann with the weight 1 + χ(q), a character mod D; the local
  factor at D counts the D − 2 admissible residue pairs, hence (D−2), not (D−1)). Numerics in
  `results/goldbach_class_parity_{23,31}.txt`: G/(𝔖₂ li₂) and G_C/main within 1–2 % of 1 (Chebyshev bias of the principal
  primes), and the ratio of the split-pair counts for χ(N) = +1 against −1 matches (D−3)/(D−1), not (D−2)/D.
* **Remark 3.14:** Assing–Blomer–Li's Theorem 1.3 uses the trivial nebentypus only through their Voronoi formula
  (Lemma 7.4); for ω = χ_{−D} the Kloosterman sums acquire ω(c) (cusp pair (∞, 0)) or become twisted (D | q), so their
  Theorem 2.4 must be re-proved with the Kuznetsov formula for that cusp pair with character. Not done; T_ψ stays a hypothesis.

## Reproducing

```
cd code
python linnik_heegner.py            # ~2 min; writes results/linnik_heegner.txt
python linnik_class.py 23 4e6       # single classes, D = 23 (or 31); ~2 min; results/linnik_class_23.txt
python linnik_class_variance.py 23 4e6    # mean square vs diagonal, Wilton supremum
python linnik_class_arcs.py 23 3e5        # depletion profile rho(q) and the major-arc model
python linnik_class_smooth.py 23 3e4 1e5  # smoothed variance vs the beta-integrated model (~30 min)
python linnik_class_split.py 23 2e5       # V split into major/minor arcs; env KS, CQS select K and Q/sqrt X
python linnik_class_h5.py 47 4e6          # class number five: two modes, ~1 min; results/linnik_class_h5_47.txt (or 79)
python goldbach_class_parity.py 23 4e6    # sieve weights, Omega-strata, Goldbach pairs in a class vs the main term, ~4 min; results/goldbach_class_parity_23.txt (or 31)
```

Python ≥ 3.10 with numpy. Exact integer counts; floating point only in the Euler products. The paper compiles with
`pdflatex` (twice) from `paper/`.
