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

## Reproducing

```
cd code
python linnik_heegner.py            # ~2 min; writes results/linnik_heegner.txt
```

Python ≥ 3.10 with numpy. Exact integer counts; floating point only in the Euler products. The paper compiles with
`pdflatex` (twice) from `paper/`.
