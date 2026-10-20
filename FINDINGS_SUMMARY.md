# Chowla's Cosine Problem: Findings from Two Days of Computation

*Noah Puckett, October 2026. Informal research summary.*

## 1. The Problem

Pick some positive whole numbers, say 1, 3 and 4. Each number gives a wave, and we add the waves up:

  f(x) = cos(x) + cos(3x) + cos(4x).

At x = 0 every wave is at its peak, so f(0) = 3, the number of numbers picked. Each wave also spends as much time below zero as above it, so f averages out to zero over a full period and has to dip below zero somewhere. For this example the lowest point is about −1.37.

Chowla's question is how shallow that dip can be made. The right yardstick turns out to be the square root of the number of numbers chosen, so for a set A of n numbers we use the **score**

  score(A) = −(lowest value of f) / √n.

Lower is better. The example above scores 1.37 / √3 ≈ 0.79.

**Chowla's conjecture** says the score cannot be pushed to zero: there is a fixed constant c > 0 such that every set of n positive integers dips to −c√n or below. It has been open for decades. Epoch AI lists it among its FrontierMath Open Problems, where a construction pushing the score toward zero would count as a disproof.

What was known before this work:

| Result | Value | Reference |
|---|---|---|
| Classical construction (differences of a Sidon set) | score → 1/√2 ≈ 0.707 | folklore |
| Two-sided Sidon construction | score → 1/2 | Bu et al., arXiv:2610.17160 (Oct 2026) |
| Every set dips at least this deep | order n^{1/5} | Bedert (2025); Jin, Milojević, Tomon, Zhang (2025); Shankar (2026) |

So the true answer is somewhere between n^{1/5} and (1/2)√n. Chowla's conjecture says the √n end is correct.

## 2. What I Did

Bu et al. start from a **Sidon set**, a set in which every sum of two elements is different. They take a random subset of one and build the final set from it with a short formula. A random subset is enough for their limit theorem, but at any fixed size some subsets do better than others.

I replaced the random choice with a search. The program swaps one element of the subset for another and keeps the swap when the dip gets shallower, sometimes accepting a worse swap so it does not get stuck. The searches ran for about two days on a 16-core PC and a MacBook. Every score below was computed on a grid at least 64 times finer than the largest number in the set, refined with Newton's method at the lowest points, and recomputed on the second machine.

## 3. Results

**Table 1.** Scores at the sizes reported by Bu et al. Lower is better.

| n | Classical | Bu et al. (random subset) | This work (optimized subset) |
|---|---|---|---|
| 1,560 | 0.714 | 0.725 | **0.643** |
| 3,540 | 0.711 | 0.708 | **0.624** |
| 6,320 | 0.710 | 0.674 | **0.613** |
| 14,280 | 0.709 | 0.655 | **0.598** |
| 25,440 | 0.709 | 0.635 | **0.590** |

The optimized subsets improve on every published value, by 0.045 to 0.08. This is what one would expect from optimizing a step that was random. The value is in having concrete, checkable sets with better constants. Results at intermediate sizes (2,450 to 11,990) all fall between 0.59 and 0.64.

**A separate result at n = 1,560: 0.594.** A different method treats every candidate number as a dial between 0 and 1, moves the dials to reduce the dip, and slowly forces each to 0 or 1. Starting from the 0.643 set, two rounds of this reached 0.594. It did not help at larger sizes.

**Table 2.** Best scores found for small n. These come from a search, so they are not proven optimal.

| n | 6 | 8 | 10 | 12 | 14 | 16 | 18 |
|---|---|---|---|---|---|---|---|
| best found | 0.650 | 0.692 | 0.651 | 0.639 | 0.620 | 0.648 | 0.646 |

At n = 6 the search found {2, 4, 8, 12, 14, 16}, which is Mercer's optimal set {1, 2, 4, 6, 7, 8} doubled. Doubling every element does not change the score. No result in this project went below 0.5.

## 4. The Thesis

Optimizing the subset in the two-sided Sidon construction gives better constants than any published value at every size tested. The construction itself, however, appears to be limited to 1/2. Section 5 proves a version of this: for starting sets with no additive coincidences up to a given order, the best possible score approaches 1/2 as that order grows. The case of ordinary Sidon sets, which is what everyone actually uses, is still open. If that case goes the same way, beating 1/2 will need a construction of a different shape.

## 5. The Theorem

**The idea.** In the two-sided construction, with k starting numbers and about k² final numbers, the sum of waves can be written exactly as

  f = ρ² − Q,

where ρ is a sum of k cosines and Q is the sum of their squares. Since ρ² is never negative, f can only be very negative where ρ is close to zero, and there f ≈ −Q. Each cos² averages 1/2, so Q is usually about k/2. Dividing by √(k²) = k gives a score near 1/2. The theorem turns this into a proof for a large class of starting sets.

**Definitions.** A set S is a **Bₕ set** if all sums of h elements of S are different, apart from reordering the same elements. Sidon sets are B₂ sets.

**Theorem.** Let m ≥ 1 and let S be a B_{2m+1} set of k integers, all multiples of some integer d > (2m+1)/3. Build the set A from S as in Bu et al. Then

  score(A) ≥ (1 − λ_m)/2 − O(1/k),

where the constants λ_m are explicit and tend to 0 as m grows:

| Condition on S | Lower bound on score(A) |
|---|---|
| B₃, even integers | (√6 − 2)/2 ≈ 0.2247 |
| B₅, multiples of 2 | 0.3098 |
| B₇, multiples of 3 | 0.3547 |
| B₁₁, multiples of 4 | 0.4013 |
| B₁₅, multiples of 6 | 0.4252 |
| B_h with h → ∞ | → 1/2 |

**Proof sketch.**

1. *Weights.* The lowest value of f is at most its average against any weight w ≥ 0. We take w = q(ρ²)², where q is a polynomial of degree m chosen so that w is concentrated where ρ is near 0. For m = 1 this is w = (ρ² − b)².
2. *Exact averages.* The waves have frequencies 3s + 1/2. A product of N of them has a nonzero average only if Σ±(3s + 1/2) = 0, which means 3·Σ±s = −(P − M)/2, where P and M are the numbers of plus and minus signs. Since every s is a multiple of d and N/6 < d, this forces P = M. The Bₕ property then forces the plus terms to match the minus terms. So all averages of products of up to 2h of the cosines are exactly the same as for independent random waves, and they can be computed in closed form.
3. *Optimizing.* For m = 1 the averages give E[f·w] / E[w] = (k/2)(12 − 4β)/(3 − 2β + β²), up to lower-order terms, where b = βk/2. The minimum is at β = 3 + √6, where the ratio equals 2 − √6 ≈ −0.4495. Dividing by √(k(k−1)) gives the bound. For larger m the best q comes from a small eigenvalue problem, which produces the rest of the table.

**What this does not cover.** The Erdős–Turán and Singer sets used by Bu et al. and in Table 1 are Sidon sets, not B₃ sets, so the theorem does not apply to them directly. Extending it to Sidon sets is the main open question. The theorem also gives bounds below 1/2 for every fixed h. The value 1/2 appears only in the limit.

**Sanity check.** Random B₃ sets made of even numbers, with k from 10 to 26, score between 0.81 and 0.99. This only confirms that nothing contradicts the bound; it is not evidence that the bound is sharp.

## 6. What Did Not Work

I also tried five other approaches at n ≈ 1,560, where the bar to beat was 0.643.

| Approach | Best score | Why it failed |
|---|---|---|
| Rounding a smooth ramp of weights into a set | 0.93 | rounding errors act like noise and create deep dips |
| Sidon sets with a few deliberate coincidences | 0.755 | merging duplicate numbers costs more than it gains |
| Pairs of neighbouring numbers (a and a + 1) | 2.3 | thousands of duplicates |
| Two constructions at very different scales | about 0.91 | the two sets of dips line up |
| Evenly spaced sets | 8.6 | regular spacing produces one very deep dip |

A short calculation suggests why. For any construction of the form f = Ψ² − (average of Ψ²), the dip equals that average, and comparing with the variance of f gives score = 1/√(2(κ − 1)), where κ is the kurtosis of Ψ, a measure of how spiky it is. A score below 1/2 needs κ > 3. However, κ > 3 happens exactly when the frequencies satisfy extra additive relations, and those relations are what make numbers in the final set coincide. This is a heuristic explanation, not a proof.

## 7. How Reliable Is This?

- **Tables 1 and 2.** Computed, refined and recomputed on two machines. The actual sets are included with this summary, so anyone can check them directly.
- **Bu et al.'s values.** Taken from the table in arXiv:2610.17160 and recomputed from their formulas.
- **The theorem.** The B₃ case has a complete short proof, including the correction in step 2 that requires the elements of S to be multiples of d. It should still be checked by a professional mathematician. The rows for larger h use the same argument together with a numerical eigenvalue computation.
- **The kurtosis argument in Section 6** is a heuristic.
- **Tools.** All computation ran on my own machines. I used an AI assistant (Claude) to help write the code and work through the mathematics, and I checked the results by computation and by hand.

## 8. Open Questions

1. Does the theorem extend from Bₕ sets to ordinary Sidon sets?
2. Is there a construction that is not of the form "square minus its average" and scores below 1/2?
3. For small n, do the best scores stay near 0.62 to 0.65, or do they eventually drop?

## 9. Code and Data

Included alongside this summary: the optimized sets for every n in Table 1, the 0.594 set, the small-n sets, all search and verification scripts, and the result logs from both machines.
