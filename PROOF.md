# Barrier theorem for square-type (two-sided Sidon) constructions

## Setting
Let S = {s_1, ..., s_k} be a set of integers, theta_j(x) = (3 s_j + 1/2) x, c_j = cos theta_j,
rho = sum_j c_j, Q = sum_j c_j^2, and

  A(S) = { 3|s_i - s_j| : i < j }  U  { |3(s_i + s_j) + 1| : i < j }.

If S is Sidon (B_2), then |A| = k(k-1) and, for every real x,

  f_A(x) = rho(x)^2 - Q(x)     (Bu et al., arXiv 2610.17160, Prop. 2.1).

S is a **B_h set** if every equation a_1 + ... + a_h = b_1 + ... + b_h with a's, b's in S has only the trivial
solutions (the b's are a permutation of the a's).

## Lemma (moments are exactly Gaussian-model)
Average over x in [0, 4π). For a product of N <= 2h factors cos(theta_{j_1} x) ... cos(theta_{j_N} x), the average
equals 2^{-N} times the number of sign patterns with sum_t eps_t (3 s_{j_t} + 1/2) = 0. The 1/2-terms force equally many
+ and - signs; the B_h property then forces the + multiset to equal the - multiset. So the average equals the value for
**independent** uniform phases. Every moment of (rho, Q) of total degree <= 2h therefore equals its value in the
independent model, where E c^2 = 1/2, E c^4 = 3/8 and E c^6 = 5/16.

## Theorem 1 (B_3)
If S is a B_3 set of size k, then min_x f_A(x) <= -(0.4495 - O(1/k)) * k/2, so

  -min f_A / sqrt|A|  >=  0.2247 - O(1/k).

**Proof.** For any weight w >= 0 with E w > 0 we have min f <= E[f w] / E[w]. Take w = (rho^2 - b)^2, so that f w has
degree 6 in the c's and the Lemma applies with h = 3. Exactly, in the independent model:

- E rho^2 = E Q = k/2
- E rho^4 = 3k^2/4 - 3k/8
- E Q rho^2 = k^2/4 + k/8
- E rho^6 = 5k/16 + 45k(k-1)/16 + 15k(k-1)(k-2)/8
- E Q rho^4 = k [5/16 + 15(k-1)/16 + 3(k-1)^2/8]

E[f w] = E rho^6 - E Q rho^4 - 2b (E rho^4 - E Q rho^2) + b^2 (E rho^2 - E Q). The last term vanishes. With
b = beta k/2, to leading order:

  E[f w] = (k/2)^3 (12 - 4 beta) + O(k^2),   E[w] = (k/2)^2 (3 - 2 beta + beta^2) + O(k).

The ratio (12 - 4 beta) / (3 - 2 beta + beta^2) is minimized at beta ~ 5.45, giving -0.4495. Divide by
sqrt|A| = sqrt(k(k-1)) to get the claim. ∎

## Theorem 2 (B_{2m+1}, general m)
Take w = q(rho^2)^2 with deg q = m. Then f w has degree 4m + 2, so B_{2m+1} suffices. The best weight gives
min f <= -(1 - lam_m + O(1/k)) * k/2, where lam_m = min_q E[z^2 q(z^2)^2] / E[q(z^2)^2] for a standard Gaussian z
(a Gaussian-quadrature / generalized eigenvalue problem). Computed:

| m | condition | dip >= (·) k/2 | score >= |
|---|---|---|---|
| 1 | B_3 | 0.4495 | 0.2247 |
| 2 | B_5 | 0.6197 | 0.3098 |
| 3 | B_7 | 0.7094 | 0.3547 |
| 4 | B_9 | 0.7648 | 0.3824 |
| 5 | B_11 | 0.8025 | 0.4013 |
| 6 | B_13 | 0.8298 | 0.4149 |
| 7 | B_15 | 0.8504 | 0.4252 |

lam_m -> 0 as m -> infinity (concentrating q at z = 0), so for B_h sets with h -> infinity the bound tends to **1/2**,
matching Bu et al.'s upper bound. **Square-type constructions over highly dissociated sets cannot beat 1/2.**

## Scope and open gap
- The constructions used in practice (Erdos-Turan, Singer) are Sidon (B_2) but not B_3, so the theorem does not apply to
  them directly. Their higher moments deviate from the Gaussian model; numerically they behave the same way (score > 1/2).
- Closing the gap (from B_h down to Sidon only) is the main open step. It is the natural place for a professional co-author.
- The O(1/k) terms are explicit polynomials and can be made fully explicit for each k.
