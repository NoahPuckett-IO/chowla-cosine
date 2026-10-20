# Chowla findings (2026-10-17)

Score = -min f_A / sqrt(|A|), lower is better. All "ours" values dense-verified (64x grid).

## Headline: optimized Sidon selection beats Bu et al. (arXiv 2610.17160) Table 1
Method (bu_opt.py): keep Bu et al.'s two-sided identity f_A = rho^2 - k/2 - g/2, but instead of a random
Bernoulli selection of S, run simulated annealing over k-subsets S of an Erdos-Turan Sidon pool
{2pt + (t^2 mod p)} (centred), maximizing min f on a grid with O(grid) incremental updates. 40 min/job.

| n | Bu et al. | Sidon diffs | ours | k | pool p |
|---|---|---|---|---|---|
| 1560 | 0.725 | ~0.713 | 0.6461 | 40 | 47 |
| 3540 | 0.708 | ~0.711 | 0.6264 | 60 | 67 |
| 6320 | 0.674 | ~0.710 | 0.6167 (Mac) | 80 | 89 |
| 14280 | 0.655 | ~0.709 | 0.6042 | 120 | 131 |
| 25440 | 0.635 | ~0.708 | 0.7744 (lost) | 160 | 163 |

Lessons: pool p ~ 1.2-1.4 k works; p ~ k is too tight (n=25440 loss, also k60 p61, k80 p83).
Sets saved: runs/buopt_k*_p*_s1.json (S and A).

## Earlier (superseded)
- LM rounding of Fejer weights: ~0.75-0.80 for n <= 1000. Beaten by classical Sidon (0.707). See LITERATURE.md.

## Caveats / next
- Numerical only, finite n; not a proof. Bu's asymptotic 1/2 still stands; our values are above 1/2.
- Next: k=160 with p ~ 200-230; shift parameter; Singer-set pools; polish A directly.

## Round 2 combined best (PC + Mac, dense-verified), 2026-10-17 ~17:15
| n | 1560 | 2450 | 3540 | 4830 | 6320 | 8010 | 9900 | 11990 | 14280 | 25440 |
|---|---|---|---|---|---|---|---|---|---|---|
| ours | 0.6450 | 0.6395 | 0.6242 | 0.6210 | 0.6140 | 0.6072 | 0.6050 | 0.5975 | 0.5984 | 0.5897 |
| Bu | 0.725 | - | 0.708 | - | 0.674 | - | - | - | 0.655 | 0.635 |
Beats Bu at every n in their table, incl. n=25440 (bigger pools fixed it). Trend: slowly falling, ~0.59 at n=25k, still above 1/2.

## Overnight 2026-10-18

| n | Our Innovative Try (CORR) | Our Old Try (BUOPT) | Bu et al.'s Try |
|---|---|---|---|
| 1560 | 0.6750 @ n=1576 (n0=1560); also 0.8080 @ 1638, 0.8232 @ 1607, 1.0828 @ 1872 | 0.6426 | 0.725 |
| 2450 | — | 0.6339 | — |
| 3540 | — | 0.6242 | 0.708 |
| 4830 | — | 0.6156 | — |
| 6320 | — | 0.6131 | 0.674 |
| 8010 | — | 0.6068 | — |
| 9900 | — | 0.6050 | — |
| 11990 | — | 0.5975 | — |
| 14280 | — | 0.5984 | 0.655 |
| 25440 | — | 0.5897 | 0.635 |

Verdict: No CORR score is below 0.5 (best 0.6750 at n=1576, from n0=1560); CORR runs exist only near n≈1560 and scores get worse as mfrac grows.
CORR does not beat BUOPT at comparable n (0.6750 vs 0.6426 at n≈1560); BUOPT still beats Bu et al. at every n they report.

## Overnight 2026-10-18
| n | Our Innovative Try (combo) | Our Old Try (bu_opt) | Bu et al.'s Try |
|---|---|---|---|
| ~1,560 | 0.675 (n=1576, +1%) | 0.643 | 0.725 |
| 2,450 | - | 0.634 | - |
| 3,540 | - | 0.624 | 0.708 |
| 6,320 | - | 0.613 | 0.674 |
| 14,280 | - | 0.598 | 0.655 |
| 25,440 | - | 0.590 | 0.635 |
Combo with +3/5/10/20% extra numbers: 0.81-1.08 (worse). Verdict: the combo did NOT work as built.
The LP's fractional gain (0.43 on its active points) was lost: cutting-plane rounds too slow and
pipage rounding too crude. No score below 0.5. Old Try improved slightly (n=1560 0.643, 6320 0.613).

## 2026-10-18 morning: attempts to break 1/2 (none did)
- Melt-down (continuous relaxation + crystallize), from Bu set: n=1560 -> 0.6315 (new best at that n). Fejer start: >=1.04.
- Breakout swaps (24 variants): Bu set is a strict local optimum for single moves (0.6426 unchanged).
- Puckett hunt (anneal over all Sidon sets): 0.644 (n=2256), 0.655 (n=1560). Dense ramp sets: 0.93-0.97.
- Small-n map (anneal, A within [1,3n]): n=6 0.650 (= Mercer exact), 8 0.692, 10 0.651, 12 0.639, 14 0.620, 16 0.648,
  18 0.646, 20 0.702, 24 0.686, 32 0.813, 48 0.857, 100 0.984, 200 1.085 (search weak for n>=20).
- Theory note: realizable constructions give Hermite polynomials in z = rho/sigma (z^2-1 -> 1/2 wall; z^4-6z^2+3 -> 0.87).
  The one-sided shape z^4-3 (~0.22) needs the degree-2 block with multiplicity ~k, which a set cannot have.
- Kurtosis reformulation: for f = Psi^2 - mean(Psi^2), score = 1/sqrt(2(kappa-1)), kappa = kurtosis of Psi.
  Sidon <=> minimal additive energy <=> kappa ~ 3 <=> 1/2. Sub-1/2 needs kappa > 3 WITHOUT sum/difference collisions.
  Test (kurt.py): clustered Sidon S0+{0,d,...} raises kappa to 6-17 but merged collisions push scores to 2.1-8.5 (vs 0.77).
- Re-melt cycle (from 0.6315 set): n=1560 -> 0.5944 (record). n=6320 0.6191, n=14280 1.41 (melt fails at large L).

## 2026-10-18 noon: five alternative routes, all ruled out at n~1560 (bar: Bu 0.6426, record 0.5944)
| Idea | Best | Note |
|---|---|---|
| 1 one-sided rounding of Fejer ramp (shape.py) | 0.93-0.97 | random ramp 1.50-1.84 |
| 2 near-Sidon, collisions allowed, true merged set (nearsidon.py) | 0.755 | most 0.85-1.18 |
| 3 local multiplicity (clusters A+{0,d}) | 2.28-3.37 | |
| 4 multi-scale gating (union/product at scale L) | ~0.91 analytic | dips align; run OOM |
| 5 repulsive sets (DPP; jittered lattice proxy) | 8.59 | true DPP sampler too slow |
Re-melt cycle 2 from 0.5944: best 0.6035 (plateau).
