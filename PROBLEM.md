# Chowla's Cosine Problem (Epoch FrontierMath: Open Problems)

Source: https://epoch.ai/frontiermath/open-problems/chowla-cosine
Status (2026-10-17): Unsolved. Tier "Moderately interesting". Category: Analysis, Construction - Finite, Improved Bound.

## Statement (from Epoch's prompt)
For a finite SET A of positive integers, f_A(x) = sum_{a in A} cos(a x).
Average of f_A is 0 and f_A(0) = |A|, so f_A dips below 0. How negative must it be, relative to sqrt(|A|)?

Task: write `solution(c: float, n: int) -> str` that, for EVERY c > 0 and n, returns a set A with
|A| >= n and f_A(x) >= -c*sqrt(|A|) for all real x.
Output: distinct positive integers, whitespace/comma separated. Keep max(A) well under 10^6
(lower ceiling for small c). Verifier = numerical sampling check.

## Known state
- Constructions since the 1960s achieve the bound for c = 1.
- NO known example for c <= 1/20.
- Solving it fully = DISPROVING Chowla's conjecture (that -sqrt(|A|) is the right order).
  Epoch notes the conjecture may be true, in which case the problem is unsolvable.
- Recent lower-bound papers (2025): arXiv 2509.05260, arXiv 2509.03490 (polynomial lower bounds
  on how negative f_A must be; they do NOT rule out c -> 0 at the sqrt scale ... check this first).

## Reality check
This is not "beat a number". A real solve needs a new construction idea that works for all c.
Realistic wins, in order:
1. Explicit large sets with c well below 1 (e.g. 0.5, 0.2) -> notable partial progress, reportable.
2. A pattern in optimized sets that generalizes -> candidate algorithm.
3. Full algorithm for every c -> the actual solve (would be a big deal).
