# Literature check (2026-10-17, quick search, not exhaustive)

Normalized score = -min f_A / sqrt(|A|). Lower is better.

| Source | What it gives |
|---|---|
| Epoch problem page | c = 1 achievable since the 1960s; no example known for c <= 1/20 |
| Mercer, arXiv 1709.06612 | Exact/experimental optima for tiny n: n=2..6 -> 0.80, 0.76, 0.76, 0.73, 0.65. Folklore construction: lambda(n) <= sqrt(2n)+1/2 (constant ~1.41). No large-n construction below 1. |
| Bedert, arXiv 2509.05260 / 2509.03490 (2025) | Lower bound: min f_A <= -n^(1/7 - o(1)) |
| Shankar, arXiv 2609.05338 (Sep 2026) | Log-free lower bound K(N) >> N^(1/5). Only upper construction cited: Sidon differences, O(sqrt N), no explicit constant. No computer search. |

Takeaways
- Tiny optimal sets (n <= 6) already sit below 1, so "below 1" alone is not news; what matters is
  LARGE n. We found no published large-n sets below 1.
- Our LM rounding: 0.7519 at n=150 and 0.7596 at n=250 (dense-verified). If this holds at n in the
  thousands, it is (as far as this check shows) a new numerical record for large sets.
- Not checked yet: AI-search papers (AlphaEvolve-style, 2025-26) and Erdos-problems forum threads.
  Do that before claiming novelty.

## UPDATE (2026-10-17 afternoon) - our ~0.78 is BEATEN
- Classical Sidon-difference sets (Erdos-Turan B = {2pk + k^2 mod p}, A = positive differences) score
  sqrt(k/(2(k-1))) -> 1/sqrt(2) = 0.707. Verified here: p=41 n=820 -> 0.7159; p=47 n=1081 -> 0.7147.
  Epoch's "c = 1" was loose accounting.
- Bu, Khan, Yi et al., arXiv 2610.17160 (Oct 2026): proves limsup lambda(n)/sqrt(n) <= 1/2
  (Singer difference sets + Fejer majorant + parity deletion + random selection), numerics to n=25,440.
- Chowla is not in the AlphaEvolve problem repo; Epoch still lists it unsolved.
- Verdict: our LM sweeps (0.75-0.80) are not a record. Current bar = Bu et al. finite-n values (~0.5).
