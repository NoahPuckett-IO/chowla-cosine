"""Concrete check of Theorem 1: random greedy B_3 sets S (k = 10..24), build A(S), measure the score.
The theorem says score >= 0.2247 - O(1/k); we expect ~1/2 or more in practice."""
import itertools, math, os, time
import numpy as np
import chowla as c
from bu_opt import build_A
rng = np.random.default_rng(0)
for k in (10, 14, 18, 22, 26):
    for trial in range(2):
        S, sums3 = [], set()
        while len(S) < k:
            s = int(rng.integers(-k ** 3 * 4, k ** 3 * 4))
            if s in S: continue
            cand = S + [s]; new = set()
            ok = True
            for t in itertools.combinations_with_replacement(cand, 3):
                if s not in t: continue
                v = sum(t)
                if v in sums3 or v in new: ok = False; break
                new.add(v)
            if ok:
                # also need B_2: pair sums distinct
                ps = [a + b for a, b in itertools.combinations_with_replacement(cand, 2)]
                if len(ps) == len(set(ps)): S, sums3 = cand, sums3 | new
        S = [2 * s for s in S]; A = build_A(S); sc = -c.fmin(A, oversample=32) / math.sqrt(len(A))
        line = f"B3CHECK2 d=2 k={k} n={len(A)} score={sc:.4f} (theorem: >= 0.2247 - O(1/k))"
        print(line, flush=True)
        with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")
