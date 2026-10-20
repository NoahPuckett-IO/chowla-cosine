"""Puckett-set test: clustered Sidon sets S = {s0 + c : s0 in S0, c in C} plugged into Bu's identity
A = {3|a-b|} U {|3(a+b)+1|} (a<b in S), duplicates merged (a SET). Kurtosis of Psi rises with clustering;
merging collisions costs. Measures the net score."""
import math, sys
import numpy as np
import chowla as c
def kurt(S, G=1 << 16):
    x = np.pi * np.arange(G) / G; P = np.cos(np.outer(3 * np.array(S) + 0.5, x)).sum(0)
    return (P ** 4).mean() / (P ** 2).mean() ** 2
def build(S):
    A = set()
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            A.add(3 * abs(S[i] - S[j])); A.add(abs(3 * (S[i] + S[j]) + 1))
    return sorted(A)
for p in (23, 29, 41):
    S0 = [2 * p * t + (t * t % p) for t in range(p)]; S0 = [s - S0[p // 2] for s in S0]
    for C in ([0], [0, 1], [0, 2], [0, 7], [0, 1, 3], [0, 1, 2], [0, 5, 11], [0, 1, 2, 3]):
        S = sorted({s + d for s in S0 for d in C}); A = build(S)
        print(f"p={p} C={C} k={len(S)} n={len(A)} kurt={kurt(S):.2f} score={c.score(A):.4f}", flush=True)
