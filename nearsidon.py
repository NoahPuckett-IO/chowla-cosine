"""Idea 2: near-Sidon S (3-APs / collisions allowed) in Bu's identity, scored on the TRUE merged set A.
Frequency counts tracked so f = sum over distinct frequencies.   python nearsidon.py K M SECONDS SEED"""
import json, math, os, sys, time
import numpy as np
import chowla as c
def freqs(s, rest): return [3 * abs(s - t) for t in rest] + [abs(3 * (s + t) + 1) for t in rest]
def run(k, M, seconds, seed):
    rng = np.random.default_rng(seed); t0 = time.time()
    G = 1 << int(math.ceil(math.log2(8 * (6 * M + 2)))); x = np.pi * np.arange(G + 1) / G
    S = [int(v) for v in rng.choice(np.arange(-M, M + 1), k, replace=False)]
    cnt = {}; f = np.zeros(G + 1)
    def add(fr, sgn):
        nonlocal f
        for q in fr:
            if q == 0: continue
            old = cnt.get(q, 0); new = old + sgn; cnt[q] = new
            if old == 0 and new == 1: f += np.cos(q * x)
            elif old == 1 and new == 0: f -= np.cos(q * x)
    for i, s in enumerate(S): add(freqs(s, S[:i]), 1)
    n = sum(1 for v in cnt.values() if v > 0); cur = -f.min() / math.sqrt(n); best, bestS = cur, list(S)
    while time.time() - t0 < seconds:
        T = 0.05 * (1 - (time.time() - t0) / seconds) + 1e-4
        i = rng.integers(k); out = S[i]; rest = S[:i] + S[i + 1:]; new = int(rng.integers(-M, M + 1))
        if new in S: continue
        f0 = f.copy(); c0 = dict(cnt)
        add(freqs(out, rest), -1); add(freqs(new, rest), 1)
        n2 = sum(1 for v in cnt.values() if v > 0); v = -f.min() / math.sqrt(n2)
        if v <= cur or rng.random() < math.exp((cur - v) / T):
            cur = v; S = rest + [new]
            if v < best: best, bestS = v, list(S)
        else: f, cnt = f0, c0
    A = sorted({q for i, s in enumerate(bestS) for q in freqs(s, bestS[:i]) if q > 0}); n = len(A)
    dense = -c.fmin(A, oversample=32) / math.sqrt(n)
    line = f"NEARSIDON k={k} M={M} n={n} dense={dense:.4f}"
    print(line, flush=True)
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")
if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
