"""Map the true best score for moderate n (beyond Mercer's exact n <= 6): annealing with restarts over
A within [1, M], incremental f updates.   python smallsets.py n M SECONDS SEED"""
import json, math, os, sys, time
import numpy as np
import chowla as c

def run(n, M, seconds, seed):
    rng = np.random.default_rng(seed); t0 = time.time()
    G = 1 << int(math.ceil(math.log2(16 * M))); x = np.pi * np.arange(G + 1) / G
    C = np.cos(np.outer(np.arange(M + 1), x))
    best, bestA = 9e9, None
    while time.time() - t0 < seconds:
        A = list(rng.choice(np.arange(1, M + 1), n, replace=False)); inA = np.zeros(M + 1, bool); inA[A] = True
        f = C[A].sum(0); cur = f.min(); steps = int(4000 * n)
        for t in range(steps):
            T = 0.5 * math.sqrt(n) * (1 - t / steps) + 1e-3
            i = rng.integers(n); a = A[i]; b = int(rng.integers(1, M + 1))
            if inA[b]: continue
            f2 = f - C[a] + C[b]; v = f2.min()
            if v >= cur or rng.random() < math.exp((v - cur) / T):
                f, cur = f2, v; inA[a] = False; inA[b] = True; A[i] = b
        sc = -cur / math.sqrt(n)
        if sc < best: best, bestA = sc, sorted(int(t) for t in A)
    dense = -c.fmin(bestA, oversample=256) / math.sqrt(n)
    line = f"SMALL n={n} M={M} score={dense:.4f} A={bestA if n <= 24 else '...'}"
    print(line, flush=True)
    json.dump({"A": bestA, "dense": dense}, open(os.path.join(c.HERE, "runs", f"small_n{n}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
