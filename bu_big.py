"""bu_opt for large n: same annealing over Sidon subsets inside Bu's identity, but cosines are computed
on the fly (no p x G tables), so memory is tiny and every core can run one.
  python bu_big.py K P SECONDS SEED"""
import json, math, os, sys, time
import numpy as np
import chowla as c
from bu_opt import build_A

def run(k, p, seconds, seed):
    rng = np.random.default_rng(seed); t0 = time.time()
    B = np.array([2 * p * t + (t * t % p) for t in range(p)]); B = B - B[p // 2]
    mx = 6 * int(np.abs(B).max()) + 2
    G = 1 << int(math.ceil(math.log2(3 * mx))); x = np.pi * np.arange(G + 1) / G
    R = lambda s: np.cos((3 * s + 0.5) * x); Gm = lambda s: np.cos((6 * s + 1) * x)
    inS = np.zeros(p, bool); inS[rng.choice(p, k, replace=False)] = True
    rho = sum(R(s) for s in B[inS]); g = sum(Gm(s) for s in B[inS])
    obj = lambda r, gg: (r * r - k / 2 - gg / 2).min()
    cur = obj(rho, g); best, bestS = cur, inS.copy(); it = 0
    while time.time() - t0 < seconds:
        it += 1
        T = 0.5 * (1 - (time.time() - t0) / seconds) + 1e-3
        a = rng.choice(np.flatnonzero(inS)); b = rng.choice(np.flatnonzero(~inS))
        r2 = rho - R(B[a]) + R(B[b]); g2 = g - Gm(B[a]) + Gm(B[b]); v = obj(r2, g2)
        if v >= cur or rng.random() < math.exp((v - cur) / T):
            rho, g, cur = r2, g2, v; inS[a], inS[b] = False, True
            if v > best: best, bestS = v, inS.copy()
    del rho, g, x
    S = [int(s) for s in B[bestS]]; A = build_A(S); n = len(A)
    dense = -c.fmin(A, oversample=16) / math.sqrt(n)
    line = f"BUBIG k{k} p{p} n{n} dense={dense:.4f} iters={it}"
    print(line, flush=True)
    json.dump({"k": k, "p": p, "S": S, "dense": dense}, open(os.path.join(c.HERE, "runs", f"bubig_k{k}_p{p}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
