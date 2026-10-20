"""Optimize the Sidon set S inside Bu et al.'s two-sided identity instead of sampling it.
For Sidon S (size k): A = {3|si-sj|} U {|3(si+sj)+1|}, |A| = k(k-1), and
  f_A(x) = rho(x)^2 - k/2 - g(x)/2,  rho = sum cos((3s+1/2)x),  g = sum cos((6s+1)x).
Pool = Erdos-Turan Sidon set {2pk + (k^2 mod p)}, centred. Annealing swaps one element of S with an
unused pool element; rho and g update in O(grid) per move.
  python bu_opt.py K P SECONDS SEED
"""
import json, math, os, sys, time
import numpy as np
import chowla as c

def build_A(S):
    S = sorted(S); A = set()
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            A.add(3 * abs(S[i] - S[j])); A.add(abs(3 * (S[i] + S[j]) + 1))
    return sorted(A)

def run(k, p, seconds, seed):
    rng = np.random.default_rng(seed)
    B = np.array([2 * p * t + (t * t % p) for t in range(p)]); B = B - B[p // 2]
    mx = 6 * int(np.abs(B).max()) + 2
    G = 1 << int(math.ceil(math.log2(6 * mx)))
    x = np.pi * np.arange(G + 1) / G                      # f even, period 2pi -> [0, pi]
    R = np.empty((p, len(x)), np.float32); Gm = np.empty((p, len(x)), np.float32)
    for t in range(p): R[t] = np.cos((3 * B[t] + 0.5) * x); Gm[t] = np.cos((6 * B[t] + 1) * x)
    idx = rng.choice(p, k, replace=False); inS = np.zeros(p, bool); inS[idx] = True
    rho = R[inS].sum(0); g = Gm[inS].sum(0)
    obj = lambda r, gg: (r * r - k / 2 - gg / 2).min()
    cur = obj(rho, g); best, bestS = cur, inS.copy(); t0 = time.time(); it = 0
    while time.time() - t0 < seconds:
        it += 1
        T = 0.5 * (1 - (time.time() - t0) / seconds) + 1e-3
        a = rng.choice(np.flatnonzero(inS)); b = rng.choice(np.flatnonzero(~inS))
        r2 = rho - R[a] + R[b]; g2 = g - Gm[a] + Gm[b]; v = obj(r2, g2)
        if v >= cur or rng.random() < math.exp((v - cur) / T):
            rho, g, cur = r2, g2, v; inS[a], inS[b] = False, True
            if v > best: best, bestS = v, inS.copy()
    S = [int(s) for s in B[bestS]]; A = build_A(S); n = len(A)
    sc = c.score(A); dense = -c.fmin(A, oversample=64) / math.sqrt(n)
    line = f"BUOPT k{k} p{p} n{n} grid={-best/math.sqrt(n):.4f} score={sc:.4f} dense={dense:.4f} iters={it}"
    print(line, flush=True)
    os.makedirs(os.path.join(c.HERE, "runs"), exist_ok=True)
    json.dump({"k": k, "p": p, "S": S, "A": A, "score": sc, "dense": dense},
              open(os.path.join(c.HERE, "runs", f"buopt_k{k}_p{p}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh:
        fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    k, p, sec, seed = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    run(k, p, sec, seed)
