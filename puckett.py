"""Hunt for PUCKETT SETS: Sidon S where cosines only cancel where they are also small, i.e.
Sum cos^2(theta_s) < (1/2 - delta) k wherever Sum cos(theta_s) ~ 0. Plugged into Bu's identity
f = rho^2 - Q, (rho = Sum cos theta_s, Q = Sum cos^2, theta_s = (3s+1/2)x) that beats the 1/2 wall.
Search: anneal over ALL Sidon sets S in [-M, M] (any element can move anywhere, Sidon kept).
  python puckett.py K M SECONDS SEED"""
import json, math, os, sys, time
import numpy as np
import chowla as c
from bu_opt import build_A

def run(k, M, seconds, seed):
    rng = np.random.default_rng(seed); t0 = time.time()
    G = 1 << int(math.ceil(math.log2(8 * (6 * M + 2)))); x = np.pi * np.arange(G + 1) / G
    S = []
    diffs = set()
    while len(S) < k:                                   # random greedy Sidon start
        s = int(rng.integers(-M, M + 1))
        nd = [abs(s - t) for t in S]
        if s in S or len(set(nd)) < len(nd) or any(d in diffs for d in nd): continue
        S.append(s); diffs.update(nd)
    R = {s: np.cos((3 * s + 0.5) * x) for s in S}
    rho = sum(R.values()); Q = sum(r * r for r in R.values())
    cur = (rho * rho - Q).min(); best, bestS = cur, list(S); it = 0
    while time.time() - t0 < seconds:
        it += 1
        T = 1.5 * (1 - (time.time() - t0) / seconds) + 0.02
        out = S[rng.integers(k)]; rest = [t for t in S if t != out]
        dr = {abs(a - b) for i, a in enumerate(rest) for b in rest[i + 1:]} if it % 50 == 0 or True else None
        new = int(rng.integers(-M, M + 1))
        if new in rest: continue
        nd = [abs(new - t) for t in rest]
        if len(set(nd)) < len(nd) or any(d in dr for d in nd): continue
        rn = np.cos((3 * new + 0.5) * x)
        rho2 = rho - R[out] + rn; Q2 = Q - R[out] ** 2 + rn * rn
        v = (rho2 * rho2 - Q2).min()
        if v >= cur or rng.random() < math.exp((v - cur) / T):
            rho, Q, cur = rho2, Q2, v; S = rest + [new]; del R[out]; R[new] = rn
            if v > best: best, bestS = v, list(S)
    A = build_A(bestS); n = len(A)
    sc = c.score(A); dense = -c.fmin(A, oversample=64) / math.sqrt(n)
    line = f"PUCKETT k={k} M={M} n={n} score={sc:.4f} dense={dense:.4f} iters={it}"
    print(line, flush=True)
    json.dump({"S": sorted(bestS), "A": A, "dense": dense}, open(os.path.join(c.HERE, "runs", f"puckett_k{k}_M{M}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
