"""Ideas 3-5 quick constructions near n=1560."""
import json, math, os, time
import numpy as np
import chowla as c
from bu_opt import build_A
def log(name, A):
    A = sorted(set(int(a) for a in A if a > 0)); line = f"IDEA {name} n={len(A)} score={c.score(A):.4f}"
    print(line, flush=True)
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")
A = json.load(open(os.path.join(c.HERE, "runs", "buopt_k40_p47_s7.json")))["A"]
log("baseline-bu", A)
for d in (1, 2, 5): log(f"3-cluster+{d}", A + [a + d for a in A])                       # local multiplicity
log("3-cluster+0,1,2", [a + d for a in A for d in (0, 1, 2)])
def bu_small(p):
    S = [2 * p * t + (t * t % p) for t in range(p)]; return build_A([s - S[p // 2] for s in S])
B = bu_small(7)                                                                          # n = 42
for L in (10 ** 4, 10 ** 5):
    log(f"4-union-dilate L={L}", A + [L * a for a in A])
    prod = [b1 + L * b2 for b1 in B for b2 in B] + [abs(b1 - L * b2) for b1 in B for b2 in B] + B + [L * b for b in B]
    log(f"4-product-gate L={L}", prod)
rng = np.random.default_rng(0); N = 3120                                                 # DPP, sine kernel, density 1/2
d = np.subtract.outer(np.arange(N), np.arange(N)); K = np.where(d == 0, 0.5, np.sin(np.pi * d / 2) / (np.pi * np.where(d == 0, 1, d)))
w, V = np.linalg.eigh(K)
for t in range(3):
    keep = rng.random(N) < w; Vk = V[:, keep]; pts = []
    while Vk.shape[1] > 0:
        pr = (Vk ** 2).sum(1); pr /= pr.sum(); i = rng.choice(N, p=pr); pts.append(i + 1)
        j = np.argmax(np.abs(Vk[i])); vj = Vk[:, j]; Vk = Vk - np.outer(vj, Vk[i] / vj[i]); Vk = np.delete(Vk, j, 1)
        if Vk.shape[1]: Vk, _ = np.linalg.qr(Vk)
    log(f"5-dpp-sine sample{t}", pts)
