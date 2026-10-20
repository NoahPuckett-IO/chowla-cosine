"""One sweep job: python sweep_lm.py K. Tries eta/tc/seeds, verifies best with a dense grid, saves it."""
import json, os, sys, time, chowla as c
K = int(sys.argv[1]); best = (9e9, None, None)
for eta in (0.2, 0.3, 0.5):
    for tc in (0.3, 0.4, 0.5):
        for s in range(3):
            A = c.lm_round(K, tc=tc, gamma=0.1, eta=eta, seed=s)
            sc = c.score(A)
            if sc < best[0]: best = (sc, A, (eta, tc, s))
sc, A, prm = best
dense = -c.fmin(A, oversample=64) / len(A) ** 0.5
json.dump({"K": K, "n": len(A), "score": sc, "dense": dense, "params": prm, "A": A},
          open(os.path.join(c.HERE, "runs", f"lm_K{K}.json"), "w"))
with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh:
    fh.write(f"{time.strftime('%F %T')} LM K{K} n{len(A)} score={sc:.4f} dense={dense:.4f} params={prm}\n")
