"""Dense ramp ('shape') sets: start from a Fejer-ramp rounding of [1, 2n] (about half the numbers, thinning
linearly), then dip-targeted soft-min swaps (breakout) to push the leftover error into a one-sided shape.
  python shape.py n SECONDS SEED TAU MODE"""
import json, os, sys
import chowla as c, breakout
n, sec, seed, tau, mode = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
A = c.pipage_round(2 * n, seed=seed)
p = os.path.join(c.HERE, "runs", f"shape_start_n{n}_s{seed}.json"); json.dump({"A": A}, open(p, "w"))
print("start score", round(c.score(A), 4), len(A), flush=True)
breakout.main(p, sec, seed, tau, mode)
