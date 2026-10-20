"""Break out of the Bu family: dip-targeted swaps on the WHOLE set A.
Worst grid point x*; remove a in A with cos(a x*) most negative, add b not in A with cos(b x*) most
positive; keep a move if the global min improves (small annealing allowance).
  python breakout.py START.json SECONDS SEED"""
import json, math, os, sys, time
import numpy as np
import chowla as c

def main(start, seconds, seed, tau=1.5, mode='swap'):
    rng = np.random.default_rng(seed); t0 = time.time()
    A = sorted(json.load(open(start))["A"]); n = len(A); L = int(1.3 * max(A))
    G = 1 << int(math.ceil(math.log2(8 * L))); x = np.pi * np.arange(G + 1) / G
    inA = np.zeros(L + 1, bool); inA[A] = True
    ind = np.zeros(2 * G); ind[A] = 1; f = np.fft.rfft(ind).real[:G + 1]
    sm = lambda v: v.min() - tau * np.log(np.exp(-(v - v.min()) / tau).sum())
    scr = lambda v, m: -v.min() / math.sqrt(m)
    m = n; cur = sm(f); best = scr(f, m); bestA = inA.copy(); it = 0
    allf = np.arange(1, L + 1)
    while time.time() - t0 < seconds:
        it += 1
        js = np.argsort(f)[:3]                           # current worst points
        j = js[rng.integers(3)]
        cw = np.cos(allf * x[j])
        ins = allf[inA[1:]]; outs = allf[~inA[1:]]
        ca = cw[ins - 1]; cb = cw[outs - 1]
        a_c = ins[np.argsort(ca)[:20]]; b_c = outs[np.argsort(-cb)[:20]]
        a = int(a_c[rng.integers(20)]); b = int(b_c[rng.integers(20)])
        mv = mode if mode != 'mix' else ('swap', 'add', 'rem')[rng.integers(3)]
        if mv == 'grow': mv = 'add' if rng.random() < 0.7 else 'swap'
        if mv == 'swap': f2 = f - np.cos(a * x) + np.cos(b * x); m2 = m
        elif mv == 'add': f2 = f + np.cos(b * x); m2 = m + 1
        else: f2 = f - np.cos(a * x); m2 = m - 1
        v = sm(f2) / math.sqrt(m2) * math.sqrt(m)
        if v > cur - 0.05 * rng.random():
            f, cur, m = f2, sm(f2), m2
            if mv != 'add': inA[a] = False
            if mv != 'rem': inA[b] = True
            if scr(f, m) < best: best, bestA = scr(f, m), inA.copy()
    A = [int(t) for t in np.flatnonzero(bestA)]
    sc = c.score(A); dense = -c.fmin(A, oversample=64) / math.sqrt(len(A))
    line = f"BREAK {mode} start={os.path.basename(start)} tau={tau} n={len(A)} start={os.path.basename(start)} score={sc:.4f} dense={dense:.4f} iters={it}"
    print(line, flush=True)
    json.dump({"A": A, "dense": dense}, open(os.path.join(c.HERE, "runs", f"break_n{len(A)}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
